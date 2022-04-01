# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
import pandas as pd
import os, re
import functools
import shutil
from biocluster.config import Config
from bson.objectid import ObjectId
import types
from biocluster.api.database.base import Base
from collections import defaultdict

# 文件移动、运行时间相关方法

def clean_mongo(search={}, targets=[], all=False):
    proj = "metabolome"
    db_client = Config().get_mongo_client(mtype=proj)
    db = db_client[Config().get_mongo_dbname(proj)]
    col = db["sg_table_relation"]
    main_tables = defaultdict(set)
    for one in col.find_one()["target"]:
        if one[0] in targets or all:
            main_ids = remove_detail(db, search, *one)
            main_tables[one[0]].update(main_ids)

    for m in main_tables:
        remove_main(db, m, list(main_tables[m]))


def remove_detail(db, search, main_col, detail_col, main_key):
    print("开始删除 {} 详情表 {}".format(detail_col, search))
    main_ids = set()
    for one in db[main_col].find(search):
        if "main_id" in one:
            main_ids.add(one["main_id"])
        else:
            main_ids.add(one["_id"])
    db[detail_col].remove({main_key: {"$in": list(main_ids)}})
    print("完成删除 {} 详情表内容 {}".format(detail_col, main_ids))
    return main_ids

def remove_main(db, main_col, main_ids):
    print("开始删除 {} 主表 {}".format(main_col, main_ids))
    db[main_col].remove({"main_id": {"$in": main_ids}})
    db[main_col].remove({"_id": {"$in": main_ids}})
    db["sg_status"].update({"table_id": {"$in": main_ids}}, {"$set": {"status": "delete_relation"}})
    print("完成删除 {} 主表".format(main_col, main_ids))

def link_file(oldfile, newfile):
    """
    hard link file from oldfile to newfile
    :param oldfile:
    :param newfile:
    :return:
    """
    if not os.path.isfile(oldfile):
        raise Exception("不存在文件：%s" % oldfile)
    if os.path.exists(newfile):
        os.remove(newfile)
    os.link(oldfile, newfile)

def link_dir(olddir, newdir):
    """
    hard link directory from olddir to newdir
    :param olddir:
    :param newdir:
    :return:
    """
    if not os.path.isdir(olddir):
        raise Exception("不存在路径: %s" % olddir)
    allfiles = os.listdir(olddir)
    if not os.path.exists(newdir):
        os.makedirs(newdir)
    # else:
    #     shutil.rmtree(newdir)
    #     os.makedirs(newdir)
    for file in allfiles:
        oldfile = os.path.join(olddir, file)
        newfile = os.path.join(newdir, file)
        if os.path.isdir(newfile):
            shutil.rmtree(newfile)
        if os.path.isfile(oldfile):
            link_file(oldfile, newfile)
        else:
            link_dir(oldfile, newfile)

def time_count(func):  # 统计函数运行时间，作为方法的装饰器
    @functools.wraps(func)
    def wrapper(*args, **kw):
        start = time.time()
        func_name = func.__name__
        start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))
        print('Run %s at %s' % (func_name, start_time))
        func(*args, **kw)
        end = time.time()
        end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end))
        print('End %s at %s' % (func_name, end_time))
        print("{}函数执行完毕，共运行{}s".format(func_name, end - start))
    return wrapper

def wait_file(path, times=10, sleep=10):
    '''
    等待某个文件生成
    :param path: 文件名称
    :param wait_times: 最大等待次数
    :param sleep: 每次等待时间，单位秒
    :return:
    '''

    while times > 0:
        if not os.path.isfile(path):
            time.sleep(sleep)
            times -= 1
            self.wait_file(path, times=times, sleep=sleep)
        return path
    raise Exception("超过文件等待次数，需检查文件%s" % path)

# 代谢组处理相关方法

def check_info(info):
    '''
    检查原始表中的注释内容，将代表无注释的字符统一处理成'-'
    :param info: 注释信息
    :return:  str
    '''
    if info in ['_', '-', ' ', '', '\r']:
        return '-'
    else:
        return info

def check_metab(info):
    '''
    检查代谢物名称中是否有特殊字符，返回修改后的代谢物名称
    :param info: 代谢物或其他有特殊字符的名称
    :return: info
    '''
    info = check_info(info)
    # if info == '-':
    #     return info
    # info = info.replace("α", "|alpha|").replace("β", "|beta|").replace("ω", "|omega|").replace("->", "|right|")
    return info

def check_metab_list(info_list):
    """
    对一个代谢物列表进行处理
    :param info_list:
    :return:
    """
    new_list = []
    for info in info_list:
        new_info = check_metab(info)
        new_list.append(new_info)
    return new_list

def check_head(table, head_name):
    """
    查看head_name是否在table的表头中出现
    :param table: 数据表
    :param head_name: 查看的表头名称
    :return: bool
    """
    data = pd.read_table(table)
    try:
        test = data[head_name]
        return True
    except:
        return False

def check_index(table, index_name, index_col_name, allowzero=True, getid=True):
    """
    查看行名是否在table的行名中出现
    :param table: 数据表
    :param index_name: 查看的行名称
    :param index_col_name: 查看的列名称
    :param allowzero: 是否允许含有0
    :return:bool
    """
    data = pd.read_table(table, index_col=index_col_name)
    try:
        if 0 in data.loc[index_name] and not allowzero:
            return False
        else:
            if getid:
                return data['metab_id'][index_name]
            return True
    except:
        return False

def combine_table(metab_table, output_table, metabset=None, group=None, ignore_control=False):
    '''
    根据代谢集提取代谢物，根据分组表提取样品
    分组表中可能含有Control组的样品，可以通过ignore_control判断是否获取Control组的数据
    :param metab_table: metab_abun 或metab_desc以及其他已metab_id列开头的数据表，使用group参数时，必须包含group中的所有样本
    :param output_table: 输出新的结果表
    :param metabset: 代谢集文件，不能为mul_metabset格式
    :param group: 分组文件，两列（只能包含一种分组），表头可自定义，但必须要有
    :param ignore_control: 是否保留control样品，默认保留
    :return:
    '''
    if type(ignore_control) != bool:
        raise Exception("ignore_control param must be True or False")
    if not os.path.isfile(metab_table):
        raise Exception("no file exists: %s" % metab_table)
    if not os.path.exists(os.path.dirname(output_table)):
        os.makedirs(os.path.dirname(output_table))
    table_data = pd.read_table(metab_table, index_col='metab_id')
    if metabset:
        if not os.path.isfile(metabset):
            raise Exception("no file exists: %s" % metabset)
        set_data = pd.read_table(metabset, header=None, index_col=0)
        table_data = pd.merge(table_data, set_data, right_index=True, left_index=True, how='inner')
        if len(table_data) == 0:
            raise Exception("no metab_id in both metab_table and metaset")
    if group:
        if not os.path.isfile(group):
            raise Exception("no file exists: %s" % group)
        group_data = pd.read_table(group)
        group_data.index = group_data[group_data.columns[0]]
        if ignore_control:
            group_name = group_data.columns[1]
            group_set = set(group_data[group_name])
            if "Control" in group_set:
                group_set.remove("Control")
                group_data = group_data[group_data[group_name].isin(group_set)]
        table_data = table_data[group_data.index]
    table_data.index.name = 'metab_id'
    table_data.to_csv(output_table, sep="\t", quoting=3)


class Relation(Base):
    '''
    提供id与详情互转功能： convert_to_id/convert_id
    id与详情map文件生成:  save_map
    id与详情对应关系变量: return_dic {info: matb_id}
                          return_id_dic {metab_id: info}
    '''

    def __init__(self):
        """

        :param mongo: 是否从mongo里取数据
        :return:
        """
        super(Relation, self).__init__()
        self.dic_name_list = []  # [列名, 列名, 列名]
        self.dic_info = {}  # {列名: {该列详情: id}}
        self.id_info = {}  # {列名：{id: 该列详情}}
        self._project_type = "metabolome"

    def get_mongo_info(self, coll, main_name, main_id):
        """
        从mongo库获取数据
        :param coll: collection名称
        :param main_name: 主表字段
        :param main_id: 主表_id
        :return:
        """
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                raise Exception('main_id必须为ObjectId对象或其对应的字符串！')
        result = self.db[coll].find({main_name: main_id})
        print "从mongo库读取%s个记录" % result.count()
        name_list = []
        for one in result:
            for name in one.keys():
                if name in ["_id", "metab_id"]:
                    continue
                name_list.append(name)
                if self.dic_info.has_key(name):
                    self.dic_info[name][one[name]] = one["metab_id"]
                    self.id_info[name][one['metab_id']] = one[name]
                else:
                    self.dic_info[name] = {one[name]: one["metab_id"]}
                    self.id_info[name] = {one["metab_id"]: one[name]}
        self.dic_name_list = list(set(name_list))
        print self.dic_name_list

    def get_table_info(self, table, id_name=None):
        """
        从本地表格读取各列与id对应关系
        :param table: 表格路径
        :return:
        """
        if id_name:
            id_name =id_name
        else:
            id_name = "metab_id"
        data = pd.read_table(table) #pd.read_table
        data.drop_duplicates([id_name],keep='first',inplace=True)
        data = data.set_index([id_name])
        self.dic_name_list = data.columns.tolist()
        print self.dic_name_list
        for i in data.index.tolist():
            for col_name in self.dic_name_list:
                if self.dic_info.has_key(col_name):
                    self.dic_info[col_name][data[col_name][i]] = i
                    self.id_info[col_name][i] = data[col_name][i]
                else:
                    self.dic_info[col_name] = {data[col_name][i]: i}
                    self.id_info[col_name] = {i: data[col_name][i]}

    def save_map(self, dic_name, outfile):
        """
        将读取到的数据存成本地文件，两列，"metab_id\tdic_name"
        :param dic_name: 此详情列在读取的表格或mongo库中的名称
        :param outfile: 输出的文件名称
        :return:
        """
        if dic_name not in self.dic_name_list:
            raise Exception("%s不在提取的内容中，提取的内容：%s" % (dic_name, self.dic_name_list))
        with open(outfile, 'w') as f:
            f.write("metab_id\t" + dic_name + "\n")
            for metab_id in self.id_info[dic_name].keys():
                f.write("%s\t%s\n" % (metab_id, self.id_info[dic_name][metab_id]))

    def return_dic(self, dic_name):
        """
        返回一个详情数据{metab_detail: metab_id}
        :param dic_name: 此详情数据在mongo库或详情表中的名称
        :return: diction
        """
        if type(dic_name) == list:
            tmp_list = []
            for i in dic_name:
                if i in self.dic_name_list:
                    tmp_list.append(self.dic_info[i])
                else:
                    raise Exception("%s不在提取的内容中，提取的内容：%s" % (i, self.dic_name_list))
            return tmp_list
        else:
            if dic_name in self.dic_name_list:
                return self.dic_info[dic_name]
            else:
                raise Exception("%s不在提取的内容中， 提取的内容：%s" % (dic_name, self.dic_name_list))

    def return_id_dic(self, dic_name):
        """
        返回一个详情数据{metab_id: metab_detail}
        :param dic_name: 此详情数据在mongo库或详情表中的名称
        :return:diction
        """
        if type(dic_name) == list:
            tmp_list = []
            for i in dic_name:
                if i in self.dic_name_list:
                    tmp_list.append(self.id_info[i])
                else:
                    raise Exception("%s不在提取的内容中，提取的内容：%s" % (i, self.dic_name_list))
            return tmp_list
        else:
            if dic_name in self.dic_name_list:
                return self.id_info[dic_name]
            else:
                raise Exception("%s不在提取的内容中， 提取的内容：%s" % (dic_name, self.dic_name_list))

    def convert_to_id(self, dic_name, info):
        '''
        将详情转换成id
        :param dic_name:  详情所在的字段
        :param info:  详情信息
        :return:
        '''
        if dic_name not in self.dic_name_list:
            raise Exception("%s不在提取的内容中，提取的内容：%s" % (dic_name, self.dic_name_list))
        if self.dic_info[dic_name].has_key(info):
            return self.dic_info[dic_name][info]
        else:
            raise Exception("%s在读取的数据中未发现，检查%s内容" % (info, dic_name))

    def convert_id(self, dic_name, metab_id):
        '''
        将id转换成详情
        :param dic_name:  详情所在的字段
        :param metab_id: 代谢物id
        :return:
        '''
        if dic_name not in self.dic_name_list:
            raise Exception("%s不在提取的内容中，提取的内容：%s" % (dic_name, self.dic_name_list))
        if self.id_info[dic_name].has_key(metab_id):
            return self.id_info[dic_name][metab_id]
        else:
            raise Exception("未发现metabid：%s" % metab_id)

    def get_trans_file(self, map_table, id_name_dict, oldfile, newfile, myrow=True, mycol=False,add_id=False,drop_metab_id=True):  # add by shaohua.yuan
        """
        根据map_file中metab_id,替换为metab名称,生成新file,map_file和old_file第一列内容为metab_id
        """
        #map_file = "id_map.txt"
        #metab_trans.Metabolite("Metabolite")
        #metab_trans.save_map("Metabolite", "id_map.txt")
        #map_table = pd.read_table(map_file, sep="\t", header=0)
        #map_table.rename(columns={map_table.columns[0]: "metab_id"}, inplace=True)
        old_table = pd.read_table(oldfile, sep="\t", header=0)
        old_table.rename(columns={old_table.columns[0]: "metab_id"}, inplace=True)
        merge = old_table
        if myrow:
            tmp_merge = pd.merge(merge, map_table, how="left", on="metab_id")
            columns_order = ["Metabolite"] + merge.columns.tolist()
            merge = tmp_merge[columns_order]
            print "row"
            print merge.head()
        if mycol:
            new_columns = []
            id_name_dict["metab_id"] = "metab_id"
            id_name_dict["Metabolite"] = "Metabolite"
            [new_columns.append(id_name_dict[x]) for x in merge.columns]
            merge.columns = [new_columns]   # 20190729
        if drop_metab_id:
            def _tmp_fun(x):
                if pd.isna(x['Metabolite']):
                    x['Metabolite'] = x['metab_id']
                return x
            #merge = merge.apply(_tmp_fun,axis=1)
            merge = merge.drop('metab_id', 1)
            #merge = merge.drop(['metab_id'], axis = 1)
        merge.to_csv(newfile, sep="\t", index=False, quoting=3)

    def get_dataframe_and_dict(self, map_path, mode=1, metab_name=None,id_name=None):  ## add by shaohua.yuan
        '''
        返回Metabolite,metab_id的dataframe
        :param map_path:
        :param mode: 1 is {metab:Metabolite},2 is {Metabolite:metab}
        '''
        metab_trans = Relation()
        id_name_dict = {}
        if metab_name:
            metab_name = metab_name
        else:
            metab_name = "Metabolite"
        if id_name:
            id_name =id_name
        else:
            id_name = "metab_id"
        metab_trans.get_table_info(map_path,id_name=id_name)
        if mode == 2:
            id_name_dict = metab_trans.return_dic(metab_name)
            map_table = pd.DataFrame.from_dict(id_name_dict, orient='index')
            map_table[metab_name] = map_table.index
            map_table.columns = [id_name, metab_name]
        elif mode == 1:
            id_name_dict = metab_trans.return_id_dic(metab_name)
            map_table = pd.DataFrame.from_dict(id_name_dict, orient='index')
            map_table[id_name] = map_table.index
            map_table.columns = [metab_name, id_name]
        return map_table, id_name_dict

    def get_metab_to_idfile(self, map_table, id_name_dict, oldfile, newfile, myrow=True, mycol=False):
        """
        根据map_file中metab_id,替换为metab名称,生成新file,map_file和old_file第一列内容为metab_id
        """
        old_table = pd.read_table(oldfile, sep="\t", header=0)
        old_table.rename(columns={old_table.columns[0]: "Metabolite"}, inplace=True)
        merge = old_table
        #print merge.head()
        if myrow:
            tmp_merge = pd.merge(merge, map_table, how="left", on="Metabolite")
            columns_order = ["metab_id"] + merge.columns.tolist()
            merge = tmp_merge[columns_order]
        if mycol:
            new_columns = []
            id_name_dict["metab_id"] = "metab_id"
            id_name_dict["Metabolite"] = "Metabolite"
            [new_columns.append(id_name_dict[x]) for x in merge.columns]
            merge.columns = [new_columns]

        def _tmp_fun(x):
            if pd.isna(x['metab_id']):
                x['metab_id'] = x['Metabolite']
            return x
        #merge = merge.apply(_tmp_fun,axis=1)
        merge = merge.drop('Metabolite', 1)
        merge.to_csv(newfile, sep="\t", index=False, quoting=3)

    def trans_col_data(self, id_name_dict, infile, outfile):
        with open(infile, "r") as f, open(outfile, "w") as outf:
            for line in f:
                line = line.strip().split("\t")
                name = line[0]
                ids = line[1]
                id_list = ids.split(";")
                data_list = []
                for each in id_list:
                    if id_name_dict.has_key(each):
                        metab = id_name_dict[each]
                        data_list.append(metab)
                    else:
                        data_list.append(each)
                metabs = ";".join(data_list)
                outf.write(name + "\t" + metabs + "\n")

    def add_metabolites_column(self, id_name_dict, oldfile, newfile, id_name=None):  # add by ysh 20190409
        with open(oldfile, "r") as f, open(newfile, "w") as outf:
            head = f.readline().strip().split("\t")
            if id_name:
                id_name =id_name
            else:
                id_name = "metab_id"
            if head:
                metab_index = head.index(id_name)
                outf.write("\t".join(head) + "\t" + "Metabolite" + "\n")
            for line in f:
                line = line.strip()
                line1 = line.split("\t")
                metab_ids = line1[metab_index]
                id_list = re.split('[;|]',metab_ids)
                metabolite_list =[]
                for x in id_list:
                    if x in id_name_dict.keys():
                        metabolite_list.append(id_name_dict[x])
                    else:
                        metabolite_list.append(x)
                metabolites = ";".join(metabolite_list)
                outf.write(line + "\t" + metabolites + "\n")

    ## 从desc_file 中获取o_id 和 metab_id 或 Metabolite 的对应关系。oldfile中增加一列o_id
    def add_oid_fun(self,oldfile,id_desc_file,oid='ID',link_k='Metabolite',oldfile_link_id=0):
        mdata = pd.read_table(id_desc_file,sep='\t',header=0)
        if oid not in mdata.columns:
            return
        oids = mdata[oid].tolist()
        link_ks = mdata[link_k].tolist()
        rets = dict(zip(link_ks,oids))
        with open(oldfile+'_tmp_add_oid','w') as fw, open(oldfile) as fr:
            head = fr.readline()
            if 'ID' in head:
                return
            fw.write('ID\t'+head)
            for line in fr:
                name = line.split('\t')[oldfile_link_id]
                if name in rets.keys():
                    fw.write(rets[name]+'\t'+line)
                else:
                    print(name + ' not in :')
                    print(line)
                    fw.write('-\t'+line)
                    #fw.write(' \t'+line)
        os.remove(oldfile)
        os.rename(oldfile+'_tmp_add_oid', oldfile)












if __name__ == "__main__":  # test Ralation Class
    import time

    table_path = "metab_desc.txt"
    collection = "exp_mix"
    main_name = "exp_id"
    main_id = "5b026d8da4e1af4146c58520"
    obj = Relation()
    obj.get_table_info(table_path)  # 以本地文件的方式
    obj.save_map("KEGG Compound ID", "compound_table.save")  # 输出map文件测试
    mongo_obj = Relation()
    mongo_obj.get_mongo_info(collection, main_name, main_id)  # 以mongo库的方式
    mongo_obj.save_map("kegg_id", "compound_mongo.save")  # 输出map文件测试
    # 输出diction测试
    print "read table in ok"
    time.sleep(5)
    f1 = open("return_dic_table.log", "w")
    str_ = obj.return_dic("KEGG Compound ID")
    str_ = str(str_)
    f1.write(str_)
    f1.close()
    f2 = open("return_dic_mongo.log", "w")
    str_ = mongo_obj.return_dic("kegg_id")
    str_ = str(str_)
    f2.write(str_)
    f2.close()

    # 输出diction测试
    f3 = open("return_id_table.log", "w")
    str_ = obj.return_id_dic("KEGG Compound ID")
    str_ = str(str_)
    f3.write(str_)
    f3.close()
    f4 = open("return_id_mongo.log", "w")
    str_ = mongo_obj.return_id_dic("kegg_id")
    str_ = str(str_)
    f4.write(str_)
    f4.close()

    # 输出转换后的id测试
    try:
        print "convert to id success: %s" % obj.convert_to_id("KEGG Compound ID", "C12859")
    except:
        raise Exception(e)
    try:
        print "convert to id success %s" % mongo_obj.convert_to_id("kegg_id", "C12859")
    except:
        raise Eception(e)

    # 输出id对应的测试
    try:
        print "convert id success: %s" % obj.convert_id("KEGG Compound ID", "metab_44")
    except:
        raise Exception(e)
    try:
        print "convert id success: %s" % mongo_obj.convert_id("kegg_id", "metab_44")
    except:
        raise Exception(e)
