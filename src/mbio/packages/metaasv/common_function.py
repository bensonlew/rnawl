# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
## @20200404


import pandas as pd
import os
import re
import functools
import time
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import json
import pandas as pd
import math
from bson import ObjectId


def link_dir(olddir, newdir):
    """
    hard link directory from olddir to newdir
    :param olddir: old dir
    :param newdir: new dir
    :return:
    """
    if not os.path.isdir(olddir):
        raise Exception("不存在路径: %s" % olddir)
    allfiles = os.listdir(olddir)
    if not os.path.exists(newdir):
        os.makedirs(newdir)
    newfiles = [os.path.join(newdir,i) for i in allfiles]
    for newfile in newfiles:
        if os.path.exists(newfile):
            if os.path.isfile(newfile):
                os.remove(newfile)
            elif os.path.isdir(newfile):
                shutil.rmtree(newfile)
    if len(allfiles) >= 1:
        for i in allfiles:
            if os.path.isfile(os.path.join(olddir, i)):
                os.link(os.path.join(olddir, i),os.path.join(newdir, i))
            elif os.path.isdir(os.path.join(olddir, i)):
                link_dir(os.path.join(olddir, i),os.path.join(newdir, i))
    else:
        raise Exception("结果文件夹为空:%s" % allfiles)


def link_file(oldfile, newfile):
    """
    hard link file from oldfile to newfile
    :param oldfile: oldfile
    :param newfile: newfile
    :return:
    """
    if not os.path.isfile(oldfile):
        raise Exception("不存在文件：%s" % oldfile)
    if os.path.exists(newfile):
        os.remove(newfile)
    os.link(oldfile, newfile)


def time_count(func):
    # 统计函数运行时间，作为方法的装饰器
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
            wait_file(path, times=times, sleep=sleep)
        return path
    raise Exception("超过文件等待次数，需检查文件%s" % path)


def check_option(options, property, list=None, min=None, max=None):
    """
    参数二次检查
    :param options: 需要检查的options名称，有四种属性infile、string、int、float
    :param property: 属性
    :param list: 是否传过来list属性
    :param min: 最小值
    :param max: 最大值
    :return: True
    """
    db = Config().get_mongo_client(mtype="metaasv", ref=False)[Config().get_mongo_dbname("metaasv", ref=False)]
    if property in ['infile']:
        if not options.is_set:
            raise OptionError("{}不存在文件：{}！".format(property, options.prop['path']))

    elif property in ['string']:
        if list:
            if options not in list:
                raise OptionError("{}不在：{}中！".format(options, list))
        else:
            if not options:
                raise OptionError("{}不存在！".format(options))

    elif property in ['float']:
        if list:
            if options not in list:
                raise OptionError("{}不在：{}中！".format(options, list))
        else:
            if not options:
                raise OptionError("{}不存在！".format(options))
        if min and max:
            if options > max and options < min:
                raise OptionError("{} 必须大于{}，小于{}".format(options.name, min, max))

    elif property in ['int']:
        if list:
            if options not in list:
                raise OptionError("{}不在：{}中！".format(options, list))
        else:
            if not options:
                raise OptionError("{}不存在！".format(options))
        if min and max:
            if options > max and options < min:
                raise OptionError("{} 必须大于{}，小于{}".format(options.name, min, max))


def calculate_abundance(asv_table, group=None, method=None, top=None, abundance_method=None):
    """
    功能：计算丰度
    :param asv_table: asv丰度表
    :param group: group表
    :param method: 是否按照分组合并，合并的方法
    :param top: 是否挑选top物种
    :param abundance_method: 计算绝对丰度还是相对丰度
    :return:
    """
    db = Config().get_mongo_client(mtype="metaasv", ref=False)[Config().get_mongo_dbname("metaasv", ref=False)]
    table = pd.read_table(asv_table, sep='\t', index_col=0, header=0)
    table.columns = [str(i) for i in table.columns]
    table_name = table.index.name
    if group:
        group_table = pd.read_table(group, sep='\t', index_col=0, header=0)
        group_table["sample"] = [str(i) for i in group_table.index]
        group_sample = ""
        if method in ["sum"]:
            group_sample = group_table.join(table.T, on="sample").groupby(group_table.columns[0]).sum()  # 求和
        elif method in ["average"]:
            group_sample = group_table.join(table.T, on="sample").groupby(group_table.columns[0]).mean()  # 求均值
        elif method in ["middle"]:
            group_sample = group_table.join(table.T, on="sample").groupby(group_table.columns[0]).median()  # 中位数
        elif method not in ["average", "sum", "middle"]:
            group_sample = group_table.join(table.T, on="sample")
            group_sample.drop(group_sample.columns[:2], axis=1, inplace=True)
        abund = group_sample.T
        abund.index.name = table_name
    else:
        abund = table
    abund['all_sum'] = abund.apply(lambda x: x.sum(), axis=1)
    abund_table = abund.sort_values(by=['all_sum'], ascending=0)
    del abund_table["all_sum"]
    abund_table = abund_table.ix[list((abund_table > 0).any(axis=1))] #去除都为0的物种
    if len(abund_table) < 1:
        raise Exception('在所选参数下数据为空，请重新设置水平或分组方案参数!')
    dir_name = os.path.dirname(asv_table)
    if top not in ["", None]:
        top_num = int(top)
        abund_table = abund_table.head(top_num)
    if abundance_method in ['absolute']:
        abund_table_path = os.path.join(dir_name, "abundance.xls")
        abund_table.to_csv(abund_table_path, sep="\t", encoding="utf-8")
        return abund_table
    elif abundance_method in ['relative']:
        abund_table.columns = [str(i) for i in abund_table.columns]
        abund_table.loc['row_sum'] = abund_table.apply(lambda x: x.sum(), axis=0)
        abund_table_percent = abund_table.apply(lambda x: x / abund_table.loc['row_sum'], axis=1).drop('row_sum')
        abund_table_path = os.path.join(dir_name, "abundance.percents.xls")
        abund_table_percent.to_csv(abund_table_path, sep="\t", encoding="utf-8")
        return abund_table_path

def calculate_asv(asv_table, group=None, method=None, top=None, abundance_method=None):
    """
    功能：计算asv的相对丰度
    :param asv_table: asv丰度表
    :param group: group表
    :param method: 是否按照分组合并，合并的方法
    :param top: 是否挑选top物种
    :param abundance_method: 计算绝对丰度还是相对丰度
    :return:
    """
    table = pd.read_table(asv_table, sep='\t', index_col=0, header=0)
    table.columns = [str(i) for i in table.columns]
    abund = table
    abund['Total'] = abund.apply(lambda x: x.sum(), axis=1)
    abund_table = abund.sort_values(by=['Total'], ascending=0)
    # asv_number = 0
    # for num in list(abund['Total'].values):
    #     asv_number += num
    # del abund_table["Total"]
    # abund['Percent'] = abund_table.apply(lambda x: float(x) / asv_number, axis=1)
    abund_table = abund_table.ix[list((abund_table > 0).any(axis=1))] #去除都为0的物种
    if len(abund_table) < 1:
        raise Exception('在所选参数下数据为空，请重新设置水平或分组方案参数!')
    dir_name = os.path.dirname(asv_table)
    if top not in ["", None]:
        top_num = int(top)
        abund_table = abund_table.head(top_num)
    if abundance_method in ['absolute']:
        abund_table_path = os.path.join(dir_name, "abundance.xls")
        abund_table.to_csv(abund_table_path, sep="\t", encoding="utf-8")
        return abund_table
    elif abundance_method in ['relative']:
        abund_table.columns = [str(i) for i in abund_table.columns[0:-1]] + ["Percent"]
        abund_table.loc['row_sum'] = abund_table.apply(lambda x: x.sum(), axis=0)
        abund_table_percent = abund_table.apply(lambda x: x / abund_table.loc['row_sum'], axis=1).drop('row_sum')
        abund_table_path = os.path.join(dir_name, "abundance.percents.xls")
        abund_table_percent.to_csv(abund_table_path, sep="\t", encoding="utf-8")
        print(abund_table.columns)
        return abund_table_path


def filter_zero_and_replace(asv_table,output):
    """
    将丰度表中的所有丰度为0的在此基础上加上整张丰度表最小值的十分之一
    :param asv_table:
    :return:
    """
    min_list = []
    data_table = pd.read_table(asv_table, sep="\t", header=0,index_col=0)
    # data_table = asv_table
    columns = [str(i) for i in list(data_table.columns)[1:]]
    for column in columns:
        data_list = data_table[data_table[column] != 0.0][column]
        if len(data_list) != 0:
            min_num = min(data_list)
            min_list.append(min_num)
    all_min_num = min(min_list)
    data_table = data_table.replace(0.0, all_min_num)
    data_table.to_csv(output, sep="\t")
    return data_table

def filter_zero(asv_table):
    """
    将丰度表中的所有丰度为0的去掉
    :param asv_table:
    :return:
    """
    min_list = []
    dir_name = os.path.dirname(asv_table)
    new_path = os.path.join(dir_name, "input_asv.xls")
    data_table = pd.read_table(asv_table, sep="\t", header=0,index_col=0)
    columns = [str(i) for i in data_table.columns]
    for column in columns:
        data_list = [data_table[column] != 0.0]
        min_num = min(data_list)
        min_list.append(min_num)
    all_min_num = min(min_list)
    data_table = data_table.replace(0.0, all_min_num)
    data_table.to_csv(new_path, sep="\t")
    return new_path


def filter_asv_set(asv_table, my_json,asv_id, level, output):
    """
    过滤和保留基因集
    :param asv_table: 输入的ASV_table表
    :param my_json: 过滤的json信息
    :param asv_id: asv_id的表
    :param level: 选择的分类水平
    :param output: 输出文件信息
    :return:
    """
    db = Config().get_mongo_client(mtype="metaasv", ref=False)[Config().get_mongo_dbname("metaasv", ref=False)]
    level_dict = {1:"d__", 2:"k__",3:"p__",4:"c__",5:"o__",6:"f__",7:"g__",8:"s__",9:"asv",}
    json_info = my_json
    set_id = ObjectId(json_info)
    collection = db['asv_set']
    new_asv_id = ObjectId(asv_id)
    result = collection.find_one({"_id": set_id})
    if result:
        asv_list = set(result["asv_list"])
    else:
        raise Exception("未能成功的找到对应的基因集！")
    asv_collection = db['asv_detail_level']
    if level in level_dict:
        find_level = level_dict[level]
    else:
        find_level = "asv"
    results = asv_collection.find({"asv_id": new_asv_id, "level_id": 9})
    choose_list = []
    if result:
        for res in results:
            asv_level = res['asv']
            if asv_level in asv_list:
                choose_level = str(res[find_level])
                if choose_level not in choose_list:
                    choose_list.append(choose_level)
    else:
        raise Exception("未能成功查找到对应的level水平，请检查asv_id是否正确或者数据注释信息是否正确！")

    data = pd.read_table(asv_table, sep='\t', header=0)
    asv_table_list = data['ASV ID'].values
    origin_asv_list = []
    origin_asv_dict = {}
    for asv in asv_table_list:
        asv_name = re.split(r';', asv)[-1].strip()
        if asv_name not in origin_asv_list:
            origin_asv_list.append(asv_name)
        if asv_name not in origin_asv_dict:
            origin_asv_dict[asv_name] = asv
    choose_list = set(choose_list)
    print(choose_list)
    print(origin_asv_list)
    origin_asv_list = set(origin_asv_list)
    origin_asv_list2 = origin_asv_list.intersection(choose_list)
    print("++++++++++++++++",origin_asv_list2, "----------------")
    if len(origin_asv_list) == 0:
        raise Exception("经过筛选ASV集的结果为空，请切换ASV集！")
    asv_full_name_list = []
    for old_asv in origin_asv_list2:
        if old_asv in origin_asv_dict:
            full_asv =  origin_asv_dict[old_asv]
            if full_asv not in asv_full_name_list:
                asv_full_name_list.append(full_asv)
    print(asv_full_name_list)
    new_data = data[data['ASV ID'].isin(asv_full_name_list)]

    new_data.to_csv(output, sep="\t", index=0)

def filter_asv(set_id, mongo_version=None):
    """
    返回asv集
    :param my_json:
    :return:
    """
    if mongo_version:
        Config.DBVersion = mongo_version
    db = Config().get_mongo_client(mtype="metaasv", ref=False)[Config().get_mongo_dbname("metaasv", ref=False)]
    set_id = ObjectId(set_id)
    collection = db['asv_set']
    result = collection.find_one({"_id": set_id})
    if result:
        asv_list = set(result["asv_list"])
        return(asv_list)
    else:
        raise OptionError("ASV集查不到，请传入正确的set_id")


def normalize_data(asv_table,outfile, abundance_method=None, type=None):
    """
    功能：计算丰度、排序
    :param asv_table: asv丰度表
    :param abundance_method: 计算方法
    :param type: 类型 环境因子或者asv表(asv,env)
    :return:
    """
    table = pd.read_table(asv_table, sep='\t', index_col=0, header=0)
    table.columns = [str(i) for i in table.columns]
    abund = table
    abund['all_sum'] = abund.apply(lambda x: x.sum(), axis=1)
    abund_table = abund.sort_values(by=['all_sum'], ascending=0)
    del abund_table["all_sum"]
    # abund_table = abund_table.ix[list((abund_table > 0).any(axis=1))] #去除都为0的物种
    ##暂时去掉上面那一行，因为存在预测样本的时候随机森林报错
    if type in ["asv"]:
        abund_table.reset_index(inplace=True)
        abund_table.drop_duplicates(subset="ASV ID",keep="first", inplace=True) #去除重复的ASV_ID的物种
        abund_table.set_index("ASV ID", inplace=True)
    if len(abund_table) < 1:
        raise Exception('在所选参数下数据为空，请重新设置水平或分组方案参数!')
    dir_name = os.path.dirname(asv_table)
    if abundance_method in ['absolute']:
        abund_table.to_csv(outfile, sep="\t", encoding="utf-8")
        return outfile
    elif abundance_method in ['Relative']:
        abund_table.columns = [str(i) for i in abund_table.columns]
        abund_table.loc['row_sum'] = abund_table.apply(lambda x: x.sum(), axis=0)
        abund_table_percent = abund_table.apply(lambda x: x / abund_table.loc['row_sum'], axis=1).drop('row_sum')
        abund_table_percent.to_csv(outfile, sep="\t", encoding="utf-8")
        return outfile
    elif abundance_method in ['Min-Max']:
        abund_table.columns = [str(i) for i in abund_table.columns]
        for column in list(abund_table.columns):
            min_number = min(list(abund_table[column].values))
            max_number = max(list(abund_table[column].values))
            abund_table[column] = abund_table[column].apply(lambda x: float(x-min_number) /(max_number-min_number))
        abund_table.to_csv(outfile, sep="\t", encoding="utf-8")
        return outfile
    elif abundance_method in ['log10']:
        abund_table.columns = [str(i) for i in abund_table.columns]
        for column in list(abund_table.columns):
            max_number = max(list(abund_table[column]))
            abund_table[column] = abund_table[column].apply(lambda x: x if x == 0.0 else math.log(x) / math.log(max_number))
        abund_table.to_csv(outfile, sep="\t", encoding="utf-8")
        return outfile
    elif abundance_method in ['Z-score']:
        abund_table.columns = [str(i) for i in abund_table.columns]
        for column in abund_table.columns:
            abund_table[column] = (abund_table[column] - abund_table[column].mean()) / abund_table[column].std()
            # abund_table[column] = z_score.abs() > 2.2
        # abund_table_path = os.path.join(dir_name, "abundance.xls")
        abund_table.to_csv(outfile, sep="\t", encoding="utf-8")
        return outfile


def copy_task(task_id,from_task,insert_id, collection,main_collection, key_relation=None,query=None):
    """
    根据指定的任务copy详情表数据
    :param task_id: 任务id
    :param from_task: 来源任务id
    :param insert_id: 要插入的主表id
    :param collection: 详情表
    :param main_collection: 主表
    :param key_relation: 主表和详情表关联字段
    :return:
    """
    db = Config().get_mongo_client(mtype="metaasv", ref=False)[Config().get_mongo_dbname("metaasv", ref=False)]
    collection_detail = db[collection]
    main_collection_choose = db[main_collection]
    if query:
        main_result = main_collection_choose.find_one({"task_id": from_task, "query_id": query, "status" : "end"})
    else:
        main_result = main_collection_choose.find_one({"task_id": from_task, "status" : "end"})
    if main_result:
        data_list = []
        main_id = main_result["_id"]
        results = collection_detail.find({key_relation: main_id})
        if results:
            for res in results:
                del res['_id']
                res[key_relation] = ObjectId(insert_id)
                if res not in data_list:
                    data_list.append(res)
        else:
            raise OptionError("在表{}中查不到对应的task_id:{}的任务".format(collection, from_task))
        try:
            collection_detail.insert_many(data_list)
        except Exception,e:
            raise Exception("copy 详情表:%s数据失败%s"%(collection,e))
        try:
            resul = main_collection_choose.find_one({"task_id": task_id, "query_id": query})
            if resul:
                new_main_id = resul["_id"]
                main_collection_choose.update_one({"task_id": task_id, "query_id": query}, {"$set": {"main_id": new_main_id}})
        except Exception,e:
            raise Exception("更新主表%s失败！%s"%(collection,e))
    else:
        raise OptionError("在表{}中查不到对应的task_id:{}".format(main_collection, from_task))


def check_file(file_type, file_path):
    """
    检查文件类型是否是压缩格式
    :param file_type: 文件或者文件夹（file dir）
    :param file_path: 输入文件的路径
    :return:
    """
    if file_type in ["file"]:
        file_name = os.path.basename(file_path)
        if re.search(r"\.gz$", file_name):
            return True
        elif re.search(r"tar\.gz$", file_name):
            return True
        elif re.search(r"\.zip$", file_name):
            return True
        else:
            return False
    elif file_type in ["dir"]:
        list_path = os.path.join(file_path, "list.txt")
        if not os.path.exists(list_path):
            dir_list = os.listdir(file_path)
            with open(os.path.join(file_path, "list.txt"), "w") as w:
                for file in dir_list:
                    if re.search(r".gz", file):
                        sample_name = file.strip(".gz")
                    elif re.search(r"tar.gz", file):
                        sample_name = file.strip(".tar.gz")
                    elif re.search(r".zip", file):
                        sample_name = file.strip(".zip")
                    else:
                        sample_name = file
                    w.write("{}\t{}\t{}\n".format(file, sample_name, "s"))
            for file in dir_list:
                if re.search(r".gz|tar.gz|.zip", file):
                    return True
        else:
            i = 1
            with open(list_path, 'r') as f:##读第一行的文件名称，有标题
                lines = f.readlines()
                for line in lines[1:]:
                    i += 1
                    if i == 3:
                        break
                    line = line.strip().split("\t")
                    file_name = line[0]
                    if re.search(r"\.gz$", file_name):
                        break
                    elif re.search(r"tar\.gz$", file_name):
                        return True
                    elif re.search(r"\.zip$", file_name):
                        return True
                    else:
                        return False


def get_sample_check(task_id, query_id):
    """
    从MongoDBdb中根据task_id从sample_check中查到样本信息
    :param task_id:
    :param query_id: 文件对应的位置
    :return:
    """
    db = Config().get_mongo_client(mtype="metaasv", ref=False)[Config().get_mongo_dbname("metaasv", ref=False)]
    file_sample_dict = {}
    origin_new_dict = {}
    file_name_list = []
    sample_check = db["sample_check"]
    sample_check_detail = db["sample_check_detail"]
    main_results = sample_check.find_one({"task_id": task_id, "query_id" : query_id, "status": "end"})
    if main_results:
        main_id = main_results["main_id"]
        results = sample_check_detail.find({"check_id": main_id})
        if results:
            for result in results:
                origin_specimen_dict = {}
                is_check = result["is_check"]
                origin_specimen = result["origin_specimen"]
                specimen = result["specimen"]
                file_name = result['file_name']
                if "file_name" in result:
                    true_file_name = file_name
                else:
                    true_file_name = result['true_file_name']
                if is_check in ["true"]:
                    origin_specimen_dict[origin_specimen] = specimen
                    # origin_new_dict[file_name][origin_specimen] = specimen
                if true_file_name in file_sample_dict.keys():
                    file_name_list = file_sample_dict[true_file_name]
                else:
                    file_name_list = []
                if origin_specimen_dict not in file_name_list:
                    if origin_specimen_dict != {}:
                        file_name_list.append(origin_specimen_dict)
                    file_sample_dict[true_file_name] = file_name_list
    return file_sample_dict


def get_group_from_table(asv_table,outfile):
    """
    根据asv表获得group表，设置默认为all的情况
    :param asv_table:
    :param outfile: 输出结果文件
    :return:
    """
    with open(asv_table, 'r') as f, open(outfile, 'w') as w:
        w.write("#sample\tAll\n")
        line = f.readline()
        line = line.strip().split("\t")
        for name in line[1:]:
            w.write("{}\t{}\n".format(name.strip(), "All"))
    return outfile

def find_group_name(task_id):
    """
    从MongoDB中查找分组名称，并返回group_detail
    :param task_id:
    :return:
    """
    db = Config().get_mongo_client(mtype="metaasv", ref=False)[Config().get_mongo_dbname("metaasv", ref=False)]
    collection = db["specimen_group"]
    tempfind = collection.find_one({'task_id': task_id})
    group_detail = {}
    if tempfind:
        group_names = tempfind['category_names']
        specimen_names = tempfind["specimen_names"]
        n_group_names = range(len(group_names))
        for i in n_group_names:
            group = group_names[i]
            specimen = specimen_names[i]
            if group not in group_detail:
                group_detail[group] = [str(x) for x in specimen.keys()]
        return group_detail