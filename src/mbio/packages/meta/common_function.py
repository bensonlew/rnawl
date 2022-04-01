# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
## @2020922


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
from collections import defaultdict

def get_mongo():
    client = Config().get_mongo_client(mtype="meta")
    db = client[Config().get_mongo_dbname("meta")]
    return client, db

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


def get_sg_sample_check(task_id, query_id):
    """
    从MongoDBdb中根据task_id从sample_check中查到样本信息
    :param task_id:
    :param query_id: 文件对应的位置
    :return:
    """
    mongo,db = get_mongo()
    #db = Config().get_mongo_client(mtype="meta", ref=False)[Config().get_mongo_dbname("meta", ref=False)]
    #client = Config().get_mongo_client(mtype="meta", ref=True) # by zzg 20210223 关闭mongo连接
    #db = client[Config().get_mongo_dbname("meta", ref=True)]
    file_sample_dict = {}
    sample_check = db["sg_sample_check"]
    sample_check_detail = db["sg_sample_check_detail"]
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
                if file_name:
                    true_file_name = file_name
                else:
                    true_file_name = result['true_file_name']
                if is_check in ["true"]:
                    origin_specimen_dict[origin_specimen] = specimen
                if true_file_name in file_sample_dict.keys():
                    file_name_list = file_sample_dict[true_file_name]
                else:
                    file_name_list = []
                if origin_specimen_dict not in file_name_list:
                    if origin_specimen_dict != {}:
                        file_name_list.append(origin_specimen_dict)
                    file_sample_dict[true_file_name] = file_name_list
    #client.close()
    return file_sample_dict


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
    db = Config().get_mongo_client(mtype="meta", ref=False)[Config().get_mongo_dbname("meta", ref=False)]
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


def envname_info(infile):
    direction, col, table, lable = 'row', 0, True, False
    for postfix in ['Maaslin.xls', 'correlation.xls', 'pvalue.xls',
                    'pearsons_pvalue.xls', 'pearsons_correlation.xls']:
        if infile.endswith(postfix):
            return 'row', 0, True, True
    for postfix in ['rda_vif.txt']:
        if infile.endswith(postfix):
            return 'row', None, True, True
    for postfix in ['_biplot.xls', '_envfit.xls', 'env.R2adj.xls',
                    'Maaslin.txt', 'randomForest_imptance_table.xls',
                    'envfit_vector.xls', 'envfit_vector_scores.xls',
                    'randomForest_top3_vimp.xls', 'permanova.xls']:
        if infile.endswith(postfix):
            return 'col', 0, True, True
    for postfix in ['corr_network_clustering.txt', 'corr_network_by_cut.txt',
                    'corr_network_centrality.txt', 'corr_network_node_degree.txt']:
        if infile.endswith(postfix):
            return 'col', 1, False, True
    return 'row', 0, True, False


def envname_hander(env_id, infile, newfile, direction='row', col=0, table=False, db_type='meta'):
    """在上传的结果文件中恢复真实的环境因子名称

    Args:
        env_id (str): 环境因子主表id
        infile (str): 输入文件，包含待转换的环境因子名称
        newfile (str): 改名后的输出文件路径
        direction (str, optional): 环境因子是行(row)或列(col)上. Defaults to 'row'.
        col (int, optional): 指定环境因子在第几列或行. Defaults to 0.
        other (bool, optional): 表明输入文件是否是有行列名成的二维表 . Defaults to False.
    """
    db = Config().get_mongo_client(mtype=db_type)[Config().get_mongo_dbname(db_type)]
    if db_type == 'metagenomic':
        col_name = 'env'
    else:
        col_name = 'sg_env'
    env_info = db[col_name].find_one({'_id': ObjectId(env_id)})
    print(db)
    print(env_id)
    print(db_type)
    print(col_name)
    ana_envnames = env_info['env_names'].split(',')
    if 'env_true_name' in env_info:
        true_envnames = env_info['env_true_name']
    else:
        true_envnames = ana_envnames
    names_dict = defaultdict(str)
    names_dict.update(dict(zip(ana_envnames, true_envnames)))
    if table:
        tb = pd.read_csv(infile, sep='\t')
        if direction == 'row':
            tb.columns = map(lambda x: names_dict[x] or x, tb.columns)
        elif direction == 'col':
            col = tb.columns[col]
            tb[col] = map(lambda x: names_dict[x] or x, tb[col])
        columns = ["" if "Unnamed: " in c else c for c in tb.columns]
        tb.columns = columns
        tb.to_csv(newfile, sep='\t', index=False)
    else:
        with open(infile, 'r') as fin, open(newfile, 'w') as fout:
            for line in fin:
                li = line.strip().split('\t')
                if len(li) > col:
                    if li[col] in names_dict:
                        li[col] = names_dict[li[col]] or li[col]
                    if infile.endswith('corr_network_by_cut.txt'):
                        li[0] = names_dict[li[0]] or li[0]
                fout.write('\t'.join(li) + '\n')


def env_link(obj, env_id, olddir, newdir, db_type='meta'):
    if not os.path.isdir(olddir):
        raise Exception("不存在路径: %s" % olddir)
    allfiles = os.listdir(olddir)
    if not os.path.exists(newdir):
        os.makedirs(newdir)
    newfiles = [os.path.join(newdir, i) for i in allfiles]
    for newfile in newfiles:
        if os.path.exists(newfile):
            if os.path.isfile(newfile):
                os.remove(newfile)
            elif os.path.isdir(newfile):
                shutil.rmtree(newfile)
    if len(allfiles) >= 1:
        for i in allfiles:
            infile = os.path.join(olddir, i)
            newfile = os.path.join(newdir, i)
            if os.path.isfile(os.path.join(olddir, i)):
                direction, col, table, lable = envname_info(infile)
                if not lable:
                    os.link(infile, newfile)
                else:
                    envname_hander(env_id, infile, newfile,
                                   direction, col, table, db_type)
            elif os.path.isdir(os.path.join(olddir, i)):
                env_link(obj, env_id, infile, newfile, db_type)
    else:
        raise Exception("结果文件夹为空:%s" % allfiles)


def envname_restore(func):
    @functools.wraps(func)
    def wrapper(self, *arg, **kwargs):
        env_file = 'env_file'
        if 'envtable' in self._options:
            env_file = 'envtable'
        if self.option(env_file).is_set:
            env_id = self.option('env_id')
            db_type = self._sheet.db_type
            if env_id:
                self.logger.info('修改结果文件中环境因子名称为真实名称')
                outpath = os.path.join(self.work_dir, 'output_upload')
                env_link(self, env_id, self.output_dir, outpath, db_type)
                self.logger.info('完成名称修改')
                self.output_dir = outpath
        func(self, *arg, **kwargs)
    return wrapper

def filter_otu_set(otu_table, my_json,otu_id, level, output):
    """
    过滤和保留基因集
    :param otu_table: 输入的otu_table表
    :param my_json: 过滤的json信息
    :param otu_id: otu_id的表
    :param level: 选择的分类水平
    :param output: 输出文件信息
    :return:
    """
    mongo, db = get_mongo()
    #db = Config().get_mongo_client(mtype="meta", ref=False)[Config().get_mongo_dbname("meta", ref=False)]
    level_dict = {1:"d__", 2:"k__",3:"p__",4:"c__",5:"o__",6:"f__",7:"g__",8:"s__",9:"otu",}
    json_info = my_json
    set_id = ObjectId(json_info)
    collection = db['sg_otuset']
    new_otu_id = ObjectId(otu_id)
    result = collection.find_one({"_id": set_id})
    if result:
        otu_list = set(result["otuset_list"])
    else:
        raise Exception("未能成功的找到对应的基因集！")
    otu_collection = db['sg_otu_detail_level']
    if level in level_dict:
        find_level = level_dict[level]
    else:
        find_level = "otu"
    results = otu_collection.find({"otu_id": new_otu_id, "level_id": 9})
    choose_list = []
    if result:
        for res in results:
            otu_level = res['otu']
            if otu_level in otu_list:
                choose_level = str(res[find_level])
                if choose_level not in choose_list:
                    choose_list.append(choose_level)
    else:
        raise Exception("未能成功查找到对应的level水平，请检查otu_id是否正确或者数据注释信息是否正确！")

    data = pd.read_table(otu_table, sep='\t', header=0)
    otu_table_list = data['OTU ID'].values
    origin_otu_list = []
    origin_otu_dict = {}
    for otu in otu_table_list:
        otu_name = re.split(r';', otu)[-1].strip()
        if otu_name not in origin_otu_list:
            origin_otu_list.append(otu_name)
        if otu_name not in origin_otu_dict:
            origin_otu_dict[otu_name] = otu
    choose_list = set(choose_list)
    print(choose_list)
    print(origin_otu_list)
    origin_otu_list = set(origin_otu_list)
    origin_otu_list2 = origin_otu_list.intersection(choose_list)
    print("++++++++++++++++",origin_otu_list2, "----------------")
    if len(origin_otu_list) == 0:
        raise Exception("经过筛选otu集的结果为空，请切换otu集！")
    otu_full_name_list = []
    for old_otu in origin_otu_list2:
        if old_otu in origin_otu_dict:
            full_otu =  origin_otu_dict[old_otu]
            if full_otu not in otu_full_name_list:
                otu_full_name_list.append(full_otu)
    print(otu_full_name_list)
    new_data = data[data['OTU ID'].isin(otu_full_name_list)]

    new_data.to_csv(output, sep="\t", index=0)

def get_save_pdf_status(task_id=None):
    """
    查询任务save_pdf字段
    :param task_id: task_id
    :return:
    """
    """
    mongo, db = get_mongo()
    #db = Config().get_mongo_client(mtype="meta", ref=False)[Config().get_mongo_dbname("meta", ref=False)]

    collection = db['sg_task']
    result = collection.find_one({"task_id": task_id})
    if result:
        if "save_pdf" in result and result["save_pdf"] == 1:
            save_pdf = True
        else:
            save_pdf = False
    else:
        save_pdf = False
    """
    return False

def get_name(table_id=None,table=None,name="name"):
    """
        查询主表name字段
        :param table_id: table_id
        :return:
        """
    mongo, db = get_mongo()
    # db = Config().get_mongo_client(mtype="meta", ref=False)[Config().get_mongo_dbname("meta", ref=False)]

    collection = db[table]
    result = collection.find_one({"_id": ObjectId(table_id)})
    return result[name]

def get_pan_id(table_id=None,table=None,name="name"):
    """
        查询主表name字段
        :param table_id: table_id
        :return:
        """
    mongo, db = get_mongo()
    # db = Config().get_mongo_client(mtype="meta", ref=False)[Config().get_mongo_dbname("meta", ref=False)]

    collection = db[table]
    result = collection.find_one({"_id": ObjectId(table_id)})
    return result[name]

def group_file_spilt(group_file, output_dir):
    group_detail_list = []
    group_n = []
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)
    df = pd.read_table(group_file,header=0,index_col=None)
    for x in df.columns.values[1:4]:
        df_single = df[[df.columns.values[0],x]]
        df_dropna = df_single.dropna(subset=[x])
        df_dropna.to_csv(output_dir+"/"+x+".group.txt",sep="\t",index=None)

def get_level_id(table_name,table_id,keys):
    """
    查询任务save_pdf字段
    :param task_id: task_id
    :return:
    """
    mongo, db = get_mongo()
    #db = Config().get_mongo_client(mtype="meta", ref=False)[Config().get_mongo_dbname("meta", ref=False)]

    collection = db[table_name]
    result = collection.find_one({"_id": ObjectId(table_id)})
    level_id = result[keys]
    return level_id