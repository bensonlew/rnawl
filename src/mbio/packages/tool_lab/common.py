# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify: 20210706

import sys,os
from biocluster.config import Config
import re
from biocluster.file import download,exists
import shutil


def down_seq_files(db, db_version, download_dir, task_id, sample_list):
    """
    功能是根据任务下载组装序列文件
    :param download_dir:
    :param task_id:
    :return:
    """
    if db in ['bacgenome']:
        client = Config().get_mongo_client(mtype=db, db_version=db_version)
        db = client[Config().get_mongo_dbname(db, db_version=db_version)]
        if not os.path.exists(download_dir):
            os.mkdir(download_dir)
        else:
            shutil.rmtree(download_dir)
            os.mkdir(download_dir)
        assemble = db.assemble
        assemble_seq = db.assemble_seq
        assemble_id = assemble.find_one({"task_id": task_id})["main_id"]
        analysis_type = assemble.find_one({"task_id": task_id})["analysis_type"]
        if analysis_type == "uncomplete":
            for each in sample_list.keys():
                path = assemble_seq.find_one({"assemble_id": assemble_id, "specimen_id": each})["seq_path"]
                s3_path = path + "/" + each + "/assembly_predict/assembly/" + each + "_scaf.fna"
                new_path = os.path.join(download_dir, sample_list[each] + ".fna")
                download(s3_path, new_path)
        else:
            for sample in sample_list.keys():
                complete_list = []
                complete_list_uniq = []
                if not os.path.exists(download_dir + '/' + sample_list[sample]):
                    os.mkdir(download_dir + '/' + sample_list[sample])
                path = assemble_seq.find_one({"assemble_id": assemble_id, "specimen_id": sample})["seq_path"]
                sample_data = assemble_seq.find({"assemble_id": assemble_id, "specimen_id": sample})
                for i in sample_data:
                    complete_list.append(i)
                for each in complete_list:
                    if each:
                        if each["type"] not in complete_list_uniq:
                            complete_list_uniq.append(each["type"])
                for a in complete_list_uniq:
                    s3_path = path + "/" + sample + "/assembly_predict/assembly/seq_dir/" + a + ".fasta"
                    new_path = download_dir + '/' + sample_list[sample] + '/' + sample_list[sample] + '_' + a + '.fna'
                    download(s3_path, new_path)
        return (download_dir, analysis_type)


def down_gbk_files(db, db_version, download_dir, task_id, sample_list):
    """
    功能是根据任务下载gbk文件
    :param download_dir:
    :param task_id:
    :return:
    """
    if db in ['fungi','fungigenome']:
        db = "fungigenome"
        client = Config().get_mongo_client(mtype=db, db_version=db_version)
        db = client[Config().get_mongo_dbname(db, db_version=db_version)]
        if not os.path.exists(download_dir):
            os.mkdir(download_dir)
        else:
            shutil.rmtree(download_dir)
            os.mkdir(download_dir)
        assemble = db.assemble
        assemble_seq = db.assemble_seq
        assemble_id = assemble.find_one({"task_id": task_id})["_id"]
        analysis_type = assemble.find_one({"task_id": task_id})["analysis_type"]
        for sample in sample_list.keys():
            complete_list = []
            complete_list_uniq = []
            if not os.path.exists(download_dir + '/' + sample_list[sample]):
                os.mkdir(download_dir + '/' + sample_list[sample])
            path = assemble_seq.find_one({"assemble_id": assemble_id, "specimen_id": sample})["seq_path"]
            sample_data = assemble_seq.find({"assemble_id": assemble_id, "specimen_id": sample})
            for i in sample_data:
                complete_list.append(i)
            for each in complete_list:
                if each:
                    if each["type"] not in complete_list_uniq:
                        complete_list_uniq.append(each["type"])
            for a in complete_list_uniq:
                s3_path = path + "/" + sample + "/project_overview/All_predict_info/" + sample + ".gbk"
                new_path = download_dir  + '/' + sample_list[sample] + '.gbk'
                download(s3_path, new_path)
        return (download_dir, analysis_type)