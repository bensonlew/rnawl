# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last modify: 2018.3.5

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re,os
from types import StringTypes
from biocluster.config import Config
from bson.son import SON
import shutil
import glob


class Gbk(Base):
    def __init__(self, bind_object):
        super(Gbk, self).__init__(bind_object)
        self._project_type = "bacgenome"


    @report_check
    def add_gbk(self,params=None, project_sn=None, task_id=None):
        task_id = task_id if task_id else self.bind_object.sheet.id
        project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
        collection = self.db["gbk"]
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "faile",
            "desc": "gbk分析",
            "params": {'ananlysis':'gbk'},
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "version" : "3.0"
        }
        main_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': main_id},{'$set':{'main_id':main_id}}) #guanqing 20180813
        return main_id

    @report_check
    def add_gbk_detail(self,update=True,path=None,main_id =None,sample_name=None,gbk_path =None,new_path =None,task_id=None,gbk_detail_id=None):
        if update:
            files = os.listdir(path)
            data_list = []
            path_up_1 =  os.path.dirname(path.rstrip('/'))
            for file in files:
                if file in ['Scaffold']:
                    data = {
                        "gbk_id": ObjectId(main_id),
                        "specimen_id": sample_name,
                        "gbk_path": sample_name+'.gbk',
                        "ptt_path" : sample_name +'.ptt',
                        "gff_path" : sample_name + '.gff',
                        "path":gbk_path,
                    }
                    ori_params = self.add_gbk_pre_info(path+'/Scaffold/'+sample_name+'.gbk')
                    data['ori_params'] = ori_params
                    data_son = SON(data)
                    data_list.append(data_son)
                else:
                    seq_type = file
                    if 'Chromosome' in file:
                        genome_type = 'chromosome'
                    elif 'Plasmid' in file:
                        genome_type = 'plasmid'
                    else:
                        genome_type = re.sub("\d*$","",file)
                    data = {
                        "gbk_id": ObjectId(main_id),
                        "specimen_id": sample_name,
                        #"gbk_path": file,
                        "path": gbk_path,
                        "genome_type" : genome_type,
                        "seq_type" : seq_type
                    }
                    gbks = glob.glob(path + '/' + file + '/*.gbk')
                    for g in gbks:
                        g_name = os.path.basename(g)
                        data['gbk_path'] = g_name
                        ori_params = self.add_gbk_pre_info(g)
                        data['ori_params'] = ori_params
                        g_pre = g_name.rstrip('.gbk')
                        if os.path.exists(path_up_1+'/'+ g_pre+'.gff'):
                            data['gff_path'] = g_pre+'.gff'
                        else:
                            data['gff_path'] = '-'

                        if os.path.exists(path_up_1+'/'+ g_pre+'.ptt'):
                            data['ptt_path'] = g_pre+'.ptt'
                        else:
                            data['ptt_path'] = '-'

                        data_son = SON(data)
                        data_list.append(data_son)



            # ## 增加ptt和gff文件路径
            # ptt_data = {
            #     "gbk_id": ObjectId(main_id),
            #     "specimen_id": sample_name,
            #     "path": gbk_path,
            #     "ptt_path" : sample_name + '.ptt'
            # }
            #
            # gff_data = {
            #     "gbk_id": ObjectId(main_id),
            #     "specimen_id": sample_name,
            #     "path": gbk_path,
            #     "gff_path" : sample_name+'.gff'
            # }
            #
            # data_list.append(SON(ptt_data))
            # data_list.append(SON(gff_data))
            try:
                collection = self.db["gbk_detail"]
                collection.insert_many(data_list)
                main_collection = self.db["gbk"]
                main_collection.update({"_id": ObjectId(main_id)},
                                       {"$set": {"status": 'end'}})
            except Exception, e:
                self.bind_object.logger.info("导入%s结果表出错:%s" % (path, e))
            else:
                self.bind_object.logger.info("导入%s结果表成功!" % path)
        else:
            collection = self.db["gbk_detail"]
            collection.update_one({"_id": ObjectId(gbk_detail_id)},
                                        {'$set': {'new_path': new_path}})


    def add_gbk_pre_info(self,gbk_file):
        with open(gbk_file) as f:
            lines = f.readlines()
        ori_params = {
        "line1" : lines[0].rstrip('\n'),
        "line3" : lines[2].rstrip('\n'),
        "line4" : lines[3].rstrip('\n'),
        "line5" : lines[4].rstrip('\n'),
        "line8" : lines[7].rstrip('\n'),
        "line9" : lines[8].rstrip('\n'),
        "des" : lines[1].rstrip('\n'),
        "source" : lines[5].rstrip('\n'),
        "organism" : lines[6].rstrip('\n'),
        "author" : lines[9].rstrip('\n'),
        "title" : lines[10].rstrip('\n'),
        "journal" : lines[11].rstrip('\n')
        }
        return ori_params
