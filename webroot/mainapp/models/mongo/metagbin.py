# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'

import types
import json
import random, os, re
from bson import SON
from bson.objectid import ObjectId
from biocluster.config import Config
from mainapp.models.workflow import Workflow
from mainapp.models.mongo.meta import Meta
from biocluster.file import exists,list_dir


class Metagbin(Meta):
    """
    获取mongo字段
    """
    def __init__(self, bind_object=None):
        self._bind_object = bind_object
        super(Metagbin, self).__init__(self._bind_object)
        self._project_type = 'metagbin'

    def get_sample_genefile(self, task_id, genome, predict="cds", type='fnn'):
        """
        获取预测基因的路径或16s_rRNA的gene路径
        :param task_id:
        :param type:
        :return:
        """
        if predict not in ["cds", "rrna"]:
            raise Exception("输入predict必须为gene, trna, rrna!")
        if type not in ['fnn', 'faa']:
            raise Exception("输入type必须为fnn, faa!")
        collection = self.db['predict_' + predict]
        result = collection.find_one({"task_id": task_id, "genome_id": genome})
        path = result[type + '_path']
        return path

    def get_genome(self, task_id, genome):
        """
        获取基因组对象存储路径
        :param task_id:
        :param genome:
        :return:
        """
        collection = self.db["assembly"]
        result = collection.find_one({"task_id": task_id, "genome_id": genome})
        path = result["genome_path"]
        return path

    def get_projectsn(self, task_id):
        """
        获取project_sn
        :param task_id:
        :return:
        """
        collection = self.db['sg_task']
        result = collection.find_one({"task_id": task_id})
        if not result:
            raise Exception('{}没有响应的id:{}!'.format(collection, task_id))
        project_sn = result["project_sn"]
        return project_sn

    def get_insert(self,task_id, genome_id):
        """
        根据task_id和genome_id去mongo中查assembly主表的插入片段大小和最大读长
        :return:
        """
        main_collection = self.db['assembly']
        result = main_collection.find_one({'task_id':task_id, 'genome_id': genome_id})
        if not result:
            raise Exception('{}没有响应的id：{}'.format(main_collection, task_id))
        insert_size = int(result['insert'])
        max_length = int(result['max'])
        return insert_size, max_length

    def get_taxon(self, task_id, bin_id):
        """
        根据task_id查找main_id，然后根据main_id和bin_id提取出来是古菌还是细菌
        :param task_id:
        :param bin_id:
        :return:
        """
        main_collection = self.db['bin']
        main = main_collection.find_one({"task_id": task_id})
        main_id = main['main_id']
        bin_detail = self.db['bin_detail']
        new_bin_id = str(bin_id).split("_")[1]
        result = bin_detail.find_one({"bi_id":main_id, "bin_id": new_bin_id})
        domain = result['domain']
        return domain











