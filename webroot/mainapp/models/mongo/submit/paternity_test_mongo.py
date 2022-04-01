# -*- coding: utf-8 -*-
# __author__ = 'moli.zhou'
import datetime
from mainapp.libs.param_pack import param_pack
from bson import SON
from bson import ObjectId
from biocluster.config import Config
from mainapp.models.mongo.core.base import Base


class PaternityTest(Base):
    """
    写入流程主表，从任务主表中提取父母本编号信息等
    """

    def __init__(self):
        # self.mongo_client = Config().mongo_client
        # self.database = self.mongo_client[Config().MONGODB+'_paternity_test']
        # self.mongo_client_get_tab = Config().biodb_mongo_client
        # self.database_tab = self.mongo_client_get_tab['sanger_paternity_test_ref']
        super(PaternityTest, self).__init__()
        self._project_type = "pt"
        # self.client = Config().get_mongo_client(mtype=self._project_type)
        # self.database = self.client[Config().get_mongo_dbname(self._project_type)]
        # self.client_ref = Config().get_mongo_client(mtype=self._project_type, ref=True)
        # self.database_ref = self.client_ref[Config().get_mongo_dbname(self._project_type, ref=True)]

    def add_pt_father(self, father_id, err_min, dedup):
        params = dict()
        params['err_min'] = err_min
        params['dedup'] = dedup
        name = 'err-' + str(err_min) + '_dedup-' + str(dedup)
        insert_data = {
            "father_id": father_id,
            "name": name,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        # collection = self.database['sg_pt_father']
        collection = self.db['sg_pt_father']
        new_params = param_pack(params)
        insert_data["params"] = new_params
        pt_father_id = collection.insert_one(insert_data).inserted_id
        return pt_father_id

    def get_query_info(self, task):
        if task is None:
            raise Exception('未获取到father_id')
        # collection = self.database['sg_father']
        collection = self.db['sg_father']
        task_info = collection.find({"_id": ObjectId(task)})
        for info in task_info:
            # 	print i
            # #新建一个字典包含所需信息，然后返回这个字典
            # 	info = dict(
            # 		dad_id=i[u'dad_id'],
            # 		mom_id=i[u'mom_id'],
            # 		preg_id=i[u'preg_id'],
            # 		ref_fasta=i[u'ref_fasta'],
            # 		targets_bedfile=i[u'targets_bedfile'],
            # 		ref_point=i[u'ref_point'],
            # 		project_sn=i[u'project_sn'],
            # 		fastq_path=i[u'fastq_path']
            # 	)
            return info

    def get_ref_info(self, father_id):
        if father_id is None:
            raise Exception('未获取到father_id')
        # collection = self.database['sg_pt_ref_file']
        collection = self.db['sg_pt_ref_file']
        task_info = collection.find_one({"father_id": ObjectId(father_id)})
        return task_info

    def insert_main_table(self, collection, data):
        # return self.database[collection].insert_one(SON(data)).inserted_id
        return self.db[collection].insert_one(SON(data)).inserted_id
