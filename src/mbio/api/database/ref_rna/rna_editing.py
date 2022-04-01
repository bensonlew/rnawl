# -*- coding: utf-8 -*-
# __author__ = 'moli.zhou'
import datetime
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
from mainapp.libs.param_pack import param_pack


class RnaEditing(Base):
    '''
    将前端需要调用的结果文件导入mongo数据库，之结果保存的tsanger collection
    '''

    def __init__(self, bind_object):
        super(RnaEditing, self).__init__(bind_object)
        # self._db_name = Config().MONGODB
        self._project_type = 'ref_rna'
        #self.mongo_client = Config().mongo_client
        #self.database = self.mongo_client['tsanger_ref_rna']

    # "status": "end",
    @report_check
    def add_editing(self, task_id, db):
        params = dict()
        params['database'] = db
        name = str(db)
        insert_data = {
            "task_id": task_id,
            "name": name,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }

        collection = self.db['sg_tf']
        new_params = param_pack(params)
        insert_data["params"] = new_params
        # collection.insert_data["params"] = params
        try:
            flow_id = collection.insert_one(insert_data).inserted_id
            # collection.insert_one(insert_data)
        except Exception as e:
            self.bind_object.logger.error('导入转录因子主表出错：{}'.format(e))
        else:
            self.bind_object.logger.info("导入转录因子主表成功")
        return flow_id

    @report_check
    def add_editing_overview(self, file, flow_id):
        self.bind_object.logger.info("开始转录因子细节表")
        sg_tf = list()
        with open(file, 'r') as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                insert_data = {
                    "task_id": self.bind_object.id,
                    "family_id": flow_id,
                    "amino_id": line[0],
                    "tf_family": line[1],
                    "DEG": line[2],
                }
                sg_tf.append(insert_data)
            try:
                collection = self.db['sg_tf_detail']
                collection.insert_many(sg_tf)
            except Exception as e:
                self.bind_object.logger.error('导入转录因子细节表出错：{}'.format(e))
            else:
                self.bind_object.logger.info("导入转录因子细节表成功")

    @report_check
    def add_editing_detail(self, file, flow_id):
        self.bind_object.logger.info("开始转录因子细节表")
        sg_tf = list()
