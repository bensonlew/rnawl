# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modify:20171127
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import unittest
import json

from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from mbio.api.database.denovo_rna_v2.api_assembly import ApiBase


class DenovoassAlign(ApiBase):
    def __init__(self, bind_object):
        super(DenovoassAlign, self).__init__(bind_object)
        self.result_file = {}
        self.result_dir = ''

    def set_result_dir(self, align_mudule_dir):
        '''
        根据组装模块结果导入结果路径
        '''
        self.result_dir = align_mudule_dir
        self.bind_object.logger.info("导入**** {}".format(align_mudule_dir))
        self.result_file['align_stat'] = os.path.join(align_mudule_dir, "alignment_rate.txt")

        for key, value in self.result_file.items():
            if os.path.exists(value):
                pass
            else:
                self.bind_object.set_error('结果文件%s 不存在，请检查', variables=(value), code="52002101")

    def check_id(self, object_id):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('onject_id必须为ObjectId对象或其对应的字符串！', code="52002102")
        return object_id

    @report_check
    def add_align(self):
        '''
        插入主表
        '''
        params =  {"method":"Trinity"}

        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': 'Align_salmon_' + datetime.datetime.now().strftime('%Y%m%d_%H%M%S'),
            'desc': 'denovo比对主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            # 'result_dir': self.result_dir,
            # 'type': 'workflow',
            'status': 'end',
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
        }
        collection = self.db['sg_align']
        assemble_id = collection.insert_one(insert_data).inserted_id
        return assemble_id

    @report_check
    def add_align_stat(self, align_id):
        '''
        插入比对详情表
        '''
        align_id = self.check_id(align_id)
        data_list = list()

        with open(self.result_file['align_stat'], 'rb') as f:
            lines = f.readlines()
            '''
            '''
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = [
                    ('align_id', align_id),
                    ('sample', line[0]),
                    ('total_reads', int(line[1])),
                    ('aligned_reads', int(line[2])),
                    ('aligned_rate', format(float(line[3]), '.2%')),
                ]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db["sg_align_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入%s信息出错:%s" % ("sg_align_detail", e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % ("sg_align_detail"))

    def run(self, result_dir):
        self.bind_object.logger.info("开始导表、数据路径为 {}".format(result_dir))
        self.set_result_dir(result_dir)
        obj_id = self.add_align()
        self.add_align_stat(obj_id)
        self.update_db_record('sg_align', obj_id, status="end", main_id=obj_id)

class TestFunction(unittest.TestCase):
    """
    测试导表函数
    """
    from mbio.workflows.denovo_rna_v2.denovo_test_api import DenovoTestApiWorkflow
    from biocluster.wsheet import Sheet
    import random

    def test_mongo(self):
        print("test12212112*")
        data = {
            "id": "denovo_assemble",
            "project_sn": "denovo_assemble",
            "type": "workflow",
            "name": "denovo_rna_v2.denovo_test_api",
            "options": {
            },
        }
        wsheet = self.Sheet(data=data)
        wf = self.DenovoTestApiWorkflow(wsheet)

        test_dir = '/mnt/ilustre/users/sanger-dev/workspace/20190618/DenovoAssemble_denovo_ass6741/Quant/'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        print("TEST fUNCTION")
        wf.test_api = wf.api.api("denovo_rna_v2.denovoass_align")
        wf.test_api.run(test_dir)

if __name__ == '__main__':
    unittest.main()
