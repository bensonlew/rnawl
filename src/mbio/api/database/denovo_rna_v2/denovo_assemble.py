# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modify:20171123
from biocluster.api.database.base import Base, report_check
import re
import os
import datetime
import types
import unittest
import json

from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from mbio.api.database.denovo_rna_v2.api_base import ApiBase


class DenovoAssemble(ApiBase):
    def __init__(self, bind_object):
        super(DenovoAssemble, self).__init__(bind_object)
        self._project_type = 'denovo_rna_v2'
        self.result_file = {}
        self.result_dir = ''
        self.filtered = 'yes'

    def set_result_dir(self, result_dir, result_unigene_dir, evolution_dir, evolution_unigene_dir):
        '''
        根据组装模块结果导入结果路径
        '''
        self.result_dir = result_dir
        self.bind_object.logger.info("导入**** {}".format(result_dir))
        self.result_file['assemble_stat'] = os.path.join(result_dir, "trinity_stat/Trinity_stat.xls")
        self.result_file['assemble_stat_unigene'] = os.path.join(result_unigene_dir, "trinity_stat/Trinity_stat.xls")
        self.result_file['assemble_length'] = os.path.join(result_dir, "trinity_stat/trans_count_stat_500.txt")
        self.result_file['assemble_length_gene'] = os.path.join(result_dir, "trinity_stat/unigene_count_stat_500.txt")
        if self.filtered == 'yes':
            self.result_file['filter_assemble_stat'] = os.path.join(evolution_dir, "trinity_stat/Trinity_stat.xls")
            self.result_file['filter_assemble_stat_unigene'] = os.path.join(evolution_unigene_dir, "trinity_stat/Trinity_stat.xls")
            self.result_file['filter_assemble_length'] = os.path.join(evolution_dir, "trinity_stat/trans_count_stat_500.txt")
            self.result_file['filter_assemble_length_gene'] = os.path.join(evolution_dir, "trinity_stat/unigene_count_stat_500.txt")
            #assemble_length = os.path.join(result_dir, "/trinity_stat/Trinity_stat.xls")
            #assemble_length_gene = os.path.join(result_dir, "/trinity_stat/Trinity_stat.xls")

        for key, value in self.result_file.items():
            if os.path.exists(value):
                pass
            else:
                self.bind_object.set_error('结果文件%s 不存在，请检查', variables=(value), code="52000601")

    def check_id(self, object_id):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('assemble_id必须为ObjectId对象或其对应的字符串！', code="52000602")
        return object_id

    @report_check
    def add_assemble(self):
        '''
        插入主表
        '''
        params =  {"method":"Trinity"}
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': 'Assemble_' + datetime.datetime.now().strftime('%Y%m%d_%H%M%S'),
            'desc': 'denovo组装拼接主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            #'result_dir': self.result_dir,
            #'type': 'workflow',
            'status': 'end',
            'filter': self.filtered,
            'version': 'v2',
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
        }
        collection = self.db['sg_assembly']
        assemble_id = collection.insert_one(insert_data).inserted_id
        return assemble_id

    @report_check
    def add_assemble_stat(self, assemble_id, table_type):
        '''
        插入组装统计详情表
        '''
        assemble_id = self.check_id(assemble_id)
        data_list = []
        stat_table = self.result_file['assemble_stat']
        stat_table_G = self.result_file['assemble_stat_unigene']
        if table_type == "filter":
            stat_table = self.result_file['filter_assemble_stat']
            stat_table_G = self.result_file['filter_assemble_stat_unigene']

        with open(stat_table, 'rb') as f, open(stat_table_G, 'rb') as f2:
            lines = f.readlines()
            f2.readline()
            '''
            Type    Resource
            Total transcripts num   2
            Total unigenes num      2
            Total sequence base     7196
            Largest 6908
            Smallest        288
            Average length  3598.00
            N50     6908
            E90N50  0
            GC percent      43.09
            Mean mapped reads      0
            TransRate score 0.17604
            Busco score     0.0%(0.0%)
            '''
            for line in lines[1:]:
                line = line.strip().split('\t')
                line_G = f2.readline().strip().split('\t')
                if line[0] == "BUSCO score":
                    if line_G[1].startswith("C:"):
                        stat_match = re.search(r'C:(.*)\[S:(.*),D:(.*)\],F:(.*),M:(.*),n:(.*)', line_G[1])
                        stat_c = stat_match.group(1)
                        stat_d = stat_match.group(3)
                        line_G[1] = stat_c + '(' + stat_d + ')'
                    if line[1].startswith("C:"):
                        stat_match = re.search(r'C:(.*)\[S:(.*),D:(.*)\],F:(.*),M:(.*),n:(.*)', line_G[1])
                        stat_c = stat_match.group(1)
                        stat_d = stat_match.group(3)
                        line[1] = stat_c + '(' + stat_d + ')'

                data = [
                    ('assembly_id', assemble_id),
                    ('type', table_type),
                    ('name', line[0]),
                    ('unigene', line_G[1]),
                    ('transcript', line[1])
                ]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db["sg_assembly_stat"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入%s信息出错:%s" % ("sg_assemble_stat", e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % ("sg_assemble_stat"))

    def add_length_stat(self, assemble_id, seq_type):
        '''
        插入长度统计详情表
        '''
        assemble_id = self.check_id(assemble_id)
        data_list = []
        step_num_list = []
        step_pct_list = []
        stat_table = ""

        if self.filtered == 'no':
            stat_table = self.result_file['assemble_length']
            if seq_type == "G":
                stat_table=self.result_file['assemble_length_gene']
            else:
                pass
        else:
            stat_table = self.result_file['filter_assemble_length']
            if seq_type == "G":
                stat_table=self.result_file['filter_assemble_length_gene']
            else:
                pass

        with open(stat_table, 'rb') as f:
            lines = f.readlines()
            num_sum = float(lines[-1].strip().split("\t")[1])
            for line in lines[1:-1]:
                step_num = dict()
                step_pct = dict()
                step_range = line.strip().split("\t")[0]
                num = line.strip().split("\t")[1]
                pct = format(float(num)/num_sum, '.0%')
                step_num[step_range] = num
                step_pct[step_range] = pct
                step_num_list.append(step_num)
                step_pct_list.append(step_pct)

            data = [
                ('assembly_id', assemble_id),
                ('seq_type', seq_type),
                ('step', 500),
                ('step_data', step_num_list),
                ('step_pct', step_pct_list)
            ]
            data = SON(data)
            data_list.append(data)

        try:
            collection = self.db["sg_assembly_len"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入%s信息出错:%s" % ("sg_assemble_len", e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % ("sg_assemble_len"))

    def run(self, result_dir, result_unigene_dir, evolution_dir=None, evolution_unigene_dir=None):
        if evolution_dir:
            pass
        else:
            self.filtered = "no"
        self.bind_object.logger.info("开始到表情数据路径为 {}".format(result_dir))
        self.set_result_dir(result_dir, result_unigene_dir, evolution_dir, evolution_unigene_dir)
        obj_id = self.add_assemble()
        self.add_assemble_stat(obj_id,  'assemble')
        if self.filtered == "yes":
            self.add_assemble_stat(obj_id, 'filter')
        self.add_length_stat(obj_id, 'G')
        self.add_length_stat(obj_id, 'T')
        self.update_db_record('sg_assembly', obj_id, status="end", main_id=obj_id)


class TestFunction(unittest.TestCase):
    """
    测试导表函数
    """
    from mbio.workflows.denovo_rna_v2.denovo_test_api import DenovoTestApiWorkflow
    from biocluster.wsheet import Sheet
    import random

    def test_mongo(self):
        data = {
            # "id": "denovo_rna_v2" + str(random.randint(1,10000)),
            "id": "denovo_rna_v2_upgrade",
            "project_sn": "denovo_rna_v2_upgrade",
            "type": "workflow",
            "name": "denovo_rna_v2.denovo_test_api",
            "options": {
            },
        }
        wsheet = self.Sheet(data=data)
        wf = self.DenovoTestApiWorkflow(wsheet)

        test_dir = '/mnt/ilustre/users/sanger-dev/workspace/20190801/DenovoAssemble_denovo_ass4855/DenovoAssemble2Filter/output'
        evolution_dir = '/mnt/ilustre/users/sanger-dev/workspace/20190801/DenovoAssemble_denovo_ass4855/AssembleEvalution/output'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("denovo_rna_v2.denovo_assemble")
        # wf.test_api.run(test_dir)
        wf.test_api.run(test_dir, test_dir, evolution_dir, evolution_dir)
        '''
        test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_denovo/test_data1/output'
        evolution_dir = '/mnt/ilustre/users/sanger-dev/workspace/20171227/Single_evolution_assemble6135/AssembleEvalution/output'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("denovo_rna_v2.denovo_assemble")
        wf.test_api.run(test_dir)
        wf.test_api.run(test_dir, evolution_dir)
        '''

if __name__ == '__main__':
    unittest.main()
