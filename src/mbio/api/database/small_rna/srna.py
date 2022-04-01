# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

from mbio.api.database.small_rna.api_base import ApiBase
import json
import pandas as pd
import datetime
from bson.objectid import ObjectId
import os
import numpy as np
from collections import OrderedDict
import math
from scipy import stats
import unittest

class Srna(ApiBase):
    def __init__(self, bind_object):
        super(Srna, self).__init__(bind_object)

    def add_known_mirna(self, known_mirna_detail, project_sn="small_rna", task_id="small_rna", params=None, pdfs=None):
        self._bind_object.logger.info('start to add known mirna info to mongo')
        self._bind_object.logger.debug('incoming known mirna in add_known_mirna is {}'.format(known_mirna_detail))
        known_mirna_pd = pd.read_table(known_mirna_detail)

        # check params type
        if params == None:
            self._bind_object.logger.debug('type of incoming params is {}, make an empty str'.format(type(params)))
            params = str()
        elif type(params) == dict:
            self._bind_object.logger.debug('type of incoming params is {}, dumps it to str'.format(type(params)))
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        elif type(params) == str:
            self._bind_object.logger.debug('type of incoming params is {}, pass'.format(type(params)))
        else:
            self._bind_object.set_error('type of incoming params is {}, abord'.format(type(params)))

        # prepare main_info value
        time_now = datetime.datetime.now()
        name = 'KnownMirna_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        desc = 'known_mirna_main_table'
        # insert one document to sg_known_mirna, return _id as main_id
        main_info = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name,
            'created_ts': created_ts,
            'desc': desc,
            'params': params,
            'pdf_dir': pdfs,
            'status': 'start'
        }
        main_id = self.create_db_table('sg_known_mirna', [main_info])
        self._bind_object.logger.info('succeed in creating document in sg_known_mirna')

        # prepare record_dict
        known_mirna_pd['known_mirna_id'] = ObjectId(main_id)
        row_dict_list = known_mirna_pd.to_dict('records')

        # insert many documents to sg_known_mirna_detail
        self.create_db_table('sg_known_mirna_detail', row_dict_list)
        self._bind_object.logger.info('succeed in creating documents in sg_known_mirna_detail')

        # upsert status and main_id
        record_dict = {'_id': main_id, 'task_id': task_id}
        self.update_db_record_by_dict('sg_known_mirna', record_dict, status='end', main_id=main_id)
        return main_id

    def add_ncrna_stat(self, ncrna_stat, project_sn="small_rna", task_id="small_rna", params=None):
        self._bind_object.logger.info('start to add ncrna info to mongo')
        self._bind_object.logger.debug('incoming known mirna in add_ncrna_stat is {}'.format(ncrna_stat))

        # check params type
        if params == None:
            self._bind_object.logger.debug('type of incoming params is {}, make an empty str'.format(type(params)))
            params = str()
        elif type(params) == dict:
            self._bind_object.logger.debug('type of incoming params is {}, dumps it to str'.format(type(params)))
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        elif type(params) == str:
            self._bind_object.logger.debug('type of incoming params is {}, pass'.format(type(params)))
        else:
            self._bind_object.set_error('type of incoming params is {}, abord'.format(type(params)))

        # prepare main_info value
        time_now = datetime.datetime.now()
        name = 'NcrnaStat_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        desc = 'ncrna_stat_main_table'
        # insert one document to sg_known_mirna, return _id as main_id
        with open(ncrna_stat, "r") as f:
            head = f.readline()
            self.dict_ncrna = {}
            self.samples = head.strip().split("\t")[1:]
            self.sample_order = ",".join(self.samples)
            for sample in self.samples:
                if sample not in self.dict_ncrna.keys():
                    self.dict_ncrna[sample] = {}
                    self.dict_ncrna[sample]["rRNA"] = 0
                    self.dict_ncrna[sample]["snRNA"] = 0
                    self.dict_ncrna[sample]["snoRNA"] = 0
                    self.dict_ncrna[sample]["tRNA"] = 0
                    self.dict_ncrna[sample]["repbase"] = 0
                    self.dict_ncrna[sample]["exon"] = 0
                    self.dict_ncrna[sample]["intron"] = 0
                    self.dict_ncrna[sample]["unannotated"] = 0
            for line in f:
                items = line.strip().split("\t")
                if items[0] == "rRNA":
                    for index, value in enumerate(items[1:]):
                        self.dict_ncrna[self.samples[index]]["rRNA"] += int(value.split("(")[0])
                if items[0] == "tRNA":
                    for index, value in enumerate(items[1:]):
                        self.dict_ncrna[self.samples[index]]["tRNA"] += int(value.split("(")[0])
                if items[0] == "snRNA":
                    for index, value in enumerate(items[1:]):
                        self.dict_ncrna[self.samples[index]]["snRNA"] += int(value.split("(")[0])
                if items[0] == "snoRNA":
                    for index, value in enumerate(items[1:]):
                        self.dict_ncrna[self.samples[index]]["snoRNA"] += int(value.split("(")[0])
                if items[0] == "repbase":
                    for index, value in enumerate(items[1:]):
                        self.dict_ncrna[self.samples[index]]["repbase"] += int(value.split("(")[0])
                if items[0] == "exon":
                    for index, value in enumerate(items[1:]):
                        self.dict_ncrna[self.samples[index]]["exon"] += int(value.split("(")[0])
                if items[0] == "intron":
                    for index, value in enumerate(items[1:]):
                        self.dict_ncrna[self.samples[index]]["intron"] += int(value.split("(")[0])
                if items[0] == "unannotated":
                    for index, value in enumerate(items[1:]):
                        self.dict_ncrna[self.samples[index]]["unannotated"] += int(value.split("(")[0])

        main_info = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name,
            'created_ts': created_ts,
            'desc': desc,
            'params': params,
            'status': 'start',
            'sample_order': self.sample_order
        }
        main_id = self.create_db_table('sg_ncrna_stat', [main_info])
        self._bind_object.logger.info('succeed in creating document in sg_ncrna_stat')

        # prepare record_dict
        ncrna_stat = pd.read_table(ncrna_stat)
        ncrna_stat['ncrna_stat_id'] = ObjectId(main_id)
        row_dict_list = ncrna_stat.to_dict('records')

        # insert many documents to sg_ncrna_stat_detail
        self.create_db_table('sg_ncrna_stat_detail', row_dict_list)
        self._bind_object.logger.info('succeed in creating documents in sg_ncrna_stat_detail')

        # insert documents to sg_ncrna_stat_graph
        data_list = []
        for sample in self.samples:
            data = {
                "sample": sample,
                "rRNA" : self.dict_ncrna[sample]['rRNA'],
                "snRNA": self.dict_ncrna[sample]['snRNA'],
                "snoRNA": self.dict_ncrna[sample]['snoRNA'],
                "tRNA": self.dict_ncrna[sample]['tRNA'],
                "repbase": self.dict_ncrna[sample]['repbase'],
                "exon": self.dict_ncrna[sample]['exon'],
                "intron": self.dict_ncrna[sample]['intron'],
                "unannotated": self.dict_ncrna[sample]['unannotated'],
                "ncrna_stat_id": ObjectId(main_id)
                }
            data_list.append(data)
        try:
            collection = self.db["sg_ncrna_stat_graph"]
            result = collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入ncrna统计画图信息数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入ncrna统计画图信息数据成功")

        # upsert status and main_id
        record_dict = {'_id': main_id, 'task_id': task_id}
        self.update_db_record_by_dict('sg_ncrna_stat', record_dict, status='end', main_id=main_id)
        return main_id

    def add_novel_mirna(self, novel_mirna_detail, project_sn="small_rna", task_id="small_rna", params=None, category=None, pdfs=None):
        self._bind_object.logger.info('start to add novel mirna info to mongo')
        self._bind_object.logger.debug('incoming novel mirna in add_novel_mirna is {}'.format(novel_mirna_detail))
        novel_mirna_pd = pd.read_table(novel_mirna_detail)

        # check params type
        if params == None:
            self._bind_object.logger.debug('type of incoming params is {}, make an empty str'.format(type(params)))
            params = str()
        elif type(params) == dict:
            self._bind_object.logger.debug('type of incoming params is {}, dumps it to str'.format(type(params)))
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        elif type(params) == str:
            self._bind_object.logger.debug('type of incoming params is {}, pass'.format(type(params)))
        else:
            self._bind_object.set_error('type of incoming params is {}, abord'.format(type(params)))

        # prepare main_info value
        time_now = datetime.datetime.now()
        if category.lower() == "plant":
            method = "mireap"
            score_energy = 'energy'
        else:
            method = "mirdeep"
            score_energy = 'score'
        name = 'NovelMirna_{}_{}'.format(method, time_now.strftime('%Y%m%d_%H%M%S'))
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        desc = 'novel_mirna_main_table'
        # insert one document to sg_known_mirna, return _id as main_id
        main_info = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name,
            'created_ts': created_ts,
            'desc': desc,
            'params': params,
            'pdf_dir': pdfs,
            'status': 'start',
            'score_energy': score_energy
        }
        main_id = self.create_db_table('sg_novel_mirna', [main_info])
        self._bind_object.logger.info('succeed in creating document in sg_novel_mirna')

        # prepare record_dict
        novel_mirna_pd['novel_mirna_id'] = ObjectId(main_id)
        row_dict_list = novel_mirna_pd.to_dict('records')

        # insert many documents to sg_novel_mirna_detail
        self.create_db_table('sg_novel_mirna_detail', row_dict_list)
        self._bind_object.logger.info('succeed in creating documents in sg_novel_mirna_detail')

        # upsert status and main_id
        record_dict = {'_id': main_id, 'task_id': task_id}
        self.update_db_record_by_dict('sg_novel_mirna', record_dict, status='end', main_id=main_id)
        return main_id

    def add_mirna_stat(self, mirna_stat, project_sn="small_rna", task_id="small_rna", params=None, group=None):
        self._bind_object.logger.info('start to add mirna info to mongo')
        self._bind_object.logger.debug('incoming mirna stat in add_mirna_stat is {}'.format(mirna_stat))
        if group:
            sample_list = list()
            with open(group, 'r') as g:
                for line in g.readlines():
                    if line.startswith('#'):
                        continue
                    sample_list.append(line.strip().split('\t')[0])
        mirna_stat = pd.read_table(mirna_stat)

        # check params type
        if params == None:
            self._bind_object.logger.debug('type of incoming params is {}, make an empty str'.format(type(params)))
            params = str()
        elif type(params) == dict:
            self._bind_object.logger.debug('type of incoming params is {}, dumps it to str'.format(type(params)))
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        elif type(params) == str:
            self._bind_object.logger.debug('type of incoming params is {}, pass'.format(type(params)))
        else:
            self._bind_object.set_error('type of incoming params is {}, abord'.format(type(params)))

        # prepare main_info value
        time_now = datetime.datetime.now()
        name = 'MirnaStat_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        desc = 'mirna_stat_main_table'
        # insert one document to sg_known_mirna, return _id as main_id
        main_info = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name,
            'created_ts': created_ts,
            'desc': desc,
            'params': params,
            'status': 'start'
        }
        main_id = self.create_db_table('sg_mirna_stat', [main_info])
        self._bind_object.logger.info('succeed in creating document in sg_mirna_stat')

        # prepare record_dict
        mirna_stat['mirna_stat_id'] = ObjectId(main_id)
        row_dict_list = mirna_stat.to_dict('records')
        if group:
            row_dict_list.sort(key=lambda x: sample_list.index(x['sample']))
        # insert many documents to sg_mirna_stat_detail
        self.create_db_table('sg_mirna_stat_detail', row_dict_list)
        self._bind_object.logger.info('succeed in creating documents in sg_mirna_stat_detail')

        # upsert status and main_id
        record_dict = {'_id': main_id, 'task_id': task_id}
        self.update_db_record_by_dict('sg_mirna_stat', record_dict, status='end', main_id=main_id)
        return main_id

    def add_srna_stat(self, srna_stat, srna_stat_for_graph, project_sn="small_rna", task_id="small_rna", params=None, group=None):
        self._bind_object.logger.info('start to add srna info to mongo')
        self._bind_object.logger.debug('incoming srna stat in add_srna_stat is {}'.format(srna_stat))
        if group:
            sample_list = list()
            with open(group, 'r') as g:
                for line in g.readlines():
                    if line.startswith('#'):
                        continue
                    sample_list.append(line.strip().split('\t')[0])
        srna_stat = pd.read_table(srna_stat)

        # check params type
        if params == None:
            self._bind_object.logger.debug('type of incoming params is {}, make an empty str'.format(type(params)))
            params = str()
        elif type(params) == dict:
            self._bind_object.logger.debug('type of incoming params is {}, dumps it to str'.format(type(params)))
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        elif type(params) == str:
            self._bind_object.logger.debug('type of incoming params is {}, pass'.format(type(params)))
        else:
            self._bind_object.set_error('type of incoming params is {}, abord'.format(type(params)))

        # prepare main_info value
        time_now = datetime.datetime.now()
        name = 'SrnaStat_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        desc = 'srna_stat_main_table'
        # insert one document to sg_known_mirna, return _id as main_id
        main_info = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name,
            'created_ts': created_ts,
            'desc': desc,
            'params': params,
            'status': 'start',
            'sample_order': sample_list
        }
        main_id = self.create_db_table('sg_srna_stat', [main_info])
        self._bind_object.logger.info('succeed in creating document in sg_srna_stat')

        # prepare record_dict
        srna_stat['srna_stat_id'] = ObjectId(main_id)
        row_dict_list = srna_stat.to_dict('records')

        # insert many documents to sg_srna_stat_detail
        self.create_db_table('sg_srna_stat_detail', row_dict_list)
        self._bind_object.logger.info('succeed in creating documents in sg_srna_stat_detail')

        # insert stat graph
        with open(srna_stat_for_graph, "r") as f:
            data_list = []
            first_line = f.readline()
            for line in f:
                line = line.strip().split("\t")
                data = {
                    "srna_stat_id" : ObjectId(main_id),
                    "sample": line[0],
                    "genome": {"repbase" : int(line[7]), "exon" : int(line[8]), "intron" : int(line[9])},
                    "miRNA": {"known_mirna": int(line[1]), "novel_mirna": int(line[2])},
                    "sncRNA": {"rRNA" : int(line[3]), "snoRNA" : int(line[4]), "snRNA" : int(line[5]), "tRNA" : int(line[6])},
                    "unknown": {"unknown": int(line[10])},
                    }
                data_list.append(data)
            try:
                collection = self.db["sg_srna_stat_graph"]
                result = collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入srna统计画图信息数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入srna统计画图信息数据成功")

        # upsert status and main_id
        record_dict = {'_id': main_id, 'task_id': task_id}
        self.update_db_record_by_dict('sg_srna_stat', record_dict, status='end', main_id=main_id)
        return main_id

    ####################################################################################################

class TestFunction(unittest.TestCase):
    """
    测试导表函数
    """

    def test_mongo(test):
        from mbio.workflows.small_rna.small_rna_test_api import SmallRnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        data = {
        "id": "small_rna",
        "project_sn": "small_rna",
        "type": "workflow",
        "name": "small_rna.small_rna_test_api",
        "options": {},
        }
        wsheet = Sheet(data=data)
        wf = SmallRnaTestApiWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        known_mirna = "/mnt/ilustre/users/sanger-dev/workspace/20181127/Single_srna17-28-32/Srna/output/known_mirna/known_mirna_detail.xls"
        novel_mirna = "/mnt/ilustre/users/sanger-dev/workspace/20181127/Single_srna17-28-32/Srna/output/novel_mirna/novel_mirna_detail.xls"
        pdfs_known = "/mnt/ilustre/users/sanger-dev/workspace/20181127/Single_srna17-28-32/Srna/output/known_mirna/structure_pdf"
        pdfs_novel = "/mnt/ilustre/users/sanger-dev/workspace/20181127/Single_srna17-28-32/Srna/output/novel_mirna/structure_pdf"
        ncrna_stat = "/mnt/ilustre/users/sanger-dev/workspace/20181127/Single_srna17-28-32/Srna/SrnaStat/output/ncrna_stat.xls"
        srna_stat = "/mnt/ilustre/users/sanger-dev/workspace/20181127/Single_srna17-28-32/Srna/SrnaStat/output/srna_stat.xls"
        mirna_stat = "/mnt/ilustre/users/sanger-dev/workspace/20181127/Single_srna17-28-32/Srna/SrnaStat/output/mirna_stat.xls"
        srna_stat_for_graph = "/mnt/ilustre/users/sanger-dev/workspace/20181127/Single_srna17-28-32/Srna/SrnaStat/output/srna_stat_for_graph.xls"
        category = "animal"
        wf.api_srna = wf.api.api("small_rna.srna")
        wf.api_srna.add_known_mirna(known_mirna, project_sn="small_rna", task_id="small_rna", params={"method": "quantifier"}, pdfs=pdfs_known)
        wf.api_srna.add_ncrna_stat(ncrna_stat, project_sn="small_rna", task_id="small_rna", params={"method": "statistics"})
        wf.api_srna.add_novel_mirna(novel_mirna, project_sn="small_rna", task_id="small_rna", params={"method": "mirdeep"}, category=category ,pdfs=pdfs_novel)
        wf.api_srna.add_mirna_stat(mirna_stat, project_sn="small_rna", task_id="small_rna", params={"method": "mirdeep"})
        wf.api_srna.add_srna_stat(srna_stat, srna_stat_for_graph, project_sn="small_rna", task_id="small_rna", params={"method": "mirdeep"})

if __name__ == '__main__':
    unittest.main()