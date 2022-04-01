# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
from bson.objectid import ObjectId
import datetime
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.ref_rna_v2.api_base import ApiBase
from biocluster.config import Config
import re
import types
from mbio.packages.rna.annot_config import AnnotConfig

class AddAnnotation(ApiBase):
    def __init__(self, bind_object):
        super(AddAnnotation, self).__init__(bind_object)
        self._project_type = 'ref_rna_v2'
        self.kegg2020 = Config().SOFTWARE_DIR + '/database/Annotation/other2020/kegg202003/ko_des'
        self.kegg2019 = Config().SOFTWARE_DIR + '/database/Annotation/other2019/kegg201909/ko_des'

    def add_name_description(self, annotation_file, annot_type, annot_method, main_id=None, task_id='annotation', project_sn='annotation'):
        if main_id is None:
            name = "Add_annotation" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            # if type(params) == dict:
            #     params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v3",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='Add annotation main table',
                params=[],
                annot_type=annot_type,
                annot_method=annot_method,
                database = 'custom',
                status="start",
            )
            main_id = self.create_db_table('sg_add_annotation', [main_info])
        else:
            if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
                main_id = ObjectId(main_id)
        db = Config().get_mongo_client(mtype='ref_rna_v2')[Config().get_mongo_dbname('ref_rna_v2')]
        connect_main = db['sg_annotation_query']
        record = connect_main.find_one({'task_id': task_id})
        record_main_id = record['main_id']
        connect_detail = db['sg_annotation_query_detail']
        with open(annotation_file) as anno:
            for line in anno.readlines():
                gene_id, gene_name1, description1 = line.strip().split('\t')
                record2 = connect_detail.find({'query_id': record_main_id, 'gene_id':gene_id})
                for rc2 in record2:
                    if rc2 is None:
                        self.bind_object.set_error("{}不在query表中".format(gene_id))
                        continue
                    _id = rc2['_id']
                    gene_name_old = rc2['gene_name']
                    description_old = rc2['description']
                    if annot_type == 'replace':
                        # gene_name = gene_name
                        # description = description
                        if annot_method == 'all':
                            self.update_db_record_new('sg_annotation_query_detail', _id, gene_name=gene_name1,
                                                  description=description1)
                        if annot_method == 'name':
                            self.update_db_record_new('sg_annotation_query_detail', _id, gene_name=gene_name1)
                        if annot_method == 'description':
                            self.update_db_record_new('sg_annotation_query_detail', _id, description=description1)
                    if annot_type == 'supplement':
                        if rc2['gene_name'] == "":
                            gene_name = gene_name1
                        else:
                            gene_name = rc2['gene_name']
                        if rc2['description'] == "":
                            description = description1
                        else:
                            description = rc2['description']
                        if annot_method == 'all':
                            self.update_db_record_new('sg_annotation_query_detail', _id, gene_name=gene_name,
                                                  description=description)
                        if annot_method == 'name':
                            self.update_db_record_new('sg_annotation_query_detail', _id, gene_name=gene_name)
                        if annot_method == 'description':
                            self.update_db_record_new('sg_annotation_query_detail', _id, description=description)
                    try:
                        record['flag']
                    except:
                        self.update_db_record_new('sg_annotation_query_detail', _id, gene_name_old=gene_name_old, description_old=description_old)
        self.update_db_record('sg_add_annotation', main_id, status="end")
        self.update_db_record('sg_annotation_query', record_main_id, flag='add')
        self.delete_db_record_some_key('sg_annotation_query', record_main_id, flaging='running')


    def add_name_description_database(self, change_database,annot_type, annot_method, main_id=None, task_id='annotation',
                             project_sn='annotation'):
        if main_id is None:
            name = "Add_annotation" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            # if type(params) == dict:
            #     params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v3",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='Add annotation main table',
                params=[],
                # annot_type=annot_type,
                # annot_method=annot_method,
                # database=change_database,
                status="start",
            )
            main_id = self.create_db_table('sg_add_annotation', [main_info])
        else:
            if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
                main_id = ObjectId(main_id)
        db = Config().get_mongo_client(mtype=self._project_type)[Config().get_mongo_dbname(self._project_type)]

        connect_main = db['sg_annotation_query']
        record = connect_main.find_one({'task_id': task_id})
        record_main_id = record['main_id']
        connect_detail = db['sg_annotation_query_detail']
        record_detail = connect_detail.find({'query_id': record_main_id})
        try:
            if change_database.lower() == 'swiss-prot':
                for detail in record_detail:
                    _id = detail['_id']
                    gene_name_old = detail['gene_name']
                    description_old = detail['description']
                    try:
                        record['flag']
                    except:
                        self.update_db_record_new('sg_annotation_query_detail', _id, gene_name_old=gene_name_old)
                        self.update_db_record_new('sg_annotation_query_detail', _id,
                                              description_old=description_old)
                    swissprot = detail['swissprot']
                    if swissprot == "" or swissprot is None:
                        continue
                    else:
                        info = swissprot.strip().split('|')[2]
                    description_s = re.findall("\((.*?)\sOS=", info)[0]
                    name_s = re.findall("(.*?)\_",info)[0]
                    if annot_type == 'replace':
                        name = name_s
                        description = description_s
                    if annot_type == 'supplement':
                        if detail['gene_name'] == "":
                            name = name_s
                        else:
                            name = detail['gene_name']
                        if detail['description'] == "":
                            description = description_s
                        else:
                            description = detail['description']
                    if annot_method == 'all':
                        self.update_db_record_new('sg_annotation_query_detail', _id, gene_name=name
                                              )
                        self.update_db_record_new('sg_annotation_query_detail', _id,
                                              description=description)
                    if annot_method == 'name':
                        self.update_db_record_new('sg_annotation_query_detail', _id, gene_name=name)
                    if annot_method == 'description':
                        self.update_db_record_new('sg_annotation_query_detail', _id, description=description)
                    # try:
                    #     record['flag']
                    # except:
                    #     self.update_db_record('sg_annotation_query_detail', _id, gene_name_old=gene_name_old, description_old=description_old)
            if change_database.lower() == 'nr':
                for detail in record_detail:
                    _id = detail['_id']
                    gene_name_old = detail['gene_name']
                    description_old = detail['description']
                    try:
                        record['flag']
                    except:
                        self.update_db_record_new('sg_annotation_query_detail', _id, gene_name_old=gene_name_old, description_old=description_old)
                    nr = detail['nr']
                    if nr == "" or nr is None :
                        continue
                    if '[' in nr:
                        description_s = re.findall("\((.*?)\s\[", nr)[0]
                    else:
                        description_s = re.findall("\((.*?)\)", nr)[0]
                    if annot_type == 'replace':
                        description = description_s
                    if annot_type == 'supplement':
                        if detail['description'] == "":
                            description = description_s
                        else:
                            description = detail['description']
                    if annot_method == 'description':
                        self.update_db_record_new('sg_annotation_query_detail', _id, description=description)
                    # try:
                    #     record['flag']
                    # except:
                    #     self.update_db_record('sg_annotation_query_detail', _id, gene_name_old=gene_name_old, description_old=description_old)

            if change_database.lower() == "kegg":
                task = db['sg_task']
                task_record = task.find_one({'task_id': task_id})
                try:
                    database_version = task_record['database_version']
                    kegg_version = database_version['kegg']
                    kegg_files_dict = AnnotConfig().get_file_dict(db="kegg", version=kegg_version)
                    kegg = kegg_files_dict["ko_des"]

                    # if database_version['kegg'] == '2017':
                    #     kegg = self.kegg2019
                    # if database_version['kegg'] == '202003':
                    #     kegg = self.kegg2020

                except:
                    kegg = self.kegg2019
                ko_id_dict = dict()

                with open(kegg, 'r') as keggs:
                    for line in keggs.readlines():
                        ko_id, ko_info = line.strip().split('\t')
                        # name = re.findall("(.*?)\;", ko_info)
                        if ';' not in ko_info:
                            description = ko_info
                        else:
                            if '[E' in ko_info:
                                description = re.findall("\;\s(.*?)\s\[", ko_info)[0]
                            else:
                                description = re.findall("\;\s(.*?)$", ko_info)[0]
                        ko_id_dict[ko_id]=description
                for detail in record_detail:
                    _id = detail['_id']
                    gene_name_old = detail['gene_name']
                    description_old = detail['description']
                    ko_id = detail['ko_id']
                    ko_name = detail['ko_name']
                    if ko_id == "" or ko_name == "" or ko_id is None:
                        try:
                            record['flag']
                        except:
                            self.update_db_record_new('sg_annotation_query_detail', _id, gene_name_old=gene_name_old,
                                                  description_old=description_old)
                        continue
                    try:
                        record['flag']
                    except:
                        self.update_db_record_new('sg_annotation_query_detail', _id, gene_name_old=gene_name_old, description_old=description_old)
                    raw_ko_id = ko_id
                    raw_ko_name = ko_name
                    if ';' in raw_ko_id:
                        ko_id = ""
                        ko_name = ""
                        for n,i in enumerate(raw_ko_id.split(';')):
                            if i.startswith("K"):
                                ko_id = i
                                ko_name = raw_ko_name.split(';')[n]
                                break
                    try:
                        description_s = ko_id_dict['ko:{}'.format(ko_id)]
                    except:
                        description_s = ""
                    if annot_type == 'replace':
                        name = ko_name
                        description = description_s
                    if annot_type == 'supplement':
                        if detail['gene_name'] == "":
                            name = ko_name
                        else:
                            name = detail['gene_name']
                        if detail['description'] == "":
                            description = description_s
                        else:
                            description = detail['description']
                    if annot_method == 'all':
                        self.update_db_record_new('sg_annotation_query_detail', _id, gene_name=name,
                                              description=description)
                    if annot_method == 'name':
                        self.update_db_record_new('sg_annotation_query_detail', _id, gene_name=name)
                    if annot_method == 'description':
                        self.update_db_record_new('sg_annotation_query_detail', _id, description=description)
        except Exception, e:
            print e
            self.delete_db_record_some_key('sg_annotation_query', record_main_id, flaging='running')
            self.bind_object.set_error('更新gene_name和description出现错误，请运行一键还原')

                # try:
                #     record['flag']
                # except:
                #     self.update_db_record('sg_annotation_query_detail', _id, gene_name_old=gene_name_old, description_old=description_old)
        self.update_db_record('sg_add_annotation', main_id, status="end")
        self.update_db_record('sg_annotation_query', record_main_id, flag='add')
        self.delete_db_record_some_key('sg_annotation_query', record_main_id, flaging='running')
    def ghost_add_annotation(self, main_id=None, task_id='annotation', project_sn='annotation'):
        if main_id is None:
            name = "Ghost_Add_annotation" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            # if type(params) == dict:
            #     params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v3",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='Ghost Add annotation main table',
                params=[],
                # annot_type=annot_type,
                # annot_method=annot_method,
                # database=change_database,
                status="start",
            )
            main_id = self.create_db_table('sg_ghost_add_annotation', [main_info])
        else:
            if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
                main_id = ObjectId(main_id)
        db = Config().get_mongo_client(mtype='ref_rna_v2')[Config().get_mongo_dbname('ref_rna_v2')]
        connect_main = db['sg_annotation_query']
        record = connect_main.find_one({'task_id': task_id})
        record_main_id = record['main_id']
        connect_detail = db['sg_annotation_query_detail']
        record_detail = connect_detail.find({'query_id': record_main_id})

        for detail in record_detail:
            _id = detail['_id']
            try:
                gene_name_old = detail['gene_name_old']
                description_old = detail['description_old']
            except:
                continue
            self.update_db_record_new('sg_annotation_query_detail', _id, gene_name=gene_name_old,
                                  description=description_old)
            self.delete_db_record_some_key('sg_annotation_query_detail', _id, gene_name_old=gene_name_old,
                                           description_old=description_old)

        self.delete_db_record_some_key('sg_annotation_query', record_main_id, flag='add')
        self.delete_db_record_some_key('sg_annotation_query', record_main_id, flaging='ghosting')
        self.update_db_record('sg_ghost_add_annotation', main_id, status="end")
        self.remove_db_record('sg_add_annotation', task_id=task_id)

def update_db_record(self, table_name, record_id=None, query_dict=None, insert_dict=None, **kwargs):
    if record_id is not None:
        if isinstance(record_id, types.StringTypes):
            record_id = ObjectId(record_id)
        elif isinstance(record_id, ObjectId):
           record_id = record_id
        else:
            self.bind_object.set_error('type of main id must be String or ObjectId', code="53701107")
    conn = self.db[table_name]
    if query_dict:
        if record_id is not None:
            query_dict.update({'_id': record_id})
    else:
        if record_id is not None:
            query_dict = {'_id': record_id}
        else:
            self.bind_object.set_error('query dict must be provided while record id is None', code="53701108")
    if insert_dict:
        kwargs.update(insert_dict)
    conn.update(query_dict, {'$set': kwargs}, upsert=True)



class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.ref_rna_v2.refrna_test_api import RefrnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'add_annotation_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'ref_rna_v2.refrna_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = RefrnaTestApiWorkflow(wheet)
        wf.sheet.id = 'ref_rna_v3'
        wf.sheet.project_sn = 'ref_rna_v3'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('ref_rna_v3.add_annotation')
        wf.test_api.ghost_add_annotation(
            # annotation_file='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/add_annotation/annotation',
            task_id="tsg_37259",
            # change_database="custom",
            # annot_type='replace',
            # annot_method='all'
        )

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)