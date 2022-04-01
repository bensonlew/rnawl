# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson import SON
import pandas as pd
from bson.objectid import ObjectId
from types import StringTypes


class ArghubAnalysis(Base):
    def __init__(self, bind_object=None):
        super(ArghubAnalysis, self).__init__(bind_object)
        self._project_type = 'arghub'

    def add_database_main(self, data):
        return self._run_add('arghub', data, main=True)

    def add_database_file(self, filepath, main_id, table_name):
        main_id = self._check_id(main_id)
        df = pd.read_csv(filepath, header=0, index_col=False, sep='\t')
        df['hub_id'] = main_id
        data = self.df_to_mongo(df)
        self._run_add(table_name, data)

    def update_main(self, main_id=None):
        data = {}
        data = self.main_data(data)
        self.bind_object.logger.info("main_id: {},update_main: {}".format(main_id, data))
        self._run_add('analysis', data, main_id=ObjectId(main_id))

    def add_result(self, filepath, main_id, colloction_name):
        main_id = self._check_id(main_id)
        df = pd.read_csv(filepath, header=0, index_col=False, sep='\t')
        df['analysis_id'] = main_id
        data = self.df_to_mongo(df)
        self._run_add(colloction_name, data)

    def df_to_mongo(self, df):
        keys = map(lambda x: x.lower(), df.columns)
        mongo_data = []
        df.apply(lambda x: mongo_data.append(SON(dict(zip(keys, x)))), axis=1)
        return mongo_data

    def _check_id(self, objectid):
        if not isinstance(objectid, ObjectId):
            if isinstance(objectid, StringTypes):
                objectid = ObjectId(objectid)
            else:
                self.bind_object.set_error(
                    "{}必须为ObjectID对象或其字符串形式".format(objectid))
        return objectid

    def _run_add(self, name, data, main_id=None):
        try:
            collection = self.db[name]
            if main_id:
                collection.update_one(
                    {'_id': main_id},
                    {'$set': data},
                )
            else:
                collection.insert_many(data)
        except Exception, e:
            self.bind_object.set_error('导入表%s出错：%s' % (name, e))
            #print('导入表%s出错：%s' % (name, e))
        else:
            self.bind_object.logger.info('导入表%s成功！' % name)
            print('导入表%s成功！' % name)

    def main_data(self, data):
        data['anti_column_data'] = {
            "name": "anti_type"
        }
        data['antibiotic_column_data'] = {
            "name": "anti_class",
        }
        if self.bind_object.option('input_type') == 'read':
            data['abundance_column_data'] = {
                "name": "read_num",
            }
            data['read_table_data'] = {
                "column": [
                    {
                        "field": "arghub_id",
                        "title": "ARG Hub ID",
                        "filter": "true",
                        "sort": "true",
                        "type": "string"
                    },
                    {
                        "field": "gene_name",
                        "title": "Gene Name",
                        "filter": "true",
                        "sort": "true",
                        "type": "string"
                    },
                    {
                        "field": "anti_name",
                        "title": "Anti-bacterial Name",
                        "filter": "true",
                        "sort": "false",
                        "type": "string"
                    },
                    {
                        "field": "anti_type",
                        "title": "Anti-bacterial type",
                        "filter": "true",
                        "sort": "false",
                        "type": "string"
                    },
                    {
                        "field": "resistance_mechanism",
                        "title": "Resistance Mechanism",
                        "filter": "true",
                        "sort": "false",
                        "type": "string"
                    },
                    {
                        "field": "variant_type",
                        "title": "Varient Type",
                        "filter": "true",
                        "sort": "false",
                        "type": "string"
                    },
                    {
                        "field": "variant_info",
                        "title": "Varient Info",
                        "filter": "true",
                        "sort": "false",
                        "type": "string"
                    },
                    {
                        "field": "reads_num",
                        "title": "Reads Number",
                        "filter": "true",
                        "sort": "false",
                        "type": "int"
                    },
                    {
                        "field": "gene_depth",
                        "title": "Gene Depth",
                        "filter": "true",
                        "sort": "false",
                        "type": "float"
                    },
                ],
                "condition": {}
            }
        else:
            data['gene_table_data'] = {
                "column": [
                    {
                        "field": "gene_id",
                        "title": "Gene ID",
                        "filter": "true",
                        "sort": "true",
                        "type": "string"
                    },
                    {
                        "field": "arghub_id",
                        "title": "ARG Hub ID",
                        "filter": "true",
                        "sort": "true",
                        "type": "string"
                    },
                    {
                        "field": "gene_name",
                        "title": "Gene Name",
                        "filter": "true",
                        "sort": "true",
                        "type": "string"
                    },
                    {
                        "field": "length",
                        "title": "Length (bp)",
                        "filter": "true",
                        "sort": "true",
                        "type": "int"
                    },
                    {
                        "field": "anti_name",
                        "title": "Anti-bacterial Name",
                        "filter": "true",
                        "sort": "false",
                        "type": "string"
                    },
                    {
                        "field": "anti_type",
                        "title": "Anti-bacterial type",
                        "filter": "true",
                        "sort": "false",
                        "type": "string"
                    },
                    {
                        "field": "resistance_mechanism",
                        "title": "Resistance Mechanism",
                        "filter": "true",
                        "sort": "false",
                        "type": "string"
                    },
                    {
                        "field": "variant_type",
                        "title": "Varient Type",
                        "filter": "true",
                        "sort": "false",
                        "type": "string"
                    },
                    {
                        "field": "variant_info",
                        "title": "Varient Info",
                        "filter": "true",
                        "sort": "false",
                        "type": "string"
                    },
                    {
                        "field": "align_len",
                        "title": "Align Length (bp)",
                        "filter": "true",
                        "sort": "false",
                        "type": "int"
                    },
                    {
                        "field": "q_coverage",
                        "title": "Coverage(％)",
                        "filter": "true",
                        "sort": "false",
                        "type": "float"
                    },
                    {
                        "field": "evalue",
                        "title": "Evalue",
                        "filter": "true",
                        "sort": "false",
                        "type": "float"
                    },
                    {
                        "field": "score",
                        "title": "Score",
                        "filter": "true",
                        "sort": "false",
                        "type": "float"
                    },
                ],
                "condition": {}
            }
            if self.bind_object.option('input_type') in ['single', 'meta']:
                data['gene_table_data']["column"].extend([
                    {
                        "field": "location",
                        "title": "Location",
                        "filter": "true",
                        "sort": "true",
                        "type": "string"
                    },
                    {
                        "field": "start",
                        "title": "Start",
                        "filter": "true",
                        "sort": "true",
                        "type": "int"
                    },
                    {
                        "field": "end",
                        "title": "End",
                        "filter": "true",
                        "sort": "true",
                        "type": "float"
                    },
                    {
                        "field": "strand",
                        "title": "Strand",
                        "filter": "true",
                        "sort": "false",
                        "type": "string"
                    },
                ])
            if self.bind_object.option('aligner') != 'hmmscan':
                data['gene_table_data']["column"].extend([
                    {
                        "field": "identity",
                        "title": "Identity (％)",
                        "filter": "true",
                        "sort": "false",
                        "type": "float"
                    },
                    {
                        "field": "evaluation",
                        "title": "evaluation",
                        "filter": "true",
                        "sort": "false",
                        "type": "int"
                    },
                ])
        if self.bind_object.option('input_type') in ['meta', 'single']:
            data['gene_neary_mge'] = {
                "column": [
                    {
                        "field": "mge_name",
                        "title": "MGE Name",
                        "filter": "false",
                        "sort": "false",
                        "type": "string"
                    },
                    {
                        "field": "mge_arglist",
                        "title": "Resistance Gene No.",
                        "filter": "false",
                        "sort": "false",
                        "type": "int"
                    },
                    {
                        "field": "type",
                        "title": "MGE Type",
                        "filter": "false",
                        "sort": "false",
                        "type": "string"
                    },
                    {
                        "field": "location",
                        "title": "Query Location",
                        "filter": "false",
                        "sort": "false",
                        "type": "string"
                    },
                    {
                        "field": "start",
                        "title": "Query Start",
                        "filter": "false",
                        "sort": "true",
                        "type": "int"
                    },
                    {
                        "field": "end",
                        "title": "Query End",
                        "filter": "false",
                        "sort": "false",
                        "type": "int"
                    },
                    {
                        "field": "length",
                        "title": "Query Length （bp）",
                        "filter": "false",
                        "sort": "false",
                        "type": "int"
                    },
                ],
                "condition": {}
            }
            data['gene_neary_ele'] = {
                "column": [{
                    "field": "location",
                    "title": "Query Location",
                    "filter": "true",
                    "sort": "false",
                    "type": "string"
                }, {
                    "field": "melem_name",
                    "title": "Component",
                    "filter": "true",
                    "sort": "true",
                    "type": "string"
                }, {
                    "field": "melem_type",
                    "title": "Product",
                    "filter": "true",
                    "sort": "true",
                    "type": "string"
                }, {
                    "field": "strand",
                    "title": "Strand",
                    "filter": "true",
                    "sort": "false",
                    "type": "string"
                }, {
                    "field": "start",
                    "title": "Start",
                    "filter": "true",
                    "sort": "false",
                    "type": "int"
                }, {
                    "field": "end",
                    "title": "End",
                    "filter": "true",
                    "sort": "false",
                    "type": "int"
                }, {
                    "field": "length",
                    "title": "Length(bp)",
                    "filter": "true",
                    "sort": "false",
                    "type": "int"
                }],
                "condition": {}
            }
            data['genome_mge_data'] = {
                    "name": "mge_name",
                    "id": "location",
                    "start": "start",
                    "end": "end",
                    "y": "location",
                    "group": "type",
                    "condition": {"data_type": "mge"}
                }
            data['genome_melem_data'] = {
                    "name": "melem_name",
                    "id": "location",
                    "start": "start",
                    "end": "end",
                    "y": "location",
                    "group": "melem_type",
                    "condition": {"type": "elem"}
                }
            data['mge_melem_data'] = {
                    "name": "melem_name",
                    "id": "location",
                    "start": "start",
                    "end": "end",
                    "y": "location",
                    "group": "melem_type",
                    "condition": {"type": "elem"}
                }
            data['mge_data'] = {
                "column": [{
                    "field": "mge_name",
                    "title": "MGE Name",
                    "filter": "true",
                    "sort": "false",
                    "type": "string"
                }, {
                    "field": "type",
                    "title": "MGE Type",
                    "filter": "true",
                    "sort": "false",
                    "type": "string"
                }, {
                    "field": "location",
                    "title": "Location",
                    "filter": "true",
                    "sort": "false",
                    "type": "string"
                }, {
                    "field": "start",
                    "title": "Start",
                    "filter": "true",
                    "sort": "false",
                    "type": "int"
                }, {
                    "field": "end",
                    "title": "End",
                    "filter": "true",
                    "sort": "false",
                    "type": "int"
                }, {
                    "field": "strand",
                    "title": "Strand",
                    "filter": "true",
                    "sort": "false",
                    "type": "string"
                }, {
                    "field": "length",
                    "title": "Length （bp）",
                    "filter": "true",
                    "sort": "false",
                    "type": "int"
                }, {
                    "field": "seq",
                    "title": "Sequence",
                    "filter": "true",
                    "sort": "false",
                    "type": "string"
                }, {
                }],
                "condition": {"data_type": "mge"}
            }
            data['mge_anti_data'] = {
                "column": [{
                    "field": "gene_id",
                    "title": "Gene ID",
                    "filter": "true",
                    "sort": "false",
                    "type": "string"
                }, {
                    "field": "arghub_id",
                    "title": "ARG Hub ID",
                    "filter": "true",
                    "sort": "false",
                    "type": "string"
                }, {
                    "field": "gene_name",
                    "title": "Gene Name",
                    "filter": "true",
                    "sort": "false",
                    "type": "string"
                }, {
                    "field": "location",
                    "title": "Location",
                    "filter": "true",
                    "sort": "false",
                    "type": "string"
                }, {
                    "field": "start",
                    "title": "Start",
                    "filter": "true",
                    "sort": "false",
                    "type": "int"
                }, {
                    "field": "end",
                    "title": "End",
                    "filter": "true",
                    "sort": "false",
                    "type": "int"
                }, {
                    "field": "strand",
                    "title": "Strand",
                    "filter": "true",
                    "sort": "false",
                    "type": "string"
                }, {
                    "field": "length",
                    "title": "Length （bp）",
                    "filter": "true",
                    "sort": "false",
                    "type": "int"
                }, {
                    "field": "align_len",
                    "title": "Align Length（bp）",
                    "filter": "true",
                    "sort": "false",
                    "type": "string"
                }],
                "condition": {}
            }
            if self.bind_object.option('aligner') != 'hmmscan':
                data['mge_anti_data']["column"].append(
                    {
                    "field": "identity",
                    "title": "Identity（％）",
                    "filter": "true",
                    "sort": "false",
                    "type": "float"
                })
            data['melem_data'] = {
                "column": [{
                    "field": "location",
                    "title": "Location",
                    "filter": "true",
                    "sort": "false",
                    "type": "string"
                }, {
                    "field": "melem_name",
                    "title": "Component",
                    "filter": "true",
                    "sort": "false",
                    "type": "string"
                }, {
                    "field": "melem_type",
                    "title": "Product",
                    "filter": "true",
                    "sort": "false",
                    "type": "string"
                }, {
                    "field": "start",
                    "title": "Start",
                    "filter": "true",
                    "sort": "false",
                    "type": "int"
                }, {
                    "field": "end",
                    "title": "End",
                    "filter": "true",
                    "sort": "false",
                    "type": "int"
                }, {
                    "field": "strand",
                    "title": "Strand",
                    "filter": "true",
                    "sort": "false",
                    "type": "string"
                }, {
                    "field": "length",
                    "title": "Length （bp）",
                    "filter": "true",
                    "sort": "false",
                    "type": "int"
                }, {
                    "field": "seq",
                    "title": "Sequence",
                    "filter": "true",
                    "sort": "false",
                    "type": "string"
                }],
                "condition": {"type": "melem"}
            }
        return data
