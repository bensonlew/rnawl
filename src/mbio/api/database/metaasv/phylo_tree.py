# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'


from biocluster.api.database.base import Base, report_check
import os
from collections import defaultdict
import json
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import re



class PhyloTree(Base):
    def __init__(self, bind_object):
        super(PhyloTree, self).__init__(bind_object)
        self._project_type = 'metaasv'

    @report_check
    def add_phylo_tree_info_for_meta(self, asv_id, params=None):
        collection = self.db["phylo_tree"]
        # params = {}
        main_data = [('project_sn', self.bind_object.sheet.project_sn),
                     ('task_id', self.bind_object.sheet.id),
                     ('asv_id', asv_id),
                     ('name', 'Tree_Origin'),
                     ('status', 'end'),
                     ('params', ""),
                     ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                     ('newicktree', open(self.bind_object.output_dir + '/ASV/ASV_phylo.tre').read())]
        self.main_id = collection.insert_one(SON(main_data)).inserted_id
        try:
            # if self.bind_object.sheet.option('color_level_id'):
                # self._add_species_categories()
            # self.specimen_categories = self._add_format_otu()
            collection.update_one({'_id': self.main_id}, {'$set': {'main_id': self.main_id}})
        except Exception as e:
            self.bind_object.logger.error("Phylo tree导入数据失败: %s" % e)
            self.bind_object.set_error("Phylo tree导入数据失败")
        self.bind_object.logger.info('Phylo tree导入数据库完成。')
        return str(self.main_id)

    @report_check
    def add_phylo_tree_for_report(self, otu_id,task_id,tree_type,tree_file):
        collection = self.db["sg_phylo_tree"]
        main_data = [('project_sn', self.bind_object.sheet.project_sn),
                     ('task_id',task_id),
                     ('table_id', ObjectId(otu_id)),
                     ('name', 'Tree_Origin'),
                     ('status', 'end'),
                     ('params', 'null'),
                     ('tree_type', tree_type ),
                     ('level_id', 9),
                     ("table_type" ,"otu"),
                     ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                     ('value', open(tree_file).read())]
        self.main_id = collection.insert_one(SON(main_data)).inserted_id
        self.bind_object.logger.info('Phylo tree导入数据库完成。')
        return str(self.main_id)

    # 专门用于即时计算导表的模块，不可放在metabase中使用。对应的worflow为metabase.report.plot_tree
    @report_check
    def add_phylo_tree_info(self, main_id,seq_num=None):
        self.main_id = ObjectId(main_id)
        tree = open(self.bind_object.output_dir + '/phylo_tree.tre')
        tree_str = tree.read()
        collection = self.db['phylo_tree_detail']
        insert_data = {
            "name": "",
            "phylo_tree_id" : self.main_id,
            "data": tree_str
        }
        self.detail_id = collection.insert_one(insert_data).inserted_id
        self.specimen_categories = self._add_format_otu()
        tree_str_list = re.sub('([\(,])([dkpcofgs])_[^_]',r'\1\2__',tree_str)
        main_collection = self.db['phylo_tree']
        if self.bind_object.sheet.option('color_level_id'):
            self._add_species_categories()
        tree_data = {
                "tree_data": {"name":"name",
                        "categories": ""}
                }
        tree_data_json = json.dumps(tree_data, sort_keys=True, separators=(',', ':'))
        if seq_num:
            main_collection.update_one({'_id': self.main_id},
                                   {'$set': {
                                       "main_id": self.main_id,
                                       "tree_data": tree_data_json,
                                       'newicktree': tree_str_list,'seq_num': seq_num,
                                       'specimen_categories': self.specimen_categories}})
        else:
            main_collection.update_one({'_id': self.main_id},
                                   {'$set': {
                                       "main_id": self.main_id,
                                       "tree_data": tree_data_json,
                                       'newicktree': tree_str_list,
                                       'specimen_categories': self.specimen_categories}})
        tree.close()
        self.bind_object.logger.info('Phylo tree导入数据库完成。')
        return str(self.main_id)

    def add_phylo_tree_info_2(self, main_id,tree_file, seq_num=None, software=None):
        main_id = ObjectId(main_id)
        tree = open(tree_file)
        main_collection = self.db['personal_phylo_tree']
        collection = self.db['personal_phylo_tree_detail']
        insert_data = {
            "name" : "",
            "phylo_tree_id" : main_id,
            'data': tree.read()
        }
        collection.insert_one(insert_data)
        tree_data = {
                "tree_data": {"name":"name",
                        "categories": ""}
                }
        if software:
            settled_params = {"software" : software}
        else:
            settled_params = {"software" : "IQ-TREE"}
        settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))

        tree_data_json = json.dumps(tree_data, sort_keys=True, separators=(',', ':'))
        if seq_num:
            main_collection.update_one({'_id': main_id},{'$set': {'newicktree': tree.read(),
                                                                  'main_id': main_id,
                                                                  'settled_params': settled_params_json,
                                                                  'tree_data': tree_data_json,
                                                                  'seq_num':seq_num}})
        else:
            main_collection.update_one({'_id': main_id},{'$set': {'newicktree': tree.read(),
                                                                  'settled_params': settled_params_json,
                                                                  'tree_data': tree_data_json,
                                                                  'main_id': main_id}})
        tree.close()
        self.bind_object.logger.info('Phylo tree导入数据库完成!')

    @report_check
    def _add_species_categories(self):
        species_group = self.bind_object.output_dir + '/species_group.xls'
        if not os.path.isfile(species_group):
            self.categories = None
            return self.categories
        collection = self.db['phylo_tree_detail']
        with open(species_group) as f:
            f.readline()
            group = defaultdict(list)
            for i in f:
                line_sp = i.strip().split('\t')
                group[line_sp[1]].append(line_sp[0])
            categories = json.dumps(group.keys())
            species = json.dumps(group.values())
            try:
                collection.update_one({"_id": self.detail_id}, {"$set":{"species": species,"categories": categories}})
            except:
                self.set_error("更新详情表phylo_tree_detail失败!")

    @report_check
    def _add_format_otu(self):
        collection = self.db['phylo_tree_bar']
        with open(self.bind_object.output_dir + '/species_table.xls') as f:
            origin_categories = f.readline().rstrip().split('\t')[1:]
            categories = origin_categories
            categories.insert(0, 'species_name')
            insert_data = []
            for i in f:
                sp_list=i.strip().split('\t')
                one = zip(categories, i.strip().split('\t'))
                one.insert(0, ('phylo_tree_id', self.main_id))
                insert_data.append(SON(one))
        collection.insert_many(insert_data)
        columns_data = {
            "columns_data": {"name":"species_name",
                "data": origin_categories}
            }
        columns_data_json = json.dumps(columns_data, sort_keys=True, separators=(',', ':'))
        try:
            main_collection = self.db['phylo_tree']
            main_collection.update_one({'_id': self.main_id}, {'$set': {'columns_data': columns_data_json}})
        except Exception as e:
            self.bind_object.logger.error("Phylo tree更新失败: %s" % e)
            self.bind_object.set_error("Phylo tree更新失败")
        return categories[1:]

    @report_check
    def add_phylo_tree_info_workflow(self, main_id,tree_path,seq_num=None):
        self.main_id = ObjectId(main_id)
        tree = open(tree_path)
        tree_str = tree.read()
        collection = self.db['phylo_tree_detail']
        insert_data = {
            "name": "",
            "phylo_tree_id" : self.main_id,
            "data": tree_str
        }
        self.detail_id = collection.insert_one(insert_data).inserted_id
        # self.specimen_categories = self._add_format_otu()
        tree_str_list = re.sub('([\(,])([dkpcofgs])_[^_]',r'\1\2__',tree_str)
        main_collection = self.db['phylo_tree']
        # if self.bind_object.sheet.option('color_level_id'):
        #     self._add_species_categories()
        tree_data = {
                "tree_data": {"name":"name",
                        "categories": ""}
                }
        tree_data_json = json.dumps(tree_data, sort_keys=True, separators=(',', ':'))

        if seq_num:
            main_collection.update_one({'_id': self.main_id},
                                   {'$set': {
                                       "main_id": self.main_id,
                                       "tree_data": tree_data_json,
                                       'newicktree': tree_str_list,'seq_num': seq_num,}})
        else:
            main_collection.update_one({'_id': self.main_id},
                                   {'$set': {
                                       "main_id": self.main_id,
                                       "tree_data": tree_data_json,
                                       'newicktree': tree_str_list,}})
        tree.close()
        self.bind_object.logger.info('Phylo tree导入数据库完成。')
        return str(self.main_id)