# -*- coding: utf-8 -*-

from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
import re
# from biocluster.config import Config
from bson.son import SON
from mbio.api.database.whole_transcriptome.api_base import ApiBase
from bson.objectid import ObjectId


class Enterotyping(Base):
    def __init__(self, bind_object):
        super(Enterotyping, self).__init__(bind_object)
        self._project_type = "tool_lab"
        self.table_type = {
            "BCA_circle": "BCA_circle",
            "BCA_point": "BCA_point",
            "pcoa_circle": "PCoA_circle",
            "pcoa_point": "PCoA_point"
        }

    @report_check
    def add_enterotype(self, main_id=None, main=False, task_id=None, params=None,
                       cluster_name=None, spe_name=None, name=None, stats=False):
        self._tables = []
        _main_collection = self.db['enterotype']
        if not main_id:
            insert_mongo_json = {
                'project_sn': self.bind_object.sheet.project_sn,
                'task_id': task_id,
                'name': name if name else 'Enterotypes_Origin',
                'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
                'status': 'end',
                'desc': '',
                'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            }
            main_id = _main_collection.insert_one(insert_mongo_json).inserted_id
            _main_collection.update({"_id": main_id}, {"$set": {"main_id": main_id}})
        else:
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
        result = _main_collection.find_one({'_id': main_id})
        if not result:
            self.bind_object.set_error('找不到主表id对应的表')
        _main_collection.update_one({'_id': main_id}, {'$set': {'cluster_name': cluster_name}})
        _main_collection.update_one({'_id': main_id}, {'$set': {'spe_name': spe_name}})
        if stats:
            point_list = ["PCoA_point", "BCA_point"]
        else:
            point_list = ["PCoA_point"]
        point_list_data_json = json.dumps(point_list, sort_keys=True, separators=(',', ':'))
        _main_collection.update_one({'_id': main_id}, {'$set': {"main_id": main_id, "point_list": point_list_data_json}})
        self.bind_object.logger.info("主表导入成功")
        return main_id

    @report_check
    def add_enterotype_detail(self, file_path, table_type, update_id, coll_name='enterotype_detail',
                              main_coll='enterotype', update_column=True, db=None,group_table=None):
        """
        导入所有画图的数据
        """

        sample_group_dict = {}
        with open(group_table, 'r') as grp:
            lines = grp.readlines()
            group_list = []
            for line in lines[1:]:
                line = line.strip().split("\t")
                group_name = line[1]
                if group_name not in group_list:
                    group_list.append(group_name)
                sample_group_dict[line[0]] = group_name
        group_sample_dict = {}  ##获得group和样本list对应关系的一个dict
        for group in group_list:
            new_list = []
            for sample in sample_group_dict.keys():
                if sample_group_dict[sample] in [group]:
                    new_list.append(sample)
            group_sample_dict[group] = new_list
        self._tables.append(table_type)
        if not update_id:
            self.bind_object.set_error('需要提供主表ID')
        if not isinstance(update_id, ObjectId):
            update_id = ObjectId(update_id)
        if not db:
            db = self.db
        collection = db["enterotype_detail"]
        sample_enterotype = {}
        type_list = []
        self.bind_object.logger.info("file_path: {}".format(file_path))
        with open(file_path, 'rb') as f:
            all_lines = f.readlines()
            columns = all_lines[0].rstrip().split('\t')[1:]
            for line in all_lines[1:]:
                values = line.rstrip().split('\t')
                sample = str(values[0])
                enterotype = values[1]
                if enterotype not in type_list:
                    type_list.append(enterotype)
                sample_enterotype[sample] = enterotype
        data_temp = []
        for group in group_list:
            self.bind_object.logger.info("group: {}".format(group))
            new_group_list = group_sample_dict[group]
            for type in type_list:
                first_list = []
                for sample in new_group_list:
                    if sample_enterotype[sample] in [type]:
                        first_list.append(sample)
                number = len(first_list)
                insert_data = {
                    "enterotype_id": ObjectId(update_id),
                    "enterotype": type,
                    "group_name": group,
                    "num": number,
                    "group": ",".join(first_list)
                }
                data_temp.append(insert_data)
        try:
            collection.insert_many(data_temp)
            self.bind_object.logger.info("enterotype_detail导入数据成功!")
        except:
            self.bind_object.logger.error("enterotype_detail导入数据失败!")
            self.bind_object.set_error("enterotype_detail导入数据失败")
        try:
            main_collection = self.db["enterotype"]
            summary_table_data = {"table_data": ["enterotype"] + group_list}
            summary_table_data_json = json.dumps(summary_table_data, sort_keys=True, separators=(',', ':'))
            main_collection.update_one({'_id': update_id}, {'$set': {"main_id": ObjectId(update_id),
                                                                     "summary": group_list,
                                                                     "summary_table_data": summary_table_data_json}})
        except:
            self.bind_object.set_error("enterotype更新数据失败")

    @report_check
    def add_enterotype_bar(self, file_path, table_type, update_id, coll_name='enterotype_bar',
                           main_coll='enterotype', update_column=True, db=None,top1_name=None):
        """
        导入柱形图的数据表，类型分为ch和enterotype
        """
        self._tables.append(table_type)
        if not update_id:
            self.bind_object.set_error('需要提供主表ID')
        if not isinstance(update_id, ObjectId):
            update_id = ObjectId(update_id)
        if not db:
            db = self.db
        collection = db[coll_name]
        main_collection = self.db[main_coll]
        data_temp = []
        with open(file_path, 'rb') as f:
            all_lines = f.readlines()
            columns = all_lines[0].rstrip().split('\t')[1:]
            for line in all_lines[1:]:
                values = line.rstrip().split('\t')
                name = values[0]
                if table_type in ['ch']:
                    insert_data = {
                        'enterotype_id': update_id,
                        'type': table_type,
                        'cluster': str(name),
                        "ch_index": float(values[1])
                    }
                    data_temp.append(insert_data)
                else:
                    """
                    insert_data = {
                        'enterotype_id': update_id,
                        'type': table_type,
                        'enterotype': name,
                    }
                    values_dict = dict(zip(columns, map(lambda x: round(float(x), 4), values[1:])))
                    data_temp.append(dict(insert_data, **values_dict))
                    """
                    for x in range(len(columns)):
                        insert_data = {
                            'enterotype_id': update_id,
                            'type': table_type,
                            'enterotype': "type" + name,
                            'name': columns[x],
                            'value': values[x+1],
                            "top1":eval(top1_name)["type" + name]  # + "_" + name
                        }
                        data_temp.append(insert_data)
        if data_temp:
            collection.insert_many(data_temp)
            self.bind_object.logger.info("enterotype_bar导入数据成功!")
        else:
            return None
        if update_column:
            default_column = ['ch', 'enterotype']
            if table_type in ['ch']:
                ch_column_data = {"name": "cluster", "data": "ch_index", "category": "category", "condition": {"type": "ch"}}
                ch_column_data_json = json.dumps(ch_column_data, sort_keys=True, separators=(',', ':'))
                main_collection.update_one({'_id': update_id},
                                           {'$set': {"main_id": update_id,
                                                     "ch_column_data": ch_column_data_json}})
            elif table_type in ['enterotype']:
                enterotype_column_data = {"name": "name", "data": "value", "category": "enterotype","condition": {"type": "enterotype"}}
                enterotype_column_data_json = json.dumps(enterotype_column_data, sort_keys=True, separators=(',', ':'))
                main_collection.update_one({'_id': update_id},
                                           {'$set': {"main_id": update_id,
                                                     "enterotype_column_data": enterotype_column_data_json}})
            else:
                self.bind_object.set_error('错误的表格类型：%s不能在主表中插入相应表头', variables=(table_type))

    @report_check
    def add_enterotype_scatter(self, file_path, table_type, update_id, coll_name='enterotype_scatter',
                               main_coll='enterotype', update_column=True, db=None,group_data=None,type_data=None,top1_name=None):
        """
        导入散点图的数据表
        """
        self._tables.append(table_type)
        self.bind_object.logger.info("table_type111111:{}".format(table_type))
        if not update_id:
            self.bind_object.set_error('需要提供主表ID')
        if not isinstance(update_id, ObjectId):
            update_id = ObjectId(update_id)
        if not db:
            db = self.db
        collection = db[coll_name]
        main_collection = self.db[main_coll]

        group_dict = {}
        type_dict = {}
        with open(group_data) as v,open(type_data) as g:
            data1 = v.readlines()
            data2 = g.readlines()
            for i in data1[1:]:
                group_dict[i.strip().split("\t")[0]] = i.strip().split("\t")[1]
            for x in data2[1:]:
                type_dict[x.strip().split("\t")[0]] = "type" + x.strip().split("\t")[1]
        with open(file_path, 'rb') as f:
            all_lines = f.readlines()
            columns = all_lines[0].rstrip().split('\t')[1:]
            data_temp = []
            for line in all_lines[1:]:
                values = line.rstrip().split('\t')
                name = values[0]
                insert_data = {
                    'enterotype_id': update_id,
                    'name': name
                }
                if table_type in ['pcoa_point']:
                    insert_data['category'] = group_dict[name]
                    insert_data['group'] = type_dict[name]
                    insert_data['type'] = "scatter"
                    insert_data['subtype'] = "PCoA_point"
                    values_dict = dict(zip(columns, values[1:]))
                    """
                    # arrow
                    insert_data2 = {
                        'enterotype_id': update_id,
                        'type': "pcoa_arrow",
                        'name': name,
                        'x': "x",
                        'y': "y",
                        'x1': 0,
                        'y1': 0,
                        'x2': float(values[1]),
                        'y2': float(values[2])
                    }
                    data_temp.append(insert_data2)
                    """
                    data_temp.append(dict(insert_data, **values_dict))
                elif table_type in ['pcoa_circle']:
                    for x in range(len(columns)):
                        insert_data = {
                            'enterotype_id': update_id,
                            'subtype': 'PCoA_circle',
                            'type': "ellipse",
                            'name': eval(top1_name)["type" + str(columns[x])], # +"_"+str(columns[x]),
                            'group': "type" + str(columns[x]),
                            'data':  values[x+1],
                            'x':'x',
                            'y':'y',
                            'method':'trend',
                            'method_type':'method_type'
                        }
                        data_temp.append(insert_data)
                elif table_type in ['BCA_circle']:
                    for x in range(len(columns)):
                        insert_data = {
                            'enterotype_id': update_id,
                            'subtype': 'BCA_circle',
                            'type': "ellipse",
                            'name': eval(top1_name)["type" + str(columns[x])]+"_"+str(columns[x]),
                            'group': "type" + str(columns[x]),
                            'data': values[x + 1],
                            'x': 'x',
                            'y': 'y',
                            'method': 'trend',
                            'method_type': 'method_type'
                        }
                        data_temp.append(insert_data)
                else:  ##BCA_point
                    if len(values) != 2:
                        insert_data['category'] = group_dict[name]
                        insert_data['group'] = type_dict[name]
                        insert_data['type'] = "scatter"
                        insert_data['subtype'] = "BCA_point"
                        values_dict = dict(zip(columns, map(lambda x: round(float(x), 4), values[1:])))
                        data_temp.append(dict(insert_data, **values_dict))
            if data_temp:
                collection.insert_many(data_temp)
                self.bind_object.logger.info("enterotype_scatter导入数据成功!")
            else:
                return None
        if update_column:
            default_column = ['BCA_circle', 'BCA_point', 'pcoa_circle', 'pcoa_point']
            self.bind_object.logger.info("table_type:{}".format(table_type))
            if table_type in ['pcoa_point']:
                pcoa_scatter_data = {"name": "name", "group":"group", "category": "category", "data": ["x","y"], "condition": {"type": "scatter"}}
                pcoa_scatter_data_json = json.dumps(pcoa_scatter_data, sort_keys=True, separators=(',', ':'))
                main_collection.update_one({'_id': update_id},
                                           {'$set': {"main_id": update_id,
                                                     "pcoa_scatter_data": pcoa_scatter_data_json}})
                #arrow_data = {"name": "name", "condition": {"type": "pcoa_arrow"}}
                #arrow_data_json = json.dumps(arrow_data, sort_keys=True, separators=(',', ':'))
                #main_collection.update_one({'_id': update_id},{'$set': {"main_id": update_id,"pcoa_arrow_data": arrow_data_json}})
            elif table_type in ['pcoa_circle']:
                pcoa_ellipse_data = {"data": "data", "group": "group", "name": "name", "condition": {"type": "ellipse"}}
                pcoa_ellipse_data_json = json.dumps(pcoa_ellipse_data, sort_keys=True, separators=(',', ':'))
                main_collection.update_one({'_id': update_id},
                                           {'$set': {"main_id": update_id,
                                                     "pcoa_ellipse_data": pcoa_ellipse_data_json}})

            elif table_type in ['BCA_point']:
                pass
            #    if data_temp:
            #        bca_scatter_data = {"name": "name", "category": "category", "data": ["x","y"], "condition": {"type": "BCA_point"}}
            #        bca_scatter_data_json = json.dumps(bca_scatter_data, sort_keys=True, separators=(',', ':'))
            #        main_collection.update_one({'_id': update_id},
            #                                    {'$set': {"main_id": update_id,
            #                                            "bca_scatter_data": bca_scatter_data_json}})
            elif table_type in ['BCA_circle']:
                pass
            #    bca_ellipse_data = {"data": "data", "group": "group", "name": "name", "condition": {"type": "BCA_circle"}}
            #    bca_ellipse_data_json = json.dumps(bca_ellipse_data, sort_keys=True, separators=(',', ':'))
            #    main_collection.update_one({'_id': update_id},
            #                               {'$set': {"main_id": update_id,
            #                                         "bca_ellipse_data": bca_ellipse_data_json}})
            else:
                self.bind_object.set_error('错误的表格类型：%s不能在主表中插入相应表头', variables=(table_type))

    @report_check
    def add_enterotype_detail_cluster(self, id=None, file_path=None, name=None, update_column=True, update_type=None):
        insert_data = list()
        main_collection = self.db['enterotype']
        with open(file_path, 'rb') as r:
            all_lines = r.readlines()
            head = all_lines[0].strip()
            head = re.split('\t', head)
            sample = head[:-2]
            for line in all_lines[1:101]:
                line = line.strip()
                line = re.split('\t', line)
                otu_detail = dict()
                otu_detail['species_name'] = line[0].split("; ")[-1]  # 保留最低层级的物种名称
                otu_detail['enterotype_id'] = ObjectId(id)
                otu_detail['type'] = str(name)
                otu_detail["sum"] = line[-2]
                otu_detail["percent"] = line[-1]
                insert_data.append(otu_detail)
        try:
            collection = self.db['enterotype_cluster']
            collection.insert_many(insert_data)
        except Exception as e:
            self.bind_object.logger.error("导入enterotype_cluster表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入enterotype_cluster表格成功")
        if update_column:
            settled_params = {"software": "R-3.3.1 (cluster,clusterSim)"}
            settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
            if update_type:
                type_table_data = {"table_data": ["species_name", "sum", "percent"],
                                   "condition": {"type": update_type.split(",")}}
                table_data_json = json.dumps(type_table_data, sort_keys=True, separators=(',', ':'))
                main_collection.update_one({'_id': ObjectId(id)}, {'$set': {"settled_params": settled_params_json,
                                                                            "type_table_data": table_data_json}})
