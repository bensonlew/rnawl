# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
#last modify:qc.zhang-20181107:当组成分析为bubble或者heatmap时增加composition_detail表字段level_list;当组成分析类型为bar或者bubble时增加主表字段spe_sort用以区分新老项目和排序；


from biocluster.api.database.base import Base, report_check
import re
import json
import datetime
from bson.objectid import ObjectId
from mbio.packages.metagenomic.id_convert import name2id


class Composition(Base):
    def __init__(self, bind_object):
        super(Composition, self).__init__(bind_object)
        self._project_type = "metagenomic"
        self.sample_2_id = ""

    @report_check
    def add_composition(self, graphic_type, anno_type, species_tree=None, specimen_tree=None, name=None, params=None,
                        specimen_list=None, species_list=None):
        task_id = self.bind_object.sheet.id
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": task_id,
            #"geneset_id": geneset_id,
            "anno_type": anno_type,
            #"anno_id": anno_id,
            "name": name,
            #"level_id": level_id,
            "status": "end",
            "specimen_tree": specimen_tree,
            "species_tree": species_tree,
            "species_list": species_list,
            "specimen_list": specimen_list,
            #"group_id": group_id,
            "graphic_type": graphic_type,
            "desc": "组成分析",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["composition"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_composition_detail(self, file_path, composition_id, species_tree=None, specimen_tree=None,
                               group_method=None, level_color=None, main_col="composition"):
        task_collection = self.db[main_col]
        compositioninfo= task_collection.find_one({'_id': ObjectId(composition_id)})
        task_name2id = compositioninfo["task_id"]
        graphic_type = compositioninfo["graphic_type"]
        anno_type = compositioninfo["anno_type"]
        self.sample_2_id = name2id(task_name2id, type="task")
        self.bind_object.logger.info(self.sample_2_id)
        self.bind_object.logger.info("成功进入composition_detail函数")
        species_list = ""
        specimens_list = [] ##add by qingchen.zhang
        new_sample_list = ''
        if specimen_tree != "":
            with open(specimen_tree, 'r') as f:
                specimen_tree = f.readline()
                speciment = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', specimen_tree)  #zouguanqing 20181228
                #new_sample_list =','.join([i[1] for i in speciment])   #zouguanqing 20181228
                for i in range(len(speciment)):
                    sp = speciment[i][1]
                    if sp in self.sample_2_id.keys() and group_method == "":
                        sp_id = self.sample_2_id[sp]
                        specimen_tree = specimen_tree.replace("(" + sp + ":", "(" + str(sp_id) + ":")
                        specimen_tree = specimen_tree.replace("," + sp + ":", "," + str(sp_id) + ":")
                for i in range(len(speciment)):
                    e = speciment[i][1]
                    if e in self.sample_2_id.keys() and group_method == "":
                        speciem = self.sample_2_id[e]
                        if specimens_list and str(speciem) not in specimens_list:
                           specimens_list.append(str(speciem))##add by qingchen.zhang
                        elif specimens_list == []:##add by qingchen.zhang
                             specimens_list.append(str(speciem))##add by qingchen.zhang
                    else:
                        if specimens_list and e not in specimens_list:
                            specimens_list.append(str(e))##add by qingchen.zhang
                        elif specimens_list == []:##add by qingchen.zhang
                            specimens_list.append(str(e))##add by qingchen.zhang
        if species_tree != "":
            with open(species_tree, 'r') as f:
                species_tree = f.readline()
        data_list = []
        with open(file_path, "r") as f:
            self.bind_object.logger.info("开始导丰度表")
            if graphic_type=="bubble" or graphic_type=="heatmap":
                if level_color:
                    specimens = f.readline().strip().split("\t")[1:-1]
                else:
                    specimens = f.readline().strip().split("\t")[1:]
            else:
                specimens = f.readline().strip().split("\t")[1:]
            for line in f:
                line = line.strip().split("\t")
                if graphic_type=="bubble" :
                    if level_color:
                        if line[0] == "others":
                            data = {
                            "composition_id": ObjectId(composition_id),
                            "species_name": line[0],
                            "level_color": "others",
                            }
                        else:
                            data = {
                            "composition_id": ObjectId(composition_id),
                            "species_name": line[0],
                            "level_color": line[-1],
                            }
                    else:
                        data = {
                            "composition_id": ObjectId(composition_id),
                            "species_name": line[0],
                        }
                elif graphic_type=="heatmap":
                    if level_color:
                        data = {
                        "composition_id": ObjectId(composition_id),
                        "species_name": line[0],
                        "level_color": line[-1],
                        }
                    else:
                        data = {
                            "composition_id": ObjectId(composition_id),
                            "species_name": line[0],
                        }
                else:
                    data = {
                        "composition_id": ObjectId(composition_id),
                        "species_name": line[0],
                    }
                species_list = species_list + ";" + line[0]
                for n, e in enumerate(specimens):
                    if e in self.sample_2_id.keys() and group_method == "":
                        speciem = self.sample_2_id[e]
                        data[str(speciem)] = line[n + 1]
                        if specimen_tree == "":
                            if specimens_list and str(speciem) not in specimens_list:
                                specimens_list.append(str(speciem))##add by qingchen.zhang
                            elif specimens_list == []:##add by qingchen.zhang
                                specimens_list.append(str(speciem))##add by qingchen.zhang
                    else:
                        data[e] = line[n + 1]
                        if specimen_tree == "":
                            if specimens_list and e not in specimens_list:
                                specimens_list.append(str(e))##add by qingchen.zhang
                            elif specimens_list == []:##add by qingchen.zhang
                                specimens_list.append(str(e))##add by qingchen.zhang
                data_list.append(data)
            species_list = species_list.lstrip(';')
        try:
            collection = self.db[main_col + "_detail"]
            collection.insert_many(data_list)
            self.bind_object.logger.info("开始刷新主表写树")
            main_collection = self.db[main_col]
            if graphic_type=="heatmap":
                main_collection.update({"_id": ObjectId(composition_id)}, {"$set": {"specimen_tree": specimen_tree,
                                                                            "species_tree": species_tree,
                                                                            "specimen_list": ','.join(specimens_list),
                                                                            "species_list": species_list,
                                                                            "version" :1}})
            else:
                main_collection.update({"_id": ObjectId(composition_id)}, {"$set": {"specimen_tree": specimen_tree,
                                                                            "species_tree": species_tree,
                                                                            "specimen_list": ','.join(specimens_list),
                                                                            "species_list": species_list}})
        except Exception, e:
            self.bind_object.logger.info("导入丰度表格出错:%s" % e)
            self.bind_object.set_error("导入丰度表格出错", code="52800401")
        else:
            self.bind_object.logger.info("导入丰度表格成功")

