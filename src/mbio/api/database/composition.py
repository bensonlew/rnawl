# -*- coding: utf-8 -*-
# __author__ = 'xuting'

from biocluster.api.database.base import Base, report_check
import re
import datetime
from bson.objectid import ObjectId
from types import StringTypes
# from biocluster.config import Config


class Composition(Base):

    def __init__(self, bind_object):
        super(Composition, self).__init__(bind_object)
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB
        self.task_id = ""
        self.name_id = dict()

    @report_check
    def add_sg_otu(self, params, from_otu_table, name=None, newick_id=None):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!", code="51000501")
        collection = self.db["sg_composition"]
        result = collection.find_one({"_id": from_otu_table})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}在sg_otu表里找到相应的记录".format(str(from_otu_table)))
            self.bind_object.set_error("在sg_otu表中找不到相应记录", code="51000502")
        project_sn = result['project_sn']
        self.task_id = result['task_id']
        if not name:
            name = "community_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        insert_data = {
            "project_sn": project_sn,
            'task_id': self.task_id,
            'from_id': str(from_otu_table),
            'name': self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else name,
            "params": params,
            "newick_id": newick_id,
            'status': 'end',
            'desc': 'otu table after Cluster Analysis',
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "show": 0,
            "type": "circos"
        }
        collection = self.db["sg_composition"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        #collection.update({"_id": inserted_id}, {"$set": {"main_id": inserted_id}})
        return inserted_id

    @report_check
    def add_sg_composition(self, params, from_otu_table, name=None, newick_id=None):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!", code="51000501")
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": from_otu_table})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}在sg_otu表里找到相应的记录".format(str(from_otu_table)))
            self.bind_object.set_error("在sg_otu表中找不到相应记录", code="51000502")
        project_sn = result['project_sn']
        self.task_id = result['task_id']
        if not name:
            name = "community_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        insert_data = {
            "project_sn": project_sn,
            'task_id': self.task_id,
            'from_id': str(from_otu_table),
            'otu_id': from_otu_table,
            'name': self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else name,
            "params": params,
            "newick_id": newick_id,
            'status': 'end',
            'desc': 'otu table after Cluster Analysis',
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "show": 0,
            "type": "circos"
        }
        collection = self.db["sg_composition"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        # collection.update({"_id": inserted_id}, {"$set": {"main_id": inserted_id}})
        return inserted_id

    @report_check
    def add_sg_otu_detail(self, file_path, new_otu_id, from_otu_id, type):
        if from_otu_id != 0 and not isinstance(from_otu_id, ObjectId):
            if isinstance(from_otu_id, StringTypes):
                from_otu_id = ObjectId(from_otu_id)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!", code="51000503")
        self.bind_object.logger.info("开始导入sg_composition_detail表")
        #self._get_name_id(from_otu_id)
        find_otu = self.db['sg_composition'].find_one({"_id": ObjectId(new_otu_id)})
        if find_otu:
            self.task_id = find_otu['task_id']
        else:
            self.bind_object.set_error("COMP_ID没有找到相关的主表信息", code="51000504")
        insert_data = list()
        spe_str=""  #guanqing.zou 物种按丰度排序，以|分割的字符串 20180411
        
        ##guanqing.zou 20180508 begin
        with open(file_path, 'rb') as r:
            lines = r.readlines()
            head = lines[0].strip('\r\n')
            head = re.split('\t', head)
            new_head = head[1:]
            if type in ["ternary","Ternary"] :
                oth = lines[1].split('\t')
                oth_com = re.split(r"\s*;\s*", oth[0])
                othlen = len(oth_com)

            for line in lines[1:]:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                sample_num = line[1:]
                classify_list = re.split(r"\s*;\s*", line[0])
                spe_str+=classify_list[-1]+'|'       #guanqing 20180411
                    
                otu_detail = dict()
                otu_detail['comp_id'] = ObjectId(new_otu_id)

                ##guanqing 20180423
                #for cf in classify_list:
                if type in ["circos","Circos"] :
                    cf = classify_list[-1]
                elif type in ["ternary","Ternary"] :
                    if classify_list[-1] not in ['other','others'] :
                        for ele in range(othlen):
                            if classify_list[ele] != oth_com[ele] :
                                oth_com[ele] = 'others'
                        cf = ';'.join([aa.strip() for aa in classify_list])
                    else:
                        cf = ';'.join([aa.strip() for aa in oth_com])    
                #otu_detail[cf[0:3].lower()] = cf
                otu_detail["species_name"]= cf
        ###guanqing.zou 20180508 end
                for i in range(0, len(sample_num)):
                    otu_detail[new_head[i]] = sample_num[i]
                otu_detail['task_id'] = self.task_id
                insert_data.append(otu_detail)
        try:
            collection = self.db['sg_composition_detail']
            collection.insert_many(insert_data)
            main_collection = self.db['sg_composition']  #guanqing.zou
            main_collection.update({"_id":ObjectId(new_otu_id)},{"$set":{"spe_sort":spe_str[0:-1]}})  #guanqing 20180411
            main_collection.update({"_id":ObjectId(new_otu_id)},{"$set":{"specimen_list":new_head}})
            #main_collection.update({"_id": ObjectId(new_otu_id)}, {"$set": {"main_id": ObjectId(new_otu_id)}})
        except Exception as e:
            self.bind_object.logger.error("导入sg_composition_detail表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入sg_composition_detail表格成功")

