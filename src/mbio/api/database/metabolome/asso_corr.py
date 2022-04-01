# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:20180718
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
# from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from mbio.packages.metabolome.common import check_metab


class AssoCorr(Base):
    def __init__(self, bind_object):
        super(AssoCorr, self).__init__(bind_object)
        self._project_type = "metabolome"

    @report_check
    def add_association_corr(self,  name=None, main_id=None, params =None, assoc_tree=None,metab_tree=None,listfile=None):
        if not main_id:
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else 'AssoCorr_Origin',
                'params': params if params else '',
                'status': 'end',
                'assoc_list': '',
                'assoc_tree': '',
                'metab_tree': '',
                'main_id': ''
            }
            try:
                collection = self.db['association_corr']
                association_corr_id = collection.insert_one(insert_data).inserted_id
                self.update_table("main_id", association_corr_id, association_corr_id)
            except Exception, e:
                self.bind_object.set_error('导入association_corr主表异常:%s', variables=(e), code="54700101")
        else:
            self.update_table("main_id", main_id, main_id)
            association_corr_id = main_id
        if assoc_tree:
            if not os.path.exists(assoc_tree):
                self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(assoc_tree), code="54700102")
            with open(assoc_tree,"r") as f:
                assoc_tree = f.readline().strip()
                asso_list = f.readline().strip().split(";")
            self.update_table("assoc_tree", assoc_tree, association_corr_id)
            self.update_table("assoc_list", asso_list, association_corr_id)
        if metab_tree:
            if not os.path.exists(metab_tree):
                self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(metab_tree), code="54700103")
            with open(metab_tree,"r") as f2:
                metabtree = f2.readline().strip()
                metabtree = check_metab(metabtree)
            self.update_table("metab_tree", metabtree, association_corr_id)
        if listfile:
            if not os.path.exists(listfile):
                self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(listfile), code="54700104")
            with open(listfile,"r") as f2:
                asso_list = f2.readline().strip().split("\t")
                asso_list = asso_list[1:len(asso_list)]
            self.update_table("assoc_list", asso_list, association_corr_id)
        return association_corr_id

    @report_check
    def add_association_corr_detail(self, association_corr_id, corr_file,p_file):
        metabset_corr_id = self.check_id(association_corr_id,"association_corr_id")
        if not os.path.exists(corr_file):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(corr_file), code="54700105")
        data_list = []
        result = self.db['association_corr'].find_one({'_id': association_corr_id})
        if not result:
            self.bind_object.set_error('找不到%s对应的主表id', variables=(association_corr_id), code="54700106")
        else:
            task_id = result['task_id']
            #samples_dic = name2id(task_id, type="task")
        with open(corr_file, 'rb') as f:
            head = f.next()
            sams = head.strip().split("\t")
            for line in f:
                line = line.strip().split('\t')
                name = check_metab(line[0])
                insert_data = {
                    'corr_id': association_corr_id,
                    'metab':name,
                    'type': 'corr'
                }
                for i in range(1, len(sams)):
                    sam_corr = float(line[i])
                    insert_data[sams[i]] = sam_corr
                #data_list.append(insert_data)
                try:
                    collection = self.db['association_corr_detail']
                    #collection.insert_many(data_list)
                    collection.insert(insert_data, check_keys=False)
                except Exception as e:
                    self.bind_object.set_error("导入表格%s信息出错:%s" , variables=(corr_file, e), code="54700107")
                else:
                    pass
        self.bind_object.logger.info("导入表格%s信息成功!" % corr_file)
        data_list = []
        with open(p_file, 'rb') as f:
            head = f.next()
            sams = head.strip().split("\t")
            for line in f:
                line = line.strip().split('\t')
                name = check_metab(line[0])
                insert_data = {
                    'corr_id': metabset_corr_id,
                    'metab':name,
                    'type': 'pvalue'
                }
                for i in range(1, len(sams)):
                    sam_corr = float(line[i])
                    insert_data[sams[i]] = sam_corr
                #data_list.append(insert_data)
                try:
                    collection = self.db['association_corr_detail']
                    #collection.insert_many(data_list)
                    collection.insert(insert_data, check_keys=False)
                except Exception as e:
                    self.bind_object.set_error("导入表格%s信息出错:%s" , variables=(p_file, e), code="54700108")
                else:
                    pass
        self.bind_object.logger.info("导入表格%s信息成功!" % p_file)


    @report_check
    def update_table(self, str, name, main_table_id):
        try:
            self.db['association_corr'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {str: name}})
        except Exception as e:
            self.bind_object.set_error('association_corr%s字段出错:%s', variables=(str,e), code="54700109")

    @report_check
    def check_id(self, object_id, type):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('%s必须为ObjectId对象或其对应的字符串！', variables=(type), code="54700110")
        return object_id
