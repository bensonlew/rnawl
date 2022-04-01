# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:20180523
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
# from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from mbio.packages.metabolome.common import check_metab, check_metab_list


class MetabsetCorr(Base):
    def __init__(self, bind_object):
        super(MetabsetCorr, self).__init__(bind_object)
        self._project_type = "metabolome"

    @report_check
    def add_metabset_corr(self, name=None, main_id=None, params =None, tree_file=None, metab=None, list_file=None,table_type=None,metab_set=None):
        if not main_id:
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else 'MetabCorr_Origin',
                'params': params if params else '',
                'status': 'end',
                'metab': metab,
                'tree': '',
                'main_id': ''
            }
            if table_type:
                insert_data['table_type'] = table_type
            if metab_set:
                insert_data['metab_set'] = metab_set

            try:
                collection = self.db['metabset_corr']
                metabset_corr_id = collection.insert_one(insert_data).inserted_id
                self.update_table("main_id", metabset_corr_id, metabset_corr_id)
            except Exception, e:
                self.bind_object.set_error('导入metabset_corr主表异常:%s', variables=(e), code="54701401")
        else:
            self.update_table("main_id", main_id, main_id)
            metabset_corr_id = main_id
        if tree_file:
            if not os.path.exists(tree_file):
                self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(tree_file), code="54701402")
            with open(tree_file, "r") as f:
                metabtree = f.readline().strip()
                metabtree = check_metab(metabtree)
                metab = f.readline().strip().split(";")
                metab = check_metab_list(metab)
            self.update_table("tree", metabtree, metabset_corr_id)
            self.update_table("metab", metab, metabset_corr_id)
        if list_file:
            if not os.path.exists(list_file):
                self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(list_file), code="54701403")
            with open(list_file, "r") as f2:
                metab = f2.readline().strip().split("\t")
                metab = metab[1:len(metab)]
                metab = check_metab_list(metab)
                if metab[0] == 'Metabolite':
                    metab.pop(0)
            self.update_table("metab", metab, metabset_corr_id)
        return metabset_corr_id

    @report_check
    def add_metabset_corr_detail(self, metabset_corr_id, corr_file, p_file):
        metabset_corr_id = self.check_id(metabset_corr_id, "metabset_corr_id")
        if not os.path.exists(corr_file):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(corr_file), code="54701404")
        data_list = []
        result = self.db['metabset_corr'].find_one({'_id': metabset_corr_id})
        if not result:
            self.bind_object.set_error('找不到%s对应的主表id', variables=(metabset_corr_id), code="54701405")
        else:
            task_id = result['task_id']
            #samples_dic = name2id(task_id, type="task")

        with open(corr_file, 'rb') as f:
            head = f.next()
            sams = head.strip().split("\t")
            if 'ID' in sams:
                mid = 1
            else:
                mid = 0

            for line in f:
                line = line.strip().split('\t')
                name = line[mid]
                insert_data = {
                    'corr_id': metabset_corr_id,
                    'metab': name,
                    'type': 'corr'
                }
                for i in range(mid+1, len(sams)):
                    sam_corr = float(line[i])
                    insert_data[sams[i]] = sam_corr
                    if mid == 1:
                        insert_data['o_id'] = line[0]  #20190617
                #data_list.append(insert_data)
                try:
                    collection = self.db['metabset_corr_detail']
                    #collection.insert_many(data_list)
                    collection.insert(insert_data, check_keys=False)
                except Exception as e:
                    self.bind_object.set_error("导入表格%s信息出错:%s" , variables=(corr_file, e), code="54701406")
                else:
                    pass
        self.bind_object.logger.info("导入表格%s信息成功!" % corr_file)
        data_list = []
        with open(p_file, 'rb') as f:
            head = f.next()
            sams = head.strip().split("\t")
            if 'ID' in sams:
                mid = 1
            else:
                mid = 0

            for line in f:
                line = line.strip().split('\t')
                name = line[mid]
                insert_data = {
                    'corr_id': metabset_corr_id,
                    'metab': name,
                    'type': 'pvalue'
                }
                for i in range(mid+1, len(sams)):
                    try:
                        sam_corr = float(line[i])
                    except ValueError:
                        sam_corr = 1.0
                    except IndexError:
                        sam_corr = 1.0
                    insert_data[sams[i]] = sam_corr
                    if mid == 1:
                        insert_data['o_id'] = line[0]  #20190617
                #data_list.append(insert_data)
                try:
                    collection = self.db['metabset_corr_detail']
                    #collection.insert_many(data_list)
                    collection.insert(insert_data, check_keys=False)
                except Exception as e:
                    self.bind_object.set_error("导入表格%s信息出错:%s" , variables=(corr_file, e), code="54701407")
                else:
                    pass
        self.bind_object.logger.info("导入表格%s信息成功!" % corr_file)

        if mid ==1:
            main_collection = self.db['metabset_corr']
            main_collection.update({"_id":metabset_corr_id},{"$set":{"var_head":"o_id"}})  #20190617




    @report_check
    def update_table(self, str, name, main_table_id):
        try:
            self.db['metabset_corr'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {str: name}})
        except Exception as e:
            self.bind_object.set_error('更新metabset_corr主表%s字段出错:%s', variables=(str, e), code="54701408")

    @report_check
    def check_id(self, object_id, type):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('%s必须为ObjectId对象或其对应的字符串！', variables=(type), code="54701409")
        return object_id
