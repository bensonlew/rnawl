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
from mbio.packages.metabolome.common import check_metab
import glob


class MetabsetCluster(Base):
    def __init__(self, bind_object):
        super(MetabsetCluster, self).__init__(bind_object)
        self._project_type = "metabolome"

    @report_check
    def add_metabset_cluster(self, specimen=None, name=None, main_id=None, params =None, sam_tree=None,
                             matab_tree=None,list_file=None,table_type=None,metab_set=None):
        if not main_id:
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else 'MetabCluster_Origin',
                'params': params if params else '',
                'status': 'end',
                'specimen': specimen,
                'sample_tree': '',
                'metab_tree': '',
                'main_id': ''
            }
            if table_type:
                insert_data['table_type'] = table_type
            if metab_set:
                insert_data['metab_set'] = metab_set
            try:
                collection = self.db['metabset_cluster']
                metabset_cluster_id = collection.insert_one(insert_data).inserted_id
                self.update_table("main_id", metabset_cluster_id, metabset_cluster_id)
            except Exception as e:
                self.bind_object.set_error('导入metabset_cluster主表异常:%s', variables=(e), code="54700801")
        else:
            self.update_table("main_id", main_id, main_id)
            metabset_cluster_id = main_id
        if sam_tree:
            if not os.path.exists(sam_tree):
                self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(sam_tree), code="54700802")
            with open(sam_tree,"r") as f:
                specimen_tree = f.readline().strip()
                specimen = f.readline().strip().split(";")
            self.update_table("sample_tree", specimen_tree, metabset_cluster_id)
            self.update_table("specimen", specimen, metabset_cluster_id)
        if matab_tree:
            if not os.path.exists(matab_tree):
                self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(matab_tree), code="54700803")
            with open(matab_tree,"r") as f:
                metab_tree = f.readline().strip()
                metab_tree = check_metab(metab_tree)
            self.update_table("metab_tree", metab_tree, metabset_cluster_id)
        if list_file:
            if not os.path.exists(list_file):
                self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(list_file), code="54700804")
            with open(list_file,"r") as f2:
                specimen = f2.readline().strip().split("\t")
                specimen = specimen[1:len(specimen)]
            self.update_table("specimen", specimen, metabset_cluster_id)
        return metabset_cluster_id

    @report_check
    def add_metabset_cluster_detail(self, metabset_cluster_id, expression_file):
        metabset_cluster_id = self.check_id(metabset_cluster_id, "metabset_cluster_id")
        if not os.path.exists(expression_file):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(expression_file), code="54700805")
        data_list = []
        result = self.db['metabset_cluster'].find_one({'_id': metabset_cluster_id})
        if not result:
            self.bind_object.set_error('找不到%s对应的主表id', variables=(metabset_cluster_id), code="54700806")
        else:
            task_id = result['task_id']
            #samples_dic = name2id(task_id, type="task")
        #20190618
        fdir = os.path.dirname(expression_file)
        cluster_files = glob.glob(fdir+'/metab.*_cluster.xls')
        if len(cluster_files) != 1:
            #self.bind_object.set_error('metab.*_cluster.xls num is not 1')
            has_sub = 'F'

        else:
            has_sub = 'T'
            cluster_file = cluster_files[0]
            cluster_dic = {}
            with open(cluster_file) as fr:
                for line in fr:
                    line = line.rstrip('\n')
                    spline = line.split('\t')
                    for e in spline[1].split(';'):
                        if e in cluster_dic.keys():
                            self.bind_object.set_error('metabolite name is repeat')
                        cluster_dic[e] = spline[0]

        with open(expression_file, 'rb') as f:
            head = f.next()
            sams = head.strip().split("\t")
            if 'ID' in sams:
                sample_start = 3  ###2
            else:
                sample_start = 2  ##1
            for line in f:
                line = line.strip().split('\t')
                metab = check_metab(line[sample_start-2])   #样本前2位是代谢物名称
                insert_data = {
                    'cluster_id': metabset_cluster_id,
                    'metab': metab,
                    'metab_id' : line[sample_start-1]  # 添加metab_id
                }
                for i in range(sample_start, len(sams)):
                    sam_abu = float(line[i])
                    insert_data[sams[i]] = sam_abu
                if sample_start == 3:
                    insert_data['o_id'] = line[0]

                if len(cluster_files) == 1:
                    if metab not in cluster_dic.keys():
                        self.bind_object.set_error(metab + ' not in cluster_dic')
                    insert_data['sub_cluster'] = int(cluster_dic[metab])   #20190618
                data_list.append(insert_data)
            try:
                collection = self.db['metabset_cluster_detail']
                collection.insert_many(data_list)
                main_collection = self.db['metabset_cluster']
                if sample_start == 3:  #20190618
                    main_collection.update({"_id":metabset_cluster_id},{"$set":{"var_head":'o_id',"has_sub":has_sub}})
                else:
                    main_collection.update({"_id":metabset_cluster_id},{"$set":{"has_sub": has_sub}})


            except Exception as e:
                self.bind_object.set_error("导入表格%s信息出错:%s" , variables=(expression_file, e), code="54700807")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % expression_file)

    @report_check
    def update_table(self, str, name, main_table_id):
        try:
            self.db['metabset_cluster'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {str: name}})
        except Exception as e:
            self.bind_object.set_error('更新metabset_cluster主表%s字段出错:%s', variables=(str, e), code="54700808")

    @report_check
    def check_id(self, object_id, type):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('%s必须为ObjectId对象或其对应的字符串！', variables=(type), code="54700809")
        return object_id
