# -*- coding: utf-8 -*-

from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
# from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from mbio.packages.metabolome.common import check_metab
import glob
import pandas as pd


class EnrichCluster(Base):
    def __init__(self, bind_object):
        super(EnrichCluster, self).__init__(bind_object)
        self._project_type = "metabolome"

    @report_check
    def add_enrich_cluster(self, specimen=None, name=None, main_id=None, params =None, metabset_tree=None,
                             pathway_tree=None):
        if not main_id:
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else 'MetabEnrichCluster_Origin',
                'params': params if params else '',
                'status': 'end',
                'specimen': specimen,
                'sample_tree': '',
                'metab_tree': '',
                'main_id': ''
            }
            try:
                collection = self.db['metabset_kegg_heatmap']
                metabset_cluster_id = collection.insert_one(insert_data).inserted_id
                self.update_table("main_id", metabset_cluster_id, metabset_cluster_id)
            except Exception as  e:
                self.bind_object.set_error('导入metabset_kegg_heatmap 主表异常:%s', variables=(e), code="54700801")
        else:
            self.update_table("main_id", main_id, main_id)
            metabset_cluster_id = main_id
        if metabset_tree:
            if not os.path.exists(metabset_tree):
                self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(metabset_tree), code="54700802")
            with open(metabset_tree,"r") as f:
                specimen_tree = f.readline().strip()
            self.update_table("metabset_tree", specimen_tree, metabset_cluster_id)
        if pathway_tree:
            if not os.path.exists(pathway_tree):
                self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(pathway_tree), code="54700803")
            with open(pathway_tree,"r") as f:
                pathway_tree = f.readline().strip()
            self.update_table("pathway_tree", pathway_tree, metabset_cluster_id)

        return metabset_cluster_id

    @report_check
    def add_enrich_cluster_detail(self, heatmap_id, expression_file, enrich_pvalue=None, metab_ids_file=None):
        heatmap_id = self.check_id(heatmap_id, "heatmap_id")
        if not os.path.exists(expression_file):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(expression_file), code="54700805")
        data_list = []
        result = self.db['metabset_kegg_heatmap'].find_one({'_id': heatmap_id})
        if not result:
            self.bind_object.set_error('找不到%s对应的主表id', variables=(heatmap_id), code="54700806")
        else:
            task_id = result['task_id']

        fdir = os.path.dirname(expression_file)
        cluster_file = fdir+'/pathway.cluster.txt'
        cluster_dic = {}
        if os.path.exists(cluster_file):
            with open(cluster_file) as fr:
                for line in fr:
                    line = line.rstrip('\n')
                    spline = line.split('\t')
                    for e in spline[1].split(';'):
                        if e in cluster_dic.keys():
                            self.bind_object.set_error('%s name is repeat'% e)
                        cluster_dic[e] = spline[0]

        with open(expression_file, 'rb') as f:
            head = f.next()
            sams = head.rstrip().split("\t")[3:]
            metabset_name = []
            metabset_name_2 =[]
            for id,name in enumerate(sams,1):
                metabset_name.append(('metabset'+str(id),name))
                metabset_name_2.append((name, 'metabset'+str(id)))
            metabset_name_dic = dict(metabset_name)
            metabset_name_2_dic = dict(metabset_name_2)

            for line in f:
                line = line.strip().split('\t')
                pathway_id = line[0]
                insert_data = {
                    'heatmap_id': heatmap_id,
                    'pathway_id': pathway_id,
                    'pathway_desc' : line[1],
                    'database' :  line[2]
                }
                for col_id, keys in enumerate(metabset_name,3):
                    sam_abu = float(line[col_id])
                    insert_data[keys[0]] = sam_abu
                if os.path.exists(cluster_file):
                    if pathway_id not in cluster_dic.keys():
                        self.bind_object.set_error(pathway_id + ' not in cluster_dic')
                    insert_data['ncluster'] = int(cluster_dic[pathway_id])   #20190618
                data_list.append(insert_data)

        if enrich_pvalue:
            p_info = {}
            with open(enrich_pvalue) as fr:
                head = fr.readline()
                sams = head.rstrip().split('\t')[1:]
                for line in fr:
                    spline = line.strip().split('\t')
                    id = spline[0]
                    p_info[id] ={}
                    for index,s in enumerate(sams,1):
                        mongo_k = metabset_name_2_dic[s] + '_pvalue'
                        p_info[id][mongo_k] = float(spline[index])

            for index in range(len(data_list)):
                pathway_id = data_list[index]['pathway_id']
                pvalue = p_info[pathway_id]
                data_list[index].update(pvalue)

        if metab_ids_file:
            metab_ids_info = {}
            metab_ids = pd.read_table(metab_ids_file, sep='\t',index_col=0)
            metab_ids.fillna("",inplace=True)  #20200629
            metab_ids = metab_ids.to_dict('index')

            for map_id in metab_ids:
                metab_ids_info[map_id] = {}
                for name in metabset_name_2_dic.keys():
                    ref_id = metabset_name_2_dic[name]
                    metab_ids_info[map_id].update({ref_id+'_num': metab_ids[map_id][name+'_Study_num'],ref_id+'_ids':metab_ids[map_id][name+'_Metab_ids']})

            for index in range(len(data_list)):
                pathway_id = data_list[index]['pathway_id']
                data_list[index].update(metab_ids_info[pathway_id])

        try:
            collection = self.db['metabset_kegg_heatmap_detail']
            collection.insert_many(data_list)

            self.update_table('metabset_name', metabset_name_dic, heatmap_id)

        except Exception as e:
            self.bind_object.set_error("导入表格%s信息出错:%s" , variables=(expression_file, e), code="54700807")
        else:
            self.bind_object.logger.info("导入表格%s信息成功!" % expression_file)

    @report_check
    def update_table(self, str, name, main_table_id):
        main_table_id = self.check_id(main_table_id, 'heatmap_id')
        try:
            self.db['metabset_kegg_heatmap'].update_one({'_id':main_table_id}, {'$set': {str: name}})
        except Exception as e:
            self.bind_object.set_error('更新metabset_kegg_heatmap主表%s字段出错:%s', variables=(str, e), code="54700808")

    @report_check
    def check_id(self, object_id, type):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('%s必须为ObjectId对象或其对应的字符串！', variables=(type), code="54700809")
        return object_id
