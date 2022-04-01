# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'

import json
import os  # ?
import gevent
import datetime
import time  # 统计拷贝时间
from bson import ObjectId
from biocluster.api.database.base import Base
from gevent import Greenlet
# from gevent.monkey import patch_all  # delete by GHD @ 20180321
from mainapp.libs.param_pack import group_detail_sort
from id_convert import name2id
from biocluster.config import Config
from mbio.packages.metagenomic.common import get_old_mongo


class CopyDemo(Base):
    def __init__(self, old_task_id, new_task_id, new_project_sn,
                 new_member_id, new_member_type,
                 old_db_version=1):
        '''
        old_db_version 表明old_task_id来源的项目是否在老的db_version
        '''
        super(CopyDemo, self).__init__()
        self._project_type = 'metagenomic'
        self._old_task_id = old_task_id
        self.old_db_version = old_db_version
        if old_db_version == 0:
            self.old_db = get_old_mongo()[1]
        else:
            self.old_db = self.db
        self._new_task_id = new_task_id
        self._new_project_sn = new_project_sn
        self._new_member_id = new_member_id
        self._new_member_type = new_member_type
        self.specimen_id_dict = {}
        self.group_id_dict = {}
        self.env_id_dict = {}
        self.geneset_id_dict = {}
        self.anno_id_dict = {}
        self.group_detail_dict = {}  # 是否需要？
        self.hcluster_tree_id_dict = {}
        # self.otu_id_dict = {}  # 改成注释表的
        self.all_greenlets = []
        self._exchange_dict = {
            'specimen_id': self.specimen_id_dict,
            'specimen_name': self.specimen_id_dict,  # 一些数据表中名称为specimen_name
            'group_id': self.group_id_dict,
            'env_id': self.env_id_dict,
            'geneset_id': self.geneset_id_dict,
            'anno_id': self.anno_id_dict,
            'group_detail': self.group_detail_dict,
            'hcluster_tree_id': self.hcluster_tree_id_dict,
            # 'otu_id': self.otu_id_dict,
        }

    def copy_member_id(self):
        """
        复制sg_task的数据
        :return:
        """
        coll = self.db['sg_task']
        find = coll.find_one({'task_id': self._old_task_id})
        if not find:
            raise Exception('运行错误：找不到demo任务相关信息')
        find['task_id'] = self._new_task_id
        find['member_id'] = self._new_member_id
        find['member_type'] = self._new_member_type
        find.pop('_id')
        find['project_sn'] = self._new_project_sn
        find['is_demo'] = 1
        try:
            find['demo_id'] = self._old_task_id
        except:
            find.pop('demo_id')
            find['demo_id'] = self._old_task_id
        coll.insert_one(find)

    def copy_specimen_group(self):
        """
        复制group表
        :return:
        """
        old_coll = self.old_db['specimen_group']
        coll = self.db['specimen_group']
        finds = old_coll.find({'task_id': self._old_task_id})
        news = []
        old_group_ids = []
        for i in finds:
            i['task_id'] = self._new_task_id
            i['project_sn'] = self._new_project_sn
            old_group_ids.append(str(i.pop('_id')))
            for one in i['specimen_names']:
                for sp in one:
                    try:
                        one[one.index(sp)] = self.specimen_id_dict[sp]
                    except:
                        print "can't find sp: %s" % sp
                        print "specimen_id_dict is : %s" % self.specimen_id_dict
                """
                for sp in one.copy():  # use this if specimen_names is hash
                    try:
                        one[self.specimen_id_dict[sp]] = one[sp]
                    except:
                        print "sp is :" + sp
                        print "one is :"
                        print one
                    one.pop(sp)
                """
            news.append(i)
        if news:
            result = coll.insert_many(news)
            self.group_id_dict = dict(zip(old_group_ids, [str(one) for one in result.inserted_ids]))
        else:
            self.group_id_dict = {}
        self.group_id_dict['all'] = 'all'  # 特殊ID
        self.group_id_dict[None] = None  # 特殊ID
        self.group_id_dict[''] = None  # 特殊ID
        self._exchange_dict['group_id'] = self.group_id_dict
        return self.group_id_dict

    def copy_specimen_graphic(self):
        """
        拷贝样本原始序列作图表格
        :return:
        """
        old_coll = self.old_db['specimen_graphic']
        coll = self.db['specimen_graphic']
        finds = old_coll.find({'task_id': self._old_task_id})
        news = []
        for i in finds:
            i.pop('_id')
            i["task_id"] = self._new_task_id
            if 'specimen_name' in i:
                i["specimen_name"] = self.specimen_id_dict[str(i['specimen_name'])]
            news.append(i)
        if news:
            coll.insert_many(news)
        else:
            print "specimen_graphic表没有找到task_id:%s相关信息" % self._old_task_id
            # self.bind_object.logger.info("specimen_graphic表没有找到task_id:%s相关信息" % self._old_task_id)

    def get_specimen_id_dict(self):
         """
         获得新旧的样品id，需复制data_stat_detail表之后获得
         :return:
         """
         old_id = name2id(self._old_task_id, "task", self.old_db_version)
         new_id = name2id(self._new_task_id, "task")
         for name in old_id.keys():
             self.specimen_id_dict[old_id[name]] = new_id[name]
         # print "TAKE CARE !!!"
         # print self.specimen_id_dict
         # return self.specimen_id_dict

    def copy_datastat_details(self, collection, main_field, change_dict):
        """
        用于copy deta_stat_detail,specimen_graphic表,并提取样品id
        params collection: detail表名称
        params main_field: 主表字段名称
        params change_dict: 主表新旧替换字典
        """
        time_start = datetime.datetime.now()
        old_coll = self.old_db[collection]
        coll = self.db[collection]
        type_raw = []  # type = raw 的data_stat_detail
        type_other = []  # type != raw 的data_stat_detail
        for old, new in change_dict.items():
            if (not isinstance(old, ObjectId)) and len(old) != 24:
                print "%s表的%s字段不是ObjectId类型：%s len is %s" % (collection, main_field, old, len(old))
                continue
            finds = old_coll.find({main_field: ObjectId(old)})
            data_main = self.old_db['data_stat'].find_one({'_id': ObjectId(old)})
            if data_main['type'] == 'raw':
                raw = True
            else:
                raw = False
            for i in finds:
                i.pop('_id')
                i[main_field] = ObjectId(new)
                if self.old_db_version == 0 and 'origin_id' in i:  # 从不同mongo服务器复制时需要去掉origin_id字段
                    i.pop('origin_id')
                if raw:
                    type_raw.append(i)
                else:
                    type_other.append(i)
        # a 先导入raw，b 然后获取新老specimen_id的对应关系
        # c 接着修改非raw的specimen_id为新的, d 最后导入非raw的data_stat_detail
        if type_raw:  # a
            coll.insert_many(type_raw)
        self.get_specimen_id_dict()  # b
        for one in type_other:  # c
            one['specimen_name'] = self.specimen_id_dict.get(one['specimen_name'], one['specimen_name'])
        coll.insert_many(type_other)  # d
        self.copy_specimen_graphic()
        time_end = datetime.datetime.now()
        run_time = (time_end - time_start).seconds
        print "{}复制运行时间: {}s".format(collection, run_time)

    def _copy_main_details(self, collection, main_field, change_dict, others_position=[]):
        """
        公共模块，一般用于更新detail表，根据提供的主表id字段名，和主表新旧ID字典，进行查找，再复制替换，others_position用于更新主表ID之外其他需要更新的ID

        params collection: detail表名称
        params main_field: 主表字段名称
        params change_dict: 主表新旧替换字典，一般来源于 copy_collection_with_change 的返回字典
        params others_position: detail表中除了主表还需要更新的字段，
            只能是 specimen_id,group_id,env_id,otu_id,alpha_diversity_id,newick_tree_id,specimen_distance_id
        """
        time_start = datetime.datetime.now()
        old_coll = self.old_db[collection]
        coll = self.db[collection]
        if type(change_dict) == tuple:
            print "ERROR: change_dict is tuple in _copy_main_details"
            print "collection is %s, main_field is %s, change_dict: " % (collection, main_field)
            print change_dict
            return
        for old, new in change_dict.items():
            if (not isinstance(old, ObjectId)) and len(old) != 24:
                print "%s表的%s字段不是ObjectId类型：%s len is %s" % (collection, main_field, old, len(old))
                continue
            finds = old_coll.find({main_field: ObjectId(old)})
            news = []
            for i in finds:
                i.pop('_id')
                i[main_field] = ObjectId(new)
                for position in others_position:
                    i[position] = self.exchange_ObjectId(position, i[position])
                if self.old_db_version == 0:
                    if 'specimen_name' in i:  # 老数据库的 specimen id 换成新数据库的
                        i['specimen_name'] = self.specimen_id_dict.get(i['specimen_name'], i['specimen_name'])
                    new_i = {}
                    for k in i:
                        if k in self.specimen_id_dict:  # 表头是老数据库的 specimen id 换成新的
                            n_k = self.specimen_id_dict[i]
                            new_i[n_k] = i.pop(k)
                    i.update(new_i)
                news.append(i)
            if news:
                coll.insert_many(news)
            else:
                print 'WARNING: 主表:{}没有detail表信息，请注意数据合理性,collection:{}'.format(old, collection)
        time_end = datetime.datetime.now()
        run_time = (time_end - time_start).seconds
        print "{}复制运行时间: {}s".format(collection, run_time)

    def copy_main_details(self, collection, main_field, change_dict, others_position=[], join=True):
        greenlet = Greenlet(self._copy_main_details, collection, main_field, change_dict, others_position)
        greenlet.start()
        if join is True:
            greenlet.join()
            return greenlet.value
        self.all_greenlets.append(greenlet)
        return greenlet

    def copy_collection_with_change(self, collection, change_positions=[], update_sg_status=False, targetcoll=None):
        """
        公共模块，一般用于导入主表数据，依靠task_id进行查询，修改change_positions提供的字段，相应修改ID为新的，同时更新params中的数据ID
        params collection: 主表名称
        params change_positions: 需要替换的ID,可用为specimen_id,group_id...
        params update_sg_status: 更新sg_status表
        params targetcoll: 更新到特定集合， 默认与collection参数相同
        """
        # ref_rna, 没有运行Greenlet?
        coll = self.old_db[collection]
        if targetcoll:
            targetcoll = self.db[targetcoll]
        else:
            targetcoll = self.db[collection]
        finds = coll.find({'task_id': self._old_task_id})
        news = []
        olds = []
        for i in finds:
            i['task_id'] = self._new_task_id
            if collection == 'geneset' and i['name'] != 'GENESET_Origin':
                continue
            if 'params' in i and i['params'] == "":
                continue  # 这种情况属于在页面上删除了该记录，所以不进行复制 @20180301
            if 'project_sn' in i:
                i['project_sn'] = self._new_project_sn
            olds.append(str(i.pop('_id')))
            for position in change_positions:
                if position in i:
                    i[position] = self.exchange_ObjectId(position, i[position])
            if 'params' in i:
                print "进行 %s 表的复制param" % collection
                i['params'] = self.params_exchange(i['params'])
                print "表 %s 复制结束" % collection
            news.append(i)
        if news:
            result = targetcoll.insert_many(news)
            if update_sg_status:
                self.insert_new_status(collection, news, result.inserted_ids)
            return dict(zip(olds, [str(one) for one in result.inserted_ids]))
        else:
            return {}

    def exchange_ObjectId(self, key, thisObjectId):
        """
        用于替换id，key是该ID的字段名，thisObjectId是旧的ID(ObjectId类型)
        """
        if not self._exchange_dict[key].keys():
            return None
        elif not str(thisObjectId) in self._exchange_dict[key].keys():  # mark
            return str(thisObjectId)
        if isinstance(thisObjectId, ObjectId):
            return ObjectId(self._exchange_dict[key][str(thisObjectId)])
        else:
            return self._exchange_dict[key][thisObjectId]  # 不是ObjectId时直接返回也是字符串

    def params_exchange(self, params_str):
        """
        专门用于params的数据ID替换, 注：metagenomic 需替换geneset_id,anno_id,group_id
        :return:
        """
        if params_str == None:
            return None
        try:
            params = json.loads(params_str)
        except Exception:
            print("WARNNING: 非json格式的params: {}".format(params_str))
            return params_str
        if not params:
            return None
        if 'group_detail' in params:
            for one_group in params['group_detail']:
                try:  # mark
                    params['group_detail'][one_group] = [self.specimen_id_dict[one_sp] for one_sp in
                                                     params['group_detail'][one_group]]
                except Exception:  # mark
                    print "group_detail have something wrong: "
                    print self.specimen_id_dict
            params['group_detail'] = group_detail_sort(params['group_detail'])
            if 'second_group_detail' in params:
                if params['second_group_detail'] and params['second_group_detail'] != "null":
                    for one_group in params['second_group_detail']:
                        params['second_group_detail'][one_group] = [self.specimen_id_dict[one_sp] for one_sp in
                                                                    params['second_group_detail'][one_group]]
                    params['second_group_detail'] = group_detail_sort(params['second_group_detail'])
                    if 'second_group_id' in params:
                        params['second_group_id'] = self.group_id_dict[params['second_group_id']]
        # 添加其他可能的字段
        if 'anno_id' in params and params['anno_id']:
            if str(params['anno_id']) in self._exchange_dict['anno_id'].keys():
                params['anno_id'] = self._exchange_dict['anno_id'][str(params['anno_id'])]
            else:
                print "anno_id  %s is not in dict: " % params['anno_id']
                print self._exchange_dict['anno_id']
                print params_str
        if 'geneset_id' in params:
            if str(params['geneset_id']) in self.geneset_id_dict.keys():
                params['geneset_id'] = self.geneset_id_dict[str(params['geneset_id'])]
            else:
                print "geneset_id %s is not in dict: " % params['geneset_id']
                print self.geneset_id_dict
        if 'group_id' in params:
            if str(params['group_id']) in self.group_id_dict:
                params['group_id'] = self.group_id_dict[str(params['group_id'])]
            else:
                print "group_id %s is not in dict: " % params['group_id']
                print self.group_id_dict
        if 'specimen_group' in params:
            if str(params['specimen_group']) in self.group_id_dict:
                params['specimen_group'] = self.group_id_dict[str(params['specimen_group'])]
            else:
                print "specimen_group %s is not in dict: " % params['specimen_group']
                print self.group_id_dict
        if 'env_id' in params:
            if str(params['env_id']) in self.env_id_dict:
                params['env_id'] = self.env_id_dict[str(params['env_id'])]
            else:
                print "env_id %s is not in dict: " % params['env_id']
                print self.env_id_dict
        return json.dumps(params, sort_keys=True, separators=(',', ':'))

    def insert_new_status(self, collection, main_docs, ids):
        """
        导入mongo表sg_status数据信息
        :param collection:
        :param main_docs:
        :param ids:
        :return:
        """
        coll = self.db['sg_status']
        news = []
        for index, doc in enumerate(main_docs):
            try:
                submit_location = json.loads(doc['params'])['submit_location']
            except Exception:
                print("WARNING: params参数没有submit_location字段, Doc:{}".format(doc))
                submit_location = None
            status = {
                "status": doc['status'],
                "table_id": ids[index],
                "time": doc['created_ts'],
                "task_id": self._new_task_id,
                "params": doc['params'],
                "table_name": doc['name'],
                "submit_location": submit_location,
                "type_name": collection,
                "is_new": "new",
                "desc": doc['desc'] if 'desc' in doc else None
            }
            news.append(status.copy())
        if news:
            coll.insert_many(news)

    # def abund_table_path(self):
    #     abund_table_dict = self.copy_collection_with_change("abund_table_path")  # 此表已不用

    def anno_table(self):
        anno_list = ["ardb", "card", "cazy", "cog", "kegg", "nr", "vfdb", "overview"]
        anno_detail = {
            "ardb": ["arg", "class", "type"],
            "card": ["aro", "class"],
            "cazy": ["class", "family"],
            "cog": ["category", "function", "nog"],
            "kegg": ["enzyme", "gene", "module", "orthology", "pathway"],
            "nr": ["detail"],
            "vfdb": ["pie", "vfs"],
            "overview": [],
        }
        for anno in anno_list:
            collection = "anno_" + anno
            # 　anno_main_dict = self.copy_collection_with_change(collection, change_positions=['group_id', 'group_detail'])
            anno_main_dict = self.copy_collection_with_change(collection, change_positions=['geneset_id'])
            self._exchange_dict['anno_id'].update(anno_main_dict)
            for one in anno_detail[anno]:
                collection_detail = collection + "_" + one
                # position = anno + "_id"
                position = "anno_overview_id" if anno == "overview" else anno + "_id"
                self.copy_main_details(collection_detail, position, anno_main_dict, join=False)

    def anno_personal(self):
        anno_list = ["cyps", "mvir", "pfam", "phi", "probio", "qs", "sec", "tcdb","ttss", "go"]
        for anno in anno_list:
            collection = "anno_" + anno
            # 　anno_main_dict = self.copy_collection_with_change(collection, change_positions=['group_id', 'group_detail'])
            anno_main_dict = self.copy_collection_with_change(collection, change_positions=['geneset_id'])
            self._exchange_dict['anno_id'].update(anno_main_dict)

    def anosim(self):
        #　anosim_id_dict = self.copy_collection_with_change("anosim",
        #                                                    change_positions=['anno_id', 'geneset_id', 'group_id'])
        anosim_id_dict = self.copy_collection_with_change("anosim")
        self.copy_main_details('anosim_detail', 'anosim_id', anosim_id_dict, join=False)

    def assem_stat(self):
        assem_id_dict = self.copy_collection_with_change("assemble_stat")
        self.copy_main_details('assemble_stat_bar', 'assem_id', assem_id_dict, join=False)
        # self.copy_main_details('assemble_stat_bar', 'assem_id', assem_id_dict, others_position=['specimen_name'], join=False)
        self.copy_main_details('assemble_stat_detail', 'assem_id', assem_id_dict, join=False)
        # self.copy_main_details('assemble_stat_detail', 'assem_id', assem_id_dict, others_position=['specimen_name'], join=False)

    def beta_deversity(self):
        # 　beta_id_dict = self.copy_collection_with_change("beta_diversity", change_positions=['anno_id', 'geneset_id', 'group_id'])
        beta_id_dict = self.copy_collection_with_change("beta_diversity", change_positions=['env_id'])
        self.copy_main_details('beta_diversity_detail', 'beta_diversity_id', beta_id_dict, join=False)

    def composition(self):
        # composition_id_dict = self.copy_collection_with_change("composition", change_positions=['anno_id', 'geneset_id', 'group_id'])
        # specimen_list字段是样品id以逗号分割
        composition_id_dict = self.copy_collection_with_change("composition")
        self.copy_main_details('composition_detail', 'composition_id', composition_id_dict, join=False)

    def data_stat(self):
        data_stat_id_dict = self.copy_collection_with_change("data_stat")
        self.copy_datastat_details('data_stat_detail', 'data_stat_id', data_stat_id_dict)
        # self.copy_main_details('data_stat_detail', 'data_stat_id', data_stat_id_dict, join=False)
        # self.get_specimen_id_dict()
        # self.copy_specimen_graphic()

    def enterotype(self):
        enterotype_id_dict = self.copy_collection_with_change("enterotype")
        self.copy_main_details('enterotype_detail', 'enterotype_id', enterotype_id_dict, join=False)
        self.copy_main_details('enterotype_detail_cluster', 'enterotype_id', enterotype_id_dict, join=False)

    def env_vif(self):
        vif_id_dict = self.copy_collection_with_change("env_vif")
        self.copy_main_details('env_vif_detail', 'vif_id', vif_id_dict, join=False)

    def geneset(self):
        self.geneset_id_dict = self.copy_collection_with_change("geneset")
        self._exchange_dict['geneset_id'] = self.geneset_id_dict
        self.copy_main_details('geneset_bar', 'geneset_id', self.geneset_id_dict, join=False)
        self.copy_main_details('geneset_readsn', 'geneset_id', self.geneset_id_dict, join=False)
        self.copy_main_details('geneset_readsr', 'geneset_id', self.geneset_id_dict, join=False)
        # readsn/readsr含有以样品id为字段的字段

    def hcluster_tree(self):
        self.hcluster_tree_id_dict = self.copy_collection_with_change("hcluster_tree")
        self._exchange_dict['hcluster_tree_id'] = self.hcluster_tree_id_dict

    def heatmap_cor(self):
        heatmap_cor_id_dict = self.copy_collection_with_change('heatmap_cor', change_positions=['env_id'])
        self.copy_main_details('heatmap_cor_detail', 'heatmap_cor_id', heatmap_cor_id_dict, join=False)

    def lefse(self):
        lefse_id_dict = self.copy_collection_with_change('lefse', change_positions=['geneset_id', 'group_id', 'anno_id'])
        self.copy_main_details('lefse_detail', 'species_lefse_id', lefse_id_dict, join=False)

    def mantel_test(self):
        mantel_id_dict = self.copy_collection_with_change('mantel_test')
        self.copy_main_details('mantel_test_detail', 'mantel_id', mantel_id_dict, join=False)

    def metastat(self):
        # "metastat, metastat_plot, metastat_detail-metastat_id"
        metastat_id_dict = self.copy_collection_with_change('metastat')
        self.copy_main_details('metastat_detail', 'metastat_id', metastat_id_dict, join=False)
        self.copy_main_details('metastat_plot', 'metastat_id', metastat_id_dict, join=False)

    def network(self):
        # "network, network_degree, network_link, network_node-network_id"
        network_id_dict = self.copy_collection_with_change('network', change_positions=["geneset_id"])
        self.copy_main_details('network_degree', 'network_id', network_id_dict, join=False)
        self.copy_main_details('network_link', 'network_id', network_id_dict, join=False)
        self.copy_main_details('network_node', 'network_id', network_id_dict, join=False)

    def network_cor(self):
        # "net_work_cor, network_cor_degree, network_cor_link, network_cor_code-network_cor_id"
        network_cor_id_dict = self.copy_collection_with_change('network_cor', change_positions=["geneset_id"])
        self.copy_main_details('network_cor_degree', 'network_cor_id', network_cor_id_dict, join=False)
        self.copy_main_details('network_cor_link', 'network_cor_id', network_cor_id_dict, join=False)
        self.copy_main_details('network_cor_node', 'network_cor_id', network_cor_id_dict, join=False)

    def permanova(self):
        # "permanova, permanova_detail-permanova_id"
        permanova_id_dict = self.copy_collection_with_change('permanova')
        self.copy_main_details('permanova_detail', 'permanova_id', permanova_id_dict, join=False)

    def predict_gene(self):
        predict_gene_id_dict = self.copy_collection_with_change("predict_gene")
        # self.copy_main_details("predict_gene_bar", "predict_gene_id", predict_gene_id_dict, others_position=['specimen_name'], join=False)
        # self.copy_main_details("predict_gene_detail", "predict_gene_id", predict_gene_id_dict, others_position=['specimen_name'], join=False)
        self.copy_main_details("predict_gene_bar", "predict_gene_id", predict_gene_id_dict, join=False)
        self.copy_main_details("predict_gene_detail", "predict_gene_id", predict_gene_id_dict, join=False)
        predict_total_dict = self.copy_collection_with_change("predict_gene_total")
        # self.copy_main_details("predict_gene_total", "predict_gene_id", predict_gene_id_dict, join=False)

    def regression(self):
        # "regression, regression_curve, regression_line-regression_id"
        regression_id_dict = self.copy_collection_with_change("regression")
        self.copy_main_details("regression_curve", "regression_id", regression_id_dict, join=False)
        self.copy_main_details("regression_line", "regression_id", regression_id_dict, join=False)

    def specimen_distance(self):
        distance_id_dict = self.copy_collection_with_change('specimen_distance', change_positions=['hcluster_tree_id'])
        self.copy_main_details('specimen_distance_detail', 'specimen_distance_id', distance_id_dict, join=False)

    def venn(self):
        venn_id_dict = self.copy_collection_with_change('venn')
        self.copy_main_details('venn_detail', 'venn_id', venn_id_dict, join=False)
        self.copy_main_details('venn_graph', 'venn_id', venn_id_dict, join=False)
        self.copy_main_details('venn_pie_data', 'venn_id', venn_id_dict, join=False)

    def contribute(self):
        contribute_id_dict = self.copy_collection_with_change('contribute', change_positions=['geneset_id'])
        self.copy_main_details('contribute_detail', 'contribute_id', contribute_id_dict, join=False)

    def gl(self, func):
        greenlet = Greenlet(func)
        greenlet.start()
        self.all_greenlets.append(greenlet)

    def check_copied(self):
        geneset = self.db['geneset']
        find = geneset.find_one({'task_id': self._new_task_id})
        return find and True

    def from_task(self):
        if self.check_copied():
            print "skip the copy step, because the mongo data have been copied"
        else:
            self.data_stat()  # 需要首先复制data_stat，用于样本id的获取
            self.copy_specimen_group()  # 样本分组表
            self.env_id_dict = self.copy_collection_with_change('env')
            self.copy_main_details('env_detail', 'env_id', self.env_id_dict, others_position=['specimen_id'], join=False)
            self._exchange_dict['env_id'] = self.env_id_dict
            for func in [self.assem_stat, self.predict_gene, self.geneset]:
                self.gl(func)
            gevent.joinall(self.all_greenlets)
        return True

    def run(self, from_task=False):
        # patch_all()  # delete by GHD @ 20180321
        if from_task:
            return self.from_task()
        print "STRART COPYMEMBERID"
        self.copy_member_id()
        print "START datastat"
        self.data_stat()
        print "specimen_id is "
        print self.specimen_id_dict
        print "START GROUP!!!!"
        self.copy_specimen_group()
        print "GROUP END!!!!"
        print "group_id is: "
        print self.group_id_dict
        self.env_id_dict = self.copy_collection_with_change('env')
        print "env_id is: "
        print self.env_id_dict
        self.copy_main_details('env_detail', 'env_id', self.env_id_dict, others_position=['specimen_id'], join=False)
        self._exchange_dict['env_id'] = self.env_id_dict
        self.assem_stat()
        self.predict_gene()
        self.geneset()
        self.anno_table()
        self.hcluster_tree()
        self.anno_personal()
        greenlet = Greenlet(self.anosim)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.beta_deversity)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.composition)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.enterotype)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.env_vif)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.heatmap_cor)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.lefse)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.mantel_test)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.metastat)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.network)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.network_cor)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.permanova)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.regression)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.specimen_distance)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.venn)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.contribute)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        gevent.joinall(self.all_greenlets)
        import socket
        reload(socket)


if __name__ == '__main__':
    start_time = time.time()
    copy_task = CopyDemo('tsg_26035_test2', 'demotest_metagenome41', '111112', '111113', '2')  # input params from task existed
    copy_task.run()
    end_time = time.time()
    print "total time: {}".format(end_time - start_time)
