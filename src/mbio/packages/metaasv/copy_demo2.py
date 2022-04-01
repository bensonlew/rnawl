# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'@20200523

import json
import datetime
import gevent
import time  # 统计拷贝时间
from bson import ObjectId
from biocluster.api.database.base import Base
from gevent import Greenlet
from gevent.monkey import patch_all


class CopyDemo(Base):
    def __init__(self, old_task_id, new_task_id, new_project_sn, new_member_id, new_member_type):
        super(CopyDemo, self).__init__()
        self._project_type = 'metaasv'
        self._old_task_id = old_task_id
        self._new_task_id = new_task_id
        self._new_project_sn = new_project_sn
        self._new_member_id = new_member_id
        self._new_member_type = new_member_type
        self.specimen_id_dict = {}
        self.all_greenlets = []

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
        print '进行{}表的复制'.format(collection)
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
            finds = coll.find({main_field: ObjectId(old)})
            news = []
            for i in finds:
                i.pop('_id')
                i[main_field] = ObjectId(new)
                for position in others_position:
                    i[position] = self.exchange_ObjectId(position, i[position])
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
        coll = self.db[collection]
        if targetcoll:
            targetcoll = self.db[targetcoll]
        else:
            targetcoll = self.db[collection]
        finds = coll.find({'task_id': self._old_task_id, 'status': 'end'})
        news = []
        olds = []
        for i in finds:
            i['task_id'] = self._new_task_id
            if 'project_sn' in i:
                i['project_sn'] = self._new_project_sn
            olds.append(str(i.pop('_id')))
            for position in change_positions:
                if position in i:
                    i[position] = self.exchange_ObjectId(position, i[position])
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

    def alpha_diversity(self):
        """
        alpha_diversity多样性
        :return:
        """
        self.copy_collection_with_change('alpha_diversity')

    def alpha_diversity_diff(self):
        """
        alpha_diversity多样性组间差异检验
        :return:
        """
        self.copy_collection_with_change("alpha_diversity_diff")

    def asv(self):
        """
        asv分类学分析
        :return:
        """
        self.copy_collection_with_change("asv")
        details = self.db['sg_otu_specimen'].find({"asv_id":old_main_id},{"_id":0})
        insert_data = []
        for detail in details:
            detail['otu_id'] = main_table_id
            old_sample_id = str(detail['specimen_id'])
            detail['specimen_id'] = ObjectId(self.sample_map[old_sample_id])
            insert_data.append(detail)
        self.db['sg_otu_specimen'].insert_many(insert_data)

        self.copy_detail(old_main_id, main_table_id, 'sg_otu_detail','otu_id')
        self.copy_detail(old_main_id, main_table_id, 'sg_otu_detail_level','otu_id')
        self.copy_detail(old_main_id, main_table_id, 'sg_otu_seq','otu_id')
        self.copy_detail(old_main_id, main_table_id, 'sg_otu_summary','otu_id')


    def asv_set(self):
        """
        asv_set基因集
        :return:
        """
        self.copy_collection_with_change('asv_set')

    def barpie(self):
        """
        群落组成分析--barpie
        :return:
        """
        self.copy_collection_with_change('barpie')

    def circos(self):
        """
        circos图
        :return:
        """
        self.copy_collection_with_change("circos")

    def cor_heatmap(self):
        """
        相关性heatmap图
        :return:
        """
        self.copy_collection_with_change("cor_heatmap")

    def cor_network(self):
        """
        单因素网络图
        :return:
        """
        self.copy_collection_with_change("cor_network")

    def data_stat(self):
        """
        样本信息统计
        :return:
        """
        data_stat_dict = self.copy_collection_with_change("data_stat")
        self.copy_datastat_details('data_stat_detail', 'stat_id', data_stat_dict)


    def dbrda(self):
        """
        dbrda分析
        :return:
        """
        self.copy_collection_with_change('dbrda')

    def enterotype(self):
        """
        样本菌群分析
        :return:
        """
        self.copy_collection_with_change('enterotype')

    def env(self):
        """
        环境因子
        :return:
        """
        self.copy_collection_with_change('env')

    def env_group(self):
        """
        环境因子分组
        :return:
        """
        self.copy_collection_with_change('env_group')

    def funguild(self):
        """
        funguild功能预测分析
        :return:
        """
        self.copy_collection_with_change('funguild')

    def hcluster(self):
        """
        层级聚类分析
        :return:
        """
        self.copy_collection_with_change("hcluster")

    def heatmap(self):
        """
        群落heatmap分析
        :return:
        """
        self.copy_collection_with_change("heatmap")

    def lefse(self):
        """
        LEFse分析
        :return:
        """
        self.copy_collection_with_change("lefse")

    def mantel(self):
        """
        mantel_test分析
        :return:
        """
        self.copy_collection_with_change("mantel")

    def multiple_group(self):
        """
        多组比较
        :return:
        """
        self.copy_collection_with_change("multiple_group")

    def two_group(self):
        """
        两组比较
        :return:
        """
        self.copy_collection_with_change("two_group")

    def two_sample(self):
        """
        两样本比较
        :return:
        """
        self.copy_collection_with_change("two_sample")

    def nmds(self):
        """
        NMDS分析
        :return:
        """
        self.copy_collection_with_change("nmds")

    def pan_core(self):
        """
        pan_core分析
        :return:
        """
        self.copy_collection_with_change("pan_core")

    def pca(self):
        """
        PCA分析
        :return:tnf
        """
        self.copy_collection_with_change("pca")

    def pcoa(self):
        """
        PCoA分析
        :return:
        """
        self.copy_collection_with_change("pcoa")

    def permanova(self):
        """
        PERMANOVA分析
        :return:
        """
        self.copy_collection_with_change("permanova")

    def personal_phylo_tree(self):
        """
        个性化系统发生树
        :return:
        """
        self.copy_collection_with_change("personal_phylo_tree")

    def phylo_tree(self):
        """
        系统发生树
        :return:
        """
        self.copy_collection_with_change("phylo_tree")

    def picrust2(self):
        """
        picrust2功能预测
        :return:
        """
        self.copy_collection_with_change("picrust2")

    def randomForest(self):
        """
        随机森林分析
        :return:tnf
        """
        self.copy_collection_with_change("randomForest")

    def rank_abundance(self):
        """
        rank_abundance分析
        :return:
        """
        self.copy_collection_with_change("rank_abundance")

    def rarefaction(self):
        """
        rarefaction分析
        :return:
        """
        self.copy_collection_with_change("rarefaction")

    def rda_cca(self):
        """
        rda_cca分析
        :return:tnf
        """
        self.copy_collection_with_change("rda_cca")

    def regression(self):
        """
        regression回归分析
        :return:
        """
        self.copy_collection_with_change("regression")

    def roc(self):
        """
        roc分析
        :return:
        """
        self.copy_collection_with_change("roc")

    def sample_check(self):
        """
        样本检测
        :return:
        """
        self.copy_collection_with_change("sample_check")

    def specimen(self):
        """
        样本记录表
        :return:
        """
        specimen_dict = self.copy_collection_with_change("specimen")
        self.copy_datastat_details('specimen_detail', 'specimen_id', specimen_dict)

    def tax4fun(self):
        """
        tax4fun分析
        :return:
        """
        self.copy_collection_with_change("tax4fun")

    def two_corr_network(self):
        """
        双因素网络图分析
        :return:tnf
        """
        self.copy_collection_with_change("two_corr_network")

    def venn(self):
        """
        venn分析
        :return:
        """
        self.copy_collection_with_change("venn")

    def vif(self):
        """
        vif分析
        :return:
        """
        self.copy_collection_with_change("vif")

    def vpa(self):
        """
        PCA分析
        :return:tnf
        """
        self.copy_collection_with_change("vpa")

    def run(self):
        print "STRART COPYMEMBERID"
        self.copy_member_id()
        greenlet = Greenlet(self.alpha_diversity)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.alpha_diversity_diff)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.asv)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.asv_set)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.barpie)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.circos)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.cor_heatmap)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.cor_network)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.data_stat)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.dbrda)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.enterotype)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.env)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.env_group)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.funguild)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.hcluster)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.heatmap)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.lefse)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.mantel)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.multiple_group)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.two_group)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.two_sample)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.nmds)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.pan_core)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.pca)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.pcoa)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.permanova)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.personal_phylo_tree)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.phylo_tree)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.picrust2)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.randomForest)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.rank_abundance)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.rarefaction)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.rda_cca)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.regression)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.roc)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.sample_check)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.specimen)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.tax4fun)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.two_corr_network)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.venn)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.vif)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.vpa)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        gevent.joinall(self.all_greenlets)
        import socket
        reload(socket)


if __name__ == '__main__':
    start_time = time.time()
    copy_task = CopyDemo('tsg_33191', '188_5c26cc136502c ', '111112', '111113', '2')  # input params from task existed
    copy_task.run()
    end_time = time.time()
    print "total time: {}".format(end_time - start_time)
