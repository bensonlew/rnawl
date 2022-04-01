# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'@20191104

import json
import datetime
import gevent
import time  # 统计拷贝时间
from bson import ObjectId
from biocluster.api.database.base import Base
from gevent import Greenlet
from gevent.monkey import patch_all


class BacComparativeCopyDemo(Base):
    def __init__(self, old_task_id, new_task_id, new_project_sn, new_member_id, new_member_type):
        super(BacComparativeCopyDemo, self).__init__()
        self._project_type = 'bac_comparative'
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
        # ref_rna, 没有运行Greenlet?
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

    def aai(self):
        """
        aai分析
        :return:
        """
        self.copy_collection_with_change('aai')

    def ani(self):
        """
        ani分析
        :return:
        """
        self.copy_collection_with_change('ani')

    def card(self):
        """
        card
        :return:
        """
        self.copy_collection_with_change('anno_card')

    def cazy(self):
        """
        :return:
        """
        self.copy_collection_with_change('anno_cazy')

    def cog(self):
        """
        :return:
        """
        self.copy_collection_with_change('anno_cog')

    def kegg(self):
        """
        :return:
        """
        self.copy_collection_with_change('anno_kegg')

    def phi(self):
        """
        :return:
        """
        self.copy_collection_with_change('anno_phi')

    def secretory(self):
        """
        :return:
        """
        self.copy_collection_with_change('anno_secretory')

    def signalp(self):
        """
        :return:
        """
        self.copy_collection_with_change('anno_signalp')

    def summary(self):
        """
        :return:
        """
        self.copy_collection_with_change('anno_summary')

    def tcdb(self):
        """
        :return:
        """
        self.copy_collection_with_change('anno_tcdb')

    def tmhmm(self):
        """
        :return:
        """
        self.copy_collection_with_change('anno_tmhmm')

    def antismash(self):
        """
        :return:
        """
        self.copy_collection_with_change('antismash')

    def vfdb(self):
        """
        :return:
        """
        self.copy_collection_with_change('anno_vfdb')

    def cog_diff(self):
        """
        genome管理
        :return:
        """
        self.copy_collection_with_change('cog_diff')

    def core_gene(self):
        """
        :return:
        """
        self.copy_collection_with_change('core_gene')

    def cor_pan_comp(self):
        """
        :return:
        """
        self.copy_collection_with_change('cor_pan_comp')

    def core_pan_diff(self):
        """
        :return:
        """
        self.copy_collection_with_change('core_pan_diff')

    def correlation(self):
        """
        :return:
        """
        self.copy_collection_with_change('correlation')

    def data(self):
        """
        :return:
        """
        data_dict = self.copy_collection_with_change('data')
        self.copy_main_details('data_attributes', 'data_id', data_dict, join=False)
        self.copy_main_details('data_gene', 'data_id', data_dict, join=False)
        self.copy_main_details('data_specimen', 'data_id', data_dict, join=False)

    def gene_predict(self):
        """
        :return:
        """
        self.copy_collection_with_change('gene_predict')

    def gene_stat(self):
        """
        :return:
        """
        self.copy_collection_with_change('gene_stat')

    def hcluster(self):
        """
        :return:
        """
        self.copy_collection_with_change('hcluster')

    def island(self):
        """
        :return:
        """
        self.copy_collection_with_change('island')

    def kegg_diff(self):
        """
        :return:
        """
        self.copy_collection_with_change('kegg_diff')

    def kegg_enrichment(self):
        """
        :return:
        """
        self.copy_collection_with_change('kegg_enrichment')

    def nmds(self):
        """
        :return:
        """
        self.copy_collection_with_change('nmds')

    def pan(self):
        """
        :return:
        """
        self.copy_collection_with_change('pan')

    def pan_genome(self):
        """
        :return:
        """
        self.copy_collection_with_change('pan_genome')

    def pan_group(self):
        """
        :return:
        """
        self.copy_collection_with_change('pan_group')

    def pan_venn(self):
        """
        :return:
        """
        self.copy_collection_with_change('pan_venn')

    def pca(self):
        """
        :return:
        """
        self.copy_collection_with_change('pca')

    def pcoa(self):
        """
        :return:
        """
        self.copy_collection_with_change('pcoa')

    def prephage(self):
        """
        :return:
        """
        self.copy_collection_with_change('prephage')

    def specimen_group(self):
        """
        :return:
        """
        self.copy_collection_with_change('specimen_group')

    def synteny(self):
        """
        :return:
        """
        self.copy_collection_with_change('synteny')

    def tree(self):
        """
        :return:
        """
        self.copy_collection_with_change('tree')

    def run(self):
        print "STRART COPYMEMBERID"
        self.copy_member_id()
        greenlet = Greenlet(self.data)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.tree)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.synteny)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.specimen_group)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.prephage)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.pcoa)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.pca)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.nmds)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.pan_venn)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.pan_group)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.pan_genome)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.pan)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.cog)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.ani)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.aai)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.signalp)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.tcdb)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.phi)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.cazy)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.card)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.summary)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.kegg_enrichment)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.island)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.vfdb)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.kegg_diff)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.hcluster)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.tmhmm)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.correlation)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.antismash)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.cog_diff)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.kegg)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.secretory)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.core_gene)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.gene_predict)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.gene_stat)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.cor_pan_comp)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.core_pan_diff)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        gevent.joinall(self.all_greenlets)
        import socket
        reload(socket)


if __name__ == '__main__':
    start_time = time.time()
    copy_task = BacComparativeCopyDemo('tsg_33191', '188_5c26cc136502c ', '111112', '111113', '2')
    copy_task.run()
    end_time = time.time()
    print "total time: {}".format(end_time - start_time)
