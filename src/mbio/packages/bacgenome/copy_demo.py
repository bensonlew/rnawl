# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import json
import gevent
import datetime
import time  # 统计拷贝时间
from bson import ObjectId
from biocluster.api.database.base import Base
from gevent import Greenlet
from gevent.monkey import patch_all



class CopyDemo(Base):
    def __init__(self, old_task_id, new_task_id, new_project_sn, new_member_id, new_member_type):
        super(CopyDemo, self).__init__()
        self._project_type = 'bacgenome'
        self._old_task_id = old_task_id
        self._new_task_id = new_task_id
        self._new_project_sn = new_project_sn
        self._new_member_id = new_member_id
        self._new_member_type = new_member_type
        self.specimen_id_dict = {}
        self.all_greenlets = []


    def copy_member_id(self):
        """
        复制bac_task的数据
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
        finds = coll.find({'task_id': self._old_task_id})
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

    def datastat(self):
        datastat_id_dict = self.copy_collection_with_change("datastat")
        self.copy_main_details('datastat_gene', 'datastat_id', datastat_id_dict, join=False)
        self.copy_main_details('datastat_specimen', 'datastat_id', datastat_id_dict, join=False)

    def gbk(self):
        gbk_id_dict = self.copy_collection_with_change("gbk")
        self.copy_main_details('gbk_detail', 'gbk_id', gbk_id_dict, join=False)

    def assemble(self):
        ass_dict = self.copy_collection_with_change("assemble")
        self.copy_main_details('assemble_seq', 'assemble_id', ass_dict, join=False)

    def predict_gene(self):
        gene_predict_dict= self.copy_collection_with_change("gene_predict")
        self.copy_main_details('gene_predict_seq', 'predict_id', gene_predict_dict, join=False)

    def predict_rrna(self):
        rrna_dict = self.copy_collection_with_change("rrna_predict")
        self.copy_main_details('rrna_predict_detail', 'predict_id', rrna_dict, join=False)

    def predict_trna(self):
        trna_dict = self.copy_collection_with_change("trna_predict")
        self.copy_main_details('trna_predict_detail', 'predict_id', trna_dict, join=False)

    def predict_repeat(self):
        self.copy_collection_with_change("repeat_predict")

    def anno_nr(self):
        self.copy_collection_with_change("anno_nr")

    def anno_cog(self):
        self.copy_collection_with_change("anno_cog")

    def anno_kegg(self):
        self.copy_collection_with_change("anno_kegg")

    def anno_go(self):
        self.copy_collection_with_change("anno_go")

    def anno_pfam(self):
        self.copy_collection_with_change("anno_pfam")

    def anno_cazy(self):
        self.copy_collection_with_change("anno_cazy")

    def anno_card(self):
        self.copy_collection_with_change("anno_card")

    def anno_phi(self):
        self.copy_collection_with_change("anno_phi")

    def anno_tcdb(self):
        self.copy_collection_with_change("anno_tcdb")

    def anno_summary(self):
        summary_dict= self.copy_collection_with_change("anno_summary")
        self.copy_main_details('anno_summary_detail', 'summary_id', summary_dict, join=False)

    def anno_swissprot(self):
        self.copy_collection_with_change("anno_swissprot")

    def anno_vfdb(self):
        self.copy_collection_with_change("anno_vfdb")

    def anno_tmhmm(self):
        self.copy_collection_with_change("anno_tmhmm")

    def anno_secretory(self):
        self.copy_collection_with_change("anno_secretory")

    def anno_antismash(self):
        self.copy_collection_with_change("anno_antismash")

    def assess_gc(self):
        self.copy_collection_with_change("assess_gc")

    def assess_kmer(self):
        self.copy_collection_with_change("assess_kmer")

    def assess_size(self):
        self.copy_collection_with_change("assess_size")

    def island(self):
       island_dict = self.copy_collection_with_change("island")
       self.copy_main_details('island_detail', 'island_id', island_dict, join=False)

    def criprs(self):
        self.copy_collection_with_change("crispr")

    def prephage(self):
        pre_dict = self.copy_collection_with_change("prephage")
        self.copy_main_details('prephage_detail', 'prephage_id', pre_dict, join=False)

    def gene_graph(self):
        self.copy_collection_with_change("gene_graph")

    def promote(self):
        self.copy_collection_with_change("promote")

    def cgview(self):
        self.copy_collection_with_change("cgview")

    def circos(self):
        self.copy_collection_with_change("circos")
        circos_id = self.copy_collection_with_change("circos_table")
        self.copy_main_details('circos_table_detail', 'circos_table_id', circos_id, join=False)


    # v2 new add
    def regulator(self):
        self.copy_collection_with_change("anno_regulator")

    # v2 new add
    def methylation(self):
        self.copy_collection_with_change("methylation")

    def signalp(self):
        self.copy_collection_with_change("anno_signalp")

    def run(self):
        self.copy_member_id()
        greenlet = Greenlet(self.datastat)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.gbk)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.circos)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.predict_gene)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.anno_summary)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.assemble)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.predict_rrna)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.predict_trna)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.island)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.prephage)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        gevent.joinall(self.all_greenlets)
        self.gene_graph()
        self.predict_repeat()
        self.assess_gc()
        self.assess_kmer()
        self.assess_size()
        self.anno_nr()
        self.anno_kegg()
        self.anno_cog()
        self.anno_go()
        self.anno_cazy()
        self.anno_card()
        self.anno_tmhmm()
        self.anno_tcdb()
        self.anno_phi()
        self.anno_vfdb()
        self.anno_pfam()
        self.anno_swissprot()
        self.anno_secretory()
        self.promote()
        self.anno_antismash()
        self.criprs()
        self.cgview()
        ## add v2  zouguanqing
        self.methylation()
        self.regulator()
        self.signalp()

        import socket
        reload(socket)


if __name__ == '__main__':
    start_time = time.time()
    #copy_task = CopyDemo('tsg_26035_test2', 'demotest_metagenome41', '111112', '111113', '2')  # input params from task existed
    #copy_task.run()
    end_time = time.time()
    print "total time: {}".format(end_time - start_time)