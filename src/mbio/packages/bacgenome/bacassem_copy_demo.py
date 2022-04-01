# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'@20190103

import json
import datetime
import gevent
import time  # 统计拷贝时间
from bson import ObjectId
from biocluster.api.database.base import Base
from gevent import Greenlet
from gevent.monkey import patch_all


class BacassemCopyDemo(Base):
    def __init__(self, old_task_id, new_task_id, new_project_sn, new_member_id, new_member_type):
        super(BacassemCopyDemo, self).__init__()
        self._project_type = 'bac_assem'
        self._old_task_id = old_task_id
        self._new_task_id = new_task_id
        self._new_project_sn = new_project_sn
        self._new_member_id = new_member_id
        self._new_member_type = new_member_type
        self.all_greenlets = []
        self._exchange_dict = {
            "genome_id": {},
            "gap_id": {}
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

    def params_exchange(self, params_str):
        """
        专门用于params的数据ID替换, 注：metagenomic 需替换geneset_id,anno_id,group_id
        :return:
        """
        if params_str == None:
            return params_str
        try:
            params = json.loads(params_str)
        except Exception:
            print("WARNNING: 非json格式的params: {}".format(params_str))
            return params_str
        if not params:
            return params_str
        if 'genome' in params:
            if str(params['genome']) in self._exchange_dict["genome_id"].keys():
                params['genome'] = self._exchange_dict["genome_id"][str(params['genome'])]
            else:
                print "genome %s is not in dict: " % params['genome']
                print self._exchange_dict["genome_id"]
        return json.dumps(params, sort_keys=True, separators=(',', ':'))

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

    def sample_info(self):
        self.copy_collection_with_change("sample_info")

    def genome(self):
        coll = self.db["genome"]
        finds = coll.find({'task_id': self._old_task_id})
        news = []
        olds = []
        for i in finds:
            i['task_id'] = self._new_task_id
            olds.append(str(i.pop('_id')))
            i["gap_id"] = self.exchange_ObjectId("gap_id", i["gap_id"])
            news.append(i)
        if news:
            result = coll.insert_many(news)
            genome_id =  dict(zip(olds, [str(one) for one in result.inserted_ids]))
            self._exchange_dict["genome_id"] = genome_id

    def draft(self):
        self.copy_collection_with_change("draft")

    def draft_assess(self):
        self.copy_collection_with_change("draft_assess")

    def chr_assess(self):
        self.copy_collection_with_change("kmer_dist")
        self.copy_collection_with_change("gc_depth")
        self.copy_collection_with_change("blast_nt")
        self.copy_collection_with_change("complete_stat", change_positions=["genome_id"], update_sg_status=True)
        self.copy_collection_with_change("kmer_pca", update_sg_status=True)
        self.copy_collection_with_change("organism", update_sg_status=True)

    def gap_fill(self):
        gap_fill_id = self.copy_collection_with_change("gap_fill", update_sg_status=True)
        self._exchange_dict["gap"] = gap_fill_id

    def run(self):
        print "STRART COPYMEMBERID"
        self.copy_member_id()
        self.sample_info()
        self.draft()
        self.draft_assess()
        self.gap_fill()
        self.genome()
        greenlet = Greenlet(self.chr_assess)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        # greenlet = Greenlet(self.bin)
        # greenlet.start()
        # self.all_greenlets.append(greenlet)

        gevent.joinall(self.all_greenlets)
        import socket
        reload(socket)


if __name__ == '__main__':
    start_time = time.time()
    copy_task = BacassemCopyDemo('tsg_33191', '188_5c26cc136502c ', '111112', '111113', '2')  # input params from task existed
    copy_task.run()
    end_time = time.time()
    print "total time: {}".format(end_time - start_time)
