# -*- coding: utf-8 -*-
# __author__ = 'sheng.he'
# lastmodied: 20160921
import json
import gevent
import datetime
from bson import ObjectId
from gevent import Greenlet
#from gevent.monkey import patch_all
from biocluster.config import Config
from mainapp.libs.param_pack import group_detail_sort


class CopyMongo(object):
    """
    需要回滚操作。。。,暂时未提供
    """
    def __init__(self, old_task_id, new_task_id, new_project_sn, new_member_id, db='tsanger_ref_rna'):
        self.config = Config()
        self.db = self.config.mongo_client[db]
        self._old_task_id = old_task_id
        self._new_task_id = new_task_id
        self._new_project_sn = new_project_sn
        self._new_member_id = new_member_id
        self.specimen_id_dict = {}
        self.group_id_dict = {}
        # self.env_id_dict = {}
        # self.otu_id_dict = {}
        # self.alpha_diversity_id_dict = {}
        # self.newick_tree_id_dict = {}
        # self.specimen_distance_id_dict = {}
        self.all_greenlets = []
        self._exchange_dict = {  # 根据特定字段名称，进行特定的ID新旧替换
            'specimen_id': self.specimen_id_dict,
            'group_id': self.group_id_dict,
            'env_id': self.env_id_dict,
            'otu_id': self.otu_id_dict,
            'alpha_diversity_id': self.alpha_diversity_id_dict,
            'newick_tree_id': self.newick_tree_id_dict,
            'specimen_distance_id': self.specimen_distance_id_dict
            }

    def run(self):
        """
        运行执行复制特定ID数据的操作，如果有新的分析请参照下面的写法添加代码，不同分析表结构不同，所有需要手动添加。
        """
        #patch_all()
        self.copy_member_id()
        self.copy_sg_specimen()
        # self.copy_sg_specimen_group()
        # self.copy_collection_with_change('sg_specimen_sequence', change_positions=['specimen_id', ], join=False)
        # self.copy_collection_with_change('sg_specimen_step', change_positions=['specimen_id', ], join=False)
        # self.copy_sg_otu()
        # self.copy_collection_with_change('sg_otu_detail', change_positions=['otu_id', ], join=False)
        # self.copy_collection_with_change('sg_otu_detail_level', change_positions=['otu_id', ], join=False)
        # self.env_id_dict = self.copy_collection_with_change('sg_env')
        # self.copy_main_details('sg_env_detail', 'env_id', self.env_id_dict, others_position=['specimen_id'], join=True)
        # self.env_id_dict[None] = None
        # self.env_id_dict[''] = None
        # self._exchange_dict['env_id'] = self.env_id_dict
        # self.copy_main_details('sg_otu_specimen', 'otu_id', self.otu_id_dict, others_position=['specimen_id'], join=False)
        # self.copy_sg_newick_tree()
        # self.recopy_update_otu()
        greenlet = Greenlet(self.species_env_correlation)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.species_mantel_check)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.otu_venn)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.otu_pan_core)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.alpha_rarefaction_curve)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.beta_specimen_distance)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.beta_multi_analysis)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.beta_multi_anosim)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.species_difference_check)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.species_difference_lefse)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.network)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.randomforest)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.phylo_tree)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.hc_heatmap)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.sequence_info)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.corr_network)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.alpha_diversity)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.sixteen_function_predict)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.roc)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        gevent.joinall(self.all_greenlets)
        gevent.joinall(self.all_greenlets)
        import socket
        reload(socket)


    def _copy_collection_with_change(self, collection, change_positions=[], update_sg_status=False):
        """
        公共模块，一般用于导入主表数据，依靠task_id进行查询，修改change_positions提供的字段，相应修改ID为新的，同时更新params中的数据ID

        params collection: 主表名称
        params change_positions: 需要替换的ID,
            可用为specimen_id,group_id,env_id,otu_id,alpha_diversity_id,newick_tree_id,specimen_distance_id
        params update_sg_status: 更新 sg_status表
        """
        time_start = datetime.datetime.now()
        coll = self.db[collection]
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
            if 'params' in i:
                i['params'] = self.params_exchange(i['params'])
            news.append(i)
        if news:
            result = coll.insert_many(news)
            if update_sg_status:
                self.insert_new_status(collection, news, result.inserted_ids)
            time_end = datetime.datetime.now()
            run_time = (time_end - time_start).seconds
            print "{}复制运行时间: {}s".format(collection, run_time)
            return dict(zip(olds, [str(one) for one in result.inserted_ids]))
        else:
            time_end = datetime.datetime.now()
            run_time = (time_end - time_start).seconds
            print "{}复制运行时间: {}s".format(collection, run_time)
            return {}

    def copy_collection_with_change(self, collection, change_positions=[], update_sg_status=False, join=True):
        greenlet = Greenlet(self._copy_collection_with_change, collection, change_positions, update_sg_status)
        greenlet.start()
        if join is True:
            greenlet.join()
            return greenlet.value
        self.all_greenlets.append(greenlet)
        return greenlet

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
        coll = self.db[collection]
        for old, new in change_dict.items():
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

    def insert_new_status(self, collection, main_docs, ids):
        """
        导入mongo表sg_status数据信息
        """
        coll = self.db.sg_status
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

    def copy_sg_newick_tree(self):
        """
        复制sg_newick_tree的数据，返回新旧ID字典
        """
        finds = self.db.sg_newick_tree.find({'task_id': self._old_task_id})
        olds = []
        news = []
        for i in finds:
            olds.append(str(i.pop('_id')))
            i['task_id'] = self._new_task_id
            if 'project_sn' in i:
                i['project_sn'] = self._new_project_sn
            if 'params' in i:
                i['params'] = self.params_exchange(i['params'])
            if 'table_id' in i:
                i['table_id'] = ObjectId(self.otu_id_dict[str(i['table_id'])])
            news.append(i)
        if news:
            result = self.db.sg_newick_tree.insert_many(news)
            self.newick_tree_id_dict = dict(zip(olds, [str(one) for one in result.inserted_ids]))
            news = zip(result.inserted_ids, news)
            news = [i for i in news if i[1]["tree_type"] == "cluster"]  # 只用cluster是beta分析中的结果，需要进行更新sg_status
            docs = [i[1] for i in news]
            ids = [i[0] for i in news]
            self.insert_new_status('sg_newick_tree', docs, ids)
        else:
            self.newick_tree_id_dict = {}
        self.newick_tree_id_dict[''] = None
        self.newick_tree_id_dict[None] = None
        self._exchange_dict['newick_tree_id'] = self.newick_tree_id_dict
        return self.newick_tree_id_dict

    def recopy_update_otu(self):
        """
        otu表与newick tree有相互重复引用ID，在导入otu表和newicktree数据后，需要返回重新更新otu表中的newick_tree_id(newick_id)数据
        """
        for i in self.otu_id_dict.values():
            find = self.db.sg_otu.find_one({'_id': ObjectId(i)})
            if 'newick_id' in find:
                self.db.sg_otu.update_one({'_id': ObjectId(i)}, {'$set': {'newick_id': self.exchange_ObjectId('newick_tree_id', find['newick_id'])}})
            elif 'newick_tree_id' in find:
                self.db.sg_otu.update_one({'_id': ObjectId(i)}, {'$set': {'newick_tree_id': self.exchange_ObjectId('newick_tree_id', find['newick_id'])}})
        pass

    def copy_sg_env_detail(self):  # 未使用
        for old, new in self.env_id_dict.items():
            finds = self.db.sg_env_detail.find({'env_id': ObjectId(old)})
            news = []
            for i in finds:
                i.pop('_id')
                i['env_id'] = ObjectId(new)
                i['specimen_id'] = ObjectId(self.specimen_id_dict[str(i['specimen_id'])])
                news.append(i)
            if news:
                self.db.sg_env_detail.insert_many(news)
            else:
                print 'WARNING: 环境因子主表:{}没有detail表信息，请注意数据合理性'.format(old)

    def copy_sg_otu(self):
        """
        复制otu主表
        """
        finds = self.db.sg_otu.find({"task_id": self._old_task_id})
        news = []
        old_otu_ids = []
        for i in finds:
            i['task_id'] = self._new_task_id
            i['project_sn'] = self._new_project_sn
            old_otu_ids.append(str(i.pop('_id')))
            news.append(i)
        if news:
            result = self.db.sg_otu.insert_many(news)
        else:
            raise Exception('不存在任何OTU表，请确认数据的完整性')
        self.otu_id_dict = dict(zip(old_otu_ids, [str(one) for one in result.inserted_ids]))
        self._exchange_dict['otu_id'] = self.otu_id_dict
        update_sg_status_docs = []
        for i in result.inserted_ids:
            find = self.db.sg_otu.find_one({'_id': i})
            update_dict = {}
            if 'from_id' in find:
                update_dict['from_id'] = self.otu_id_dict[find['from_id']]
            if ('params' in find) and find['params']:
                params = json.loads(find['params'])
                if 'group_detail' in params:
                    for one_group in params['group_detail']:
                        params['group_detail'][one_group] = [self.specimen_id_dict[one_sp] for one_sp in params['group_detail'][one_group]]
                    params['group_detail'] = group_detail_sort(params['group_detail'])
                if 'group_id' in params:
                    params['group_id'] = self.group_id_dict[params['group_id']]
                if 'otu_id' in params:
                    params['otu_id'] = self.otu_id_dict[params['otu_id']]
                update_dict['params'] = json.dumps(params, sort_keys=True, separators=(',', ':'))
            find.update(update_dict)
            update_sg_status_docs.append(find)
            self.db.sg_otu.update_one({'_id': i}, {'$set': update_dict})
        update_sg_status_docs = [i for i in update_sg_status_docs if i['type'] in ['otu_statistic', "otu_statistic,otu_venn,otu_group_analysis", 'otu_group_analyse']]
        ids = [i['_id'] for i in update_sg_status_docs]
        self.insert_new_status('sg_otu', update_sg_status_docs, ids)
        return self.otu_id_dict

    def copy_sg_specimen_sequence(self):  # 未使用
        finds = self.db.sg_specimen_sequence.find({'task_id': self._old_task_id})
        news = []
        for i in finds:
            i['task_id'] = self._new_task_id
            i['project_sn'] = self._new_project_sn
            i.pop('_id')
            i['specimen_id'] = ObjectId(self.specimen_id_dict[str(i['specimen_id'])])
            news.append(i)
        if news:
            result = self.db.sg_specimen_sequence.insert_many(news)
            return result.inserted_ids
        else:
            raise Exception('不存在任何样本序列， 请检查数据完整性')

    def copy_sg_specimen_group(self):
        """
        复制group表
        """
        finds = self.db.sg_specimen_group.find({"task_id": self._old_task_id})
        news = []
        old_group_ids = []
        for i in finds:
            i['task_id'] = self._new_task_id
            i['project_sn'] = self._new_project_sn
            old_group_ids.append(str(i.pop('_id')))
            for one in i['specimen_names']:
                for sp in one.copy():
                    one[self.specimen_id_dict[sp]] = one[sp]
                    one.pop(sp)
            news.append(i)
        if news:
            result = self.db.sg_specimen_group.insert_many(news)
            self.group_id_dict = dict(zip(old_group_ids, [str(one) for one in result.inserted_ids]))
        else:
            self.group_id_dict = {}
        self.group_id_dict['all'] = 'all'  # 特殊ID
        self.group_id_dict[None] = None  # 特殊ID
        self.group_id_dict[''] = None  # 特殊ID
        self._exchange_dict['group_id'] = self.group_id_dict
        return self.group_id_dict

    def copy_sg_specimen(self):
        """
        复制样本表
        """
        finds = self.db.sg_specimen.find({"task_id": self._old_task_id})
        news = []
        old_specimen_ids = []
        for i in finds:
            old_specimen_ids.append(str(i.pop('_id')))
            i['task_id'] = self._new_task_id
            i['project_sn'] = self._new_project_sn
            news.append(i)
        if news:
            result = self.db.sg_specimen.insert_many(news)
        else:
            raise Exception('没有任何样本信息，请核对任务结果是否完整')
        self.specimen_id_dict = dict(zip(old_specimen_ids, [str(one) for one in result.inserted_ids]))
        self._exchange_dict['specimen_id'] = self.specimen_id_dict
        return self.specimen_id_dict

    @staticmethod
    def pop_id(doc):  # 未使用
        """
        """
        pop_one = doc.pop('_id')
        return str(pop_one)

    def exchange_ObjectId(self, key, thisObjectId):
        """
        用于替换id，key是该ID的字段名，thisObjectId是旧的ID(ObjectId类型)
        """
        if isinstance(thisObjectId, ObjectId):
            return ObjectId(self._exchange_dict[key][str(thisObjectId)])
        else:
            return self._exchange_dict[key][thisObjectId]  # 不是ObjectId时直接返回也是字符串

    def exchange_strId(self, key, strId):  # 未使用
        """
        """
        return self._exchange_dict[key][strId]

    def copy_member_id(self):
        """
        复制sg_task的数据
        """
        coll = self.db.sg_task
        find = coll.find_one({'task_id': self._old_task_id})
        if not find:
            raise Exception('运行错误：找不到demo任务相关信息')
        find['task_id'] = self._new_task_id
        find['member_id'] = self._new_member_id
        find.pop('_id')
        find['project_sn'] = self._new_project_sn
        coll.insert_one(find)

    def params_exchange(self, params_str):
        """
        专门用于params的数据ID替换
        """
        try:
            params = json.loads(params_str)
        except Exception:
            print("WRANNING：非json格式的params：{}".format(params_str))
            return params_str
        if not params:
            return None
        if 'group_detail' in params:
            for one_group in params['group_detail']:
                params['group_detail'][one_group] = [self.specimen_id_dict[one_sp] for one_sp in params['group_detail'][one_group]]
            params['group_detail'] = group_detail_sort(params['group_detail'])
            if 'second_group_detail' in params:
                if params['second_group_detail']:
                    for one_group in params['second_group_detail']:
                        params['second_group_detail'][one_group] = [self.specimen_id_dict[one_sp] for one_sp in params['second_group_detail'][one_group]]
                    params['second_group_detail'] = group_detail_sort(params['second_group_detail'])
                    if 'second_group_id' in params:
                        params['second_group_id'] = self.group_id_dict[params['second_group_id']]
        if 'group_id' in params:
            params['group_id'] = self.group_id_dict[params['group_id']]
        if 'otu_id' in params:
            params['otu_id'] = self.otu_id_dict[params['otu_id']]
        if 'env_id' in params:
            params['env_id'] = self.env_id_dict[params['env_id']]
        if 'alpha_diversity_id' in params:
            params['alpha_diversity_id'] = self.alpha_diversity_id_dict[params['alpha_diversity_id']]
        return json.dumps(params, sort_keys=True, separators=(',', ':'))

if __name__ == '__main__':
    copy_task = CopyMongo('tsanger_7186', 'tsg_7186_22', '10004002_1', 'shenghe_test')
    copy_task.run()

