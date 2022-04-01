# !/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "hongdongxuan"
# lastmodified = "guhaidong"
# lastmodifieddate = "20171109"

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
from mainapp.libs.param_pack import group_detail_sort, filter_json_sort, sub_group_detail_sort
from socket import error as SocketError
import urllib2
import urllib
import httplib
import json
from bson.objectid import ObjectId
import gevent
import time
import hashlib
from gevent.event import AsyncResult
import pymongo
import copy
import random
import datetime
import re
from types import StringTypes


class PipeSubmitAgent(Agent):
    """
    用于submit子分析的接口,针对的是所有的分析，后面还会有医口与农口
    version v1.0
    author: hongdongxuan
    last_modify: 2017.03.09
    """

    def __init__(self, parent):
        super(PipeSubmitAgent, self).__init__(parent)
        options = [
            {"name": "data", "type": "string"},
            {"name": "pipe_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "batch_task_id", "type": "string"}
        ]
        self.add_option(options)
        self.step.add_steps("piple_submit")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.piple_submit.start()
        self.step.update()

    def stepfinish(self):
        self.step.piple_submit.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("data"):
            raise OptionError("必须输入接口传进来的data值", code="32705901")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(PipeSubmitAgent, self).end()


class PipeSubmitTool(Tool):

    def __init__(self, config):
        super(PipeSubmitTool, self).__init__(config)
        self._version = '1.0.1'
        # self._config = Config()
        # self._client = self._config.mongo_client
        # self.db = self._client[self._config.MONGODB]
        self._project_type = "meta"
        self.client = Config().get_mongo_client(mtype=self._project_type)
        self.db = self.client[Config().get_mongo_dbname(self._project_type)]
        self.max_instant_running = 18
        self.max_submit_running = 40
        self.count_submit_running = 0
        self.count_instant_running = 0

    def rely_error(self, ana, error_rely):
        ana.is_end = True
        self.logger.error('依赖对象运行错误: {}'.format(error_rely.api))
        if ana.instant:
            self.count_instant_running += 1
        else:
            self.count_submit_running += 1
        self.one_end(ana)
        # raise Exception('ERROR: 依赖对象运行错误: {}'.format(error_rely.api))

    def get_origin_task_id(self):
        split_id = self.sheet.id.split('_')
        split_id.pop()
        split_id.pop()
        self.task_id = '_'.join(split_id)
        return self.task_id

    def update_mongo_ends_count(self, ana):

        # if ana.instant:
        #     self.logger.info("api: {} END_COUNT+1".format(ana.api))
        #     self.db['sg_pipe_batch'].find_one_and_update({'_id': ObjectId(self.option('pipe_id'))},
        #                                                  {'$inc': {"ends_count": 1}})
        # if ana._params_check_end or not ana.success:
        if ana._params['submit_location'] == "otu_pan_core":  # pancore 分析默认两个主表，一次结束 加 2
            inc = 2
        else:
            inc = 1
        self.logger.info("api: {} END_COUNT+{}".format(ana.api, inc))
        # self.db['sg_pipe_batch'].find_one_and_update({'_id': ObjectId(self.option('pipe_id'))},
        #                                              {'$inc': {"ends_count": inc}})
        self.db['sg_pipe_batch'].update({'_id': ObjectId(self.option('pipe_id'))}, {'$inc': {"ends_count": inc}})


    def one_end(self, ana):
        if ana._params['submit_location'] == "otu_pan_core":  # pancore 分析默认两个主表，一次结束 加 2
            self.count_ends += 2
        else:
            self.count_ends += 1
        if not ana.instant:
            self.logger.info('投递任务投递成功:{} 任务参数检查获取结果: {}'.format(ana.api, ana._params_check_end))
        if ana._params_check_end or not ana.success:  # 用于当参数都相同的时候没有进度条
            self.update_mongo_ends_count(ana)
        if ana.instant:
            self.count_instant_running -= 1
            self.logger.info("当前运行中的即时任务数为: {}".format(
                self.count_instant_running))
        else:
            self.count_submit_running -= 1
            self.logger.info("当前运行中的投递任务数为: {}, which analysis:{}".format(
                self.count_submit_running, ana._params['submit_location']))
        self.logger.info("END COUNT: {}".format(self.count_ends))

        # if not ana.success and self.count_ends == self.all_count:
        #     self.final_end(self.all_count)

        if self.count_ends == self.all_count:
            self.final_end(self.all_count)   # 最后一个分析要重置下all_count
            self.all_end.set()
        pass

    def final_end(self, all_count):
        """
        用于最后一个分析因为参数不合适，计算失败, 不投递出去计算，这样导致进度条更新不了
        :return:
        """
        self.logger.info("all_count:{}".format(all_count))
        real_all_count = self.db['sg_pipe_detail'].find({'pipe_batch_id': ObjectId(self.option('pipe_id')),
                                                         'status': {'$in': ['end', 'start', "failed"]}}).count()
        self.logger.info("real_all_count:{}".format(real_all_count))
        failed = self.db['sg_pipe_detail'].find({'pipe_batch_id': ObjectId(self.option('pipe_id')),
                                                 "submit_location": {'$in': ["alpha_diversity_index", "otu_pan_core"]},
                                                 'status': "failed"}).count()
        self.logger.info("第一次查询多样性指数与pan_core失败的个数:{}".format(failed))
        real_all_count += failed
        if real_all_count < all_count - 1:
            # self.db['sg_pipe_batch'].find_one_and_update({'_id': ObjectId(self.option('pipe_id'))},
            #                                              {'$set': {"all_count": real_all_count}}, upsert=True)
            self.db['sg_pipe_batch'].update({'_id': ObjectId(self.option('pipe_id'))},
                                                         {'$set': {"all_count": real_all_count}}, upsert=True)
            all_counts = real_all_count
        else:
            all_counts = all_count - 1  # 这里没有考虑将抽平放进去
        for i in xrange(10):
            ends_counts = self.db['sg_pipe_detail'].find({'pipe_batch_id': ObjectId(self.option('pipe_id')),
                                                          'status': {'$in': ['end', "failed"]}}).count()
            self.logger.info("第{}次查询ends_count值: {}".format(i + 1, ends_counts))
            failed_ = self.db['sg_pipe_detail'].find({'pipe_batch_id': ObjectId(self.option('pipe_id')),
                                                     "submit_location": {
                                                         '$in': ["alpha_diversity_index", "otu_pan_core"]},
                                                      'status': "failed"}).count()
            self.logger.info("第二次查询多样性指数与pan_core失败的个数:{}".format(failed_))
            ends_counts += failed_
            if ends_counts >= all_counts:
                # self.db['sg_pipe_batch'].find_one_and_update({'_id': ObjectId(self.option('pipe_id'))},
                #                                              {'$set': {"ends_count": ends_counts}}, upsert=True)
                self.db['sg_pipe_batch'].update({'_id': ObjectId(self.option('pipe_id'))},
                                                             {'$set': {"ends_count": ends_counts}}, upsert=True)
                self.logger.info("Tool中更新表格已经完成了，last_end_count:{}".format(ends_counts))
                break
            gevent.sleep(10)

    def get_params(self, config_name):
        """
        获取参数, 需要处理的参数请先特殊处理
        """
        config = self.analysis_params[config_name]
        params = {}
        for i in config['main']:
            if i in self.web_data:
                params[i] = self.web_data[i]
            else:
                self.logger.info("没有提供参数: {}".format(i))
        sub_params = self.web_data['sub_analysis'][config_name]
        for i in config['others']:
            if i == "second_group_detail" and sub_params[i]:
                params[i] = group_detail_sort(sub_params[i])  # lefse接口中使用的是group_detail_sort
            else:
                params[i] = sub_params[i]
        api = '/' + '/'.join(sub_params['api'].split('|'))
        instant = config['instant']
        collection_name = config['collection_name']
        return (api, instant, collection_name, params)

    def format_special_params(self):
        """
        预处理 参数， 如传入的filter_json和group_detail
        """
        pass

    def get_class(self, name):
        class_name = ''.join([i.capitalize() for i in name.split('_')])
        if class_name in globals():
            return globals()[class_name]
        else:
            self.logger.error("没有找到相应的类模块， 请添加：{}".format(class_name))
            self.set_error("模块加载错误", code="32705902")

    def get_otu_subsample_params(self, group_detail):
        """
        """
        params = {
            'otu_id': self.web_data['otu_id'],
            'filter_json': filter_json_sort(self.web_data['filter_json']),
            'group_id': 'all',
            'group_detail': group_detail_sort(group_detail),
            'submit_location': 'otu_statistic',
            'task_type': 'reportTask',
            'size': self.web_data['size']
        }
        return '/meta/otu_subsample', False, 'sg_otu', params

    def signature(self):
        timestamp = str(int(time.time()))
        nonce = str(random.randint(1000, 10000))
        web_key = self.mysql_client_key
        sha1 = hashlib.sha1()
        key_list = [web_key, nonce, timestamp]
        key_list.sort()
        map(sha1.update, key_list)
        hashkey = sha1.hexdigest()
        signature = {
            "client": self.task_client,
            "nonce": nonce,
            "timestamp": timestamp,
            "signature": hashkey
        }
        return urllib.urlencode(signature)

    def pipe_main_mongo_insert(self, level, group_info):
        insert_data = {
            'pipe_batch_id': ObjectId(self.option('pipe_id')),
            'level_id': int(level),
            'group_id': "all" if str(group_info['group_id']) == "all" else ObjectId(group_info['group_id']),
            'task_id': self.option('task_id')
        }
        inserted_id = self.db['sg_pipe_main'].insert_one(insert_data).inserted_id
        if int(level) == self.min_level and str(group_info['group_id']) == self.first_group_id:
            self.get_picture_id(pipe_main_id=inserted_id, main_table_id=self.option('pipe_id'))
        return inserted_id
        pass

    def get_picture_id(self, pipe_main_id, main_table_id):
        """
        批量导图片时使用， 目前选取的是最低分类水平，分组方案中all后面的第一个分组
        :param pipe_main_id:
        :param main_table_id:
        :return:
        """
        if pipe_main_id != 0 and not isinstance(pipe_main_id, ObjectId):
            if isinstance(pipe_main_id, StringTypes):
                pipe_main_id = ObjectId(pipe_main_id)
            else:
                self.set_error("pipe_main_id必须为ObjectId对象或其对应的字符串!", code="32705903")
        collection = self.db["sg_pipe_batch"]
        data = {
            "pipe_main_id": pipe_main_id
        }
        try:
            collection.update({"_id": ObjectId(main_table_id)}, {'$set': data}, upsert=False)
        except:
            self.set_error("sg_pipe_batch中pipe_main_id更新失败，请检查！", code="32705904")

    def run_webapitest(self):
        """
        进行参数判断后投递所有的接口，如果后面要添加新的分析，主要要添加以下几个地方的内容：
        1）analysis_params中保存的是每个分析的submit_location命名的字典，instant用于区分投递或者即时的；waits中包含的是分析之间
        的依赖关系，是基于哪些分析的结果进行下一步计算的；main保存的是主要的参数如："env_id", "env_labs"等，other中保存的是每个
        分析中特有的一些参数；collection_name是这个子分析对应的mongo中的主表名
        2）定义每个分析的类，如：CorrNetworkAnalyse（该分析的submit_location驼峰命名）
        3）一键化的所有的分析，都是以submit_location来区分，所以在对接的时候要着重考虑与前端的submit_location是不是一致的。
        4)所有任务都被设置为投递任务，但是，被依赖的任务必须设置为即时的。
        :return:
        """
        self.analysis_params = {
            "randomforest_analyse": {"instant": False, "waits": ["otu_subsample"], "main": [], "others": ["task_type", "ntree_id", "submit_location"], "collection_name": "sg_randomforest"},
            "sixteens_prediction": {"instant": False, "waits": ["otu_subsample"], "main": [], "others": ["group_method", "task_type", "submit_location"], "collection_name": "sg_16s"},
            "species_lefse_analyse": {"instant": False, "waits": ["otu_subsample"], "main": [], "others": ["task_type", "second_group_id", "lda_filter", "second_group_detail", "submit_location", "start_level", "strict", "end_level"], "collection_name": "sg_species_difference_lefse"},
            "alpha_rarefaction_curve": {"instant": False, "waits": ["otu_subsample"], "main": [], "others": ["index_type", "freq", "submit_location", "task_type"], "collection_name": "sg_alpha_rarefaction_curve"},
            "otunetwork_analyse": {"instant": False, "waits": ["otu_subsample"], "main": [], "others": ["task_type", "add_Algorithm", "submit_location"], "collection_name": "sg_network"},
            "roc_analyse": {"instant": False, "waits": ["otu_subsample"], "main": [], "others": ["task_type", "top_n_id", "method_type", "submit_location"], "collection_name": "sg_roc"},
            "alpha_diversity_index": {"instant": False, "waits": ["otu_subsample"], "main": [], "others": ["index_type", "submit_location", "task_type"], "collection_name": "sg_alpha_diversity"},
            "alpha_ttest": {"instant": False, "waits": ["alpha_diversity_index", "otu_subsample"], "main": [], "others": ["task_type", "submit_location", "test_method"], "collection_name": "sg_alpha_ttest"},
            "beta_multi_analysis_plsda": {"instant": False, "waits": ["otu_subsample"], "main": [], "others": ["analysis_type", "task_type", "submit_location"], "collection_name": "sg_beta_multi_analysis"},
            "beta_sample_distance_hcluster_tree": {"instant": False, "waits": ["otu_subsample"], "main": [], "others": ["distance_algorithm", "hcluster_method", "task_type", "submit_location"], "collection_name": "sg_newick_tree"},
            "beta_multi_analysis_pearson_correlation": {"instant": False, "waits": ["otu_subsample"], "main": ["env_id", "env_labs"], "others": ["task_type", "species_cluster", "submit_location", "method", "top_species", "env_cluster"], "collection_name": "sg_species_env_correlation"},
            "beta_multi_analysis_rda_cca": {"instant": False, "waits": ["otu_subsample"], "main": ["env_id", "env_labs"], "others": ["analysis_type", "task_type", "submit_location"], "collection_name": "sg_beta_multi_analysis"},
            "hc_heatmap": {"instant": False, "waits": ["otu_subsample"], "main": [], "others": ["task_type", "add_Algorithm", "submit_location", "sample_method", "species_number", "method"], "collection_name": "sg_hc_heatmap"},
            "beta_multi_analysis_anosim": {"instant": False, "waits": ["otu_subsample"], "main": [], "others": ["distance_algorithm", "permutations", "submit_location", "task_type"], "collection_name": "sg_beta_multi_anosim"},
            "species_difference_multiple": {"instant": False, "waits": ["otu_subsample"], "main": [], "others": ["task_type", "correction", "methor", "submit_location", "coverage", "test"], "collection_name": "sg_species_difference_check"},
            "beta_multi_analysis_results": {"instant": False, "waits": ["otu_subsample"], "main": ["env_id", "env_labs"], "others": ["otu_method", "env_method", "submit_location", "task_type"], "collection_name": "sg_species_mantel_check"},
            "beta_multi_analysis_pcoa": {"instant": False, "waits": ["otu_subsample"], "main": [], "others": ["analysis_type", "distance_algorithm", "submit_location", "task_type"], "collection_name": "sg_beta_multi_analysis"},
            "beta_multi_analysis_dbrda": {"instant": False, "waits": ["otu_subsample"], "main": ["env_id", "env_labs"], "others": ["analysis_type", "distance_algorithm", "submit_location", "task_type"], "collection_name": "sg_beta_multi_analysis"},
            "species_difference_two_group": {"instant": False, "waits": ["otu_subsample"], "main": [], "others": ["ci", "task_type", "correction", "methor", "submit_location", "coverage", "test", "type"], "collection_name": "sg_species_difference_check"},
            "beta_multi_analysis_nmds": {"instant": False, "waits": ["otu_subsample"], "main": [], "others": ["analysis_type", "distance_algorithm", "submit_location", "task_type"], "collection_name": "sg_beta_multi_analysis"},
            "otu_group_analyse": {"instant": False, "waits": ["otu_subsample"], "main": [], "others": ["task_type", "submit_location"], "collection_name": "sg_otu"},
            "otu_venn": {"instant": False, "waits": ["otu_subsample"], "main": [], "others": ["task_type", "submit_location"], "collection_name": "sg_otu_venn"},
            "beta_multi_analysis_pca": {"instant": False, "waits": ["otu_subsample"], "main": ["env_id", "env_labs"], "others": ["analysis_type", "task_type", "submit_location"], "collection_name": "sg_beta_multi_analysis"},
            "corr_network_analyse": {"instant": False, "waits": ["otu_subsample"], "main": [], "others": ["abundance", "coefficient", "ratio_method", "lable", "submit_location", "task_type"], "collection_name": "sg_corr_network"},
            "otu_pan_core": {"instant": False, "waits": ["otu_subsample"], "main": [], "others": ["task_type", "submit_location"], "collection_name": "sg_otu_pan_core"},
            "plot_tree": {"instant": False, "waits": ["otu_subsample"], "main": [], "others": ["task_type", "color_level_id", "submit_location", "topn"], "collection_name": "sg_phylo_tree"},
            "enterotyping": {"instant": False, "waits": ["otu_subsample"], "main": [], "others": ["task_type", "submit_location"], "collection_name": "sg_enterotyping"},
        }
        self.web_data = json.loads(self.option('data'))
        self.web_data['sub_analysis'] = json.loads(
            self.web_data['sub_analysis'])
        self.sub_analysis_len = len(self.web_data['sub_analysis'])
        self.otu_id = self.web_data['otu_id']
        self.group_infos = json.loads(self.web_data['group_info'])
        self.levels = [int(i) for i in self.web_data[
            'level_id'].strip().split(',')]
        self.min_level = int(sorted(self.levels)[-1])
        self.first_group_id = str(self.web_data['first_group_id'])
        self.count_ends = 0
        self.task_client = self.web_data['client']
        # self.mysql_client_key = worker_client().get_key(self.task_client)
        # self.mysql_client_key = worker_client().add_task({})
        self.mysql_client_key = 'mykey'
        # monkey.patch_socket(aggressive=False, dns=False)
        # monkey.patch_ssl()
        self.signature = self.signature()
        self.task_id = self.option("task_id")
        self.url = "http://bcl.i-sanger.com" if self.task_client == "client01" else "http://bcl.tsanger.com"
        self.all = {}
        self.all_count = 0
        sixteens_prediction_flag = False  # 16s功能预测分析特殊性，没有分类水平参数
        lefse_flag = False  # lefse特殊性，没有分类水平参数
        api, instant, collection_name, params = self.get_otu_subsample_params(self.web_data['group_detail'])
        otu_subsample = SubmitOtuSubsample(self, collection_name, params, api, instant)
        pipe_count = 0
        for level in self.levels:
            self.all[level] = {}
            for group_info in self.group_infos:
                pipe_main_id = self.pipe_main_mongo_insert(level, group_info)
                pipe = {}
                for i in self.web_data['sub_analysis']:
                    if i == 'sixteens_prediction' and sixteens_prediction_flag:
                        continue
                    elif i == "species_lefse_analyse" and lefse_flag:
                        continue
                    api, instant, collection_name, params = self.get_params(i)
                    params['group_id'] = group_info['group_id']
                    params['group_detail'] = group_detail_sort(group_info['group_detail'])
                    params['level_id'] = level
                    pipe[i] = self.get_class(i)(
                        self, collection_name, params, api, instant, pipe_main_id, pipe_count, self.min_level)
                pipe['otu_subsample'] = otu_subsample
                for analysis, submit in pipe.iteritems():
                    # self.logger.info("等待5s")
                    # time.sleep(5)
                    if analysis == 'otu_subsample':
                        continue
                    waits = [pipe[i]
                             for i in self.analysis_params[analysis]["waits"]]
                    submit.start(waits, timeout=15000)
                self.all[level][group_info["group_id"]] = pipe
                self.all_count += (len(pipe) - 1)
                pipe_count += 1
            sixteens_prediction_flag = True
            lefse_flag = True
        self.all_count += 1  # otu subsample 任务
        self.update_all_count()
        otu_subsample.start([], timeout=8000)
        self.all_end = AsyncResult()
        try:
            self.logger.info("所有任务计数: {}".format(self.all_count))
            self.all_end.get()
        except:
            self.logger.info('所有任务已经投递结束，计数:{}，总数:{}'.format(self.count_ends, self.all_count))
            for i in self.all:
                for one in self.all[i].values():
                    for key in one.values():
                        if not key.is_end:
                            self.logger.info("没有结束的submit: {}, api: {}, params: {}".format(key, key.api, key._params))
        self.end()

    def update_all_count(self):
        # self.db['sg_pipe_batch'].find_one_and_update({'_id': ObjectId(self.option('pipe_id'))},
        #                                              {'$set': {"all_count": self.all_count - 1}})
        self.db['sg_pipe_batch'].update({'_id': ObjectId(self.option('pipe_id'))},
                                                     {'$set': {"all_count": self.all_count - 1}})

    def run(self):
        super(PipeSubmitTool, self).run()
        self.run_webapitest()


class Submit(object):
    """投递对象"""

    def __init__(self, bind_object, collection, params, api, instant, pipe_main_id=None, pipe_count=0, min_level=0):
        """
        :params bind_object:
        :params collection:
        :params params:
        :params api:
        :params instant:
        """
        self.workflow_id = ''  # workflow的ID
        self.task_id = bind_object.task_id
        self.pipe_main_id = pipe_main_id
        self.api = api  # 接口url
        self.main_table_id = ''  # 主表_id
        self._params = params  # 参数
        self.success = False  # 是否成功
        self.is_end = False  # 是否结束
        self.error_info = None  # 错误信息
        self.mongo_collection = collection  # mongo表名称
        self.instant = instant  # 即时任务
        self._end_event = AsyncResult()  # end事件
        self.db = bind_object.db  # mongodb对象
        self.bind_object = bind_object  # tool对象
        self.out_params = {}  # 输出参数，单一个分析需要给其他分析提供参数时提供
        self.result = {}  # 返回结果
        self._params_check_end = False  # 参数检查完成的任务
        self.pipe_count = pipe_count  #
        self.min_level = min_level  # 所有level的数组中分类水平最低的

    def params_pack(self, dict_params):
        """
        参数打包，用于比对mongo数据库
        """
        return json.dumps(dict_params, sort_keys=True, separators=(',', ':'))

    def get_workflow_id(self):
        """
        获取workflow_id
        """
        # self.workflow_id = self.result['workflow_id']
        if self.workflow_id:
            return self.workflow_id
        result = self.db.workflowid2analysisid.find_one(
            {"main_id": ObjectId(self.main_table_id)}, {'workflow_id': 1})
        if result:
            self.workflow_id = result['workflow_id']
            return self.workflow_id
        else:
            self.bind_object.logger.error(
                "没有通过主表ID: {} 在workflowid2analysisid找到任务ID信息".format(self.main_table_id))
            self.set_error("没有找到任务ID信息", code="32705905")

    def start(self, waits, timeout=10000):
        self.waits = waits
        gevent.spawn(self._submit, waits, timeout=timeout)

    def end_fire(self):
        """
        分析结束后的后续工作
        """
        if "success" in self.result and self.result['success']:
            self.success = True
            self.set_out_params()
        self.is_end = True
        self.bind_object.one_end(self)

    def _submit(self, waits, timeout):
        """投递任务"""
        for i in waits:
            i.end_event.get(timeout=timeout)
            if not i.success:
                self._end_event.set()
                self.bind_object.rely_error(self, i)
                return
        self.run_permission()
        self.post_to_webapi()
        if 'success' in self.result and self.result['success']:
            self.main_table_id = self.result['content']['ids']['id']
            if not self.instant:
                status = 'end' if self._params_check_end else 'start'
                self.insert_pipe_detail(self.result['content']["ids"]['name'], self.main_table_id, status)
                # self.check_end(ObjectId(self.main_table_id))
            else:
                self.insert_pipe_detail(self.result['content']['ids']['name'],
                                        ObjectId(self.result['content']['ids']['id']), 'end')
        elif 'success' in self.result and not self.result['success']:
            if 'content' in self.result:
                self.insert_pipe_detail(self.result['content']['ids']['name'],
                                        ObjectId(self.result['content']['ids']['id']), 'failed')
            else:
                self.insert_pipe_detail(None, None, 'failed')
        else:
            self.bind_object.logger.error("任务接口返回值不规范: {}".format(self.result))
        self.end_fire()
        self._end_event.set()

    def post_to_webapi(self):
        """
        投递接口
        """
        # if not self.instant:
        if self._params['submit_location'] != "otu_statistic":
            self.bind_object.logger.info("submit_location: %s" % (self._params['submit_location']))
            wait_time = self.pipe_count * (self.bind_object.sub_analysis_len + 15) + random.randint(0, self.bind_object.sub_analysis_len * 5)
            self.bind_object.logger.info("等待时间%s" % (wait_time))
            gevent.sleep(wait_time)

        if self.check_params():
            self._params_check_end = True
            return self.result
        self.result = self.post()

    def insert_pipe_detail(self, table_name, table_id, status):
        """
        导入细节表的函数，先判断下是不是otu_statistic抽平分析，如果不是才会进行导入细节表
        :param table_name:表的名字
        :param table_id:主表的id
        :param status:分析的状态
        :return:
        """
        if not re.match(r'^OTUTaxonAnalysis', str(table_name)):
            if self._params['group_id'] == 'all':
                group_id = 'all'
                group_name = 'All'
            else:
                group_id = ObjectId(self._params['group_id'])
                group_name = self.db['sg_specimen_group'].find_one({'_id': ObjectId(self._params['group_id'])},
                                                                   {'group_name': 1})['group_name']
            if "level_id" in self._params:
                level_id = self._params['level_id']
                level_name = {1: 'Domain', 2: 'Kingdom', 3: 'Phylum', 4: 'Class', 5: 'Order', 6: 'Family', 7: 'Genus',
                              8: 'Species', 9: "OTU"}[level_id]
            else:
                level_id = self.min_level
                level_name = {1: 'Domain', 2: 'Kingdom', 3: 'Phylum', 4: 'Class', 5: 'Order', 6: 'Family', 7: 'Genus',
                              8: 'Species', 9: "OTU"}[self.min_level]

            if re.match(r'^16sFunctionPrediction|^LEfSe', str(table_name)):
                pipe_main_id = self.get_main_id(self.bind_object.option('pipe_id'), group_id, level_id)
            else:
                pipe_main_id = self.pipe_main_id

            insert_data = {
                "task_id": self.task_id,
                "otu_id": ObjectId(self.bind_object.otu_id),
                'group_name': group_name,
                'level_name': level_name,
                "submit_location": self._params['submit_location'],
                'params': self.json_params,
                'pipe_main_id': ObjectId(pipe_main_id),
                'pipe_batch_id': ObjectId(self.bind_object.option('pipe_id')),
                'table_id': ObjectId(table_id),
                'status': status,
                'desc': "" if status == "end" else self.result['info'],
                'level_id': str(level_id),
                "group_id": group_id,
                'type_name': self.mongo_collection,
                'table_name': table_name,
                'created_ts': datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            }
            self._pipe_detail_id = self.db['sg_pipe_detail'].insert_one(insert_data).inserted_id
        elif str(status) == "failed":
            insert_data = {
                "task_id": self.task_id,
                "otu_id": ObjectId(self.bind_object.otu_id),
                'group_name': "",
                'level_name': "",
                "submit_location": self._params['submit_location'],
                'params': self.json_params,
                'pipe_main_id': ObjectId(self.pipe_main_id),
                'pipe_batch_id': ObjectId(self.bind_object.option('pipe_id')),
                'table_id': ObjectId(table_id),
                'status': "failed",
                'desc': "因为OtuSubsample分析计算失败，后面的依赖分析都不能进行，请重新设定基本参数，再次尝试!",
                'level_id': "",
                "group_id": "",
                'type_name': self.mongo_collection,
                'table_name': table_name,
                'created_ts': datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            }
            self._pipe_detail_id = self.db['sg_pipe_detail'].insert_one(insert_data).inserted_id
        else:
            pass

    def get_main_id(self, pipe_id, group_id, level_id):
        """
        特殊用途，16s与lefse两个分析的pipe_main_id,要去查表获取
        :return:
        """
        if group_id == 'all':
            group_id = 'all'
        else:
            group_id = ObjectId(group_id)
        mongo_result = self.db['sg_pipe_main'].find_one({"pipe_batch_id": ObjectId(pipe_id),
                                                         "group_id": group_id, "level_id": int(level_id)})["_id"]
        return mongo_result

    def run_permission(self):
        if self.instant:
            for i in xrange(1000):
                if self.bind_object.count_instant_running < self.bind_object.max_instant_running:
                    break
                gevent.sleep(3)
            else:
                self.bind_object.logger.warn(
                    "任务等待时间达到上限3000s,直接运行: {}".format(self.api))
            self.bind_object.count_instant_running += 1
        else:
            for i in xrange(1000):
                if self.bind_object.count_submit_running < self.bind_object.max_submit_running:
                    break
                gevent.sleep(3)
            else:
                self.bind_object.logger.warn(
                    "任务等待时间达到上限3000s,直接运行: {}".format(self.api))
            self.bind_object.count_submit_running += 1
        self.bind_object.logger.info("任务开始投递: {}, 当前即时任务数: {}, 投递任务数: {}".format(
            self.api, self.bind_object.count_instant_running, self.bind_object.count_submit_running))

    def url_params_format(self):
        temp_params = copy.deepcopy(self._params)
        dumps_list = ["group_detail", "second_group_detail", "filter_json"]
        for i in dumps_list:
            if i in temp_params and temp_params[i]:
                temp_params[i] = json.dumps(
                    temp_params[i], sort_keys=True, separators=(',', ':'))
        # if not self.instant:
            # temp_params['pipe_id'] = str(self.pipe_main_id)
        temp_params['batch_id'] = str(self.bind_object.option('pipe_id'))
        temp_params['batch_task_id'] = str(self.bind_object.option('batch_task_id'))
        return temp_params

    def post(self):
        """
        执行post
        """
        self.url = self.bind_object.url + self.api
        httpHandler = urllib2.HTTPHandler(debuglevel=1)
        httpsHandler = urllib2.HTTPSHandler(debuglevel=1)
        opener = urllib2.build_opener(httpHandler, httpsHandler)
        urllib2.install_opener(opener)
        params = self.url_params_format()
        data = urllib.urlencode(params)
        url = "%s?%s" % (self.url, self.bind_object.signature)
        request = urllib2.Request(url, data)
        for i in range(3):
            try:
                response = gevent_url_fetch(request)[1]
                # response = urllib2.urlopen(request)
            except (urllib2.HTTPError, urllib2.URLError, httplib.HTTPException, SocketError) as e:
                self.bind_object.logger.warn(
                    "接口投递失败, 第{}次(最多三次尝试), url: {}, 错误:{}".format(i + 1, self.api, e))
                gevent.sleep(50)  # 休眠50s 目的是为了等端口重启好
            else:
                the_page = response.read()
                return json.loads(the_page)
        return {"info": "分析项: {} 投递失败！".format(self.api), "success": False}

    @property
    def end_event(self):
        return self._end_event

    def check_end(self, id, timeout=3600):
        """
        检查投递任务是否结束  目前只用于特殊情况， 即参数相同，但是不由 一键化提交
        """
        for i in xrange(1000):
            self.bind_object.logger.info("第 {} 次 查询数据库，检测任务({}, {})是否完成".format(
                i + 1, self.mongo_collection, self.main_table_id))
            mongo_result = self.db[self.mongo_collection].find_one(
                {"_id": ObjectId(self.main_table_id)}, {"status": 1, "desc": 1, "name": 1})
            self.result = {'success': True, "info": mongo_result['desc'],
                           'content': {"ids": {'id': str(mongo_result['_id']), 'name': mongo_result["name"]}}}
            if mongo_result["status"] == "end":
                break
            elif mongo_result["status"] != "start":
                self.result["success"] = False
                break
            gevent.sleep(20)
        return
        # self.get_workflow_id()
        # try:
        #     self.bind_object.logger.info("开始向WPM服务请求等待任务: {} 结束".format(self.workflow_id))
        #     get = gevent_url_fetch(self.bind_object.url + "/report/" + self.workflow_id)[1]
        #     self.bind_object.logger.info("向WPM服务请求等待结果为: {}".format(get.read()))
        #     # wait(self.workflow_id, timeout=timeout)
        # except Exception as e:
        #     self.bind_object.logger.info("ERROR:{}".format(traceback.format_exc()))
        #     self.bind_object.logger.info("任务等待超时: workflow_id: {}, url: {}".format(self.workflow_id, self.api))

    def check_params(self):
        """
        检查参数，发现end状态，直接放回计算完成，发现start状态，直接监控直到结束
        """
        self.waits_params_get()  # 依赖分析的参数获取
        self.set_params_type()  # 根据每个分析的接口中对参数的打包格式，进行对应处理
        self.json_params = self.params_pack(self._params)
        result = self.db[self.mongo_collection].find({'task_id': self.task_id, 'params': self.json_params,
                                                      'status': {'$in': ['end', 'start', "failed"]}})
        if not result.count():
            self.bind_object.logger.info("参数比对没有找到相关结果: 任务: {}, collection: {}".format(
                self.api, self.mongo_collection))
            return False
        else:
            lastone = result.sort('created_ts', pymongo.DESCENDING)[0]
            self.bind_object.logger.info("参数比对找到已经运行的结果: 任务: {}, 状态: {}, collection: {}, _id: {}".format(
                self.api, lastone['status'], self.mongo_collection, lastone["_id"]))
            self.main_table_id = lastone['_id']
            if lastone['status'] == 'end':
                self.result = {'success': True, "info": lastone['desc'],
                               'content': {"ids": {'id': str(lastone['_id']), 'name': lastone["name"]}}}
                return self.result
            elif lastone['status'] == 'start':
                self.check_end(lastone['_id'])
                self.result = {'success': True, "info": lastone['desc'],
                               'content': {"ids": {'id': str(lastone['_id']), 'name': lastone["name"]}}}
                result = self.db[self.mongo_collection].find_one(
                    {'_id': lastone['_id']}, {'status': 1, 'desc': 1, "name": 1})
                if result['status'] != 'end':  # 任务完成后检查状态
                    self.result['success'] = False
                return self.result
        return False

    def waits_params_get(self):
        """
        在子类中重写，获取等待的分析结果
        """
        pass
        return self._params

    def set_out_params(self):
        """
        在子类中重写，设置结果到out_params,方便其他分析使用
        """
        pass

    def set_params_type(self):
        """
        在子类中重写，设置参数的格式，与每个接口中类型对应
        :return:
        """
        pass
        return self._params


class BetaSampleDistanceHclusterTree(Submit):
    """
    """

    def waits_params_get(self):
        self._params['otu_id'] = self.waits[0].out_params['otu_id']
        return self._params


class SubmitOtuSubsample(Submit):
    # def set_out_params(self):
    #     self.out_params['otu_id'] = self.result["content"]['ids']['id']
    def set_out_params(self):
        for i in xrange(1000):
            mongo_result = self.db['sg_otu'].find_one({"_id": ObjectId(self.result["content"]['ids']['id'])})["status"]
            self.bind_object.logger.info("第{}次查询抽平结果表！ new_otu_id:{}, status:{}".format
                                         (i + 1, self.result["content"]['ids']['id'], mongo_result))
            if mongo_result == "end":
                self.out_params['otu_id'] = self.result["content"]['ids']['id']
                break
            elif mongo_result == "failed":
                self.success = False
                break
            gevent.sleep(15)


class AlphaDiversityIndex(Submit):

    # def set_out_params(self):new_
    #     self.out_params['alpha_diversity_id'] = self.result[
    #         "content"]['ids']['id']
    def set_out_params(self):
        for i in xrange(1000):
            mongo_result = self.db['sg_alpha_diversity'].find_one({"_id": ObjectId(self.result["content"]['ids']['id'])})["status"]
            self.bind_object.logger.info("第{}次查询多样性指数！ alpha_diversity_id:{}, status:{}".format
                                         (i + 1, self.result["content"]['ids']['id'], mongo_result))
            if mongo_result == "end":
                self.out_params['alpha_diversity_id'] = self.result["content"]['ids']['id']
                break
            elif mongo_result == "failed":
                self.success = False
                break
            gevent.sleep(15)

    def waits_params_get(self):
        self._params['otu_id'] = self.waits[0].out_params['otu_id']
        return self._params


class AlphaTtest(Submit):

    def waits_params_get(self):
        self._params['alpha_diversity_id'] = self.waits[0].out_params['alpha_diversity_id']
        self._params['otu_id'] = self.waits[1].out_params['otu_id']
        # del self._params['level_id']  # 如果要差异检验也是按照不同分类水平来的话就不要删除level_id
        return self._params


class SixteensPrediction(BetaSampleDistanceHclusterTree):
    def set_params_type(self):
        del self._params['level_id']
        return self._params


class SpeciesLefseAnalyse(BetaSampleDistanceHclusterTree):
    def set_params_type(self):
        self._params['end_level'] = int(self._params['end_level'])
        self._params['start_level'] = int(self._params['start_level'])
        self._params['strict'] = int(self._params['strict'])
        if re.search(r'\.0$', self._params['lda_filter']):
            self._params['lda_filter'] = int(float(self._params['lda_filter']))
        elif re.search(r'\..*$', self._params['lda_filter']):
            self._params['lda_filter'] = float(self._params['lda_filter'])
        else:
            self._params['lda_filter'] = int(self._params['lda_filter'])
        del self._params['level_id']
        return self._params


class OtunetworkAnalyse(BetaSampleDistanceHclusterTree):
    pass


class AlphaRarefactionCurve(BetaSampleDistanceHclusterTree):
    def set_params_type(self):
        self._params['freq'] = int(self._params['freq'])
        sort_index = self._params['index_type'].split(',')
        sort_index.sort()
        sort_index = ','.join(sort_index)
        self._params['index_type'] = sort_index
        return self._params


class RandomforestAnalyse(BetaSampleDistanceHclusterTree):
    pass


class RocAnalyse(BetaSampleDistanceHclusterTree):
    pass


class BetaMultiAnalysisPlsda(BetaSampleDistanceHclusterTree):
    pass


class BetaMultiAnalysisPearsonCorrelation(BetaSampleDistanceHclusterTree):
    pass


class BetaMultiAnalysisRdaCca(BetaSampleDistanceHclusterTree):
    pass


class HcHeatmap(BetaSampleDistanceHclusterTree):
    pass


class BetaMultiAnalysisAnosim(BetaSampleDistanceHclusterTree):
    pass


class SpeciesDifferenceMultiple(BetaSampleDistanceHclusterTree):
    def set_params_type(self):
        self._params['coverage'] = float(self._params['coverage'])
        return self._params


class BetaMultiAnalysisResults(BetaSampleDistanceHclusterTree):
    pass


class BetaMultiAnalysisPcoa(BetaSampleDistanceHclusterTree):
    pass


class BetaMultiAnalysisDbrda(BetaSampleDistanceHclusterTree):
    pass


class SpeciesDifferenceTwoGroup(BetaSampleDistanceHclusterTree):
    def set_params_type(self):
        self._params['ci'] = float(self._params['ci'])
        self._params['coverage'] = float(self._params['coverage'])
        return self._params


class BetaMultiAnalysisNmds(BetaSampleDistanceHclusterTree):
    pass


class OtuGroupAnalyse(BetaSampleDistanceHclusterTree):
    pass


class OtuVenn(BetaSampleDistanceHclusterTree):
    pass


class BetaMultiAnalysisPca(BetaSampleDistanceHclusterTree):
    pass


class CorrNetworkAnalyse(BetaSampleDistanceHclusterTree):
    def set_params_type(self):
        self._params['abundance'] = int(self._params['abundance'])
        self._params['coefficient'] = float(self._params['coefficient'])
        self._params['lable'] = float(self._params['lable'])
        return self._params


class OtuPanCore(BetaSampleDistanceHclusterTree):

    def __init__(self, *args, **kwargs):
        """
        pancore有两个主表，需要做额外的处理，例如，算作两个分析，all_counts + 1
        """
        super(BetaSampleDistanceHclusterTree, self).__init__(*args, **kwargs)
        self.bind_object.all_count += 1

    def _submit(self, waits, timeout):
        """投递任务"""
        for i in waits:
            i.end_event.get(timeout=timeout)
            if not i.success:
                self.bind_object.rely_error(self, i)
                return
        self.run_permission()
        self.post_to_webapi()
        if 'success' in self.result and self.result['success']:
            for i in self.result['content']['ids']:
                self.main_table_id = ObjectId(i['id'])
                self.insert_pipe_detail(i['name'], ObjectId(i['id']), 'end')
        elif 'success' in self.result and not self.result['success']:
            if 'content' in self.result:
                for i in self.result['content']['ids']:
                    self.main_table_id = ObjectId(i['id'])
                    self.insert_pipe_detail(i['name'], ObjectId(i['id']), 'failed')
            else:
                self.insert_pipe_detail(None, None, 'failed')
        else:
            self.bind_object.logger.error("任务接口返回值不规范: {}".format(self.result))
        self.end_fire()
        self._end_event.set()


    def check_params(self):
        """
        检查参数，发现end状态，直接放回计算完成，发现start状态，直接监控直到结束
        """
        self.waits_params_get()  # 依赖分析的参数获取
        self.set_params_type()  # 根据每个分析的接口中对参数的打包格式，进行对应处理
        self.json_params = self.params_pack(self._params)
        result = self.db[self.mongo_collection].find({'task_id': self.task_id, 'params': self.json_params,
                                                      'status': {'$in': ['end', 'start', "failed"]}})
        if not result.count():
            self.bind_object.logger.info("参数比对没有找到相关结果: 任务: {}, collection: {}".format(
                self.api, self.mongo_collection, ))
            return False
        else:
            lastone = result.sort('created_ts', pymongo.DESCENDING)[0]
            self.bind_object.logger.info("参数比对找到已经运行的结果: 任务: {}, 状态: {}, collection: {}, _id: {}".format(
                self.api, lastone['status'], self.mongo_collection, lastone["_id"]))
            self.main_table_id = lastone['_id']
            if lastone['status'] == 'end':
                ids = self.find_pan_core_ids(main_table_id=self.main_table_id)
                self.result = {'success': True, "info": lastone['desc'],
                               'content': {"ids": [{'id': str(lastone['_id']), 'name': lastone["name"]},
                                                   {'id': ids["id"], 'name': ids['name']}]}}
                return self.result
            elif lastone['status'] == 'start':
                self.check_end(lastone['_id'])
                self.result = {'success': True, "info": lastone['desc'],
                               'content': {"ids": {'id': str(lastone['_id']), 'name': lastone["name"]}}}
                result = self.db[self.mongo_collection].find_one(
                    {'_id': lastone['_id']}, {'status': 1, 'desc': 1, "name": 1})
                if result['status'] != 'end':  # 任务完成后检查状态
                    self.result['success'] = False
                return self.result
        return False

    def find_pan_core_ids(self, main_table_id):
        """
        pan_core分析有两个主表，找core的主表id
        :param main_table_id:
        :return:
        """
        result = self.db[self.mongo_collection].find_one({"_id": ObjectId(main_table_id)})
        result1 = self.db[self.mongo_collection].find({"unique_id": str(result['unique_id'])})
        for m in result1:
            if str(m["_id"]) != str(main_table_id):
                ids = {"id": str(m["_id"]), "name": m["name"]}
        return ids


class PlotTree(BetaSampleDistanceHclusterTree):
    def set_params_type(self):
        self._params['color_level_id'] = int(self._params['color_level_id'])
        return self._params


class Enterotyping(BetaSampleDistanceHclusterTree):
    pass


class gevent_HTTPConnection(httplib.HTTPConnection):

    def connect(self):
        import socket
        from gevent import socket as cosocket
        if self.timeout is socket._GLOBAL_DEFAULT_TIMEOUT:
            timeout = cosocket._GLOBAL_DEFAULT_TIMEOUT
        else:
            timeout = self.timeout
        self.sock = cosocket.create_connection((self.host, self.port), timeout)


class gevent_HTTPHandler(urllib2.HTTPHandler):

    def http_open(self, request):
        return self.do_open(gevent_HTTPConnection, request)


def gevent_url_fetch(url):
    opener = urllib2.build_opener(gevent_HTTPHandler)
    resp = opener.open(url)
    return resp.headers, resp

# if __name__ == "__main__":
#     PipeSubmitTool().run_webapitest()
