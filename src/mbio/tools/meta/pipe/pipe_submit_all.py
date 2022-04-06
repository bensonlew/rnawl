## !/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/python
# -*- coding: utf-8 -*-
# __author__ = "hongdongxuan"

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mainapp.libs.signature import CreateSignature
from biocluster.config import Config
from mainapp.libs.param_pack import *
from socket import error as SocketError
import urllib2
import urllib
import httplib
import sys
import os
import json
from bson.objectid import ObjectId
import gevent
import time
import re


class PipeSubmitAllAgent(Agent):
    """
    用于submit子分析的接口,针对的是所有的分析，后面还会有医口与农口
    version v1.0
    author: hongdongxuan
    last_modify: 2017.03.09
    """
    def __init__(self, parent):
        super(PipeSubmitAllAgent, self).__init__(parent)
        options = [
            {"name": "data", "type": "string"},
            {"name": "pipe_id", "type": "string"}
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
            raise OptionError("必须输入接口传进来的data值")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        self._memory = '50G'

    def end(self):

        super(PipeSubmitAllAgent, self).end()


class PipeSubmitAllTool(Tool):

    def __init__(self, config):
        super(PipeSubmitAllTool, self).__init__(config)
        self._version = '1.0.1'
        self._config = Config()
        self._client = self._config.mongo_client
        self.db = self._client['tsanger']


    def run_webapitest(self):
        """
        进行参数判断后投递所有的接口，如果后面要添加新的分析，主要要添加以下几个地方的内容：
        1）analysis_table：analysis_table中保存的是每个分析的submit_location与每个分析对应的主表名字（mongo或者controller中可以找到）
        2）子分析分投递型的与即时型的，投递的分析要在前面重写投递的参数，即时的只要在sub_analysis这个字典中添加就ok（注意了要有每个分析的api）
        3）参数判断的那个函数中，要根据每个子分析的controller中的参数重组格式一模一样，只有这样才能匹配上对应的params，同时要在这个函数中对应的地方写上submit_location
        4）一键化的所有的分析，都是以submit_location来区分，所以在对接的时候要着重考虑与前端的submit_location是不是一致的。
        :return:
        """
        analysis_table = {"randomforest_analyse": "sg_randomforest", "sixteens_prediction": "sg_16s",
                          "species_lefse_analyse": "sg_species_difference_lefse", "alpha_rarefaction_curve":"sg_alpha_rarefaction_curve",
                          "otunetwork_analyse": "sg_network", "roc_analyse": "sg_roc",
                          "alpha_diversity_index":"sg_alpha_diversity","alpha_ttest":"sg_alpha_ttest",
                          "beta_multi_analysis_plsda": "sg_beta_multi_analysis","beta_sample_distance_hcluster_tree": "sg_newick_tree",
                          "beta_multi_analysis_pearson_correlation":"sg_species_env_correlation","beta_multi_analysis_rda_cca": "sg_beta_multi_analysis","hc_heatmap": "sg_hc_heatmap",
                          "beta_multi_analysis_anosim": "sg_beta_multi_anosim", "species_difference_multiple": "sg_species_difference_check",
                          "beta_multi_analysis_results": "sg_species_mantel_check", "beta_multi_analysis_pcoa": "sg_beta_multi_analysis",
                          "beta_multi_analysis_dbrda": "sg_beta_multi_analysis", "species_difference_two_group": "sg_species_difference_check",
                          "beta_multi_analysis_nmds":"sg_beta_multi_analysis","otu_group_analyse":"sg_otu",
                          "otu_venn":"sg_otu_venn","beta_multi_analysis_pca": "sg_beta_multi_analysis",
                          "corr_network_analyse": "sg_corr_network", "otu_pan_core": "sg_otu_pan_core",
                          "plot_tree": "sg_phylo_tree", "enterotyping": "sg_enterotyping"}
        data = self.option('data')
        data = json.loads(data)
        group_infos = data['group_info']
        group_infos = json.loads(group_infos)  #这里group_infos才是dict
        level_ids = []  #line 109 ~112将level_id字符串进行分割存入列表中，用于后面组合
        level = str(data['level_id']).strip().split(",")
        for m in level:
            level_ids.append(m)
        sub_analysis = data['sub_analysis']
        sub_analysis = json.loads(sub_analysis)
        key_list = [] #line 116 ~ 122全局计算有多少个分析，这里有个特例，pancore，一个分析会产生两个主表
        for key in sub_analysis:
            key_list.append(str(key))
        if "otu_pan_core" in key_list:
            analysis_num = len(sub_analysis) + 1
        else:
            analysis_num = len(sub_analysis)
        client = data['client']
        method = "post"
        if client == "client01":
            base_url = 'http://bcl.i-sanger.com'
        elif client == "client03":
            # base_url = "http://10.101.203.193:9100"
            base_url = "http://192.168.12.102:9090"
        else:
            raise Exception("***client必须是client01或者client03***")
        otu_id = self.otu_subsample(data=data, client=client, base_url=base_url, method=method) #计算抽平，后面分析基于抽平后的id进行计算
        list2 = [] #用于存储分类水平与分组方案的所有的组合
        for level in level_ids:
            for group in group_infos:
                if 'env_id' in data.keys() and 'env_labs' in data.keys():
                    m = {"otu_id": str(otu_id), "level_id": str(level), "group_id": group['group_id'],
                         "group_detail": json.dumps(group['group_detail']), "env_id": str(data['env_id']),
                         "env_labs": str(data['env_labs'])}
                else:
                    m = {"otu_id": str(otu_id), "level_id": str(level), "group_id": group['group_id'],
                         "group_detail": json.dumps(group['group_detail'])}
                list2.append(m)
        all_analysis_num = analysis_num * len(list2) #计算总共有多少个分析
        all_results = []  # 存储所有子分析的ids
        ready_analysis_num = 0 #定义已经完成的分析的个数
        ready_analysis_list = [] #定时刷新一次mongo表，将每次刷新投递运行成功的结果保存到这个列表中，然后计算出最大值，用于后面即时的分析叠加
        #下面定义多样性指数与指数间差异分析的接口投递，所有的投递的接口都要在这里自行定义
        anaylsis_names = []
        alpha_diversity_index_data = {}
        alpha_ttest_data = {}
        species_lefse_analyse_data = {}
        sixteens_prediction_data = {}
        alpha_rarefaction_curve_data = {}
        otunetwork_analyse_data = {}
        roc_analyse_data = {}
        randomforest_analyse_data = {}
        for key in sub_analysis:
            anaylsis_names.append(str(key))
            if key == "alpha_diversity_index":
                alpha_diversity_index_data = sub_analysis[key]
            elif key == "alpha_ttest":
                alpha_ttest_data = sub_analysis[key]
            elif key == "sixteens_prediction":
                sixteens_prediction_data = sub_analysis[key]
            elif key == "species_lefse_analyse":
                species_lefse_analyse_data = sub_analysis[key]
            elif key == "alpha_rarefaction_curve":
                alpha_rarefaction_curve_data = sub_analysis[key]
            elif key == "otunetwork_analyse":
                otunetwork_analyse_data = sub_analysis[key]
            elif key == "roc_analyse":
                roc_analyse_data = sub_analysis[key]
            elif key == "randomforest_analyse":
                randomforest_analyse_data = sub_analysis[key]
            else:
                pass
        if "randomforest_analyse" in anaylsis_names:
            ready_analysis_num = self.randomforest_analyse(list2, randomforest_analyse_data, analysis_table,
                                                           ready_analysis_num, all_analysis_num, client, base_url,
                                                           method, all_results)
            ready_analysis_list.append(ready_analysis_num)
        else:
            pass
        if "sixteens_prediction" in anaylsis_names:
            ready_analysis_num = self.sixteens_prediction(list2, sixteens_prediction_data,
                                                          analysis_table, ready_analysis_num, all_analysis_num, client,
                                                          base_url, method, all_results)
            ready_analysis_list.append(ready_analysis_num)
        else:
            pass
        if "species_lefse_analyse" in anaylsis_names:
            ready_analysis_num = self.species_lefse_analyse(list2, species_lefse_analyse_data,
                                                            analysis_table, ready_analysis_num, all_analysis_num,
                                                            client, base_url, method, all_results)
            ready_analysis_list.append(ready_analysis_num)
        else:
            pass
        if "alpha_rarefaction_curve" in anaylsis_names:
            ready_analysis_num = self.alpha_rarefaction_curve(list2, alpha_rarefaction_curve_data,
                                                              analysis_table, ready_analysis_num, all_analysis_num,
                                                              client, base_url, method, all_results)
            ready_analysis_list.append(ready_analysis_num)
        else:
            pass

        if "otunetwork_analyse" in anaylsis_names:
            ready_analysis_num = self.otunetwork_analyse(list2, otunetwork_analyse_data, analysis_table,
                                                         ready_analysis_num, all_analysis_num, client, base_url,
                                                         method, all_results)
            ready_analysis_list.append(ready_analysis_num)
        else:
            pass
        if "roc_analyse" in anaylsis_names:
            ready_analysis_num = self.roc_analyse(list2, roc_analyse_data, analysis_table,
                                                  ready_analysis_num, all_analysis_num, client, base_url,
                                                  method, all_results)
            ready_analysis_list.append(ready_analysis_num)
        else:
            pass
        #用于判断是否有投递的分析，有投递的分析，则求最后一次刷新mongo表中完成的最大个数，否则就初始化为0
        if len(ready_analysis_list) == 0:
            pass
        else:
            ready_analysis_num = max(ready_analysis_list)

        if "alpha_ttest" in anaylsis_names and "alpha_diversity_index" in anaylsis_names:
            for n in list2:
                alpha_diversity_index_data['otu_id'] = n['otu_id']
                alpha_diversity_index_data['level_id'] = n['level_id']
                alpha_diversity_index_data['group_id'] = n['group_id']
                alpha_diversity_index_data['group_detail'] = n['group_detail']
                if len(self.params_check(alpha_diversity_index_data,
                                         analysis_table[alpha_diversity_index_data['submit_location']],
                                         alpha_diversity_index_data['submit_location'])) == 0:
                    results_alpha_diversity_index = self.run_controllers(api="/".join(str(alpha_diversity_index_data['api']).strip().split("|")),
                                                                         client=client, base_url=base_url,
                                                                         params=alpha_diversity_index_data, method=method)
                    results_alpha_diversity_index = json.loads(results_alpha_diversity_index)
                    all_results.append(results_alpha_diversity_index)
                    ready_analysis_num += 1
                    self.update_status(all_analysis_num, ready_analysis_num)
                    if 'sub_anaylsis_id' in results_alpha_diversity_index.keys():
                        alpha_diversity_id = results_alpha_diversity_index['sub_anaylsis_id']['id']
                    else:
                        alpha_diversity_id = ""
                else:
                    all_results.append(self.params_check(alpha_diversity_index_data,
                                                         analysis_table[alpha_diversity_index_data['submit_location']],
                                                         alpha_diversity_index_data['submit_location']))
                    ready_analysis_num += 1
                    self.update_status(all_analysis_num, ready_analysis_num)
                    results_alpha_diversity_index = self.params_check(alpha_diversity_index_data,
                                                                      analysis_table[alpha_diversity_index_data['submit_location']],
                                                                      alpha_diversity_index_data['submit_location'])
                    alpha_diversity_id = results_alpha_diversity_index['sub_anaylsis_id']['id']
                    print "打印出alpha_diversity_id"
                    print alpha_diversity_id
                if alpha_diversity_id:
                    alpha_ttest_data['otu_id'] = n['otu_id']
                    alpha_ttest_data['level_id'] = n['level_id']
                    alpha_ttest_data['group_id'] = n['group_id']
                    alpha_ttest_data['alpha_diversity_id'] = str(alpha_diversity_id)
                    alpha_ttest_data['group_detail'] = n['group_detail']
                    if len(self.params_check(alpha_ttest_data, analysis_table[alpha_ttest_data['submit_location']],
                                             alpha_ttest_data['submit_location'])) == 0:
                        results_alpha_ttest = self.run_controllers(api="/".join(str(alpha_ttest_data['api']).strip().split("|")),
                                                                   client=client, base_url=base_url,
                                                                   params=alpha_ttest_data, method=method)
                        results_alpha_ttest = json.loads(results_alpha_ttest)
                        all_results.append(results_alpha_ttest)
                        ready_analysis_num += 1
                        self.update_status(all_analysis_num, ready_analysis_num)
                    else:
                        all_results.append(self.params_check(alpha_ttest_data,
                                                             analysis_table[alpha_ttest_data['submit_location']],
                                                             alpha_ttest_data['submit_location']))
                        ready_analysis_num += 1
                        self.update_status(all_analysis_num, ready_analysis_num)
                else:
                    alpha_diversity_results = {"level_id": n['level_id'], "group_id": n['group_id'],
                                               "info": "多样性指数分析出错，无法进行该分析！", "success": False,
                                               "submit_location": "alpha_ttest"}
                    all_results.append(alpha_diversity_results)
                    ready_analysis_num += 1
                    self.update_status(all_analysis_num, ready_analysis_num)
        elif "alpha_diversity_index" in anaylsis_names and "alpha_ttest" not in anaylsis_names:
            for n in list2:
                alpha_diversity_index_data['otu_id'] = n['otu_id']
                alpha_diversity_index_data['level_id'] = n['level_id']
                alpha_diversity_index_data['group_id'] = n['group_id']
                alpha_diversity_index_data['group_detail'] = n['group_detail']
                if len(self.params_check(alpha_diversity_index_data,
                                         analysis_table[alpha_diversity_index_data['submit_location']],
                                         alpha_diversity_index_data['submit_location'])) == 0:
                    results_alpha_diversity_index = self.run_controllers(
                        api="/".join(str(alpha_diversity_index_data['api']).strip().split("|")),
                        client=client, base_url=base_url, params=alpha_diversity_index_data, method=method)
                    results_alpha_diversity_index = json.loads(results_alpha_diversity_index)
                    all_results.append(results_alpha_diversity_index)
                    ready_analysis_num += 1
                    self.update_status(all_analysis_num, ready_analysis_num)
                else:
                    all_results.append(self.params_check(alpha_diversity_index_data,
                                                         analysis_table[alpha_diversity_index_data['submit_location']],
                                                         alpha_diversity_index_data['submit_location']))
                    ready_analysis_num += 1
                    self.update_status(all_analysis_num, ready_analysis_num)
        else:
            pass
        #删除子分析字典中alpha_diversity_index，alpha_ttest两个分析,重构sub_analysis数组
        if "alpha_diversity_index" in anaylsis_names:
            del sub_analysis['alpha_diversity_index']
        else:
            pass
        if "alpha_ttest" in anaylsis_names:
            del sub_analysis['alpha_ttest']
        else:
            pass
        if "sixteens_prediction" in anaylsis_names:
            del sub_analysis['sixteens_prediction']
        else:
            pass
        if "species_lefse_analyse" in anaylsis_names:
            del sub_analysis['species_lefse_analyse']
        else:
            pass
        if "alpha_rarefaction_curve" in anaylsis_names:
            del sub_analysis['alpha_rarefaction_curve']
        else:
            pass
        if "otunetwork_analyse" in anaylsis_names:
            del sub_analysis['otunetwork_analyse']
        else:
            pass
        if "roc_analyse" in anaylsis_names:
            del sub_analysis['roc_analyse']
        else:
            pass
        if "randomforest_analyse" in anaylsis_names:
            del sub_analysis['randomforest_analyse']
        else:
            pass
        #定义没有依赖的分析的通用接口投递
        for info in list2:
            for anaylsis in sub_analysis:
                sub_analysis[anaylsis]['otu_id'] = info['otu_id']
                sub_analysis[anaylsis]['level_id'] = info['level_id']
                sub_analysis[anaylsis]['group_id'] = info['group_id']
                sub_analysis[anaylsis]['group_detail'] = info['group_detail']
                if anaylsis in ['beta_multi_analysis_pca', 'beta_multi_analysis_rda_cca', 'beta_multi_analysis_dbrda',
                                'beta_multi_analysis_pearson_correlation', 'beta_multi_analysis_results'] \
                        and 'env_id' in info.keys() and 'env_labs' in info.keys():
                    sub_analysis[anaylsis]['env_id'] = info['env_id']
                    sub_analysis[anaylsis]['env_labs'] = info['env_labs']
                else:
                    pass
            gevent_data = []
            submit_location_dict = {}
            for key in sub_analysis:
                name = []
                data = []  # 用于存储子分析参数的具体数值
                submit_info = {}
                if len(self.params_check(sub_analysis[key], analysis_table[sub_analysis[key]['submit_location']],
                                         sub_analysis[key]['submit_location'])) == 0:
                    for key1 in sub_analysis[key]:
                        name.append(key1)
                        data.append(sub_analysis[key][key1])
                        submit_info = {"api": "/".join(str(sub_analysis[key]['api']).strip().split('|')),
                                       "submit_location": str(sub_analysis[key]['submit_location'])}
                    name = ";".join(name)
                    data = ";".join(data)
                    api = submit_info['api']
                    method = "post"
                    # time.sleep(60)
                    submit_location_dict[str(gevent.spawn(self.webapitest, method=method, api=api, name=name, data=data, client=client, base_url=base_url))] = submit_info['submit_location']
                    gevent_data.append(gevent.spawn(self.webapitest, method=method, api=api, name=name, data=data, client=client, base_url=base_url))
                else:
                    all_results.append(
                        self.params_check(sub_analysis[key], analysis_table[sub_analysis[key]['submit_location']],
                                          sub_analysis[key]['submit_location']))
                    ready_analysis_num += 1
                    self.update_status(all_analysis_num, ready_analysis_num)
            print "gevent_data", gevent_data
            print "submit_location_dict", submit_location_dict
            gevent.joinall(gevent_data)
            for result_page in gevent_data:
                a = str(result_page)
                print "test", submit_location_dict[a]
                result = json.loads(json.loads(result_page.value))
                print result['success']
                ready_analysis_num += 1
                self.update_status(all_analysis_num, ready_analysis_num)
                if result['success'] == True:
                    results = {"level_id": info['level_id'], "group_id": info['group_id'],
                               "sub_anaylsis_id": result['content']['ids'], "success": True,
                               "submit_location": submit_location_dict[str(result_page)]}
                    all_results.append(results)
                else:
                    results = {"level_id": info['level_id'], "group_id": info['group_id'],
                               "info": result['info'], "success": False,
                               "submit_location": submit_location_dict[str(result_page)]}
                    all_results.append(results)
                    # return_page = self.webapitest(method, api, name, data, client, base_url)
                    # result = json.loads(return_page)
                    # result = json.loads(result)
                    # ready_analysis_num += 1
                    # self.update_status(all_analysis_num, ready_analysis_num)
                    # if result['success'] == True:
                    #     results = {"level_id": info['level_id'], "group_id": info['group_id'],
                    #                "sub_anaylsis_id": result['content']['ids'], "success": True,
                    #                "submit_location": submit_info['submit_location']}
                    #     all_results.append(results)
                    # else:
                    #     results = {"level_id": info['level_id'], "group_id": info['group_id'],
                    #                "info": result['info'], "success": False,
                    #                "submit_location": submit_info['submit_location']}
                    #     all_results.append(results)
                # else:
                #     all_results.append(self.params_check(sub_analysis[key], analysis_table[sub_analysis[key]['submit_location']], sub_analysis[key]['submit_location']))
                #     ready_analysis_num += 1
                #     self.update_status(all_analysis_num, ready_analysis_num)
        output_table = os.path.join(self.work_dir, "ids.txt")
        with open(output_table, "w") as w:
            w.write(str(all_results))


    def params_check(self, ever_analysis_params, analysis_table_name, submit_location):
        """
        这里是进行每个分析的参数判断，如果参数一样就不进行投递。
        :param ever_analysis_params: 每个分析的所有的参数
        :param analysis_table_name: 这个分析的主表
        :param submit_location: 该分析的submit_location，注意这个字段一定要与前端对应，所有的分析都是以此字段进行区别。
        :return:
        """
        if not isinstance(ever_analysis_params, dict):
            raise Exception("传入的params不是一个字典！")
        if submit_location in ["beta_multi_analysis_rda_cca", "beta_multi_analysis_pca", "beta_multi_analysis_pcoa",
                               "beta_multi_analysis_nmds", "beta_multi_analysis_plsda", "beta_multi_analysis_dbrda", "otunetwork_analyse",
                               "randomforest_analyse", "otu_pan_core", "roc_analyse", "otu_group_analyse", "hc_heatmap",
                               "beta_multi_analysis_pearson_correlation", "otu_venn",
                               "beta_multi_analysis_anosim", "beta_sample_distance_hcluster_tree", "enterotyping",
                               "beta_multi_analysis_results"]:
            my_param = dict()
            for key in ever_analysis_params:
                if key == "level_id":
                    my_param[key] = int(ever_analysis_params[key])
                elif key == "group_detail":
                    my_param[key] = group_detail_sort(ever_analysis_params[key])
                elif key == "api":
                    pass
                else:
                    my_param[key] = ever_analysis_params[key]
            new_params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))
        elif submit_location in ["species_lefse_analyse"]:
            my_param = dict()
            for key in ever_analysis_params:
                if key == "group_detail":
                    my_param[key] = group_detail_sort(ever_analysis_params[key])
                elif key == "second_group_detail" and ever_analysis_params['second_group_detail']:
                    my_param[key] = group_detail_sort(json.loads(ever_analysis_params[key]))
                elif key == "second_group_detail" and not ever_analysis_params['second_group_detail']:
                    my_param[key] = ever_analysis_params[key]
                elif key == "strict" or key == "start_level" or key == "end_level":
                    my_param[key] = int(ever_analysis_params[key])
                elif key == "api":
                    pass
                elif key == "level_id":
                    pass
                elif key == "lda_filter" and re.search(r'\.0$', ever_analysis_params['lda_filter']):
                    my_param[key] = int(float(ever_analysis_params[key]))
                elif key == "lda_filter" and re.search(r'\..*$', ever_analysis_params['lda_filter']):
                    my_param[key] = float(ever_analysis_params[key])
                elif key == "lda_filter" and not re.search(r'\.0$', ever_analysis_params['lda_filter']) and not \
                        re.search(r'\..*$', ever_analysis_params['lda_filter']):
                    my_param[key] = int(ever_analysis_params[key])
                else:
                    my_param[key] = ever_analysis_params[key]
            new_params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))
        elif submit_location in ["alpha_rarefaction_curve"]:
            my_param = dict()
            for key in ever_analysis_params:
                if key == "level_id" or key == "freq":
                    my_param[key] = int(ever_analysis_params[key])
                elif key == "group_detail":
                    my_param[key] = group_detail_sort(ever_analysis_params[key])
                elif key == "index_type":
                    sort_index = ever_analysis_params[key].split(',')
                    sort_index.sort()
                    sort_index = ','.join(sort_index)
                    my_param[key] = sort_index
                elif key == "api":
                    pass
                else:
                    my_param[key] = ever_analysis_params[key]
            new_params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))
        elif submit_location in ["sixteens_prediction", "alpha_ttest"]:
            my_param = dict()
            for key in ever_analysis_params:
                if key == "group_detail":
                    my_param[key] = group_detail_sort(ever_analysis_params[key])
                elif key == "api":
                    pass
                elif key == "level_id":
                    pass
                else:
                    my_param[key] = ever_analysis_params[key]
            new_params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))
        elif submit_location in ["corr_network_analyse", "plot_tree", "species_difference_two_group", "species_difference_multiple"]:
            my_param = dict()
            for key in ever_analysis_params:
                if key == "level_id" or key == "abundance" or key == "color_level_id":
                    my_param[key] = int(ever_analysis_params[key])
                elif key == "group_detail":
                    my_param[key] = group_detail_sort(ever_analysis_params[key])
                elif key == "api":
                    pass
                elif key == "lable" or key == "coefficient" or key == "ci" or key == "coverage":
                    my_param[key] = float(ever_analysis_params[key])
                else:
                    my_param[key] = ever_analysis_params[key]
            new_params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))
        elif submit_location in ["alpha_diversity_index"]:
            my_param = dict()
            for key in ever_analysis_params:
                if key == "group_detail":
                    my_param[key] = group_detail_sort(ever_analysis_params[key])
                elif key == "api":
                    pass
                elif key == "index_type":
                    sort_index = ever_analysis_params[key].split(',')
                    sort_index = ','.join(sort_index)
                    my_param[key] = sort_index
                else:
                    my_param[key] = ever_analysis_params[key]
            new_params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))
        else:
            my_param = dict()
            for key in ever_analysis_params:
                my_param[key] = ever_analysis_params[key]
            new_params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))
        collection_sg_otu = self.db["sg_otu"]
        result = collection_sg_otu.find_one({"_id": ObjectId(str(ever_analysis_params['otu_id']))})
        task_id = result['task_id']
        collections = self.db[analysis_table_name]
        # print json.dumps(new_params)
        results = {}
        if submit_location == "otu_pan_core":
            result1 = collections.find_one({"task_id": task_id, "params": new_params, "status": "end"})
            analysis_id = []
            if result1:
                ids = {"id": str(result1['_id']), "name": result1['name']}
                analysis_id.append(ids)
                another_id = collections.find({"task_id": task_id, "created_ts": result1['created_ts']})
                for m in another_id:
                    if str(m['_id']) != str(result1['_id']) and json.dumps(m['params']) == json.dumps(new_params):
                        ids1 = {"id": str(m['_id']), "name": m['name']}
                        analysis_id.append(ids1)
                        results = {"level_id": ever_analysis_params['level_id'],
                                   "group_id": ever_analysis_params['group_id'],
                                   "sub_anaylsis_id": analysis_id, "success": "unknown",
                                   "submit_location": submit_location}
                        break
        else:
            if collections.find_one({"task_id": task_id, "params": new_params, "status": "start"}):
                m = collections.find_one({"task_id": task_id, "params": new_params, "status": "start"})
                results = {"level_id": ever_analysis_params['level_id'],
                           "group_id": ever_analysis_params['group_id'],
                           "sub_anaylsis_id": {"id": str(m['_id']), "name": m['name']}, "success": "unknown",
                           "submit_location": submit_location}
            elif collections.find_one({"task_id": task_id, "params": new_params, "status": "end"}):
                m = collections.find_one({"task_id": task_id, "params": new_params, "status": "end"})
                results = {"level_id": ever_analysis_params['level_id'],
                           "group_id": ever_analysis_params['group_id'],
                           "sub_anaylsis_id": {"id": str(m['_id']), "name": m['name']}, "success": True,
                           "submit_location": submit_location}
        # print "打印出不需要计算的分析"
        # print results
        return results

    def run_controllers(self, api, client, base_url, params, method, header=None):
        """
        用于重构每个子分析的输入参数，（参考sub_anaylsis中的参数）并进行投递计算,注输入的参数必须是字典格式
        params样例：{\"submit_location\": \"corrnetwork_analyse\", \"api\": \"meta/corr_network\",
        \"task_type\": \"reportTask\",\"lable\": \"0.03\", \"ratio_method\": \"pearson\", \"coefficient\":
        \"0.08\", \"abundance\": \"150\"}  注这里面要添加level_id与group_id
        :return:
        """
        if not isinstance(params, dict):
            raise Exception("传入的params不是一个字典！")
        name = []
        data = []
        results = {}
        for key in params:
            name.append(key)
            data.append(params[key])
        name = ";".join(name)
        data = ";".join(data)
        time.sleep(0.5)
        return_page = self.webapitest(method, api, name, data, client, base_url)
        result = json.loads(return_page)
        result = json.loads(result)
        if result['success'] == True:
            results = {"level_id": params['level_id'], "group_id": params['group_id'],
                       "sub_anaylsis_id": result['content']['ids'], "success": True,
                       "submit_location": params['submit_location']}
        else:
            results = {"level_id": params['level_id'], "group_id": params['group_id'],
                       "info": result['info'], "success": False,
                       "submit_location": params['submit_location']}

        return json.dumps(results)

    def webapitest(self, method, api, name, data, client, base_url, header=None):
        """
        :param method: [get,post]
        :param api: the api address and url
        :param name: names for data, split by \";\"
        :param data: data or files for names, split by \";\"
        :param client: client name
        :param base_url: the base url of api, http://192.168.12.102:9010
        :param header: use header to submit signature info
        :return:
        """
        if data and not name:
            print("Error:must give the name option when the data is given!")
            sys.exit(1)
        if name and not data:
            print("Error:must give the data option when the name is given!")
            sys.exit(1)
        datas = {}
        if name:
            names_list = name.split(";")
            data_list = data.split(";")
            for index in range(len(names_list)):
                if index < len(data_list):
                    if os.path.isfile(data_list[index]):
                        with open(data_list[index], "r") as f:
                            content = f.readlines()
                            content = "".join(content)
                    else:
                        content = data_list[index]
                    datas[names_list[index]] = content
                else:
                    datas[names_list[index]] = ""
        httpHandler = urllib2.HTTPHandler(debuglevel=1)
        httpsHandler = urllib2.HTTPSHandler(debuglevel=1)
        opener = urllib2.build_opener(httpHandler, httpsHandler)

        urllib2.install_opener(opener)
        data = urllib.urlencode(datas)

        signature_obj = CreateSignature(client)
        signature = {
            "client": signature_obj.client,
            "nonce": signature_obj.nonce,
            "timestamp": signature_obj.timestamp,
            "signature": signature_obj.signature
        }

        signature = urllib.urlencode(signature)
        url = "%s/%s" % (base_url, api)
        if not header:
            if "?" in url:
                url = "%s&%s" % (url, signature)
            else:
                url = "%s?%s" % (url, signature)

        if method == "post":
            # print("post data to url %s ...\n\n" % url)
            # print("post data:\n%s\n" % data)
            request = urllib2.Request(url, data)
            # print "request"
            # print request
        else:
            if data:
                if "?" in url:
                    url = "%s&%s" % (url, data)
                else:
                    url = "%s?%s" % (url, data)
            else:
                url = "%s%s" % (base_url, api)
            # print("get url %s ..." % url)
            request = urllib2.Request(url)

        if header:
            request.add_header('client', signature_obj.client)
            request.add_header('nonce', signature_obj.nonce)
            request.add_header('timestamp', signature_obj.timestamp)
            request.add_header('signature', signature_obj.signature)

        try:
            response = urllib2.urlopen(request)
        except (urllib2.HTTPError, urllib2.URLError, httplib.HTTPException, SocketError), e:
            print("%s-- \n" % e)
            time.sleep(50)    #这一步一定要休眠50s 目的是为了等端口重启好，再测试过程中出现接口连续重启3次的情况
            #这一步骤是用于处理接口投递时出现异常导致接口投递失败，当失败的时候休眠5s然后重新投递2次，badstatsline属于httplib.HTTPException异常
            for n in range(2):
                time.sleep(10)
                try:
                    response = urllib2.urlopen(request)
                except (urllib2.HTTPError, urllib2.URLError, httplib.HTTPException, SocketError), e:
                    print("%s \n" % e)
                    the_page = {"info": "submit failed!", "success": False}
                    the_page = json.dumps(the_page)
                else:
                    the_page = response.read()
                    print("Return page:\n%s" % the_page)
                    break
        else:
            the_page = response.read()
            print("Return page:\n%s" % the_page)
            print '\n'
            # print type(the_page)
        return json.dumps(the_page)


    def update_status(self, all_analysis_num, ready_analysis_num, analysis_type=None, all_results=None,
                      analysis_table=None):
        """
        用于进度条显示
        all_analysis_num:参数是指所有的分析的个数（分类水平*分组方案*分析个数）
        ready_analysis_num:参数是指已经完成的分析的个数
        all_results: 针对的是投递的分析，需要根据all_results中的id去查找状态
        analysis_table：每个分析submit_location对应的主表
        analysis_type:submit or instant
        pipe_id:批次表的主表id
        :return:
        """
        main_table_id = self.option("pipe_id")
        collection_pipe = self.db["sg_pipe_batch"]
        if analysis_type and analysis_type == "submit":
            for m in all_results:
                collections_sub_analysis = self.db[analysis_table[m["submit_location"]]]
                if 'sub_anaylsis_id' in m.keys():
                    if "id" in m["sub_anaylsis_id"].keys():
                        sub_anaylsis_main_id = str(m['sub_anaylsis_id']['id'])
                        result = collections_sub_analysis.find_one({"_id": ObjectId(sub_anaylsis_main_id)})
                        if result and result['status'] != 'start':
                            ready_analysis_num += 1
                        else:
                            pass
                    elif isinstance(m['sub_anaylsis_id'], list):
                        for n in m['sub_anaylsis_id']:
                            sub_anaylsis_main_id = str(n['id'])
                            result = collections_sub_analysis.find_one({"_id": ObjectId(sub_anaylsis_main_id)})
                            if result and result['status'] != 'start':
                                ready_analysis_num += 1
                            else:
                                pass
                else:
                    ready_analysis_num += 1

        percent = str(ready_analysis_num) + "/" + str(all_analysis_num)
        data = {
            "percent": percent
        }
        try:
            collection_pipe.update({"_id": ObjectId(main_table_id)}, {'$set': data}, upsert=False)
        except:
            raise Exception("进度条更新失败了！")
        return ready_analysis_num

    def run(self):
        super(PipeSubmitAllTool, self).run()
        self.run_webapitest()
        self.end()


    def otu_subsample(self, data, client, base_url, method):
        """
        进行抽平计算，返回的值是抽平后的otu_id
        :param data:
        :param client:
        :param base_url:
        :param method:
        :return:
        """
        if not isinstance(data, dict):
            raise Exception("传入的params不是一个字典！")
        collection_sg_otu = self.db["sg_otu"]
        result = collection_sg_otu.find_one({"_id": ObjectId(str(data['otu_id']))})
        task_id = result['task_id']
        from_id = result['from_id']
        my_param = dict()
        my_param["group_id"] = data['group_id']
        my_param['otu_id'] = from_id
        my_param["submit_location"] = "otu_statistic"
        my_param["size"] = data['size']
        my_param["filter_json"] = filter_json_sort(data['filter_json'])
        my_param["group_detail"] = group_detail_sort(data['group_detail'])
        my_param["task_type"] = "reportTask"
        params_otu = param_pack(my_param)
        if collection_sg_otu.find_one({"task_id": task_id, "params": params_otu, "status": "end"}):
            otu_id = str(collection_sg_otu.find_one({"task_id": task_id, "params": params_otu, "status": "end"})["_id"])
        else:
            otu_id = "none"
        if otu_id == "none":
            params = {"otu_id": str(data['otu_id']), "group_id": data['group_id'], "group_detail": data['group_detail'],
                      "submit_location": "otu_statistic",
                      "filter_json": data['filter_json'], "task_type": "reportTask", "size": data['size'],
                      "level_id": "9"}
            api_statistic = "meta/otu_subsample"
            results_statistic = self.run_controllers(api=api_statistic, client=client, base_url=base_url, params=params,
                                                     method=method)
            results_statistic = json.loads(results_statistic)
            if 'sub_anaylsis_id' in results_statistic.keys() and results_statistic['success'] == True:
                otu_id = results_statistic['sub_anaylsis_id']['id']
            else:
                raise Exception("进行抽平或者样本筛选的过程失败，程序被终止！")
        return otu_id

    def otunetwork_analyse(self, list2, otunetwork_analyse_data, analysis_table, ready_analysis_num, all_analysis_num,
                           client, base_url, method, all_results):
        """
        :param list2: 这是所有的分组方案与分类水平进行随机组合的列表
        :param otunetwork_analyse_data: 每个子分析需要除了otu_id，group_id，group_detail参数之外的其它高级参数
        :param analysis_table: 这个在run_webapitest开始地方定义，将所有分析的submit_location与main_table_name对应的字典
        :param client:
        :param base_url:
        :param method:
        :param all_results:用于存储所有分析的返回的id以及表的名字，或者返回的错误信息。该列表是个全局列表
        :return: 返回的值是用于进度条计算的，返回的是投递的分析已经完成的个数
        """
        for n in list2:
            ready_analysis_num = 0
            otunetwork_analyse_data['otu_id'] = n['otu_id']
            otunetwork_analyse_data['level_id'] = n['level_id']
            otunetwork_analyse_data['group_id'] = n['group_id']
            otunetwork_analyse_data['group_detail'] = n['group_detail']
            if len(self.params_check(otunetwork_analyse_data,
                                     analysis_table[otunetwork_analyse_data['submit_location']],
                                     otunetwork_analyse_data['submit_location'])) == 0:
                result_otunetwork_analyse = self.run_controllers(
                    api="/".join(str(otunetwork_analyse_data['api']).strip().split("|")),
                    client=client, base_url=base_url, params=otunetwork_analyse_data, method=method)
                result_otunetwork_analyse = json.loads(result_otunetwork_analyse)
                all_results.append(result_otunetwork_analyse)
                ready_analysis_num = self.update_status(all_analysis_num=all_analysis_num, analysis_type="submit",
                                                        ready_analysis_num=ready_analysis_num, all_results=all_results,
                                                        analysis_table=analysis_table)
            else:
                all_results.append(self.params_check(otunetwork_analyse_data,
                                                     analysis_table[otunetwork_analyse_data['submit_location']],
                                                     otunetwork_analyse_data['submit_location']))
                ready_analysis_num = self.update_status(all_analysis_num=all_analysis_num, analysis_type="submit",
                                                        ready_analysis_num=ready_analysis_num, all_results=all_results,
                                                        analysis_table=analysis_table)
        return ready_analysis_num


    def roc_analyse(self, list2, roc_analyse_data, analysis_table, ready_analysis_num, all_analysis_num, client,
                    base_url, method, all_results):
        """
        计算ROC的程序
        """
        for n in list2:
            ready_analysis_num = 0
            roc_analyse_data['otu_id'] = n['otu_id']
            roc_analyse_data['level_id'] = n['level_id']
            roc_analyse_data['group_id'] = n['group_id']
            roc_analyse_data['group_detail'] = n['group_detail']
            if len(self.params_check(roc_analyse_data,
                                     analysis_table[roc_analyse_data['submit_location']],
                                     roc_analyse_data['submit_location'])) == 0:
                result_roc_analyse = self.run_controllers(api="/".join(str(roc_analyse_data['api']).strip().split("|")),
                                                          client=client, base_url=base_url,
                                                          params=roc_analyse_data, method=method)
                result_roc_analyse = json.loads(result_roc_analyse)
                all_results.append(result_roc_analyse)
                ready_analysis_num = self.update_status(all_analysis_num=all_analysis_num,
                                                        ready_analysis_num=ready_analysis_num, analysis_type="submit",
                                                        all_results=all_results, analysis_table=analysis_table)
            else:
                all_results.append(self.params_check(roc_analyse_data,
                                                     analysis_table[roc_analyse_data['submit_location']],
                                                     roc_analyse_data['submit_location']))
                ready_analysis_num = self.update_status(all_analysis_num=all_analysis_num,
                                                        ready_analysis_num=ready_analysis_num, analysis_type="submit",
                                                        all_results=all_results, analysis_table=analysis_table)
        return ready_analysis_num

    def randomforest_analyse(self, list2, randomforest_analyse_data, analysis_table, ready_analysis_num,
                             all_analysis_num, client, base_url, method, all_results):
        """
        计算随机森林的程序
        """
        for n in list2:
            ready_analysis_num = 0
            randomforest_analyse_data['otu_id'] = n['otu_id']
            randomforest_analyse_data['level_id'] = n['level_id']
            randomforest_analyse_data['group_id'] = n['group_id']
            randomforest_analyse_data['group_detail'] = n['group_detail']
            if len(self.params_check(randomforest_analyse_data,
                                     analysis_table[randomforest_analyse_data['submit_location']],
                                     randomforest_analyse_data['submit_location'])) == 0:
                result_randomforest_analyse = self.run_controllers(
                    api="/".join(str(randomforest_analyse_data['api']).strip().split("|")),
                    client=client, base_url=base_url,
                    params=randomforest_analyse_data, method=method)
                result_randomforest_analyse = json.loads(result_randomforest_analyse)
                all_results.append(result_randomforest_analyse)
                ready_analysis_num = self.update_status(all_analysis_num=all_analysis_num,
                                                        ready_analysis_num=ready_analysis_num, analysis_type="submit",
                                                        all_results=all_results, analysis_table=analysis_table)
            else:
                all_results.append(self.params_check(randomforest_analyse_data,
                                                     analysis_table[randomforest_analyse_data['submit_location']],
                                                     randomforest_analyse_data['submit_location']))
                ready_analysis_num = self.update_status(all_analysis_num=all_analysis_num,
                                                        ready_analysis_num=ready_analysis_num, analysis_type="submit",
                                                        all_results=all_results, analysis_table=analysis_table)
        return ready_analysis_num


    def sixteens_prediction(self, list2, sixteens_prediction_data, analysis_table, ready_analysis_num,
                            all_analysis_num, client, base_url, method, all_results):
        """
        计算16s功能预测的程序
        """
        for n in list2:
            ready_analysis_num = 0
            sixteens_prediction_data['otu_id'] = n['otu_id']
            sixteens_prediction_data['level_id'] = n['level_id']
            sixteens_prediction_data['group_id'] = n['group_id']
            sixteens_prediction_data['group_detail'] = n['group_detail']
            if len(self.params_check(sixteens_prediction_data,
                                     analysis_table[sixteens_prediction_data['submit_location']],
                                     sixteens_prediction_data['submit_location'])) == 0:
                result_sixteens_prediction = self.run_controllers(
                    api="/".join(str(sixteens_prediction_data['api']).strip().split("|")),
                    client=client, base_url=base_url,
                    params=sixteens_prediction_data, method=method)
                result_sixteens_prediction = json.loads(result_sixteens_prediction)
                all_results.append(result_sixteens_prediction)
                ready_analysis_num = self.update_status(all_analysis_num=all_analysis_num,
                                                        ready_analysis_num=ready_analysis_num, analysis_type="submit",
                                                        all_results=all_results, analysis_table=analysis_table)
            else:
                all_results.append(self.params_check(sixteens_prediction_data,
                                                     analysis_table[sixteens_prediction_data['submit_location']],
                                                     sixteens_prediction_data['submit_location']))
                ready_analysis_num = self.update_status(all_analysis_num=all_analysis_num,
                                                        ready_analysis_num=ready_analysis_num, analysis_type="submit",
                                                        all_results=all_results, analysis_table=analysis_table)
        return ready_analysis_num

    def species_lefse_analyse(self, list2, species_lefse_analyse_data, analysis_table, ready_analysis_num,
                              all_analysis_num, client, base_url, method, all_results):
        """
        计算lefse的程序
        """
        for n in list2:
            ready_analysis_num = 0
            species_lefse_analyse_data['otu_id'] = n['otu_id']
            species_lefse_analyse_data['level_id'] = n['level_id']
            species_lefse_analyse_data['group_id'] = n['group_id']
            species_lefse_analyse_data['group_detail'] = n['group_detail']
            if len(self.params_check(species_lefse_analyse_data,
                                     analysis_table[species_lefse_analyse_data['submit_location']],
                                     species_lefse_analyse_data['submit_location'])) == 0:
                result_species_lefse_analyse = self.run_controllers(
                    api="/".join(str(species_lefse_analyse_data['api']).strip().split("|")),
                    client=client, base_url=base_url,
                    params=species_lefse_analyse_data, method=method)
                result_species_lefse_analyse = json.loads(result_species_lefse_analyse)
                all_results.append(result_species_lefse_analyse)
                ready_analysis_num = self.update_status(all_analysis_num=all_analysis_num,
                                                        ready_analysis_num=ready_analysis_num, analysis_type="submit",
                                                        all_results=all_results, analysis_table=analysis_table)
            else:
                all_results.append(self.params_check(species_lefse_analyse_data,
                                                     analysis_table[species_lefse_analyse_data['submit_location']],
                                                     species_lefse_analyse_data['submit_location']))
                ready_analysis_num = self.update_status(all_analysis_num=all_analysis_num,
                                                        ready_analysis_num=ready_analysis_num, analysis_type="submit",
                                                        all_results=all_results, analysis_table=analysis_table)
        return ready_analysis_num

    def alpha_rarefaction_curve(self, list2, alpha_rarefaction_curve_data, analysis_table, ready_analysis_num,
                                all_analysis_num, client, base_url, method, all_results):
        """
        计算多样性指数的程序
        """
        for n in list2:
            ready_analysis_num = 0  #这个要放在每个循环里面，每循环一次都要初始化
            alpha_rarefaction_curve_data['otu_id'] = n['otu_id']
            alpha_rarefaction_curve_data['level_id'] = n['level_id']
            alpha_rarefaction_curve_data['group_id'] = n['group_id']
            alpha_rarefaction_curve_data['group_detail'] = n['group_detail']
            if len(self.params_check(alpha_rarefaction_curve_data,
                                     analysis_table[alpha_rarefaction_curve_data['submit_location']],
                                     alpha_rarefaction_curve_data['submit_location'])) == 0:
                result_alpha_rarefaction_curve = self.run_controllers(
                    api="/".join(str(alpha_rarefaction_curve_data['api']).strip().split("|")),
                    client=client, base_url=base_url,
                    params=alpha_rarefaction_curve_data, method=method)
                result_alpha_rarefaction_curve = json.loads(result_alpha_rarefaction_curve)
                all_results.append(result_alpha_rarefaction_curve)
                ready_analysis_num = self.update_status(all_analysis_num=all_analysis_num,
                                                        ready_analysis_num=ready_analysis_num, analysis_type="submit",
                                                        all_results=all_results, analysis_table=analysis_table)
            else:
                all_results.append(self.params_check(alpha_rarefaction_curve_data,
                                                     analysis_table[alpha_rarefaction_curve_data['submit_location']],
                                                     alpha_rarefaction_curve_data['submit_location']))
                ready_analysis_num = self.update_status(all_analysis_num=all_analysis_num,
                                                        ready_analysis_num=ready_analysis_num, analysis_type="submit",
                                                        all_results=all_results, analysis_table=analysis_table)
        return ready_analysis_num
