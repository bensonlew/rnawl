# !/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python
# -*- coding: utf-8 -*-
# __auther__ = 'yuguo'
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
from collections import defaultdict
import datetime
import gevent
import json
import codecs
import re


class MedManageModule(Module):
    """
    该module用于管理医学样本报告分析投递顺序。
    lasted modified by hongdong@20180820
    """
    def __init__(self, work_id):
        super(MedManageModule, self).__init__(work_id)
        options = [
            {"name": "sample_tab", "type": "infile", "format": "medical.tab"},  # 拆分表
            {"name": "sample_dir", "type": "outfile", "format": "medical.sample_dir"},  # 拆分后的fastq文件路径
            {"name": "batch_id", "type": "string"},  # 数据拆分主表id
            {"name": "board_batch", "type": "string"},
            {"name": "member_id", "type": "string"},  # 用户ID
            {"name": "project_types", "type": "string", "default": "all"}  # 项目类型，用于区分拆分完之后激发的工作流
        ]
        self.add_option(options)
        self.Sanmple_dispath = None
        self.task_logs = defaultdict(list)  # {samples:[work_key, dispath_time, samples, info]}
        self.task1 = []  # 第一批投递任务
        self.task2 = []  # 第二批投递任务，必须等第一批运行完成后触发提交
        self.task3 = []  # 第三批投递任务，如果存在不能组建家系的S
        self.tid = 0
        self.nipt_n = 0
        self.pt_n = 0
        self.split_type = None

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("sample_tab"):
            raise OptionError("缺少必要参数split_tab：样本拆分表")
        return True

    def creat_tasks(self, project_types):
        self.Sanmple_dispath = SanmpleDispath(self, self.option("sample_tab").prop['path'])
        self.Sanmple_dispath.get_works(project_types)
        self.logger.info("works:{}".format(self.Sanmple_dispath.works))
        self.logger.info("task1:{}".format(self.Sanmple_dispath.task1))
        self.logger.info("task2:{}".format(self.Sanmple_dispath.task2))
        self.logger.info("task3:{}".format(self.Sanmple_dispath.task3))
        self.task1 = self.add_task_tools(self.Sanmple_dispath.task1)
        self.task2 = self.add_task_tools(self.Sanmple_dispath.task2)
        self.task3 = self.add_task_tools(self.Sanmple_dispath.task3)
        # self.logger.info(self.task1)
        # self.logger.info(self.task2)

    def add_task_tools(self, task):
        num = len(task)
        if num > 0:
            tools = [self.add_tool("medical.datasplit_v3.submit_work") for i in range(num)]
            return tools
        else:                       # modified by HD at 20180106 解决当没有亲子样本的时候，报错
            return []
        # self.logger.info(self.task1)       

    def run_pt_snp(self):
        """
        提交task1：PT所有的F、M跑call snp。
        """
        self.task1[0].set_options({
            "api_name": "s/med/pt_v3/pt_analysis",
            "params": json.dumps({
                'batch_id': self.option('batch_id'),
                'board_batch': self.option('board_batch'),
                'member_id': self.option('member_id'),
                'samples': json.dumps(self.Sanmple_dispath.works['pt-fm-snp']),
                'fastq_path': self.option('sample_dir').prop['path'] + "PT/",
                'ana_type': "F-M"  # 用于区分是FM call snp 还是S call snp
            }),
            "failed_to_stop": "true"
        })
        self.task1[0].run()
        dispath_time = datetime.datetime.now().strftime("%Y%m%d:%H%M%S")
        self.Sanmple_dispath.add_task_time('task1', 0, dispath_time)
        self.logger.info('开始运行task1：run_pt_snp\n[batch_id:{},\nsamples:{},\nall_snp_num:{},\nall_family_num:{}]'
                         .format(self.option('batch_id'), json.dumps(self.Sanmple_dispath.works['pt-fm-snp']),
                                 self.Sanmple_dispath.pt_counts[0], self.Sanmple_dispath.pt_counts[1]))

    def run_pt_family(self):
        """
        提交task2：按照家系单位投递PT的家系分析，如果有son，则流程里先做son的snp
        """
        for tid in sorted(self.Sanmple_dispath.task2.keys()):
            if re.search(r'^pt', self.Sanmple_dispath.task2[tid][2]):
                family = self.Sanmple_dispath.task2[tid][3]
                f = re.sub(r'\*', '', family)
                son, mom, dad = f.split(',')
                self.task2[tid].set_options({
                    "api_name": "s/med/pt_v3/pt_analysis",
                    "params": json.dumps({
                        'batch_id': self.option('batch_id'),
                        'board_batch': self.option('board_batch'),
                        'member_id': self.option('member_id'),
                        'mom_id': mom,
                        'dad_id': dad,
                        'son_id': son,
                        'fastq_path': self.option('sample_dir').prop['path'] + "PT/",
                        'err_min_num': '11'
                    })
                })
                gevent.sleep(1)
                self.task2[tid].run()
                dispath_time = datetime.datetime.now().strftime("%Y%m%d:%H%M%S")
                self.Sanmple_dispath.add_task_time('task2', tid, dispath_time)
                self.logger.info('开始运行task2：run_pt_family\n[batch_id:{},\nfamily:{},\nall_snp_num:{},\nal'
                                 'l_family_num:{}]'.format(self.option('batch_id'), family,
                                                           self.Sanmple_dispath.pt_counts[0],
                                                           self.Sanmple_dispath.pt_counts[1]))

    def run_nipt(self):
        """
        提交task2：按紧急批次投递nipt workflow
        """
        for tid in sorted(self.Sanmple_dispath.task2.keys()):
            if re.search('nipt', self.Sanmple_dispath.task2[tid][2]):
                self.task2[tid].set_options({    # modified by hd 20180106 self.nipt_n改为tid
                    "api_name": "s/med/nipt",    # modified by hd 20180103 140-148行
                    "params": json.dumps({
                        'batch_id': self.option('batch_id'),
                        'member_id': self.option('member_id'),
                        'sample_list': self.Sanmple_dispath.task2[tid][3],
                        'fastq_path': self.option('sample_dir').prop['path'] + "NIPT/",
                        'urgent': self.Sanmple_dispath.task2[tid][2],
                        'single': self.split_type
                        # 'all_nipt_num': self.Sanmple_dispath.nipt_counts   # modified by HD 20180204
                    })
                })
                self.task2[tid].run()
                dispath_time = datetime.datetime.now().strftime("%Y%m%d:%H%M%S")
                self.Sanmple_dispath.add_task_time('task2', 0, dispath_time)
                self.logger.info('开始运行task2：run_{}\nnipt_counts:{},'
                                 .format(self.Sanmple_dispath.task2[tid][2], self.Sanmple_dispath.nipt_counts))

    def run_pt_s_snp(self):
        """
        提交task3：将组建不了家系胎儿样本S进行第三批投递，运行callsnp，以备在以后家系完整的板里能调出snp结果
        """
        self.task3[0].set_options({
            "api_name": "s/med/pt_v3/pt_analysis",
            "params": json.dumps({
                'batch_id': self.option('batch_id'),
                'board_batch': self.option('board_batch'),
                'member_id': self.option('member_id'),
                'samples': json.dumps(self.Sanmple_dispath.works['pt-s-snp']),
                'fastq_path': self.option('sample_dir').prop['path'] + "PT/",
                'all_snp_num': self.Sanmple_dispath.pt_counts[0],
                'all_family_num': self.Sanmple_dispath.pt_counts[1],
                'ana_type': "S"  # 用于区分是FM call snp 还是S call snp
            }),
            "failed_to_stop": "true"
        })
        self.task3[0].run()
        dispath_time = datetime.datetime.now().strftime("%Y%m%d:%H%M%S")
        self.Sanmple_dispath.add_task_time('task3', 0, dispath_time)
        self.logger.info('开始运行task3：run_s_snp\n[batch_id:{},\nsamples:{},\nall_snp_num:{},\nall_family_num:{}]'
                         .format(self.option('batch_id'), json.dumps(self.Sanmple_dispath.works['pt-s-snp']),
                                 self.Sanmple_dispath.pt_counts[0], self.Sanmple_dispath.pt_counts[1]))

    def write_work_log(self, logfile, tasks):
        """
        记录任务列表
        self.work_dir + "/workflow_stats.log"
        :return:
        """
        self.logger.info("开始写worklog：{}".format(logfile))
        with open(logfile, "w") as f:
            for t in tasks:
                if len(t) > 0:
                    for key in sorted(t.keys()):
                            f.write('\t'.join(t[key]) + "\n")

    def submit_task1(self):
        if self.option("project_types") not in ["nipt"]:
            if self.Sanmple_dispath.works['pt-fm-snp']:  # modified by hd 20180104 snp样本可能为空
                self.run_pt_snp()

    def submit_task2(self):
        if self.option("project_types") == 'pt':
            self.run_pt_family()
        elif self.option("project_types") == 'nipt':
            # pass
            self.run_nipt()
        else:
            self.run_pt_family()
            self.run_nipt()

    def submit_task3(self):
        if self.option("project_types") not in ['nipt']:
            self.run_pt_s_snp()

    def run(self):
        """
        修改运行逻辑modified by hongdong 20171212
        :return:
        """
        super(MedManageModule, self).run()
        self.split_type = self.api.api('medical.paternity_test_v3.paternity_test_v3').find_split_type(
            self.option("batch_id"))
        self.creat_tasks(self.option("project_types"))
        self.write_work_log(self.work_dir + "/tasks_plans.log", self.Sanmple_dispath.tasks)
        if len(self.task1) > 0:
            if len(self.task2) > 0:
                self.task1[0].on('end', self.submit_task2)
                if len(self.task3) > 0:
                    self.on_rely(self.task2, self.submit_task3)
                    self.on_rely(self.task3, self.end)
                else:
                    self.on_rely(self.task2, self.end)
                self.submit_task1()
            else:
                if len(self.task3) > 0:
                    self.task1[0].on('end', self.submit_task3)
                    self.on_rely(self.task3, self.end)
                else:
                    self.task1[0].on('end', self.end)
                self.submit_task1()
        elif len(self.task2) > 0:
            if len(self.task3) > 0:
                self.on_rely(self.task2, self.submit_task3)
                self.on_rely(self.task3, self.end)
            else:
                self.on_rely(self.task2, self.end)
            self.submit_task2()
        elif len(self.task3) > 0:    # add by hd 20180110 如果这块板子全部是亲子胎儿
            self.on_rely(self.task3, self.end)
            self.submit_task3()
        else:
            self.end()

    def end(self):
        self.write_work_log(self.work_dir + "/tasks_stats.log", self.Sanmple_dispath.tasks)
        super(MedManageModule, self).end()


class SanmpleDispath(object):
    """
    作者：kefei.huang
    时间：20171115
    last modified by hongdong at 20180820
    使用datasplit的结果，根据特定的规则，对样品的投递顺序进行排序。
    目前主要有两个workflow，亲子和产筛。
    目前来看，将这两个workflow设定为并行状态，不排除以后根据特别的情况，做成workflow串行的状态。
    所以目前，根据投递时间还有投递的资源，粗略的控制投递顺序。
    emergency={
        "加急1天"：4,
        "加急2天"：3,
        "加急3天"：2,
        "不加急"：1
    }
    """

    def __init__(self, bind_object, sample_tab):  # sample_tab path
        self.bind_object = bind_object
        self.sample_tab = sample_tab
        self.works = {
            'pt-fm-snp': [],  # F和M callsnp的样本列表
            'pt-emergency': [],  # 2  s_m_f的family列表
            'pt-hurry': [],  # 1
            "pt-hurryhurry": [],
            'pt-normal': [],
            'pt-s-snp': [],  # 未组家系的S单独call snp
            'pt-outboard': [],  # 按case组成家系但F\M\S都不在该板的家系
            'nipt-emergency': [],
            'nipt-hurry': [],
            'nipt-normal': [],
            'other': []
        }
        self.emergency = {
            "特加急": 2,
            "加急": 1
        }
        self.pt_emergency = {
            "加急1天": 4,
            "加急2天": 3,
            "加急3天": 2,
            "不加急": 1
        }
        self.emergency_encode()
        self.pt_caseDic = defaultdict(list)  # pt样本家系case  {case_name: [f1,m,s,f2]}
        self.pt_statuDic = {}  # pt样本紧急程度 {sample_id: pt-特加急}
        self.pt_familys = defaultdict(list)
        self.task1 = defaultdict(list)
        self.task2 = defaultdict(list)
        self.task3 = defaultdict(list)
        self.notask = defaultdict(list)
        self.t2 = 0
        self.nt = 0
        self.pt_counts = [0, 0]  # call snp num, family num
        self.nipt_counts = 0
        self.tasks = []

    def emergency_encode(self):
        for k, v in self.emergency.items():
            new = k.encode("utf-8").decode("utf-8")
            self.emergency.pop(k)
            self.emergency[new] = v
        for k, v in self.pt_emergency.items():
            new = k.encode("utf-8").decode("utf-8")
            self.pt_emergency.pop(k)
            self.pt_emergency[new] = v

    def add_task_time(self, task, tid, dispath_time):
        task = getattr(self, task)
        task[tid][1] = dispath_time

    def get_works(self, project_types):
        with codecs.open(self.sample_tab, "r", "utf-8") as SPLIT:
            titles = SPLIT.readline().split(',')
            sample_idx, analysis_idx, emergency_idx = [titles.index(i) for i in ['sample_name', 'analysis_type',
                                                                                 'emergency']]
            for line in SPLIT:
                linetemp = line.split(",")
                pattern = linetemp[emergency_idx]
                # self.bind_object.logger.info(linetemp)
                if linetemp[analysis_idx] in ('pt', 'dcpt', "wqcf"):
                    self.pt_statuDic[linetemp[sample_idx]] = pattern
                    case = linetemp[sample_idx].split("-")[0]
                    self.pt_caseDic[case].append(linetemp[sample_idx])
                    result = self.bind_object.api.api('medical.paternity_test_v3.paternity_test_v3').get_case_sampls(
                        case)
                    if result:
                        self.pt_caseDic[case].extend([temp['sample_id'] for temp in result])
                elif linetemp[analysis_idx] == 'nipt':
                    if pattern in self.emergency.keys():
                        if self.emergency[pattern] == 2:
                            self.works['nipt-emergency'].append(linetemp[sample_idx])
                        elif self.emergency[pattern] == 1:
                            self.works['nipt-hurry'].append(linetemp[sample_idx])
                    else:
                        self.works['nipt-normal'].append(linetemp[sample_idx])
                else:
                    self.works['other'].append(linetemp[sample_idx])
                    self.notask[self.nt] = ['notask-' + str(self.nt), '-', 'other', linetemp[sample_idx]]
                    self.nt += 1
        self.get_pt_works(project_types)
        # 暂时注释掉产筛
        for emerge in ['nipt-emergency', 'nipt-hurry', 'nipt-normal']:
            if self.works[emerge]:   # 过滤掉任务为空的情况，不再添加到tool中 add by hongdong 20180103
                if project_types in ['all', 'nipt']:
                    self.task2[self.t2] = ['task2-' + str(self.t2), '-', emerge, ', '.join(self.works[emerge])]
                    self.t2 += 1
                    self.nipt_counts += len(self.works[emerge])
                else:
                    self.notask[self.nt] = ['notask-' + str(self.nt), '-', emerge, ', '.join(self.works[emerge])]
                    self.nt += 1
        self.tasks = [self.task1, self.task2, self.task3, self.notask]

    def get_pt_works(self, project_types):
        self.bind_object.logger.info(self.pt_caseDic)
        for k in self.pt_caseDic.keys():
            fathers = []
            mothers = []
            sons = []
            for each in set(self.pt_caseDic[k]):
                if re.search("-S", each) and re.match(r'^TEST.*', each):   # 添加胎儿的测试样本
                    self.works['pt-s-snp'].append(each)
                if re.search("-S", each):
                    sons.append(each)
                elif re.search("-M", each):
                    mothers.append(each)
                    if each in self.pt_statuDic.keys():
                        self.works['pt-fm-snp'].append(each)
                elif re.search("-F", each):
                    fathers.append(each)
                    if each in self.pt_statuDic.keys():
                        self.works['pt-fm-snp'].append(each)
                else:
                    self.notask[self.nt] = ['notask' + str(self.nt), '-', 'pt-unknow', each]
                    self.nt += 1
                    self.bind_object.logger.info(each + " is not one of father,mother and son, check files")
            if len(fathers) > 0 and len(mothers) > 0 and len(sons) > 0:
                for son_each in sons:
                    for mother_each in mothers:
                        for father_each in fathers:
                            family = "_".join([son_each, mother_each, father_each])
                            self.check_pt_family(family)
            else:  # 无法组建家系时, S需要单独call snp
                if len(sons) > 0:
                    for s in sons:
                        if s in self.pt_statuDic.keys():
                            if re.match(r'^TEST.*', s):   # add by hd 20180205 过滤掉测试样本
                                continue
                            self.works['pt-s-snp'].append(s)
        if project_types in ['all', 'pt']:
            if len(self.works['pt-fm-snp']) > 0:
                self.task1[0] = ['task1-0', '-', 'pt-fm-snp', ','.join(self.works['pt-fm-snp'])]
            if len(self.works['pt-s-snp']) > 0:
                self.task3[0] = ['task3-0', '-', 'pt-s-snp', ','.join(self.works['pt-s-snp'])]
        else:
            self.notask[self.nt] = ['notask-' + str(self.nt), '-', 'pt-fm-snp', ','.join(self.works['pt-fm-snp'])]
            self.nt += 1
            self.notask[self.nt] = ['notask-' + str(self.nt), '-', 'pt-s-snp', ','.join(self.works['pt-s-snp'])]
            self.nt += 1
        for emerge in ['pt-emergency', 'pt-hurryhurry', 'pt-hurry', 'pt-normal']:
            if len(self.works[emerge]) > 0:   # modified by hd 20180104 解决亲子家系列表为空的时候运行出错
                for family in self.works[emerge]:
                    if project_types in ['all', 'pt']:
                        self.task2[self.t2] = ['task2-' + str(self.t2), '-', emerge, family]
                        self.t2 += 1
                        self.pt_counts[1] += 1
                    else:
                        self.notask[self.nt] = ['notask-' + str(self.nt), '-', emerge, family]
                        self.nt += 1
        self.pt_counts[0] = len(self.pt_statuDic.keys())
        self.bind_object.logger.info("pt_statuDic:{}".format(self.pt_statuDic.keys()))
        self.bind_object.logger.info("pt_counts:{}".format(self.pt_counts))
        # return tid

    def check_pt_family(self, family):
        """
        检查pt一个家系组合的紧急类型
        emergency={
        "加急1天"：4,
        "加急2天"：3,
        "加急3天"：2,
        "不加急"：1
        }
        """
        family_tag = []
        report_type = []  # 0：normal，1：hurry，2：emergercy，10：样本不在该板
        for temp in family.split('_'):
            if temp in self.pt_statuDic.keys():
                family_tag.append(temp)
                if self.pt_statuDic[temp] in self.pt_emergency.keys():
                    report_type.append(self.pt_emergency[self.pt_statuDic[temp]])
                else:
                    report_type.append(0)  # 非加急：普通
            else:
                report_type.append(10)  # ##不在这份列表里面
                family_tag.append(temp + "*")
        family = ','.join(family_tag)
        if 4 in report_type:
            self.works['pt-emergency'].append(family)   # 加急1天
        elif 3 in report_type:
            self.works['pt-hurryhurry'].append(family)   # 加急2天
        elif 2 in report_type:
            self.works['pt-hurry'].append(family)   # 加急3天
        elif sum(report_type) == 30:  # F\M\S都不在该板
            self.works['pt-outboard'].append(family)
        else:
            self.works['pt-normal'].append(family)
