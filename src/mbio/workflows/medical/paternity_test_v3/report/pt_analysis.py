# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from biocluster.api.file.lib.transfer import MultiFileTransfer
from bson import ObjectId
import gevent
import json
import os
import re


class PtAnalysisWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        用途及运行逻辑：
        1)当option中有samples这个参数的时候，该流程就是进行call snp，call完就结束
        2)当option中没有samples时，该流程用于家系的分析出报告，如果家系中的3个样本都有tab文件，就直接进行家系分析与查重模块，
        如果有样本的tab文件不存在就进行call SNP然后在进行家系分析与查重
        auther: hongdong
        last modified by hongdong@20180821
        :param wsheet_object:
        """
        self._sheet = wsheet_object
        super(PtAnalysisWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fastq_path", "type": "infile", "format": "sequence.fastq_dir"},  # fastq所在路径(文件夹
            {"name": "err_min_num", "type": "int", "default": 9},  # 该选项用于要循环多少个错配，默认是2-8
            {"name": "dad_id", "type": "string"},
            {"name": "mom_id", "type": "string"},
            {"name": "old_mom_id", "type": "string"},
            {"name": "son_id", "type": "string"},
            {"name": "batch_id", "type": "string"},  # 每拆分一次就有一个batch_id，仅用于区分批次
            {"name": "update_info", "type": "string"},  # 用于后面更新主表
            {"name": "member_id", "type": "string"},  # 会员ID
            {"name": "family_id", "type": "string"},
            {"name": "samples", "type": "string"}  # samples是用于样本的call snp，如果有该字段就只进行call snp，
            # 否则就是进行家系分析
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.fastq_call_snp_modules = []
        self.modules_analysis = []
        self.modules_analysis_homohybrid = []
        self.tools_result = []
        self.tools_dedup = []
        self.similarity_tools = []
        self.other_family_workflow = []
        self.rdata = []
        self.path = None  # 用于存储报告中的图片
        self.father_err_id, self.father_err_id_homohybrid = None, None
        self.dad_list = []
        self.board_batch = None
        self.api_pt = self.api.api('medical.paternity_test_v3.paternity_test_v3')
        self.father_id = ''
        # self.father_id = json.loads(self.option('update_info')).keys()[0]
        self.ref_data = self.config.SOFTWARE_DIR + "/database/human/pt_ref/tab_v3_fm/"
        self.targets_bedfile = self.config.SOFTWARE_DIR + "/database/human/pt_ref/dcpt.test.1010.bed"  # 3000个SNP位点
        self.ref_fasta = self.config.SOFTWARE_DIR + "/database/human/hg38.chromosomal_assembly/ref.fa"  # 参考基因组
        self.ref_point = self.config.SOFTWARE_DIR + "/database/human/pt_ref/targets.bed.rda"  # 每个位点的基因型概率
        self.is_update = "true"  # 在进度条更新的过程中用于区分是胎儿样本是家系分析中进行call snp
        # 还是单独流程中进行call snp的
        self.err_list = [5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40, 45, 50]
        self.err_list_homohybrid = [5, 6, 7, 8, 9, 10, 15, 20]  # 胎儿纯合与杂合的运行错配个数
        # self.err_list = [5]
        # self.err_list_homohybrid = [5]
        self.batch_id = ''
        self.ms_match = None
        self.dp_correction = None
        self.msismatch = True
        self.ms_son_id = ''   # S降噪后的样本id
        self.son_is_dcpt = True  # 胎儿是否是多重还是杂补的
        self.ms_is_match = True  # 家系分析结束后判断M与S的位点是否匹配
        self.split_type = None  # 拆分类型，PE或者SE
        self.end_by_workflow = True   # 判断是正常结束还是因为缺少样本手动终止
        self.run_first = True  # 用于判断是不是又重新组建了家系，在该工作流中进行投递
        self.bucket = self.config.get_project_region_bucket(project_type="pt_v3").rstrip("/")

    def check_options(self):
        """检查参数设置"""
        if not self.option('fastq_path').is_set:
            raise OptionError('必须提供fastq文件所在的路径')
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def fastq_call_snp_run(self, samples):
        """
        用于激发样本的call snp module
        :param samples:
        :return:
        """
        n = 0
        for sample_name in samples:
            fastq_call_snp = self.add_module("medical.paternity_test_v3.fastq_call_snp")
            analysis_type = self.api_pt.get_sample_analysis_type(sample_name)
            self.step.add_steps('fastq2tab{}'.format(n))
            fastq_call_snp.set_options({
                "sample_id": sample_name,
                "fastq_path": self.option("fastq_path"),
                "ref_fasta": self.ref_fasta,
                "targets_bedfile": self.targets_bedfile,
                "batch_id": self.batch_id,
                "analysis_type": analysis_type,
                "board_batch": self.board_batch,
                "is_update": self.is_update,  # "true" or "False"
                "split_type": self.split_type
            })
            step = getattr(self.step, 'fastq2tab{}'.format(n))
            step.start()
            fastq_call_snp.on('end', self.finish_update, 'fastq2tab{}'.format(n))
            self.fastq_call_snp_modules.append(fastq_call_snp)
            n += 1
        for j in range(len(self.fastq_call_snp_modules)):
            self.fastq_call_snp_modules[j].on('end', self.set_output, 'fastq2tab')

    def family_analysis_run(self):
        """
        家系分析生成父权制与调试表,如果只进行一个错配的计算的话，要将self.result_info_run启动 放到link到output之后。
        """
        task = []
        if self.check_tab():   # 用于检测胎儿call snp完成后，父本母本是否存在
            self.get_dad_no_analysis()   # 用于检测后面是否需要查重
            if not self.son_is_dcpt:
                preg_tab = self.output_dir + '/' + self.ms_son_id + '.tab'
            else:
                preg_tab = self.dp_correction.output_dir + "/" + self.ms_son_id + '.tab'
            # 胎儿杂合的流程
            for m in self.err_list:
                family_analysis = self.add_module("medical.paternity_test_v3.family_analysis")
                family_analysis.set_options({
                    "dad_tab": self.output_dir + '/' + self.option("dad_id") + '.tab',
                    "mom_tab": self.output_dir + '/' + self.option("mom_id") + '.tab',
                    "preg_tab": preg_tab,
                    "ref_point": self.ref_point,
                    "err_min": m
                })
                family_analysis.on('end', self.set_output, 'pt_analysis')
                task.append(family_analysis)
            # 增加胎儿纯合与杂合的流程
            for n in self.err_list_homohybrid:
                family_analysis_homohybrid = self.add_module("medical.paternity_test_v3.family_analysis_homohybrid")
                family_analysis_homohybrid.set_options({
                    "dad_tab": self.output_dir + '/' + self.option("dad_id") + '.tab',
                    "mom_tab": self.output_dir + '/' + self.option("mom_id") + '.tab',
                    "preg_tab": preg_tab,
                    "ref_point": self.ref_point,
                    "err_min": n
                })
                family_analysis_homohybrid.on('end', self.set_output, 'pt_analysis_homohybrid')
                task.append(family_analysis_homohybrid)
            self.on_rely(task, self.check_ms_is_match)
            for module in task:
                gevent.sleep(1)
                module.run()
        else:
            self.end_by_workflow = False
            self.end()

    def check_ms_is_match(self):
        """
        在家系分析运行完成后检查下ms是否匹配  --WQ181495-M-1_WQ181495-S-1_info_show.txt
        :return:
        """
        if self.api_pt.check_is_match(self.output_dir + '/pt_result_5/{}'.format(
                "_".join([self.option("mom_id"), self.ms_son_id, "info_show.txt"]))):
            self.logger.info("ms匹配，将进行后面的查重或者个体识别分析！")
            if len(self.dad_list) == 0:
                self.sample_similarity_run()
            else:
                self.dedup_run()
        else:
            self.msismatch = False
            self.ms_match_run()

    def dp_correction_run(self):
        """
        降噪tool运行
        :return:
        """
        self.logger.info("开始MS降噪")
        if self.run_first:
            self.export_tab(self.option("mom_id"))
            mom_tab = self.output_dir + '/' + self.option("mom_id") + '.tab'
            preg_tab = self.output_dir + '/' + self.option("son_id") + '.tab'
            self.export_tab(self.option("son_id"))
        else:
            self.export_tab(self.option("old_mom_id"))
            mom_tab = self.output_dir + '/' + self.option("old_mom_id") + '.tab'
            preg_tab = self.output_dir + '/' + self.option("son_id").split("_")[0] + '.tab'
            self.export_tab(self.option("son_id").split("_")[0])
        self.dp_correction = self.add_tool("medical.paternity_test_v3.dp_correction")
        self.dp_correction.set_options({
            "mom_tab": mom_tab,
            "preg_tab": preg_tab,
            "new_id": self.ms_son_id
        })
        self.dp_correction.on("end", self.set_output, "dp_correction")
        self.dp_correction.on("end", self.family_analysis_run)
        self.dp_correction.run()

    def find_sample_id(self, sample_id):
        """
        根据正则表达式去匹配查找对应的样本id，这里要注意的是'WQ181495-2-M-1' 与'WQ181495-M-1' 其实是同一个家系
        re.match("(WQ.*)-(\d)-(M|F|S)(.*)", a)
        :return:
        """
        m = re.match("(WQ.*)-(\d)-(M|F|S)(.*)", sample_id)
        if m:
            new_id = m.group(3) + m.group(4)
        else:
            new_id = "-".join(sample_id.split('-')[1:])
        return new_id

    def ms_match_run(self):
        """
        当该家系母本与胎儿不匹配的时候，会将该模本的相关批次的样本拿来与胎儿进行检查--这里只做错配为5的结果
        :return:
        """
        self.ms_match = self.add_module("medical.paternity_test_v3.ms_match")
        self.ms_match.set_options({
            "dad_tab": self.output_dir + '/' + self.option("dad_id") + '.tab',
            "mom_tab": self.output_dir + '/' + self.option("mom_id") + '.tab',
            "preg_tab": self.output_dir + '/' + self.ms_son_id + '.tab'
            if self.son_is_dcpt else self.output_dir + '/' + self.option("son_id") + '.tab',
            "ref_point": self.ref_point,
            "err_min": 5
        })
        self.ms_match.on("end", self.set_output, "ms_match")
        self.ms_match.on("end", self.end)
        self.ms_match.run()

    def get_dad_no_analysis(self):
        """
        获取需要查重的样品：
        1、4个月内的F样品
        2、批次中的母本
        :return:
        """
        ref_list = self.api_pt.get_ref_list()                              # 查询四个月内的所有样品
        ref_list = [i for i in ref_list if "-F" in i]                      # 筛选出父本ID
        batch_list = self.api_pt.get_all_batch_sample(self.option("dad_id"), self.option("mom_id"),
                                                      self.option("son_id"))
        ref_list.extend(batch_list)
        ref_list = list(set(ref_list))
        # self.logger.info("FMS所有的批次样本列表：{}".format(ref_list))
        father_dedup_id = self.api_pt.find_dedup_err_id(5, self.option("mom_id"), self.ms_son_id)
        if father_dedup_id:
            dad_list_ = set(self.api_pt.find_dedup_father_id(father_dedup_id))   # 已经进行查重的样本，这里包含了FM
            for m in ref_list:
                if m not in dad_list_ and os.path.exists(os.path.join(self.ref_data, "{}.tab".format(m))):
                    self.dad_list.append(m)
        else:
            for m in ref_list:
                if os.path.exists(os.path.join(self.ref_data, "{}.tab".format(m))):
                    self.dad_list.append(m)
            self.logger.info("库中未查找到该M-S的查重结果，后面将进行全库查重分析！")
        if len(self.dad_list) == 0:
            self.logger.info("查重父本列表为空，将跳过查重分析模块！")
        # self.logger.info(self.dad_list)

    def dedup_run(self):
        """
        引进新的查重机制，首先检测出所有没有进行家系分析的父本，然后导表的时候要将样本的批次信息导入
        dad_list: "WQ123-F1,WQ1234-F1,WQ5678-F1"
        :return:
        """
        n = 0
        for p in self.err_list:
            pt_analysis_dedup = self.add_tool("medical.paternity_test_v3.dedup_v2")
            self.step.add_steps('dedup_{}'.format(n))
            pt_analysis_dedup.set_options({
                "dad_id": self.option("dad_id"),
                "mom_tab": self.output_dir + '/' + self.option("mom_id") + '.tab',
                "preg_tab": self.output_dir + '/' + self.ms_son_id + '.tab' if self.son_is_dcpt else
                self.output_dir + '/' + self.option("son_id") + '.tab',
                "ref_point": self.ref_point,
                "err_min": p,
                "dad_list": ','.join(self.dad_list),
                "father_path": self.ref_data
            })
            step = getattr(self.step, 'dedup_{}'.format(n))
            step.start()
            pt_analysis_dedup.on('end', self.finish_update, 'dedup_{}'.format(n))
            self.tools_dedup.append(pt_analysis_dedup)
            n += 1
        for j in range(len(self.tools_dedup)):
            self.tools_dedup[j].on('end', self.set_output, 'dedup')
        if len(self.tools_dedup) > 1:
            self.on_rely(self.tools_dedup, self.sample_similarity_run)
        else:
            self.tools_dedup[0].on("end", self.sample_similarity_run)
        for tool in self.tools_dedup:
            gevent.sleep(1)
            tool.run()

    def sample_similarity_run(self):
        """个体识别"""
        for id_ in [self.option("dad_id"), self.option('mom_id')]:
            similarity = self.add_tool("medical.paternity_test_v3.sample_similarity")
            similarity.set_options({
                "sample_id": id_
            })
            self.similarity_tools.append(similarity)
        if self.similarity_tools:
            if len(self.similarity_tools) > 1:
                self.on_rely(self.similarity_tools, self.end)
            else:
                self.similarity_tools[0].on("end", self.end)
        else:
            raise Exception("similarity_tools列表为空！")
        for tool in self.similarity_tools:
            gevent.sleep(1)
            tool.run()

    def submit_family_analysis_run(self):
        """
        根据MS降噪后的tab文件，检索这个家系中的其余的父本与母本，然后重新组建家系运行亲子鉴定流程--这里只有第一次运行的时候才会回调
        ---这一需求暂时不用
        :return:
        """
        familys = self.api_pt.make_familys(self.option("dad_id").split('-')[0])
        if familys and self.run_first and self.son_is_dcpt:
            for family in familys:
                if family[0] == self.option("dad_id") and family[1] == self.option('mom_id'):
                    pass
                else:
                    submit_tool = self.add_tool("medical.datasplit_v3.submit_work")
                    submit_tool.set_options({
                        "api_name": "s/med/pt_v3/pt_analysis",
                        "params": json.dumps({
                            'batch_id': self.option('batch_id'),
                            # 'board_batch': self.board_batch,
                            'member_id': self.option('member_id'),
                            'mom_id': family[1],
                            'dad_id': family[0],
                            'son_id': self.ms_son_id,
                            'fastq_path': self.option('fastq_path').prop['path'],
                            'err_min_num': '11',
                            'old_mom_id': self.option("mom_id")
                        })
                    })
                    self.other_family_workflow.append(submit_tool)
            if len(self.other_family_workflow) > 1:
                self.on_rely(self.other_family_workflow, self.end)
            else:
                self.other_family_workflow[0].on("end", self.end)
            for tool in self.other_family_workflow:
                gevent.sleep(1)
                tool.run()
        else:
            self.logger.info("没有新的家系组建成功，流程直接结束！")
            self.end()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_output(self, event):
        obj = event["bind_object"]
        if event['data'] == "fastq2tab":
            if os.path.exists(obj.output_dir + '/bam2tab'):
                self.linkdir(obj.output_dir + '/bam2tab', self.output_dir)
            elif os.path.exists(obj.output_dir + '/bam2tab_dc'):
                self.linkdir(obj.output_dir + '/bam2tab_dc', self.output_dir)
        if event['data'] == "pt_analysis":
            self.linkdir(obj.output_dir + '/family_analysis', self.output_dir + '/pt_result_' +
                         str(obj.option('err_min')))
            self.linkdir(obj.output_dir + '/family_merge', self.output_dir + '/pt_result_' +
                         str(obj.option('err_min')))
            self.linkdir(obj.output_dir + '/result_info', self.output_dir + '/pt_result_' +
                         str(obj.option('err_min')))
        if event['data'] == "pt_analysis_homohybrid":
            self.linkdir(obj.output_dir + '/family_analysis', self.output_dir + '/pt_homohybrid_result_' +
                         str(obj.option('err_min')))
            self.linkdir(obj.output_dir + '/family_merge', self.output_dir + '/pt_homohybrid_result_' +
                         str(obj.option('err_min')))
            self.linkdir(obj.output_dir + '/result_info', self.output_dir + '/pt_homohybrid_result_' +
                         str(obj.option('err_min')))
        if event['data'] == "dedup":
            self.linkdir(obj.output_dir, self.output_dir + '/pt_result_' + str(obj.option('err_min')))
        if event['data'] == "dp_correction":
            self.linkdir(obj.output_dir, self.output_dir)

    def check_tab(self):
        """
        用于胎儿call snp完成后，用于实时检查下该家系的父本与母本tab文件是否存在，如果不存在就终止流程，存在后，
        再检查下是否这个家系已经分析过了，如果检查都通过，要先将父本，母本与胎儿tab文件重新导出来到output目录中
        :return:
        """
        params_json = {
            'dad_id': self.option("dad_id"),
            'mom_id': self.option("mom_id"),
            'son_id': self.ms_son_id
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        if not self.api_pt.tab_exist(self.option("dad_id")):
            if self.run_first:
                self.api_pt.add_sample_pt([self.option("son_id")])
            self.logger.info('该家系中胎儿call snp完成，'
                             '由于父本{}的tab文件不存在，流程终止--可以忽略！'.format(self.option("dad_id")))
            return False
        elif not self.api_pt.tab_exist(self.option("mom_id")):
            if self.run_first:
                self.api_pt.add_sample_pt([self.option("son_id")])
            self.logger.info('该家系中胎儿call snp完成，由于母本{}的tab文件不存在，'
                             '流程终止--可以忽略！'.format(self.option("mom_id")))
            return False
        if self.api_pt.find_params_result(params):
            if self.run_first:
                self.api_pt.analysis_status(self.option("dad_id"), "3")  # 更新sg_family表中的analysis_status为3
            self.logger.info('该家系以前已经分析过了，不再进行重复分析，流程终止--可以忽略！')
            return False
        self.export_tab(self.option("dad_id"))
        self.export_tab(self.option("mom_id"))
        if self.run_first:
            self.export_tab(self.option("son_id"))
        else:
            self.export_tab(self.option("son_id").split('_')[0])
        return True

    def export_tab(self, sample_id):
        """
        判断文件是否存在，不存在的时候，导出对应文件
        :param sample_id:
        :return:
        """
        if not os.path.exists(self.output_dir + '/' + sample_id + '.tab'):
            self.api_pt.export_tab_file(sample_id, self.output_dir)
        else:
            self.logger.info("样本{}tab文件已经存在！".format(sample_id))

    def get_sample_no_analysis(self, samples_list):
        """
         用于获取option中的samples中没有进行过call snp的样本列表用于后面的call snp分析
        :param samples_list:
        :return:
        """
        samples = []
        for m in samples_list:
            if not self.api_pt.tab_exist(m):
                samples.append(m)
            else:
                self.logger.info("样本{}已经call snp完成了，不在进行重复call snp".format(m))
                self.api_pt.update_analysis_status(self.option("batch_id"), "snp", "true")
        return samples

    def get_father_id(self):
        """
        根据update_info获取到sg_father主表id，进行家系分析的update_info格式如下：
        {"5a76c030a4e1af070b53b0a9": "sg_father", "batch_id": "5a684bdc89dde4200e6b7f34"}
        :return:
        """
        for i in json.loads(self.option('update_info')):
            if i == "batch_id":
                self.batch_id = json.loads(self.option('update_info'))[i]
                continue
            self.father_id = i
        self.logger.info("self.father_id 为{}, self.batch_id 为{}".format(self.father_id, self.batch_id))

    def run(self):
        """
        1)samples列表存在，则将所有的样本进行call snp
            1.1 call snp完成后，更新sg_analysis_status主表status为end
            1.2 每个样本call snp完成后，更新下sg_analysis_status中的end_count个数加1
        2)samples列表不存在
            2.1 判断old_mom_id是否存在，
                2.1.1 如果存在的话，就执行M-S降噪后(那old_mom_id与胎儿进行)的结果进行家系分析
                2.1.2 如果不存在的话，就去判断，胎儿时多重还是杂补的
                    2.1.2.1 如果是杂补的，判断，如果tab文件存在就进行家系分析，不存在的时候，先call snp 然后在进行家系分析
                    2.1.2.2 如果是多重的，先M-S降噪，然后进行家系分析
                2.1.3 家系分析完成后，检查胎儿与母本是否匹配，
                    2.1.3.1 匹配的话就直接进行后续查重与个体识别，以及回调，进行其它家系分析
                    2.1.3.2 不匹配的话，就运行，将该母本同一批次的所有样本，进行胎儿与母本是否匹配计算，看能够找到与胎儿匹配的母本
            2.2 传入的son的tab文件存在，检查父本与母本tab文件是否存在，存在的话就进行家系分析，父本母本都不存在就停止分析
            2.3 传入的son的tab文件不存在
                2.3.1 该son样本在该板子中继续，首先将son样本进行call snp，call snp成功后，在检查该家系的父本与母本tab文件
                是否都存在，如果都存在的话，检查该家系分析是否已经分析完成了（状态为end或者start）。
                2.3.2 该son样本不在该板子中终止流程。
        所有流程运行完之后就进行该样本的qc+批次信息的整合汇总
        :return:
        """
        samples = []
        self.board_batch, self.split_type = self.api_pt.get_board_name(self.option("batch_id"))  # 获取 拆分类型 PE or SE
        if self.option('samples'):
            samples = self.get_sample_no_analysis(json.loads(self.option('samples')))
            if len(samples) > 0:
                self.batch_id = json.loads(self.option('update_info')).keys()[0]
                self.fastq_call_snp_run(samples)
                if self.fastq_call_snp_modules:
                    if len(self.fastq_call_snp_modules) > 1:
                        self.on_rely(self.fastq_call_snp_modules, self.end)
                    else:
                        self.fastq_call_snp_modules[0].on('end', self.end)
                else:
                    raise Exception("fastq2tab_tools对象列表为空，无法进行call snp分析--流程正常终止！")
                for module in self.fastq_call_snp_modules:
                    gevent.sleep(1)
                    module.run()
            else:
                self.logger.info('该批次的样本都进行过了call snp分析，没有新样本可以进行分析！')
                self.api_pt.update_snp_status(self.option("batch_id"))
                self.end_by_workflow = False
                gevent.spawn_later(5, self.end)
        else:
            self.get_father_id()  # 初始化一下father_id值
            if self.option("old_mom_id"):
                self.run_first = False   # 不是第一次运行
                self.ms_son_id = self.option("son_id")  # M S 合并后的样本id
            else:
                if self.api_pt.get_sample_analysis_type(self.option("son_id")) != 'wqcf':  # 初始化一下胎儿是否是多重的
                    self.son_is_dcpt = False
                    self.ms_son_id = self.option("son_id")
                else:
                    self.ms_son_id = "{}_{}".format(self.option("son_id"), self.find_sample_id(self.option("mom_id")))
            self.logger.info("ms_son_id:{}".format(self.ms_son_id))
            self.api_pt.analysis_status(self.option("dad_id"), "2")   # 更新sg_family表中的analysis_status为2
            if self.api_pt.tab_exist(self.option("son_id")) or not self.run_first:
                self.logger.info("开始检查胎儿是否存在，或者是不是第一次运行！")
                if self.api_pt.tab_exist(self.option("dad_id")) and self.api_pt.tab_exist(self.option("mom_id")):
                    if self.son_is_dcpt:
                        self.dp_correction_run()
                    else:
                        self.family_analysis_run()
                else:
                    self.api_pt.analysis_status(self.option("dad_id"), "1")  # 该情况默认为缺样本analysis_status为1
                    self.logger.info('由于父本或者母本的tab文件不存在，流程终止--可以忽略！')
                    self.end_by_workflow = False
                    gevent.spawn_later(5, self.end)
            else:
                if self.api_pt.check_sample_in_board(self.option("batch_id"), self.option("son_id")):
                    samples.append(self.option("son_id"))
                    self.is_update = 'false'  # 在家系分析中 胎儿的call snp  不用更新状态
                    self.fastq_call_snp_run(samples)
                    if self.son_is_dcpt:
                        self.fastq_call_snp_modules[0].on('end', self.dp_correction_run)
                    else:
                        self.fastq_call_snp_modules[0].on('end', self.family_analysis_run)
                    self.fastq_call_snp_modules[0].run()
                else:
                    self.api_pt.analysis_status(self.option("dad_id"), "1")  # 该情况默认为缺样本analysis_status为1
                    self.logger.info('由于样本{}不在该板子中，无法家系分析，'
                                     '流程终止--可以忽略！'.format(self.option("son_id")))
                    self.end_by_workflow = False
                    gevent.spawn_later(5, self.end)
        super(PtAnalysisWorkflow, self).run()

    def set_db(self):
        """结果文件批量入库"""
        api_pt = self.api.api('medical.paternity_test_v3.paternity_test_v3')
        result_name = self.option("dad_id") + "_" + self.option("mom_id") + "_" + self.ms_son_id
        if self.msismatch:
            # 胎儿杂合导表部分
            for i in self.err_list:
                dir_path = self.output_dir + '/pt_result_' + str(i)
                if not os.path.exists(dir_path):
                    continue
                name = self.option("dad_id") + "_" + self.option("mom_id") + "_" + self.ms_son_id + "_error-" + str(i)
                self.father_err_id = api_pt.add_sg_father_err(self.option("batch_id"), self.option("family_id"),
                                                              self.father_id, i, "1", name, self.option('member_id'))
                results = os.listdir(dir_path)
                for f in results:
                    if f == result_name + "_family_joined_tab.txt":
                        api_pt.add_sg_pt_father_detail(dir_path + '/' + f, self.father_err_id)
                    elif f == self.option("mom_id") + "_" + self.ms_son_id + '_info_show.txt':
                        api_pt.add_info_detail(dir_path + '/' + f, self.father_err_id, self.option("dad_id"),
                                               self.option("mom_id"), self.ms_son_id)
                    elif f == result_name + '.txt':
                        results = api_pt.find_dedup_err_id(i, self.option("mom_id"), self.ms_son_id)
                        if results:
                            father_dedup_id = results
                        else:
                            father_dedup_id = api_pt.add_dedup_main(i, self.option("mom_id"), self.ms_son_id)
                        api_pt.update_to_father_err(father_dedup_id, self.father_err_id)
                        api_pt.import_dedup_data(dir_path + '/' + f, father_dedup_id)
                    elif f == result_name + "_family.png":
                        if self.config.RGW_ENABLE:
                            pic_path = "{}/pt_result_pic/".format(self.bucket) + \
                                       self.father_id + '/pt_result_' + str(i) + '/'
                        else:
                            pic_path = "rerewrweset/MEDfiles/pt_result_pic/" + \
                                       self.father_id + '/pt_result_' + str(i) + '/'
                        api_pt.add_result_pic(pic_path, result_name, self.father_id, self.father_err_id)
                    elif f == result_name + '_test_pos.txt':
                        api_pt.add_test_pos(dir_path + '/' + f, self.father_err_id)
                    elif f == result_name + "_family_analysis.txt":
                        api_pt.add_analysis_tab(dir_path + '/' + f, self.father_err_id, self.father_id, i)
                if len(self.dad_list) == 0:
                    father_dedup_id_ = api_pt.find_dedup_err_id(i, self.option("mom_id"), self.ms_son_id)
                    if father_dedup_id_:
                        api_pt.update_to_father_err(father_dedup_id_, self.father_err_id)
                api_pt.update_father_err_status(self.father_err_id)
            # 胎儿纯合与杂合导表部分
            for i in self.err_list_homohybrid:
                dir_path = self.output_dir + '/pt_homohybrid_result_' + str(i)
                if not os.path.exists(dir_path):
                    continue
                name = self.option("dad_id") + "_" + self.option("mom_id") + "_" + self.ms_son_id + "_error-" + str(i)
                self.father_err_id_homohybrid = api_pt.add_sg_father_err(self.option("batch_id"),
                                                                         self.option("family_id"), self.father_id, i,
                                                                         "1", name, self.option('member_id'),
                                                                         postype='homohybrid')
                results = os.listdir(dir_path)
                for f in results:
                    if f == result_name + "_family_match_update.png":
                        if self.config.RGW_ENABLE:
                            pic_path = "{}/pt_result_pic/".format(self.bucket) + self.father_id +\
                                       '/pt_homohybrid_result_' + str(i) + '/'
                        else:
                            pic_path = "rerewrweset/MEDfiles/pt_result_pic/" + \
                                       self.father_id + '/pt_homohybrid_result_' + str(i) + '/'
                        api_pt.add_result_pic(pic_path, result_name, self.father_id, self.father_err_id_homohybrid,
                                              postype='homohybrid')
                    elif f == result_name + '_test_pos_update.txt':
                        api_pt.add_test_pos(dir_path + '/' + f, self.father_err_id_homohybrid, postype='homohybrid')
                    elif f == result_name + "_family_analysis.txt":
                        api_pt.add_analysis_tab(dir_path + '/' + f, self.father_err_id_homohybrid, self.father_id, i,
                                                postype='homohybrid')
                    elif f == result_name + "_num_pos_map.txt":
                        api_pt.add_num_pos_map(dir_path + '/' + f, self.father_err_id_homohybrid, self.father_id)
                api_pt.update_father_err_status(self.father_err_id_homohybrid, postype='homohybrid')
        else:
            api_pt.update_db_record("sg_father", {"_id": ObjectId(self.father_id)}, {"ms_match": "NO"})  # MS不匹配
            api_pt.add_sg_father_msmatch(self.ms_match.output_dir + "/ms_match.txt", self.father_id, self.ms_son_id)
        api_pt.update_father_tab_batch(self.father_id, self.option("dad_id"))  # 更新父本的建库批次与抽提批次
        api_pt.update_report_status(self.option("dad_id"))  # 用于更新family表中报告状态为2
        for m in self.err_list:
            api_pt.update_dedup_father(m, self.option("mom_id"), self.ms_son_id, self.father_id)  # 用于更新查重
            # 批次到father表,最后面失败不会影响整个流程
        # api_pt.update_family_analysis_process(self.option("batch_id"))  # 页面用于更新整个家系分析进度
        # api_pt.add_sg_family_info(self.father_id, self.option("dad_id"), self.option("mom_id"))
        # 用于出报告的数据汇总，以及后期他们想看家系的出报告的情况
        self.api_pt.analysis_status(self.option("dad_id"), "3")  # 更新sg_family表中的analysis_status为3

    def end(self):
        if self.end_by_workflow:
            if not self.option('samples'):
                self.set_db()
                self.api_pt.add_sample_pt([self.option("dad_id"), self.option("mom_id"), self.option("son_id")])
            else:
                self.api_pt.add_sample_pt(json.loads(self.option('samples')))
        else:
            if not self.option('samples'):
                self.api_pt.update_db_record("sg_father", {"_id": ObjectId(self.father_id)}, {"is_show": "2"})
            # 要将sg_father中is_show更新为2 页面家系中只展示is_show为1的记录
        self.upload_target_path()
        super(PtAnalysisWorkflow, self).end()

    def upload_target_path(self):
        """
        将结果link到指定磁盘路径下
        :return:
        """
        self.logger.info("开始设置结果目录！")
        if self.config.RGW_ENABLE:
            transfer = MultiFileTransfer()
            transfer.add_upload(self.output_dir + "/", "{}/pt_result_pic/{}/".format(self.bucket, self.father_id),
                                base_path=self.output_dir + "/")
            transfer.perform()
        else:
            if self._sheet.client == 'client01':
                self.path = '/mnt/ilustre/data/rerewrweset/MEDfiles/pt_result_pic/{}'.format(self.father_id)
            else:
                self.path = '/mnt/ilustre/tsanger-data/rerewrweset/MEDfiles/pt_result_pic/{}'.format(self.father_id)
            if not os.path.exists(self.path):
                os.makedirs(self.path)
            self.linkdir(self.output_dir, self.path)
        self.logger.info("设置结果目录成功！")
