# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
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
        last modified by hongdong@20171204
        :param wsheet_object:
        """
        self._sheet = wsheet_object
        super(PtAnalysisWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fastq_path", "type": "infile", "format": "sequence.fastq_dir"},  # fastq所在路径(文件夹
            {"name": "err_min_num", "type": "int", "default": 9},  # 该选项用于要循环多少个错配，默认是2-8
            {"name": "dad_id", "type": "string"},
            {"name": "mom_id", "type": "string"},
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
        self.tools_result = []
        self.tools_dedup = []
        self.similarity_tools = []
        self.rdata = []
        self.path = None  # 用于存储报告中的图片
        self.father_err_id = None
        self.dad_list = []
        self.board_batch = None
        self.api_pt = self.api.api('medical.paternity_test_v2')
        self.father_id = ''
        # self.father_id = json.loads(self.option('update_info')).keys()[0]
        self.ref_data = self.config.SOFTWARE_DIR + "/database/human/pt_ref/tab_data/"
        self.targets_bedfile = self.config.SOFTWARE_DIR + "/database/human/pt_ref/snp.chr.sort.3.bed"  # 3000个SNP位点
        self.ref_fasta = self.config.SOFTWARE_DIR + "/database/human/hg38.chromosomal_assembly/ref.fa"  # 参考基因组
        self.ref_point = self.config.SOFTWARE_DIR + "/database/human/pt_ref/targets.bed.rda"  # 每个位点的基因型概率
        self.is_update = "true"  # 在进度条更新的过程中用于区分是胎儿样本是家系分析中进行call snp 还是单独流程中进行call snp的
        self.batch_id = ''

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
        # fastq_dir = os.listdir(self.option('fastq_path').prop['path'])
        n = 0
        for sample_name in samples:
            # if "{}_R1.fastq.gz".format(sample_name) not in fastq_dir:
            #     raise Exception("样本：{}不在该下机板子中！".format(sample_name))
            fastq_call_snp = self.add_module("medical.paternity_test_v2.fastq_call_snp")
            self.step.add_steps('fastq2tab{}'.format(n))
            fastq_call_snp.set_options({
                "sample_id": sample_name,
                "fastq_path": self.option("fastq_path"),
                "ref_fasta": self.ref_fasta,
                "targets_bedfile": self.targets_bedfile,
                "batch_id": self.batch_id,
                # "batch_id": self.option('batch_id'),
                "analysis_type": "pt" if '-S' in sample_name else "dcpt",
                "board_batch": self.board_batch,
                "is_update": self.is_update  # "true" or "False"
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
        self.check_tab()   # 用于检测胎儿call snp完成后，父本母本是否存在
        self.get_dad_no_analysis()   # 用于检测后面是否需要查重
        n = 0
        for m in range(2, self.option('err_min_num')):
            family_analysis = self.add_module("medical.paternity_test_v2.family_analysis")
            self.step.add_steps('family_analysis{}'.format(n))
            result_dir = os.path.join(self.output_dir, 'pt_result_' + str(m))
            if not os.path.exists(result_dir):
                os.mkdir(result_dir)
            family_analysis.set_options({
                "dad_tab": self.output_dir + '/' + self.option("dad_id") + '.tab',
                "mom_tab": self.output_dir + '/' + self.option("mom_id") + '.tab',
                "preg_tab": self.output_dir + '/' + self.option("son_id") + '.tab',
                "ref_point": self.ref_point,
                "err_min": m
            })
            step = getattr(self.step, 'family_analysis{}'.format(n))
            step.start()
            family_analysis.on('end', self.finish_update, 'family_analysis{}'.format(n))
            self.modules_analysis.append(family_analysis)
            n += 1
        self.logger.info(self.modules_analysis)
        for j in range(len(self.modules_analysis)):
            self.modules_analysis[j].on('end', self.set_output, 'pt_analysis')
        if self.modules_analysis:
            if len(self.modules_analysis) > 1:
                self.on_rely(self.modules_analysis, self.result_info_run)
            elif len(self.modules_analysis) == 1:
                self.modules_analysis[0].on('end', self.result_info_run)
        else:
            raise Exception("modules_analysis列表为空！")
        for module in self.modules_analysis:
            gevent.sleep(1)
            module.run()

    def result_info_run(self):
        """画出相关结果图SNP分型图等"""
        n = 0
        for l in range(2, self.option('err_min_num')):
            result_dir = os.path.join(self.output_dir, 'pt_result_' + str(l))
            results = os.listdir(result_dir)
            self.logger.info(results)
            for f in results:
                if re.match(r'.*family_joined_tab\.Rdata$', f):
                    result_info = self.add_tool("medical.paternity_test_v2.result_info")
                    self.step.add_steps('result_info{}'.format(n))
                    self.rdata = os.path.join(result_dir, f)
                    result_info.set_options({
                        "tab_merged": self.rdata
                    })
                    step = getattr(self.step, 'result_info{}'.format(n))
                    step.start()
                    self.tools_result.append(result_info)
                    n += 1
                else:
                    pass
        for j in range(len(self.tools_result)):
            self.tools_result[j].on('end', self.set_output, 'result_info')
        if self.tools_result:
            if len(self.tools_result) > 1:
                if len(self.dad_list) == 0:
                    self.logger.info("111111111--")
                    self.on_rely(self.tools_result, self.end)   # 当需要查重的父本列表为空的时候，直接跳过查重部分
                else:
                    self.logger.info("111111111---")
                    self.on_rely(self.tools_result, self.dedup_run)
            elif len(self.tools_result) == 1:
                if len(self.dad_list) == 0:
                    self.tools_result[0].on('end', self.end)   # 当需要查重的父本列表为空的时候，直接跳过查重部分
                    self.logger.info("222222222--")
                else:
                    self.tools_result[0].on('end', self.dedup_run)
                    self.logger.info("222222222---")
        else:
            raise Exception("tools_result列表为空！")
        for tool in self.tools_result:
            gevent.sleep(1)
            tool.run()

    def get_dad_no_analysis(self):
        """
        获取没有进行过查重分析的父本
        :return:
        """
        ref_list = self.get_ref_list()
        father_dedup_id = self.api_pt.find_dedup_err_id(2, self.option("mom_id"), self.option("son_id"))
        if father_dedup_id:
            dad_list_ = self.api_pt.find_dedup_father_id(father_dedup_id)
            for m in ref_list:
                if m not in list(set(dad_list_)):
                    self.dad_list.append(m)
        else:
            self.dad_list = ref_list
            self.logger.info("库中未查找到该M-S的查重结果，后面将进行全库查重分析！")
        if len(self.dad_list) == 0:
            self.logger.info("查重父本列表为空，将跳过查重分析模块！")
        self.logger.info(self.dad_list)

    def get_ref_list(self):
        """
        用于获取查重参考库中父本的所有的样本id
        :return:
        """
        ref_list = []
        ref_dad = os.listdir(self.ref_data)
        for m in ref_dad:
            ref_list.append(m.strip().split('.')[0])
        return ref_list

    def dedup_run(self):
        """
        引进新的查重机制，首先检测出所有没有进行家系分析的父本，然后导表的时候要将样本的批次信息导入
        dad_list: "WQ123-F1,WQ1234-F1,WQ5678-F1"
        :return:
        """
        n = 0
        for p in range(2, self.option('err_min_num')):
            pt_analysis_dedup = self.add_tool("medical.paternity_test_v2.dedup_v2")
            self.step.add_steps('dedup_{}'.format(n))
            pt_analysis_dedup.set_options({
                "dad_id": self.option("dad_id"),
                "mom_tab": self.output_dir + '/' + self.option("mom_id") + '.tab',
                "preg_tab": self.output_dir + '/' + self.option("son_id") + '.tab',
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
            similarity = self.add_tool("medical.paternity_test_v2.sample_similarity")
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
            tool.run()

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
        if event['data'] == "result_info":
            dir_path = os.path.dirname(obj.option("tab_merged").prop['path'])
            if self._sheet.client == 'client01':
                self.path = '/mnt/ilustre/data/rerewrweset/MEDfiles/pt_result_pic/' + self.father_id
            else:
                self.path = '/mnt/ilustre/tsanger-data/rerewrweset/MEDfiles/pt_result_pic/' + self.father_id
            if not os.path.exists(self.path):
                os.makedirs(self.path)
            self.linkdir(obj.output_dir, dir_path)
            target_dir = os.path.basename(dir_path.rstrip())
            self.linkdir(obj.output_dir, self.path + '/' + target_dir)
        if event['data'] == "dedup":
            self.linkdir(obj.output_dir, self.output_dir + '/pt_result_' + str(obj.option('err_min')))

    def check_tab(self):
        """
        用于胎儿call snp完成后，用于实时检查下该家系的父本与母本tab文件是否存在，如果不存在就终止流程，存在后，
        再检查下是否这个家系已经分析过了，如果检查都通过，要先将父本，母本与胎儿tab文件重新导出来到output目录中
        :return:
        """
        params_json = {
            'dad_id': self.option("dad_id"),
            'mom_id': self.option("mom_id"),
            'son_id': self.option("son_id")
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        if not self.api_pt.tab_exist(self.option("dad_id")):
            self.api_pt.add_sample_pt([self.option("son_id")])
            self.logger.info('该家系中胎儿call snp完成，'
                             '由于父本{}的tab文件不存在，流程终止--可以忽略！'.format(self.option("dad_id")))
            super(PtAnalysisWorkflow, self).end()
        elif not self.api_pt.tab_exist(self.option("mom_id")):
            self.api_pt.add_sample_pt([self.option("son_id")])
            self.logger.info('该家系中胎儿call snp完成，由于母本{}的tab文件不存在，'
                             '流程终止--可以忽略！'.format(self.option("mom_id")))
            super(PtAnalysisWorkflow, self).end()
        if self.api_pt.find_params_result(params):
            self.api_pt.analysis_status(self.option("dad_id"), "3")  # 更新sg_family表中的analysis_status为3
            self.logger.info('该家系以前已经分析过了，不再进行重复分析，流程终止--可以忽略！')
            super(PtAnalysisWorkflow, self).end()
        if not os.path.exists(self.output_dir + '/' + self.option("dad_id") + '.tab'):
            self.api_pt.export_tab_file(self.option("dad_id"), self.output_dir)
        else:
            self.logger.info("样本{}tab文件已经存在！".format(self.option("dad_id")))
        if not os.path.exists(self.output_dir + '/' + self.option("mom_id") + '.tab'):
            self.api_pt.export_tab_file(self.option("mom_id"), self.output_dir)
        else:
            self.logger.info("样本{}tab文件已经存在！".format(self.option("mom_id")))
        if not os.path.exists(self.output_dir + '/' + self.option("son_id") + '.tab'):
            self.api_pt.export_tab_file(self.option("son_id"), self.output_dir)
        else:
            self.logger.info("样本{}tab文件已经存在！".format(self.option("son_id")))

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
            self.logger.info("self.father_id 为{}, self.batch_id 为{}".format(i, self.batch_id))

    def run(self):
        """
        1)samples列表存在，则将所有的样本进行call snp
            1.1 call snp完成后，更新sg_analysis_status主表status为end
            1.2 每个样本call snp完成后，更新下sg_analysis_status中的end_count个数加1
        2)samples列表不存在
            2.1 传入的son的tab文件存在，检查父本与母本tab文件是否存在，存在的话就进行家系分析，父本母本都不存在就停止分析
            2.2 传入的son的tab文件不存在
                2.2.1 该son样本在该板子中继续，首先将son样本进行call snp，call snp成功后，在检查该家系的父本与母本tab文件
                是否都存在，如果都存在的话，检查该家系分析是否已经分析完成了（状态为end或者start）。
                2.2.2 该son样本不在该板子中终止流程。
        所有流程运行完之后就进行该样本的qc+批次信息的整合汇总
        :return:
        """
        samples = []
        self.board_batch = self.api_pt.get_board_name(self.option("batch_id"))
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
                super(PtAnalysisWorkflow, self).run()
            else:
                self.logger.info('该批次的样本都进行过了call snp分析，没有新样本可以进行分析！')
                self.api_pt.update_snp_status(self.option("batch_id"))
                self.start_listener()
                self.fire("start")
                self.end()
        else:
            self.get_father_id()  # 初始化一下father_id值
            self.api_pt.analysis_status(self.option("dad_id"), "2")   # 更新sg_family表中的analysis_status为2
            if self.api_pt.tab_exist(self.option("son_id")):
                if self.api_pt.tab_exist(self.option("dad_id")) and self.api_pt.tab_exist(self.option("mom_id")):
                    self.family_analysis_run()
                else:
                    self.logger.info('由于父本或者母本的tab文件不存在，流程终止--可以忽略！')
                    super(PtAnalysisWorkflow, self).end()
            else:
                if self.api_pt.check_sample_in_board(self.option("batch_id"), self.option("son_id")):
                    samples.append(self.option("son_id"))
                    self.is_update = 'false'  # 在家系分析中 胎儿的call snp  不用更新状态
                    self.fastq_call_snp_run(samples)
                    if self.fastq_call_snp_modules:
                        if len(self.fastq_call_snp_modules) > 1:
                            self.on_rely(self.fastq_call_snp_modules, self.family_analysis_run)
                        else:
                            self.fastq_call_snp_modules[0].on('end', self.family_analysis_run)
                    else:
                        raise Exception("fastq2tab_tools对象列表为空，无法进行call snp分析--流程正常终止！")
                    for module in self.fastq_call_snp_modules:
                        module.run()
                else:
                    self.api_pt.analysis_status(self.option("dad_id"), "1")  # 该情况默认为缺样本analysis_status为1
                    self.logger.info('由于样本{}不在该板子中，无法家系分析，'
                                     '流程终止--可以忽略！'.format(self.option("son_id")))
                    super(PtAnalysisWorkflow, self).end()
            super(PtAnalysisWorkflow, self).run()

    def set_db(self):
        """结果文件批量入库"""
        api_pt = self.api.api('medical.paternity_test_v2')
        for i in range(2, self.option('err_min_num')):
            dir_path = self.output_dir + '/pt_result_' + str(i)
            if not os.path.exists(dir_path):
                continue
            name = self.option("dad_id") + "_" + self.option("mom_id") + "_" + self.option("son_id") + "_err" \
                                                                                                       "or-" + str(i)
            self.father_err_id = api_pt.add_sg_father_err(self.option("batch_id"), self.option("family_id"),
                                                          self.father_id, i, "1", name, self.option('member_id'))
            results = os.listdir(dir_path)
            result_name = self.option("dad_id") + "_" + self.option("mom_id") + "_" + self.option("son_id")
            for f in results:
                if f == result_name + "_family_joined_tab.txt":
                    api_pt.add_sg_pt_father_detail(dir_path + '/' + f, self.father_err_id)
                elif f == self.option("mom_id") + "_" + self.option("son_id") + '_info_show.txt':
                    api_pt.add_info_detail(dir_path + '/' + f, self.father_err_id)
                elif f == result_name + '.txt':
                    results = api_pt.find_dedup_err_id(i, self.option("mom_id"), self.option("son_id"))
                    if results:
                        father_dedup_id = results
                    else:
                        father_dedup_id = api_pt.add_dedup_main(i, self.option("mom_id"), self.option("son_id"))
                    api_pt.update_to_father_err(father_dedup_id, self.father_err_id)
                    api_pt.import_dedup_data(dir_path + '/' + f, father_dedup_id)
                elif f == result_name + "_family.png":
                    api_pt.add_result_pic("pt_result_pic/" + self.father_id + '/pt_result_' + str(i) + '/', result_name,
                                          self.father_id, self.father_err_id)
                elif f == result_name + '_test_pos.txt':
                    api_pt.add_test_pos(dir_path + '/' + f, self.father_err_id)
                elif f == result_name + "_family_analysis.txt":
                    api_pt.add_analysis_tab(dir_path + '/' + f, self.father_err_id, self.father_id, i)
            if len(self.dad_list) == 0:
                father_dedup_id_ = api_pt.find_dedup_err_id(i, self.option("mom_id"), self.option("son_id"))
                if father_dedup_id_:
                    api_pt.update_to_father_err(father_dedup_id_, self.father_err_id)
            api_pt.update_father_err_status(self.father_err_id)
        api_pt.update_father_tab_batch(self.father_id, self.option("dad_id"))  # 更新父本的建库批次与抽提批次
        api_pt.update_report_status(self.option("dad_id"))  # 用于更新family表中报告状态为2
        for m in range(2, self.option('err_min_num')):
            api_pt.update_dedup_father(m, self.option("mom_id"), self.option("son_id"), self.father_id)  # 用于更新查重
            # 批次到father表,最后面失败不会影响整个流程
        # api_pt.update_family_analysis_process(self.option("batch_id"))  # 页面用于更新整个家系分析进度
        api_pt.add_sg_family_info(self.father_id, self.option("dad_id"), self.option("mom_id"))
        # 用于出报告的数据汇总，以及后期他们想看家系的出报告的情况
        self.api_pt.analysis_status(self.option("dad_id"), "3")  # 更新sg_family表中的analysis_status为3

    def end(self):
        if not self.option('samples'):
            self.set_db()
            self.api_pt.add_sample_pt([self.option("dad_id"), self.option("mom_id"), self.option("son_id")])
        else:
            self.api_pt.add_sample_pt(json.loads(self.option('samples')))
        super(PtAnalysisWorkflow, self).end()
