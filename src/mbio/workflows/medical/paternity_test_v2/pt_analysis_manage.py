# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
from biocluster.core.exceptions import OptionError
from biocluster.workflow import Workflow
import os


class PtAnalysisManageWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        用途及运行逻辑：该工作流紧接着拆分流程，将构建出运行列表，加急等信息，生成家系列表，然后进行样本的call snp，
        以及运行family_submit tool，将所有的家系投递出去
        auther: hongdongxuan
        time: 20171126
        :param wsheet_object:
        """
        self._sheet = wsheet_object
        super(PtAnalysisManageWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fastq_path", "type": "infile", "format": "sequence.fastq_dir"},  # 拆分后的亲子鉴定的fastq文件路径
            {"name": "datasplit_id", "type": "string"},  # 每拆分就有一个batch_id，仅用于区分批次，也用于后面的进度统计
            {"name": "update_info", "type": "string"},  # 用于后面更新主表，该主表是在拆分流程中进行导表
            {"name": "member_id", "type": "string"},  # 会员ID
            {"name": "types", "type": "string"},
            {"name": "board_batch", "type": "string"}  # 板号名字
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.fastq2tab_tools = []
        self.family_list = []
        self.familysubmit = self.add_tool("medical.paternity_test_v2.family_submit")
        self.ref_data = self.config.SOFTWARE_DIR + "/database/human/pt_ref/tab_data/"
        self.targets_bedfile = self.config.SOFTWARE_DIR + "/database/human/pt_ref/snp.chr.sort.3.bed"  # 3000个SNP位点
        self.ref_fasta = self.config.SOFTWARE_DIR + "/database/human/hg38.chromosomal_assembly/ref.fa"  # 参考基因组
        self.ref_point = self.config.SOFTWARE_DIR + "/database/human/pt_ref/targets.bed.rda"  # 每个位点的基因型概率

    def check_options(self):
        """检查参数设置"""
        if not self.option('fastq_path').is_set:
            raise OptionError('必须提供fastq文件所在的路径')
        if not self.option("member_id"):
            raise OptionError('必须提供member_id!')
        return True

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def fastq2tab_run(self, samples):
        """
        用于激发样本的call snp module
        :param samples:
        :return:
        """
        n = 0
        for sample_name in samples:
            if "{}_R1.fastq.gz".format(sample_name) not in os.listdir(self.option('fastq_path').prop['path']):
                raise Exception("样本：{}不在{}下机板子中！".format(sample_name, self.option("board_batch")))
            fastq2tab = self.add_module("medical.paternity_test_v2.fastq_call_snp")
            self.step.add_steps('fastq2tab{}'.format(n))
            fastq2tab.set_options({
                "sample_id": sample_name,
                "fastq_path": self.option("fastq_path"),
                "cpu_number": 4,
                "ref_fasta": self.ref_fasta,
                "targets_bedfile": self.targets_bedfile,
                "batch_id": self.option('batch_id'),
                "type": "pt" if '-S' in sample_name else "dcpt"
            })
            step = getattr(self.step, 'fastq2tab{}'.format(n))
            step.start()
            fastq2tab.on('end', self.finish_update, 'fastq2tab{}'.format(n))
            self.fastq2tab_tools.append(fastq2tab)
            n += 1
        for j in range(len(self.fastq2tab_tools)):
            self.fastq2tab_tools[j].on('end', self.set_output, 'fastq2tab')
        if self.fastq2tab_tools:
            if len(self.fastq2tab_tools) > 1:
                self.on_rely(self.fastq2tab_tools, self.family_submit_run)
            else:
                self.fastq2tab_tools[0].on('end', self.family_submit_run)
        for tool in self.fastq2tab_tools:
            tool.run()

    def family_submit_run(self):
        """
        用于在父本与母本call snp结束后进行家系接口投递分析
        :return:
        """
        if len(self.family_list) == 0:
            self.logger.info("家系列表为空，不进行后面的家系分析，流程到此结束，bye！")
            self.end()
        self.familysubmit.set_option({
            "family_list": self.family_list,
            "datasplit_id": self.option("datasplit_id"),
            "types": self.option("types"),
            "member_id": self.option("member_id"),
            "fastq_path": self.option("fastq_path"),
            "err_min_num": 3
        })
        self.familysubmit.on('end', self.end)

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
        obj = event['bind_object']
        if event['data'] == 'fastq2tab':
            self.linkdir(obj.output_dir, 'fastq2tab')

    def run(self):
        """
        检查拆分的主表id查询的板号与self.option('board_name')是否一致，然后获取到board_name对应的所有的亲子的样本，
        然后去sg_sample_tab，检查下是否已经存在，存在说明已经分析过了，就不再进行分析，否则就添加到分析队列samples中，除此之外
        还要进行样本信息导入到sg_family表格中。当Samples为空，直接进行家系分析，
        组建好的家系列表family_list为空的时候就直接self.end, 不然就进行家系投递分析
        :return:
        """
        super(PtAnalysisManageWorkflow, self).run()
        pt_api = self.api.api("medical.paternity_test_v2")
        if pt_api.get_board_name(self.option("datasplit_id")) != self.option("board_batch"):
            raise Exception("该工作流中的拆分主表id与board_batch不一致！")
        samples = pt_api.get_sample_list(self.option("board_batch"))
        pt_api.add_sg_family(samples)  # 导入样本信息到sg_family表格中
        self.family_list = pt_api.get_family_list(samples)
        if len(samples) != 0:
            self.fastq2tab_run(samples)
        else:
            self.family_submit_run()

    def end(self):
        super(PtAnalysisManageWorkflow, self).end()
