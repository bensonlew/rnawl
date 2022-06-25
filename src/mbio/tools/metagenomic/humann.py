# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# __last_modify__ = '20191224'


from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil

class HumannAgent(Agent):
    """
    humann2基于reads进行丰度统计和功能计算
    """

    def __init__(self, parent):
        super(HumannAgent, self).__init__(parent)
        options = [
            {"name": "read1", "type": "infile", "format": "sequence.fastq"},
            {"name": "read2", "type": "infile", "format": "sequence.fastq"},
            {"name": "sample", "type": "string"},
            {"name": "search_mode", "type": "string", "default": "uniref50"}, #uniref50 or uniref90
            {"name": "prescreen_threshold", "type": 'float', "default": 0.01},
            {"name": "identity_threshold", "type": 'int', "default": 50},
            {"name": "subject_coverage", "type": 'int', "default": 50},
            {"name": "query_coverage", "type": 'int', "default": 90},
            {"name": "pathways", "type": "string", "default": "metacyc"},#metacyc,unipathway
            {"name": "translated_alignment", "type": "string", "default": "diamond"},  #usearch,rapsearch,diamond
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("read1").is_set:
            raise OptionError("请输入read1文件！")
        if not self.option("read2").is_set:
            raise OptionError("请输入read2文件！")
        if not self.option("sample"):
            raise OptionError("请输入sample的样品名称！")
        if self.option("search_mode") not in ["uniref50", "uniref90"]:
            raise OptionError("选择的基因数据库不是uniref50或者uniref90！")
        if self.option("pathways") not in ["metacyc", "unipathway"]:
            raise OptionError("选择的基因数据库不是metacyc或者unipathway！")
        if self.option("translated_alignment") not in ["usearch", "rapsearch", "diamond"]:
            raise OptionError("选择的比对方法不是diamond，rapsearch、usearch其中一个！")

    def set_resource(self):
        """
        所需资源
        """
        num = os.path.getsize(self.option("read1").prop['path'])/1000000000
        self._cpu = 10
        self._memory = str(num*10) + 'G'


class HumannTool(Tool):
    def __init__(self, config):
        super(HumannTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/bioinfo/meta/usearch-v7.0:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/RAPSearch2.24_64bits/bin:" + self.config.SOFTWARE_DIR + "/miniconda2/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/metaphlan2:" + self.config.SOFTWARE_DIR + "/bioinfo/align/diamond-0.8.35:"
        self.bowtie2 =self.config.SOFTWARE_DIR + "/bioinfo/align/bowtie2-2.3.4.3-linux-x86_64"
        self.set_environ(PATH=self.path)
        self.humann = "/miniconda2/bin/humann2"
        self.out =self.work_dir + "/result"
        self.mergefq = "/bioinfo/metaGenomic/sortmerna-2.1b/scripts/merge-paired-reads.sh"
        self.reads = self.work_dir + "/merge.fq"

    def merge_fq(self):
        cmd = "{} {} {} {}".format(self.mergefq, self.option('read1').prop['path'], self.option('read2').prop['path'], self.reads)
        command = self.add_command("merge_fq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("merge_fq运行完成！")
        else:
            self.set_error("merge_fq运行完成运行出错!")

    def run_humann2(self):
        """
        description
        :return:
        """
        cmd = "{} --threads 8 --input {} --output {} --search-mode {} --prescreen-threshold {} --identity-threshold {} --translated-subject-coverage-threshold {} --translated-query-coverage-threshold {} --pathways {} --translated-alignment {} --remove-temp-output".format(self.humann, self.reads, self.out, self.option("search_mode"), self.option("prescreen_threshold"), self.option("identity_threshold"), self.option("subject_coverage"), self.option("query_coverage"), self.option("pathways"), self.option("translated_alignment"))
        cmd +=" --bowtie2={}".format(self.bowtie2)
        command = self.add_command("run_humann2", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_humann2运行完成！")
        else:
            self.set_error("run_humann2运行完成运行出错!")


    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        if os.path.exists(self.output_dir + "/" + self.option("sample")):
            shutil.rmtree(self.output_dir + "/" + self.option("sample"))
        os.mkdir(self.output_dir + "/" + self.option("sample"))
        for i in ["merge_genefamilies.tsv","merge_pathabundance.tsv","merge_pathcoverage.tsv"]:
            os.link(self.out + "/" +i, self.output_dir + "/" + self.option("sample") + "/" + i)


    def run(self):
        super(HumannTool, self).run()
        self.merge_fq()
        self.run_humann2()
        self.set_output()
        self.end()