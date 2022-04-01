# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# __last_modify__ = '20191226'


from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil

class MetaphlanAgent(Agent):
    """
    metaphlan2基于reads进行物种分类以及丰度统计
    """

    def __init__(self, parent):
        super(MetaphlanAgent, self).__init__(parent)
        options = [
            {"name": "read1", "type": "infile", "format": "sequence.fastq"},
            {"name": "read2", "type": "infile", "format": "sequence.fastq"},
            {"name": "sample", "type": "string"},
            {"name": "read_min_len", "type": "int", "default": 70},
            {"name": "bt2_ps", "type": "string", "default": "very-sensitive"},#sensitive,very-sensitive,sensitive-local,very-sensitive-local
            {"name": "tax_lev", "type": "string", "default": "a"},#a,k,p,c,o,f,g,s
            {"name": "min_cu_len", "type": "int", "default": 2000},
            {"name": "stat", "type": "string", "default": "avg_g"},#avg_g,avg_l,tavg_g,tavg_l,wavg_g,wavg_l,med
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

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 8
        self._memory = '50G'


class MetaphlanTool(Tool):
    def __init__(self, config):
        super(MetaphlanTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/metaphlan2:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/metaphlan2/utils:" + self.config.SOFTWARE_DIR + "/bioinfo/align/diamond-0.8.35:"
        self.bowtie2 =self.config.SOFTWARE_DIR + "/bioinfo/align/bowtie2-2.3.4.3-linux-x86_64/bowtie2"
        self.set_environ(PATH=self.path)
        self.metaphlan = "/bioinfo/metaGenomic/metaphlan2/metaphlan2.py"
        self.out =self.work_dir + "/" + self.option("sample") + ".profiled.txt"
        self.path = self.option("read1").prop['path'] + "," + self.option("read2").prop['path']
        self.bowtie2out = self.work_dir + "/" + self.option("sample") + ".bowtie2out.txt"

    def run_metaphlan2(self):
        """
        description
        :return:
        """
        if os.path.exists(self.bowtie2out):
            os.remove(self.bowtie2out)
        cmd = "{} {} --nproc 8 -o {} --input_type fastq --sample_id {} --min_cu_len {} --bt2_ps {} --tax_lev {} --min_cu_len {} --stat {} --bowtie2out {}".format(self.metaphlan, self.path, self.out, self.option("sample"), self.option("min_cu_len"), self.option("bt2_ps"), self.option("tax_lev"), self.option("min_cu_len"), self.option("stat"), self.bowtie2out)
        cmd +=" --bowtie2_exe={}".format(self.bowtie2)
        command = self.add_command("run_metaphlan2", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_metaphlan2运行完成！")
        else:
            self.set_error("run_metaphlan2运行完成运行出错!")

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        if len(os.listdir(self.output_dir)) >=1:
            shutil.rmtree(self.output_dir)
        os.link(self.out, self.output_dir + "/" + self.option("sample") + ".profiled.txt")
        os.link(self.bowtie2out, self.output_dir + "/" + self.option("sample") + ".bowtie2out.txt")

    def run(self):
        super(MetaphlanTool, self).run()
        self.run_metaphlan2()
        self.set_output()
        self.end()