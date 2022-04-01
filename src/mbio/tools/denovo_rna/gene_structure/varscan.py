# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from mbio.packages.denovo_rna.gene_structure.snp_position import snp_stat


class VarscanAgent(Agent):
    """
    varscan:SNP calling软件
    version 1.0
    author: qindanhua
    last_modify: 2016.07.11
    """

    Method = ["pileup2snp", "mpileup2snp", "pileup2indel", "mpileup2indel"]

    def __init__(self, parent):
        super(VarscanAgent, self).__init__(parent)
        options = [
            {"name": "pileup", "type": "infile", "format": "denovo_rna.gene_structure.pileup"},  # mpileup 输出格式
            {"name": "method", "type": "string", "default": "pileup2snp"},  # mpileup 输出格式
            {"name": "bed", "type": "infile", "format": "denovo_rna.gene_structure.bed"},  # bed格式文件
            # {"name": "vcf", "type": "outfile", "format": "vcf"}     # Variant Call Format
        ]
        self.add_option(options)
        self.step.add_steps('varscan')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.varscan.start()
        self.step.update()

    def step_end(self):
        self.step.varscan.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        """
        if not self.option("pileup").is_set:
            raise OptionError("请传入pileup文件")
        if self.option("method") not in self.Method:
            raise OptionError("选择正确的工具")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 10
        self._memory = '11G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            # ["./estimators.xls", "xls", "alpha多样性指数表"]
        ])
        if self.option("bed").is_set:
            result_dir.add_relpath_rules([
                ["./snp_position_stat.xls", "xls", "snp编码位置信息统计表"],
                ["./snp_type_stat.xls", "xls", "snp类型统计表"],
                ["./snp.xls", "xls", "snp信息表"]
            ])
        else:
            result_dir.add_relpath_rules([
                ["./pileup_out.xls", "xls", "snp信息表"]
            ])
        # print self.get_upload_files()
        super(VarscanAgent, self).end()


class VarscanTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(VarscanTool, self).__init__(config)
        self.varscan_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/VarScan.v2.3.9.jar"
        self.java_path = "program/sun_jdk1.8.0/bin/"

    def pileup2snp(self):
        cmd = "{}java -jar {} pileup2snp {} --min-coverage 8 --min-reads2 3 --min-strands2 2 --min-avg-qual 30 " \
              "--min-var-freq 0.30".format(self.java_path, self.varscan_path, self.option("pileup").prop["path"])
        self.logger.info("开始运行pileup2snp")
        command = self.add_command("pileup2snp", cmd)
        self.logger.info(cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行pileup2snp结束")
        else:
            self.set_error("运行pileup2snp出错")

    def set_output(self):
        self.logger.info("set out put")
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        if self.option("bed").is_set:
            snp_stat(self.work_dir+'/pileup2snp.o', self.option("bed").prop["path"])
            os.link(self.work_dir + "/snp.xls", self.output_dir + "/snp.xls")
            os.link(self.work_dir + "/snp.type.stat.xls", self.output_dir + "/snp.type.stat.xls")
            os.link(self.work_dir + "/snp.position.stat.xls", self.output_dir + "/snp.position.stat.xls")
        else:
            os.link(self.work_dir+'/pileup2snp.o', self.output_dir+'/pileup_out.xls')
        self.logger.info("done")

    def run(self):
        """
        运行
        """
        super(VarscanTool, self).run()
        if self.option("method") == "pileup2snp":
            self.pileup2snp()
        self.set_output()
        self.end()
