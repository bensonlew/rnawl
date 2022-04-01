from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.files.align.bwa.bam import BamFile
import glob
import os

class CashBamAgent(Agent):
    def __init__(self, parent):
        super(CashBamAgent, self).__init__(parent)
        options = [
            {"name": "case_name", "type": "string", "default": None},
            {"name": "Control_name", "type": "string", "default": None},
            {"name": "case_bam", "type": "string", "default": None},
            {"name": "control_bam", "type": "string", "default": None},
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},
            {"name": "outprefix", "type": "string", "default": None}
        ]
        self.add_option(options)
        self.step.add_steps('cash_bam')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def check_options(self):
        if self.option("case_bam") is None:
            raise OptionError("case_bam is not set")
        if self.option("control_bam") is None:
            raise OptionError("control_bam is not set")
        for bam_A in self.option('case_bam').strip().split(","):
            file = BamFile()
            file.set_path(bam_A)
            if not file.check():
                raise Exception("%s not exist" % bam_A)
        for bam_B in self.option('control_bam').strip().split(","):
            file = BamFile()
            file.set_path(bam_B)
            if not file.check():
                raise Exception("%s not exist" % bam_B)

    def set_resource(self):
        self._memory = '20G'
        self._cpu = 8

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["alldiff.statistics.txt","","AS事件统计表"]
            ["alldiff.txt","","AS事件详细表"]
        ])

        super(CashBamAgent,self).end()

class CashBamTool(Tool):
    def __init__(self,config):
        super(CashBamTool.self).__init__(config)
        self._version = 'v2.2.1' #cash last release
        self.java_path = 'program/sun_jdk1.8.0/bin/java'
        self.soft_path = 'bioinfo/rna/cash_v2.2.1/cash.jar'

    def run_cash(self):
        cmd = self.java_path
        cmd += ' -jar {} '.format(self.soft_path)
        cmd += '--Case:{} '.format(self.option("case_name"))
        cmd += '{} '.format(self.option("case_bam"))
        cmd += '--Control:{} '.format(self.option("control_name"))
        cmd += '{} '.format(self.option("control_bam"))
        cmd += '--GTF {} '.format(self.option("ref_gtf"))
        cmd += '--Output {} '.format(self.option("outprefix"))
        self.logger.info('开始运行cash')
        cash_cmd = self.add_command("cash_cmd",cmd).run()
        if cash_cmd.return_code == 0:
            self.logger.info("cash运行完成")
        else:
            self.set_error("cash运行出错!")

    def set_output(self):
        as_detail = glob.glob(self.option("outprefix") + '/*alldiff.txt')
        as_stat = glob.glob(self.option("outprefix") + '/*alldiff.txt')
        all_file = as_detail + as_stat
        for each in all_file:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir,fname)
            if os.path.exist(link):
                os.remove(link)
            os.link(each,link)

    def run(self):
        super(CashBamTool,self).run()
        self.run_cash()
        self.set_output()

