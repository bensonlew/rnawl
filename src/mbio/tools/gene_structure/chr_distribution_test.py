# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import glob
import subprocess


class ChrDistributionTestAgent(Agent):
    """
    Rseqc-2.3.6:RNA测序分析工具
    version 1.0
    author: qindanhua
    last_modify: 2016.07.27
    """

    def __init__(self, parent):
        super(ChrDistributionTestAgent, self).__init__(parent)
        options = [
            {"name": "bam", "type": "infile", "format": "align.bwa.bam"},  # bam格式文件,排序过的
            {"name": "freq", "type": "int", "default": 100000},  # 统计窗口
            {"name": "chr_info", "type": "string", "default": ''},  # 统计窗口

        ]
        self.add_option(options)
        self.step.add_steps('dup')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.dup.start()
        self.step.update()

    def step_end(self):
        self.step.dup.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        """
        if not self.option("bam").is_set:
            raise OptionError("请传入比对结果bam格式文件")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 10
        self._memory = '15G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        super(ChrDistributionTestAgent, self).end()


class ChrDistributionTestTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(ChrDistributionTestTool, self).__init__(config)
        self.python_path = "program/Python/bin/"
        self.python_full_path = self.config.SOFTWARE_DIR + "/program/Python/bin/"
        self.samtools_path = self.config.SOFTWARE_DIR + "/bioinfo/align/samtools-1.3.1/"
        self.circos_script = self.config.SOFTWARE_DIR + '/bioinfo/plot/scripts/circos_distribution.py'
        self.bam_name = ""

    def bam_to_sam(self):
        bam_path = self.option("bam").prop["path"]
        self.bam_name = bam_path.split("/")[-1]
        self.logger.info(self.work_dir + "/" + self.bam_name + ".sam")
        cmd = "{}samtools view {} > {}".format(self.samtools_path, bam_path, self.work_dir + "/" + self.bam_name + ".sam")
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError:
            self.logger.info("bam_to_sam出错")
            return False

    def chr_distribution(self):
        # bam_path = self.option("bam").prop["path"]
        # self.bam_name = bam_path.split("/")[-1]
        cmd = "cat {} | cut -f 3 |uniq -c > {}".format(self.bam_name + ".sam", self.bam_name + "_chr_stat.xls")
        self.logger.info(cmd)
        self.logger.info("开始运行统计reads在染色体的分布情况")
        try:
            result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
            with open("chr_stat.out", "w") as w:
                w.write(result)
            self.logger.info("开始运行统计reads在染色体的分布情况结束")
            return True
        except subprocess.CalledProcessError:
            self.logger.info("cat统计reads在染色体的分布出错")
            return False

    def circos_stat(self):
        # circos_stst = self.work_dir + "/" + self.bam_name + "/circos_stat"
        circos_stst = self.work_dir + "/circos_stat"
        sam_path = self.work_dir + "/" + self.bam_name + ".sam"
        self.logger.info(self.option("chr_info"))
        if not os.path.exists(circos_stst):
            os.mkdir(circos_stst)
        cmd = "{}python {} -i {} -l {} -f {} -o {}".format(self.python_path, self.circos_script, sam_path,
                                                           self.option("chr_info"), self.option("freq"), circos_stst)
        command = self.add_command("circos_stat", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行circos_stat结束")
        else:
            self.set_error("运行circos_stat出错")

    def set_output(self):
        self.logger.info("set out put")
        out_put = self.output_dir + "/" + self.bam_name + "_chr_stat.xls"
        self.logger.info(out_put)
        # self.logger.info("sed \'$d\' {}".format(self.bam_name + "_chr_stat.xls"))
        # self.logger.info("awk \'BEGIN{print \"#chr\tread_num\"};{print $2\" \"$1}\'")
        cmd = "sed \'$d\' %s |awk \'BEGIN{print \"#chr\tread_num\"};{print $2\" \"$1}\' | awk \'NR==1||$2>100000||length($1)<8\' |sort -n -k 1 > %s" % (self.bam_name + "_chr_stat.xls", out_put)
        self.logger.info(cmd)
        os.system(cmd)
        self.logger.info("set done")
        self.end()

    def run(self):
        """
        运行
        """
        super(ChrDistributionTestTool, self).run()
        self.bam_to_sam()
        self.chr_distribution()
        if not self.option("chr_info") == "":
            self.circos_stat()
        self.set_output()
