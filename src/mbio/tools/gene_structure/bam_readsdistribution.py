# -*- coding: utf-8 -*-
# __author__ = 'zengjing'

"""reads 区域分布"""

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import subprocess
import glob
import shutil

class BamReadsdistributionAgent(Agent):
    """
    Rseqc-2.3.6:RNA测序分析工具
    reads区域分布
    last_modify: 2016.09.09
    """
    version = ""

    def __init__(self, parent):
        super(BamReadsdistributionAgent, self).__init__(parent)
        options = [
            {'name': 'bam', 'type': 'infile','format':'align.bwa.bam,align.bwa.bam_dir'},  # 排序过的bam格式文件
            {'name': 'bed', 'type': 'infile','format':'gene_structure.bed'}    # bed格式文件
            #{'name': 'output', 'type': 'outfile','format':'ref_rna.txt'}
        ]
        self.add_option(options)
        self._memory_increase_step = 50
        self.step.add_steps('reads_distribution')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.reads_distribution.start()
        self.step.update()

    def step_end(self):
        self.step.reads_distribution.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        """
        if not self.option('bam').is_set:
            raise OptionError("请传入比对结果bam格式文件", code = "32200501")
        if not self.option('bed').is_set:
            raise OptionError("请传入参考基因组结构注释bed格式文件", code = "32200502")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        self._memory = '20G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        super(BamReadsdistributionAgent, self).end()


class BamReadsdistributionTool(Tool):
    def __init__(self, config):
        super(BamReadsdistributionTool, self).__init__(config)
        self.cmd_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/"
        self.shell_path = "bioinfo/rna/scripts"
        python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/'
        self.set_environ(PATH=python_path)

    def run_reads_distribution(self):
        """
        运行read_distribution
        """
        cmd = "{}/read_distribution.sh {} {} {} {}"\
            .format(self.shell_path, self.cmd_path, self.option("bam").prop["path"], self.option("bed").prop["path"],
                    os.path.basename(self.option("bam").prop["path"]).split(".")[0])
        self.logger.info("开始运行read_distribution")
        self.logger.info(cmd)
        # try:
        #     subprocess.check_output(cmd, shell=True)
        # except subprocess.CalledProcessError:
        #     raise Exception("运行出错!")
        cmd_obj = self.add_command("cmd", cmd, ignore_error=True)
        cmd_obj.run()
#        distribution = os.path.join(self.work_dir, "reads_distribution.txt")
 #       shutil.copyfile(distribution, "../ReadsDistribution/output/reads_distribution.txt")
        self.wait(cmd_obj)
        if cmd_obj.return_code == 0:
            self.logger.info("read_distribution运行完成")
        elif cmd_obj.return_code in [1, -9, 137, 139]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            cmd_obj.rerun()
            self.wait(cmd_obj)
            if cmd_obj.return_code == 0:
                self.logger.info("read_distribution重新运行完成")
            else:
                self.set_error("read_distribution重新运行失败")

    def set_out(self):
        self.logger.info("set out put")
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        file_path = glob.glob(r"*reads_distribution.txt")
        print(file_path)
        for f in file_path:
            output_dir = os.path.join(self.output_dir, f)
            os.link(os.path.join(self.work_dir, f), output_dir)
        self.wait()

    def run(self):
        """
        运行
        """
        super(BamReadsdistributionTool, self).run()
        self.run_reads_distribution()
        self.set_out()
        self.end()
