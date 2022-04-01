## !/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "hongdongxuan"
#last_modify:20161125

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import shutil


class ReBam2tabAgent(Agent):
    """
    调用re_bam2tab.sh脚本，完成将合并bam文件，并转换成*.mem.sort.hit.vcf.tab文件
    version v1.0
    author: hongdongxuan
    last_modify: 2016.11.25
    """
    def __init__(self, parent):
        super(ReBam2tabAgent, self).__init__(parent)
        options = [
            {"name": "bam_file1", "type": "string"},  # bam文件1
            {"name": "bam_file2", "type": "string"},  # bam文件2
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},  # 参考序列
            {"name": "targets_bedfile", "type": "string"} #位点信息
        ]
        self.add_option(options)
        self.step.add_steps("re_Bam2tab")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.re_Bam2tab.start()
        self.step.update()

    def stepfinish(self):
        self.step.re_Bam2tab.finish()
        self.step.update()


    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("bam_file1"):
            raise OptionError("必须输入bam1文件的所在路径")
        if not self.option("bam_file2"):
            raise OptionError("必须输入bam2文件的所在路径")
        if not self.option("ref_fasta").is_set:
            raise OptionError("必须输入参考基因组序列fasta文件")
        if not self.option('targets_bedfile'):
            raise OptionError('必须提供target_bedfile文件')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        self._memory = '200G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
                    ])
        result_dir.add_regexp_rules([
            [r".mem.sort.hit.vcf.tab", "tab", "所有位点的信息"],
        ])
        super(ReBam2tabAgent, self).end()


class ReBam2tabTool(Tool):
    """
    运行脚本：re_bam2tab.sh sample_id bam_file1 bam_file2 ref targets_bedfile
    """
    def __init__(self, config):
        super(ReBam2tabTool, self).__init__(config)
        self._version = '1.0.1'
        self.cmd_path = "bioinfo/medical/scripts/re_bam2tab.sh"

    def run_reBam2tab(self):
        reBam2tab_cmd = self.cmd_path + " 0 %s %s %s %s" % (self.option("bam_file1"), self.option("bam_file2"),
                                                     self.option("ref_fasta").prop["path"], self.option("targets_bedfile"))
        print reBam2tab_cmd
        self.logger.info(reBam2tab_cmd)
        self.logger.info("开始运行cmd")
        cmd = self.add_command("cmd", reBam2tab_cmd).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("运行Bam2tab成功")
        else:
            self.logger.info("运行Bam2tab出错")

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        results = os.listdir(self.work_dir + "/")
        for f in results:
            if re.search(r'.*mem\.sort\.hit\.vcf\.tab$', f):
                shutil.copy(self.work_dir + "/" + f, self.output_dir)
            else:
                pass
        self.logger.info('设置文件夹路径成功')

    def run(self):
        super(ReBam2tabTool, self).run()
        self.run_reBam2tab()
        self.set_output()
        self.end()