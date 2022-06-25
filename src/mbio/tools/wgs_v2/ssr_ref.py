# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190320

import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class SsrRefAgent(Agent):
    """
    参考基因组设计SSR,需将结果ref.fa.misa的路径写入mongo的sg_task表，用于样本的SSR设计
    """
    def __init__(self, parent):
        super(SsrRefAgent, self).__init__(parent)
        options = [
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta", "required": True},  # 参考基因组的ref.fa文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("ref_fa").is_set:
            raise OptionError("请设置参考基因组ref.fa")

    def set_resource(self):
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(SsrRefAgent, self).end()


class SsrRefTool(Tool):
    def __init__(self, config):
        super(SsrRefTool, self).__init__(config)
        self.perl_path = "miniconda2/bin/perl"
        self.python_path = "miniconda2/bin/python"
        self.misa_path = self.config.PACKAGE_DIR + "/wgs_v2/misa.pl"
        self.misa_ini = self.config.PACKAGE_DIR + "/wgs_v2/misa.ini"
        self.ssr_stat_path = self.config.PACKAGE_DIR + "/wgs_v2/ssr_misa_stat.py"

    def run_misa(self):
        """
        运行misa.pl,对ref.fa进行SSR设计,输出结果ref.fa.ssrfa、ref.fa.newmisa、ref.fa.misa、ref.fa.statistics
        """
        ref_fa = os.path.join(self.work_dir, "ref.fa")
        if os.path.exists(ref_fa):
            os.remove(ref_fa)
        os.link(self.option("ref_fa").prop["path"], ref_fa)
        cmd = "{} {} {} {}".format(self.perl_path, self.misa_path, ref_fa, self.misa_ini)
        command = self.add_command("ref_misa", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("对ref.fa设计SSR成功")
        else:
            self.set_error("对ref.fa设计SSR失败，请检查")

    def run_ssr_stat(self):
        """
        对misa的结果ref.fa.misa进行统计，得到SSR类型分布进行统计
        """
        if os.path.exists(os.path.join(self.output_dir, "ref.fa.misa")):
            os.remove(os.path.join(self.output_dir, "ref.fa.misa"))
        os.link(os.path.join(self.work_dir, "ref.fa.misa"), os.path.join(self.output_dir, "ref.fa.misa"))
        cmd = "{} {} -i {} ".format(self.python_path, self.ssr_stat_path, os.path.join(self.work_dir, "ref.fa.misa"))
        cmd += "-o {} -s {}".format(os.path.join(self.output_dir, "ssr.stat.xls"), "Reference")
        command = self.add_command("ssr_stat", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("SSR统计成功")
        else:
            self.set_error("SSR统计失败，请检查")

    def run(self):
        super(SsrRefTool, self).run()
        self.run_misa()
        self.run_ssr_stat()
        self.end()
