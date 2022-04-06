# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190321

import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class SsrMergeAgent(Agent):
    """
    合并参考基因组和各样本的SSR结果，组成新的vcf文件
    """
    def __init__(self, parent):
        super(SsrMergeAgent, self).__init__(parent)
        options = [
            # {"name": "ref_misa", "type": "infile", "format": "wgs_v2.bcf", "required": True},  # 参考基因组的SSR:ref.fa.newmisa
            {"name": "ref_misa", "type": "string", "required": True},  # 参考基因组的SSR:ref.fa.newmisa
            {"name": "ssr_list", "type": "infile", "format": "wgs_v2.bam_list", "required": True},  # 样本对应的ssr.result结果文件list
        ]
        self.add_option(options)

    def check_options(self):
        # if not self.option("ref_misa").is_set:
        #     raise OptionError("请设置参考基因组ref.fa.newmisa")
        if not self.option("ssr_list").is_set:
            raise OptionError("请设置样本的ssr.result结果文件list")

    def set_resource(self):
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(SsrMergeAgent, self).end()


class SsrMergeTool(Tool):
    def __init__(self, config):
        super(SsrMergeTool, self).__init__(config)
        self.python_path = "miniconda2/bin/python"
        self.ssr_merge_path = self.config.PACKAGE_DIR + "/wgs_v2/ssr_merge.py"
        self.ref_misa_path = self.config.SOFTWARE_DIR + "/database/dna_geneome/" + self.option("ref_misa")
        if not os.path.exists(self.ref_misa_path):
            raise OptionError("请设置参考基因组ref.fa.newmisa")

    def run_ssr_merge(self):
        """
        合并参考基因组和各样本的SSR结果，组成新的vcf文件
        """
        # cmd = "{} {} -r {} ".format(self.python_path, self.ssr_merge_path, self.option("ref_misa").prop["path"])
        cmd = "{} {} -r {} ".format(self.python_path, self.ssr_merge_path, self.ref_misa_path)
        cmd += "-l {} -o {}".format(self.option("ssr_list").prop["path"], os.path.join(self.output_dir, "final.ssr.result.vcf"))
        command = self.add_command("ssr_merge", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("SSR合并成功")
        else:
            self.set_error("SSR合并失败，请检查")

    def run(self):
        super(SsrMergeTool, self).run()
        self.run_ssr_merge()
        self.end()
