# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190320

import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class BamSsrAgent(Agent):
    """
    通过diff_ssr.pl对每个bam进行SSR设计，得到每个bam的ssr.result
    """
    def __init__(self, parent):
        super(BamSsrAgent, self).__init__(parent)
        options = [
            {"name": "bam_file", "type": "infile", "format": "align.bwa.sam", "required": True},  # 样本对应的bam文件
            # {"name": "ref_misa", "type": "infile", "format": "wgs_v2.bcf", "required": True},  # 参考基因组的SSR:ref.fa.misa
            {"name": "ref_misa", "type": "string", "required": True},  # 参考基因组的SSR:ref.fa.misa
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("bam_file").is_set:
            raise OptionError("请设置样本的bam文件")
        # if not self.option("ref_misa").is_set:
        #     raise OptionError("请设置参考基因组的ref.fa.misa")

    def set_resource(self):
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(BamSsrAgent, self).end()


class BamSsrTool(Tool):
    def __init__(self, config):
        super(BamSsrTool, self).__init__(config)
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.diff_ssr = self.config.PACKAGE_DIR + "/wgs_v2/diff_ssr.pl"
        self.ref_misa_path = self.config.SOFTWARE_DIR + "/database/dna_geneome/" + self.option("ref_misa")
        if not os.path.exists(self.ref_misa_path):
            raise OptionError("请设置参考基因组的ref.fa.misa")

    def run_diff_ssr(self):
        """
        用diff_ssr.pl对bam进行SSR设计
        """
        outfile = os.path.join(self.output_dir, os.path.basename(self.option("bam_file").prop["path"]).split(".bam")[0])
        # cmd = "{} {} -input {} ".format(self.perl_path, self.diff_ssr, self.option("ref_misa").prop["path"])
        cmd = "{} {} -input {} ".format(self.perl_path, self.diff_ssr, self.ref_misa_path)
        cmd += "-bam {} -output {}".format(self.option("bam_file").prop["path"], outfile)
        command = self.add_command("diff_ssr", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("bam设计SSR运行完成")
        else:
            self.set_error("bam设计SSR运行失败")

    def run(self):
        super(BamSsrTool, self).run()
        self.run_diff_ssr()
        self.end()
