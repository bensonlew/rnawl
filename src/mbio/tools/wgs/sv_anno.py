# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.10

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class SvAnnoAgent(Agent):
    """
    sv 注释
    """
    def __init__(self, parent):
        super(SvAnnoAgent, self).__init__(parent)
        options = [
            {"name": "sv_file", "type": "string"},  # call sv文件
            {"name": "ref_gff", "type": "string"},  # ref.gff
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sv_file"):
            raise OptionError("请设置sv文件", code="34506801")
        if not self.option("ref_gff"):
            raise OptionError("请设置ref.gff文件", code="34506802")
        if not os.path.exists(self.option("sv_file")):
            raise OptionError("sv文件:%s不存在，请检查",variables=(self.option("sv_file")), code="34506803")
        if not os.path.exists(self.option("ref_gff")):
            raise OptionError("ref_gff文件:%s不存在，请检查",variables=(self.option("ref_gff")), code="34506804")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(SvAnnoAgent, self).end()


class SvAnnoTool(Tool):
    def __init__(self, config):
        super(SvAnnoTool, self).__init__(config)
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.sv_anno = self.config.PACKAGE_DIR + "/wgs/sv_anno.pl"

    def run_sv_anno(self):
        """
        sv_anno.pl sv注释
        """
        anno_file = os.path.basename(self.option("sv_file")).split(".sv")[0] + ".sv.anno.xls"
        cmd = "{} {} -i {} ".format(self.perl_path, self.sv_anno, self.option("sv_file"))
        cmd += "-g {} -o {}".format(self.option("ref_gff"), self.work_dir + "/" + anno_file)
        command = self.add_command("sv_anno", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("sv_anno.pl运行完成")
        else:
            self.set_error("sv_anno.pl运行失败", code="34506801")
        if os.path.exists(self.output_dir + "/" + anno_file):
            os.remove(self.output_dir + "/" + anno_file)
        os.link(self.work_dir + "/" + anno_file, self.output_dir + "/" + anno_file)

    def run(self):
        super(SvAnnoTool, self).run()
        self.run_sv_anno()
        self.end()
