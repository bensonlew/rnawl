# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.10

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class CnvAnnoAgent(Agent):
    """
    cnv 注释
    """
    def __init__(self, parent):
        super(CnvAnnoAgent, self).__init__(parent)
        options = [
            {"name": "cnv_file", "type": "string"},  # call cnv文件
            {"name": "ref_gff", "type": "string"},  # ref.gff
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("cnv_file"):
            raise OptionError("请设置cnv文件", code="34501401")
        if not self.option("ref_gff"):
            raise OptionError("请设置ref.gff文件", code="34501402")
        if not os.path.exists(self.option("cnv_file")):
            raise OptionError("cnv文件:%s不存在，请检查",variables=(self.option("cnv_file")), code="34501403")
        if not os.path.exists(self.option("ref_gff")):
            raise OptionError("ref_gff文件:%s不存在，请检查",variables=(self.option("ref_gff")), code="34501404")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(CnvAnnoAgent, self).end()


class CnvAnnoTool(Tool):
    def __init__(self, config):
        super(CnvAnnoTool, self).__init__(config)
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.cnv_anno = self.config.PACKAGE_DIR + "/wgs/cnv_anno.pl"

    def run_cnv_anno(self):
        """
        cnv_anno.pl cnv注释
        """
        anno_file = os.path.basename(self.option("cnv_file")).split(".cnv")[0] + ".cnv.anno.xls"
        cmd = "{} {} -i {} ".format(self.perl_path, self.cnv_anno, self.option("cnv_file"))
        cmd += "-g {} -o {}".format(self.option("ref_gff"), self.work_dir + "/" + anno_file)
        command = self.add_command("cnv_anno", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("cnv_anno.pl运行完成")
        else:
            self.set_error("cnv_anno.pl运行失败", code="34501401")
        if os.path.exists(self.output_dir + "/" + anno_file):
            os.remove(self.output_dir + "/" + anno_file)
        os.link(self.work_dir + "/" + anno_file, self.output_dir + "/" + anno_file)

    def run(self):
        super(CnvAnnoTool, self).run()
        self.run_cnv_anno()
        self.end()
