# -*- coding: utf-8 -*-
# __author__ = 'liuwentian'
# modified 2018.06.14

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class SplitbylgAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(SplitbylgAgent, self).__init__(parent)
        options = [
            {"name": "bin", "type": "string"},  # 判断是否为bin，传入参数为yes/no
            {"name": "poptype", "type": "string"},  # 传入poptype参数，例：CP,F2等
            {"name": "marker_path", "type": "infile", "format": "dna_gmap.marker"},
            # 传入marker文件，没bin：pop.filtered.marker；有bin：Total.bin.marker
            {"name": "lg_path", "type": "infile", "format": "dna_gmap.lg"},  # 传入Total.lg
            {"name": "file_path", "type": "outfile", "format": "dna_gmap.marker_path"}  # 将output的路径作为tool：smooth的输入文件。
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("bin"):
            raise OptionError("必须输入bin", code="34801301")
        if not self.option("poptype"):
            raise OptionError("必须输入poptype", code="34801302")
        if not self.option('marker_path').is_set:
            raise OptionError("必须输入marker文件！", code="34801303")
        if not self.option("lg_path"):
            raise OptionError("必须输入Total.lg", code="34801304")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        super(SplitbylgAgent, self).end()


class SplitbylgTool(Tool):
    def __init__(self, config):
        super(SplitbylgTool, self).__init__(config)
        self.perl_path = 'miniconda2/bin/perl'
        self.splitbyLG_CP_pl = self.config.PACKAGE_DIR + '/dna_gmap/splitbyLG-CP.pl'
        self.splitbyLG_NOCP_pl = self.config.PACKAGE_DIR + '/dna_gmap/splitbyLG-NOCP.pl'

    def clean_output(self):
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))

    def run_splitbylg_cp(self):
        """
        poptype:CP
        splitbyLG_CP
        """
        cmd = "{} {} -l {} -i {} -d {} -t {}".format(self.perl_path, self.splitbyLG_CP_pl,
                                                     self.option('lg_path').prop["path"],
                                                     self.option('marker_path').prop["path"],
                                                     self.output_dir, self.option('poptype'))
        command = self.add_command("splitbylg_cp", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("splitbyLG_CP完成")
        else:
            self.set_error("splitbyLG_CP失败", code="34801301")
        if os.path.exists(self.output_dir):
            self.option("file_path").set_path(self.output_dir)

    def run_splitbylg_nocp(self):
        """
        poptype:NOCP
        splitbyLG_NOCP
        """
        cmd = "{} {} -l {} -i {} -d {} -t {}".format(self.perl_path, self.splitbyLG_NOCP_pl,
                                                     self.option('lg_path').prop["path"],
                                                     self.option('marker_path').prop["path"],
                                                     self.output_dir, self.option('poptype'))
        command = self.add_command("splitbylg_nocp", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("splitbyLG_NOCP完成")
        else:
            self.set_error("splitbyLG_NOCP失败", code="34801302")
        if os.path.exists(self.output_dir):
            self.option("file_path").set_path(self.output_dir)

    def run(self):
        super(SplitbylgTool, self).run()
        self.clean_output()
        if self.option('poptype') == 'CP':
            self.run_splitbylg_cp()
        else:
            self.run_splitbylg_nocp()
        self.end()
