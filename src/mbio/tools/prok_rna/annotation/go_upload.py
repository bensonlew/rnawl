# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2017.04.13
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess


class GoUploadAgent(Agent):
    """
    用GO
    """
    def __init__(self, parent):
        super(GoUploadAgent, self).__init__(parent)
        options = [
            {"name": "gos_list_upload", "type": "infile", "format": "prok_rna.anno_upload"},
            {"name": "go2level_out", "type": "outfile", "format": "prok_rna.level2"},
            {"name": "golist_out", "type": "outfile", "format": "prok_rna.go_list"}
        ]
        self.add_option(options)
        self.step.add_steps('go_annotation')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self.queue = 'BLAST2GO'  # 投递到指定的队列BLAST2GO

    def step_start(self):
        self.step.go_annotation.start()
        self.step.update()

    def step_end(self):
        self.step.go_annotation.finish()
        self.step.update()

    def check_options(self):
        if self.option("gos_list_upload").is_set:
            pass
        else:
            raise OptionError("必须提供BLAST结果文件", code = "35000801")

    def set_resource(self):
        self._cpu = 10
        self._memory = '25G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["./blast2go.annot", "annot", "Go annotation based on blast output"],
            ["./query_gos.list", "list", "Merged Go annotation"],
            ["./go1234level_statistics.xls", "xls", "Go annotation on 4 levels"],
            ["./go123level_statistics.xls", "xls", "Go annotation on 3 levels"],
            ["./go12level_statistics.xls", "xls", "Go annotation on 2 levels"],
            # ["./go2level.xls", "xls", "Go annotation on level 2"],
            # ["./go3level.xls", "xls", "Go annotation on level 3"],
            # ["./go4level.xls", "xls", "Go annotation on level 4"]
        ])
        super(GoUploadAgent, self).end()


class GoUploadTool(Tool):

    def __init__(self, config):
        super(GoUploadTool, self).__init__(config)
        self._version = "1.0"
        self.b2g_user = "biocluster102"
        self.b2g_password = "sanger-dev-123"
        self.python = self.config.SOFTWARE_DIR + "/program/Python/bin/python"
        self.goAnnot = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/goAnnot.py"
        self.goSplit = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/goSplit.py"

    def run(self):
        super(GoUploadTool, self).run()
        self.run_annotation()

    def run_annotation(self):
        self.option("gos_list_upload").get_transcript_anno(outdir=self.work_dir + "/query_gos.list")
        self.option("golist_out", self.work_dir + "/query_gos.list")
        cmd = '%s %s %s %s %s %s' % (self.python, self.goAnnot, self.work_dir + "/query_gos.list", 'localhost', self.b2g_user, self.b2g_password)  # 10.100.203.193
        if os.path.exists(self.output_dir + '/query_gos.list'):
            os.remove(self.output_dir + '/query_gos.list')
        os.link(self.work_dir + '/query_gos.list', self.output_dir + '/query_gos.list')
        self.logger.info("运行goAnnot.py")
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info("运行goAnnot.py完成")
            if os.path.exists(self.output_dir + '/go1234level_statistics.xls'):
                os.remove(self.output_dir + '/go1234level_statistics.xls')
            os.link(self.work_dir + '/go1234level_statistics.xls',
                    self.output_dir + '/go1234level_statistics.xls')
            if os.path.exists(self.output_dir + '/go123level_statistics.xls'):
                os.remove(self.output_dir + '/go123level_statistics.xls')
            os.link(self.work_dir + '/go123level_statistics.xls',
                    self.output_dir + '/go123level_statistics.xls')
            if os.path.exists(self.output_dir + '/go12level_statistics.xls'):
                os.remove(self.output_dir + '/go12level_statistics.xls')
            os.link(self.work_dir + '/go12level_statistics.xls',
                    self.output_dir + '/go12level_statistics.xls')
        except subprocess.CalledProcessError:
            self.set_error("运行goAnnot.py出错", code = "35000802")
        # self.run_gosplit()
        self.end()

    def run_gosplit(self):
        cmd = '{} {}'.format(self.python, self.goSplit)
        cmd += ' %s' % self.work_dir + '/go_detail.xls'
        self.logger.info("运行goSplit.py")
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info("运行goSplit.py完成")
            outfiles = ['go2level.xls', 'go3level.xls', 'go4level.xls']
            for item in outfiles:
                linkfile = self.output_dir + '/' + item
                if os.path.exists(linkfile):
                    os.remove(linkfile)
                os.link(self.work_dir + '/' + item, linkfile)
            self.option('go2level_out', self.output_dir + '/go2level.xls')
        except subprocess.CalledProcessError:
            self.set_error("运行goSplit.py出错")
        self.end()
