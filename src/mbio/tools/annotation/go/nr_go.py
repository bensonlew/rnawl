# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os,re
from biocluster.core.exceptions import OptionError
import subprocess


class NrGoAgent(Agent):
    """
    to perform Gene Ontology Annotation
    author: liubinxu
    last_modified: 20181112
    """

    def __init__(self, parent):
        super(NrGoAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "sequence.profile_table"},
            {"name": "blast2go_annot", "type": "outfile", "format": "annotation.go.blast2go_annot"},
            {"name": "p", "type": "int", "default": 10},
        ]
        self.add_option(options)
        self.step.add_steps('go_annotation')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self._memory_increase_step = 60  # 每次重运行增加内存60G by qingchen.zhang @ 20201204

    def step_start(self):
        self.step.go_annotation.start()
        self.step.update()

    def step_end(self):
        self.step.go_annotation.finish()
        self.step.update()

    def check_options(self):
        if self.option("blastout").is_set:
            pass
        else:
            raise OptionError("必须提供BLAST结果文件", code = "31205001")

    def set_resource(self):
        self._cpu = self.option("p") + 1
        self._memory = '120G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["./blast2go.annot", "annot", "Go annotation based on blast output"],
            ["./query_gos.list", "list", "Merged Go annotation"],
        ])
        super(NrGoAgent, self).end()


class NrGoTool(Tool):

    def __init__(self, config):
        super(NrGoTool, self).__init__(config)
        self.set_environ(JAVA_HOME=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0")
        self.set_environ(JRE_HOME=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/jre")
        self.set_environ(CLASSPATH=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/lib")
        self.set_environ(CLASSPATH=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/jre/lib")
        self._version = "1.0"
        self.b2g_user = "biocluster102"
        self.b2g_password = "sanger-dev-123"
        # self.idmapping_db = self.config.SOFTWARE_DIR + "/database/Annotation/other/idmapping.tb"
        self.idmapping_db = self.config.SOFTWARE_DIR + "/database/Annotation/all/PIR/version_20200628/idmapping.tb" # 升级数据库 20201111
        self.nr2go_script = self.config.PACKAGE_DIR + "/prok_rna/get_GO_from_blast_by_nr.py"

    def run(self):
        super(NrGoTool, self).run()
        self.nr2go()
        self.set_out()
        self.end()

    def nr2go(self):
        cmd1 = 'miniconda2/bin/python {}'.format(self.nr2go_script)
        cmd1 += ' %s %s %s %s' % (self.option('blastout').prop['path'],
                                  self.idmapping_db,
                                  self.option("p"),
                                  "go_annot.xls")
        # Config.DB_HOST,Config.DB_USER,Config.DB_PASSWD
        self.logger.info("运行NR2GO")
        self.logger.info(cmd1)
        command = self.add_command("nr2go", cmd1)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行nr2go结束")
        elif command.return_code in [1, '1']:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by qingchen.zhang @20201204
        else:
            self.set_error("运行nr2go出错", code = "31205001")

    def set_out(self):
        outfiles = ["go_annot.xls"]
        for item in outfiles:
            linkfile = self.output_dir + '/' + item
            if os.path.exists(linkfile):
                os.remove(linkfile)
            os.link(self.work_dir + '/' + item, linkfile)
        self.option("blast2go_annot", self.output_dir + "/go_annot.xls")
