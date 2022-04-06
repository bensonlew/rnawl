# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess
from mbio.packages.rna.annot_config import AnnotConfig


class Nr2goAgent(Agent):
    """
    to perform Gene Ontology Annotation
    author: liubinxu
    last_modified: 20180509
    """

    def __init__(self, parent):
        super(Nr2goAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "blastout_go", "type": "outfile", "format": "ref_rna_v2.blast_xml"},
            {"name": "blast2go_annot", "type": "outfile", "format": "ref_rna_v2.blast2go_annot"},
            {"name": "pir_version", "type": "string", "default": "2019"},
            {"name": "known_go", "type": "string", "default": None},
            {"name": "p", "type": "int", "default": 10},
        ]
        self.add_option(options)
        self.step.add_steps('go_annotation')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        # self.queue = 'BLAST2GO'  # 投递到指定的队列BLAST2GO

    def step_start(self):
        self.step.go_annotation.start()
        self.step.update()

    def step_end(self):
        self.step.go_annotation.finish()
        self.step.update()

    def check_options(self):
        if self.option("blastout").is_set:
            '''
            document = ET.__parse_details(self.option("blastout").prop['path'])
            root = document.getroot()
            db = root.find('BlastOutput_db')
            if db.text == 'nr':
                pass
            else:
                raise OptionError("BLAST比对数据库不支持")
            '''
            pass
        else:
            raise OptionError("必须提供BLAST结果文件", code = "33710902")

    def set_resource(self):
        self._cpu = self.option("p") + 1
        self._memory = '75G'

    def end(self):
        super(Nr2goAgent, self).end()


class Nr2goTool(Tool):

    def __init__(self, config):
        super(Nr2goTool, self).__init__(config)
        self.set_environ(JAVA_HOME=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0")
        self.set_environ(JRE_HOME=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/jre")
        self.set_environ(CLASSPATH=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/lib")
        self.set_environ(CLASSPATH=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/jre/lib")
        self._version = "1.0"
        self.b2g_user = "biocluster102"
        self.b2g_password = "sanger-dev-123"
        self.idmapping_db = self.idmapping_db = AnnotConfig().get_file_path(
            file ="idmapping.tb",
            db = "pir",
            version = self.option("pir_version"))
        self.nr2go_script = self.config.PACKAGE_DIR + "/prok_rna/get_GO_from_blast_by_nr2.py"


    def run(self):
        super(Nr2goTool, self).run()
        self.convert_xml2table()
        self.nr2go()
        self.set_out()
        self.end()

    def convert_xml2table(self):
        self.logger.info('转换xml结果为表格格式')
        self.option("blastout").convert2table("blast_table.xls")

    def nr2go(self):
        cmd1 = 'miniconda2/bin/python {}'.format(self.nr2go_script)
        cmd1 += ' %s %s %s %s' % ("blast_table.xls",
                                  self.idmapping_db,
                                  self.option("p"),
                                  "blast2go_annot.xls")
        # Config.DB_HOST,Config.DB_USER,Config.DB_PASSWD
        self.logger.info("运行NR2GO")
        self.logger.info(cmd1)
        command = self.add_command("nr2go", cmd1)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行nr2go结束")
        else:
            self.set_error("运行nr2go出错", code = "33710902")

    def set_out(self):
        outfiles = ["blast2go_annot.xls"]
        for item in outfiles:
            linkfile = self.output_dir + '/' + item
            if os.path.exists(linkfile):
                os.remove(linkfile)
            os.link(self.work_dir + '/' + item, linkfile)
        self.option("blast2go_annot", self.output_dir + "/blast2go_annot.xls")
