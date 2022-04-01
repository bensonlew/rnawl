# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError


class Blast2goAgent(Agent):
    """
    Running Blast2go
    last_modified: 20180927
    """

    def __init__(self, parent):
        super(Blast2goAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "align.blast.blast_xml"},
            {"name": "blast2go_annot", "type": "outfile", "format": "annotation.go.blast2go_annot"}
        ]
        self.add_option(options)
        self.step.add_steps('blast2go')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self.queue = 'BLAST2GO'  # 投递到指定的队列BLAST2GO

    def step_start(self):
        self.step.blast2go.start()
        self.step.update()

    def step_end(self):
        self.step.blast2go.finish()
        self.step.update()

    def check_options(self):
        if self.option("blastout").is_set:
            pass
        else:
            raise OptionError("BLAST result file with xml type must be provide!")

    def set_resource(self):
        self._cpu = 5
        self._memory = '35G'

    def end(self):
        super(Blast2goAgent, self).end()


class Blast2goTool(Tool):

    def __init__(self, config):
        super(Blast2goTool, self).__init__(config)
        self.set_environ(JAVA_HOME=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0")
        self.set_environ(JRE_HOME=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/jre")
        self.set_environ(CLASSPATH=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/lib")
        self.set_environ(CLASSPATH=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/jre/lib")
        self._version = "1.0"
        self.b2g_user = "biocluster102"
        self.b2g_password = "sanger-dev-123"
        self.xml_file = self.option('blastout').prop['path']
        self.pre = self.xml_file.split('/')[-1]

    def run(self):
        super(Blast2goTool, self).run()
        self.run_b2g()
        self.set_output()

    def run_b2g(self):
        cmd = '/program/sun_jdk1.8.0/bin/java -Xmx30g -cp ' + self.config.SOFTWARE_DIR + '/bioinfo/annotation/b2g4pipe_v2.5/*:'
        cmd += self.config.SOFTWARE_DIR + '/bioinfo/annotation/b2g4pipe_v2.5/ext/*: es.blast2go.prog.B2GAnnotPipe'
        cmd += ' -in {} -prop {}/bioinfo/annotation/b2g4pipe_v2.5/b2gPipe.properties -annot -out {}'.format(self.xml_file, self.config.SOFTWARE_DIR, self.work_dir + '/' + self.pre + '.blast2go')
        self.logger.info('运行b2g程序输入为{}'.format(self.xml_file))
        self.logger.info(cmd)
        b2g = self.add_command('blast2go', cmd)
        b2g.run()
        self.wait('blast2go')
        if b2g.return_code == 0:
            self.logger.info('运行blast2go完成')
        else:
            self.set_error('running b2g error')

    def set_output(self):
        if os.path.exists(self.output_dir + '/' + self.pre + '.blast2go.annot'):
            os.remove(self.output_dir + '/' + self.pre + '.blast2go.annot')
        os.link(self.work_dir + '/' + self.pre + '.blast2go.annot',self.output_dir + '/' + self.pre + '.blast2go.annot')
