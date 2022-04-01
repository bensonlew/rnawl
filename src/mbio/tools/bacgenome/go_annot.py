# -*- coding: utf-8 -*-
# __author__ = 'qignchen.zhang'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess


class GoAnnotAgent(Agent):
    """
    """

    def __init__(self, parent):
        super(GoAnnotAgent, self).__init__(parent)
        options = [
            {"name": "go2level_infile", "type": "infile", "format": "annotation.go.blast2go_annot"},
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
        if self.option("go2level_infile").is_set:
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
            raise OptionError("必须提供BLAST结果文件")

    def set_resource(self):
        self._cpu = 10
        self._memory = '60G'

    def end(self):
        super(GoAnnotAgent, self).end()


class GoAnnotTool(Tool):

    def __init__(self, config):
        super(GoAnnotTool, self).__init__(config)
        self.set_environ(JAVA_HOME=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0")
        self.set_environ(JRE_HOME=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/jre")
        self.set_environ(CLASSPATH=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/lib")
        self.set_environ(CLASSPATH=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/jre/lib")
        self._version = "1.0"
        self.b2g_user = "biocluster102"
        self.b2g_password = "sanger-dev-123"

    def run(self):
        super(GoAnnotTool, self).run()
        self.format_go()

    def format_go(self):
        """
        根据输入文件整理一下
        :return:
        """
        go_file = self.work_dir  + "/go_statistics2.xls"
        with open(self.option("go2level_infile").prop['path'], 'r') as f, open(go_file, 'w') as w:
            for line in f:
                spline=line.strip().split("\t")
                gene = spline[0]
                go_name = spline[1].split(";")
                go_des = spline[3]
                for go in go_name:
                    w.write("{}\t{}\t{}\n".format(gene, go, go_des))
        self.run_not_level()

    def run_not_level(self):
        cmd = '{}/program/Python/bin/python {}/annotation/go/go_desc.py'.format(self.config.SOFTWARE_DIR, self.config.PACKAGE_DIR)
        cmd += ' {} {}'.format(self.work_dir  + "/go_statistics2.xls", self.work_dir + '/go_statistics.xls' )
        self.logger.info('运行go_desc.py')
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd,shell=True)
            self.logger.info('运行go_desc.py 完成')
            if os.path.exists(self.output_dir + '/go_statistics.xls'):
                os.remove(self.output_dir + '/go_statistics.xls')
            os.link(self.work_dir + '/go_statistics.xls', self.output_dir + '/go_statistics.xls')
        except:
            self.set_error("运行go_desc出错")
        self.end()

