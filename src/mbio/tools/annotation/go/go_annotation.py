# -*- coding: utf-8 -*-
# __author__ = 'wangbixuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess


class GoAnnotationAgent(Agent):
    """
    to perform Gene Ontology Annotation
    author: wangbixuan
    last_modified: 20160728
    """

    def __init__(self, parent):
        super(GoAnnotationAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "align.blast.blast_xml"},
            {"name": "go2level_out", "type": "outfile", "format": "annotation.go.level2"},
            {"name": "golist_out", "type": "outfile", "format": "annotation.go.go_list"},
            {"name": "blast2go_annot", "type": "outfile", "format": "annotation.go.blast2go_annot"}
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
            raise OptionError("必须提供BLAST结果文件" ,code="31201501")

    def set_resource(self):
        self._cpu = 10
        self._memory = '60G'

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
        super(GoAnnotationAgent, self).end()


class GoAnnotationTool(Tool):

    def __init__(self, config):
        super(GoAnnotationTool, self).__init__(config)
        self.set_environ(JAVA_HOME=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0")
        self.set_environ(JRE_HOME=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/jre")
        self.set_environ(CLASSPATH=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/lib")
        self.set_environ(CLASSPATH=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/jre/lib")
        self._version = "1.0"
        self.b2g_user = "biocluster102"
        self.b2g_password = "sanger-dev-123"

    def run(self):
        super(GoAnnotationTool, self).run()
        self.run_b2g()

    def run_b2g(self):
        self.blast_nr_out = self.work_dir + '/temp_blast_nr.xml'
        self.option("blastout").change_blast_version(self.blast_nr_out)
        cmd = '/program/sun_jdk1.8.0/bin/java -Xmx30g -cp ' + self.config.SOFTWARE_DIR + '/bioinfo/annotation/b2g4pipe_v2.5/*:'
        cmd += self.config.SOFTWARE_DIR + '/bioinfo/annotation/b2g4pipe_v2.5/ext/*: es.blast2go.prog.B2GAnnotPipe'
        cmd += ' -in {} -prop {}/bioinfo/annotation/b2g4pipe_v2.5/b2gPipe.properties -annot -out {}'.format(self.blast_nr_out, self.config.SOFTWARE_DIR, self.work_dir + '/blast2go')
        self.logger.info('运行b2g程序')
        self.logger.info(cmd)
        b2g = self.add_command('b2g', cmd)
        b2g.run()
        self.wait('b2g')
        if b2g.return_code == 0:
            self.logger.info('运行b2g完成')
            linkfile = self.output_dir + '/blast2go.annot'
            if os.path.exists(linkfile):
                os.remove(linkfile)
            os.link(self.work_dir + '/blast2go.annot', linkfile)
            self.option('blast2go_annot', linkfile)
            self.logger.debug("b2g end")
            if os.path.exists(self.work_dir + '/blast2go.annot'):
                if os.path.getsize(self.work_dir + '/blast2go.annot') != 0:
                    self.run_gomerge()
                else:
                    self.end()
            else:
                self.end()
        else:
            self.set_error('运行b2g出错', code="31201501")

    def run_gomerge(self):
        cmd1 = '{}/miniconda2/bin/python {}/bioinfo/annotation/scripts/goMerge.py'.format(
            self.config.SOFTWARE_DIR, self.config.SOFTWARE_DIR)
        cmd1 += ' %s %s' % (
            self.work_dir + '/blast2go.annot', 'GO.list')
        # Config.DB_HOST,Config.DB_USER,Config.DB_PASSWD
        self.logger.info("运行mergeGO.py")
        self.logger.info(cmd1)
        try:
            subprocess.check_output(cmd1, shell=True)
            if os.path.exists(self.output_dir + '/GO.list'):
                os.remove(self.output_dir + '/GO.list')
            if os.path.exists(self.output_dir + '/query_gos.list'):
                os.remove(self.output_dir + '/query_gos.list')
            os.link(self.work_dir + '/GO.list',
                    self.output_dir + '/query_gos.list')
            self.option('golist_out', self.output_dir + '/query_gos.list')
        except subprocess.CalledProcessError:
            self.set_error('运行mergeGO.py出错', code="31201502")
        self.run_annotation()

    def run_annotation(self):
        cmd2 = '{}/miniconda2/bin/python {}/bioinfo/annotation/scripts/goAnnot.py'.format(
            self.config.SOFTWARE_DIR, self.config.SOFTWARE_DIR)
        cmd2 += ' %s %s %s %s' % (
            self.work_dir + '/GO.list', 'localhost', self.b2g_user, self.b2g_password)  # 10.100.203.193
        self.logger.info("运行goAnnot.py")
        self.logger.info(cmd2)
        try:
            subprocess.check_output(cmd2, shell=True)
            self.logger.info("运行goAnnot.py完成")
            outfiles = ['go1234level_statistics.xls', 'go123level_statistics.xls', 'go12level_statistics.xls']
            for item in outfiles:
                linkfile = self.output_dir + '/' + item
                if os.path.exists(linkfile):
                    os.remove(linkfile)
                os.link(self.work_dir + '/' + item, linkfile)
        except subprocess.CalledProcessError:
            self.set_error("运行goAnnot.py出错", code="31201503")
        # self.run_gosplit()
        self.run_not_level()

    def run_gosplit(self):
        cmd3 = '{}/miniconda2/bin/python {}/bioinfo/annotation/scripts/goSplit.py'.format(
            self.config.SOFTWARE_DIR, self.config.SOFTWARE_DIR)
        cmd3 += ' %s' % self.work_dir + '/go_detail.xls'
        self.logger.info("运行goSplit.py")
        self.logger.info(cmd3)
        try:
            subprocess.check_output(cmd3, shell=True)
            self.logger.info("运行goSplit.py完成")
            outfiles = ['go2level.xls', 'go3level.xls', 'go4level.xls']
            for item in outfiles:
                linkfile = self.output_dir + '/' + item
                if os.path.exists(linkfile):
                    os.remove(linkfile)
                os.link(self.work_dir + '/' + item, linkfile)
            self.option('go2level_out', self.output_dir + '/go2level.xls')
        except subprocess.CalledProcessError:
            self.set_error("运行goSplit.py出错", code="31201504")
        self.end()

    def run_not_level(self):
        cmd = '{}/miniconda2/bin/python {}/annotation/go/go_desc.py'.format(self.config.SOFTWARE_DIR, self.config.PACKAGE_DIR)
        cmd += ' {} {}'.format(self.work_dir + '/blast2go.annot', self.work_dir + '/go_statistics.xls' )
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

