# -*- coding: utf-8 -*-
# __author__ = 'wangbixuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess
from mbio.packages.rna.annot_config import AnnotConfig


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
            {"name": "blastout2", "type": "infile", "format": "denovo_rna_v2.blast_xml"},
            {"name": "go2level_out", "type": "outfile", "format": "annotation.go.level2"},
            {"name": "golist_out", "type": "outfile", "format": "annotation.go.go_list"},
            {"name": "blast2go_annot", "type": "outfile", "format": "annotation.go.blast2go_annot"},
            {"name": "pir_version", "type": "string", "default": "2019"}, #pir database version
            {"name": "go_version", "type": "string", "default": "2019"}, #pir database version
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
            self.option("blastout2", self.option("blastout").prop['path'])
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
            raise OptionError("必须提供BLAST结果文件", code = "32000701")

    def set_resource(self):
        self._cpu = 10
        file_size = float(os.path.getsize(self.option('blastout').prop['path'])) / 1024 / 1024
        mem = int(float(file_size)/1024 * 18) + 2
        mem = max(mem, 50)
        self._memory = "{}G".format(str(mem))

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
        # self.idmapping_db = self.config.SOFTWARE_DIR + "/database/Annotation/other/idmapping.tb"
        self.nr2go_script = self.config.PACKAGE_DIR + "/prok_rna/get_GO_from_blast_by_nr2.py"
        self.idmapping_db = AnnotConfig().get_file_path(
            file ="idmapping.tb",
            db = "pir",
            version = self.option("pir_version"))
        self.go_obo = AnnotConfig().get_file_dict(db="go", version=self.option("go_version"))['go']

    def run(self):
        super(GoAnnotationTool, self).run()
        self.run_nr2go()

    def convert_xml2table(self):
        self.logger.info('转换xml结果为表格格式')
        self.option("blastout").convert2table("blast_table.xls")

    def run_nr2go(self):
        self.convert_xml2table()
        cmd1 = 'program/Python/bin/python {}'.format(self.nr2go_script)
        cmd1 += ' %s %s %s %s' % ("blast_table.xls",
                                  self.idmapping_db,
                                  10,
                                  "blast2go.annot")
        # Config.DB_HOST,Config.DB_USER,Config.DB_PASSWD
        self.logger.info("运行NR2GO")
        self.logger.info(cmd1)
        command = self.add_command("nr2go", cmd1)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行nr2go结束")
        else:
            self.set_error('运行b2g出错', code = "32000702")

        linkfile = self.output_dir + '/blast2go.annot'
        if os.path.exists(linkfile):
            os.remove(linkfile)

        os.link(self.work_dir + '/blast2go.annot', linkfile)
        self.option('blast2go_annot', linkfile)
        self.logger.debug("b2g end")
        self.run_gomerge()


    def run_b2g(self):
        self.blast_nr_out = self.work_dir + '/temp_blast_nr'
        self.blast_nr_out_xml = self.work_dir + '/temp_blast_nr.xml'

        self.option("blastout2").convert_xml2go(self.blast_nr_out_xml)
        self.option("blastout2", self.blast_nr_out_xml)

        split_file=self.option("blastout2").change_blast_version2(self.blast_nr_out, sub_num=20000)
        a = 0
        for xml_file in split_file:
            a += 1
            cmd = '/program/sun_jdk1.8.0/bin/java -Xmx30g -cp ' + self.config.SOFTWARE_DIR + '/bioinfo/annotation/b2g4pipe_v2.5/*:'
            cmd += self.config.SOFTWARE_DIR + '/bioinfo/annotation/b2g4pipe_v2.5/ext/*: es.blast2go.prog.B2GAnnotPipe'
            cmd += ' -in {} -prop {}/bioinfo/annotation/b2g4pipe_v2.5/b2gPipe.properties -annot -out {}'.format(xml_file, self.config.SOFTWARE_DIR, xml_file + '.blast2go')
            self.logger.info('运行b2g程序输入为{}'.format(xml_file))
            self.logger.info(cmd)
            b2g = self.add_command('b2g' + str(a) , cmd)
            b2g.run()
            self.wait('b2g' + str(a))
            if b2g.return_code == 0:
                self.logger.info('运行b2g完成')
            else:
                self.set_error('运行b2g出错', code = "32000702")
        for xml_file in split_file:
            os.system("cat " + xml_file + ".blast2go.annot >> blast2go.annot")
        linkfile = self.output_dir + '/blast2go.annot'
        if os.path.exists(linkfile):
            os.remove(linkfile)

        os.link(self.work_dir + '/blast2go.annot', linkfile)
        self.option('blast2go_annot', linkfile)
        self.logger.debug("b2g end")
        self.run_gomerge()

    def run_gomerge(self):
        cmd1 = '{}/program/Python/bin/python {}/bioinfo/annotation/scripts/goMerge.py'.format(
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
            self.set_error('运行mergeGO.py出错', code = "32000703")
        self.run_annotation()

    def run_annotation(self):
        cmd2 = '{}/program/Python/bin/python {}/bioinfo/annotation/scripts/goAnnot.py'.format(
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
            self.set_error("运行goAnnot.py出错", code = "32000704")
        # self.run_gosplit()
        self.end()

    def run_gosplit(self):
        cmd3 = '{}/program/Python/bin/python {}/bioinfo/annotation/scripts/goSplit.py'.format(
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
            self.set_error("运行goSplit.py出错", code = "32000705")
        self.end()
