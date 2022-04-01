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
            {"name": "go2level_out", "type": "outfile", "format": "annotation.go.level2"},
            {"name": "golist_out", "type": "outfile", "format": "annotation.go.go_list"},
            {"name": "blast2go_annot", "type": "outfile", "format": "annotation.go.blast2go_annot"},
            {"name": "protein_fasta", "type": "infile", "format": "itraq_and_tmt.common"},
            {"name": "merge_known", "type": "bool", "default": True},
            {"name": "version", "type": "string", "default": "2019"},
            {"name": "go_version", "type": "string", "default": "2019"}, #pir database version
        ]
        self.add_option(options)
        self._memory_increase_step = 50
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
        elif self.option("blast2go_annot").is_set:
            pass
        else:
            raise OptionError("必须提供BLAST结果文件", code = "33701101")

    def set_resource(self):
        self._cpu = 2
        self._memory = '100G'

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
        self.merge_go = self.config.PACKAGE_DIR + "/ref_rna_v2/goMerge.py"
        self._version = "1.0"
        self.b2g_user = "biocluster102"
        self.b2g_password = "sanger-dev-123"
        self.idmapping_db = AnnotConfig().get_file_path(
            file ="idmapping.tb",
            db = "pir",
            version = self.option("version"))
        self.go_obo = AnnotConfig().get_file_dict(db="go", version=self.option("go_version"))['go']

        self.nr2go_script = self.config.PACKAGE_DIR + "/prok_rna/get_GO_from_blast_by_nr2.py"
        self.go_annotation_py = os.path.join(self.config.PACKAGE_DIR, 'rna/annotation/go_annotation2.py')

    def run(self):
        super(GoAnnotationTool, self).run()
        self.convert_xml2table()
        self.nr2go()
        self.run_gomerge()

    def convert_xml2table(self):
        self.logger.info('转换xml结果为表格格式')
        self.option("blastout").convert2table("blast_table.xls")

    def nr2go(self):
        cmd1 = 'program/Python/bin/python {}'.format(self.nr2go_script)
        cmd1 += ' %s %s %s %s' % ("blast_table.xls",
                                  self.idmapping_db,
                                  10,
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

    def run_gomerge(self):
        self.option('blast2go_annot', os.path.abspath("blast2go_annot.xls"))
        cmd1 = 'program/Python/bin/python {}'.format(self.merge_go)
        cmd1 += ' %s %s' % (self.option('blast2go_annot').prop['path'], 'GO.list')
        # cmd1 += ' %s %s' % ("blast2go_annot.xls", 'GO.list')
        # Config.DB_HOST,Config.DB_USER,Config.DB_PASSWD
        self.logger.info("运行mergeGO.py")
        self.logger.info(cmd1)
        command = self.add_command("go_merge", cmd1, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行mergeGO成功")
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("运行mergeGO.py出错", code = "33701102")
        '''
        try:
            subprocess.check_output(cmd1, shell=True)
        except subprocess.CalledProcessError:
            self.set_error('运行mergeGO.py出错')
        '''
        if os.path.exists(self.output_dir + '/GO.list'):
            os.remove(self.output_dir + '/GO.list')
        if os.path.exists(self.output_dir + '/query_gos.list'):
            os.remove(self.output_dir + '/query_gos.list')
        os.link(self.work_dir + '/GO.list',
                self.output_dir + '/query_gos.list')
        self.option('golist_out', self.output_dir + '/query_gos.list')
        self.run_annotation()

    def run_annotation(self):
        cmd2 = 'program/Python/bin/python {}'.format(self.go_annotation_py)
        cmd2 += ' %s %s %s' % (
            self.work_dir + '/GO.list', self.work_dir, self.go_obo)  # 10.100.203.193
        self.logger.info("运行go_annotation")
        self.logger.info(cmd2)
        command = self.add_command("go_annot", cmd2, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行go_annot成功")
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("运行go_annot出错", code = "33701103")
        '''
        try:
            subprocess.check_output(cmd2, shell=True)
            self.logger.info("运行goAnnot.py完成")
        except subprocess.CalledProcessError:
            self.set_error("运行goAnnot.py出错")
        # self.run_gosplit()
        '''
        infiles = ["level4.stat.tsv", 'level3.stat.tsv', 'level2.stat.tsv']
        outfiles = ['go1234level_statistics.xls', 'go123level_statistics.xls', 'go12level_statistics.xls']
        for inf,item in zip(infiles, outfiles):
            linkfile = self.output_dir + '/' + item
            if os.path.exists(linkfile):
                os.remove(linkfile)
            os.link(self.work_dir + '/' + inf, linkfile)
        self.end()

    def run_gosplit(self):
        cmd3 = 'program/Python/bin/python {}/bioinfo/annotation/scripts/goSplit.py'.format(self.config.SOFTWARE_DIR)
        cmd3 += ' %s' % self.work_dir + '/go_detail.xls'
        self.logger.info("运行goSplit.py")
        self.logger.info(cmd3)
        command = self.add_command("go_annot", cmd3, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行go_split成功")
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("运行go_split出错", code = "33701104")
        outfiles = ['go2level.xls', 'go3level.xls', 'go4level.xls']
        for item in outfiles:
            linkfile = self.output_dir + '/' + item
            if os.path.exists(linkfile):
                os.remove(linkfile)
            os.link(self.work_dir + '/' + item, linkfile)
        self.option('go2level_out', self.output_dir + '/go2level.xls')
        '''
        try:
            subprocess.check_output(cmd3, shell=True)
            self.logger.info("运行goSplit.py完成")
        except subprocess.CalledProcessError:
            self.set_error("运行goSplit.py出错")
        '''
        self.end()
