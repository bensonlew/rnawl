# -*- coding: utf-8 -*-
# __author__ = 'wangbixuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess
import pandas as pd
from mbio.packages.rna.annot_config import AnnotConfig

class GoAnnotationAgent(Agent):
    """
    to perform Gene Ontology Annotation
    author: wangbixuan
    last_modified: liubinxu 20180807
    """

    def __init__(self, parent):
        super(GoAnnotationAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "prok_rna.blast_xml"},
            {"name": "go2level_out", "type": "outfile", "format": "prok_rna.level2"},
            {"name": "go_detail", "type": "outfile", "format": "prok_rna.common"},
            {"name": "golist_out", "type": "outfile", "format": "prok_rna.go_list"},
            {"name": "blast2go_annot", "type": "infile", "format": "prok_rna.blast2go_annot"},
            {"name": "pir_version", "type": "string", "default": "2019"}, #pir database version
            {"name": "go_version", "type": "string", "default": "2019"}, #pir database version
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
        elif self.option("blast2go_annot").is_set:
            pass
        else:
            raise OptionError("必须提供BLAST结果文件", code = "35000701")

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

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
        self.merge_go = self.config.PACKAGE_DIR + "/prok_rna/goMerge.py"

        self.idmapping_db = AnnotConfig().get_file_path(
            file ="idmapping.tb",
            db = "pir",
            version = self.option("pir_version"))
        self.go_obo = AnnotConfig().get_file_dict(db="go", version=self.option("go_version"))['go']
        self.annot_go = os.path.join(self.config.PACKAGE_DIR, 'rna/annotation/go_annotation2.py')
        self.annot_go2 = os.path.join(self.config.PACKAGE_DIR, 'prok_rna/go_last_level.py')
        # self._version = "1.0"
        # self.b2g_user = "biocluster102"
        # self.b2g_password = "sanger-dev-123"

    def run(self):
        super(GoAnnotationTool, self).run()
        self.run_gomerge()

    def run_gomerge(self):
        cmd1 = 'program/Python/bin/python {}'.format(self.merge_go)
        cmd1 += ' %s %s' % (self.option('blast2go_annot').prop['path'], 'GO.list')
        # Config.DB_HOST,Config.DB_USER,Config.DB_PASSWD
        self.logger.info("运行mergeGO.py")
        self.logger.info(cmd1)
        command = self.add_command("merge_go", cmd1)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行merge_go结束")
        else:
            self.set_error("运行merge_go出错")

        self.run_annotation()

    def run_annotation(self):
        cmd2 = 'program/Python/bin/python {}'.format(self.annot_go)
        cmd2 += ' %s  %s %s' % (
            self.work_dir + '/GO.list', self.work_dir, self.go_obo)

        self.logger.info("运行goAnnot.py")
        self.logger.info(cmd2)
        command = self.add_command("annot_go", cmd2)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行annot_go结束")
        else:
            self.set_error("运行annot_go出错")
        infiles = ["level4.stat.tsv", 'level3.stat.tsv', 'level2.stat.tsv']
        outfiles = ['go1234level_statistics.xls', 'go123level_statistics.xls', 'go12level_statistics.xls']
        for inf,item in zip(infiles, outfiles):
            linkfile = self.output_dir + '/' + item
            if os.path.exists(linkfile):
                os.remove(linkfile)
            os.link(self.work_dir + '/' + inf, linkfile)
        self.run_annotation2()
        self.run_go_stat()
        self.set_out()

    def run_annotation2(self):
        cmd2 = 'program/Python/bin/python {}'.format(self.annot_go2)
        cmd2 += ' %s  %s %s' % (
            self.work_dir + '/GO.list', self.work_dir, self.go_obo)

        self.logger.info("运行goAnnot2")
        self.logger.info(cmd2)
        command = self.add_command("annot_go2", cmd2)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行annot_go2结束")
        else:
            self.set_error("运行annot_go2出错")
        infiles = ["go_detail_stat.tsv", 'go_detail_des.tsv']
        outfiles = ["go_detail_stat.tsv", 'go_detail_des.tsv']
        for inf,item in zip(infiles, outfiles):
            linkfile = self.output_dir + '/' + item
            if os.path.exists(linkfile):
                os.remove(linkfile)
            os.link(self.work_dir + '/' + inf, linkfile)


    def run_go_stat(self):
        '''
        统计GO数量
        '''
        with open(self.work_dir + '/GO.list', 'r') as go_f:
            genes = [line.strip().split("\t")[0] for line in go_f.readlines()]
            gene_num = len(set(genes))
        go_class = pd.read_table(self.output_dir + '/go12level_statistics.xls')
        type_num = len(set(go_class[go_class.columns[0]]))
        class_gene_num = {}
        for gotype in ['molecular_function', 'cellular_component', 'biological_process']:
            seqs_list = list(go_class[go_class[go_class.columns[0]]==gotype]['Seq List'])
            all_seqs = ';'.join(seqs_list).split(';')
            class_gene_num[gotype] = len(set(all_seqs))
        with open(self.work_dir + '/GO.stat.xls', 'w') as gostat_f:
            gostat_f.write("GO Category No.\t{}\nGene No. of CC\t{}\nGene No. of MF\t{}\nGene No. of BP\t{}\nGene No.\t{}\n".format(type_num, class_gene_num['cellular_component'], class_gene_num['molecular_function'], class_gene_num['biological_process'], gene_num))

    def set_out(self):
        # '''
        outfiles = [ 'GO.list', 'go_detail.xls', 'GO.stat.xls']
        for item in outfiles:
            linkfile = self.output_dir + '/' + item
            if os.path.exists(linkfile):
                os.remove(linkfile)
            os.link(self.work_dir + '/' + item, linkfile)
        # '''
        self.option('go2level_out', self.output_dir + '/go12level_statistics.xls')
        self.option('go_detail', self.output_dir + '/go_detail.xls')
        self.option('golist_out', self.output_dir + '/GO.list')
        if os.path.exists(self.output_dir + '/query_gos.list'):
            os.remove(self.output_dir + '/query_gos.list')
        os.link(self.work_dir + '/GO.list',
                self.output_dir + '/query_gos.list')
        self.option('golist_out', self.output_dir + '/query_gos.list')
        self.end()
