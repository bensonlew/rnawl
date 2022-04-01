# -*- coding: utf-8 -*-
# __author__ = 'wangbixuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess
from Bio import SeqIO
import sqlite3
from mbio.packages.annotation.go.go_terms import Terms

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
            raise OptionError("必须提供BLAST结果文件", code = "32500801")
        if self.option("merge_known"):
            if not self.option("protein_fasta").is_set:
                raise OptionError("合并已知注释需要提供fasta文件", code = "32500802")
            else:
                pass

    def set_resource(self):
        self._cpu = 10
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
        self._version = "1.0"
        self.b2g_user = "biocluster102"
        self.b2g_password = "sanger-dev-123"
        self.acc2go_db = self.config.SOFTWARE_DIR + "/database/Annotation/latest_sqlite/acc2go.db"

    def run(self):
        super(GoAnnotationTool, self).run()
        self.run_b2g()

    def run_b2g(self):
        self.blast_nr_out = self.work_dir + '/temp_blast_nr.xml'
        self.option("blastout").change_blast_version(self.blast_nr_out)
        cmd = '/program/sun_jdk1.8.0/bin/java -Xmx100g -cp ' + self.config.SOFTWARE_DIR + '/bioinfo/annotation/b2g4pipe_v2.5/*:'
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
            self.run_gomerge()
        else:
            self.set_error('运行b2g出错', code = "32500803")

    def get_set(self, fa_file):
        '''
        获取序列名列表
        in: fa_file fasta path
        return: seq_id list
        '''
        seq_list = list()
        for seq in SeqIO.parse(fa_file, "fasta"):
            seq_list.append(seq.id)
        return seq_list

    def merge_known(self, seq_list, annot_file, annot_out_file):
        '''
        合并go注释与数据库已存注释
        in: seq_list seq_id list
        in: annot_file GO annotation file
        out: annot_out_file merged file
        '''
        anno_dict = dict()
        with open(annot_file, 'rb') as anno_f:
            for line in anno_f.readlines():
                cols = line.strip().split("\t")
                if anno_dict.has_key(cols[0]):
                    anno_dict[cols[0]] += ";" + cols[1]
                else:
                    anno_dict[cols[0]] = cols[1]
        conn = sqlite3.connect(self.acc2go_db)
        cursor = conn.cursor()
        with open(annot_out_file, 'wb') as anno_out_f:
            for seq in seq_list:
                try:
                    cursor.execute('select * from acc2go where acc="{}"'.format(seq))
                    go = cursor.fetchall()[0][1]
                except:
                    go = ""
                    self.logger.info("没找到seq_id:{}".format(seq))
                if anno_dict.has_key(seq):
                    if go == "":
                        go = anno_dict[seq]
                    else:
                        go = go + ";" + anno_dict[seq]
                if go == "":
                    pass
                else:
                    anno_out_f.write("{}\t{}\n".format(seq, ";".join(list(set(go.split(";"))))))

    def go_merge_parent(self):
        '''
        合并GO父节点注释
        '''
        GO_db = Terms(obo_fp=self.config.SOFTWARE_DIR + "/database/GO/go-basic.obo")
        with open("GO.list", 'r') as go_in, open("GO_P.list", "w") as go_out:
            for line in go_in:
                cols = line.strip().split("\t")
                terms = cols[1].split(";")
                terms_parents = GO_db.get_list_parent_id(terms)
                go_out.write("{}\t{}\n".format(cols[0], ";".join(terms_parents)))

        if os.path.exists('GO.list'):
            os.remove('GO.list')
        os.system('cp GO_P.list GO.list')

    def run_gomerge(self):
        cmd1 = '{}/program/Python/bin/python {}/bioinfo/annotation/scripts/goMerge.py'.format(
            self.config.SOFTWARE_DIR, self.config.SOFTWARE_DIR)
        cmd1 += ' %s %s' % (
            self.work_dir + '/blast2go.annot', 'GO_annot.list')

        # Config.DB_HOST,Config.DB_USER,Config.DB_PASSWD
        self.logger.info("运行mergeGO.py")
        self.logger.info(cmd1)
        try:
            subprocess.check_output(cmd1, shell=True)

        except subprocess.CalledProcessError:
            self.set_error('运行mergeGO.py出错', code = "32500804")

        if self.option("merge_known"):
            seq_list = self.get_set(self.option('protein_fasta').prop['path'])
            self.merge_known(seq_list, "GO_annot.list", "GO.list")
        else:
            if os.path.exists('GO.list'):
                os.remove('GO.list')
            os.link('GO_annot.list', 'GO.list')
        self.go_merge_parent()

        if os.path.exists(self.output_dir + '/GO.list'):
            os.remove(self.output_dir + '/GO.list')
        if os.path.exists(self.output_dir + '/query_gos.list'):
            os.remove(self.output_dir + '/query_gos.list')
        os.link(self.work_dir + '/GO.list',
            self.output_dir + '/query_gos.list')
        self.option('golist_out', self.output_dir + '/query_gos.list')
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
            self.set_error("运行goAnnot.py出错", code = "32500805")
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
            self.set_error("运行goSplit.py出错", code = "32500806")
        self.end()
