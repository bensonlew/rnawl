# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess
from mbio.packages.rna.annot_config import AnnotConfig


class Blast2goAgent(Agent):
    """
    to perform Gene Ontology Annotation
    author: liubinxu
    last_modified: 20180509
    """

    def __init__(self, parent):
        super(Blast2goAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "blastout_go", "type": "outfile", "format": "ref_rna_v2.blast_xml"},
            {"name": "blast2go_annot", "type": "outfile", "format": "ref_rna_v2.blast2go_annot"},
            {"name": "known_go", "type": "string", "default": None},
            {"name": "pir_version", "type": "string", "default": "2019"}, #pir database version
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
        else:
            raise OptionError("必须提供BLAST结果文件", code = "33700201")

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
        ])
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

    def run(self):
        super(Blast2goTool, self).run()
        self.run_b2g()
        self.convert_go()
        self.end()

    def run_b2g(self):
        self.blast_nr_out = self.work_dir + '/temp_blast_nr.xml'
        self.option("blastout").convert_xml2go(self.blast_nr_out)
        self.option("blastout_go", self.blast_nr_out)
        cmd = '/program/sun_jdk1.8.0/bin/java -Xmx50g -cp ' + self.config.SOFTWARE_DIR + '/bioinfo/annotation/b2g4pipe_v2.5/*:'
        cmd += self.config.SOFTWARE_DIR + '/bioinfo/annotation/b2g4pipe_v2.5/ext/*: es.blast2go.prog.B2GAnnotPipe'
        cmd += ' -in {} -prop {}/bioinfo/annotation/b2g4pipe_v2.5/b2gPipe.properties -annot -out {}'.format(self.blast_nr_out, self.config.SOFTWARE_DIR, self.work_dir + '/blast2go')
        self.logger.info('运行b2g程序')
        self.logger.info(cmd)
        b2g = self.add_command('b2g', cmd, ignore_error=True)
        b2g.run()
        self.wait('b2g')
        if b2g.return_code == 0:
            self.logger.info('运行b2g完成')
            linkfile = self.output_dir + '/blast2go.annot'
            if os.path.exists(linkfile):
                os.remove(linkfile)
            os.link(self.work_dir + '/blast2go.annot', linkfile)
            self.logger.debug("b2g end")
        elif b2g.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error('运行b2g出错', code = "33700202")

    def get_known_go(self, go_file):
        '''
        获取已知的GO注释,
        ; 返回dict {tran_id: go1;go2}
        '''
        tran2go = dict()
        with open(go_file, 'rb') as go:
            for line in go.readlines():
                cols = line.strip().split("\t")
                if len(cols) >= 3 and cols[2]:
                    tran2go[cols[1]] = cols[2]
        return tran2go

    def convert_go(self):
        # 输出注释，对应blast相关evalue, identity, similarity
        self.logger.info('转换go结果格式')
        self.option("blastout_go").convert2table("blast_table.xls")
        blast_dict = dict()
        with open("blast_table.xls", 'rb') as blast_file:
            blast_file.readline()
            for line in blast_file.readlines():
                cols = line.strip().split("\t")
                query_hit = cols[5]
                blast_dict[query_hit] = cols

        trans_list = list()
        with open(self.work_dir + '/blast2go.annot', 'rb') as b2g, open(self.work_dir + '/blast2go_annot.xls', 'w') as b2gout:
            for line in b2g.readlines():
                line = line.strip().split("\t")
                # trans_list.append(line[0])
                blast_list = blast_dict[line[0]]
                if len(line) == 3:
                    b2gout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        line[0].split("__id__")[0],
                        line[1],
                        line[2],
                        line[0].split("__id__")[1],
                        blast_list[1],
                        blast_list[3],
                        blast_list[4],
                    ))
                elif len(line) == 2:
                    b2gout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        line[0].split("__id__")[0],
                        line[1],
                        "",
                        line[0].split("__id__")[1],
                        blast_list[1],
                        blast_list[3],
                        blast_list[4],
                    ))
                else:
                    pass
            trans_set = set(trans_list)
            if self.option("known_go"):
                tran2go = self.get_known_go(self.option("known_go"))
                for tran,gos in tran2go.items():
                    for go in gos.split(";"):
                        b2gout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            tran,go,"known","known",0,100,100
                        ))

        self.logger.info('转换格式完成')
        linkfile = self.output_dir + '/blast2go_annot.xls'
        if os.path.exists(linkfile):
            os.remove(linkfile)
        os.link(self.work_dir + '/blast2go_annot.xls', linkfile)
        self.option('blast2go_annot', linkfile)
