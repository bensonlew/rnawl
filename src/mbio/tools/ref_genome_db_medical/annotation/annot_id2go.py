# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess
from mbio.packages.rna.annot_config import AnnotConfig
import csv
import unittest
import sqlite3
import gzip


class AnnotId2goAgent(Agent):
    """
    to perform Gene Ontology Annotation
    author: liubinxu
    last_modified: 20180509
    """

    def __init__(self, parent):
        super(AnnotId2goAgent, self).__init__(parent)
        options = [
            {"name": "ids_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "blast2go_annot", "type": "infile", "format": "ref_rna_v2.blast2go_annot"},
            {"name": "pir_version", "type": "string", "default": "2019"},
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
        if self.option("ids_file").is_set:
            pass
        else:
            raise OptionError("必须提供id结果文件", code = "33710902")

    def set_resource(self):
        self._cpu = self.option("p") + 1
        self._memory = '60G'

    def end(self):
        super(AnnotId2goAgent, self).end()


class AnnotId2goTool(Tool):

    def __init__(self, config):
        super(AnnotId2goTool, self).__init__(config)
        self.idmapping_db = self.idmapping_db = AnnotConfig().get_file_path(
            file ="idmapping.tb",
            db = "pir",
            version = self.option("pir_version"))
        self.nr2go_script = self.config.PACKAGE_DIR + "/ref_genome_db_medical//get_go_byid.py"
        self.blast_gos = dict()

    def run(self):
        super(AnnotId2goTool, self).run()
        self.get_blast2go()
        self.convert_id2table()
        self.id2go()
        self.merge_idgo_blastgo()
        self.set_out()
        self.end()

    def get_blast2go(self):
        with open(self.option("blast2go_annot").prop['path'], 'r') as f:
            for line in f:
                cols = line.strip().split("\t")
                if cols[0] in self.blast_gos:
                    self.blast_gos[cols[0]] += ";" + cols[1]
                self.blast_gos[cols[0]] = cols[1]

    def convert_id2table(self):
        self.logger.info('转换xml结果为表格格式')
        with open(self.option("ids_file").prop['path'], 'r') as f, open("temp_ids.tsv", 'w') as fo:
            for id_dict in csv.DictReader(f, delimiter='\t'):
                fo.write("\t".join([id_dict['transcript_id'], id_dict['uniprot_gn_id']]) + "\n")

    def id2go(self):
        cmd1 = 'program/Python/bin/python {}'.format(self.nr2go_script)
        cmd1 += ' %s %s %s %s' % ("temp_ids.tsv",
                                  self.idmapping_db,
                                  self.option("p"),
                                  "id2go_annot.xls")
        # Config.DB_HOST,Config.DB_USER,Config.DB_PASSWD
        self.logger.info("运行id2GO")
        self.logger.info(cmd1)
        command = self.add_command("nr2go", cmd1)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行nr2go结束")
        else:
            self.set_error("运行nr2go出错", code = "33710902")

    def merge_idgo_blastgo(self):
        with open('id2go_annot.xls', 'r') as f, open('all2go_annot.xls', 'w') as fo:
            fo.write('#Seq_id\tGos\tsource\n')
            idmapping_genes = list()
            for line in f:
                fo.write('{}\tid_mapping\n'.format(line.strip()))
                idmapping_genes.append(line.split("\t")[0])
            for k,v in self.blast_gos.items():
                if k in idmapping_genes:
                    pass
                else:
                    fo.write('{}\t{}\tblast\n'.format(k, v))

    def set_out(self):
        outfiles = ["id2go_annot.xls", 'all2go_annot.xls']
        for item in outfiles:
            linkfile = self.output_dir + '/' + item
            if os.path.exists(linkfile):
                os.remove(linkfile)
            os.link(self.work_dir + '/' + item, linkfile)


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'annot_id2go_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db_medical.annotation.annot_id2go',
            'instant': False,
            'options': dict(
                ids_file = "/mnt/ilustre/users/sanger-dev/workspace/20200811/Single_annot_getid_5420_4347/AnnotGetid/ids.tsv"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
