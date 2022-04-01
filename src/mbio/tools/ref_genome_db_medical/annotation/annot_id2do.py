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


class AnnotId2doAgent(Agent):
    """
    to perform Gene Ontology Annotation
    author: liubinxu
    last_modified: 20200811
    """

    def __init__(self, parent):
        super(AnnotId2doAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "ids_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "known_do", "type": "string", "default": None},
            {"name": "do_version", "type": "string", "default": "202008"},
        ]
        self.add_option(options)
        self.step.add_steps('do_annotation')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        # self.queue = 'BLAST2GO'  # 投递到指定的队列BLAST2GO

    def step_start(self):
        self.step.do_annotation.start()
        self.step.update()

    def step_end(self):
        self.step.do_annotation.finish()
        self.step.update()

    def check_options(self):
        if self.option("ids_file").is_set:
            pass
        else:
            raise OptionError("必须提供id结果文件", code = "33710902")

    def set_resource(self):
        self._cpu = '2'
        self._memory = '20G'

    def end(self):
        super(AnnotId2doAgent, self).end()


class AnnotId2doTool(Tool):

    def __init__(self, config):
        super(AnnotId2doTool, self).__init__(config)

        self.do_tax = AnnotConfig().get_file_path(
            file ="enterz2do.tsv",
            db = "do",
            version = self.option("do_version"))
        self.id_do_dict = dict()
        self.gene2do = dict()

    def run(self):
        super(AnnotId2doTool, self).run()
        # self.convert_xml2table()
        self.get_do_db()
        self.id2do()
        self.export_do()
        self.set_out()
        self.end()

    def convert_xml2table(self):
        self.logger.info('转换xml结果为表格格式')
        self.option("blastout").convert2table("blast_table.xls")

    def get_do_db(self):
        """
        从参考库DO_sequence中找到gi号对应的description
        """

        with open(self.do_tax, 'r') as f:
            f.readline()
            for line in f:
                cols = line.strip().split("\t")
                gene_id = cols[0]
                if cols[1] != "NA":
                    self.gene2do[gene_id] = cols[1].split("|")

    def id2do(self):
        with open(self.option("ids_file").prop['path'], 'r') as f:
            for id_dict in csv.DictReader(f, delimiter='\t'):
                if id_dict['transcript_id'] in self.id_do_dict:
                    pass
                else:
                    if id_dict['entrezgene_id'] != "\\N" and  id_dict['entrezgene_id'] in self.gene2do:
                        self.id_do_dict[id_dict['transcript_id']] = self.gene2do[id_dict['entrezgene_id']]


    def export_do(self):
        with open("do_annot.tsv", 'w') as fo:
            fo.write("transcript_id\tdo\tdescription\n")
            for k, v in self.id_do_dict.items():
                fo.write("\t".join([k, ";".join(v)]) + "\n")

    def set_out(self):
        outfiles = ["do_annot.tsv"]
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
            'id': 'annot_id2do_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db_medical.annotation.annot_id2do',
            'instant': False,
            'options': dict(
                ids_file = "/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/ids.tsv"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
