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


class AnnotId2keggAgent(Agent):
    """
    to perform Gene Ontology Annotation
    author: liubinxu
    last_modified: 20200811
    """

    def __init__(self, parent):
        super(AnnotId2keggAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "ids_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "species", "type": "string", "default": None},
            {"name": "kegg_version", "type": "string", "default": "202007"},
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
        self._cpu = '2'
        self._memory = '20G'

    def end(self):
        super(AnnotId2keggAgent, self).end()


class AnnotId2keggTool(Tool):
    def __init__(self, config):
        super(AnnotId2keggTool, self).__init__(config)

        self.abr_dict = dict({
            "Homo_sapiens": "hsa",
            "Mus_musculus": "mmu",
            "Rattus_norvegicus": "rno"
        })

        # self.id2kegg = self.config.PACKAGE_DIR + "/ref_genome_db_medical//id2kegg.py"


    def run(self):
        super(AnnotId2keggTool, self).run()
        # self.convert_xml2table()
        self.get_kegg_db()
        self.export_kegg()
        self.set_out()
        self.end()

    def convert_xml2table(self):
        self.logger.info('转换xml结果为表格格式')
        self.option("blastout").convert2table("blast_table.xls")

    def get_kegg_db(self):
        """
        从参考库KEGG_sequence中找到gi号对应的description
        """
        self.gene2ko = dict()
        keggdb = AnnotConfig().get_file_path(
            file="data",
            db="kegg",
            db_type="file",
            version=self.option("kegg_version"))

        kegg_faa = "{}/{}.faa".format(keggdb, self.abr_dict[self.option("species")])

        with open(kegg_faa, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    cols = line.strip().lstrip(">").split(":")
                    self.gene2ko[cols[-1] + ":" + cols[0]] = cols[1]
        return self.gene2ko

    def export_kegg(self):
        with open(self.option("ids_file").prop['path'], 'r') as f, open("kegg_annot.tsv", 'w') as fo:
            fo.write("transcript_id\tkegg_gene\tko\n")
            for id_dict in csv.DictReader(f, delimiter='\t'):
                if id_dict['kegg'] in self.gene2ko:
                    kegg_gene = id_dict['kegg']
                    ko = self.gene2ko[kegg_gene]
                    fo.write("\t".join([id_dict["transcript_id"], kegg_gene, ko]) + "\n")
                elif self.abr_dict[self.option('species')] + ":" + id_dict['entrezgene_id'] in self.gene2ko:
                    kegg_gene = self.abr_dict[self.option('species')] + ":" + id_dict['entrezgene_id']
                    ko = self.gene2ko[kegg_gene]
                    fo.write("\t".join([id_dict["transcript_id"], kegg_gene, ko]) + "\n")
                else:
                    pass

    def set_out(self):
        outfiles = ["kegg_annot.tsv"]
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
            'id': 'annot_id2kegg_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db_medical.annotation.annot_id2kegg',
            'instant': False,
            'options': dict(
                ids_file = "/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/ids.tsv",
                species = "Homo_sapiens"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
