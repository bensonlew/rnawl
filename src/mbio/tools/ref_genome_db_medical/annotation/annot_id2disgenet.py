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


class AnnotId2disgenetAgent(Agent):
    """
    to perform Gene Ontology Annotation
    author: liubinxu
    last_modified: 20200811
    """

    def __init__(self, parent):
        super(AnnotId2disgenetAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "ids_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "known_disgenet", "type": "string", "default": None},
            {"name": "disgenet_version", "type": "string", "default": None},
        ]
        self.add_option(options)
        self.step.add_steps('disgenet_annotation')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        # self.queue = 'BLAST2GO'  # 投递到指定的队列BLAST2GO

    def step_start(self):
        self.step.disgenet_annotation.start()
        self.step.update()

    def step_end(self):
        self.step.disgenet_annotation.finish()
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
        super(AnnotId2disgenetAgent, self).end()


class AnnotId2disgenetTool(Tool):

    def __init__(self, config):
        super(AnnotId2disgenetTool, self).__init__(config)

        self.disgenet_tax = self.config.SOFTWARE_DIR + "/database/DisGeNET/7.0/browser_source_summary_gda.tsv"
        self.id_disgenet_dict = dict()
        self.gene2disgenet = dict()
        self.gene2disid = dict()
        self.gene2disname = dict()

    def run(self):
        super(AnnotId2disgenetTool, self).run()
        # self.convert_xml2table()
        self.get_disgenet_db()
        self.id2disgenet()
        self.export_disgenet()
        self.set_out()
        self.end()


    def get_disgenet_db(self):
        """
        从参考库DISGENET_sequence中找到gi号对应的description
        """

        with open(self.disgenet_tax, 'r') as f:
            f.readline()
            for line in f:
                cols = line.strip().split("\t")
                gene_id = cols[1]
                if gene_id in self.gene2disgenet:
                    if cols[0] not in self.gene2disgenet[gene_id]:
                        self.gene2disgenet[gene_id].append(cols[0])
                    else:
                        pass
                    self.gene2disid[gene_id].append(cols[7])
                    self.gene2disname[gene_id].append(cols[6])
                else:
                    self.gene2disgenet[gene_id] = [cols[0]]
                    self.gene2disid[gene_id] = [cols[7]]
                    self.gene2disname[gene_id] = [cols[6]]

    def id2disgenet(self):
        with open(self.option("ids_file").prop['path'], 'r') as f:
            for id_dict in csv.DictReader(f, delimiter='\t'):
                if id_dict['transcript_id'] in self.id_disgenet_dict:
                    pass
                else:
                    if id_dict['entrezgene_id'] != "\\N" and  id_dict['entrezgene_id'] in self.gene2disgenet:
                        enterz_id = id_dict['entrezgene_id']
                        self.id_disgenet_dict[id_dict['transcript_id']] = [self.gene2disgenet[enterz_id], self.gene2disid[enterz_id], self.gene2disname[enterz_id], enterz_id]


    def export_disgenet(self):
        with open("disgenet_annot.tsv", 'w') as fo:
            fo.write("transcript_id\tdisgenet\tdisease_id\tdisease_name\tentrez_id\n")
            for k, v in self.id_disgenet_dict.items():
                fo.write("\t".join([k,
                                    "; ".join(v[0]),
                                    "; ".join(v[1]),
                                    "; ".join(v[2]),
                                    v[3]]) + "\n")

    def set_out(self):
        outfiles = ["disgenet_annot.tsv"]
        for item in outfiles:
            linkfile = self.output_dir + '/' + item
            if os.path.exists(linkfile):
                os.remove(linkfile)
            os.link(self.work_dir + '/' + item, linkfile)


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to disgenet test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'annot_id2disgenet_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db_medical.annotation.annot_id2disgenet',
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
