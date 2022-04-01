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


class AnnotId2eggnogAgent(Agent):
    """
    to perform Gene Ontology Annotation
    author: liubinxu
    last_modified: 20200811
    """

    def __init__(self, parent):
        super(AnnotId2eggnogAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "ids_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "known_eggnog", "type": "string", "default": None},
            {"name": "eggnog_version", "type": "string", "default": "202006"},
        ]
        self.add_option(options)
        self.step.add_steps('eggnog_annotation')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        # self.queue = 'BLAST2GO'  # 投递到指定的队列BLAST2GO

    def step_start(self):
        self.step.eggnog_annotation.start()
        self.step.update()

    def step_end(self):
        self.step.eggnog_annotation.finish()
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
        super(AnnotId2eggnogAgent, self).end()


class AnnotId2eggnogTool(Tool):

    def __init__(self, config):
        super(AnnotId2eggnogTool, self).__init__(config)

        self.eggnog_tax = AnnotConfig().get_file_path(
            file ="40674_members.tsv.gz",
            db = "eggnog",
            version = self.option("eggnog_version"))
        self.id_eggnog_dict = dict()
        self.gene2eggnog = dict()
        self.blast_id = dict()

    def run(self):
        super(AnnotId2eggnogTool, self).run()
        self.convert_xml2table()
        self.get_eggnog_tax_db()
        self.id2eggnog()
        self.export_eggnog()
        self.set_out()
        self.end()

    def convert_xml2table(self):
        self.logger.info('转换xml结果为表格格式')
        self.option("blastout").convert2table("blast_table.xls")
        with open('blast_table.xls', 'r') as f:
            for id_dict in csv.DictReader(f, delimiter='\t'):
                if id_dict['Query-Name'] in self.blast_id:
                    self.blast_id[id_dict['Query-Name']].append(id_dict['Hit-Name'])
                else:
                    self.blast_id[id_dict['Query-Name']] = [id_dict['Hit-Name']]

    def get_multi_id2eggnog(self, seq_ids):
        eggnog_list = list()
        for seq_id in seq_ids:
            id_clean = ".".join(seq_id.split(".")[1:])
            if id_clean in self.gene2eggnog:
                eggnog_list.extend(self.gene2eggnog[id_clean])
        return list(set(eggnog_list))

    def get_eggnog_tax_db(self):
        """
        从参考库EGGNOG_sequence中找到gi号对应的description
        """
        with gzip.open(self.eggnog_tax, 'r') as f:
            for line in f:
                cols = line.strip().split("\t")
                genes = cols[4].split(",")
                for gene in genes:
                    gene_id = ".".join(gene.split(".")[1:])
                    if gene_id in self.gene2eggnog:
                        self.gene2eggnog[gene_id].add(cols[1])
                    else:
                        self.gene2eggnog[gene_id] = set([cols[1]])

    def id2eggnog(self):
        with open("gene2eggnog.dict", 'w') as fo:
            fo.write("{}".format(self.gene2eggnog))
        with open(self.option("ids_file").prop['path'], 'r') as f:
            for id_dict in csv.DictReader(f, delimiter='\t'):
                if id_dict['transcript_id'] in self.id_eggnog_dict:
                    pass
                else:
                    if id_dict['protein_id'] != "\\N" and  id_dict['protein_id'] in self.gene2eggnog:
                        self.id_eggnog_dict[id_dict['transcript_id']] = self.gene2eggnog[id_dict['protein_id']]
                    elif 'refseq_peptide' in id_dict and id_dict['refseq_peptide'] != "\\N" and  id_dict['refseq_peptide'] in self.gene2eggnog:
                        self.id_eggnog_dict[id_dict['transcript_id']] = self.gene2eggnog[id_dict['refseq_peptide']]
                    elif 'refseq_peptide' in id_dict and id_dict['refseq_peptide'] != "\\N" and  id_dict['refseq_peptide_predicted'] in self.gene2eggnog:
                        self.id_eggnog_dict[id_dict['transcript_id']] = self.gene2eggnog[id_dict['refseq_peptide_predicted']]

    def export_eggnog(self):
        with open("eggnog_annot.tsv", 'w') as fo:
            fo.write("transcript_id\teggnog\tsource\n")
            for k, v in self.id_eggnog_dict.items():
                fo.write("\t".join([k, ";".join(v), "protein_id"]) + "\n")

            for k, v in self.blast_id.items():
                if k in self.id_eggnog_dict:
                    pass
                else:
                    eggnog_str = ";".join(self.get_multi_id2eggnog(v))
                    fo.write("\t".join([k, eggnog_str, "blast"]) + "\n")


    def set_out(self):
        outfiles = ["eggnog_annot.tsv"]
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
            'id': 'annot_id2eggnog_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db_medical.annotation.annot_id2eggnog',
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
