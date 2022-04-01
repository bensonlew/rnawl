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


class AnnotId2pfamAgent(Agent):
    """
    to perform Gene Ontology Annotation
    author: liubinxu
    last_modified: 20200811
    """

    def __init__(self, parent):
        super(AnnotId2pfamAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "species", "type": "string", "default": None},
            {"name": "ids_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "pfam_domain", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "pfam_version", "type": "string", "default": "33.1"},
        ]
        self.add_option(options)
        self.step.add_steps('pfam_annotation')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        # self.queue = 'BLAST2GO'  # 投递到指定的队列BLAST2GO

    def step_start(self):
        self.step.pfam_annotation.start()
        self.step.update()

    def step_end(self):
        self.step.pfam_annotation.finish()
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
        super(AnnotId2pfamAgent, self).end()


class AnnotId2pfamTool(Tool):

    def __init__(self, config):
        super(AnnotId2pfamTool, self).__init__(config)
        self.pfam_file_dict = dict({
            "Homo_sapiens": self.config.SOFTWARE_DIR + "/database/gene_db/hsapiens.all_merged.xls.domain.tsv",
            "Mus_musculus": self.config.SOFTWARE_DIR + "/database/gene_db/mmusculus.all_merged.xls.domain.tsv",
            "Rattus_norvegicus": self.config.SOFTWARE_DIR + "/database/gene_db/rnorvegicus.all_merged.xls.domain.tsv"
        })
        self.id_pfam_dict = dict()
        self.gene2pfam = dict()
        self.pfam_scan = dict()

    def run(self):
        super(AnnotId2pfamTool, self).run()
        self.get_pfam()
        self.get_pfam_tax_db()
        self.get_pfam_db()
        self.id2pfam()
        self.set_out()
        self.end()

    def get_pfam(self):
        with open(self.option("pfam_domain").prop['path'], 'r') as f:
            for id_dict in csv.DictReader(f, delimiter='\t'):
                if id_dict['Seq_id'] in self.pfam_scan:
                    self.pfam_scan[id_dict['Seq_id']].append([id_dict['Seq_id'], id_dict['Pfam_id'], id_dict['PfamStart'], id_dict['PfamEnd'], '', '', id_dict['Domain'], id_dict['DomainDescription']])
                else:
                    self.pfam_scan[id_dict['Seq_id']] = [[id_dict['Seq_id'], id_dict['Pfam_id'], id_dict['PfamStart'], id_dict['PfamEnd'], '', '', id_dict['Domain'], id_dict['DomainDescription']]]


    def convert_xml2table(self):
        self.logger.info('转换xml结果为表格格式')
        self.option("blastout").convert2table("blast_table.xls")

    def get_pfam_db(self):
        self.pfam_db_dict = dict()
        self.pfam_db = AnnotConfig().get_file_path(
            file = "Pfam-A.clans.tsv",
            db = "pfam",
            version = self.option("pfam_version")
        )
        with open(self.pfam_db, "r") as f:
            for line in f:
                cols = line.strip("\n").split("\t")
                self.pfam_db_dict[cols[0]] = cols[1:]
        # print self.pfam_db_dict

    def get_pfam_tax_db(self):
        """
        从参考库PFAM
        """
        with open(self.pfam_file_dict[self.option('species')], 'r') as f:
            for pfam_dict in csv.DictReader(f, delimiter='\t'):
                if pfam_dict['domain_db'] == 'pfam':
                    tran_id = pfam_dict['ensembl_peptide_id']
                    if tran_id in self.gene2pfam:
                        self.gene2pfam[tran_id].add((pfam_dict['domain_id'], pfam_dict['domain_start'], pfam_dict['domain_end']))
                    else:
                        self.gene2pfam[tran_id] = set([(pfam_dict['domain_id'], pfam_dict['domain_start'], pfam_dict['domain_end'])])
                else:
                    pass


    def id2pfam(self):
        with open("gene2pfam.dict", 'w') as fo:
            fo.write("{}".format(self.gene2pfam))
        with open(self.option("ids_file").prop['path'], 'r') as f,  open("pfam_annot.tsv", 'w') as fo:
            fo.write("Seq_id\tpfam\tdomain_start\tdomain_end\tclan\tclass\tDomain\tDomainDescription\tsource\n")
            for id_dict in csv.DictReader(f, delimiter='\t'):
                if id_dict['protein_id'] != "\\N" and id_dict['protein_id'] in self.gene2pfam:
                    for pfam in self.gene2pfam[id_dict['protein_id']]:
                        if pfam[0] != "\\N":
                            fo.write("\t".join([id_dict['transcript_id']] + list(pfam) + self.pfam_db_dict.get(pfam[0], ["", "", "", ""]) ) + "\tensembl\n")
                elif 'ensembl_peptide_id' in id_dict and id_dict['ensembl_peptide_id'] != "\\N" and id_dict['ensembl_peptide_id'] in self.gene2pfam:
                    for pfam in self.gene2pfam['ensembl_transcript_id']:
                        if pfam[0] != "\\N":
                            fo.write("\t".join([id_dict['transcript_id']] + list(pfam) + self.pfam_db_dict.get(pfam[0], ["", "", "", ""])) + "\tensembl\n")

            for k, v in self.pfam_scan.items():
                if k in self.gene2pfam:
                    pass
                else:
                    for pfam_range in v:
                        fo.write("{}\thmm_scan\n".format("\t".join(pfam_range)))

    def set_out(self):
        outfiles = ["pfam_annot.tsv"]
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
            'id': 'annot_id2pfam_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db_medical.annotation.annot_id2pfam',
            'instant': False,
            'options': dict(
                ids_file = "/mnt/ilustre/users/sanger-dev/workspace/20200811/Single_annot_getid_5420_4347/AnnotGetid/ids.tsv",
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
