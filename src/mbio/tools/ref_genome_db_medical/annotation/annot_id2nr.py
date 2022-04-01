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


class AnnotId2nrAgent(Agent):
    """
    to perform Gene Ontology Annotation
    author: liubinxu
    last_modified: 20200811
    """

    def __init__(self, parent):
        super(AnnotId2nrAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "blastout_go", "type": "outfile", "format": "ref_rna_v2.blast_xml"},
            {"name": "blast2go_annot", "type": "outfile", "format": "ref_rna_v2.blast2go_annot"},
            {"name": "ids_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "known_go", "type": "string", "default": None},
            {"name": "pir_version", "type": "string", "default": "2019"},
            {"name": "nr_version", "type": "string", "default": "202006"},
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
        super(AnnotId2nrAgent, self).end()


class AnnotId2nrTool(Tool):

    def __init__(self, config):
        super(AnnotId2nrTool, self).__init__(config)

        self.idmapping_db = self.idmapping_db = AnnotConfig().get_file_path(
            file ="idmapping.tb",
            db = "pir",
            version = self.option("pir_version"))
        # self.id2nr = self.config.PACKAGE_DIR + "/ref_genome_db_medical//id2nr.py"
        self.ids_refseqpep = dict()
        self.blast_id = dict()

    def run(self):
        super(AnnotId2nrTool, self).run()
        # self.convert_xml2table()
        self.convert_xml2table()
        self.id2nr()
        self.export_nr()
        self.set_out()
        self.end()

    def convert_xml2table(self):
        self.logger.info('转换xml结果为表格格式')
        self.option("blastout").convert2table("blast_table.xls")
        with open('blast_table.xls', 'r') as f:
            for id_dict in csv.DictReader(f, delimiter='\t'):
                if id_dict['Query-Name'] in self.blast_id:
                    pass
                else:
                    self.blast_id[id_dict['Query-Name']] = [id_dict['Hit-Name'], id_dict['Hit-Description']]

    def get_nr_sqlite(self):
        """
        从参考库NR_sequence中找到gi号对应的description
        """

        nracc2des = AnnotConfig().get_file_path(
            file="nr_acc2des.db",
            db="nr",
            db_type="file",
            version=self.option("nr_version"))

        conn = sqlite3.connect(nracc2des)
        cursor = conn.cursor()
        return cursor

    def get_nr_des(self):
        nracc2des = AnnotConfig().get_file_path(
            file="medical_acc2des.txt",
            db="nr",
            db_type="file",
            version=self.option("nr_version"))

        nr_acc2des_dict = dict()
        nr_acc2des_dict2 = dict()
        with open(nracc2des, 'r') as f:
            for line in f:
                acc = line.strip().split()[0]
                des = " ".join(line.strip().split()[1:])
                nr_acc2des_dict[acc] = des
                if "." in acc:
                    acc1 = acc.split(".")[0]
                    nr_acc2des_dict2[acc1] = [acc, des]

        return nr_acc2des_dict, nr_acc2des_dict2

    def id2nr(self):
        with open(self.option("ids_file").prop['path'], 'r') as f:
            for id_dict in csv.DictReader(f, delimiter='\t'):
                if id_dict['transcript_id'] in self.ids_refseqpep:
                    pass
                else:
                    if id_dict['refseq_peptide'] != "\\N" and id_dict['refseq_peptide'] != "":
                        self.ids_refseqpep[id_dict['transcript_id']] = id_dict['refseq_peptide'].split("|")[0]
                    elif id_dict['refseq_peptide_predicted'] != "\\N" and id_dict['refseq_peptide_predicted'] != "":
                        self.ids_refseqpep[id_dict['transcript_id']] = id_dict['refseq_peptide_predicted'].split("|")[0]

    def export_nr(self):
        cursor = self.get_nr_sqlite()
        nr_acc2des_dict, nr_acc2des_dict2 = self.get_nr_des()
        with open("nr_annot.tsv", 'w') as fo:
            fo.write("transcript_id\tnr\tdescription\tsource\n")
            for k, v in self.ids_refseqpep.items():
                if v in nr_acc2des_dict:
                    description = nr_acc2des_dict.get(v, "")
                else:
                    description = nr_acc2des_dict2.get(v, ["", ""])[1]

                fo.write("\t".join([k, v, description, "refseq"]) + "\n")

            for k, v in self.blast_id.items():
                if k in self.ids_refseqpep:
                    pass
                else:
                    fo.write("\t".join([k, v[0], v[1], "blast"]) + "\n")

    def set_out(self):
        outfiles = ["nr_annot.tsv"]
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
            'id': 'annot_id2nr_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db_medical.annotation.annot_id2nr',
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
