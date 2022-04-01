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


class AnnotId2swissprotAgent(Agent):
    """
    to perform Gene Ontology Annotation
    author: liubinxu
    last_modified: 20200811
    """

    def __init__(self, parent):
        super(AnnotId2swissprotAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "blastout_go", "type": "outfile", "format": "ref_rna_v2.blast_xml"},
            {"name": "blast2go_annot", "type": "outfile", "format": "ref_rna_v2.blast2go_annot"},
            {"name": "ids_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "known_go", "type": "string", "default": None},
            {"name": "pir_version", "type": "string", "default": "2019"},
            {"name": "swissprot_version", "type": "string", "default": "202006"},
            {"name": "species", "type": "string", "default": ""},
            {"name": "p", "type": "int", "default": 10},
        ]
        self.add_option(options)
        self.on('start', self.step_start)
        self.on('end', self.step_end)


    def step_start(self):
        self.step.update()

    def step_end(self):
        self.step.update()

    def check_options(self):
        if self.option("ids_file").is_set:
            pass
        else:
            raise OptionError("必须提供id结果文件")

    def set_resource(self):
        self._cpu = '2'
        self._memory = '20G'

    def end(self):
        super(AnnotId2swissprotAgent, self).end()


class AnnotId2swissprotTool(Tool):

    def __init__(self, config):
        super(AnnotId2swissprotTool, self).__init__(config)

        self.idmapping_db = self.idmapping_db = AnnotConfig().get_file_path(
            file ="idmapping.tb",
            db = "pir",
            version = self.option("pir_version"))
        # self.id2nr = self.config.PACKAGE_DIR + "/ref_genome_db_medical//id2nr.py"
        self.swissprot_dict = dict()
        self.ids_swissprot = dict()
        self.blast_id = dict()
        '''
        self.id_mapping_dict = dict({
            "Homo_sapiens": ,
            "Mus_musculus": self.config.SOFTWARE_DIR + "/database/gene_db/mmusculus.all_merged.xls",
            "Rattus_norvegicus": self.config.SOFTWARE_DIR + "/database/gene_db/rnorvegicus.all_merged.xls"
        })
        '''

    def run(self):
        super(AnnotId2swissprotTool, self).run()
        self.convert_xml2table()
        self.id2swissprot()
        self.export_swissprot()
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

    def get_swissprot_sqlite(self):
        """
        从参考库SWISSPROT_sequence中找到gi号对应的description
        """

        swissprotacc2des = AnnotConfig().get_file_path(
            file="swissprot_acc2des.db",
            db="swissprot",
            db_type="file",
            version=self.option("swissprot_version"))

        conn = sqlite3.connect(swissprotacc2des)
        cursor = conn.cursor()

        try:
            cursor.execute('select * from acc2des')
            records = cursor.fetchall()
            for record in records:
                swiss_id = record[0]
                des = record[1]
                self.swissprot_dict[swiss_id.split("|")[1]] = [swiss_id, des]

        except Exception as e:
            self.set_error("swissprot 数据库链接错误 {}".format(e))
        return self.swissprot_dict

    def id2swissprot(self):
        """
        数据库获取信息逻辑
        """
        self.get_swissprot_sqlite()
        with open(self.option("ids_file").prop['path'], 'r') as f:
            for id_dict in csv.DictReader(f, delimiter='\t'):
                if id_dict['transcript_id'] in self.ids_swissprot:
                    pass
                else:
                    if id_dict['uniprotswissprot'] != "\\N":
                        self.ids_swissprot[id_dict['transcript_id']] = self.swissprot_dict.get(id_dict['uniprotswissprot'],
                                                                                           [id_dict['uniprotswissprot'], ""])

    def export_swissprot(self):
        with open("swissprot_annot.tsv", 'w') as fo:
            fo.write("transcript_id\tswissprot\tdescription\tsource\n")
            for k, v in self.ids_swissprot.items():
                try:
                    cursor.execute('select * from acc2des where acc="{}"'.format(v))
                    desc = cursor.fetchall()[0][1]
                    description = desc.strip()
                except:
                    description = ""
                fo.write("\t".join([k, v[0], v[1], "uniprot"]) + "\n")

            for k, v in self.blast_id.items():
                if k in self.ids_swissprot:
                    pass
                else:
                    fo.write("\t".join([k, v[0], v[1], "blast"]) + "\n")

    def set_out(self):
        outfiles = ["swissprot_annot.tsv"]
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
            'id': 'annot_id2swissprot_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db_medical.annotation.annot_id2swissprot',
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
