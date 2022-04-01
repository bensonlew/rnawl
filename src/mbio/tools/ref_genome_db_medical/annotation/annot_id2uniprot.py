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


class AnnotId2uniprotAgent(Agent):
    """
    to perform Gene Ontology Annotation
    author: liubinxu
    last_modified: 20200811
    """

    def __init__(self, parent):
        super(AnnotId2uniprotAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "ids_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "pir_version", "type": "string", "default": "2019"},
            {"name": "uniprot_version", "type": "string", "default": "202009"},
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
        super(AnnotId2uniprotAgent, self).end()


class AnnotId2uniprotTool(Tool):

    def __init__(self, config):
        super(AnnotId2uniprotTool, self).__init__(config)

        self.uniprot_dict = dict()
        self.ids_uniprot = dict()
        self.blast_id = dict()

        '''
        self.id_mapping_dict = dict({
            "Homo_sapiens": ,
            "Mus_musculus": self.config.SOFTWARE_DIR + "/database/gene_db/mmusculus.all_merged.xls",
            "Rattus_norvegicus": self.config.SOFTWARE_DIR + "/database/gene_db/rnorvegicus.all_merged.xls"
        })
        '''

    def run(self):
        super(AnnotId2uniprotTool, self).run()
        self.convert_xml2table()
        # self.id_des_sqlite =  self.get_uniprot_sqlite()
        self.name2des = self.get_uniprot_file()
        self.id2uniprot()
        self.export_uniprot()
        self.set_out()
        self.end()

    def convert_xml2table(self):
        self.logger.info('转换xml结果为表格格式')
        if os.path.exists("blast_table.xls"):
            pass
        else:
            self.option("blastout").convert2table("blast_table.xls")
        with open('blast_table.xls', 'r') as f:
            for id_dict in csv.DictReader(f, delimiter='\t'):
                if id_dict['Query-Name'] in self.blast_id:
                    pass
                else:
                    uni_id = id_dict['Hit-Name'].split("|")[1]
                    self.blast_id[id_dict['Query-Name']] = [uni_id, id_dict['Hit-Description']]

    def get_uniprot_sqlite(self):
        """
        从参考库UNIPROT_sequence中找到gi号对应的description
        """

        uniprotacc2des = AnnotConfig().get_file_path(
            file="uniprot_acc2des.db",
            db="uniprot",
            db_type="file",
            version=self.option("uniprot_version"))

        conn = sqlite3.connect(uniprotacc2des)
        cursor = conn.cursor()

        return cursor

    def get_uniprot_file(self):
        """
        从参考库UNIPROT_sequence中找到gi号对应的description
        """

        uniprotacc2des = AnnotConfig().get_file_path(
            file="uniprot_acc2des.db",
            db="uniprot",
            db_type="file",
            version=self.option("uniprot_version"))

        uniprotacc2des_file = os.path.dirname(uniprotacc2des) + "/Mammalia.name_des.txt"
        name2des = dict()
        with open(uniprotacc2des_file, 'r') as f:
            for line in f:
                uniprot_name, des = line.strip().split("\t")
                uniname = uniprot_name.split("|")[1]
                name2des[uniname] = des
        return name2des


        # try:
        #     cursor.execute('select * from acc2des')
        #     records = cursor.fetchall()
        #     for record in records:
        #         swiss_id = record[0]
        #         des = record[1]
        #         self.uniprot_dict[swiss_id.split("|")[1]] = [swiss_id, des]

        # except Exception as e:
        #     self.set_error("uniprot 数据库链接错误 {}".format(e))
        # return self.uniprot_dict

    def id2uniprot(self):
        """
        数据库获取信息逻辑
        """
        self.get_uniprot_sqlite()
        with open(self.option("ids_file").prop['path'], 'r') as f:
            for id_dict in csv.DictReader(f, delimiter='\t'):
                if id_dict['transcript_id'] in self.ids_uniprot:
                    pass
                else:
                    if id_dict['uniprot_gn_id'] != "\\N":
                        uniprot_id = id_dict['uniprot_gn_id'].split("_")[0]
                        des = ""
                        if uniprot_id in self.name2des:
                            des = self.name2des[uniprot_id]
                        # try:
                        #     cursor = self.id_des_sqlite
                        #     cursor.execute('select * from acc2des where uniacc="{}"'.format(uniprot_id))
                        #     records = cursor.fetchall()
                        #     record = records[0]
                        #     swiss_id = record[1]
                        #     des = record[2]

                        # except Exception as e:
                        #     self.logger.debug("未找到注释信息 {}".format(e))
                        self.ids_uniprot[id_dict['transcript_id']] = [uniprot_id, des]

    def export_uniprot(self):
        with open("uniprot_annot.tsv", 'w') as fo:
            fo.write("transcript_id\tuniprot\tdescription\tsource\n")
            for k, v in self.ids_uniprot.items():
                try:
                    cursor = self.id_des_sqlite
                    cursor.execute('select * from acc2des where acc="{}"'.format(v))
                    desc = cursor.fetchall()[0][1]
                    description = desc.strip()
                except:
                    description = ""
                fo.write("\t".join([k, v[0], v[1], "uniprot"]) + "\n")

            for k, v in self.blast_id.items():
                if k in self.ids_uniprot:
                    pass
                else:
                    fo.write("\t".join([k, v[0], v[1], "blast"]) + "\n")

    def set_out(self):
        outfiles = ["uniprot_annot.tsv"]
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
        xml_path = '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38.p13_ensembl_98/Annotation_v2/annot_mapdb'
        data = {
            'id': 'annot_id2uniprot_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db_medical.annotation.annot_id2uniprot',
            'instant': False,
            'options': dict(
                ids_file = "/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/ids.tsv",
                blastout = xml_path + '/swissprot/blast.xml',
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
