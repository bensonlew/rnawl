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


class AnnotId2reactomeAgent(Agent):
    """
    to perform Gene Ontology Annotation
    author: liubinxu
    last_modified: 20200811
    """

    def __init__(self, parent):
        super(AnnotId2reactomeAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "ids_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "species", "type": "string", "default": None},
            {"name": "known_reactome", "type": "string", "default": None},
            {"name": "reactome_version", "type": "string", "default": "72"},
        ]
        self.add_option(options)
        self.step.add_steps('reactome_annotation')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        # self.queue = 'BLAST2GO'  # 投递到指定的队列BLAST2GO

    def step_start(self):
        self.step.reactome_annotation.start()
        self.step.update()

    def step_end(self):
        self.step.reactome_annotation.finish()
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
        super(AnnotId2reactomeAgent, self).end()


class AnnotId2reactomeTool(Tool):
    def __init__(self, config):
        super(AnnotId2reactomeTool, self).__init__(config)

        # self.reactome_tax = AnnotConfig().get_file_path(
        #     file ="UniProt2Reactome.txt",
        #     db = "reactome",
        #     version = self.option("reactome_version"))

        self.reactome_from = AnnotConfig().get_file_path(
            file ="PhysicalEntity_2_inferredFrom.tsv",
            db = "reactome",
            version = self.option("reactome_version"))

        self.gene2path = AnnotConfig().get_file_path(
            file ="UniProt2Reactome_PE_All_Levels.txt",
            db = "reactome",
            version = "72")
        self.reactome_entity2acc = AnnotConfig().get_file_path(
            file ="EntityWithAccessionedSequence.tsv",
            db = "reactome",
            version = self.option("reactome_version"))

        self.reactome_id2refentity = AnnotConfig().get_file_path(
            file ="ReferenceEntity.tsv",
            db = "reactome",
            version = self.option("reactome_version"))

        self.abr_dict = dict({
            "Homo_sapiens": "R-HSA-",
            "Mus_musculus": "R-MMU-",
            "Rattus_norvegicus": "R-RNO-"
        })
        self.id_reactome_dict = dict()
        self.uni2reactome = dict()

    def run(self):
        super(AnnotId2reactomeTool, self).run()
        # self.convert_xml2table()
        self.get_reactome_gene_list()
        self.get_reactome_from()
        self.get_reactome_db()
        self.id2reactome()
        self.export_reactome()
        self.set_out()
        self.end()

    def get_reactome_gene_list(self):
        """
        用于筛选不在数据库的id映射， 物种间可能映射错误
        """
        self.reactome_gene_list = list()
        with open(self.gene2path, 'r') as f:
            for line in f:
                cols = line.strip().split("\t")
                if cols[1] in self.reactome_gene_list:
                    pass
                else:
                    self.reactome_gene_list.append(cols[1])

    def get_reactome_db(self):
        """
        从参考库REACTOME_sequence中找到gi号对应的description
        """

        self.entity2id = dict()
        with open(self.reactome_entity2acc, 'r') as f:
            for id_dict in csv.DictReader(f, delimiter='\t'):
                if id_dict['referenceEntity'] in self.entity2id:
                    self.entity2id[id_dict['referenceEntity']].add(id_dict['DB_ID'])
                else:
                    self.entity2id[id_dict['referenceEntity']] = set([id_dict['DB_ID']])

        self.uni2reactome = dict()
        with open(self.reactome_id2refentity, 'r') as f:
            for id2_dict in csv.DictReader(f, delimiter='\t'):
                uniprot_id = id2_dict['identifier'].strip("'")
                entity_id = id2_dict['DB_ID'].strip("'")
                if entity_id in self.entity2id:
                    acc_ids = self.entity2id[entity_id]
                    # 大鼠小鼠对应方式不同
                    if self.option("species") != "Homo_sapiens":
                        acc_ids = set(self.reactome_from_dict[x] for x in acc_ids if x in self.reactome_from_dict)
                else:
                    acc_ids = set()
                if len(acc_ids) >= 1:
                    if uniprot_id in self.uni2reactome:
                        self.uni2reactome[uniprot_id].extend(list(acc_ids))
                    else:
                        self.uni2reactome[uniprot_id] = list(acc_ids)
        '''
        with open(self.reactome_tax, 'r') as f:
            f.readline()
            for line in f:
                cols = line.strip().split("\t")
                uniprot_id = cols[1].strip("'")
                if cols[2] != "-":
                    if uniprot_id in self.uni2reactome:
                        self.uni2reactome[uniprot_id].add(self.abr_dict[self.option("species")] + cols[2])
                    else:
                        self.uni2reactome[uniprot_id] = set([self.abr_dict[self.option("species")] + cols[2]])
        '''
    def get_reactome_from(self):
        self.reactome_from_dict = dict()
        with open(self.reactome_from, 'r') as f:
            for id2_dict in csv.DictReader(f, delimiter='\t'):
                db_id = id2_dict['DB_ID'].strip("'")
                from_id = id2_dict['inferredFrom'].strip("'")
                self.reactome_from_dict[db_id] = from_id

    def id2reactome(self):
        with open(self.option("ids_file").prop['path'], 'r') as f:
            for id_dict in csv.DictReader(f, delimiter='\t'):
                '''
                biomart 为通路编号
                if id_dict['reactome'].startswith("R-"):
                    self.id_reactome_dict[id_dict['transcript_id']] = set(id_dict['reactome'].split("|"))
                else:
                '''
                if id_dict['uniprot_gn_id'] in self.uni2reactome:
                    reac_ids = [self.abr_dict[self.option("species")] + acc
                                for acc in set(self.uni2reactome[id_dict['uniprot_gn_id']])
                                if self.abr_dict[self.option("species")] + acc in self.reactome_gene_list
                    ]
                    self.id_reactome_dict[id_dict['transcript_id']] = reac_ids
                else:
                    pass

    def export_reactome(self):
        with open("reactome_annot.tsv", 'w') as fo:
            fo.write("transcript_id\treactome\tdescription\n")
            for k, v in self.id_reactome_dict.items():
                fo.write("\t".join([k, ";".join(v)]) + "\n")

    def set_out(self):
        outfiles = ["reactome_annot.tsv"]
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
            'id': 'annot_id2reactome_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db_medical.annotation.annot_id2reactome',
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
