# -*- coding: utf-8 -*-
# __author__ = 'liubinxu1'
# last_modify:2018.05.08

from biocluster.module import Module
import os
import shutil
import re
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import unittest
from mbio.packages.rna.annot_config import AnnotConfig
from mbio.packages.ref_rna_v2.copy_file import CopyFile

class AnnotMergeidModule(Module):
    '''
    将序列与数据库比对 包括blast\pfam数据库和GO、COG、KEGG对应关系映射
    '''
    def __init__(self, work_id):
        super(AnnotMergeidModule, self).__init__(work_id)
        self._fasta_type = {'Protein': 'prot', 'DNA': 'nucl'}
        options = [
            {"name": "database", "type": "string", "default": "nr;kegg;uniprot;eggnog;go;pfam;do;reactome;disgenet"},

            {"name": "nr_db", "type": "string", "default": "nr"},
            {"name": "known_go", "type": "string", "default": None},
            {"name": "tax", "type": "bool", "default": False}, #要不要做物种分类
            {"name": "string_db", "type": "string", "default": "string"},
            {"name": "swissprot_db", "type": "string", "default": "swissprot"},
            {"name": "eggnog_db", "type": "string", "default": "eggnog"},
            {"name": "kegg_db", "type": "string", "default": "kegg"},
            {"name": "kegg_version", "type": "string", "default": ""},

            {"name": "species", "type": "string", "default": ""}, # 物种
            # 比对数据库 nt nr string swissprot kegg customer_mode
            {"name": "merge_type", "type": "string", "default": "partial"},
            {"name": "go_version", "type": "string", "default": ""},
            {"name": "pfam_version", "type": "string", "default": "32"},
            {'name': 'kegg_version', 'type': 'string', 'default': "2019"},
            {"name": "nr_version", "type": "string", "default": "2019"},
            {"name": "swissprot_version", "type": "string", "default": "2019"},
            {"name": "blast_nr_xml", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "blast_string_xml", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "blast_eggnog_xml", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "blast_kegg_xml", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "blast_swissprot_xml", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "blast_uniprot_xml", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "blast2go_annot", "type": "infile", "format": "ref_rna_v2.blast2go_annot"},
            {"name": "pfam_domain", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "string_version", "type": "string", "default": "2019"},
            {"name": "pir_version", "type": "string", "default": "2019"},
            {"name": "eggnog_version", "type": "string", "default": "2019"},
            {"name": "version", "type": "string", "default": "2019"},

            {"name": "ids_file", "type": "outfile", "format": "ref_rna_v2.common"},

            # 当输出格式为非5，6时，只产生文件不作为outfile
        ]
        self.add_option(options)
        self.annot_nr = self.add_tool("ref_genome_db_medical.annotation.annot_id2nr")
        self.annot_eggnog = self.add_tool("ref_genome_db_medical.annotation.annot_id2eggnog")
        self.annot_go = self.add_tool("ref_genome_db_medical.annotation.annot_id2go")
        self.annot_kegg = self.add_tool("ref_genome_db_medical.annotation.annot_id2kegg")
        self.annot_pfam = self.add_tool("ref_genome_db_medical.annotation.annot_id2pfam")
        # self.annot_swissprot = self.add_tool("ref_genome_db_medical.annotation.annot_id2swissprot")
        self.annot_uniprot = self.add_tool("ref_genome_db_medical.annotation.annot_id2uniprot")
        self.annot_do = self.add_tool("ref_genome_db_medical.annotation.annot_id2do")
        self.annot_reactome = self.add_tool("ref_genome_db_medical.annotation.annot_id2reactome")
        self.annot_disgenet = self.add_tool("ref_genome_db_medical.annotation.annot_id2disgenet")
        self.step.add_steps('annot_nr', 'annot_eggnog', 'annot_go', 'annot_pfam', 'annot_kegg', 'annot_uniprot', 'annot_do', 'annot_reactome', 'annot_disgenet')

        '''
        # blast/diamond 分布运行结果
        self.blast_nr_tools = []
        self.id2go = []
        self.blast_uniprot_tools = []
        self.blast_kegg_tools = []
        self.blast_string_tools = []
        self.blast_eggnog_tools = []

        self.hmm_pfam_tools = []

        self.catblast_tools = []
        self.string2cog_tools = []
        self.kegg2ko_tools = []
        self.nr2go_tools = []
        self.nr2ncbitax_tools = []
        self.merge_tools = {}

        self.blast_opts = {
            "query_type": self.option("query_type"),
            "outfmt": 5,
            "blast": self.option("blast"),
            "evalue": self.option("evalue"),
            "num_threads": self.option("num_threads"),
        }
        '''

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def check_options(self):
        if not self.option("ids_file").is_set:
            raise OptionError("必须设置参数ids_file")
        return True

    def run_splitfasta(self):
        '''
        fasta文件分割
        '''
        self.splitfasta.set_options({
            "fasta": self.option("query"),
            "lines": self.option("lines"),
        })
        self.splitfasta.on('start', self.set_step, {'start': self.step.split_fasta})
        self.splitfasta.on('end', self.set_step, {'end': self.step.split_fasta})
        self.splitfasta.on('end', self.set_step, {'start': self.step.blast_nr})
        self.splitfasta.on('end', self.set_step, {'start': self.step.blast_kegg})
        self.splitfasta.on('end', self.set_step, {'start': self.step.blast_string})
        self.splitfasta.on('end', self.set_step, {'start': self.step.blast_eggnog})
        self.splitfasta.on('end', self.set_step, {'start': self.step.blast_uniprot})
        self.splitfasta.on('end', self.run_blast)
        self.splitfasta.run()


    def run_id2nr(self):
        '''
        获取NR注释
        '''
        options = {
            'ids_file': self.option("ids_file"),
            'blastout': self.option("blast_nr_xml").prop['path'],
        }

        self.annot_nr.on('start', self.set_step, {'start': self.step.annot_nr})
        self.annot_nr.on('end', self.set_step, {'end': self.step.annot_nr})
        self.annot_nr.set_options(options)
        self.annot_nr.run()


    def run_id2uniprot(self):
        '''
        获取uniprot注释
        '''
        options = {
            'ids_file': self.option("ids_file"),
            'blastout': self.option("blast_uniprot_xml").prop['path'],
        }
        self.annot_uniprot.on('start', self.set_step, {'start': self.step.annot_uniprot})
        self.annot_uniprot.on('end', self.set_step, {'end': self.step.annot_uniprot})
        self.annot_uniprot.set_options(options)
        self.annot_uniprot.run()


    def run_id2kegg(self):
        '''
        获取kegg注释
        '''
        options = {
            'ids_file': self.option("ids_file"),
            'blastout': self.option("blast_kegg_xml").prop['path'],
            'species': self.option("species"),
        }
        self.annot_kegg.on('start', self.set_step, {'start': self.step.annot_kegg})
        self.annot_kegg.on('end', self.set_step, {'end': self.step.annot_kegg})
        self.annot_kegg.set_options(options)
        self.annot_kegg.run()


    def run_id2eggnog(self):
        '''
        获取eggnog注释
        '''
        options = {
            'ids_file': self.option("ids_file"),
            'blastout': self.option("blast_eggnog_xml").prop['path'],
        }

        self.annot_eggnog.on('start', self.set_step, {'start': self.step.annot_eggnog})
        self.annot_eggnog.on('end', self.set_step, {'end': self.step.annot_eggnog})
        self.annot_eggnog.set_options(options)
        self.annot_eggnog.run()

    def run_id2reactome(self):
        '''
        获取reactome注释
        '''
        options = {
            'ids_file': self.option("ids_file"),
            'species': self.option("species"),
        }

        self.annot_reactome.on('start', self.set_step, {'start': self.step.annot_reactome})
        self.annot_reactome.on('end', self.set_step, {'end': self.step.annot_reactome})
        self.annot_reactome.set_options(options)
        self.annot_reactome.run()


    def run_id2do(self):
        '''
        获取do注释
        '''
        options = {
            'ids_file': self.option("ids_file")
        }

        self.annot_do.on('start', self.set_step, {'start': self.step.annot_do})
        self.annot_do.on('end', self.set_step, {'end': self.step.annot_do})
        self.annot_do.set_options(options)
        self.annot_do.run()

    def run_id2pfam(self):
        '''
        获取pfam注释
        '''
        options = {
            'ids_file': self.option("ids_file"),
            'species': self.option("species"),
            'pfam_domain': self.option("pfam_domain").prop['path'],
        }

        self.annot_pfam.on('start', self.set_step, {'start': self.step.annot_pfam})
        self.annot_pfam.on('end', self.set_step, {'end': self.step.annot_pfam})
        self.annot_pfam.set_options(options)
        self.annot_pfam.run()


    def run_id2go(self):
        '''
        获取go注释
        '''
        options = {
            'ids_file': self.option("ids_file"),
            'blast2go_annot': self.option('blast2go_annot')
        }

        self.annot_go.on('start', self.set_step, {'start': self.step.annot_go})
        self.annot_go.on('end', self.set_step, {'end': self.step.annot_go})
        self.annot_go.set_options(options)
        self.annot_go.run()


    def run_id2disgenet(self):
        '''
        获取go注释
        '''
        options = {
            'ids_file': self.option("ids_file")
        }
        self.annot_disgenet.on('start', self.set_step, {'start': self.step.annot_disgenet})
        self.annot_disgenet.on('end', self.set_step, {'end': self.step.annot_disgenet})
        self.annot_disgenet.set_options(options)
        self.annot_disgenet.run()

    '''
    def run_blast_nr(self):
        opts = self.blast_opts.copy()
        opts.update({"database": self.option("nr_db")})
        opts.update({"version": self.option("version")})
        opts.update({"nr_version": self.option("nr_version")})
        i = 0
        for f in os.listdir(self.splitfasta.output_dir):
            opts['query'] = os.path.join(self.splitfasta.output_dir, f)
            if self.option('method') == "blast":
                blast_tool = self.add_tool('ref_genome_db_v2.annotation.blast')
            else:
                blast_tool = self.add_tool('ref_genome_db_v2.annotation.diamond')
        if len(self.blast_nr_tools) == 1:
            self.blast_nr_tools[0].on("end", self.run_id2nr, "nr")
        else:
            self.on_rely(self.blast_nr_tools, self.run_id2nr, "nr")
        for tool in self.blast_nr_tools:
            tool.run()

    '''

    def set_output(self):
        '''
        设置输出结果
        '''
        self.logger.info("event name is {}".format(self._rely.keys()))

        for db in self.option("database").split(";"):
            tool = getattr(self, "annot_" + db)
            CopyFile().linkdir(tool.output_dir, self.output_dir + "/" + db)

        self.end()

    def run(self):
        # 设置运行逻辑
        super(AnnotMergeidModule, self).run()
        self.id_tools = list()
        for db in self.option("database").split(";"):
            if db == 'nr':
                self.id_tools.append(self.annot_nr)
            elif db == 'uniprot':
                self.id_tools.append(self.annot_uniprot)
            elif db == 'kegg':
                self.id_tools.append(self.annot_kegg)
            elif db == 'eggnog':
                self.id_tools.append(self.annot_eggnog)
            elif db == 'go':
                self.id_tools.append(self.annot_go)
            elif db == 'pfam':
                self.id_tools.append(self.annot_pfam)
            elif db == 'reactome':
                self.id_tools.append(self.annot_reactome)
            elif db == 'do':
                self.id_tools.append(self.annot_do)
            elif db == 'disgenet':
                self.id_tools.append(self.annot_disgenet)
            else:
                self.logger.info("不支持该类型{}的数据库".format(db))

        self.on_rely(self.id_tools, self.set_output)

        for db in self.option("database").split(";"):
            if db == 'nr':
                self.run_id2nr()
            elif db == 'uniprot':
                self.run_id2uniprot()
            elif db == 'kegg':
                self.run_id2kegg()
            elif db == 'eggnog':
                self.run_id2eggnog()
            elif db == 'go':
                self.run_id2go()
            elif db == 'pfam':
                self.run_id2pfam()
            elif db == 'reactome':
                self.run_id2reactome()
            elif db == 'do':
                self.run_id2do()
            elif db == 'disgenet':
                self.run_id2disgenet()
            else:
                self.logger.info("不支持该类型{}的数据库".format(db))


    def end(self):
        repaths = [
            [".", "", "blast输出目录"],
            ["blast.xml", "xml", "blast xml输出结果文件"],
            ["blast_table.xls", "xls", "blast xls输出结果文件"]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        super(AnnotMergeidModule, self).end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime

        test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_ref_genome_db_v2'

        xml_path = '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38.p13_ensembl_98/Annotation_v2/annot_mapdb'
        data = {
            "id": "annot_mapdb_rerun",
            "type": "module",
            "rerun": True,
            "name": "ref_genome_db_medical.annot_mapdb",
            "skip_all_success": True,
            "instant": False,
            "options": dict(
                {
                    'ids_file' : "/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/ids.tsv",
                    'species' : "Homo_sapiens",
                    'blast_nr_xml': xml_path + '/nr/blast.xml',
                    'blast_uniprot_xml': xml_path + '/swissprot/blast.xml',
                    'blast_eggnog_xml': xml_path + '/eggnog/blast.xml',
                    'blast_kegg_xml': xml_path + '/kegg/blast.xml',

                    'blast2go_annot': xml_path + '/GO/blast2go_merge.xls',
                    'pfam_domain': xml_path + '/../annot_orfpfam/pfam_domain',
                }
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.is_skip=True
        wf.run()

if __name__ == '__main__':
    unittest.main()
