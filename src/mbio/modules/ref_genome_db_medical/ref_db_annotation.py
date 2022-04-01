# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modifiy = modified 2019.01.08

from biocluster.module import Module
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from biocluster.wsheet import Sheet
import os
import re
import glob
import json
from biocluster.config import Config
from mbio.packages.lnc_rna.copy_file import CopyFile

class RefDbAnnotationModule(Module):
    """
    交互分析进行筛选blast参数重注释时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        #self.rpc = False
        super(RefDbAnnotationModule, self).__init__(wsheet_object)
        options = [
            {"name": "nr_evalue", "type": "float", "default": 1e-3},
            {"name": "ref_fa", "type": "infile", "format": "ref_rna_v2.fasta"},
            {"name": "swissprot_evalue", "type": "float", "default": 1e-3},
            {"name": "cog_evalue", "type": "float", "default": 1e-3},
            {"name": "kegg_evalue", "type": "float", "default": 1e-3},
            {"name": "pfam_evalue", "type": "float", "default": 1e-3},
            {"name": "nr_filter_evalue", "type": "float", "default": 1e-5},
            {"name": "swissprot_filter_evalue", "type": "float", "default": 1e-5},
            {"name": "cog_filter_evalue", "type": "float", "default": 1e-5},
            {"name": "kegg_filter_evalue", "type": "float", "default": 1e-5},
            {"name": "pfam_filter_evalue", "type": "float", "default": 1e-5},
            {"name": "blast_method", "type": "string", "default": "diamond"},
            {"name": "has_new", "type": "bool", "default": True},
            {"name": "species_name", "type": "string"},
            {"name": "origin_result", "type": "infile", "format": "ref_rna_v2.common_dir"},
            {"name": "merge_type", "type": "string", "default": "partial"},
            # all 全部合并， partial 只合并未注释基因， refonly 只保留已知注释
            {"name": "trans2gene", "type": "string"},
            {"name": "gtf", "type": "string"},
            {"name": "g2t2p", "type": "string"},
            {"name": "biomart", "type": "string"},
            {"name": "biomart_type", "type": "string"},
            {"name": "enterz", "type": "string"},
            {"name": "taxonomy", "type": "string", "default": None},
            {"name": "pep", "type": "string", "default": ""},
            {"name": "trans", "type": "string", "default": ""},
            {"name": "species_class", "type": "string", "default": "Animals"},
            {"name": "species", "type": "string", "default": None},
            {"name": "known_go", "type": "string", "default": ""},
            {"name": "known_ko", "type": "string", "default": ""},

            {"name": "nr_version", "type": "string", "default": "2019"},
            {"name": "swissprot_version", "type": "string", "default": "2019"},
            {"name": "uniprot_version", "type": "string", "default": "202009"},
            {"name": "eggnog_version", "type": "string", "default": "2019"},
            {"name": "cog_version", "type": "string", "default": "2019"},
            {"name": "string_version", "type": "string", "default": "2019"},
            {"name": "pir_version", "type": "string", "default": "2019"},
            {"name": "kegg_version", "type": "string", "default": "202007"},
            {"name": "go_version", "type": "string", "default": "2019"},
            {"name": "pfam_version", "type": "string", "default": "32"},
            {"name": "genedb_version", "type": "string", "default": ""},
        ]
        self.add_option(options)

        self.orf_pfam = self.add_module("ref_genome_db_medical.annot_orfpfam")
        self.get_ids = self.add_tool("ref_genome_db_medical.annotation.annot_getid")
        self.map_db = self.add_module("ref_genome_db_medical.annot_mapdb")
        self.merge_id = self.add_module("ref_genome_db_medical.annot_mergeid")
        # self.ref_filter = self.add_module("ref_genome_db_medical.annot_filter")
        self.ref_class_medical = self.add_module("ref_genome_db_medical.annot_class_beta")
        self.ref_class = self.add_module("ref_genome_db_v2.annot_class_beta")
        self.step.add_steps("map_db", "merge_id", "get_ids", "orf_pfam", "ref_filter","ref_class","ref_class_medical")

        # self.json_path = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/annot_species.v2.json"
        # self.genome_db_api = self.api.api("ref_rna_v2.genome_db")
        self.species_list = ["Homo_sapiens", "Mus_musculus", "Rattus_norvegicus"]
        self.species2abr = {"Homo_sapiens": "hsa", "Mus_musculus": "mmu", "Rattus_norvegicus": "rno"}

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def check_options(self):
        if self.option("species_name") not in ["Homo_sapiens", "Mus_musculus", "Rattus_norvegicus"]:
            self.set_error("物种必须为 {} 中的一种".format(self.species_list))


    def get_db_annot(self):
        pass
        '''
        获取数据库注释路径信息
        ; 返回参考注释路径
        '''
        '''
        f = open(self.json_path, "r")
        json_dict = json.loads(f.read())
        genome_dict = self.genome_db_api.get_genome_dict_by_task_id(self.option("task_id"))
        if genome_dict:
            return dict(genome_dict)
        elif json_dict.has_key(self.option("species_name")):
            return json_dict[self.option("species_name")]
        else:
            self.logger.info("数据库中没有该物种{}".format(self.option("species_name")))
        '''

    def run_orf_pfam(self):
        '''
        使用蛋白序列比对pfam数据库
        '''
        options = dict(
            pep = self.option("pep"),
            lines = 5000,
            gtf = self.option("gtf"),
            g2t2p = self.option("g2t2p"),
            pfam_version = self.option("pfam_version")
        )
        
        self.orf_pfam.on('start', self.set_step, {'start': self.step.orf_pfam})
        self.orf_pfam.on('end', self.set_step, {'end': self.step.orf_pfam})
        self.orf_pfam.set_options(options)
        self.orf_pfam.run()


    def run_get_ids(self):
        '''
        获取idmapping数据
        '''
        getline_list = [
            "ensembl_gene_id",
            "external_gene_name",
            "ensembl_transcript_id",
            "external_transcript_name",
            "gene_biotype",
            "transcript_biotype",
            "uniprot_gn_id",
            "uniprotswissprot",
            "reactome",
            "kegg",
            "eggnog",
            "entrezgene_id",
            "refseq_mrna",
            "refseq_peptide",
            "refseq_ncrna",
            "refseq_peptide_predicted",
            "refseq_ncrna_predicted",
            "refseq_mrna_predicted"
        ]
        matchline_list = [
            "ensembl_transcript_id",
            "refseq_mrna",
            "refseq_ncrna",
            "refseq_ncrna_predicted",
            "refseq_mrna_predicted"
        ]
        options = dict(
            g2t2p = self.option("g2t2p"),
            species = self.option("species_name"),
            get_lines = ",".join(getline_list),
            match_lines = ",".join(matchline_list)
        )

        self.get_ids.on('start', self.set_step, {'start': self.step.get_ids})
        self.get_ids.on('end', self.set_step, {'end': self.step.get_ids})
        self.get_ids.set_options(options)
        self.get_ids.run()


    def run_map_db(self):
        '''
        使用蛋白序列比对pfam数据库
        '''
        db2nr = {
            "protists": "protist",
            "vertebrates": "metazoa",
            "metazoa": "metazoa",
            "fungi": "fungi",
            "plants": "viridiplantae",
        }
        options = dict(
            query = self.option("trans"),
            method = self.option("blast_method"),
            nr_db = db2nr[self.option("species_class")],
            known_go = self.option("known_go"),
            merge_type = self.option("merge_type"),
            nr_version = self.option("nr_version"),
            kegg_version = self.option("kegg_version"),
            eggnog_version = self.option("eggnog_version"),
            swissprot_version = self.option("swissprot_version"),
            uniprot_version = self.option("uniprot_version"),
            pir_version = self.option("pir_version"),
            lines = 100000
        )

        self.map_db.on('start', self.set_step, {'start': self.step.map_db})
        self.map_db.on('end', self.set_step, {'end': self.step.map_db})
        self.map_db.set_options(options)
        self.map_db.run()

    def run_merge_id(self):
        '''
        使用蛋白序列比对pfam数据库
        '''
        db2nr = {
            "protists": "protist",
            "vertebrates": "metazoa",
            "metazoa": "metazoa",
            "fungi": "fungi",
            "plants": "viridiplantae",
        }
        xml_path = self.map_db.output_dir
        pfam_dir = self.orf_pfam.output_dir
        options = {
            'ids_file' : self.get_ids.output_dir + '/ids.tsv',
            'species' : self.option("species_name"),
            'blast_nr_xml': xml_path + '/nr/blast.xml',
            'blast_uniprot_xml': xml_path + '/uniprot/blast.xml',
            'blast_eggnog_xml': xml_path + '/eggnog/blast.xml',
            'blast_kegg_xml': xml_path + '/kegg/blast.xml',
            'blast2go_annot': xml_path + '/GO/blast2go_merge.xls',
            'pfam_domain': pfam_dir + '/pfam_domain',
        }

        self.merge_id.on('start', self.set_step, {'start': self.step.merge_id})
        self.merge_id.on('end', self.set_step, {'end': self.step.merge_id})
        self.merge_id.set_options(options)
        self.merge_id.run()


    def run_ref_class_medical(self):
        '''
        参考注释分类
        '''
        merge_id_dir = self.merge_id.output_dir
        pfam_dir = self.orf_pfam.output_dir
        options = {
            'db': 'nr,uniprot,kegg,eggnog,pfam,go,reactome,do,disgenet',
            'taxonomy': self.option('taxonomy'),
            'nr_ids': merge_id_dir + "/" + 'nr/nr_annot.tsv',
            'uniprot_ids': merge_id_dir + "/" + 'uniprot/uniprot_annot.tsv',
            'eggnog_ids': merge_id_dir + "/" + 'eggnog/eggnog_annot.tsv',
            'kegg_ids': merge_id_dir + "/" + 'kegg/kegg_annot.tsv',
            'pfam_domain': merge_id_dir + "/" + 'pfam/pfam_annot.tsv',


            "blast2go_annot" : merge_id_dir + "/go/all2go_annot.xls",
            'go_ids': merge_id_dir + "/" + 'go/all2go_annot.xls',
            'do_ids': merge_id_dir + "/" + 'do/do_annot.tsv',
            'reactome_ids': merge_id_dir + "/" + 'reactome/reactome_annot.tsv',
            'disgenet_ids': merge_id_dir + "/" + 'disgenet/disgenet_annot.tsv',
            "gtf" : self.option("gtf"),
            "fasta": self.option("ref_fa"),
            # "gene2trans" : ref_annot_dir + "/refannot_class/all_tran2gene.txt",
            "g2t2p" : self.option("g2t2p"),
            'link_bgcolor': 'yellow',
            'png_bgcolor': 'FFFF00',
            'merge_type': self.option("merge_type"),
            "des" : self.option("biomart"),
            "des_type" : self.option("biomart_type"),
            "type": "ref",
            "nr_version": self.option("nr_version"),
            "kegg_version": self.option("kegg_version"),
            "kegg_species": self.species2abr[self.option("species_name")],
            "cog_version": self.option("eggnog_version"),
            "go_version": self.option("go_version"),
            "enterz" : self.get_ids.output_dir + '/ids.tsv',
        }

        # if os.path.exists(self.option("known_ko")):
        #     options.update({
        #         "known_ko": self.option("known_ko")
        #     })

        
        self.ref_class_medical.on('start', self.set_step, {'start': self.step.ref_class_medical})
        self.ref_class_medical.on('end', self.set_step, {'end': self.step.ref_class_medical})
        self.ref_class_medical.set_options(options)
        self.ref_class_medical.run()


    def run_ref_filter(self):
        '''
        根据筛选参数过滤参考注释
        '''
        merge_id_dir = self.merge_id.output_dir
        orf_pfam_dir = self.orf_pfam.output_dir

        options = {
            "blast_nr_xml" : merge_id_dir + "/nr/blast.xml",
            "blast_eggnog_xml" : merge_id_dir + "/eggnog/blast.xml",
            "blast_kegg_xml": merge_id_dir + "/kegg/blast.xml",
            "blast_swissprot_xml" : merge_id_dir + "/swissprot/blast.xml",
            "pfam_domain" : orf_pfam_dir + "/pfam_domain",
            "blast2go_annot" : merge_id_dir + "/GO/blast2go_merge.xls",
            'nr_evalue': self.option('nr_filter_evalue'),
            'swissprot_evalue': self.option('swissprot_filter_evalue'),
            'eggnog_evalue': self.option('cog_filter_evalue'),
            'kegg_evalue': self.option('kegg_filter_evalue'),
            'pfam_evalue': self.option('pfam_filter_evalue'),
        }
        
        self.ref_filter.on('start', self.set_step, {'start': self.step.ref_filter})
        self.ref_filter.on('end', self.set_step, {'end': self.step.ref_filter})
        self.ref_filter.set_options(options)
        self.ref_filter.run()


    def run_ref_class(self):
        '''
        参考注释分类
        '''
        ref_filter_dir = self.map_db.output_dir
        pfam_dir = self.orf_pfam.output_dir
        options = {
            'taxonomy': self.option('taxonomy'),
            'blast_nr_xml': ref_filter_dir + "/nr/blast.xml",
            'blast_kegg_xml' :ref_filter_dir + "/kegg/blast.xml",
            'blast_eggnog_xml': ref_filter_dir + "/eggnog/blast.xml",
            'blast_swissprot_xml': ref_filter_dir + "/swissprot/blast.xml",
            'pfam_domain': pfam_dir + "/pfam_domain",
            "blast2go_annot" : ref_filter_dir + "/GO/blast2go_merge.xls",
            "gtf" : self.option("gtf"),
            "fasta": self.option("ref_fa"),
            # "gene2trans" : ref_annot_dir + "/refannot_class/all_tran2gene.txt",
            "g2t2p" : self.option("g2t2p"),
            'link_bgcolor': 'yellow',
            'png_bgcolor': 'FFFF00',
            'merge_type': self.option("merge_type"),
            "des" : self.option("biomart"),
            "des_type" : self.option("biomart_type"),
            "type": "ref",
            "nr_version": self.option("nr_version"),
            "kegg_version": self.option("kegg_version"),
            "cog_version": self.option("eggnog_version"),
            "go_version": self.option("go_version"),
            "enterz" : self.option("enterz"),
        }

        if os.path.exists(self.option("known_ko")):
            options.update({
                "known_ko": self.option("known_ko")
            })

        self.ref_class.on('start', self.set_step, {'start': self.step.ref_class})
        self.ref_class.on('end', self.set_step, {'end': self.step.ref_class})
        self.ref_class.set_options(options)
        self.ref_class.run()


    def run(self):
        # self.map_db.on('end', self.run_merge_id)
        super(RefDbAnnotationModule, self).run()
        self.on_rely([self.get_ids, self.map_db, self.orf_pfam], self.run_merge_id)

        self.on_rely([self.map_db, self.orf_pfam], self.run_ref_class)
        self.merge_id.on('end', self.run_ref_class_medical)
        self.run_map_db()
        self.run_get_ids()
        self.run_orf_pfam()
        self.on_rely([self.ref_class, self.ref_class_medical], self.set_output)
        # self.ref_class.on('end', self.set_output)

    def set_output(self):
        CopyFile().linkdir(self.orf_pfam.output_dir, os.path.join(self.output_dir, "annot_orfpfam"))
        CopyFile().linkdir(self.map_db.output_dir, os.path.join(self.output_dir, "annot_map_db"))
        CopyFile().linkdir(self.ref_class.output_dir, os.path.join(self.output_dir, "annot_class2"))
        CopyFile().linkdir(self.ref_class_medical.output_dir, os.path.join(self.output_dir, "annot_class_medical"))
        CopyFile().linkdir(self.merge_id.output_dir, os.path.join(self.output_dir, "annot_merge_id"))
        self.annotdir_change()
        self.end()

    def annotdir_change(self):
        # 修改为 1 的目录结构

        v2_to_v1 = {
            "cog/cog_list_tran.xls": "cog/cog.xls",
            "cog/cog_venn_gene.txt": "anno_stat/venn/gene_cog_venn.txt",
            "cog/cog_venn_tran.txt": "anno_stat/venn/cog_venn.txt",
            "cog/summary.G.tsv": "anno_stat/cog_stat/gene_cog_summary.xls",
            "cog/summary.T.tsv": "cog/cog_summary.xls",
            "go/go_lev2_gene.stat.xls": "anno_stat/go_stat/gene_go12level_statistics.xls",
            "go/go_lev2_tran.stat.xls": "go/go12level_statistics.xls",
            "go/go_lev3_gene.stat.xls": "anno_stat/go_stat/gene_go123level_statistics.xls",
            "go/go_lev3_tran.stat.xls": "go/go123level_statistics.xls",
            "go/go_lev4_gene.stat.xls": "anno_stat/go_stat/gene_go1234level_statistics.xls",
            "go/go_lev4_tran.stat.xls": "go/go1234level_statistics.xls",
            "go/go_list_gene.xls": "anno_stat/go_stat/gene_gos.list",
            "go/go_list_tran.xls": "go/query_gos.list",
            "go/go_venn_gene.txt": "anno_stat/venn/gene_go_venn.txt",
            "go/go_venn_tran.txt": "anno_stat/venn/go_venn.txt",
            "kegg/kegg_gene_gene.xls": "anno_stat/kegg_stat/gene_kegg_table.xls",
            "kegg/kegg_gene_tran.xls": "kegg/kegg_table.xls",
            "kegg/kegg_layer_gene.xls": "anno_stat/kegg_stat/gene_kegg_layer.xls",
            "kegg/kegg_layer_tran.xls": "kegg/kegg_layer.xls",
            "kegg/kegg_pathway_gene.xls": "anno_stat/kegg_stat/gene_pathway_table.xls",
            "kegg/kegg_pathway_tran.xls": "kegg/pathway_table.xls",
            "kegg/kegg_venn_gene.txt": "anno_stat/venn/gene_kegg_venn.txt",
            "kegg/kegg_venn_tran.txt": "anno_stat/venn/kegg_venn.txt",
            "nr/nr_blast_gene.xls": "anno_stat/blast/gene_nr.xls",
            "nr/nr_blast_tran.xls": "anno_stat/blast/nr.xls",
            "nr/nr_venn_gene.txt": "anno_stat/venn/gene_nr_venn.txt",
            "nr/nr_venn_tran.txt": "anno_stat/venn/nr_venn.txt",
            "pfam/pfam_domain_gene.xls": "anno_stat/pfam_stat/gene_pfam_domain",
            "pfam/pfam_domain_tran.xls": "pfam_domain",
            "pfam/pfam_venn_gene.txt": "anno_stat/venn/gene_pfam_venn.txt",
            "pfam/pfam_venn_tran.txt": "anno_stat/venn/pfam_venn.txt",
            "swissprot/swissprot_blast_gene.xls": "anno_stat/blast/swissprot.xls",
            "swissprot/swissprot_blast_tran.xls": "anno_stat/blast/gene_swissprot.xls",
            "swissprot/swissprot_venn_gene.txt": "anno_stat/venn/gene_swissprot_venn.txt",
            "swissprot/swissprot_venn_tran.txt": "anno_stat/venn/swissprot_venn.txt",
            "kegg/kegg_pathway_gene_dir": "anno_stat/kegg_stat/gene_pathway",
            "kegg/kegg_pathway_tran_dir": "kegg/pathways",
            "all_tran2gene.txt": "tran2gene.txt",
            "all_annot.xls": "anno_stat/all_anno_detail.xls",
            "all_stat.xls": "anno_stat/all_annotation_statistics.xls"
        }
        for k,v in v2_to_v1.items():
            CopyFile().linkdir(os.path.join(self.ref_class.output_dir, k), os.path.join(os.path.join(self.output_dir, "annot_class"), v))

        v2_to_v1_xml = {
            "SwissprotAnnot/blast.xml.filter.T.xml": "blast_xml/swissprot.xml",
            "SwissprotAnnot/blast.xml.filter.G.xml": "anno_stat/blast/gene_swissprot.xml",
            "NrAnnot/blast.xml.filter.T.xml": "blast_xml/nr.xml",
            "NrAnnot/blast.xml.filter.G.xml": "anno_stat/blast/gene_nr.xml",
            "KeggAnnot/blast.xml.filter.T.xml": "blast_xml/kegg.xml",
            "KeggAnnot/blast.xml.filter.G.xml": "anno_stat/blast/gene_kegg.xml",
            "EggnogAnnot/blast.xml.filter.T.xml": "blast_xml/eggnog.xml"
        }
        for k,v in v2_to_v1_xml.items():
            CopyFile().linkdir(os.path.join(self.ref_class.work_dir, k), os.path.join(os.path.join(self.output_dir, "annot_class"), v))

    def end(self):
        super(RefDbAnnotationModule, self).end()
