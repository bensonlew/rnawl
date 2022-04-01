# -*- coding: utf-8 -*-
# __author__ = 'qiuping, liubinxu'

from __future__ import division
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import re
import shutil
import traceback
import sys
from collections import defaultdict
from mbio.packages.align.blast.xml2table import xml2table
from mbio.packages.annotation.cog_stat import cog_stat
from mbio.packages.align.blast.blastout_statistics import *
from mbio.packages.ref_rna_v2.gene2trans import gene2trans
import unittest
from mbio.packages.rna.annot_config import AnnotConfig

class AnnotStatAgent(Agent):
    '''
    last_modify: 2019.02.25
    '''
    def __init__(self, parent):
        super(AnnotStatAgent, self).__init__(parent)
        options = [
            # essential options for determining the process
            {'name': 'database', 'type': 'string', 'default': 'nr,swissprot,kegg,eggnog,pfam,go'},
            {'name': 'gene2trans', 'type': 'infile', 'format': 'lnc_rna.common'},
            # input files of each database and related options
            {'name': 'nr_xml', 'type': 'infile', 'format': 'lnc_rna.blast_xml'}, # nr
            {'name': 'swissprot_xml', 'type': 'infile', 'format': 'lnc_rna.blast_xml'}, # swissprot
            {'name': 'eggnog_cog_table', 'type': 'infile', 'format': 'lnc_rna.cog_table'}, # eggnog
            {'name': 'kegg_xml', 'type': 'infile', 'format': 'lnc_rna.blast_xml'}, # kegg
            {'name': 'kegg_anno_table', 'type': 'infile', 'format': 'lnc_rna.kegg_table'},
            {'name': 'known_ko', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'taxonomy', 'type': 'string', 'default': None},
            {'name': 'link_bgcolor', 'type': 'string', 'default': 'green'},
            {'name': 'png_bgcolor', 'type': 'string', 'default': '#00CD00'},
            {'name': 'blast2go_annot', 'type': 'infile', 'format': 'lnc_rna.common'}, # go
            {'name': 'gos_list', 'type': 'infile', 'format': 'lnc_rna.go_list'},
            {'name': 'pfam_domain', 'type': 'infile', 'format': 'lnc_rna.common'}, # pfam

            # {"name": "ncbi_taxon", "type": "infile", "format": "ref_rna_v2.nr_taxon"},
            # {"name": "gos_list_upload", "type": "infile", "format": "ref_rna_v2.anno_upload"},
            # {"name": "kos_list_upload", "type": "infile", "format": "ref_rna_v2.anno_upload"},
            # {"name": "string_xml", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            # {"name": "string_table", "type": "infile", "format": "ref_rna_v2.blast_table"},
            # {"name": "cog_summary", "type": "infile", "format": "ref_rna_v2.cog_summary"},
            # {"name": "gene_file", "type": "infile", "format": "ref_rna_v2.gene_list"},
            # {"name": "ref_genome_gtf", "type": "infile", "format": "ref_rna_v2.gtf"},
            {"name": "gene_nr_table", "type": "outfile", "format": "ref_rna_v2.blast_table"},
            # {"name": "gene_string_table", "type": "outfile", "format": "ref_rna_v2.blast_table"},
            # {"name": "gene_kegg_table", "type": "outfile", "format": "ref_rna_v2.blast_table"},
            {"name": "gene_swissprot_table", "type": "outfile", "format": "ref_rna_v2.blast_table"},
            {"name": "gene_go_level_2", "type": "outfile", "format": "ref_rna_v2.level2"},
            {"name": "gene_go_list", "type": "outfile", "format": "ref_rna_v2.go_list"},
            {"name": "gene_kegg_anno_table", "type": "outfile", "format": "ref_rna_v2.kegg_table"},

            {"name": "gene_pfam_domain", "type": "outfile", "format": "ref_rna_v2.kegg_list"},
            # {"name": "gtf", "type": "infile", "format": "ref_rna_v2.gtf"},
            {'name': 'blast_nr_table', 'type': 'outfile', 'format': 'lnc_rna.common'},
            {'name': 'blast_swissprot_table', 'type': 'outfile', 'format': 'lnc_rna.common'},
            {'name': 'kegg_version', 'type': 'string', 'default': "202003"},
            {'name': 'go_version', 'type': 'string', 'default': "20200628"},
        ]
        self.add_option(options)
        self.step.add_steps('annot_stat')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)
        self.queue = 'BLAST2GO'

    def stepstart(self):
        self.step.annot_stat.start()
        self.step.update()

    def stepfinish(self):
        self.step.annot_stat.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        self.logger.debug('{} = {}'.format('database', self.option('database')))
        if len(self.option('database').split(',')) < 1:
            raise OptionError('number of selected database must more than zero')
        if self.option('gene2trans').is_set:
            self.logger.debug('{} = {}'.format('gene2trans', self.option('gene2trans').prop['path']))
        else:
            raise OptionError('input i2u file must be provided')
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 10
        self._memory = '80G'

    def end(self):
        super(AnnotStatAgent, self).end()


class AnnotStatTool(Tool):
    def __init__(self, config):
        super(AnnotStatTool, self).__init__(config)
        self.b2g_user = "biocluster102"
        self.b2g_password = "sanger-dev-123"
        self.python_path = "/program/Python/bin/python"
        self.denovo_stat = self.config.SOFTWARE_DIR + '/bioinfo/annotation/scripts/ref_stat/'
        # self.go_annot = self.config.SOFTWARE_DIR + '/bioinfo/annotation/scripts/goAnnot.py'
        self.go_annot = self.config.PACKAGE_DIR + "/rna/annotation/goAnnot.py"
        self.go_split = self.config.SOFTWARE_DIR + '/bioinfo/annotation/scripts/goSplit.py'
        self.kegg_path = self.config.PACKAGE_DIR + "/ref_rna_v2/kegg_annotation_v2.py"
        self.map_path = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/map4.r"
        self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.3/bin/Rscript"
        self.kegg_files_dict = AnnotConfig().get_file_dict(db="kegg", version=self.option("kegg_version"))
        self.html_path = self.kegg_files_dict['html'] + '/'
        self.cog_xml = self.config.PACKAGE_DIR + "/ref_rna_v2/String2Cog.pl"
        self.perl = '/program/perl-5.24.0/bin/perl'
        self.cog_sqlite_db = self.config.SOFTWARE_DIR + "/database/COG/cog.db"
        self.cog_table = self.config.SOFTWARE_DIR + '/bioinfo/annotation/scripts/cog_annot.py'
        self.image_magick = self.config.SOFTWARE_DIR + "/program/ImageMagick/bin/convert"
        self.taxonomy_path = self.config.SOFTWARE_DIR + "/database/KEGG/species/{}.ko.txt".format(self.option("taxonomy"))
        self.gene_nr_xml = self.work_dir + '/blast/gene_nr.xml'
        self.gene_swissprot_xml = self.work_dir + '/blast/gene_swissprot.xml'
        if self.option("kegg_xml").is_set:
            self.gene_kegg_xml = self.work_dir + '/blast/gene_kegg.xml'
        dir_list = ['/blast/', '/cog_stat/', '/go_stat/', '/pfam_stat/', '/kegg_stat/', '/blast_nr_statistics/', '/blast_swissprot_statistics/']
        for i in dir_list:
            if os.path.exists(self.work_dir + i):
                shutil.rmtree(self.work_dir + i)
            os.makedirs(self.work_dir + i)
        if not os.path.exists(self.output_dir + '/venn'):
            os.makedirs(self.output_dir + '/venn')
        if self.option("gene2trans").is_set:
            tran_gene = gene2trans().get_gene_transcript(trans2gene=self.option("gene2trans").prop["path"])
        else:
            pass
        self.tran_gene = tran_gene[0]
        self.tran_list = tran_gene[1]
        self.gene_list = tran_gene[2]
        self.database = set(self.option('database').split(','))
        if 'eggnog' in self.database:
            self.database.remove('eggnog')
            self.database.add('cog')
        self.gene_anno_list = {}  # 注释到的基因序列名字
        self.anno_list = {}  # 注释到的转录本序列名字
        self.go_obo = AnnotConfig().get_file_dict(db="go", version=self.option("go_version"))['go']


    def run(self):
        super(AnnotStatTool, self).run()
        for db in self.database:
            self.logger.info("开始处理数据库{}".format(db))
            if db == 'cog':
                self.run_cog_stat()
            if db == 'nr':
                self.run_nr_stat()
            if db == 'go':
                self.run_go_stat()
            if db == 'kegg':
                self.run_kegg_stat()
            if db == 'swissprot':
                self.run_swissprot_stat()
            if db == 'pfam':
                self.run_pfam_stat()
        if 'go' in self.database:
            self.wait()
            self.logger.info('end: go stat')
        self.set_output()
        self.get_all_anno_stat()
        self.end()

    def run_nr_stat(self):
        # 筛选gene_nr.xml、gene_nr.xls
        self.logger.info("开始筛选gene_nr.xml、gene_nr.xls")
        self.option('nr_xml').sub_blast_xml(genes=self.gene_list, new_fp=self.gene_nr_xml, trinity_mode=False)
        gene2trans().get_gene_blast_xml(tran_list=self.tran_list, tran_gene=self.tran_gene, xml_path=self.gene_nr_xml, gene_xml_path=self.gene_nr_xml)
        xml2table(self.gene_nr_xml, self.work_dir + '/blast/gene_nr.xls')
        xml2table(self.option('nr_xml').prop['path'], self.work_dir + '/blast/nr.xls')
        self.option('blast_nr_table').set_path(os.path.join(self.work_dir, 'blast/nr.xls'))
        self.logger.info("完成筛选gene_nr.xml、gene_nr.xls")

    def run_swissprot_stat(self):
        self.logger.info("开始筛选gene_swissprot.xml、gene_swissprot.xls")
        self.option('swissprot_xml').sub_blast_xml(genes=self.gene_list, new_fp=self.gene_swissprot_xml, trinity_mode=False)
        gene2trans().get_gene_blast_xml(tran_list=self.tran_list, tran_gene=self.tran_gene, xml_path=self.gene_swissprot_xml, gene_xml_path=self.gene_swissprot_xml)
        xml2table(self.gene_swissprot_xml, self.work_dir + '/blast/gene_swissprot.xls')
        xml2table(self.option('swissprot_xml').prop['path'], self.work_dir + '/blast/swissprot.xls')
        self.option('blast_swissprot_table').set_path(os.path.join(self.work_dir, 'blast/swissprot.xls'))
        self.logger.info("完成筛选gene_swissprot.xml、gene_swissprot.xls")
        try:
            blastout_statistics(blast_table=self.work_dir + '/blast/gene_swissprot.xls', evalue_path=self.work_dir + '/blast_swissprot_statistics/gene_swissprot_evalue.xls', similarity_path=self.work_dir + '/blast_swissprot_statistics/gene_swissprot_similar.xls')
            blastout_statistics(blast_table=self.work_dir + '/blast/swissprot.xls', evalue_path=self.work_dir + '/blast_swissprot_statistics/swissprot_evalue.xls', similarity_path=self.work_dir + '/blast_swissprot_statistics/swissprot_similar.xls')
            self.logger.info("End: evalue,similar for gene nr blast table ")
        except Exception as e:
            self.set_error("运行swissprot evalue,similar for gene swissprot blast table出错:%s", variables = (e), code = "33702613")
            self.logger.info("Error: evalue,similar for gene swissprot blast table")

    def run_cog_stat(self):
        # 筛选gene_string.xml、gene_string.xls
        self.cog_stat_path = self.work_dir + '/cog_stat/'
        self.logger.info("开始筛选gene_string.xml、gene_string.xls")
        if self.option("eggnog_cog_table").is_set:
            self.option("eggnog_cog_table").cog_class("cog_summary.xls")
            self.option("eggnog_cog_table").cog_class("cog_stat/gene_cog_summary.xls", self.option("gene2trans").prop["path"])
        elif self.option('string_xml').is_set:
            self.option('string_xml').sub_blast_xml(genes=self.gene_list, new_fp=self.gene_string_xml, trinity_mode=False)
            gene2trans().get_gene_blast_xml(tran_list=self.tran_list, tran_gene=self.tran_gene, xml_path=self.gene_string_xml, gene_xml_path=self.gene_string_xml)
            xml2table(self.gene_string_xml, self.work_dir + '/blast/gene_string.xls')
            self.logger.info("完成筛选gene_string.xml、gene_string.xls")
        else:
            self.option('string_table').sub_blast_table(genes=self.gene_list, new_fp=self.work_dir + '/gene_string.xls')
            gene2trans().get_gene_blast_table(tran_list=self.tran_list, tran_gene=self.tran_gene, table_path=self.work_dir + '/gene_string.xls', gene_table_path=self.gene_string_table)

    def run_pfam_stat(self):
        self.pfam_stat_path = self.work_dir + '/pfam_stat/'
        def get_gene_pfam(pfam_domain, gene_list, outpath):
            """
            将pfam注释的结果文件pfam_domain筛选致函基因的结果信息
            """
            with open(pfam_domain, 'rb') as f, open(outpath, 'wb') as w:
                lines = f.readlines()
                w.write(lines[0])
                for line in lines:
                    item = line.strip().split('\t')
                    name = item[0]
                    if name in gene_list:
                        line = re.sub(r"{}".format(name), self.tran_gene[name], line)
                        w.write(line)
        get_gene_pfam(pfam_domain=self.option('pfam_domain').prop['path'], gene_list=self.gene_list, outpath=self.pfam_stat_path + 'gene_pfam_domain')
        self.option('gene_pfam_domain', self.pfam_stat_path + 'gene_pfam_domain')

    def get_known_ko_gene(self, ko_file):
        with open(ko_file, 'r') as ko_in, open(self.work_dir + "/ko_known_gene", 'w') as ko_out:
            for line in ko_in:
                cols = line.split("\t")
                if cols[1] in self.gene_list:
                    cols[1] = cols[0]
                    ko_out.write("\t".join(cols))
        return self.work_dir + "/ko_known_gene"

    def run_kegg_stat(self):
        # 筛选gene_kegg.xml、gene_kegg.xls
        self.kegg_stat_path = self.work_dir + '/kegg_stat/'
        gene_pathway = self.kegg_stat_path + '/gene_pathway/'
        if os.path.exists(gene_pathway):
            shutil.rmtree(gene_pathway)
        os.makedirs(gene_pathway)
        if self.option("kegg_xml").is_set:
            self.logger.info("开始筛选gene_kegg.xml、gene_kegg.xls")
            self.option('kegg_xml').sub_blast_xml(genes=self.gene_list, new_fp=self.gene_kegg_xml, trinity_mode=False)
            gene2trans().get_gene_blast_xml(tran_list=self.tran_list, tran_gene=self.tran_gene, xml_path=self.gene_kegg_xml, gene_xml_path=self.gene_kegg_xml)
            xml2table(self.gene_kegg_xml, self.work_dir + '/blast/gene_kegg.xls')
            self.logger.info("完成筛选gene_kegg.xml、gene_kegg.xls")
        kegg_table = self.kegg_stat_path + '/gene_kegg_table.xls'
        pidpath = self.work_dir + '/gene_pid.txt'
        pathway_table = self.kegg_stat_path + '/gene_pathway_table.xls'
        layerfile = self.kegg_stat_path + '/gene_kegg_layer.xls'
        taxonomyfile = self.kegg_stat_path + '/gene_kegg_taxonomy.xls'
        if self.option("taxonomy"):
            taxonomy = self.taxonomy_path
        else:
            taxonomy = None
        if self.option("kegg_xml").is_set:
            if self.option("known_ko").is_set:
                known_ko_gene = self.get_known_ko_gene(self.option("known_ko").prop['path'])
            else:
                known_ko_gene = "None"
            cmd = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(self.python_path, self.kegg_path, self.r_path, self.map_path, self.gene_kegg_xml, None, kegg_table, pidpath, gene_pathway, pathway_table, layerfile, taxonomy, self.option("link_bgcolor"), self.option("png_bgcolor"), self.image_magick, self.html_path, known_ko_gene, self.option('kegg_version'))
        else:
            self.option("kos_list_upload").get_gene_anno(outdir=self.work_dir + "/gene_kegg.list")
            kegg_ids = self.work_dir + "/gene_kegg.list"
            cmd = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(self.python_path, self.kegg_path, self.r_path, self.map_path, None, kegg_ids, kegg_table, pidpath, gene_pathway, pathway_table, layerfile, taxonomy, self.option("link_bgcolor"), self.option("png_bgcolor"), self.image_magick)
        self.logger.info("开始运行kegg注释脚本")
        command = self.add_command("kegg_anno", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行kegg注释脚本完成")
        else:
            self.set_error("运行kegg注释脚本出错", code = "33702611")

    def run_go_stat(self):
        self.go_stat_path = self.work_dir + '/go_stat/'
        def get_gene_go(go_result, gene_list, outpath, trinity_mode=False):
            """
            将go_annotation注释的结果文件筛选只包含基因的结果信息,保留含有基因序列的行
            go_result:go_annotation tool运行得到的blast2go.annot或query_gos.list结果文件；
            gene_list: 只包含基因序列名字的列表
            """
            with open(go_result, 'rb') as c, open(outpath, 'wb') as w:
                for line in c:
                    line = line.strip('\n').split('\t')
                    name = line[0]
                    if name in gene_list:
                        if trinity_mode:
                            name = name.split('_i')[0]
                        line[0] = name
                        w_line = '\t'.join(line)
                        w.write(w_line + '\n')
        if self.option("blast2go_annot").is_set and self.option("gos_list").is_set:
            get_gene_go(go_result=self.option('blast2go_annot').prop['path'], gene_list=self.gene_list, outpath=self.go_stat_path + '/gene_blast2go.annot')
            get_gene_go(go_result=self.option('gos_list').prop['path'], gene_list=self.gene_list, outpath=self.go_stat_path + '/gene_gos.list')
            gene2trans().get_gene_go_list(tran_list=self.tran_list, tran_gene=self.tran_gene, go_list=self.go_stat_path + '/gene_blast2go.annot', gene_go_list=self.go_stat_path + '/gene_blast2go.annot')
            gene2trans().get_gene_go_list(tran_list=self.tran_list, tran_gene=self.tran_gene, go_list=self.go_stat_path + '/gene_gos.list', gene_go_list=self.go_stat_path + '/gene_gos.list')
        else:
            self.option("gos_list_upload").get_transcript_anno(outdir=self.work_dir + "/query_gos.list")
            self.option("gos_list_upload").get_gene_anno(outdir=self.go_stat_path + '/gene_gos.list')
            self.option("gos_list", self.work_dir + "/query_gos.list")
        self.option("gene_go_list", self.go_stat_path + '/gene_gos.list')
        go_cmd1 = '{} {} {} {}'.format(self.python_path, self.go_annot, self.go_obo, self.go_stat_path + '/gene_gos.list')
        go_annot_cmd = self.add_command('go_annot_cmd', go_cmd1).run()
        self.wait(go_annot_cmd)
        if go_annot_cmd.return_code == 0:
            self.logger.info("go_annot_cmd运行完成")
        else:
            self.set_error("go_annot_cmd运行出错!", code = "33702612")

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        self.movedir2output(self.work_dir + '/blast/', 'blast')
        for db in self.database:
            if db == 'cog':
                self.movedir2output(self.cog_stat_path, 'cog_stat')
                cog_venn = self.cog_summary_list(self.work_dir + "/cog_summary.xls")
                cog_gene_venn = self.cog_summary_list(self.work_dir + '/cog_stat/' + 'gene_cog_summary.xls')
                self.get_venn(venn_list=cog_venn, output=self.output_dir + '/venn/cog_venn.txt')
                self.get_venn(venn_list=cog_gene_venn, output=self.output_dir + '/venn/gene_cog_venn.txt')
                self.anno_list['cog'] = cog_venn
                self.gene_anno_list['cog'] = cog_gene_venn
                output_dir_cog = os.path.join(self.output_dir, 'cog')
                if os.path.isdir(output_dir_cog):
                    shutil.rmtree(output_dir_cog)
                os.mkdir(output_dir_cog)
                os.link(os.path.join(self.work_dir, 'cog_summary.xls'), os.path.join(output_dir_cog, 'cog_summary.xls'))
            if db == 'nr':
                self.option('gene_nr_table', self.output_dir + '/blast/gene_nr.xls')
                self.movedir2output(self.work_dir + '/blast_nr_statistics/', 'blast_nr_statistics')
                self.option('nr_xml').get_info()
                nr_venn = self.option('nr_xml').prop['hit_query_list']
                nr_gene_venn = self.option('gene_nr_table').prop['query_list']
                self.get_venn(venn_list=nr_venn, output=self.output_dir + '/venn/nr_venn.txt')
                self.get_venn(venn_list=nr_venn, output=self.output_dir + '/venn/nr_venn.txt')
                self.get_venn(venn_list=nr_gene_venn, output=self.output_dir + '/venn/gene_nr_venn.txt')
                self.anno_list['nr'] = nr_venn
                self.gene_anno_list['nr'] = nr_gene_venn
            if db == 'kegg':
                self.movedir2output(self.kegg_stat_path, 'kegg_stat')
                self.option('gene_kegg_anno_table', self.output_dir + '/kegg_stat/gene_kegg_table.xls')
                if self.option('kegg_anno_table').is_set:
                    kegg_venn = self.option('kegg_anno_table').get_query()
                if self.option('gene_kegg_anno_table').is_set:
                    kegg_gene_venn = self.option('gene_kegg_anno_table').get_query()
                else:
                    self.option("kos_list_upload").get_transcript_anno(outdir=self.work_dir + "/kegg.list")
                    kegg_venn = self.list_num(self.work_dir + "/kegg.list")
                    kegg_gene_venn = self.list_num(self.work_dir + "/gene_kegg.list")
                self.get_venn(venn_list=kegg_venn, output=self.output_dir + '/venn/kegg_venn.txt')
                self.get_venn(venn_list=kegg_gene_venn, output=self.output_dir + '/venn/gene_kegg_venn.txt')
                self.anno_list['kegg'] = kegg_venn
                self.gene_anno_list['kegg'] = kegg_gene_venn
            if db == 'pfam':
                self.movedir2output(self.pfam_stat_path, 'pfam_stat')
                pfam_venn = self.list_num(self.option('pfam_domain').prop['path'])
                gene_pfam_venn = self.list_num(self.option('gene_pfam_domain').prop['path'])
                self.get_venn(venn_list=pfam_venn, output=self.output_dir + '/venn/pfam_venn.txt')
                self.get_venn(venn_list=gene_pfam_venn, output=self.output_dir + '/venn/gene_pfam_venn.txt')
                self.anno_list['pfam'] = pfam_venn
                self.gene_anno_list['pfam'] = gene_pfam_venn
            if db == 'go':
                self.movedir2output(self.go_stat_path, 'go_stat')
                files = os.listdir(self.work_dir)
                for f in files:
                    if re.search(r'level_statistics\.xls$', f):
                        if os.path.exists(self.output_dir + '/go_stat/gene_{}'.format(f)):
                            os.remove(self.output_dir + '/go_stat/gene_{}'.format(f))
                        os.link(self.work_dir + '/' + f, self.output_dir + '/go_stat/gene_{}'.format(f))
                    if re.search(r'level.xls$', f):
                        if os.path.exists(self.output_dir + '/go_stat/gene_{}'.format(f)):
                            os.remove(self.output_dir + '/go_stat/gene_{}'.format(f))
                        os.link(self.work_dir + '/' + f, self.output_dir + '/go_stat/gene_{}'.format(f))
                self.option('gene_go_level_2', self.output_dir + '/go_stat/gene_go12level_statistics.xls')
                go_venn = self.list_num(self.option("gos_list").prop["path"])
                go_gene_venn = self.list_num(self.option("gene_go_list").prop["path"])
                self.get_venn(venn_list=go_venn, output=self.output_dir + '/venn/go_venn.txt')
                self.get_venn(venn_list=go_gene_venn, output=self.output_dir + '/venn/gene_go_venn.txt')
                self.anno_list['go'] = go_venn
                self.gene_anno_list['go'] = go_gene_venn
            if db == 'swissprot':
                self.option('gene_swissprot_table', self.output_dir + '/blast/gene_swissprot.xls')
                self.movedir2output(self.work_dir + '/blast_swissprot_statistics/', 'blast_swissprot_statistics')
                self.option('swissprot_xml').get_info()
                swissprot_venn = self.option('swissprot_xml').prop['hit_query_list']
                swissprot_gene_venn = self.option('gene_swissprot_table').prop['query_list']
                self.get_venn(venn_list=swissprot_venn, output=self.output_dir + '/venn/swissprot_venn.txt')
                self.get_venn(venn_list=swissprot_gene_venn, output=self.output_dir + '/venn/gene_swissprot_venn.txt')
                self.anno_list['swissprot'] = swissprot_venn
                self.gene_anno_list['swissprot'] = swissprot_gene_venn
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

    def movedir2output(self, olddir, newname, mode='link'):
        if not os.path.isdir(olddir):
            self.set_error('需要移动到output目录的文件夹不存在。', code = "33702615")
        newdir = os.path.join(self.output_dir, newname)
        self.logger.info(newdir)
        if os.path.exists(newdir):
            shutil.rmtree(newdir)
        os.mkdir(newdir)
        self.logger.info(newdir)
        allfiles = os.listdir(olddir)
        oldfiles = [os.path.join(olddir, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            else:
                newdir = os.path.join(newdir, os.path.basename(oldfiles[i]))
                if not os.path.exists(newdir):
                    os.mkdir(newdir)
                for f in os.listdir(oldfiles[i]):
                    old = os.path.join(oldfiles[i], f)
                    new = os.path.join(newdir, f)
                    if os.path.exists(new):
                        os.remove(new)
                    os.link(old, new)

    def get_venn(self, venn_list, output):
        with open(output, 'wb') as w:
            for i in venn_list:
                w.write(i + '\n')

    def get_all_anno_stat(self):
        all_anno_stat = self.output_dir + '/all_annotation_statistics.xls'
        anno_num = defaultdict(dict)
        with open(all_anno_stat, 'wb') as w:
            w.write('type\ttranscripts\tgenes\ttranscripts_percent\tgenes_percent\n')
            tmp = []
            tmp_gene = []
            tmp_annot = []
            tmp_annot_gene = []
            tmp_class = []
            tmp_class_gene = []
            anno_num['total']['gene'] = len(self.gene_list)
            anno_num['total']['tran'] = len(self.tran_list)
            db_type = {
                'nr': 'NR',
                'swissprot': 'Swiss-Prot',
                'pfam': 'Pfam',
                'cog': 'COG',
                'go': 'GO',
                'kegg': 'KEGG'
            }
            for db in ['go', 'kegg', 'cog', 'nr', 'swissprot', 'pfam']:
                if anno_num['total']['tran'] == 0:
                    tran_db_percent = 0
                else:
                    tran_db_percent = '%0.4g' % (len(self.anno_list[db]) / anno_num['total']['tran'])
                if anno_num['total']['gene'] == 0:
                    gene_db_percent = 0
                else:
                    gene_db_percent = '%0.4g' % (len(self.gene_anno_list[db]) / anno_num['total']['gene'])

                w.write('{}\t{}\t{}\t{}\t{}\n'.format(db_type[db], len(self.anno_list[db]), len(self.gene_anno_list[db]),  str(tran_db_percent), str(gene_db_percent)))
                tmp += self.anno_list[db]
                tmp_gene += self.gene_anno_list[db]
                if db == 'nr' or db == 'pfam' or db == 'swissprot':
                    tmp_annot += self.anno_list[db]
                    tmp_annot_gene += self.gene_anno_list[db]
                if db == 'go' or db == 'kegg' or db == 'cog':
                    tmp_class += self.anno_list[db]
                    tmp_class_gene += self.gene_anno_list[db]

            anno_num['total_anno']['gene'] = len(set(tmp_gene))
            anno_num['total_anno']['tran'] = len(set(tmp))
            anno_num['total_anno_nsp']['gene'] = len(set(tmp_annot_gene))
            anno_num['total_anno_nsp']['tran'] = len(set(tmp_annot))
            anno_num['total_class']['gene'] = len(set(tmp_class_gene))
            anno_num['total_class']['tran'] = len(set(tmp_class))
            if anno_num['total']['tran'] == 0:
                tran_total_percent = 0
            else:
                tran_total_percent = '%0.4g' % (anno_num['total_anno']['tran'] / anno_num['total']['tran'])
            if anno_num['total']['gene'] == 0:
                gene_total_percent = 0
            else:
                gene_total_percent = '%0.4g' % (anno_num['total_anno']['gene'] / anno_num['total']['gene'])
            w.write('Total_anno\t{}\t{}\t{}\t{}\n'.format(anno_num['total_anno']['tran'], anno_num['total_anno']['gene'], str(tran_total_percent),  str(gene_total_percent)))
            if anno_num['total']['tran'] == 0:
                total_annot_nsp_tran_percent = 0
            else:
                total_annot_nsp_tran_percent = '%0.4g' % (anno_num['total_anno_nsp']['tran'] / anno_num['total']['tran'])
            if anno_num['total']['gene'] == 0:
                total_annot_nsp_gene_percent = 0
            else:
                total_annot_nsp_gene_percent = '%0.4g' % (anno_num['total_anno_nsp']['gene'] / anno_num['total']['gene'])
            if anno_num['total']['tran'] == 0:
                total_class_tran_percent = 0
            else:
                total_class_tran_percent = '%0.4g' % (anno_num['total_class']['tran'] / anno_num['total']['tran'])
            if anno_num['total']['gene'] == 0:
                total_class_gene_percent = 0
            else:
                total_class_gene_percent = '%0.4g' % (anno_num['total_class']['gene'] / anno_num['total']['gene'])
            w.write('Total\t{}\t{}\t1\t1\n'.format(anno_num['total']['tran'], anno_num['total']['gene']))

    def list_num(self, list_file):
        with open(list_file, "rb") as f:
            lines = f.readlines()
            ids = []
            for line in lines:
                line = line.strip().split("\t")
                if line[0] not in ['Query_name', 'Seq_id']:
                    if line[0] not in ids:
                        ids.append(line[0])
        return list(set(ids))

    def cog_summary_list(self, cog_summary):
        query_ids = []
        with open(cog_summary, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                try:
                    query_cog = item[4].split(";")
                    for q in query_cog:
                        if q:
                            query_ids.append(q)
                except:
                    pass
                try:
                    query_nog = item[5].split(";")
                    for q in query_nog:
                        if q:
                            query_ids.append(q)
                except:
                    pass
            query_ids = list(set(query_ids))
        return query_ids

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'annot_stat_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.annotation.annot_stat',
            'instant': False,
            'options': {
                'database': 'nr,swissprot,eggnog,go,kegg,pfam',
                'gene2trans': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_class/annot_file/i2u',
                'nr_xml': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_filter/output/nr/blast.filter.xml',
                'swissprot_xml': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_filter/output/swissprot/blast.filter.xml',
                'eggnog_cog_table': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_class/eggnog_annot/cog.xls',
                'blast2go_annot': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_filter/output/go/blast2go.filter.xls',
                'gos_list': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_class/go_annot/query_gos.list',
                'kegg_xml': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_filter/output/kegg/blast.filter.xml',
                'known_ko': '/mnt/ilustre/users/isanger/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/KEGG/Homo_sapiens.GRCh38.pathway',
                'link_bgcolor': 'yellow',
                'png_bgcolor': '#FFFF00',
                'taxonomy': 'Animals',
                'pfam_domain': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_filter/output/pfam/pfam.filter.xls',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
