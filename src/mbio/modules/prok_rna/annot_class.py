# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
from __future__ import division
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import os
import shutil
from mbio.packages.rna.annot_config import AnnotConfig

class AnnotClassModule(Module):
    """
    module for denovorna annotation
    """
    def __init__(self, work_id):
        super(AnnotClassModule, self).__init__(work_id)
        options = [
            {"name": "blast_nr_xml", "type": "infile", "format": "prok_rna.blast_xml"},
            {"name": "blast_string_xml", "type": "infile", "format": "prok_rna.blast_xml"},
            {"name": "blast_eggnog_xml", "type": "infile", "format": "prok_rna.blast_xml"},
            {"name": "blast_ncbicog_xml", "type": "infile", "format": "prok_rna.blast_xml"},
            {"name": "blast_kegg_xml", "type": "infile", "format": "prok_rna.blast_xml"},
            {"name": "blast_swissprot_xml", "type": "infile", "format": "prok_rna.blast_xml"},
            {"name": "blast_nr_table", "type": "infile", "format": "prok_rna.blast_table"},
            {"name": "blast_string_table", "type": "infile", "format": "prok_rna.blast_table"},
            {"name": "blast_kegg_table", "type": "infile", "format": "prok_rna.blast_table"},
            {"name": "blast_swissprot_table", "type": "infile", "format": "prok_rna.blast_table"},
            {"name": "pfam_domain", "type": "infile", "format": "prok_rna.kegg_list"},
            {"name": "gos_list_upload", "type": "infile", "format": "prok_rna.anno_upload"},   # 客户上传go注释文件
            {"name": "kos_list_upload", "type": "infile", "format": "prok_rna.anno_upload"},  # 客户上传kegg注释文件
            {"name": "gene_file", "type": "infile", "format": "prok_rna.gene_list"},
            {"name": "length_file", "type": "infile", "format": "prok_rna.cog_list"},  # 注释转录本序列的长度
            {"name": "new_gtf", "type": "infile", "format": "prok_rna.gtf"},  # 新转录本gtf文件
            {"name": "gtf", "type": "infile", "format": "prok_rna.gtf"},  # 参考基因组gtf文件
            {"name": "gene2trans", "type": "infile", "format": "prok_rna.common"},
            {"name": "anno_statistics", "type": "bool", "default": True},
            {"name": "go_annot", "type": "bool", "default": True},
            {"name": "nr_annot", "type": "bool", "default": False},  # 参考基因组注释不提供缺少nr_xml文件，因此将默认值改为False
            {"name": "taxonomy", "type": "string", "default": None},   # kegg数据库物种分类, Animals/Plants/Fungi/Protists/Archaea/Bacteria
            {"name": "link_bgcolor", "type": "string", "default": "yellow"},  # 通路图链接官网颜色，约定参考基因组为黄色（yellow），新序列为绿色(green), 两者共有为tomato（红）
            {"name": "png_bgcolor", "type": "string", "default": "#FFFF00"},  # 通路图静态图颜色，#00CD00(绿色)，#FFFF00（黄色）
            {"name": "gene_go_list", "type": "outfile", "format": "prok_rna.go_list"},
            {"name": "gene_kegg_table", "type": "outfile", "format": "prok_rna.kegg_table"},
            {"name": "kegg_table", "type": "outfile", "format": "prok_rna.kegg_table"},
            {"name": "gene_go_level_2", "type": "outfile", "format": "prok_rna.level2"},
            {"name": "des", "type": "infile", "format": "prok_rna.common"},
            {"name": "novel", "type": "string", "default": None},
            # 基因描述文件， 第一列为locus id, symbol  name    genomic_accession start  end strand  class  seq_type
            # {"name": "des_type", "type": "string", "default": "type3"},  # 基因描述文件类型
            # {"name": "enterz", "type": "string", "default": None},
            {"name": "annot_type", "type": "string", "default": "prok_rna" },
            # {"name": "g2t2p", "type": "string", "default": "" }, # 基因转录本蛋白对应关系列表
            {"name": "blast2go_annot", "type": "infile", "format": "prok_rna.blast2go_annot"},
            {"name": "nr_version", "type": "string", "default": "2019"},
            {"name": "go_version", "type": "string", "default": "2019"},
            {"name": "swissprot_version", "type": "string", "default": "2019"},
            {"name": "eggnog_version", "type": "string", "default": "2019"},
            {"name": "cog_version", "type": "string", "default": "2020"},
            {"name": "kegg_version", "type": "string", "default": "202003"},
            {"name": "string_version", "type": "string", "default": "2019"},
            {"name": "pir_version", "type": "string", "default":"2019"},
            {"name": "ncbi_taxonmy_version", "type": "string", "default":"2019"}

        ]
        self.add_option(options)
        self.swissprot_annot = self.add_tool('prok_rna.annotation.xml2table')
        self.nr_annot = self.add_tool('prok_rna.annotation.xml2table')
        self.eggnog_annot = self.add_tool('prok_rna.annotation.xml2table')
        self.ncbi_annot = self.add_tool('prok_rna.annotation.xml2table')

        self.eggnog_cog = self.add_tool('prok_rna.annotation.eggnog2cog')
        self.ncbi_cog = self.add_tool('prok_rna.annotation.ncbi2cog')
        self.ncbi_taxon = self.add_tool('prok_rna.annotation.ncbi_taxon')
        self.go_annot = self.add_tool('prok_rna.annotation.go_annotation')
        self.go_upload = self.add_tool('prok_rna.annotation.go_upload')
        self.string_cog = self.add_tool('prok_rna.annotation.string2cogv9')
        self.kegg_annot = self.add_tool('prok_rna.annotation.kegg_annotation')
        self.swissprot_annot = self.add_tool("prok_rna.annotation.swissprot")
        self.kegg_upload = self.add_tool('prok_rna.annotation.kegg_upload')
        self.anno_stat = self.add_tool('prok_rna.annotation.ref_anno_stat')
        self.anno_query = self.add_tool('prok_rna.annotation.ref_anno_query')
        self.step.add_steps('blast_statistics', 'nr_annot','go_annot', 'go_upload', 'kegg_annot', 'kegg_upload', 'eggnog_annot', 'eggnog_cog', 'ncbi_annot', 'ncbi_cog', 'cog_annot', 'anno_stat', 'swissprot_annot', 'anno_query')

    def check_options(self):
        if self.option('annot_type') == "denovo_rna":
            if self.option('anno_statistics'):
                if not self.option('gene2trans').is_set:
                    if self.option('gtf').is_set and self.option('g2t2p') != "":
                        tran2gene = self.work_dir + "/tran2gene.txt"
                        self.option("gtf").to_isoform2unigene(tran2gene, self.option("g2t2p"))
                        self.option('gene2trans', tran2gene)
                    else:
                        raise OptionError("进行注释统计的tool必须要设置gtf+g2t2p或者gene2trans", code = "25000101")
                if os.path.exists(self.output_dir + "/tran2gene.txt"):
                    os.remove(self.output_dir + "/tran2gene.txt")
                os.link(self.option("gene2trans").prop['path'], self.output_dir + "/tran2gene.txt")

        elif self.option('annot_type') == "ref_rna":
            if self.option('anno_statistics'):
                if not self.option('gene_file').is_set:
                    raise OptionError("进行注释统计的tool必须要设置gene_file", code = "25000102")
                if not self.option('ref_genome_gtf').is_set and not self.option('new_gtf').is_set:
                    raise OptionError("缺少gtf文件", code = "25000103")
                if not self.option("length_file").is_set:
                    raise OptionError("缺少注释转录本序列的长度文件", code = "25000104")

    def set_step(self, event):
        if 'start' in event['data']:
            event['data']['start'].start()
        if 'end' in event['data']:
            event['data']['end'].finish()
        self.step.update()

    def run_nr_anno(self):
        options = {
            'blastout': self.option('blast_nr_xml')
        }
        self.nr_annot.set_options(options)
        self.nr_annot.on('start', self.set_step, {'start': self.step.nr_annot})
        self.nr_annot.on('end', self.set_step, {'end': self.step.nr_annot})
        self.nr_annot.on('end', self.set_output, 'nr_annot')
        self.nr_annot.run()

    def run_eggnog_anno(self):
        options = {
            'blastout': self.option('blast_eggnog_xml'),


        }
        self.eggnog_annot.set_options(options)
        self.eggnog_annot.on('start', self.set_step, {'start': self.step.eggnog_annot})
        self.eggnog_annot.on('end', self.set_step, {'end': self.step.eggnog_annot})
        self.eggnog_annot.on('end', self.set_output, 'eggnog_annot')
        self.eggnog_annot.on('end', self.run_eggnog2cog)
        self.eggnog_annot.run()

    def run_eggnog2cog(self):
        option = {
            'blasttable': self.eggnog_annot.option("blast_table").prop['path'],
            'eggnog_version': self.option('eggnog_version')
        }
        self.eggnog_cog.set_options(option)
        self.eggnog_cog.on('start', self.set_step, {'start': self.step.eggnog_cog})
        self.eggnog_cog.on('end', self.set_step, {'end': self.step.eggnog_cog})
        self.eggnog_cog.on('end', self.set_output, 'eggnog_cog')
        self.eggnog_cog.run()


    def run_ncbi2cog_anno(self):
        options = {
            'blastout': self.option('blast_ncbicog_xml'),


        }
        self.ncbi_annot.set_options(options)
        self.ncbi_annot.on('start', self.set_step, {'start': self.step.ncbi_annot})
        self.ncbi_annot.on('end', self.set_step, {'end': self.step.ncbi_annot})
        self.ncbi_annot.on('end', self.set_output, 'ncbi_annot')
        self.ncbi_annot.on('end', self.run_ncbi2cog)
        self.ncbi_annot.run()

    def run_ncbi2cog(self):
        option = {
            'blasttable': self.ncbi_annot.option("blast_table").prop['path'],
            'cog_version': self.option('cog_version')
        }
        self.ncbi_cog.set_options(option)
        self.ncbi_cog.on('start', self.set_step, {'start': self.step.ncbi_cog})
        self.ncbi_cog.on('end', self.set_step, {'end': self.step.ncbi_cog})
        self.ncbi_cog.on('end', self.set_output, 'ncbi_cog')
        self.ncbi_cog.run()

    def run_ncbi_taxon(self):
        """
        注释物种分类
        """
        options = {
            'blastout': self.option('blast_nr_xml'),
            'blastdb': 'nr',
            'version': self.option('ncbi_taxonomy_version')
        }
        self.ncbi_taxon.set_options(options)
        self.ncbi_taxon.on('start', self.set_step, {'start': self.step.ncbi_taxon})
        self.ncbi_taxon.on('end', self.set_step, {'end': self.step.ncbi_taxon})
        self.ncbi_taxon.on('end', self.set_output, 'ncbi_taxon')
        self.ncbi_taxon.run()

    def run_swissprot_anno(self):
        options = {
            'blastout': self.option('blast_swissprot_xml')
        }
        self.swissprot_annot.set_options(options)
        self.swissprot_annot.on('start', self.set_step, {'start': self.step.swissprot_annot})
        self.swissprot_annot.on('end', self.set_step, {'end': self.step.swissprot_annot})
        self.swissprot_annot.on('end', self.set_output, 'swissprot_annot')
        self.swissprot_annot.run()


    def run_string2cog(self):
        options = {
            'blastout': self.option('blast_string_xml'),
            'string_table': self.option('blast_string_table')
        }
        self.string_cog.set_options(options)
        self.string_cog.on('start', self.set_step, {'start': self.step.cog_annot})
        self.string_cog.on('end', self.set_step, {'end': self.step.cog_annot})
        self.string_cog.on('end', self.set_output, 'string_cog')
        self.string_cog.run()

    def run_go_anno(self):
        """
        """
        options = {
            'blast2go_annot': self.option('blast2go_annot'),
            'go_version': self.option('go_version')
            # 'blastout': self.option('blast_nr_xml')
        }
        self.go_annot.set_options(options)
        self.go_annot.on('start', self.set_step, {'start': self.step.go_annot})
        self.go_annot.on('end', self.set_step, {'end': self.step.go_annot})
        self.go_annot.on('end', self.set_output, 'go_annot')
        self.go_annot.run()

    def run_go_upload(self):
        options = {
            'gos_list_upload': self.option('gos_list_upload')
        }
        self.go_upload.set_options(options)
        self.go_upload.on('start', self.set_step, {'start': self.step.go_upload})
        self.go_upload.on('end', self.set_step, {'end': self.step.go_upload})
        self.go_upload.on('end', self.set_output, 'go_annot')
        self.go_upload.run()

    def run_kegg_anno(self):
        """
        """
        options = {
            'blastout': self.option('blast_kegg_xml'),
            'taxonomy': self.option('taxonomy'),
            'kegg_version': self.option('kegg_version'),
            'link_bgcolor': self.option('link_bgcolor'),
            'png_bgcolor': self.option('png_bgcolor'),
            'kegg_version': self.option('kegg_version')
        }
        self.kegg_annot.set_options(options)
        self.kegg_annot.on('start', self.set_step, {'start': self.step.kegg_annot})
        self.kegg_annot.on('end', self.set_step, {'end': self.step.kegg_annot})
        self.kegg_annot.on('end', self.set_output, 'kegg_annot')
        self.kegg_annot.run()

    def run_pfam_stat(self):
        """
        pfam 统计
        """
        # self.set_output('pfam')
        pfam_base = os.path.basename(self.option("pfam_domain").prop['path'])
        if os.path.exists(self.output_dir + "/pfam_domain"):
            os.remove(self.output_dir + "/pfam_domain")
        os.link(self.option("pfam_domain").prop['path'],
                self.output_dir + "/pfam_domain")

    def run_kegg_upload(self):
        options = {
            'kos_list_upload': self.option('kos_list_upload'),
            'taxonomy': self.option('taxonomy'),
            'kegg_version': self.option('kegg_version'),
            'link_bgcolor': self.option('link_bgcolor'),
            'png_bgcolor': self.option('png_bgcolor')
        }
        self.kegg_upload.set_options(options)
        self.kegg_upload.on('start', self.set_step, {'start': self.step.kegg_upload})
        self.kegg_upload.on('end', self.set_step, {'end': self.step.kegg_upload})
        self.kegg_upload.on('end', self.set_output, 'kegg_annot')
        self.kegg_upload.run()

    def run_annot_stat(self):
        """
        """
        opts = {'gene_file': self.option('gene_file'), 'database': ','.join(self.anno_database)}
        if 'kegg' in self.anno_database:
            opts['kegg_xml'] = self.option('blast_kegg_xml')
            opts['kegg_anno_table'] = self.kegg_annot.option('kegg_table').prop['path']
            opts['kegg_version'] = self.option('kegg_version')
            opts['kos_list_upload'] = self.option('kos_list_upload')
        if 'go' in self.anno_database:
            opts['gos_list'] = self.go_annot.option('golist_out')
            opts['blast2go_annot'] = self.go_annot.option('blast2go_annot')
            opts['gos_list_upload'] = self.option('gos_list_upload')
        if 'cog' in self.anno_database:
            if self.option("blast_eggnog_xml").is_set:
                opts['eggnog_cog_table'] = self.eggnog_cog.output_dir + '/cog.xls'
            elif self.option("blast_ncbicog_xml").is_set:
                opts['eggnog_cog_table'] = self.ncbi_cog.output_dir + '/cog.xls'
            else:
                opts['string_xml'] = self.option('blast_string_xml')
                opts['string_table'] = self.option('blast_string_table')
                opts['cog_summary'] = self.string_cog.option('cog_summary')
        if 'nr' in self.anno_database:
            opts['nr_xml'] = self.option('blast_nr_xml')
            # opts['ncbi_taxon'] = self.ncbi_taxon.option("taxon_out").prop['path']
        opts['swissprot_xml'] = self.option('blast_swissprot_xml')
        opts['pfam_domain'] = self.option('pfam_domain')
        opts['des'] = self.option('des').prop['path']
        # opts['ref_genome_gtf'] = self.option('ref_genome_gtf')
        '''
        if self.option('gene2trans').is_set:
            opts['gene2trans'] = self.option('gene2trans').prop['path']
        elif self.option('new_gtf').is_set:
            opts['ref_genome_gtf'] = self.option('new_gtf').prop['path']
        elif self.option('ref_genome_gtf').is_set:
            opts['ref_genome_gtf'] = self.option('ref_genome_gtf').prop['path']
        elif self.option('gtf').is_set:
            tran2gene = self.work_dir + "/tran2gene.txt"
            self.option('gtf').to_isoform2unigene(tran2gene, self.option("g2t2p"))
            opts['gene2trans'] = tran2gene
        else:
            pass
        '''
        opts['taxonomy'] = self.option('taxonomy')
        self.anno_stat.set_options(opts)
        self.anno_stat.on('start', self.set_step, {'start': self.step.anno_stat})
        self.anno_stat.on('end', self.run_anno_query)
        self.anno_stat.on('end', self.set_step, {'end': self.step.anno_stat})
        self.anno_stat.on('end', self.set_output, 'anno_stat')
        self.anno_stat.run()


    def run_anno_query(self):
        if self.option('annot_type') == "denovo_rna":
            pass
        elif self.option('annot_type') == "prok_rna":
            opts = {
                'des': self.option('des').prop['path'],
                'kegg_version': self.option('kegg_version')
            }
            # opts = {'gene2trans': self.option('gene2trans').prop['path']}
        else:
            opts = {'length_path': self.option('length_file').prop['path']}
            if self.option('ref_genome_gtf').is_set:
                opts['ref_gtf_path'] = self.option('ref_genome_gtf').prop['path']
            if self.option('new_gtf').is_set:
                opts['new_gtf_path'] = self.option('new_gtf').prop['path']
            opts['gene_file'] = self.option('gene_file').prop['path']
        nr_path = self.anno_stat.work_dir + "/blast/nr.xls"
        swissprot_path = self.anno_stat.work_dir + "/blast/swissprot.xls"
        if self.option('novel'):
            opts['novel'] = self.option('novel')
        if os.path.exists(nr_path):
            opts['blast_nr_table'] = nr_path
        else:
            opts['blast_nr_table'] = None
        if os.path.exists(swissprot_path):
            opts['blast_swissprot_table'] = swissprot_path
        else:
            opts['blast_swissprot_table'] = None
        if self.option('pfam_domain').is_set:
            opts['pfam_domain'] = self.option('pfam_domain').prop['path']
        else:
            opts['pfam_domain'] = None
        if self.option("gos_list_upload").is_set:
            opts['gos_list'] = self.go_upload.option('golist_out').prop['path']
        else:
            opts['gos_list'] = self.go_annot.option('golist_out').prop['path']
            opts['gos_detail'] = self.go_annot.option('go_detail').prop['path']
        if self.option('blast_kegg_xml').is_set:
            opts['kegg_table'] = self.kegg_annot.option('kegg_table').prop['path']
        if self.option('kos_list_upload').is_set:
            opts['kegg_table'] = self.kegg_upload.option('kegg_table').prop['path']
        if 'cog' in self.anno_database:
            if self.option("blast_eggnog_xml").is_set:
                opts['cog_list'] = self.eggnog_cog.output_dir + '/cog.xls'
            elif self.option("blast_ncbicog_xml").is_set:
                opts['cog_list'] = self.ncbi_cog.output_dir + '/cog.xls'
        else:
            opts['cog_list'] = None
        # opts['des'] = self.option('des')
        # opts['des_type'] = self.option('des_type')
        # opts['enterz'] = self.option('enterz')

        self.anno_query.set_options(opts)
        self.anno_query.on('start', self.set_step, {'start': self.step.anno_query})
        self.anno_query.on('end', self.set_step, {'end': self.step.anno_query})
        self.anno_query.on('end', self.set_output, 'anno_query')
        self.anno_query.run()

    def run(self):
        super(AnnotClassModule, self).run()
        self.all_end_tool = []  # 所有尾部注释模块，全部结束后运行整体统计
        self.anno_database = []
        self.set_xml()
        if self.option('nr_annot'):
            self.anno_database.append('nr')
        if self.option('blast_nr_xml').is_set:
            self.all_end_tool.append(self.nr_annot)
            self.anno_database.append('nr')
            # self.all_end_tool.append(self.ncbi_taxon)
            self.run_nr_anno()
            # self.run_ncbi_taxon()
        if self.option('go_annot'):
            self.anno_database.append('go')
            self.all_end_tool.append(self.go_annot)
            if self.option("gos_list_upload").is_set:
                self.all_end_tool.append(self.go_upload)
                self.run_go_upload()
            else:
                self.run_go_anno()
        if self.option('blast_string_xml').is_set or self.option('blast_string_table').is_set:
            self.anno_database.append('cog')
            self.all_end_tool.append(self.string_cog)
            self.run_string2cog()
        if self.option('blast_eggnog_xml').is_set:
            self.anno_database.append('cog')
            self.all_end_tool.append(self.eggnog_cog)
            self.run_eggnog_anno()
        elif self.option('blast_ncbicog_xml').is_set:
            self.anno_database.append('cog')
            self.all_end_tool.append(self.ncbi_cog)
            self.run_ncbi2cog_anno()
        if self.option('blast_kegg_xml').is_set:
            self.anno_database.append('kegg')
            self.all_end_tool.append(self.kegg_annot)
            self.run_kegg_anno()
        if self.option('kos_list_upload').is_set:
            self.anno_database.append('kegg')
            self.all_end_tool.append(self.kegg_upload)
            self.run_kegg_upload()
        if self.option("blast_swissprot_xml").is_set:
            self.anno_database.append('swissprot')
            self.run_swissprot_anno()
        if self.option("pfam_domain").is_set:
            self.run_pfam_stat()
            self.anno_database.append('pfam')
        if len(self.all_end_tool) > 1:
            self.on_rely(self.all_end_tool, self.run_annot_stat)
        elif len(self.all_end_tool) == 1:
            self.all_end_tool[0].on('end', self.run_annot_stat)
        else:
            self.set_error("不需要进行任何注释工作", code = "25000105")
            self.logger.info('NEVER HERE')
    def set_xml(self):
        '''
        将xml文件链接到结果目录
        '''
        xml_dir = os.path.join(self.output_dir, "blast_xml")
        if os.path.exists(self.output_dir + '/blast_xml'):
            shutil.rmtree(self.output_dir + '/blast_xml')
        os.mkdir(xml_dir)
        # 上传xml结果
        os.link(self.option('blast_nr_xml').prop['path'], os.path.join(xml_dir, "vs_nr.xml") )
        os.link(self.option('blast_swissprot_xml').prop['path'], os.path.join(xml_dir, "vs_swissprot.xml") )
        if self.option('blast_eggnog_xml').is_set:
            os.link(self.option('blast_eggnog_xml').prop['path'], os.path.join(xml_dir, "vs_eggnog.xml") )
        elif self.option("blast_ncbicog_xml").is_set:
            os.link(self.option('blast_ncbicog_xml').prop['path'], os.path.join(xml_dir, "vs_ncbicog.xml") )
        os.link(self.option('blast_kegg_xml').prop['path'], os.path.join(xml_dir, "vs_kegg.xml") )
        os.link(self.option('pfam_domain').prop['path'], os.path.join(xml_dir, "pfam_domain") )


    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'blast_stat':
            self.linkdir(obj.output_dir, 'blast_nr_statistics')
        elif event['data'] == 'go_annot':
            self.linkdir(obj.output_dir, 'go')
        elif event['data'] == 'ncbi_taxon':
            self.linkdir(obj.output_dir, 'blast_nr_taxon')
        elif event['data'] == 'string_cog':
            self.linkdir(obj.output_dir, 'cog')
        elif event['data'] == 'eggnog_cog':
            self.linkdir(obj.output_dir, 'cog')
        elif event['data'] == 'ncbi_cog':
            self.linkdir(obj.output_dir, 'cog')
        elif event['data'] == 'kegg_annot':
            self.linkdir(obj.output_dir, 'kegg')
        elif event['data'] == 'anno_stat':
            self.linkdir(obj.output_dir, 'anno_stat')
            if os.path.exists(self.output_dir +  'cog/cog_summary.xls'):
                os.remove(self.output_dir +  'cog/cog_summary.xls')
            os.link(obj.work_dir + '/cog_summary.xls', self.output_dir +  '/cog/cog_summary.xls')
            if os.path.exists(self.output_dir +  'cog/cog_stat.xls'):
                os.remove(self.output_dir +  'cog/cog_stat.xls')
            os.link(obj.work_dir + '/cog_stat.xls', self.output_dir +  '/cog/cog_stat.xls')

            '''
            if 'kegg' in self.anno_database:
                self.option('gene_kegg_table').set_path(self.output_dir +
                                                        '/anno_stat/kegg_stat/gene_kegg_table.xls')
            if 'go' in self.anno_database:
                self.option('gene_go_list').set_path(self.output_dir +
                                                     '/anno_stat/go_stat/gene_gos.list')
                self.option('gene_go_level_2').set_path(self.output_dir +
                                                        '/anno_stat/go_stat/gene_go12level_statistics.xls')
            '''
        elif event['data'] == 'anno_query':
            if os.path.exists(self.output_dir + "/anno_stat/all_anno_detail.xls"):
                os.remove(self.output_dir + "/anno_stat/all_anno_detail.xls")
            os.link(self.anno_query.output_dir + "/all_anno_detail.xls",
                    self.output_dir + "/anno_stat/all_anno_detail.xls")
            if os.path.exists(self.output_dir + "/ref.txt"):
                os.remove(self.output_dir + "/ref.txt")
            os.link(self.option("des").prop['path'], self.output_dir + "/ref.txt")
            # if os.path.exists(self.output_dir + "/anno_stat/gene_anno_detail.xls"):
            #    os.remove(self.output_dir + "/anno_stat/gene_anno_detail.xls")
            # os.link(self.anno_query.output_dir + "/gene_anno_detail.xls",
            #         self.output_dir + "/anno_stat/gene_anno_detail.xls")
            self.end()
        else:
            pass

    def linkdir(self, olddir, newname, mode='link'):
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下
        """
        if not os.path.isdir(olddir):
            self.set_error("需要移动到output目录的文件夹不存在", code = "25000106")
        newdir = os.path.join(self.output_dir, newname)
        if os.path.exists(newdir):
            shutil.rmtree(newdir)
        os.mkdir(newdir)
        allfiles = os.listdir(olddir)
        oldfiles = [os.path.join(olddir, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            else:
                new1 = os.path.join(newdir, os.path.basename(oldfiles[i]))
                os.system("cp -r {} {}".format(oldfiles[i], new1))


    def end(self):
        repaths = [
            [".", "", "DENOVO_RNA结果文件目录"],
            ["blast_nr_statistics/output_evalue.xls", "xls", "blast结果E-value统计"],
            ["blast_nr_statistics/output_similar.xls", "xls", "blast结果similarity统计"],
            ["kegg/kegg_table.xls", "xls", "KEGG annotation table"],
            ["kegg/pathway_table.xls", "xls", "Sorted pathway table"],
            ["kegg/kegg_taxonomy.xls", "xls", "KEGG taxonomy summary"],
            ["go/blast2go.annot", "annot", "Go annotation based on blast output"],
            ["go/query_gos.list", "list", "Merged Go annotation"],
            ["go/go1234level_statistics.xls", "xls", "Go annotation on 4 levels"],
            ["go/go2level.xls", "xls", "Go annotation on level 2"],
            ["go/go3level.xls", "xls", "Go annotation on level 3"],
            ["go/go4level.xls", "xls", "Go annotation on level 4"],
            ["cog/cog_list.xls", "xls", "COG编号表"],
            ["cog/cog_summary.xls", "xls", "COG注释二级统计表"],
            ["cog/cog_table.xls", "xls", "序列COG注释详细表"],
            ["/anno_stat", "", "denovo注释统计结果输出目录"],
            ["/anno_stat/cog_stat/", "dir", "cog统计结果目录"],
            ["/anno_stat/go_stat/", "dir", "go统计结果目录"],
            ["/anno_stat/kegg_stat/", "dir", "kegg统计结果目录"],
            ["/anno_stat/blast/", "dir", "基因序列blast比对结果目录"],
            ["/anno_stat/blast_nr_statistics/", "dir", "基因序列blast比对nr库统计结果目录"],
            ["/anno_stat/blast/gene_kegg.xls", "xls", "基因序列blast比对kegg注释结果table"],
            ["/anno_stat/blast/gene_nr.xls", "xls", "基因序列blast比对nr注释结果table"],
            ["/anno_stat/blast/gene_nr.xls", "xls", "基因序列blast比对nr注释结果table"],
            ["/anno_stat/blast/gene_swissprot.xls", "xls", "基因序列blast比对到swissprot注释结果table"],
            ["/anno_stat/blast/gene_string.xml", "xml", "基因序列blast比对string注释结果xml"],
            ["/anno_stat/blast/gene_kegg.xml", "xml", "基因序列blast比对kegg注释结果xml"],
            ["/anno_stat/blast/gene_string.xml", "xml", "基因序列blast比对string注释结果xml"],
            ["/anno_stat/blast/gene_swissprot.xlm", "xml", "基因序列blast比对到swissprot注释结果xml"],
            ["/anno_stat/cog_stat/gene_cog_list.xls", "xls", "基因序列cog_list统计结果"],
            ["/anno_stat/cog_stat/gene_cog_summary.xls", "xls", "基因序列cog_summary统计结果"],
            ["/anno_stat/cog_stat/gene_cog_table.xls", "xls", "基因序列cog_table统计结果"],
            ["/anno_stat/cog_stat/gene_cog_table.xls", "xls", "基因序列cog_table统计结果"],
            ["/anno_stat/cog_stat/gene_cog_table.xls", "xls", "基因序列cog_table统计结果"],
            ["/anno_stat/cog_stat/gene_cog_table.xls", "xls", "基因序列cog_table统计结果"],
            ["/anno_stat/go_stat/gene_blast2go.annot", "annot", "Go annotation based on blast output of gene"],
            ["/anno_stat/go_stat/gene_gos.list", "list", "Merged Go annotation of gene"],
            ["/anno_stat/go_stat/gene_go1234level_statistics.xls", "xls", "Go annotation on 4 levels of gene"],
            ["/anno_stat/go_stat/gene_go2level.xls", "xls", "Go annotation on level 2 of gene"],
            ["/anno_stat/go_stat/gene_go3level.xls", "xls", "Go annotation on level 3 of gene"],
            ["/anno_stat/go_stat/gene_go4level.xls", "xls", "Go annotation on level 4 of gene"],
            ["/anno_stat/kegg_stat/gene_kegg_table.xls", "xls", "KEGG annotation table of gene"],
            ["/anno_stat/kegg_stat/gene_pathway_table.xls", "xls", "Sorted pathway table of gene"],
            ["/anno_stat/kegg_stat/gene_kegg_taxonomy.xls", "xls", "KEGG taxonomy summary of gene"],
            ["/anno_stat/kegg_stat/gene_kegg_layer.xls", "xls", "KEGG taxonomy summary of gene"],
            ["/anno_stat/kegg_stat/gene_pathway/", "dir", "基因的标红pathway图"],
            ["/anno_stat/blast_nr_statistics/gene_nr_evalue.xls", "xls", "基因序列blast结果E-value统计"],
            ["/anno_stat/blast_nr_statistics/gene_nr_similar.xls", "xls", "基因序列blast结果similarity统计"],
            ["/anno_stat/all_annotation_statistics.xls", "xls", "注释统计总览表"],
            ["/anno_stat/all_annotation.xls", "xls", "注释统计表"],
        ]
        regexps = [
            [r"kegg/pathways/ko.\d+", 'pdf', '标红pathway图'],
            [r"/blast_nr_statistics/.*_evalue\.xls", "xls", "比对结果E-value分布图"],
            [r"/blast_nr_statistics/.*_similar\.xls", "xls", "比对结果相似度分布图"],
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(AnnotClassModule, self).end()
