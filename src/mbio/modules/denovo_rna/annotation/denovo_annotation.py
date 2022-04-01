# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from __future__ import division
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.annotation.denovo_anno_stat.all_annotation_stat import AllAnnoStat
import os
import shutil


class DenovoAnnotationModule(Module):
    """
    module for denovorna annotation
    """
    def __init__(self, work_id):
        super(DenovoAnnotationModule, self).__init__(work_id)
        options = [
            {"name": "blast_nr_xml", "type": "infile", "format": "align.blast.blast_xml"},
            {"name": "blast_string_xml", "type": "infile", "format": "align.blast.blast_xml"},
            {"name": "blast_kegg_xml", "type": "infile", "format": "align.blast.blast_xml"},
            {"name": "blast_nr_table", "type": "infile", "format": "align.blast.blast_table"},
            {"name": "blast_string_table", "type": "infile", "format": "align.blast.blast_table"},
            {"name": "blast_kegg_table", "type": "infile", "format": "align.blast.blast_table"},
            {"name": "gene_file", "type": "infile", "format": "denovo_rna.express.gene_list"},
            {"name": "anno_statistics", "type": "bool", "default": True},
            {"name": "go_annot", "type": "bool", "default": True},
            {"name": "nr_annot", "type": "bool", "default": True},
            {"name": "gene_go_list", "type": "outfile", "format": "annotation.go.go_list"},
            {"name": "gene_kegg_table", "type": "outfile", "format": "annotation.kegg.kegg_table"},
            {"name": "gene_go_level_2", "type": "outfile", "format": "annotation.go.level2"},
        ]
        self.add_option(options)
        self.blast_stat_nr = self.add_tool('align.ncbi.blaststat')
        self.ncbi_taxon = self.add_tool('taxon.ncbi_taxon')
        self.go_annot = self.add_tool('annotation.go_annotation')
        self.string_cog = self.add_tool('annotation.string2cog')
        self.kegg_annot = self.add_tool('annotation.kegg_annotation')
        self.anno_stat = self.add_tool('annotation.denovo_anno_stat')
        self.step.add_steps('blast_statistics', 'go_annot', 'kegg_annot', 'cog_annot', 'taxon_annot', 'anno_stat')

    def check_options(self):
        if self.option('anno_statistics') and not self.option('gene_file').is_set:
            raise OptionError('运行注释统计的tool必须要设置gene_file')

    def set_step(self, event):
        if 'start' in event['data']:
            event['data']['start'].start()
        if 'end' in event['data']:
            event['data']['end'].finish()
        self.step.update()

    def run_blast(self):
        """
        """
        self.all_end_tool = []  # 所有尾部注释模块，全部结束后运行整体统计
        temp_options = {
            'query': self.option('query'),
            'query_type': 'nucl',
            'database': 'nr',
            'blast': 'blastx',
            'evalue': None,
            'num_threads': self.option('blast_threads'),
            'outfmt': 6
        }
        if 'nr' in self.anno_database or 'go' in self.anno_database:
            temp_options['evalue'] = self.option('nr_blast_evalue')
            self.blast_nr.set_options(temp_options)
            if 'nr' in self.anno_database:
                self.blast_nr.on('end', self.run_blast_stat)
                self.all_end_tool.append(self.blast_stat_nr)
                self.blast_nr.on('end', self.run_ncbi_taxon)
                self.all_end_tool.append(self.ncbi_taxon)
            if 'go' in self.anno_database:
                self.blast_nr.on('end', self.run_go_anno)
                self.all_end_tool.append(self.go_annot)
            self.blast_nr.on('start', self.set_step, {'start': self.step.blast_nr})
            self.blast_nr.on('end', self.set_step, {'end': self.step.blast_nr})
            self.blast_nr.on('end', self.set_output, 'nrblast')
            self.blast_nr.run()
        if 'cog' in self.anno_database:
            temp_options['database'] = 'string'
            temp_options['evalue'] = self.option('string_blast_evalue')
            self.blast_string.set_options(temp_options)
            self.blast_string.on('end', self.run_string2cog)
            self.all_end_tool.append(self.string_cog)
            self.blast_string.on('start', self.set_step, {'start': self.step.blast_string})
            self.blast_string.on('end', self.set_step, {'end': self.step.blast_string})
            self.blast_string.on('end', self.set_output, 'stringblast')
            self.blast_string.run()
        if 'kegg' in self.anno_database:
            temp_options['evalue'] = self.option('kegg_blast_evalue')
            temp_options['database'] = 'kegg'
            self.blast_kegg.set_options(temp_options)
            self.blast_kegg.on('end', self.run_kegg_anno)
            self.all_end_tool.append(self.kegg_annot)
            self.blast_kegg.on('start', self.set_step, {'start': self.step.blast_kegg})
            self.blast_kegg.on('end', self.set_step, {'end': self.step.blast_kegg})
            self.blast_kegg.on('end', self.set_output, 'keggblast')
            self.blast_kegg.run()
        if len(self.all_end_tool) > 1:
            self.on_rely(self.all_end_tool, self.run_annot_stat)
        elif len(self.all_end_tool) == 1:
            self.all_end_tool[0].on('end', self.run_annot_stat)
        else:
            self.logger.info('NEVER HERE')

    def run_annot_stat(self):
        """
        """
        opts = {'gene_file': self.option('gene_file'), 'database': ','.join(self.anno_database)}
        if 'kegg' in self.anno_database:
            opts['kegg_xml'] = self.option('blast_kegg_xml')
        if 'go' in self.anno_database:
            opts['gos_list'] = self.go_annot.option('golist_out')
            opts['blast2go_annot'] = self.go_annot.option('blast2go_annot')
        if 'cog' in self.anno_database:
            opts['string_xml'] = self.option('blast_string_xml')
            opts['cog_list'] = self.string_cog.option('cog_list')
            opts['cog_table'] = self.string_cog.option('cog_table')
        if 'nr' in self.anno_database:
            opts['nr_xml'] = self.option('blast_nr_xml')
            opts['nr_taxon_details'] = self.ncbi_taxon.option('taxon_out')
        self.anno_stat.set_options(opts)
        self.anno_stat.on('start', self.set_step, {'start': self.step.anno_stat})
        self.anno_stat.on('end', self.set_step, {'end': self.step.anno_stat})
        self.anno_stat.on('end', self.set_output, 'anno_stat')
        self.anno_stat.run()

    def run_kegg_anno(self):
        """
        """
        options = {
            'blastout': self.option('blast_kegg_xml')
        }
        self.kegg_annot.set_options(options)
        self.kegg_annot.on('start', self.set_step, {'start': self.step.kegg_annot})
        self.kegg_annot.on('end', self.set_step, {'end': self.step.kegg_annot})
        self.kegg_annot.on('end', self.set_output, 'kegg_annot')
        self.kegg_annot.run()

    def run_string2cog(self):
        options = {
            'blastout': self.option('blast_string_xml')
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
            'blastout': self.option('blast_nr_xml')
        }
        self.go_annot.set_options(options)
        self.go_annot.on('start', self.set_step, {'start': self.step.go_annot})
        self.go_annot.on('end', self.set_step, {'end': self.step.go_annot})
        self.go_annot.on('end', self.set_output, 'go_annot')
        self.go_annot.run()

    def run_blast_stat(self):
        """
        nr库比对结果统计函数
        """
        options = {
            'in_stat': self.option('blast_nr_xml')
        }
        self.blast_stat_nr.set_options(options)
        self.blast_stat_nr.on('start', self.set_step, {'start': self.step.blast_statistics})
        self.blast_stat_nr.on('end', self.set_step, {'end': self.step.blast_statistics})
        self.blast_stat_nr.on('end', self.set_output, 'blast_stat')
        self.blast_stat_nr.run()

    def run_ncbi_taxon(self):
        """
        """
        options = {
            'blastout': self.option('blast_nr_xml'),
            'blastdb': 'nr'
        }
        self.ncbi_taxon.set_options(options)
        self.ncbi_taxon.on('start', self.set_step, {'start': self.step.taxon_annot})
        self.ncbi_taxon.on('end', self.set_step, {'end': self.step.taxon_annot})
        self.ncbi_taxon.on('end', self.set_output, 'ncbi_taxon')
        self.ncbi_taxon.run()

    def run(self):
        super(DenovoAnnotationModule, self).run()
        # self.run_blast()
        self.all_end_tool = []  # 所有尾部注释模块，全部结束后运行整体统计
        self.anno_database = []
        if self.option('nr_annot'):
            self.anno_database.append('nr')
            self.all_end_tool.append(self.blast_stat_nr)
            self.all_end_tool.append(self.ncbi_taxon)
            self.run_blast_stat()
            self.run_ncbi_taxon()
        if self.option('go_annot'):
            self.anno_database.append('go')
            self.all_end_tool.append(self.go_annot)
            self.run_go_anno()
        if self.option('blast_string_xml').is_set:
            self.anno_database.append('cog')
            self.all_end_tool.append(self.string_cog)
            self.run_string2cog()
        if self.option('blast_kegg_xml').is_set:
            self.anno_database.append('kegg')
            self.all_end_tool.append(self.kegg_annot)
            self.run_kegg_anno()
        if len(self.all_end_tool) > 1:
            self.on_rely(self.all_end_tool, self.run_annot_stat)
        elif len(self.all_end_tool) == 1:
            self.all_end_tool[0].on('end', self.run_annot_stat)
        else:
            raise Exception('不需要进行任何注释工作')
            self.logger.info('NEVER HERE')

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'blast_stat':
            self.linkdir(obj.output_dir, 'blast_nr_statistics')
        elif event['data'] == 'ncbi_taxon':
            self.linkdir(obj.output_dir, 'ncbi_taxonomy')
        elif event['data'] == 'go_annot':
            self.linkdir(obj.output_dir, 'go')
        elif event['data'] == 'string_cog':
            self.linkdir(obj.output_dir, 'cog')
        elif event['data'] == 'kegg_annot':
            self.linkdir(obj.output_dir, 'kegg')
        elif event['data'] == 'anno_stat':
            self.linkdir(obj.output_dir, 'anno_stat')
            if 'kegg' in self.anno_database:
                self.option('gene_kegg_table', obj.option('gene_kegg_anno_table').prop['path'])
            if 'go' in self.anno_database:
                self.option('gene_go_list', obj.option('gene_go_list').prop['path'])
                self.option('gene_go_level_2', obj.option('gene_go_level_2').prop['path'])
            try:
                self.get_all_anno_stat(self.output_dir + '/anno_stat/all_annotation.xls')
            except Exception as e:
                self.logger.info("统计all_annotation出错：{}".format(e))
            self.end()
        else:
            pass

    def get_all_anno_stat(self, all_anno_path):
        # stat all_annotation.xls
        kwargs = {'outpath': all_anno_path, 'gene_list': self.option('gene_file').prop['gene_list']}
        for db in self.anno_database:
            if db == 'cog':
                kwargs['cog_list'] = self.string_cog.option('cog_list').prop['path']
            if db == 'go':
                kwargs['gos_list'] = self.go_annot.option('golist_out').prop['path']
            if db == 'kegg':
                kwargs['kegg_table'] = self.kegg_annot.option('kegg_table').prop['path']
            if db == 'nr':
                kwargs['blast_nr_table'] = self.option('blast_nr_table').prop['path']
                kwargs['nr_taxons'] = self.anno_stat.option('nr_taxons').prop['path']
        allstat = AllAnnoStat()
        allstat.get_anno_stat(**kwargs)

    def linkdir(self, olddir, newname, mode='link'):
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下
        """
        if not os.path.isdir(olddir):
            raise Exception('需要移动到output目录的文件夹不存在。')
        newdir = os.path.join(self.output_dir, newname)
        if not os.path.exists(newdir):
            if mode == 'link':
                shutil.copytree(olddir, newdir, symlinks=True)
            elif mode == 'copy':
                shutil.copytree(olddir, newdir)
            else:
                raise Exception('错误的移动文件方式，必须是\'copy\'或者\'link\'')
        else:
            allfiles = os.listdir(olddir)
            oldfiles = [os.path.join(olddir, i) for i in allfiles]
            newfiles = [os.path.join(newdir, i) for i in allfiles]
            for newfile in newfiles:
                if os.path.isfile(newfile) and os.path.exists(newfile):
                    os.remove(newfile)
                elif os.path.isdir(newfile) and os.path.exists(newfile):
                    shutil.rmtree(newfile)
            for i in range(len(allfiles)):
                if os.path.isfile(oldfiles[i]):
                    os.system('cp {} {}'.format(oldfiles[i], newfiles[i]))
                else:
                    os.system('cp -r {} {}'.format(oldfiles[i], newdir))

    def end(self):
        repaths = [
            [".", "", "DENOVO_RNA结果文件目录"],
            ['ncbi_taxonomy/query_taxons_detail.xls', 'xls', '序列详细物种分类文件'],
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
            ["/anno_stat/ncbi_taxonomy/", "dir", "nr统计结果目录"],
            ["/anno_stat/cog_stat/", "dir", "cog统计结果目录"],
            ["/anno_stat/go_stat/", "dir", "go统计结果目录"],
            ["/anno_stat/kegg_stat/", "dir", "kegg统计结果目录"],
            ["/anno_stat/blast/", "dir", "基因序列blast比对结果目录"],
            ["/anno_stat/blast_nr_statistics/", "dir", "基因序列blast比对nr库统计结果目录"],
            ["/anno_stat/blast/gene_kegg.xls", "xls", "基因序列blast比对kegg注释结果table"],
            ["/anno_stat/blast/gene_nr.xls", "xls", "基因序列blast比对nr注释结果table"],
            ["/anno_stat/blast/gene_nr.xls", "xls", "基因序列blast比对nr注释结果table"],
            ["/anno_stat/blast/gene_string.xml", "xml", "基因序列blast比对string注释结果xml"],
            ["/anno_stat/blast/gene_kegg.xml", "xml", "基因序列blast比对kegg注释结果xml"],
            ["/anno_stat/blast/gene_string.xml", "xml", "基因序列blast比对string注释结果xml"],
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
            ['/ncbi_taxonomy/gene_taxons_detail.xls', 'xls', '基因序列详细物种分类文件'],
            ["/anno_stat/blast_nr_statistics/gene_nr_evalue.xls", "xls", "基因序列blast结果E-value统计"],
            ["/anno_stat/blast_nr_statistics/gene_nr_similar.xls", "xls", "基因序列blast结果similarity统计"],
            ["/anno_stat/ncbi_taxonomy/gene_taxons.xls", "xls", "基因序列nr物种注释表"],
            ["/anno_stat/ncbi_taxonomy/query_taxons.xls", "xls", "nr物种注释表"],
            ["/anno_stat/all_annotation_statistics.xls", "xls", "注释统计总览表"],
            ["/anno_stat/all_annotation.xls", "xls", "注释统计表"],
        ]
        regexps = [
            # [r"nrblast/.+_vs_.+\.xml", "xml", "blast比对nr输出结果，xml格式"],
            # [r"nrblast/.+_vs_.+\.xls", "xls", "blast比对nr输出结果，表格(制表符分隔)格式"],
            # [r"stringblast/.+_vs_.+\.xml", "xml", "blast比对string输出结果，xml格式"],
            # [r"stringblast/.+_vs_.+\.xls", "xls", "blast比对string输出结果，表格(制表符分隔)格式"],
            # [r"keggblast/.+_vs_.+\.xml", "xml", "blast比对kegg输出结果，xml格式"],
            # [r"keggblast/.+_vs_.+\.xls", "xls", "blast比对kegg输出结果，表格(制表符分隔)格式"],
            [r"kegg/pathways/ko.\d+", 'pdf', '标红pathway图'],
            [r"/blast_nr_statistics/.*_evalue\.xls", "xls", "比对结果E-value分布图"],
            [r"/blast_nr_statistics/.*_similar\.xls", "xls", "比对结果相似度分布图"],
            ["^/anno_stat/ncbi_taxonomy/nr_taxon_stat", "xls", "nr物种分类统计表"],
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(DenovoAnnotationModule, self).end()
