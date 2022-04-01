# -*- coding: utf-8 -*-
# __author__ = 'zengjing'

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.annotation.denovo_anno_stat.all_annotation_stat import AllAnnoStat
import os
import shutil


class DenovoAnnotationV2Module(Module):
    """
    module for denovorna annotation
    """
    def __init__(self, work_id):
        super(DenovoAnnotationV2Module, self).__init__(work_id)
        options = [
            {"name": "blast_nr_xml", "type": "infile", "format": "align.blast.blast_xml"},
            {"name": "blast_string_xml", "type": "infile", "format": "align.blast.blast_xml"},
            {"name": "blast_kegg_xml", "type": "infile", "format": "align.blast.blast_xml"},
            {"name": "blast_swissprot_xml", "type": "infile", "format": "align.blast.blast_xml"},
            {"name": "gene_file", "type": "infile", "format": "rna.gene_list"},
            {"name": "anno_statistics", "type": "bool", "default": True},
            {"name": "go_anno", "type": "bool", "default": True},
            {"name": "nr_anno", "type": "bool", "default": True},
            {"name": "gene_go_list", "type": "outfile", "format": "annotation.go.go_list"},
            {"name": "gene_kegg_table", "type": "outfile", "format": "annotation.kegg.kegg_table"},
            {"name": "gene_go_level_2", "type": "outfile", "format": "annotation.go.level2"},
        ]
        self.add_option(options)
        self.xml2table_nr = self.add_tool('align.xml2table')
        self.xml2table_sw = self.add_tool('align.xml2table')
        self.blast_stat_nr = self.add_tool('align.blaststat')
        self.blast_stat_sw = self.add_tool('align.blaststat')
        self.ncbi_taxon = self.add_tool('taxon.ncbi_taxon')
        self.go_anno = self.add_tool('annotation.go.go_annotation')
        self.string_cog = self.add_tool('annotation.cog.string2cog')
        self.kegg_anno = self.add_tool('annotation.kegg.kegg_annotation')
        # self.anno_stat = self.add_tool('annotation.denovo_anno_stat')
        self.step.add_steps('xml2table_nr', 'xml2table_sw', 'blast_stat_nr', 'blast_stat_sw', 'go_anno', 'kegg_anno', 'cog_annot', 'taxon_annot', 'anno_stat')

    def check_options(self):
        if self.option('anno_statistics') and not self.option('gene_file').is_set:
            raise OptionError('运行注释统计的tool必须要设置gene_file')

    def set_step(self, event):
        if 'start' in event['data']:
            event['data']['start'].start()
        if 'end' in event['data']:
            event['data']['end'].finish()
        self.step.update()

    def run_xml2table_nr(self):
        """nr blast结果xml文件转为table文件"""
        options = {
            'blastout': self.option("blast_nr_xml")
        }
        self.xml2table_nr.set_options(options)
        self.xml2table_nr.on('start', self.set_step, {'start': self.step.xml2table_nr})
        self.xml2table_nr.on('end', self.set_step, {'end': self.step.xml2table_nr})
        self.xml2table_nr.on('end', self.set_output, 'xml2table_nr')
        self.xml2table_nr.run()

    def run_blast_stat_nr(self):
        """nr blast结果evalue、similarity统计"""
        options = {
            'in_stat': self.option("blast_nr_xml")
        }
        self.blast_stat_nr.set_options(options)
        self.blast_stat_nr.on('start', self.set_step, {'start': self.step.blast_stat_nr})
        self.blast_stat_nr.on('end', self.set_step, {'end': self.step.blast_stat_nr})
        self.blast_stat_nr.on('end', self.set_output, 'blast_stat_nr')
        self.blast_stat_nr.run()

    def run_ncbi_taxon(self):
        """nr物种分类"""
        options = {
            'blastout': self.option('blast_nr_xml'),
            'blastdb': 'nr'
        }
        self.ncbi_taxon.set_options(options)
        self.ncbi_taxon.on('start', self.set_step, {'start': self.step.taxon_annot})
        self.ncbi_taxon.on('end', self.set_step, {'end': self.step.taxon_annot})
        self.ncbi_taxon.on('end', self.set_output, 'ncbi_taxon')
        self.ncbi_taxon.run()

    def run_xml2table_sw(self):
        """swissprot blast结果xml文件转为table文件"""
        options = {
            'blastout': self.option("blast_swissprot_xml")
        }
        self.xml2table_sw.set_options(options)
        self.xml2table_sw.on('start', self.set_step, {'start': self.step.xml2table_sw})
        self.xml2table_sw.on('end', self.set_step, {'end': self.step.xml2table_sw})
        self.xml2table_sw.on('end', self.set_output, 'xml2table_sw')
        self.xml2table_sw.run()

    def run_blast_stat_sw(self):
        """swissprot blast结果evalue、similarity统计"""
        options = {
            'in_stat': self.option("blast_swissprot_xml")
        }
        self.blast_stat_sw.set_options(options)
        self.blast_stat_sw.on('start', self.set_step, {'start': self.step.blast_stat_sw})
        self.blast_stat_sw.on('end', self.set_step, {'end': self.step.blast_stat_sw})
        self.blast_stat_sw.on('end', self.set_output, 'blast_stat_sw')
        self.blast_stat_sw.run()

    def run_string2cog(self):
        """string2cog注释"""
        options = {
            'blastout': self.option('blast_string_xml')
        }
        self.string_cog.set_options(options)
        self.string_cog.on('start', self.set_step, {'start': self.step.cog_annot})
        self.string_cog.on('end', self.set_step, {'end': self.step.cog_annot})
        self.string_cog.on('end', self.set_output, 'string_cog')
        self.string_cog.run()

    def run_go_anno(self):
        """go注释"""
        options = {
            'blastout': self.option('blast_nr_xml')
        }
        self.go_anno.set_options(options)
        self.go_anno.on('start', self.set_step, {'start': self.step.go_anno})
        self.go_anno.on('end', self.set_step, {'end': self.step.go_anno})
        self.go_anno.on('end', self.set_output, 'go_anno')
        self.go_anno.run()

    def run_kegg_anno(self):
        """kegg注释"""
        options = {
            'blastout': self.option('blast_kegg_xml')
        }
        self.kegg_anno.set_options(options)
        self.kegg_anno.on('start', self.set_step, {'start': self.step.kegg_anno})
        self.kegg_anno.on('end', self.set_step, {'end': self.step.kegg_anno})
        self.kegg_anno.on('end', self.set_output, 'kegg_anno')
        self.kegg_anno.run()

    def run(self):
        super(DenovoAnnotationV2Module, self).run()
        self.all_end_tool = []  # 所有尾部注释模块，全部结束后运行整体统计
        self.anno_database = []
        if self.option("blast_nr_xml").is_set:
            self.anno_database.append('nr')
            self.all_end_tool.append(self.xml2table_nr)
            self.all_end_tool.append(self.blast_stat_nr)
            self.all_end_tool.append(self.ncbi_taxon)
            self.run_ncbi_taxon()
            self.run_xml2table_nr()
            self.run_blast_stat_nr()
        if self.option("blast_swissprot_xml").is_set:
            self.anno_database.append('swissprot')
            self.all_end_tool.append(self.xml2table_sw)
            self.all_end_tool.append(self.blast_stat_sw)
            self.run_xml2table_sw()
            self.run_blast_stat_sw()
        if self.option('blast_string_xml').is_set:
            self.anno_database.append('cog')
            self.all_end_tool.append(self.string_cog)
            self.run_string2cog()
        if self.option('go_anno'):
            self.anno_database.append('go')
            self.all_end_tool.append(self.go_anno)
            self.run_go_anno()
        if self.option('blast_kegg_xml').is_set:
            self.anno_database.append('kegg')
            self.all_end_tool.append(self.kegg_anno)
            self.run_kegg_anno()
        # if len(self.all_end_tool) > 1:
        #     self.on_rely(self.all_end_tool, self.run_annot_stat)
        # elif len(self.all_end_tool) == 1:
        #     self.all_end_tool[0].on('end', self.run_annot_stat)
        else:
            raise Exception('不需要进行任何注释工作')
            self.logger.info('NEVER HERE')

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'xml2table':
            self.linkdir(obj.output_dir, 'xml2table')
        if event['data'] == 'blast_stat':
            self.linkdir(obj.output_dir, 'blast_statistics')
        elif event['data'] == 'ncbi_taxon':
            self.linkdir(obj.output_dir, 'ncbi_taxonomy')
        elif event['data'] == 'go_anno':
            self.linkdir(obj.output_dir, 'go')
        elif event['data'] == 'string_cog':
            self.linkdir(obj.output_dir, 'cog')
        elif event['data'] == 'kegg_anno':
            self.linkdir(obj.output_dir, 'kegg')
        # elif event['data'] == 'anno_stat':
        #     self.linkdir(obj.output_dir, 'anno_stat')
        #     if 'kegg' in self.anno_database:
        #         self.option('gene_kegg_table', obj.option('gene_kegg_anno_table').prop['path'])
        #     if 'go' in self.anno_database:
        #         self.option('gene_go_list', obj.option('gene_go_list').prop['path'])
        #         self.option('gene_go_level_2', obj.option('gene_go_level_2').prop['path'])
        #     try:
        #         self.get_all_anno_stat(self.output_dir + '/anno_stat/all_annotation.xls')
        #     except Exception as e:
        #         self.logger.info("统计all_annotation出错：{}".format(e))
        #     self.end()
        else:
            pass
        self.end()

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
        # repaths = [
        # ]
        # regexps = [
        # ]
        # sdir = self.add_upload_dir(self.output_dir)
        # sdir.add_relpath_rules(repaths)
        # sdir.add_regexp_rules(regexps)
        super(DenovoAnnotationV2Module, self).end()
