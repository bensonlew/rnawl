# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from __future__ import division
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
from mbio.packages.annotation.denovo_anno_stat.cog_stat import cog_stat
from mbio.packages.annotation.denovo_anno_stat.nr_stat import nr_stat
from mbio.packages.align.blast.blastout_statistics import *
import os
import re
import shutil
from biocluster.config import Config
import traceback
from collections import defaultdict


class DenovoAnnoStatAgent(Agent):
    """
    统计无参转录组注释模块NR、GO、COG、KEGG的相关信息
    version v1.0
    author: qiuping
    last_modify: 2016.10.20
    """
    def __init__(self, parent):
        super(DenovoAnnoStatAgent, self).__init__(parent)
        options = [
            {"name": "nr_xml", "type": "infile", "format": "align.blast.blast_xml"},  # blast比对到nr库的xml结果文件
            {"name": "nr_taxon_details", "type": "infile", "format": "annotation.nr.nr_taxon"},
            {"name": "blast2go_annot", "type": "infile", "format": "annotation.go.blast2go_annot"},
            {"name": "gos_list", "type": "infile", "format": "annotation.go.go_list"},
            {"name": "kegg_xml", "type": "infile", "format": "align.blast.blast_xml"},
            {"name": "string_xml", "type": "infile", "format": "align.blast.blast_xml"},
            {"name": "cog_list", "type": "infile", "format": "annotation.cog.cog_list"},
            {"name": "cog_table", "type": "infile", "format": "annotation.cog.cog_table"},
            {"name": "gene_file", "type": "infile", "format": "denovo_rna.express.gene_list"},
            {"name": "database", "type": "string", "default": "nr,go,cog,kegg"},
            {"name": "gene_nr_table", "type": "outfile", "format": "align.blast.blast_table"},
            {"name": "gene_string_table", "type": "outfile", "format": "align.blast.blast_table"},
            {"name": "gene_kegg_table", "type": "outfile", "format": "align.blast.blast_table"},
            {"name": "nr_taxons", "type": "outfile", "format": "annotation.nr.nr_taxon"},
            {"name": "gene_go_list", "type": "outfile", "format": "annotation.go.go_list"},
            {"name": "gene_go_level_2", "type": "outfile", "format": "annotation.go.level2"},
            {"name": "gene_kegg_anno_table", "type": "outfile", "format": "annotation.kegg.kegg_table"},
        ]
        self.add_option(options)
        self.step.add_steps("denovo_anno_stat")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)
        self.queue = 'BLAST2GO'  # 投递到指定的队列BLAST2GO

    def stepstart(self):
        self.step.denovo_anno_stat.start()
        self.step.update()

    def stepfinish(self):
        self.step.denovo_anno_stat.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        self.anno_database = set(self.option('database').split(','))
        if len(self.anno_database) < 1:
            raise OptionError('至少选择一种注释库')
        for i in self.anno_database:
            if i not in ['nr', 'go', 'cog', 'kegg']:
                raise OptionError('需要注释的数据库不在支持范围内[nr, go, cog, kegg]:{}'.format(i))
            if i == 'go' and not self.option('blast2go_annot').is_set and not self.option('gos_list').is_set:
                raise OptionError('缺少go注释的输入文件')
            if i == 'cog' and not self.option('string_xml').is_set and not self.option('cog_list').is_set and not self.option('cog_table').is_set:
                raise OptionError('缺少cog注释的输入文件')
            if i == 'nr' and not self.option('nr_xml').is_set and not self.option('nr_taxon_details').is_set:
                raise OptionError('缺少nr注释的输入文件')
            if i == 'kegg' and not self.option('kegg_xml').is_set:
                raise OptionError('缺少kegg注释的输入文件')
        if not self.option('gene_file').is_set:
            raise OptionError('缺少gene输入文件')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        relpath = [
            [".", "", "denovo注释统计结果输出目录"],
            ["/ncbi_taxonomy/", "dir", "nr统计结果目录"],
            ["/blast_nr_statistics/", "dir", "blast比对nr库evalue等值统计目录"],
            ["/cog_stat/", "dir", "cog统计结果目录"],
            ["/go_stat/", "dir", "go统计结果目录"],
            ["/kegg_stat/", "dir", "kegg统计结果目录"],
            ["/blast/", "dir", "基因序列blast比对结果目录"],
            ["/blast/gene_kegg.xls", "xls", "基因序列blast比对kegg注释结果table"],
            ["/blast/gene_nr.xls", "xls", "基因序列blast比对nr注释结果table"],
            ["/blast/gene_string.xls", "xls", "基因序列blast比对string注释结果table"],
            ["/blast/gene_string.xml", "xml", "基因序列blast比对string注释结果xml"],
            ["/blast/gene_kegg.xml", "xml", "基因序列blast比对kegg注释结果xml"],
            ["/blast/gene_string.xml", "xml", "基因序列blast比对string注释结果xml"],
            ["/cog_stat/gene_cog_list.xls", "xls", "基因序列cog_list统计结果"],
            ["/cog_stat/gene_cog_summary.xls", "xls", "基因序列cog_summary统计结果"],
            ["/cog_stat/gene_cog_table.xls", "xls", "基因序列cog_table统计结果"],
            ["/cog_stat/gene_cog_table.xls", "xls", "基因序列cog_table统计结果"],
            ["/cog_stat/gene_cog_table.xls", "xls", "基因序列cog_table统计结果"],
            ["/cog_stat/gene_cog_table.xls", "xls", "基因序列cog_table统计结果"],
            ["/go_stat/gene_blast2go.annot", "annot", "Go annotation based on blast output of gene"],
            ["/go_stat/gene_gos.list", "list", "Merged Go annotation of gene"],
            ["/go_stat/gene_go1234level_statistics.xls", "xls", "Go annotation on 4 levels of gene"],
            ["/go_stat/gene_go2level.xls", "xls", "Go annotation on level 2 of gene"],
            ["/go_stat/gene_go3level.xls", "xls", "Go annotation on level 3 of gene"],
            ["/go_stat/gene_go4level.xls", "xls", "Go annotation on level 4 of gene"],
            ["/kegg_stat/gene_kegg_table.xls", "xls", "KEGG annotation table of gene"],
            ["/kegg_stat/gene_pathway_table.xls", "xls", "Sorted pathway table of gene"],
            ["/kegg_stat/gene_kegg_taxonomy.xls", "xls", "KEGG taxonomy summary of gene"],
            ["/kegg_stat/gene_kegg_layer.xls", "xls", "KEGG taxonomy summary of gene"],
            ["/kegg_stat/gene_pathway/", "dir", "基因的标红pathway图"],
            ['/ncbi_taxonomy/gene_taxons_detail.xls', 'xls', '基因序列详细物种分类文件'],
            ["/blast_nr_statistics/gene_nr_evalue.xls", "xls", "基因序列blast结果E-value统计"],
            ["/blast_nr_statistics/gene_nr_similar.xls", "xls", "基因序列blast结果similarity统计"],
            ["/ncbi_taxonomy/gene_taxons.xls", "xls", "基因序列nr物种注释表"],
            ["/ncbi_taxonomy/query_taxons.xls", "xls", "nr物种注释表"],
        ]
        result_dir.add_relpath_rules(relpath)
        result_dir.add_regexp_rules([
            [r"^/ncbi_taxonomy/nr_taxon_stat", "xls", "不同水平的nr物种分类统计表"],
        ])
        super(DenovoAnnoStatAgent, self).end()


class DenovoAnnoStatTool(Tool):
    """
    表达量差异检测tool
    """
    def __init__(self, config):
        super(DenovoAnnoStatTool, self).__init__(config)
        self._version = '1.0.1'
        self.b2g_user = "biocluster102"
        self.b2g_password = "sanger-dev-123"
        self.python_path = "/miniconda2/bin/python"
        self.denovo_stat = self.config.SOFTWARE_DIR + '/bioinfo/annotation/scripts/denovo_stat/'
        self.go_annot = self.config.SOFTWARE_DIR + '/bioinfo/annotation/scripts/goAnnot.py'
        self.go_split = self.config.SOFTWARE_DIR + '/bioinfo/annotation/scripts/goSplit.py'
        self.kegg_anno = self.config.SOFTWARE_DIR + '/bioinfo/annotation/scripts/kegg_annotation.py'
        self.gene_list = self.option('gene_file').prop['gene_list']
        self.gene_nr_xml = self.work_dir + '/blast/gene_nr.xml'
        self.gene_string_xml = self.work_dir + '/blast/gene_string.xml'
        self.gene_kegg_xml = self.work_dir + '/blast/gene_kegg.xml'
        dir_list = ['/blast/', '/ncbi_taxonomy/', '/cog_stat/', '/go_stat/', '/kegg_stat/', '/blast_nr_statistics/']
        for i in dir_list:
            if os.path.exists(self.work_dir + i):
                shutil.rmtree(self.work_dir + i)
            os.makedirs(self.work_dir + i)
        if not os.path.exists(self.output_dir + '/venn'):
            os.makedirs(self.output_dir + '/venn')
        self.database = self.option('database').split(',')
        self.gene_anno_list = {}  # 注释到的基因序列名字
        self.anno_list = {}  # 注释到的转录本序列名字

    def run_nr_stat(self):
        self.nr_stat_path = self.work_dir + '/ncbi_taxonomy/'
        # 筛选gene_nr.xml、gene_nr.xls
        self.logger.info("开始筛选gene_nr.xml、gene_nr.xls")
        self.option('nr_xml').sub_blast_xml(genes=self.gene_list, new_fp=self.gene_nr_xml, trinity_mode=True)
        xml2table(self.gene_nr_xml, self.work_dir + '/blast/gene_nr.xls')
        self.logger.info("完成筛选gene_nr.xml、gene_nr.xls")
        # stat gene_evalue and gene_simillar for NR
        try:
            blastout_statistics(blast_table=self.work_dir + '/blast/gene_nr.xls', evalue_path=self.work_dir + '/blast_nr_statistics/gene_nr_evalue.xls', similarity_path=self.work_dir + '/blast_nr_statistics/gene_nr_similar.xls')
            self.logger.info("End: evalue,similar for gene nr blast table ")
        except Exception as e:
            self.set_error("运行nr evalue,similar for gene nr blast table出错:{}".format(e))
            self.logger.info("Error: evalue,similar for gene nr blast table")
        # stat taxon info
        nr = nr_stat()
        try:
            nr.stats(detail_file=self.option('nr_taxon_details').prop['path'], out_dir=self.nr_stat_path, gene_list=self.gene_list)
            self.logger.info("End: stat nr ")
        except Exception as e:
            self.set_error("运行nr stat evalue出错:{}".format(e))
            self.logger.info("Error: stat nr")

    def run_cog_stat(self):
        self.cog_stat_path = self.work_dir + '/cog_stat/'
        # 筛选gene_string.xml、gene_string.xls
        self.logger.info("开始筛选gene_string.xml、gene_string.xls")
        self.option('string_xml').sub_blast_xml(genes=self.gene_list, new_fp=self.gene_string_xml, trinity_mode=True)
        xml2table(self.gene_string_xml, self.work_dir + '/blast/gene_string.xls')
        self.logger.info("完成筛选gene_string.xml、gene_string.xls")
        # cog stat
        cog = cog_stat()
        try:
            cog.stats(cog_list=self.option('cog_list').prop['path'], gene_list=self.gene_list, gene_cog_list=self.cog_stat_path + 'gene_cog_list.xls', cog_table=self.option('cog_table').prop['path'], out_dir=self.cog_stat_path)
            self.logger.info("End: stat cog ")
        except Exception as e:
            self.set_error("运行nr stat evalue出错:{}".format(e))
            self.logger.info("Error: stat cog")

    def run_kegg_stat(self):
        self.kegg_stat_path = self.work_dir + '/kegg_stat/'
        gene_pathway = self.kegg_stat_path + '/gene_pathway/'
        if os.path.exists(gene_pathway):
            shutil.rmtree(gene_pathway)
        os.makedirs(gene_pathway)
        # 筛选gene_kegg.xml、gene_kegg.xls
        self.logger.info("开始筛选gene_kegg.xml、gene_kegg.xls")
        self.option('kegg_xml').sub_blast_xml(genes=self.gene_list, new_fp=self.gene_kegg_xml, trinity_mode=True)
        xml2table(self.gene_kegg_xml, self.work_dir + '/blast/gene_kegg.xls')
        self.logger.info("完成筛选gene_kegg.xml、gene_kegg.xls")
        # kegg_stat
        try:
            kegg_anno = self.load_package('annotation.kegg.kegg_annotation')()
            kegg_anno.pathSearch(blast_xml=self.gene_kegg_xml, kegg_table=self.kegg_stat_path + '/gene_kegg_table.xls')
            kegg_anno.pathTable(kegg_table=self.kegg_stat_path + '/gene_kegg_table.xls', pathway_path=self.kegg_stat_path + '/gene_pathway_table.xls', pidpath=self.work_dir + '/gene_pid.txt')
            kegg_anno.getPic(pidpath=self.work_dir + '/gene_pid.txt', pathwaydir=gene_pathway)
            kegg_anno.keggLayer(pathway_table=self.kegg_stat_path + '/gene_pathway_table.xls', layerfile=self.kegg_stat_path + '/gene_kegg_layer.xls', taxonomyfile=self.kegg_stat_path + '/gene_kegg_taxonomy.xls')
            self.logger.info('finish: kegg stat')
        except:
            import traceback
            self.logger.info('error:{}'.format(traceback.format_exc()))
            self.set_error("运行kegg脚本出错！")

    def run_go_stat(self):
        self.go_stat_path = self.work_dir + '/go_stat/'

        def get_gene_go(go_result, gene_list, outpath, trinity_mode=True):
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
        get_gene_go(go_result=self.option('blast2go_annot').prop['path'], gene_list=self.gene_list, outpath=self.go_stat_path + '/gene_blast2go.annot')
        get_gene_go(go_result=self.option('gos_list').prop['path'], gene_list=self.gene_list, outpath=self.go_stat_path + '/gene_gos.list')
        go_cmd1 = '{} {} {} {} {} {}'.format(self.python_path, self.go_annot, self.go_stat_path + '/gene_gos.list', 'localhost', self.b2g_user, self.b2g_password)
        go_cmd2 = '{} {} {}'.format(self.python_path, self.go_split, self.work_dir + '/go_detail.xls')
        go_annot_cmd = self.add_command('go_annot_cmd', go_cmd1).run()
        self.wait(go_annot_cmd)
        if go_annot_cmd.return_code == 0:
            self.logger.info("go_annot_cmd运行完成")
            self.add_command('go_split_cmd', go_cmd2).run()
        else:
            self.set_error("go_annot_cmd运行出错!")

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        try:
            self.logger.info("设置注释统计结果目录")
            self.movedir2output(self.work_dir + '/blast/', 'blast')
            for db in self.database:
                if db == 'cog':
                    self.movedir2output(self.cog_stat_path, 'cog_stat')
                    self.option('gene_string_table', self.output_dir + '/blast/gene_string.xls')
                    # venn_stat
                    self.option('string_xml').get_info()
                    string_venn = self.option('string_xml').prop['hit_query_list']
                    string_gene_venn = self.option('gene_string_table').prop['query_list']
                    self.get_venn(venn_list=string_venn, output=self.output_dir + '/venn/string_venn.txt')
                    self.get_venn(venn_list=string_gene_venn, output=self.output_dir + '/venn/gene_string_venn.txt')
                    self.anno_list['string'] = string_venn
                    self.gene_anno_list['string'] = string_gene_venn
                if db == 'nr':
                    self.movedir2output(self.nr_stat_path, 'ncbi_taxonomy')
                    self.option('gene_nr_table', self.output_dir + '/blast/gene_nr.xls')
                    self.option('nr_taxons', self.output_dir + '/ncbi_taxonomy/query_taxons.xls')
                    self.movedir2output(self.work_dir + '/blast_nr_statistics/', 'blast_nr_statistics')
                    # venn_stat
                    self.option('nr_xml').get_info()
                    nr_venn = self.option('nr_xml').prop['hit_query_list']
                    nr_gene_venn = self.option('gene_nr_table').prop['query_list']
                    self.get_venn(venn_list=nr_venn, output=self.output_dir + '/venn/nr_venn.txt')
                    self.get_venn(venn_list=nr_gene_venn, output=self.output_dir + '/venn/gene_nr_venn.txt')
                    self.anno_list['nr'] = nr_venn
                    self.gene_anno_list['nr'] = nr_gene_venn
                if db == 'kegg':
                    self.movedir2output(self.kegg_stat_path, 'kegg_stat')
                    self.option('gene_kegg_table', self.output_dir + '/blast/gene_kegg.xls')
                    self.option('gene_kegg_anno_table', self.output_dir + '/kegg_stat/gene_kegg_table.xls')
                    # venn_stat
                    self.option('kegg_xml').get_info()
                    kegg_venn = self.option('kegg_xml').prop['hit_query_list']
                    kegg_gene_venn = self.option('gene_kegg_table').prop['query_list']
                    self.get_venn(venn_list=kegg_venn, output=self.output_dir + '/venn/kegg_venn.txt')
                    self.get_venn(venn_list=kegg_gene_venn, output=self.output_dir + '/venn/gene_kegg_venn.txt')
                    self.anno_list['kegg'] = kegg_venn
                    self.gene_anno_list['kegg'] = kegg_gene_venn
                if db == 'go':
                    self.movedir2output(self.go_stat_path, 'go_stat')
                    files = os.listdir(self.work_dir)
                    self.option('gene_go_list', self.output_dir + '/go_stat/gene_gos.list')
                    for f in files:
                        if re.search(r'level_statistics\.xls$', f):
                            if os.path.exists(self.output_dir + '/go_stat/gene_go1234level_statistics.xls'):
                                os.remove(self.output_dir + '/go_stat/gene_go1234level_statistics.xls')
                            os.link(self.work_dir + '/' + f, self.output_dir + '/go_stat/gene_go1234level_statistics.xls')
                        if re.search(r'level.xls$', f):
                            if os.path.exists(self.output_dir + '/go_stat/gene_{}'.format(f)):
                                os.remove(self.output_dir + '/go_stat/gene_{}'.format(f))
                            os.link(self.work_dir + '/' + f, self.output_dir + '/go_stat/gene_{}'.format(f))
                    self.option('gene_go_level_2', self.output_dir + '/go_stat/gene_go2level.xls')
        except Exception as e:
            print traceback.format_exc()
            self.set_error("设置注释统计分析结果目录失败{}".format(e))
            self.logger.info("设置注释统计分析结果目录失败{}".format(e))

    def movedir2output(self, olddir, newname, mode='link'):
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

    def run(self):
        super(DenovoAnnoStatTool, self).run()
        for db in self.database:
            if db == 'cog':
                self.run_cog_stat()
            if db == 'nr':
                self.run_nr_stat()
            if db == 'go':
                self.run_go_stat()
            if db == 'kegg':
                self.run_kegg_stat()
        if 'go' in self.database:
            self.wait()
            self.logger.info('end: go stat')
        self.set_output()
        self.get_all_anno_stat()
        self.end()

    def get_venn(self, venn_list, output):
        with open(output, 'wb') as w:
            for i in venn_list:
                w.write(i + '\n')

    def get_all_anno_stat(self):
        # stat all_annotation_statistics.xls
        all_anno_stat = self.output_dir + '/all_annotation_statistics.xls'
        anno_num = defaultdict(dict)
        with open(all_anno_stat, 'wb') as w:
            w.write('type\ttranscripts\tgenes\ttranscripts_percent\tgenes_percent\n')
            tmp = []
            tmp_gene = []
            anno_num['total']['gene'] = len(self.option('gene_file').prop['gene_list'])
            anno_num['total']['tran'] = self.option('string_xml').prop['query_num']
            for db in self.anno_list:
                w.write('{}\t{}\t{}\t{}\t{}\n'.format(db, len(self.anno_list[db]), len(self.gene_anno_list[db]), '%0.4g' % (len(self.anno_list[db]) / anno_num['total']['tran']), '%0.4g' % (len(self.gene_anno_list[db]) / anno_num['total']['gene'])))
                tmp += self.anno_list[db]
                tmp_gene += self.gene_anno_list[db]
            anno_num['total_anno']['gene'] = len(set(tmp_gene))
            anno_num['total_anno']['tran'] = len(set(tmp))
            w.write('total_anno\t{}\t{}\t{}\t{}\n'.format(anno_num['total_anno']['tran'], anno_num['total_anno']['gene'], '%0.4g' % (anno_num['total_anno']['tran'] / anno_num['total']['tran']), '%0.4g' % (anno_num['total_anno']['gene'] / anno_num['total']['gene'])))
            w.write('total\t{}\t{}\t1\t1\n'.format(anno_num['total']['tran'], anno_num['total']['gene']))
