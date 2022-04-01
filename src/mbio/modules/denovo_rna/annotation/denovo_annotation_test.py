# -*- coding: utf-8 -*-
from __future__ import division
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.annotation.denovo_anno_stat.all_annotation_stat import *
import os
import shutil


class DenovoAnnotationTestModule(Module):
    """
    module for denovorna annotation
    last modified:20160829
    author: wangbixuan
    """
    def __init__(self, work_id):
        super(DenovoAnnotationTestModule, self).__init__(work_id)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},
            {"name": "database", "type": "string", "default": 'nr,go,cog,kegg'},  # 默认全部四个注释
            {"name": "nr_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "string_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "kegg_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "blast_threads", "type": "int", "default": 10},
            {"name": "anno_statistics", "type": "bool", "default": True},
            {"name": "gene_file", "type": "infile", "format": "denovo_rna.express.gene_list"}
        ]
        self.add_option(options)
        self.blast_nr = self.add_tool('align.ncbi.blast')
        self.blast_string = self.add_tool('align.ncbi.blast')
        self.blast_kegg = self.add_tool('align.ncbi.blast')
        self.blast_stat_nr = self.add_tool('align.ncbi.blaststat')
        self.ncbi_taxon = self.add_tool('taxon.ncbi_taxon')
        self.go_annot = self.add_tool('annotation.go_annotation')
        self.string_cog = self.add_tool('annotation.string2cog')
        self.kegg_annot = self.add_tool('annotation.kegg_annotation')
        self.anno_stat = self.add_tool('annotation.denovo_anno_stat')
        self.step.add_steps('blast_nr', 'blast_string', 'blast_kegg', 'blast_statistics', 'go_annot', 'kegg_annot', 'cog_annot', 'taxon_annot', 'anno_stat')
        self.anno_num = dict()

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query")
        else:
            if self.option('query').prop['seq_type'] != 'DNA':
                raise OptionError('传入查询序列必须是核酸序列')
        self.anno_database = self.option('database').split(',')
        if len(self.anno_database) < 1:
            raise OptionError('至少选择一种注释库')
        for i in self.anno_database:
            if i not in ['nr', 'go', 'cog', 'kegg']:
                raise OptionError('需要注释的数据库不在支持范围内[\'nr\', \'go\', \'cog\', \'kegg\']:{}'.format(i))
        if not 1 > self.option('nr_blast_evalue') >= 0 and not 1 > self.option('string_blast_evalue') >= 0 and not 1 > self.option('kegg_blast_evalue') >= 0:
            raise OptionError('E-value值设定必须为[0-1)之间')
        if self.option('anno_statistics') and not self.option('gene_file').is_set:
            raise OptionError('运行注释统计的tool必须要设置gene_file')

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def mark_end(self, event):
        with open(self.work_dir + '/' + event['data'] + '.o', 'wb') as w:
            w.write(event['data'] + 'finish')

    def test_end(self, tp):
        if not os.listdir(tp.output_dir):
            tp.logger.debug("输出目录%s为空,你确定已经设置了输出目录?" % tp.output_dir)
        for option in tp._options.values():
            if option.type == 'outfile':
                if not option.value.is_set:
                    tp.logger.debug("输出参数%s没有设置输出文件路径,你确定此处不需要设置?" % option.name)
        tp.set_end()
        tp.fire('end')

    def test_run(self, tp):
        tp.start_listener()
        paused = False
        workflow = tp.get_workflow()
        while workflow.pause:
            if not paused:
                tp.logger.info("流程处于暂停状态，排队等待恢复运行!")
            paused = True
            workflow.is_wait = True
            gevent.sleep(1)
        tp.fire("start")

    def run_blast(self):
        """
        """
        print self.work_dir
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
            self.blast_nr.on('end', self.mark_end, 'nrblast')
            if os.path.exists(self.work_dir + '/nrblast.o'):
                print '............not run nr '
                self.blast_nr.option('outxml', self.blast_nr.work_dir + '/samll_vs_nr.xml')
                self.blast_nr.option('outtable', self.blast_nr.output_dir + '/samll_vs_nr.xls')
                # self.blast_nr._start = True
                # self.blast_nr.fire('end')
                self.test_run(self.blast_nr)
                self.test_end(self.blast_nr)
                # self.run()
            else:
                print '............run nr '
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
            self.blast_string.on('end', self.mark_end, 'stringblast')
            if os.path.exists(self.work_dir + '/stringblast.o'):
                print '............not run string '
                self.blast_string.option('outxml', self.blast_string.work_dir + '/samll_vs_string.xml')
                print '........%s' % self.blast_string.option('outxml').prop['path']
                self.blast_string.option('outtable', self.blast_string.output_dir + '/samll_vs_string.xls')
                # self.blast_string._start = True
                # self.blast_string.fire('end')
                self.test_run(self.blast_string)
                self.test_end(self.blast_string)
            else:
                print '............ run string '
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
            self.blast_kegg.on('end', self.mark_end, 'keggblast')
            if os.path.exists(self.work_dir + '/keggblast.o'):
                self.blast_kegg.option('outxml', self.blast_kegg.work_dir + '/samll_vs_kegg.xml')
                self.blast_kegg.option('outtable', self.blast_kegg.output_dir + '/samll_vs_kegg.xls')
                # self.blast_kegg._start = True
                # self.blast_kegg.fire('end')
                self.test_run(self.blast_kegg)
                self.test_end(self.blast_kegg)
            else:
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
        opts = {'gene_file': self.option('gene_file')}
        if 'kegg' in self.anno_database:
            opts['kegg_xml'] = self.blast_kegg.option('outxml')
        if 'go' in self.anno_database:
            opts['gos_list'] = self.go_annot.option('golist_out')
            opts['blast2go_annot'] = self.go_annot.option('blast2go_annot')
        if 'cog' in self.anno_database:
            opts['string_xml'] = self.blast_string.option('outxml')
            opts['cog_list'] = self.string_cog.option('cog_list')
            opts['cog_table'] = self.string_cog.option('cog_table')
        if 'nr' in self.anno_database:
            opts['nr_xml'] = self.blast_kegg.option('outxml')
            opts['nr_taxon_details'] = self.ncbi_taxon.option('taxon_out')
        self.anno_stat.set_options(opts)
        # for k in opts:
        #     # print i._name, i.prop['path']
        #     print k, opts[k]
        # for i in opts:
        #     # print i, self.anno_stat._options[i]
        #     print self.anno_stat.option(i).prop['path']
        self.anno_stat.on('start', self.set_step, {'start': self.step.anno_stat})
        self.anno_stat.on('end', self.set_step, {'end': self.step.anno_stat})
        self.anno_stat.on('end', self.set_output, 'anno_stat')
        self.anno_stat.run()

    def run_kegg_anno(self):
        """
        """
        self.anno_num['kegg'] = [self.blast_kegg.option('outtable').prop['query_num']]
        options = {
            'blastout': self.blast_kegg.option('outxml')
        }
        self.kegg_annot.set_options(options)
        self.kegg_annot.on('start', self.set_step, {'start': self.step.kegg_annot})
        self.kegg_annot.on('end', self.set_step, {'end': self.step.kegg_annot})
        self.kegg_annot.on('end', self.set_output, 'kegg_annot')
        self.kegg_annot.on('end', self.mark_end, 'kegganno')
        if os.path.exists(self.work_dir + '/kegganno.o'):
            self.kegg_annot.option('kegg_table', self.kegg_annot.output_dir + '/kegg_table.xls')
            self.test_run(self.kegg_annot)
            self.test_end(self.kegg_annot)
        else:
            self.kegg_annot.run()

    def run_string2cog(self):
        self.anno_num['string'] = [self.blast_string.option('outtable').prop['query_num']]
        options = {
            'blastout': self.blast_string.option('outxml')
        }
        self.string_cog.set_options(options)
        self.string_cog.on('start', self.set_step, {'start': self.step.cog_annot})
        self.string_cog.on('end', self.set_step, {'end': self.step.cog_annot})
        self.string_cog.on('end', self.set_output, 'string_cog')
        self.string_cog.on('end', self.mark_end, 'coganno')
        if os.path.exists(self.work_dir + '/coganno.o'):
            self.string_cog.option('cog_list', self.string_cog.output_dir + '/cog_list.xls')
            print self.string_cog.option('cog_list').prop['path']
            self.string_cog.option('cog_table', self.string_cog.output_dir + '/cog_table.xls')
            self.test_run(self.string_cog)
            self.test_end(self.string_cog)
        else:
            self.string_cog.run()

    def run_go_anno(self):
        """
        """
        options = {
            'blastout': self.blast_nr.option('outxml')
        }
        self.go_annot.set_options(options)
        self.go_annot.on('start', self.set_step, {'start': self.step.go_annot})
        self.go_annot.on('end', self.set_step, {'end': self.step.go_annot})
        self.go_annot.on('end', self.set_output, 'go_annot')
        self.go_annot.on('end', self.mark_end, 'goanno')
        if os.path.exists(self.work_dir + '/goanno.o'):
            self.go_annot.option('go2level_out', self.go_annot.output_dir + '/go2level.xls')
            self.go_annot.option('golist_out', self.go_annot.output_dir + '/query_gos.list')
            self.go_annot.option('blast2go_annot', self.go_annot.work_dir + '/blast2go.annot')
            self.test_run(self.go_annot)
            self.test_end(self.go_annot)
        else:
            self.go_annot.run()

    def run_blast_stat(self):
        """
        nr库比对结果统计函数
        """
        options = {
            'in_stat': self.blast_nr.option('outxml')
        }
        self.blast_stat_nr.set_options(options)
        self.blast_stat_nr.on('start', self.set_step, {'start': self.step.blast_statistics})
        self.blast_stat_nr.on('end', self.set_step, {'end': self.step.blast_statistics})
        self.blast_stat_nr.on('end', self.set_output, 'blast_stat')
        self.blast_stat_nr.run()

    def run_ncbi_taxon(self):
        """
        """
        self.anno_num['nr'] = [self.blast_nr.option('outtable').prop['query_num']]
        options = {
            'blastout': self.blast_nr.option('outxml'),
            'blastdb': 'nr'
        }
        self.ncbi_taxon.set_options(options)
        self.ncbi_taxon.on('start', self.set_step, {'start': self.step.taxon_annot})
        self.ncbi_taxon.on('end', self.set_step, {'end': self.step.taxon_annot})
        self.ncbi_taxon.on('end', self.set_output, 'ncbi_taxon')
        self.ncbi_taxon.on('end', self.mark_end, 'taxon')
        if os.path.exists(self.work_dir + '/taxon.o'):
            self.ncbi_taxon.option('taxon_out', self.ncbi_taxon.output_dir + '/query_taxons_detail.xls')
            self.test_run(self.ncbi_taxon)
            self.test_end(self.ncbi_taxon)
        else:
            self.ncbi_taxon.run()
        # self.ncbi_taxon.run()

    def run(self):
        super(DenovoAnnotationTestModule, self).run()
        self.run_blast()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'nrblast':
            self.linkdir(obj.output_dir, 'nrblast')
        elif event['data'] == 'stringblast':
            self.linkdir(obj.output_dir, 'stringblast')
        elif event['data'] == 'keggblast':
            self.linkdir(obj.output_dir, 'keggblast')
        elif event['data'] == 'blast_stat':
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
            print '........%s' % self.anno_num
            print '............%s' % self.option('query').prop['seq_number']
            self.anno_num['total'] = [float(self.option('query').prop['seq_number']), float(self.option('gene_file').prop['gene_num'])]
            if 'nr' in self.anno_database:
                self.anno_num['nr'].append(obj.option('gene_nr_table').prop['query_num'])
            if 'kegg' in self.anno_database:
                self.anno_num['kegg'].append(obj.option('gene_kegg_table').prop['query_num'])
            if 'cog' in self.anno_database:
                self.anno_num['string'].append(obj.option('gene_string_table').prop['query_num'])
            try:
                self.get_all_anno_stat(self.output_dir + '/anno_stat/all_annotation_statistics.xls', self.output_dir + '/anno_stat/all_annotation.xls')
            except Exception as e:
                self.logger.info("统计all_annotation出错：{}".format(e))
            self.end()
        else:
            pass

    def get_all_anno_stat(self, all_anno_stat_path, all_anno_path):
        # stat all_annotation_statistics.xls
        with open(all_anno_stat_path, 'wb') as w:
            w.write('type\ttranscripts\tgenes\ttranscripts_percent\tgenes_percent\n')
            anno_sum = [0, 0]
            dbs = self.anno_database
            if 'cog' in dbs:
                dbs.remove('cog')
                dbs.append('string')
            if 'go' in dbs:
                dbs.remove('go')
            for k in self.anno_num:
                self.anno_num[k] = [float(i) for i in self.anno_num[k]]
            for db in dbs:
                anno_sum[0] += self.anno_num[db][0]
                anno_sum[1] += self.anno_num[db][1]
                w.write('{}\t{}\t{}\t{}\t{}\n'.format(db, self.anno_num[db][0], self.anno_num[db][1], '%0.4g' % (self.anno_num[db][0] / self.anno_num['total'][0]), '%0.4g' % (self.anno_num[db][1] / self.anno_num['total'][1])))
            w.write('total_anno\t{}\t{}\t{}\t{}\n'.format(anno_sum[0], anno_sum[1], '%0.4g' % (anno_sum[0] / self.anno_num['total'][0]), '%0.4g' % (anno_sum[1] / self.anno_num['total'][1])))
            w.write('total\t{}\t{}\t1\t1\n'.format(self.anno_num['total'][0], self.anno_num['total'][1]))
            print 'aaaaaaaaaaa'
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
                kwargs['blast_nr_table'] = self.blast_nr.option('outtable').prop['path']
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
            ["/anno_stat/nr_stat/", "dir", "nr统计结果目录"],
            ["/anno_stat/cog_stat/", "dir", "cog统计结果目录"],
            ["/anno_stat/go_stat/", "dir", "go统计结果目录"],
            ["/anno_stat/kegg_stat/", "dir", "kegg统计结果目录"],
            ["/anno_stat/blast_result/", "dir", "基因序列blast比对结果目录"],
            ["/anno_stat/blast_result/gene_kegg.xls", "xls", "基因序列blast比对kegg注释结果table"],
            ["/anno_stat/blast_result/gene_nr.xls", "xls", "基因序列blast比对nr注释结果table"],
            ["/anno_stat/blast_result/gene_string.xls", "xls", "基因序列blast比对string注释结果table"],
            ["/anno_stat/blast_result/gene_string.xml", "xml", "基因序列blast比对string注释结果xml"],
            ["/anno_stat/blast_result/gene_kegg.xml", "xml", "基因序列blast比对kegg注释结果xml"],
            ["/anno_stat/blast_result/gene_string.xml", "xml", "基因序列blast比对string注释结果xml"],
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
            ['/nr_stat/gene_taxons_detail.xls', 'xls', '基因序列详细物种分类文件'],
            ["/anno_stat/nr_stat/gene_nr_evalue.xls", "xls", "基因序列blast结果E-value统计"],
            ["/anno_stat/nr_stat/gene_nr_similar.xls", "xls", "基因序列blast结果similarity统计"],
            ["/anno_stat/nr_stat/nr_taxon_stat.xls", "xls", "nr物种分类统计表"],
            ["/anno_stat/nr_stat/gene_taxons.xls", "xls", "基因序列nr物种注释表"],
            ["/anno_stat/nr_stat/query_taxons.xls", "xls", "nr物种注释表"],
        ]
        regexps = [
            [r"blast/.+_vs_.+\.xml", "xml", "blast比对输出结果，xml格式"],
            [r"blast/.+_vs_.+\.xls", "xls", "blast比对输出结果，表格(制表符分隔)格式"],
            [r"blast/.+_vs_.+\.txt", "txt", "blast比对输出结果，非xml和表格(制表符分隔)格式"],
            [r"blast/.+_vs_.+\.txt_\d+\.xml", "xml",
                "Blast比对输出多xml结果，输出格式为14的单个比对结果文件,主结果文件在txt文件中"],
            [r"blast/.+_vs_.+\.txt_\d+\.json", "json",
                "Blast比输出对多json结果，输出格式为13的单个比对结果文件,主结果文件在txt文件中"],
            [r"kegg/pathways/ko.\d+", 'pdf', '标红pathway图'],
            [r"/blast_nr_statistics/.*_evalue\.xls", "xls", "比对结果E-value分布图"],
            [r"/blast_nr_statistics/.*_similar\.xls", "xls", "比对结果相似度分布图"]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(DenovoAnnotationTestModule, self).end()
