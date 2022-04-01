# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

"""无参转录组基础分析"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import os
import shutil
import re
import gevent


class TestBaseWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        version = v1.0
        last_modify = 20160825
        """
        self._sheet = wsheet_object
        super(TestBaseWorkflow, self).__init__(wsheet_object)
        print self._parent
        options = [
            {"name": "fastq_dir", "type": "infile", 'format': "sequence.fastq,sequence.fastq_dir"},  # fastq文件夹
            {"name": "fq_type", "type": "string"},  # PE OR SE
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 有生物学重复的时候的分组文件
            {"name": "control_file", "type": "infile", "format": "denovo_rna.express.control_table"},  # 对照组文件，格式同分组文件
            {"name": "search_pfam", "type": "bool", "default": True},  # orf 是否比对Pfam数据库
            {"name": "primer", "type": "bool", "default": True},  # 是否设计SSR引物
            {"name": "kmer_size", "type": "int", "default": 25},
            {"name": "min_kmer_cov", "type": "int", "default": 2},
            {"name": "min_contig_length", "type": "int", "default": 200},  # trinity报告出的最短的contig长度。默认为200
            {"name": "SS_lib_type", "type": "string", "default": 'none'},  # reads的方向，成对的reads: RF or FR; 不成对的reads: F or R，默认情况下，不设置此参数
            {"name": "exp_way", "type": "string", "default": "fpkm"},  # edger离散值
            {"name": "diff_ci", "type": "float", "default": 0.01},  # 显著性水平
            {"name": "diff_rate", "type": "float", "default": 0.01},  # 期望的差异基因比率
            {"name": "database", "type": "string", "default": 'nr,go,cog,kegg'},  # 默认全部四个注释
            {"name": "nr_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "string_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "kegg_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "exp_analysis", "type": "string", "default": "cluster,network,kegg_rich,go_rich"},
            {"name": "gene_analysis", "type": "string", "default": "orf,snp"},
            {"name": "map_qc_analysis", "type": "string", "default": "saturation,duplication,stat,correlation"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.filecheck = self.add_tool("denovo_rna.filecheck.file_denovo")
        self.qc = self.add_module("denovo_rna.qc.quality_control")
        self.qc_stat_before = self.add_module("denovo_rna.qc.qc_stat")
        self.qc_stat_after = self.add_module("denovo_rna.qc.qc_stat")
        self.assemble = self.add_tool("denovo_rna.assemble.assemble")
        self.annotation = self.add_module('denovo_rna.annotation.denovo_annotation')
        self.orf = self.add_tool("denovo_rna.gene_structure.orf")
        self.ssr = self.add_tool("denovo_rna.gene_structure.ssr")
        self.bwa = self.add_module("denovo_rna.mapping.bwa_samtools")
        self.snp = self.add_module("denovo_rna.gene_structure.snp")
        self.map_qc = self.add_module("denovo_rna.mapping.map_assessment")
        self.exp_stat = self.add_module("denovo_rna.express.exp_analysis")
        self.exp_diff = self.add_module("denovo_rna.express.diff_analysis")
        self.orf_len = self.add_tool("meta.qc.reads_len_info")
        self.step.add_steps("qcstat", "assemble", "annotation", "express", "gene_structure", "map_stat")
        self.final_tools = list()
        self.update_status_api = self.api.denovo_update_status
        self.api_sample = self.api.denovo_rna_sample
        self.spname_spid = dict()
        self.samples = self.option('fastq_dir').samples
        self.diff_gene_id = None  # 差异基因矩阵id
        self.diff_genes = None  # 差异基因列表
        self.bam_path = ''  # rsem比对结果bam文件文件路径
        self.orf_bed = ''  # orf bed结果文件路径
        self.express_id = None
        self.api_express = self.api.denovo_express
        # self.databse = self.option('database').split(',')
        # self.option('gene_analysis') = self.option('gene_analysis').split(',')
        # self.option('exp_analysis') = self.option('exp_analysis').split(',')
        # self.map_qc_analysis = self.option('map_qc_analysis').split(',')

    def check_options(self):
        """
        检查参数设置
        """
        if not self.option('fastq_dir').is_set:
            raise OptionError('必须设置输入fastq文件夹')
        if not self.option('control_file').is_set:
            raise OptionError('必须设置输入对照方案文件')
        if not self.option('fq_type'):
            raise OptionError('必须设置测序类型：PE OR SE')
        if self.option('fq_type') not in ['PE', 'SE']:
            raise OptionError('测序类型不在所给范围内')
        if self.option('diff_ci') > 1 or self.option('diff_ci') < 0:
            raise OptionError('显著性水平不在所给范围内[0,1]')
        if self.option('diff_rate') > 1 or self.option('diff_rate') < 0:
            raise OptionError('差异基因比率不在所给范围内[0,1]')
        if self.option("fq_type") == 'SE' and self.option("SS_lib_type") not in ['F', 'R', 'none']:
            raise OptionError("SE测序时所设reads方向：{}不正确".format(self.option("SS_lib_type")))
        if self.option("fq_type") == 'PE' and self.option("SS_lib_type") not in ['FR', 'RF', 'none']:
            raise OptionError("PE测序时所设reads方向：{}不正确".format(self.option("SS_lib_type")))
        if self.option("exp_way") not in ['fpkm', 'tpm']:
            raise OptionError("所设表达量的代表指标不在范围内，请检查")
        if self.option('kmer_size') > 32 or self.option('kmer_size') < 1:
            raise OptionError("所设kmer_size不在范围内，请检查")
        if self.option('min_kmer_cov') < 1:
            raise OptionError("所设min_kmer_cov不在范围内，请检查")
        for i in self.option('exp_analysis').split(','):
            if i not in ['', 'cluster', 'network', 'kegg_rich', 'go_rich']:
                raise OptionError("差异性研究没有{}，请检查".format(i))
        for i in self.option('gene_analysis').split(','):
            if i not in ['orf', 'ssr', 'snp']:
                raise OptionError("基因结构分析没有{}，请检查".format(i))
        for i in self.option('map_qc_analysis').split(','):
            if i not in ['', 'saturation', 'coverage', 'duplication', 'correlation', 'stat']:
                raise OptionError("转录组质量评估没有{}，请检查".format(i))
        self.anno_database = self.option('database').split(',')
        if len(self.anno_database) < 1:
            raise OptionError('至少选择一种注释库')
        for i in self.anno_database:
            if i not in ['nr', 'go', 'cog', 'kegg']:
                raise OptionError('需要注释的数据库不在支持范围内[\'nr\', \'go\', \'cog\', \'kegg\']:{}'.format(i))
        if not 1 > self.option('nr_blast_evalue') >= 0 and not 1 > self.option('string_blast_evalue') >= 0 and not 1 > self.option('kegg_blast_evalue') >= 0:
            raise OptionError('E-value值设定必须为[0-1)之间')

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run_filecheck(self):
        opts = {
            'fastq_dir': self.option('fastq_dir'),
            'fq_type': self.option('fq_type'),
            'control_file': self.option('control_file')
        }
        if self.option('group_table').is_set:
            opts.update({'group_table': self.option('group_table')})
        self.filecheck.set_options(opts)
        # self.filecheck.run()
        self.filecheck.on('end', self.mark_end, 'filecheck')
        if os.path.exists(self.work_dir + '/filecheck.o'):
            print '............filecheck have already run'
            self.test_run(self.filecheck)
            self.test_end(self.filecheck)
        else:
            print '............run filecheck '
            self.filecheck.run()

    def run_qc(self):
        self.qc.set_options({
            'fastq_dir': self.option('fastq_dir'),
            'fq_type': self.option('fq_type')
        })
        self.qc.on('end', self.set_output, 'qc')
        self.qc.on('start', self.set_step, {'start': self.step.qcstat})
        self.qc.on('end', self.set_step, {'end': self.step.qcstat})
        # self.qc.run()
        self.qc.on('end', self.mark_end, 'qc')
        if os.path.exists(self.work_dir + '/qc.o'):
            print '............qc have already run'
            self.test_run(self.qc)
            self.qc.option('sickle_dir', self.qc.output_dir + '/sickle_dir')
            if self.option('fq_type') == 'SE':
                self.qc.option('fq_s', self.qc.work_dir + '/single.fq')
            else:
                self.qc.option('fq_l', self.qc.work_dir + '/left.fq')
                self.qc.option('fq_r', self.qc.work_dir + '/right.fq')
            self.test_end(self.qc)
        else:
            print '............run qc '
            self.qc.run()

    def run_qc_stat(self, event):
        if event['data']:
            print '........start run qc_stat_after'
            self.qc_stat_after.set_options({
                'fastq_dir': self.qc.option('sickle_dir'),
                'fq_type': self.option('fq_type'),
                'dup': True
            })
        else:
            self.qc_stat_before.set_options({
                'fastq_dir': self.option('fastq_dir'),
                'fq_type': self.option('fq_type')})
        if event['data']:
            self.qc_stat_after.on('end', self.set_output, 'qc_stat_after')
            # self.qc_stat_after.run()
            self.qc_stat_after.on('end', self.mark_end, 'qc_stat_after')
            if os.path.exists(self.work_dir + '/qc_stat_after.o'):
                print '............qc_stat_after have already run'
                self.test_run(self.qc_stat_after)
                self.test_end(self.qc_stat_after)
            else:
                print '............run qc_stat_after '
                self.qc_stat_after.run()
        else:
            self.qc_stat_before.on('end', self.set_output, 'qc_stat_before')
            # self.qc_stat_before.run()
            self.qc_stat_before.on('end', self.mark_end, 'qc_stat_before')
            if os.path.exists(self.work_dir + '/qc_stat_before.o'):
                print '............qc_stat_before have already run'
                self.test_run(self.qc_stat_before)
                self.test_end(self.qc_stat_before)
            else:
                print '............run qc_stat_before '
                self.qc_stat_before.run()

    def run_assemble(self):
        print '........start run run_assemble'
        opts = {
            'fq_type': self.option('fq_type'),
            'min_contig_length': self.option('min_contig_length'),
            'SS_lib_type': self.option('SS_lib_type'),
            'kmer_size': self.option('kmer_size'),
            'min_kmer_cov': self.option('min_kmer_cov'),
        }
        if self.option('fq_type') == 'SE':
            opts.update({'fq_s': self.qc.option('fq_s')})
        else:
            opts.update({
                        'fq_l': self.qc.option('fq_l'),
                        'fq_r': self.qc.option('fq_r')
                        })
        self.assemble.set_options(opts)
        self.assemble.on('end', self.set_output, 'assemble')
        self.assemble.on('start', self.set_step, {'start': self.step.assemble})
        self.assemble.on('end', self.set_step, {'end': self.step.assemble})
        # self.assemble.run()
        self.assemble.on('end', self.mark_end, 'assemble')
        if os.path.exists(self.work_dir + '/assemble.o'):
            print '............assemble have already run'
            self.test_run(self.assemble)
            self.assemble.option('gene_fa', self.assemble.work_dir + '/gene.fasta')
            self.assemble.option('trinity_fa', self.assemble.work_dir + '/trinity_out_dir/Trinity.fasta')
            self.assemble.option('gene_full_name', self.assemble.work_dir + '/gene_full_name.txt')
            self.test_end(self.assemble)
        else:
            print '............run assemble '
            self.assemble.run()

    def run_orf(self):
        orf_opts = {
            'fasta': self.assemble.option('trinity_fa'),
            'search_pfam': self.option('search_pfam')
        }
        self.orf.set_options(orf_opts)
        self.orf.on('end', self.set_output, 'orf')
        self.orf.on('start', self.set_step, {'start': self.step.gene_structure})
        # self.orf.run()
        self.orf.on('end', self.mark_end, 'orf')
        if os.path.exists(self.work_dir + '/orf.o'):
            print '............orf have already run'
            self.test_run(self.orf)
            self.orf.option('bed', self.orf.output_dir + '/Trinity.fasta.transdecoder.bed')
            self.orf.option('cds', self.orf.output_dir + '/Trinity.fasta.transdecoder.cds')
            self.orf.option('pep', self.orf.output_dir + '/Trinity.fasta.transdecoder.pep')
            self.test_end(self.orf)
        else:
            print '............run orf '
            self.orf.run()

    def run_orf_len(self):
        orf_fasta = self.orf.work_dir + '/ORF_fasta'
        self.orf_len.set_options({'fasta_path': orf_fasta})
        self.orf_len.on('end', self.set_output, 'orf_len')
        # self.orf_len.run()
        self.orf_len.on('end', self.mark_end, 'orf_len')
        if os.path.exists(self.work_dir + '/orf_len.o'):
            print '............orf_len have already run'
            self.test_run(self.orf_len)
            self.test_end(self.orf_len)
        else:
            print '............run orf_len '
            self.orf_len.run()

    def run_bwa(self):
        bwa_opts = {
            'ref_fasta': self.assemble.option('gene_fa'),
            'fq_type': self.option('fq_type'),
            'fastq_dir': self.qc.option('sickle_dir'),
        }
        self.bwa.set_options(bwa_opts)
        # self.bwa.run()
        self.bwa.on('end', self.mark_end, 'bwa')
        if os.path.exists(self.work_dir + '/bwa.o'):
            print '............bwa have already run'
            self.test_run(self.bwa)
            self.bwa.option('out_bam', self.bwa.output_dir + '/sorted_bam')
            self.test_end(self.bwa)
        else:
            print '............run bwa '
            self.bwa.run()

    def run_snp(self):
        snp_opts = {
            'bed': self.orf.option('bed'),
            'bam': self.bwa.option('out_bam'),
            'ref_fasta': self.assemble.option('gene_fa')
        }
        self.snp.set_options(snp_opts)
        self.snp.on('end', self.set_output, 'snp')
        # self.snp.run()
        self.snp.on('end', self.mark_end, 'snp')
        if os.path.exists(self.work_dir + '/snp.o'):
            print '............snp have already run'
            self.test_run(self.snp)
            self.test_end(self.snp)
        else:
            print '............run snp '
            self.snp.run()

    def run_ssr(self):
        ssr_opts = {
            'fasta': self.assemble.option('gene_fa'),
            'bed': self.orf.option('bed'),
            'primer': self.option('primer')
        }
        self.ssr.set_options(ssr_opts)
        self.ssr.on('end', self.set_output, 'ssr')
        # self.ssr.run()
        self.ssr.on('end', self.mark_end, 'ssr')
        if os.path.exists(self.work_dir + '/ssr.o'):
            print '............ssr have already run'
            self.test_run(self.ssr)
            self.test_end(self.ssr)
        else:
            print '............run ssr '
            self.ssr.run()

    def run_map_qc(self):
        map_qc_opts = {
            'bed': self.orf.option('bed'),
            'bam': self.exp_stat.option('bam_dir'),
            'fpkm': self.exp_stat.option('gene_fpkm'),
            'analysis': self.option('map_qc_analysis')
        }
        self.map_qc.set_options(map_qc_opts)
        self.map_qc.on('end', self.set_output, 'map_qc')
        self.map_qc.on('end', self.set_step, {'end': self.step.map_stat})
        # self.map_qc.run()
        self.map_qc.on('end', self.mark_end, 'map_qc')
        if os.path.exists(self.work_dir + '/map_qc.o'):
            print '............map_qc have already run'
            self.test_run(self.map_qc)
            self.test_end(self.map_qc)
        else:
            print '............run map_qc '
            self.map_qc.run()

    def run_annotation(self):
        anno_opts = {
            "query": self.assemble.option('trinity_fa'),
            "database": self.option('database'),
            "nr_blast_evalue": self.option('nr_blast_evalue'),
            "string_blast_evalue": self.option('string_blast_evalue'),
            "kegg_blast_evalue": self.option('kegg_blast_evalue'),
            "gene_file": self.assemble.option('gene_full_name'),
        }
        self.annotation.set_options(anno_opts)
        self.annotation.on('end', self.set_output, 'annotation')
        self.annotation.on('start', self.set_step, {'start': self.step.annotation})
        self.annotation.on('end', self.set_step, {'end': self.step.annotation})
        # self.annotation.run()
        self.annotation.on('end', self.mark_end, 'annotation')
        if os.path.exists(self.work_dir + '/annotation.o'):
            print '............annotation have already run'
            self.test_run(self.annotation)
            self.annotation.option('gene_go_list', self.annotation.work_dir + '/DenovoAnnoStat/go_stat/gene_gos.list')
            self.annotation.option('gene_go_level_2', self.annotation.work_dir + '/DenovoAnnoStat/go2level.xls')
            self.annotation.option('gene_kegg_table', self.annotation.work_dir + '/DenovoAnnoStat/kegg_stat/gene_kegg_table.xls')
            self.test_end(self.annotation)
        else:
            print '............run annotation '
            self.annotation.run()

    def run_exp_stat(self):
        exp_stat_opts = {
            'fq_type': self.option('fq_type'),
            'rsem_fa': self.assemble.option('trinity_fa'),
            'control_file': self.option('control_file'),
            'exp_way': self.option('exp_way'),
            'diff_ci': self.option('diff_ci'),
            'diff_rate': self.option('diff_rate')
        }
        if self.option('fq_type') == 'SE':
            exp_stat_opts.update({'fq_s': self.qc.option('sickle_dir')})
        else:
            exp_stat_opts.update({'fq_r': self.qc.option('sickle_r_dir'), 'fq_l': self.qc.option('sickle_l_dir')})
        if self.option('group_table').is_set:
            exp_stat_opts.update({
                'group_table': self.option('group_table'),
                'gname': self.option('group_table').prop['group_scheme'][0]
            })
        self.exp_stat.set_options(exp_stat_opts)
        self.exp_stat.on('end', self.set_output, 'exp_stat')
        self.exp_stat.on('start', self.set_step, {'start': self.step.express})
        # self.exp_stat.run()
        self.exp_stat.on('end', self.mark_end, 'exp_stat')
        if os.path.exists(self.work_dir + '/exp_stat.o'):
            print '............exp_stat have already run'
            self.test_run(self.exp_stat)
            self.exp_stat.option('bam_dir', self.exp_stat.work_dir + '/bowtie2_bam_dir/')
            # self.exp_stat.option('tran_count', )
            self.exp_stat.option('gene_count', self.exp_stat.merge_rsem.option('gene_count'))
            # self.exp_stat.option('tran_fpkm', )
            self.exp_stat.option('gene_fpkm', self.exp_stat.merge_rsem.option('gene_fpkm'))
            self.exp_stat.option('all_list', self.exp_stat.work_dir + '/all_list')
            self.exp_stat.diff_gene = True
            if self.exp_stat.diff_gene:
                # print self.exp_stat.diff_exp.option('diff_fpkm').prop['path']
                # print self.exp_stat.diff_exp.option('diff_fpkm').prop['path']
                self.exp_stat.option('all_list', self.exp_stat.work_dir + '/DiffExp/gene_file')
                self.exp_stat.option('diff_count', self.exp_stat.output_dir + '/diff_exp/diff_count')
                self.exp_stat.option('diff_fpkm', self.exp_stat.output_dir + '/diff_exp/diff_fpkm')
                self.exp_stat.option('diff_list_dir', self.exp_stat.work_dir + '/DiffExp/diff_list_dir')
                self.exp_stat.option('diff_list', self.exp_stat.work_dir + '/DiffExp/diff_list')
            self.test_end(self.exp_stat)
        else:
            print '............run exp_stat '
            self.exp_stat.run()

    def run_exp_diff(self):
        print '............%s' % self.exp_stat.option('diff_list_dir').prop['path']
        if self.exp_stat.diff_gene:
            exp_diff_opts = {
                'diff_fpkm': self.exp_stat.option('diff_fpkm'),
                'analysis': self.option('exp_analysis')
            }
            print '.............%sexp_diff_opts:' % exp_diff_opts
            print '.............%s' % self.option('exp_analysis')
            if 'network' in self.option('exp_analysis'):
                print 'aaaaa'
                exp_diff_opts.update({'diff_list': self.exp_stat.option('diff_list')})
            if 'kegg_rich' in self.option('exp_analysis'):
                print 'bbbb'
                exp_diff_opts.update({
                    'gene_kegg_table': self.annotation.option('gene_kegg_table'),
                    'diff_list_dir': self.exp_stat.option('diff_list_dir'),
                    'all_list': self.exp_stat.option('all_list'),
                })
            if 'go_rich' in self.option('exp_analysis'):
                print 'cccc'
                exp_diff_opts.update({
                    'gene_go_list': self.annotation.option('gene_go_list'),
                    'diff_list_dir': self.exp_stat.option('diff_list_dir'),
                    'all_list': self.exp_stat.option('all_list'),
                    'gene_go_level_2': self.annotation.option('gene_go_level_2'),
                    'diff_stat_dir': self.exp_stat.diff_exp.option('regulate_edgrstat_dir')
                })
            print '.............%sexp_diff_opts:' % exp_diff_opts
            print '.......%sexp_diff_opts:' % exp_diff_opts['diff_stat_dir']
            self.exp_diff.set_options(exp_diff_opts)
            self.exp_diff.on('end', self.set_output, 'exp_diff')
            self.exp_diff.on('end', self.set_step, {'end': self.step.express})
            self.exp_diff.run()
            self.final_tools.append(self.exp_diff)
        else:
            self.logger.info('输入文件数据量过小，没有检测到差异基因，差异基因相关分析将忽略')
        self.on_rely(self.final_tools, self.end)


    def move2outputdir(self, olddir, newname, mode='link'):
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
            self.logger.info(newfiles)
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

    def set_output(self, event):
        obj = event["bind_object"]
        # 设置qc报告文件
        if event['data'] == 'qc':
            self.move2outputdir(obj.output_dir, 'QC_stat')
        if event['data'] == 'qc_stat_before':
            self.move2outputdir(obj.output_dir, 'QC_stat/before_qc')
            # set api
            qc_stat_info = self.output_dir + '/QC_stat/before_qc'
            quality_stat = self.output_dir + '/QC_stat/before_qc/qualityStat'
            if not os.path.exists(qc_stat_info):
                raise Exception('找不到报告文件：{}'.format(qc_stat_info))
            if not os.path.exists(quality_stat):
                raise Exception('找不到报告文件：{}'.format(quality_stat))
            api_sample = self.api.denovo_rna_sample
            api_sample.add_samples_info(qc_stat=qc_stat_info, qc_adapt=None, fq_type=self.option('fq_type'))
            api_sample.add_gragh_info(quality_stat, about_qc='before')
        if event['data'] == 'qc_stat_after':
            self.move2outputdir(obj.output_dir, 'QC_stat/after_qc')
            # set api
            qc_stat_info = self.output_dir + '/QC_stat/after_qc'
            quality_stat = self.output_dir + '/QC_stat/after_qc/qualityStat'
            qc_adapt = self.output_dir + '/QC_stat/adapter.xls'
            files = [qc_stat_info, quality_stat, qc_adapt]
            for f in files:
                if not os.path.exists(f):
                    raise Exception('找不到报告文件：{}'.format(f))
            self.api_sample.add_samples_info(qc_stat=qc_stat_info, qc_adapt=qc_adapt, fq_type=self.option('fq_type'))
            self.api_sample.add_gragh_info(quality_stat, about_qc='after')
            self.spname_spid = self.api_sample.get_spname_spid()
        if event['data'] == 'assemble':
            self.move2outputdir(obj.output_dir, 'Assemble')
            # set api
            api_assemble = self.api.denovo_assemble
            trinity_path = obj.work_dir + '/trinity_out_dir/Trinity.fasta'
            gene_path = obj.work_dir + '/gene.fasta'
            stat_path = self.output_dir + '/Assemble/trinity.fasta.stat.xls'
            files = [trinity_path, gene_path, stat_path]
            for f in files:
                if not os.path.exists(f):
                    raise Exception('找不到报告文件：{}'.format(f))
            sequence_id = api_assemble.add_sequence(trinity_path, gene_path)
            api_assemble.add_sequence_detail(sequence_id, stat_path)
            for f in os.listdir(self.output_dir + '/Assemble'):
                if re.search(r'length\.distribut\.txt$', f):
                    step = f.split('_')[0]
                    file_ = self.output_dir + '/Assemble/' + f
                    api_assemble.add_sequence_step(sequence_id, file_, step)
        if event['data'] == 'map_qc':
            self.move2outputdir(obj.output_dir, 'Map_stat')
            api_map = self.api.denovo_rna_mapping
            api_map.add_mapping_stat(self.output_dir + '/Map_stat/bam_stat.xls')
            rpkm_params = {
                'analysis': 'satur',
                'quality_satur': 30,
                'low_bound': 5,
                'up_bound': 100,
                'step': 5,
                'bam_path': self.bam_path,
                'orf_bed': self.orf_bed,
            }
            rpkm_path = self.output_dir + '/Map_stat/satur/'
            rpkm_id = api_map.add_rpkm_table(rpkm_path, params=rpkm_params)
            # api_map.add_rpkm_box(rpkm_path, rpkm_id)
            api_map.add_rpkm_curve(rpkm_path, rpkm_id)
            coverage_path = self.output_dir + '/Map_stat/coverage/'
            cov_params = {
                'analysis': 'coverage',
                'min_len': 100,
                'bam_path': self.bam_path,
                'orf_bed': self.orf_bed,
            }
            api_map.add_coverage_table(coverage=coverage_path, name=None, params=cov_params)
            dup_path = self.output_dir + '/Map_stat/dup/'
            dup_params = {
                'analysis': 'dup',
                'quality_dup': 30,
                'bam_path': self.bam_path,
            }
            api_map.add_duplication_table(dup=dup_path, params=dup_params)
            corre_path = self.output_dir + '/Map_stat/correlation/'
            corre_params = {
                'analysis': 'dup',
            }
            api_map.add_correlation_table(correlation=corre_path, params=corre_params, express_id=self.express_id)
        if event['data'] == 'orf':
            self.move2outputdir(obj.output_dir, 'Gene_structure/orf')
            self.orf_bed = obj.option('bed').prop['path']
        if event['data'] == 'orf_len':
            self.move2outputdir(obj.output_dir, 'Gene_structure/orf')
            api_orf = self.api.denovo_gene_structure
            orf_path = self.output_dir + '/Gene_structure/orf/'
            if self.option('search_pfam'):
                pfam_path = orf_path + 'pfam_domain'
            else:
                pfam_path = None
            api_orf.add_orf_table(orf_bed=orf_path + 'Trinity.fasta.transdecoder.bed', reads_len_info=orf_path + 'reads_len_info', orf_domain=pfam_path, name=None, params=None)
        if event['data'] == 'ssr':
            self.move2outputdir(obj.output_dir, 'Gene_structure/ssr')
            api_ssr = self.api.denovo_gene_structure
            ssr_path = obj.work_dir
            if self.option('primer'):
                primer_path = ssr_path + '/gene.fasta.misa.results'
            else:
                primer_path = None
            ssr_params = {
                'orf_bed': self.orf_bed,
                'primer': self.option('primer')
            }
            api_ssr.add_ssr_table(ssr=ssr_path + '/gene.fasta.misa', ssr_primer=primer_path, ssr_stat=ssr_path + '/gene.fasta.statistics', name=None, params=ssr_params)
        if event['data'] == 'snp':
            self.move2outputdir(obj.output_dir, 'Gene_structure/snp')
            api_snp = self.api.denovo_gene_structure
            snp_path = self.output_dir + '/Gene_structure/snp'
            api_snp.add_snp_table(snp=snp_path)
        if event['data'] == 'exp_stat':
            self.move2outputdir(obj.output_dir, 'Express')
            self.bam_path = obj.option('bam_dir').prop['path']
            self.logger.info('%s' % self.exp_stat.diff_gene)
            # set api
            rsem_dir = self.output_dir + '/Express/rsem/'
            self.express_id = self.api_express.add_express(samples=self.samples, params=None, name=None, bam_path=self.bam_path, rsem_dir=rsem_dir)
            api_control = self.api.control
            if self.option('group_table').is_set:
                api_group = self.api.denovo_group
                group_id = api_group.add_ini_group_table(self.option('group_table').prop['path'], self.spname_spid)
                control_id = api_control.add_control(self.option('control_file').prop['path'], group_id)
                group_detail = api_group.get_group_detail(self.option('group_table').prop['path'], self.spname_spid, self.option('group_table').prop['group_scheme'][0])
            else:
                control_id = api_control.add_control(self.option('control_file').prop['path'], 'all')
            compare_column = list()
            diff_exp_dir = self.output_dir + '/Express/diff_exp/'
            diff_files = os.listdir(diff_exp_dir)
            for f in diff_files:
                if re.search(r'_edgr_stat.xls$', f):
                    con_exp = f.split('_edgr_stat.xls')[0].split('_vs_')
                    compare_column.append('|'.join(con_exp))
            diff_param = {
                'ci': self.option('diff_ci'),
                'rate': self.option('diff_rate'),
            }
            if self.option('group_table').is_set:
                express_diff_id = self.api_express.add_express_diff(params=diff_param, samples=self.samples, compare_column=compare_column, express_id=self.express_id, group_id=group_id, group_detail=group_detail, control_id=control_id, diff_exp_dir=diff_exp_dir)
            else:
                express_diff_id = self.api_express.add_express_diff(params=diff_param, samples=self.samples, compare_column=compare_column, express_id=self.express_id, group_id='all', group_detail={'all': sorted(self.api_sample.sample_ids)}, control_id=control_id, diff_exp_dir=diff_exp_dir)

            param_2 = {
                # 'express_diff_id': ,
                'compare_list': compare_column,
                'is_sum': True,
            }
            self.diff_gene_id = self.api_express.add_express(samples=self.samples, params=param_2, express_diff_id=express_diff_id, major=False)
            self.api_express.add_express_detail(self.diff_gene_id, diff_exp_dir + 'diff_count', diff_exp_dir + 'diff_fpkm', 'gene')
        if event['data'] == 'exp_diff':
            # set output
            self.move2outputdir(obj.output_dir, 'Express')
            # set api
            clust_path = self.output_dir + '/Express/cluster/hclust/'
            clust_files = os.listdir(clust_path)
            clust_params = {
                # diff_fpkm: ,
                'log': 10,
                'methor': 'hclust',
                'distance': 'euclidean',
                'sub_num': 5,
            }
            clust_id = self.api_express.add_cluster(clust_params, self.diff_gene_id, clust_path + 'samples_tree.txt', clust_path + 'genes_tree.txt', clust_path + 'hclust_heatmap.xls')
            for f in clust_files:
                if re.search(r'^subcluster_', f):
                    sub = f.split('_')[1]
                    self.api_express.add_cluster_detail(clust_id, sub, clust_path + f)
            net_path = os.path.join(self.output_dir + '/Express/network/')
            net_files = os.listdir(net_path)
            net_param = {
                # diff_fpkm: ,
                'softpower': 9,
                'similar': 0.75,
            }
            net_id = self.api_express.add_network(net_param, self.diff_gene_id, net_path + 'softPower.pdf', net_path + 'ModuleTree.pdf')
            self.api_express.add_network_detail(net_id, net_path + 'all_nodes.txt', net_path + 'all_edges.txt')
            for f in net_files:
                if re.search(r'^CytoscapeInput-edges-', f):
                    color = f.split('-edges-')[-1].split('.')[0]
                    self.api_express.add_network_module(net_id, net_path + f, color)
        if event['data'] == 'annotation':
            self.move2outputdir(obj.output_dir, 'Annotation')

    def run(self):
        self.filecheck.on('end', self.run_qc)
        self.filecheck.on('end', self.run_qc_stat, False)  # 质控前统计
        self.qc.on('end', self.run_qc_stat, True)  # 质控后统计
        self.qc.on('end', self.run_assemble)
        self.assemble.on('end', self.run_orf)
        self.assemble.on('end', self.run_exp_stat)
        self.orf.on('end', self.run_orf_len)
        self.final_tools.append(self.orf_len)
        if self.option('database'):
            self.assemble.on('end', self.run_annotation)
            self.final_tools.append(self.annotation)
        self.on_rely([self.orf, self.exp_stat], self.run_map_qc)
        self.final_tools.append(self.map_qc)
        if 'ssr' in self.option('gene_analysis'):
            self.orf.on('end', self.run_ssr)
            self.final_tools.append(self.ssr)
        if 'snp' in self.option('gene_analysis'):
            self.assemble.on('end', self.run_bwa)
            self.on_rely([self.bwa, self.orf], self.run_snp)
            self.final_tools.append(self.snp)
        if self.option('exp_analysis'):
            if ('go_rich' or 'kegg_rich') in self.option('exp_analysis'):
                self.on_rely([self.exp_stat, self.annotation], self.run_exp_diff)
            else:
                self.exp_stat.on('end', self.run_exp_diff)
        else:
            self.on_rely(self.final_tools, self.end)
        self.run_filecheck()
        super(TestBaseWorkflow, self).run()

    def end(self):
        self.send_files()
        super(TestBaseWorkflow, self).end()

    def send_files(self):
        self.logger.info('denovo_base upload files start')
        repaths = [
            ['.', "文件夹", "denovo rna 结果文件目录"],
            ["QC_stat", "文件夹", "样本数据统计文件目录"],
            ['QC_stat/before_qc', "文件夹", "质控前的样本数据统计文件目录"],
            ["QC_stat/before_qc/qualityStat/", "文件夹", "每个样本的fastq的质量统计文件的输出目录"],
            ["QC_stat/before_qc/fastq_stat.xls", "xls", "所有样本的fastq信息统计表"],
            ['QC_stat/after_qc', "文件夹", "质控后的样本数据统计文件目录"],
            ["QC_stat/after_qc/qualityStat/", "文件夹", "质控后的每个样本的fastq的质量统计文件的输出目录"],
            ["QC_stat/after_qc/fastq_stat.xls", "xls", "质控后的所有样本的fastq信息统计表"],
            ["QC_stat/after_qc/dup.xls", "xls", "所有样本的fastq序列重复统计表"],
            ["QC_stat/sickle_dir/", "文件夹", "质量剪切后的fastq文件输出目录"],
            ['Assemble', "文件夹", "Trinity拼接组装统计结果目录"],
            ['Assemble/transcript.iso.txt', "txt", "Trinity.fasta可变剪接体统计文件"],
            ['Assemble/trinity.fasta.stat.xls', "xls", "Trinity.fasta序列相关信息统计文件"],
            ['Gene_structure', "文件夹", "基因结构分析结果目录"],
            ['Gene_structure/snp', "文件夹", "snp分析结果目录"],
            ['Gene_structure/orf', "文件夹", "orf分析结果目录"],
            ["Gene_structure/orf/reads_len_info", "文件夹", "orf序列长度分布信息文件夹"],
            ['Gene_structure/ssr', "文件夹", "orf分析结果目录"],
            ["Gene_structure/ssr/misa_stat.xls", "xls", "ssr类型统计表"],
            ['Map_stat', "文件夹", "Mapping后质量统计结果目录"],
            ["Map_stat/coverage/", "文件夹", "基因覆盖度分析输出目录"],
            ["Map_stat/sorted_bam/", "文件夹", "每个样本排序后的bam文件输出目录"],
            ["Map_stat/dup/", "文件夹", "冗余序列分析输出目录"],
            ["Map_stat/satur/", "文件夹", "测序饱和度分析输出目录"],
            ["Map_stat/bam_stat.xls", "xls", "bam格式比对结果统计表"],
            ['Express/', "文件夹", "表达量分析结果目录"],
            ['Express/diff_exp', "文件夹", "表达量差异检测分析结果目录"],
            ['Express/rsem', "文件夹", "表达量计算分析结果目录"],
            ['Express/correlation', "文件夹", "表达量样本相关性分析结果目录"],
            ["Express/correlation/correlation_matrix.xls", "xls", "相关系数矩阵表"],
            ["Express/correlation/hcluster_tree_correlation_matrix.xls_average.tre", "xls", "相关系数树文件"],
            ["./Annotation", "", "DENOVO_RNA结果文件目录"],
            ['/Annotation/ncbi_taxonomy/query_taxons_detail.xls', 'xls', '序列详细物种分类文件'],
            ["/Annotation/blast_nr_statistics/output_evalue.xls", "xls", "blast结果E-value统计"],
            ["/Annotation/blast_nr_statistics/output_similar.xls", "xls", "blast结果similarity统计"],
            ["/Annotation/kegg/kegg_table.xls", "xls", "KEGG annotation table"],
            ["/Annotation/kegg/pathway_table.xls", "xls", "Sorted pathway table"],
            ["/Annotation/kegg/kegg_taxonomy.xls", "xls", "KEGG taxonomy summary"],
            ["/Annotation/go/blast2go.annot", "annot", "Go annotation based on blast output"],
            ["/Annotation/go/query_gos.list", "list", "Merged Go annotation"],
            ["/Annotation/go/go1234level_statistics.xls", "xls", "Go annotation on 4 levels"],
            ["/Annotation/go/go2level.xls", "xls", "Go annotation on level 2"],
            ["/Annotation/go/go3level.xls", "xls", "Go annotation on level 3"],
            ["/Annotation/go/go4level.xls", "xls", "Go annotation on level 4"],
            ["/Annotation/cog/cog_list.xls", "xls", "COG编号表"],
            ["/Annotation/cog/cog_summary.xls", "xls", "COG注释二级统计表"],
            ["/Annotation/cog/cog_table.xls", "xls", "序列COG注释详细表"],
            ["/Annotation/anno_stat", "", "denovo注释统计结果输出目录"],
            ["/Annotation/anno_stat/ncbi_taxonomy/", "dir", "nr统计结果目录"],
            ["/Annotation/anno_stat/cog_stat/", "dir", "cog统计结果目录"],
            ["/Annotation/anno_stat/go_stat/", "dir", "go统计结果目录"],
            ["/Annotation/anno_stat/kegg_stat/", "dir", "kegg统计结果目录"],
            ["/Annotation/anno_stat/blast/", "dir", "基因序列blast比对结果目录"],
            ["/Annotation/anno_stat/blast_nr_statistics/", "dir", "基因序列blast比对nr库统计结果目录"],
            ["/Annotation/anno_stat/blast/gene_kegg.xls", "xls", "基因序列blast比对kegg注释结果table"],
            ["/Annotation/anno_stat/blast/gene_nr.xls", "xls", "基因序列blast比对nr注释结果table"],
            ["/Annotation/anno_stat/blast/gene_nr.xls", "xls", "基因序列blast比对nr注释结果table"],
            ["/Annotation/anno_stat/blast/gene_string.xml", "xml", "基因序列blast比对string注释结果xml"],
            ["/Annotation/anno_stat/blast/gene_kegg.xml", "xml", "基因序列blast比对kegg注释结果xml"],
            ["/Annotation/anno_stat/blast/gene_string.xml", "xml", "基因序列blast比对string注释结果xml"],
            ["/Annotation/anno_stat/cog_stat/gene_cog_list.xls", "xls", "基因序列cog_list统计结果"],
            ["/Annotation/anno_stat/cog_stat/gene_cog_summary.xls", "xls", "基因序列cog_summary统计结果"],
            ["/Annotation/anno_stat/cog_stat/gene_cog_table.xls", "xls", "基因序列cog_table统计结果"],
            ["/Annotation/anno_stat/cog_stat/gene_cog_table.xls", "xls", "基因序列cog_table统计结果"],
            ["/Annotation/anno_stat/cog_stat/gene_cog_table.xls", "xls", "基因序列cog_table统计结果"],
            ["/Annotation/anno_stat/cog_stat/gene_cog_table.xls", "xls", "基因序列cog_table统计结果"],
            ["/Annotation/anno_stat/go_stat/gene_blast2go.annot", "annot", "Go annotation based on blast output of gene"],
            ["/Annotation/anno_stat/go_stat/gene_gos.list", "list", "Merged Go annotation of gene"],
            ["/Annotation/anno_stat/go_stat/gene_go1234level_statistics.xls", "xls", "Go annotation on 4 levels of gene"],
            ["/Annotation/anno_stat/go_stat/gene_go2level.xls", "xls", "Go annotation on level 2 of gene"],
            ["/Annotation/anno_stat/go_stat/gene_go3level.xls", "xls", "Go annotation on level 3 of gene"],
            ["/Annotation/anno_stat/go_stat/gene_go4level.xls", "xls", "Go annotation on level 4 of gene"],
            ["/Annotation/anno_stat/kegg_stat/gene_kegg_table.xls", "xls", "KEGG annotation table of gene"],
            ["/Annotation/anno_stat/kegg_stat/gene_pathway_table.xls", "xls", "Sorted pathway table of gene"],
            ["/Annotation/anno_stat/kegg_stat/gene_kegg_taxonomy.xls", "xls", "KEGG taxonomy summary of gene"],
            ["/Annotation/anno_stat/kegg_stat/gene_kegg_layer.xls", "xls", "KEGG taxonomy summary of gene"],
            ["/Annotation/anno_stat/kegg_stat/gene_pathway/", "dir", "基因的标红pathway图"],
            ['/ncbi_taxonomy/gene_taxons_detail.xls', 'xls', '基因序列详细物种分类文件'],
            ["/Annotation/anno_stat/blast_nr_statistics/gene_nr_evalue.xls", "xls", "基因序列blast结果E-value统计"],
            ["/Annotation/anno_stat/blast_nr_statistics/gene_nr_similar.xls", "xls", "基因序列blast结果similarity统计"],
            ["/Annotation/anno_stat/ncbi_taxonomy/gene_taxons.xls", "xls", "基因序列nr物种注释表"],
            ["/Annotation/anno_stat/ncbi_taxonomy/query_taxons.xls", "xls", "nr物种注释表"],
            ["/Annotation/anno_stat/all_annotation_statistics.xls", "xls", "注释统计总览表"],
            ["/Annotation/anno_stat/all_annotation.xls", "xls", "注释统计表"],
        ]
        regexps = [
            [r'Assemble/.*_length\.distribut\.txt$', 'txt', '长度分布信息统计文件'],
            [r"Gene_structure/snp/.*snp_position_stat\.xls", "xls", "样本snp编码位置信息统计表"],
            [r"Gene_structure/snp/.*snp_type_stat\.xls", "xls", "样本snp类型统计表"],
            [r"Gene_structure/snp/.*snp\.xls", "xls", "样本snp信息表"],
            [r"Gene_structure/ssr/.*misa$", "misa", "ssr结果"],
            [r"Gene_structure/orf/reads_len_info/.*reads_len_info\.txt$", "xls", "orf序列长度分布信息文件"],
            [r"Gene_structure/orf/transdecoder.pep$", "fasta", "蛋白质序列文件"],
            [r"Gene_structure/orf/transdecoder.cds$", "fasta", "cds序列文件"],
            [r"Gene_structure/orf/transdecoder.bed$", "bed", "orf位置信息bed格式文件"],
            [r"Map_stat/dup/.*pos\.DupRate\.xls", "xls", "比对到基因组的序列的冗余统计表"],
            [r"Map_stat/dup/.*seq\.DupRate\.xls", "xls", "所有序列的冗余统计表"],
            [r"Map_stat/satur/.*eRPKM\.xls", "xls", "RPKM表"],
            [r"Map_stat/coverage/.*cluster_percent\.xls", "xls", "饱和度作图数据"],
            [r"Express/diff_exp/.*_edgr_stat\.xls$", "xls", "edger统计结果文件"],
            [r"Express/rsem/.*results$", "xls", "单样本rsem分析结果表"],
            [r"Express/rsem/.*matrix$", "xls", "表达量矩阵"],
            [r"/Annotation/nrblast/.+_vs_.+\.xml", "xml", "blast比对nr输出结果，xml格式"],
            [r"/Annotation/nrblast/.+_vs_.+\.xls", "xls", "blast比对nr输出结果，表格(制表符分隔)格式"],
            [r"/Annotation/stringblast/.+_vs_.+\.xml", "xml", "blast比对string输出结果，xml格式"],
            [r"/Annotation/stringblast/.+_vs_.+\.xls", "xls", "blast比对string输出结果，表格(制表符分隔)格式"],
            [r"/Annotation/keggblast/.+_vs_.+\.xml", "xml", "blast比对kegg输出结果，xml格式"],
            [r"/Annotation/keggblast/.+_vs_.+\.xls", "xls", "blast比对kegg输出结果，表格(制表符分隔)格式"],
            [r"/Annotation/kegg/pathways/ko.\d+", 'pdf', '标红pathway图'],
            [r"/Annotation/blast_nr_statistics/.*_evalue\.xls", "xls", "比对结果E-value分布图"],
            [r"/Annotation/blast_nr_statistics/.*_similar\.xls", "xls", "比对结果相似度分布图"],
            ["^/Annotation/anno_stat/ncbi_taxonomy/nr_taxon_stat", "xls", "nr物种分类统计表"],
        ]
        if self.option("search_pfam"):
            repaths += [["Gene_structure/orf/pfam_domain", "", "Pfam比对蛋白域结果信息"]]
        if self.option('fq_type') == 'SE':
            repaths += [
                ['QC_stat/clip_dir', "文件夹", "SE去接头后的fastq文件输出目录"]
            ]
        else:
            repaths += [
                ["QC_stat/seqprep_dir/", "文件夹", "PE去接头后fastq文件输出目录"]
            ]
        if self.exp_stat.diff_gene:
            regexps += [
                [r"Express/cluster/hclust/subcluster_", "xls", "子聚类热图数据"],
                [r"Express/network/CytoscapeInput.*", "txt", "Cytoscape作图数据"]
            ]
            repaths += [
                ["Express/diff_exp/diff_fpkm", "xls", "差异基因表达量表"],
                ["Express/diff_exp/diff_count", "xls", "差异基因计数表"],
                ['Express/network', "文件夹", "差异基因网络共表达分析结果目录"],
                ["Express/network/all_edges.txt", "txt", "edges结果信息"],
                ["Express/network/all_nodes.txt ", "txt", "nodes结果信息"],
                ["Express/network/removeGene.xls ", "xls", "移除的基因信息"],
                ["Express/network/removeSample.xls ", "xls", "移除的样本信息"],
                ["Express/network/softPower.pdf", "pdf", "softpower相关信息"],
                ["Express/network/ModuleTree.pdf", "pdf", "ModuleTree图"],
                ["Express/network/eigengeneClustering.pdf", "pdf", "eigengeneClustering图"],
                ["Express/network/eigenGeneHeatmap.pdf", "pdf", "eigenGeneHeatmap图"],
                ["Express/network/networkHeatmap.pdf", "pdf", "networkHeatmap图"],
                ["Express/network/sampleClustering.pdf", "pdf", "sampleClustering图"],
                ["Express/cluster", "文件夹", "差异基因聚类分析分析结果目录"],
                ["Express/cluster/hclust/", "", "层级聚类分析结果目录"],
                ["Express/cluster/hclust/hc_gene_order", "txt", "按基因聚类的基因排序列表"],
                ["Express/cluster/hclust/hc_sample_order", "txt", "按样本聚类的样本排序列表"],
                ["Express/cluster/hclust/hclust_heatmap.xls", "xls", "层级聚类热图数据"]
            ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        # for i in self.get_upload_files():
        #     self.logger.info('upload file:{}'.format(str(i)))
        self.logger.info('denovo_base upload files end')

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
