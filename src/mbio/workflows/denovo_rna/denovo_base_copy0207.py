# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

"""无参转录组基础分析"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import os
import shutil
import re


class DenovoBaseWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        version = v1.0
        last_modify = 20161201
        """
        self._sheet = wsheet_object
        super(DenovoBaseWorkflow, self).__init__(wsheet_object)
        print self._parent
        options = [
            {"name": "fastq_dir", "type": "infile", 'format': "sequence.fastq,sequence.fastq_dir"},  # fastq文件夹
            {"name": "fq_type", "type": "string"},  # PE OR SE
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 有生物学重复的时候的分组文件
            {"name": "control_file", "type": "infile", "format": "denovo_rna.express.control_table"},  # 对照组文件，格式同分组文件
            {"name": "qc_quality", "type": "int", "default": 30},  # 质量剪切中保留的最小质量值
            {"name": "qc_length", "type": "int", "default": 50},   # 质量剪切中保留的最短序列长度
            {"name": "search_pfam", "type": "bool", "default": True},  # orf 是否比对Pfam数据库
            {"name": "primer", "type": "bool", "default": True},  # 是否设计SSR引物
            {"name": "kmer_size", "type": "int", "default": 25},
            {"name": "min_kmer_cov", "type": "int", "default": 2},
            {"name": "min_contig_length", "type": "int", "default": 200},  # trinity报告出的最短的contig长度。默认为200
            {"name": "SS_lib_type", "type": "string", "default": 'none'},  # reads的方向，成对的reads: RF or FR; 不成对的reads: F or R，默认情况下，不设置此参数
            {"name": "exp_way", "type": "string", "default": "fpkm"},  # edger离散值
            {"name": "diff_ci", "type": "float", "default": 0.01},  # 显著性水平
            {"name": "database", "type": "string", "default": 'go,nr,cog,kegg'},  # 默认全部四个注释
            {"name": "nr_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "string_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "kegg_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "exp_analysis", "type": "string", "default": "cluster,network,kegg_rich,go_rich,kegg_regulate"},
            {"name": "gene_analysis", "type": "string", "default": "orf,snp"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.filecheck = self.add_tool("denovo_rna.filecheck.file_denovo")
        self.qc = self.add_module("denovo_rna.qc.quality_control")
        self.qc_stat_before = self.add_module("denovo_rna.qc.qc_stat")
        self.qc_stat_after = self.add_module("denovo_rna.qc.qc_stat")
        self.assemble = self.add_tool("denovo_rna.assemble.assemble")
        self.bam_stat = self.add_tool("denovo_rna.mapping.bam_stat")
        self.annotation = self.add_module('denovo_rna.annotation.denovo_annotation')
        self.orf = self.add_tool("denovo_rna.gene_structure.orf")
        self.ssr = self.add_tool("denovo_rna.gene_structure.ssr")
        self.bwa = self.add_module("denovo_rna.mapping.bwa_samtools")
        self.snp = self.add_module("denovo_rna.gene_structure.snp")
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
        self.express_diff_id = None
        self.api_express = self.api.denovo_express
        self.api_anno = self.api.denovo_annotation

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
        if self.option('qc_quality') > 40 or self.option('diff_ci') < -15:
            raise OptionError('qc_quality不在所给范围内[-15,40]')
        if self.option('qc_length') > 150 or self.option('diff_ci') <= 0:
            raise OptionError('qc_length不在所给范围内(0,150]')
        if self.option('diff_ci') > 1 or self.option('diff_ci') < 0:
            raise OptionError('显著性水平不在所给范围内[0,1]')
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
            if i not in ['cluster', 'network', 'go_rich', 'kegg_rich', 'go_regulate', 'kegg_regulate']:
                raise OptionError("差异性研究没有{}，请检查".format(i))
        for i in self.option('gene_analysis').split(','):
            if i not in ['orf', 'ssr', 'snp']:
                raise OptionError("基因结构分析没有{}，请检查".format(i))
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
        self.filecheck.run()

    def run_qc(self):
        self.qc.set_options({
            'fastq_dir': self.option('fastq_dir'),
            'fq_type': self.option('fq_type'),
            'quality_q': self.option('qc_quality'),
            'length_q': self.option('qc_length')
        })
        self.qc.on('end', self.set_output, 'qc')
        self.qc.on('start', self.set_step, {'start': self.step.qcstat})
        self.qc.on('end', self.set_step, {'end': self.step.qcstat})
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
            self.qc_stat_after.run()
        else:
            self.qc_stat_before.on('end', self.set_output, 'qc_stat_before')
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
        self.assemble.on('end', self.set_step, {'start': self.step.annotation})
        self.assemble.run()

    def run_orf(self):
        orf_opts = {
            'fasta': self.assemble.option('trinity_fa'),
            'search_pfam': self.option('search_pfam')
        }
        self.orf.set_options(orf_opts)
        self.orf.on('end', self.set_output, 'orf')
        self.orf.on('start', self.set_step, {'start': self.step.gene_structure})
        self.orf.run()

    def run_orf_len(self):
        orf_fasta = self.orf.work_dir + '/ORF_fasta'
        self.orf_len.set_options({'fasta_path': orf_fasta})
        self.orf_len.on('end', self.set_output, 'orf_len')
        self.orf_len.run()

    def run_bwa(self):
        bwa_opts = {
            'ref_fasta': self.assemble.option('gene_fa'),
            'fq_type': self.option('fq_type'),
            'fastq_dir': self.qc.option('sickle_dir'),
        }
        self.bwa.set_options(bwa_opts)
        self.bwa.run()

    def run_snp(self):
        snp_opts = {
            'bed': self.orf.option('bed'),
            'bam': self.bwa.option('out_bam'),
            'ref_fasta': self.assemble.option('gene_fa')
        }
        self.snp.set_options(snp_opts)
        self.snp.on('end', self.set_output, 'snp')
        self.snp.run()

    def run_ssr(self):
        ssr_opts = {
            'fasta': self.assemble.option('gene_fa'),
            'bed': self.orf.option('bed'),
            'primer': self.option('primer')
        }
        self.ssr.set_options(ssr_opts)
        self.ssr.on('end', self.set_output, 'ssr')
        self.ssr.run()

    def run_blast_test(self):
        self.blast_modules = []
        self.gene_list = self.assemble.option('gene_full_name').prop['gene_list']
        blast_lines = int(self.assemble.option('trinity_fa').prop['seq_number']) / 10
        self.logger.info('.......blast_lines:%s' % blast_lines)
        blast_opts = {
            'query': self.assemble.option('trinity_fa'),
            'query_type': 'nucl',
            'database': None,
            'blast': 'blastx',
            'evalue': None,
            'outfmt': 6,
            'lines': blast_lines,
        }
        if 'go' in self.option('database') or 'nr' in self.option('database'):
            self.blast_nr = self.add_module('align.blast')
            blast_opts.update(
                {'database': 'nr', 'evalue': self.option('nr_blast_evalue')}
            )
            self.blast_nr.set_options(blast_opts)
            self.blast_modules.append(self.blast_nr)
            self.blast_nr.on('end', self.set_output, 'nrblast')
            # self.blast_nr.run()
        if 'cog' in self.option('database'):
            self.blast_string = self.add_module('align.blast')
            blast_opts.update(
                {'database': 'string', 'evalue': self.option('string_blast_evalue')}
            )
            self.blast_string.set_options(blast_opts)
            self.blast_modules.append(self.blast_string)
            self.blast_string.on('end', self.set_output, 'stringblast')
            # self.blast_string.run()
        if 'kegg' in self.option('database'):
            self.blast_kegg = self.add_module('align.blast')
            blast_opts.update(
                {'database': 'kegg', 'evalue': self.option('kegg_blast_evalue')}
            )
            self.blast_kegg.set_options(blast_opts)
            self.blast_modules.append(self.blast_kegg)
            self.blast_kegg.on('end', self.set_output, 'keggblast')
            # self.blast_kegg.run()
        self.on_rely(self.blast_modules, self.run_annotation)
        self.test_run(self.blast_kegg)
        self.test_run(self.blast_nr)
        self.test_run(self.blast_string)
        self.blast_kegg.option('outxml', '/mnt/ilustre/users/sanger-dev/workspace/20170103/DenovoBase_sg_5538/Blast2/CatBlastout1/output/blast.xml')
        self.blast_kegg.option('outtable', '/mnt/ilustre/users/sanger-dev/workspace/20170103/DenovoBase_sg_5538/Blast2/CatBlastout/output/blast_table.xls')
        self.blast_nr.option('outxml', '/mnt/ilustre/users/sanger-dev/workspace/20170103/DenovoBase_sg_5538/Blast/CatBlastout1/output/blast.xml')
        self.blast_nr.option('outtable', '/mnt/ilustre/users/sanger-dev/workspace/20170103/DenovoBase_sg_5538/Blast/CatBlastout/output/blast_table.xls')
        self.blast_string.option('outxml', '/mnt/ilustre/users/sanger-dev/workspace/20170103/DenovoBase_sg_5538/Blast1/CatBlastout1/output/blast.xml')
        self.blast_string.option('outtable', '/mnt/ilustre/users/sanger-dev/workspace/20170103/DenovoBase_sg_5538/Blast1/CatBlastout/output/blast_table.xls')
        self.test_end(self.blast_kegg)
        self.test_end(self.blast_nr)
        self.test_end(self.blast_string)

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
        self.blast_modules = []
        self.gene_list = self.assemble.option('gene_full_name').prop['gene_list']
        blast_lines = int(self.assemble.option('trinity_fa').prop['seq_number']) / 10
        self.logger.info('.......blast_lines:%s' % blast_lines)
        blast_opts = {
            'query': self.assemble.option('trinity_fa'),
            'query_type': 'nucl',
            'database': None,
            'blast': 'blastx',
            'evalue': None,
            'outfmt': 6,
            'lines': blast_lines,
        }
        if 'go' in self.option('database') or 'nr' in self.option('database'):
            self.blast_nr = self.add_module('align.blast')
            blast_opts.update(
                {'database': 'nr', 'evalue': self.option('nr_blast_evalue')}
            )
            self.blast_nr.set_options(blast_opts)
            self.blast_modules.append(self.blast_nr)
            self.blast_nr.on('end', self.set_output, 'nrblast')
            self.blast_nr.run()
        if 'cog' in self.option('database'):
            self.blast_string = self.add_module('align.blast')
            blast_opts.update(
                {'database': 'string', 'evalue': self.option('string_blast_evalue')}
            )
            self.blast_string.set_options(blast_opts)
            self.blast_modules.append(self.blast_string)
            self.blast_string.on('end', self.set_output, 'stringblast')
            self.blast_string.run()
        if 'kegg' in self.option('database'):
            self.blast_kegg = self.add_module('align.blast')
            blast_opts.update(
                {'database': 'kegg', 'evalue': self.option('kegg_blast_evalue')}
            )
            self.blast_kegg.set_options(blast_opts)
            self.blast_modules.append(self.blast_kegg)
            self.blast_kegg.on('end', self.set_output, 'keggblast')
            self.blast_kegg.run()
        self.on_rely(self.blast_modules, self.run_annotation)

    def run_annotation(self):
        anno_opts = {
            'gene_file': self.assemble.option('gene_full_name'),
        }
        if 'go' in self.option('database'):
            anno_opts.update({
                'go_annot': True,
                'blast_nr_xml': self.blast_nr.option('outxml')
            })
        else:
            anno_opts.update({'go_annot': False})
        if 'nr' in self.option('database'):
            anno_opts.update({
                'nr_annot': True,
                'blast_nr_xml': self.blast_nr.option('outxml'),
                'blast_nr_table': self.blast_nr.option('outtable')
            })
        else:
            anno_opts.update({'nr_annot': False})
        if 'kegg' in self.option('database'):
            anno_opts.update({
                'blast_kegg_xml': self.blast_kegg.option('outxml'),
                'blast_kegg_table': self.blast_kegg.option('outtable')
            })
        if 'cog' in self.option('database'):
            anno_opts.update({
                'blast_string_xml': self.blast_string.option('outxml'),
                'blast_string_table': self.blast_string.option('outtable')
            })
        self.logger.info('....anno_opts:%s' % anno_opts)
        self.annotation.set_options(anno_opts)
        self.annotation.on('end', self.set_output, 'annotation')
        self.annotation.on('end', self.set_step, {'end': self.step.annotation})
        self.annotation.run()

    def run_exp_stat(self):
        exp_stat_opts = {
            'fq_type': self.option('fq_type'),
            'rsem_fa': self.assemble.option('trinity_fa'),
            'control_file': self.option('control_file'),
            'exp_way': self.option('exp_way'),
            'diff_ci': self.option('diff_ci'),
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
        self.exp_stat.run()

    def run_bam_stat(self):
        opts = {'bam': self.exp_stat.option('bam_dir')}
        self.bam_stat.set_options(opts)
        self.bam_stat.on('end', self.set_output, 'bam_stat')
        self.bam_stat.run()

    def run_exp_diff(self):
        if self.exp_stat.diff_gene:
            exp_diff_opts = {
                'diff_fpkm': self.exp_stat.option('diff_fpkm'),
                'analysis': self.option('exp_analysis')
            }
            if 'network' in self.option('exp_analysis'):
                exp_diff_opts.update({'diff_list': self.exp_stat.option('diff_list')})
            if 'kegg_rich' in self.option('exp_analysis'):
                exp_diff_opts.update({
                    'gene_kegg_table': self.annotation.option('gene_kegg_table'),
                    'diff_list_dir': self.exp_stat.option('diff_list_dir'),
                    'all_list': self.exp_stat.option('all_list'),
                })
            if 'kegg_regulate' in self.option('exp_analysis'):
                exp_diff_opts.update({
                    'gene_kegg_table': self.annotation.option('gene_kegg_table'),
                    'diff_stat_dir': self.exp_stat.diff_exp.option('regulate_edgrstat_dir')
                })
            if 'go_rich' in self.option('exp_analysis'):
                exp_diff_opts.update({
                    'gene_go_list': self.annotation.option('gene_go_list'),
                    'diff_list_dir': self.exp_stat.option('diff_list_dir'),
                    'all_list': self.exp_stat.option('all_list')
                })
            if 'go_regulate' in self.option('exp_analysis'):
                exp_diff_opts.update({
                    'gene_go_level_2': self.annotation.option('gene_go_level_2'),
                    'diff_stat_dir': self.exp_stat.diff_exp.option('regulate_edgrstat_dir')
                })
            self.exp_diff.set_options(exp_diff_opts)
            self.exp_diff.on('end', self.set_output, 'exp_diff')
            self.exp_diff.on('end', self.set_step, {'end': self.step.express})
            self.final_tools.append(self.exp_diff)
            self.on_rely(self.final_tools, self.end)
            self.exp_diff.run()
        else:
            self.logger.info('输入文件数据量过小，没有检测到差异基因，差异基因相关分析将忽略')
            self.set_step(event={'data': {'end': self.step.express}})
            self.logger.info('......final_tools: %s' % self.final_tools)
            all_end = []
            for i in self.final_tools:
                self.logger.info('......tool: %s is_end is %s' % (i, i.is_end))
                all_end.append(i.is_end)
            if all(all_end):
                self.end()
            else:
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
            # api_sample = self.api.denovo_rna_sample
            self.api_sample.add_samples_info(qc_stat=qc_stat_info, qc_adapt=None, fq_type=self.option('fq_type'))
            self.api_sample.add_gragh_info(quality_stat, about_qc='before')
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
            self.sample_ids = self.api_sample.add_samples_info(qc_stat=qc_stat_info, qc_adapt=qc_adapt, fq_type=self.option('fq_type'))
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
        if event['data'] == 'bam_stat':
            self.move2outputdir(obj.output_dir, 'Mapping_stat')
            api_map = self.api.denovo_rna_mapping
            api_map.add_mapping_stat(stat_file=self.output_dir + '/Mapping_stat/bam_stat.xls')
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
                'primer': self.option('primer'),
                'submit_location': 'sg_denovo_ssr'
            }
            ssr_id = api_ssr.add_ssr_table(ssr=ssr_path + '/gene.fasta.misa', ssr_primer=primer_path, ssr_stat=ssr_path + '/gene.fasta.statistics', name=None, params=ssr_params)
            # update sg_status
            self.update_status_api.add_denovo_status(table_id=str(ssr_id), type_name='sg_denovo_ssr')
        if event['data'] == 'snp':
            self.move2outputdir(obj.output_dir, 'Gene_structure/snp')
            api_snp = self.api.denovo_gene_structure
            snp_path = self.output_dir + '/Gene_structure/snp'
            api_snp.add_snp_table(snp=snp_path)
        if event['data'] == 'exp_stat':
            self.move2outputdir(obj.output_dir, 'Express_stat')
            self.bam_path = obj.option('bam_dir').prop['path']
            self.logger.info('%s' % self.exp_stat.diff_gene)
            # set api
            # add express file and rsem result
            rsem_dir = self.output_dir + '/Express_stat/rsem/'
            gene_distri = self.exp_stat.merge_rsem.work_dir + '/gene_distribution.xls'
            tran_distri = self.exp_stat.merge_rsem.work_dir + '/tran_distribution.xls'
            self.express_id = self.api_express.add_express(samples=self.samples, params=None, name=None, bam_path=self.bam_path, rsem_dir=rsem_dir, gene_distri=gene_distri, tran_distri=tran_distri)
            # add control file and group_table
            api_control = self.api.control
            if self.option('group_table').is_set:
                api_group = self.api.denovo_group
                group_id = api_group.add_ini_group_table(self.option('group_table').prop['path'], self.spname_spid)
                control_id = api_control.add_control(self.option('control_file').prop['path'], group_id)
                group_detail = api_group.get_group_detail(self.option('group_table').prop['path'], self.spname_spid, self.option('group_table').prop['group_scheme'][0])
            else:
                control_id = api_control.add_control(self.option('control_file').prop['path'], 'all')
            # add diff express analysis
            compare_column = list()
            diff_exp_dir = self.output_dir + '/Express_stat/diff_exp/'
            diff_files = os.listdir(diff_exp_dir)
            for f in diff_files:
                if re.search(r'_edgr_stat.xls$', f):
                    con_exp = f.split('_edgr_stat.xls')[0].split('_vs_')
                    compare_column.append('|'.join(con_exp))
            diff_param = {
                'ci': self.option('diff_ci'),
                'submit_location': 'sg_denovo_express_diff'
            }
            if self.option('group_table').is_set:
                self.express_diff_id = self.api_express.add_express_diff(params=diff_param, samples=self.samples, compare_column=compare_column, express_id=self.express_id, group_id=group_id, group_detail=group_detail, control_id=control_id, diff_exp_dir=diff_exp_dir)
            else:
                self.express_diff_id = self.api_express.add_express_diff(params=diff_param, samples=self.samples, compare_column=compare_column, express_id=self.express_id, group_id='all', group_detail=self.sample_ids, control_id=control_id, diff_exp_dir=diff_exp_dir)
            # update sg_status
            self.update_status_api.add_denovo_status(table_id=str(self.express_diff_id), type_name='sg_denovo_express_diff')
            # add diff fpkm file
            param_2 = {
                # 'express_diff_id': ,
                'compare_list': compare_column,
                'is_sum': True,
                'submit_location': 'sg_denovo_express'
            }
            if not self.exp_stat.diff_gene:
                param_2['desc'] = '数据量太小，未检测到差异基因！'
            self.diff_gene_id = self.api_express.add_express(samples=self.samples, params=param_2, express_diff_id=self.express_diff_id, major=False)
            if self.exp_stat.diff_gene:
                self.api_express.add_express_detail(self.diff_gene_id, diff_exp_dir + 'diff_count', diff_exp_dir + 'diff_fpkm', 'gene')
            # add correlation file
            corr_api = self.api.denovo_rna_mapping
            corr_api.add_correlation_table(correlation=self.output_dir + '/Express_stat/gene_correlation/', express_id=self.express_id, detail=True, seq_type='gene')
            corr_api.add_correlation_table(correlation=self.output_dir + '/Express_stat/tran_correlation/', express_id=self.express_id, detail=True, seq_type='transcript')
        if event['data'] == 'exp_diff':
            # set output
            self.move2outputdir(obj.output_dir, 'Diff_express')
            # set api
            # set cluster
            if 'cluster' in self.option('exp_analysis'):
                clust_path = self.output_dir + '/Diff_express/cluster/hclust/'
                clust_files = os.listdir(clust_path)
                clust_params = {
                    'submit_location': 'sg_denovo_cluster',
                    'log': 10,
                    'methor': 'hclust',
                    'distance': 'euclidean',
                    'sub_num': 5,
                }
                api_clust = self.api.denovo_cluster
                with open(obj.cluster.work_dir + '/hc_sample_order', 'rb') as s:
                    samples = s.readlines()[0].strip('\n')
                with open(obj.cluster.work_dir + '/hc_gene_order', 'rb') as s:
                    genes = s.readlines()[0].strip('\n')
                clust_id = api_clust.add_cluster(params=clust_params, express_id=self.diff_gene_id, sample_tree=clust_path + 'samples_tree.txt', gene_tree=clust_path + 'genes_tree.txt', samples=samples, genes=genes)
                for f in clust_files:
                    if re.search(r'^subcluster_', f):
                        sub = f.split('_')[1]
                        api_clust.add_cluster_detail(clust_id, sub, clust_path + f)
                # update sg_status
                self.update_status_api.add_denovo_status(table_id=str(clust_id), type_name='sg_denovo_cluster')
            # set go rich
            if 'go_rich' in self.option('exp_analysis'):
                go_rich_api = self.api.denovo_go_enrich
                go_rich_path = os.path.join(self.output_dir + '/Diff_express/go_rich/')
                go_rich_dirs = os.listdir(go_rich_path)
                for d in go_rich_dirs:
                    path1 = go_rich_path + d
                    go_file, png = None, None
                    go_rich_param = {
                        'analysis_type': 'enrich',
                        'submit_location': 'sg_denovo_go_enrich',
                        'pval': 0.05,
                        'method': 'fdr',
                        'compare': ','.join(d.split('_vs_'))
                    }
                    for f in os.listdir(path1):
                        if re.match(r'go_enrich', f):
                            go_file = os.path.join(path1, f)
                        if re.search(r'png$', f):
                            png = os.path.join(path1, f)
                    go_id = go_rich_api.add_go_enrich(params=go_rich_param, go_graph_dir=png, go_enrich_dir=go_file, express_diff_id=self.express_diff_id)
                    # update sg_status
                    self.update_status_api.add_denovo_status(table_id=str(go_id), type_name='sg_denovo_go_enrich')
            # set go regulate
            if 'go_regulate' in self.option('exp_analysis'):
                go_regulate_api = self.api.denovo_go_regulate
                go_regulate_path = os.path.join(self.output_dir + '/Diff_express/go_regulate/')
                go_regulate_dirs = os.listdir(go_regulate_path)
                for d in go_regulate_dirs:
                    go_regu_param = {
                        'analysis_type': 'regulate',
                        'submit_location': 'sg_denovo_go_regulate',
                        'compare': ','.join(d.split('_vs_'))
                    }
                    path2 = go_regulate_path + d
                    go_file, png = None, None
                    f = os.path.join(path2, os.listdir(path2)[0])
                    go_regu_id = go_regulate_api.add_go_regulate(params=go_regu_param, go_regulate_dir=f, express_diff_id=self.express_diff_id)
                    # update sg_status
                    self.update_status_api.add_denovo_status(table_id=str(go_regu_id), type_name='sg_denovo_go_regulate')
            # set kegg tich
            if 'kegg_rich' in self.option('exp_analysis'):
                kegg_rich_api = self.api.denovo_kegg_rich
                kegg_rich_path = os.path.join(self.output_dir + '/Diff_express/kegg_rich/')
                kegg_rich_dirs = os.listdir(kegg_rich_path)
                for d in kegg_rich_dirs:
                    kegg_rich_param = {
                        'analysis_type': 'enrich',
                        'submit_location': 'sg_denovo_kegg_enrich',
                        'method': 'BH',
                        'compare': ','.join(d.split('_vs_'))
                    }
                    path3 = kegg_rich_path + d
                    f = os.path.join(path3, os.listdir(path3)[0])
                    kegg_id = kegg_rich_api.add_kegg_rich(params=kegg_rich_param, kegg_enrich_table=f, express_diff_id=self.express_diff_id)
                    # update sg_status
                    self.update_status_api.add_denovo_status(table_id=str(kegg_id), type_name='sg_denovo_kegg_enrich')
            # set kegg regulate
            if 'kegg_regulate' in self.option('exp_analysis'):
                kegg_regulate_api = self.api.denovo_kegg_regulate
                kegg_regulate_path = os.path.join(self.output_dir + '/Diff_express/kegg_regulate/')
                kegg_regulate_dirs = os.listdir(kegg_regulate_path)
                for d in kegg_regulate_dirs:
                    kegg_regu_param = {
                        'analysis_type': 'regulate',
                        'submit_location': 'sg_denovo_kegg_regulate',
                        'compare': ','.join(d.split('_vs_'))
                    }
                    path4 = kegg_regulate_path + d
                    stat_path = None
                    for f in os.listdir(path4):
                        if re.search(r'.xls$', f):
                            stat_path = os.path.join(path4, f)
                    pathway = path4 + '/pathways/'
                    kegg_regu_id = kegg_regulate_api.add_kegg_regulate(params=kegg_regu_param, kegg_regulate_table=stat_path, pathways_dir=pathway, express_diff_id=self.express_diff_id)
                    # update sg_status
                    self.update_status_api.add_denovo_status(table_id=str(kegg_regu_id), type_name='sg_denovo_kegg_regulate')
        if event['data'] == 'annotation':
            self.move2outputdir(obj.output_dir, 'Annotation')
            # set api
            self.api_anno.add_annotation(anno_stat_dir=obj.output_dir, databases=self.option('database'))
        if event['data'] == 'nrblast':
            self.move2outputdir(obj.output_dir, 'Annotation/nrblast')
            # blastfile = self.output_dir + '/Annotation/nrblast/' + os.listdir(self.output_dir + '/Annotation/nrblast/')[0]
            blastfile = '/mnt/ilustre/users/sanger-dev/workspace/20170103/DenovoBase_sg_5538/Blast/CatBlastout/output/blast_table.xls'
            self.api_anno.add_blast(blast_pro='blastp', blast_db='nr', e_value=self.option('nr_blast_evalue'), blast_path=blastfile, gene_list=self.gene_list)
        if event['data'] == 'keggblast':
            self.move2outputdir(obj.output_dir, 'Annotation/keggblast')
            blastfile = '/mnt/ilustre/users/sanger-dev/workspace/20170103/DenovoBase_sg_5538/Blast2/CatBlastout/output/blast_table.xls'
            # blastfile = self.output_dir + '/Annotation/keggblast/' + os.listdir(self.output_dir + '/Annotation/keggblast/')[0]
            self.api_anno.add_blast(blast_pro='blastp', blast_db='kegg', e_value=self.option('kegg_blast_evalue'), blast_path=blastfile, gene_list=self.gene_list)
        if event['data'] == 'stringblast':
            self.move2outputdir(obj.output_dir, 'Annotation/stringblast')
            blastfile = '/mnt/ilustre/users/sanger-dev/workspace/20170103/DenovoBase_sg_5538/Blast1/CatBlastout/output/blast_table.xls'
            # blastfile = self.output_dir + '/Annotation/stringblast/' + os.listdir(self.output_dir + '/Annotation/stringblast/')[0]
            self.api_anno.add_blast(blast_pro='blastp', blast_db='string', e_value=self.option('string_blast_evalue'), blast_path=blastfile, gene_list=self.gene_list)

    def run(self):
        self.logger.info('...options:%s' % self._options)
        self.filecheck.on('end', self.run_qc)
        self.filecheck.on('end', self.run_qc_stat, False)  # 质控前统计
        self.qc.on('end', self.run_qc_stat, True)  # 质控后统计
        self.qc.on('end', self.run_assemble)
        self.assemble.on('end', self.run_orf)
        self.assemble.on('end', self.run_exp_stat)
        self.exp_stat.on('end', self.run_bam_stat)
        self.orf.on('end', self.run_orf_len)
        self.final_tools.append(self.orf_len)
        self.final_tools.append(self.bam_stat)
        gene_stru = [self.orf_len]
        if self.option('database'):
            self.assemble.on('end', self.run_blast_test)
            # self.assemble.on('end', self.run_blast)
            self.final_tools.append(self.annotation)
        if 'ssr' in self.option('gene_analysis'):
            self.orf.on('end', self.run_ssr)
            self.final_tools.append(self.ssr)
            gene_stru.append(self.ssr)
        if 'snp' in self.option('gene_analysis'):
            self.assemble.on('end', self.run_bwa)
            self.on_rely([self.bwa, self.orf], self.run_snp)
            self.final_tools.append(self.snp)
            gene_stru.append(self.snp)
        self.logger.info('........gene_stru:%s' % gene_stru)
        if len(gene_stru) == 1:
            self.on('end', self.set_step, {'end': self.step.gene_structure})
        else:
            self.on_rely(gene_stru, self.set_step, {'end': self.step.gene_structure})
        if self.option('exp_analysis'):
            if 'go_rich' or 'kegg_rich' or 'go_regulate' or 'kegg_regulate' in self.option('exp_analysis'):
                self.on_rely([self.exp_stat, self.annotation], self.run_exp_diff)
            else:
                self.exp_stat.on('end', self.run_exp_diff)
        else:
            self.on_rely(self.final_tools, self.end)
        self.run_filecheck()
        super(DenovoBaseWorkflow, self).run()

    def end(self):
        self.logger.info('........start run denovo base end function')
        self.send_files()
        super(DenovoBaseWorkflow, self).end()

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
            ['Express_stat/', "文件夹", "表达量分析结果目录"],
            ['Express_stat/diff_exp', "文件夹", "差异基因统计结果目录"],
            ['Express_stat/rsem', "文件夹", "表达量计算分析结果目录"],
            ['Express_stat/diff_exp/diff_fpkm', "txt", "差异基因fpkm表达量矩阵"],
            ['Express_stat/diff_exp/diff_fpkm', "txt", "差异基因计数矩阵"],
            ['Express_stat/gene_correlation', "文件夹", "基因表达量样本相关性分析结果目录"],
            ['Express_stat/tran_correlation', "文件夹", "转录本表达量样本相关性分析结果目录"],
            ["Express_stat/gene_correlation/correlation_matrix.xls", "xls", "相关系数矩阵表"],
            ["Express_stat/gene_correlation/corr_col.tre", "树文件", "相关性分析树文件"],
            ["Express_stat/gene_correlation/corr_row.tre", "树文件", "相关性分析树文件"],
            ["Express_stat/gene_correlation/pca_importance.xls", "xls", "pca分析主成分解释度表"],
            ["Express_stat/gene_correlation/pca_rotation.xls", "xls", "pca分析主成分贡献度表"],
            ["Express_stat/gene_correlation/pca_sites.xls", "xls", "样本坐标表"],
            ["Express_stat/tran_correlation/correlation_matrix.xls", "xls", "相关系数矩阵表"],
            ["Express_stat/tran_correlation/corr_col.tre", "树文件", "相关性分析树文件"],
            ["Express_stat/tran_correlation/corr_row.tre", "树文件", "相关性分析树文件"],
            ["Express_stat/tran_correlation/pca_importance.xls", "xls", "pca分析主成分解释度表"],
            ["Express_stat/tran_correlation/pca_rotation.xls", "xls", "pca分析主成分贡献度表"],
            ["Express_stat/tran_correlation/pca_sites.xls", "xls", "样本坐标表"],
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
            ['Diff_express/network', "文件夹", "差异基因网络共表达分析结果目录"],
            ["Diff_express/network/all_edges.txt", "txt", "edges结果信息"],
            ["Diff_express/network/all_nodes.txt ", "txt", "nodes结果信息"],
            ["Diff_express/network/removeGene.xls ", "xls", "移除的基因信息"],
            ["Diff_express/network/removeSample.xls ", "xls", "移除的样本信息"],
            ["Diff_express/network/softPower.pdf", "pdf", "softpower相关信息"],
            ["Diff_express/network/ModuleTree.pdf", "pdf", "ModuleTree图"],
            ["Diff_express/network/eigengeneClustering.pdf", "pdf", "eigengeneClustering图"],
            ["Diff_express/network/eigenGeneHeatmap.pdf", "pdf", "eigenGeneHeatmap图"],
            ["Diff_express/network/networkHeatmap.pdf", "pdf", "networkHeatmap图"],
            ["Diff_express/network/sampleClustering.pdf", "pdf", "sampleClustering图"],
            ["Diff_express/cluster", "文件夹", "差异基因聚类分析分析结果目录"],
            ["Diff_express/cluster/hclust/", "", "层级聚类分析结果目录"],
            ["Diff_express/cluster/hclust/hc_gene_order", "txt", "按基因聚类的基因排序列表"],
            ["Diff_express/cluster/hclust/hclust_heatmap.xls", "xls", "层级聚类热图数据"],
            ["Diff_express/cluster/hclust/hc_sample_order", "txt", "按样本聚类的样本排序列表"],
            ["Diff_express/go_rich", "文件夹", "go富集分析结果目录"],
            ["Diff_express/go_regulate", "文件夹", "go统计分析结果目录"],
            ["Diff_express/kegg_rich", "文件夹", "kegg富集分析结果目录"],
            ["Diff_express/kegg_regulate", "文件夹", "kegg调控分析结果目录"],
            ["Gene_structure/orf/pfam_domain", "txt", "Pfam比对蛋白域结果信息"],
            ['QC_stat/clip_dir', "文件夹", "SE去接头后的fastq文件输出目录"],
            ["QC_stat/seqprep_dir/", "文件夹", "PE去接头后fastq文件输出目录"]
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
            [r"Express/diff_exp/.*_edgr_stat\.xls$", "xls", "edger统计结果文件"],
            [r"Express_stat/rsem/.*results$", "xls", "单样本rsem分析结果表"],
            [r"Express_stat/rsem/.*matrix$", "xls", "表达量矩阵"],
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
            [r"Diff_express/cluster/hclust/subcluster_", "xls", "子聚类热图数据"],
            [r"Diff_express/network/CytoscapeInput.*", "txt", "Cytoscape作图数据"],
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        # for i in self.get_upload_files():
        #     self.logger.info('upload file:{}'.format(str(i)))
        self.logger.info('denovo_base upload files end')
