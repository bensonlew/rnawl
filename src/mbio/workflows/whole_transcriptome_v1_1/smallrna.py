# -*- coding:utf-8 -*-
# __author__ = 'shicaiping,qinjincheng'

import json
import os
import shutil
import tarfile
import unittest
from biocluster.config import Config
from biocluster.workflow import Workflow
from mbio.packages.dna_evolution.send_email import SendEmail
from mbio.packages.ref_rna_v2.functions import tryforgood


class SmallrnaWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SmallrnaWorkflow, self).__init__(wsheet_object)
        options = [
            ## 分析对象
            # 分析文库选择
            {'name': 'lib_select', 'type': 'string', 'default': 'longRNA,smallRNA'},

            ## 基础参数设置
            # 测序质量 ['phred_33', 'phred_64']
            {'name': 'quality_score_system', 'type': 'string', 'default': 'phred_33'},
            # 生物学重复 [True, False]
            {'name': 'is_duplicate', 'type': 'bool', 'default': True},
            # 测序类型 ['SE', 'PE']
            {'name': 'fq_type', 'type': 'string', 'default': 'SE'},
            # 测序读长
            {'name': 'read_length', 'type': 'int', 'default': 150},
            # 建库类型 ['type I', 'type II', 'type III', 'other']
            {'name': 'lib_type', 'type': 'string', 'default': 'type I'},
            # 去除5'端前 碱基
            {'name': 'cut_5', 'type': 'int', 'default': 4},
            # 提取前 进行分析
            {'name': 'extract_length', 'type': 'int', 'default': 75},
            # 去除3'端 碱基
            {'name': 'cut_3', 'type': 'int', 'default': 4},
            # 接头序列
            {'name': 'adapter', 'type': 'string', 'default': 'TGGAATTCTCGGGTGCCAAGG'},
            # 原始序列文件
            {'name': 'fastq_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},
            # 分组方案
            {'name': 'group_table', 'type': 'infile', 'format': 'sample.group_table'},
            # 对照组文件
            {'name': 'control_file', 'type': 'infile', 'format': 'sample.control_table'},
            # 配对信息表
            {'name': 'pair_table', 'type': 'infile', 'format': 'sample.group_table'},
            # 物种类别 ['Animal', 'Plant', 'Protist', 'Fungi']
            {'name': 'taxonmy', 'type': 'string', 'default': 'Animal'},
            # 具体物种 sg_genome_db.organism_name
            {'name': 'organism_name', 'type': 'string', 'default': None},
            # 基因组编号 sg_genome_db.genome_id
            {'name': 'genome_id', 'type': 'string', 'default': None},
            # miRBase 一级分类 ['Animal', 'Plant']
            {'name': 'mirbase_category', 'type': 'string', 'default': 'Animal'},
            # miRBase 具体物种（多选时分号分隔）
            {'name': 'mirbase_specie', 'type': 'string', 'default': None},
            # 已知miRNA鉴定允许错配
            {'name': 'mismatch', 'type': 'int', 'default': 1},

            # 差异分析软件 ['DESeq2', 'edgeR', 'DEGseq']
            {'name': 'diff_method', 'type': 'string', 'default': 'DESeq2'},
            # 表达量筛选
            {'name': 'diff_filter', 'type': 'string', 'default': None},
            {'name': 'diff_threshold', 'type': 'float', 'default': 0.0},
            # FC
            {'name': 'fc', 'type': 'float', 'default': 2.0},
            # 显著性水平
            {'name': 'pvalue_padjust', 'type': 'string', 'default': 'padjust'},
            {'name': 'diff_fdr_ci', 'type': 'float', 'default': 0.05},
            # 多重验证校正方法 ['BH', 'Bonferroni', 'Holm', 'BY']
            {'name': 'padjust_way', 'type': 'string', 'default': 'BH'},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.task_id = self._sheet.id
        self.project_sn = self._sheet.project_sn
        self.modules = dict()
        self.tools = dict()

        database = Config().get_mongo_client(mtype='ref_rna_v2', dydb_forbid=True)[Config().get_mongo_dbname('ref_rna_v2', dydb_forbid=True)]
        collection = database['sg_genome_db']
        self.genome_doc = collection.find_one({'genome_id': self.option('genome_id')})
        self.db_path = os.path.join(self.config.SOFTWARE_DIR, 'database/Genome_DB_finish')

        # 用于在重运行时，删除已经导入到mongo库的表，避免数据重复
        data = os.path.join(self.work_dir, 'data.json')
        if os.path.exists(data):
            with open(data, 'r') as load_f:
                load_dict = json.load(load_f)
                if 'rerun' in load_dict and load_dict['rerun']:
                    self.logger.info("该项目重运行中，先删除mongo库中已有数据")
                    self.delete_mongo_data()

    @tryforgood
    def delete_mongo_data(self):
        self.script = os.path.join(self.config.PACKAGE_DIR, 'project_demo/delete_demo.py')
        self.program = os.path.join(self.config.SOFTWARE_DIR, 'program/Python/bin/python')
        cmd = '{} {}'.format(self.program, self.script)
        cmd += ' {} {}'.format(self.task_id, 'whole_transcriptome')
        code = os.system(cmd)
        if code == 0:
            self.logger.info("命令{}执行成功！".format(cmd))
        else:
            raise Exception("命令{}执行失败！".format(cmd))

    def check_options(self):
        for k, v in self.sheet.options().items():
            self.logger.debug('{} -> {}'.format(k, v))
        else:
            return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        self.add_steps()
        super(SmallrnaWorkflow, self).run()

    def add_steps(self):
        self.step.add_steps('file_check')
        self.step.add_steps('gzfastq2fastq')
        self.step.add_steps('mirna_qc')
        self.step.add_steps('hiseq_reads_stat_raw')
        self.step.add_steps('hiseq_reads_stat_use')
        self.step.add_steps('fasta_uniq')
        self.step.add_steps('mapper_and_stat')
        self.step.add_steps('srna')
        self.step.add_steps('diff_exp')
        self.load_libraries()

    def load_libraries(self):
        self.tools['file_check'] = self.add_tool('whole_transcriptome.smallrna.file_check')
        self.tools['gzfastq2fastq'] = self.add_tool('whole_transcriptome.smallrna.gzfastq2fastq')
        self.modules['mirna_qc'] = self.add_module('whole_transcriptome.smallrna.mirna_qc')
        self.modules['hiseq_reads_stat_raw'] = self.add_module('whole_transcriptome.smallrna.hiseq_reads_stat')
        self.modules['hiseq_reads_stat_use'] = self.add_module('whole_transcriptome.smallrna.hiseq_reads_stat')
        self.tools['fasta_uniq'] = self.add_tool('whole_transcriptome.smallrna.fasta_uniq')
        self.tools['mapper_and_stat'] = self.add_tool('whole_transcriptome.smallrna.mapper_and_stat')
        self.modules['srna'] = self.add_module('whole_transcriptome.smallrna.srna')
        self.modules['diff_exp'] = self.add_module('whole_transcriptome.diff_exp')
        self.on_rely([
            self.modules['hiseq_reads_stat_raw'],
            self.modules['hiseq_reads_stat_use'],
            self.modules['diff_exp']
        ], self.set_db)
        self.run_file_check()

    def run_file_check(self):
        fastq_dir = self.option('fastq_dir').path
        opts = {
            'fq_type': self.option('fq_type'),
            'fastq_dir': fastq_dir,
            'group_table': self.option('group_table'),
            'control_file': self.option('control_file')
        }
        self.tools['file_check'].set_options(opts)
        self.tools['file_check'].on('start', self.set_step, {'start': self.step.file_check})
        self.tools['file_check'].on('end', self.set_step, {'end': self.step.file_check})
        self.tools['file_check'].on('end', self.set_output, 'file_check')
        self.tools['file_check'].on('end', self.run_gzfastq2fastq)
        self.tools['file_check'].run()

    def run_gzfastq2fastq(self):
        fastq_path = self.option('fastq_dir').path
        opts = {'fastq_path': fastq_path}
        self.tools['gzfastq2fastq'].set_options(opts)
        self.tools['gzfastq2fastq'].on('start', self.set_step, {'start': self.step.gzfastq2fastq})
        self.tools['gzfastq2fastq'].on('end', self.set_step, {'end': self.step.gzfastq2fastq})
        self.tools['gzfastq2fastq'].on('end', self.set_output, 'gzfastq2fastq')
        self.tools['gzfastq2fastq'].on('end', self.run_mirna_qc)
        self.tools['gzfastq2fastq'].run()

    def run_mirna_qc(self):
        list_file = os.path.join(self.option('fastq_dir').path, 'list_re.txt')
        adapter = self.option('adapter')
        if self.option('quality_score_system').endswith('33'):
            fastq_format = 'Q33'
        elif self.option('quality_score_system').endswith('64'):
            fastq_format = 'Q64'
        cut_left = self.option('cut_5')
        cut_tail = self.option('cut_3')
        extract_length = self.option('extract_length')
        opts = {
            'list_file': list_file,
            'adapter': adapter,
            'fastq_format': fastq_format,
            'cut_left': cut_left,
            'cut_tail': cut_tail,
            'extract_length': extract_length
        }
        self.modules['mirna_qc'].set_options(opts)
        self.modules['mirna_qc'].on('start', self.set_step, {'start': self.step.mirna_qc})
        self.modules['mirna_qc'].on('end', self.set_step, {'end': self.step.mirna_qc})
        self.modules['mirna_qc'].on('end', self.set_output, 'mirna_qc')
        self.modules['mirna_qc'].on('end', self.run_fasta_uniq)
        self.modules['mirna_qc'].on('end', self.run_hiseq_reads_stat_raw)
        self.modules['mirna_qc'].on('end', self.run_hiseq_reads_stat_use)
        self.modules['mirna_qc'].run()

    def run_hiseq_reads_stat_raw(self):
        fastq_dir = self.modules['mirna_qc'].option('rawdata')
        fq_type = 'SE'
        if self.option('quality_score_system').endswith('33'):
            quality = 33
        elif self.option('quality_score_system').endswith('64'):
            quality = 64
        dup = False
        opts = {
            'fastq_dir': fastq_dir,
            'fq_type': fq_type,
            'quality': quality,
            'dup': dup
        }
        self.modules['hiseq_reads_stat_raw'].set_options(opts)
        self.modules['hiseq_reads_stat_raw'].on('start', self.set_step, {'start': self.step.hiseq_reads_stat_raw})
        self.modules['hiseq_reads_stat_raw'].on('end', self.set_step, {'end': self.step.hiseq_reads_stat_raw})
        self.modules['hiseq_reads_stat_raw'].on('end', self.set_output, 'hiseq_reads_stat_raw')
        self.modules['hiseq_reads_stat_raw'].run()

    def run_hiseq_reads_stat_use(self):
        fastq_dir = self.modules['mirna_qc'].option('cleandata')
        fq_type = 'SE'
        if self.option('quality_score_system').endswith('33'):
            quality = 33
        elif self.option('quality_score_system').endswith('64'):
            quality = 64
        dup = True
        opts = {
            'fastq_dir': fastq_dir,
            'fq_type': fq_type,
            'quality': quality,
            'dup': dup
        }
        self.modules['hiseq_reads_stat_use'].set_options(opts)
        self.modules['hiseq_reads_stat_use'].on('start', self.set_step, {'start': self.step.hiseq_reads_stat_use})
        self.modules['hiseq_reads_stat_use'].on('end', self.set_step, {'end': self.step.hiseq_reads_stat_use})
        self.modules['hiseq_reads_stat_use'].on('end', self.set_output, 'hiseq_reads_stat_use')
        self.modules['hiseq_reads_stat_use'].run()

    def run_fasta_uniq(self):
        config = self.modules['mirna_qc'].option('config_file').path
        opts = {'config': config}
        self.tools['fasta_uniq'].set_options(opts)
        self.tools['fasta_uniq'].on('start', self.set_step, {'start': self.step.fasta_uniq})
        self.tools['fasta_uniq'].on('end', self.set_step, {'end': self.step.fasta_uniq})
        self.tools['fasta_uniq'].on('end', self.set_output, 'fasta_uniq')
        self.tools['fasta_uniq'].on('end', self.run_mapper_and_stat)
        self.tools['fasta_uniq'].run()

    def run_mapper_and_stat(self):
        config = self.modules['mirna_qc'].option('config_file').path
        ref = os.path.join(self.db_path, self.genome_doc['dna_fa'])
        fasta = os.path.join(self.tools['fasta_uniq'].output_dir, 'uniq.fasta')
        index = os.path.join(self.db_path, self.genome_doc['dna_index'])
        gtf = os.path.join(self.db_path, self.genome_doc['gtf'])
        opts = {
            'config': config,
            'ref': ref,
            'fasta': fasta,
            'index': index,
            'gtf': gtf
        }
        self.tools['mapper_and_stat'].set_options(opts)
        self.tools['mapper_and_stat'].on('start', self.set_step, {'start': self.step.mapper_and_stat})
        self.tools['mapper_and_stat'].on('end', self.set_step, {'end': self.step.mapper_and_stat})
        self.tools['mapper_and_stat'].on('end', self.set_output, 'mapper_and_stat')
        self.tools['mapper_and_stat'].on('end', self.run_srna)
        self.tools['mapper_and_stat'].run()

    def run_srna(self):
        category = self.option('mirbase_category')
        species = self.option('mirbase_specie')
        reference = os.path.join(self.db_path, self.genome_doc['dna_fa'])
        clean_fa = os.path.join(self.tools['fasta_uniq'].output_dir, 'uniq.fasta')
        config = self.modules['mirna_qc'].option('config_file').path
        arf = os.path.join(self.tools['mapper_and_stat'].output_dir, 'reads_vs_genome.arf')
        gtf = os.path.join(self.db_path, self.genome_doc['gtf'])
        qc_output = self.modules['mirna_qc'].output_dir
        mismatch = self.option('mismatch')
        opts = {
            'category': category,
            'species': species,
            'reference': reference,
            'clean_fa': clean_fa,
            'config': config,
            'arf': arf,
            'gtf': gtf,
            'list': os.path.join(self.option('fastq_dir').path, 'list_re.txt'),
            'qc_output': qc_output,
            'mismatch': mismatch,
        }
        repeatmasker_dir = os.path.join(self.db_path, self.genome_doc['anno_path_v2'], 'repeatmasker')
        if os.path.isdir(repeatmasker_dir):
            repeat = True
            repeat_fa = os.path.join(repeatmasker_dir, 'repeatmasker_merge.repeat.fa')
            repeat_gff = os.path.join(repeatmasker_dir, 'repeatmasker_merge.SSR.gff')
            opts.update({'repeat': repeat, 'repeat_fa': repeat_fa, 'repeat_gff': repeat_gff})
        self.modules['srna'].set_options(opts)
        self.modules['srna'].on('start', self.set_step, {'start': self.step.srna})
        self.modules['srna'].on('end', self.set_step, {'end': self.step.srna})
        self.modules['srna'].on('end', self.set_output, 'srna')
        self.modules['srna'].on('end', self.run_diff_exp)
        self.modules['srna'].run()

    def run_diff_exp(self):
        program = self.option('diff_method')
        count_matrix = os.path.join(self.modules['srna'].output_dir, 'total_mirna_count.xls')
        group_table = self.option('group_table').path
        control_table = self.option('control_file').path
        exp_matrix = os.path.join(self.modules['srna'].output_dir, 'total_mirna_norm.xls')
        kind_table = os.path.join(self.modules['diff_exp'].work_dir, 'kind.txt')
        lines = ['transcript_id\tkind\n']
        lines.extend('{}\t{}\n'.format(line.split('\t')[0], 'ref') for i, line in
                     enumerate(open(os.path.join(self.modules['srna'].output_dir, 'known_mirna_count.xls'))) if i)
        lines.extend('{}\t{}\n'.format(line.split('\t')[0], 'new') for i, line in
                     enumerate(open(os.path.join(self.modules['srna'].output_dir, 'novel_mirna_count.xls'))) if i)
        open(kind_table, 'w').writelines(lines)
        threshold = self.option('diff_threshold')
        method = self.option('padjust_way').lower()
        stat_type = self.option('pvalue_padjust')
        stat_cutoff = self.option('diff_fdr_ci')
        fc = self.option('fc')
        opts = {
            'program': program,
            'count_matrix': count_matrix,
            'group_table': group_table,
            'control_table': control_table,
            'exp_matrix': exp_matrix,
            'kind_table': kind_table,
            'filter': str(self.option('diff_filter')).lower(),
            'threshold': threshold,
            'fc': fc
        }
        if program.lower() in ["degseq", "edger", "deseq2", 'limma']:
            opts.update({
                'method': method,
                'stat_type': stat_type,
                'stat_cutoff': stat_cutoff,
            })
            if self.option('pair_table').is_set:
                opts.update({'is_batch': True, 'has_batch': True, 'batch_matrix': self.option('pair_table').path})
        if program.lower() == 'noiseq':
            opts.update({
                stat_type: float(stat_cutoff),
            })
            if self.option('pair_table').is_set:
                opts.update({'is_batch': True, 'has_batch': True, 'batch_matrix': self.option('pair_table').path})
        self.modules['diff_exp'].set_options(opts)
        self.modules['diff_exp'].on('start', self.set_step, {'start': self.step.diff_exp})
        self.modules['diff_exp'].on('end', self.set_step, {'end': self.step.diff_exp})
        self.modules['diff_exp'].on('end', self.set_output, 'diff_exp')
        self.modules['diff_exp'].run()

    def set_output(self, event):
        obj = event['bind_object']
        self.move2outputdir(obj.output_dir, event['data'])

    def move2outputdir(self, olddir, newname):
        newdir = os.path.join(self.output_dir, newname)
        if not os.path.exists(newdir):
            os.makedirs(newdir)
        allfiles = os.listdir(olddir)
        oldfiles = [os.path.join(olddir, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.isfile(newfile) and os.path.exists(newfile):
                os.remove(newfile)
            elif os.path.isdir(newfile) and os.path.exists(newfile):
                shutil.rmtree(newfile)
        for i in range(len(allfiles)):
            self.move_file(oldfiles[i], newfiles[i])
        else:
            self.logger.debug('succeed in linking {} to {}'.format(olddir, newdir))

    def move_file(self, src, dst):
        if os.path.isfile(src):
            os.link(src, dst)
        else:
            os.mkdir(dst)
            for file in os.listdir(src):
                old_path = os.path.join(src, file)
                new_path = os.path.join(dst, file)
                self.move_file(old_path, new_path)

    def set_db(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.export_task_info()
        self.export_genome_info()
        self.export_qc()
        self.export_mapping()
        self.export_srna()
        self.export_diff_exp_stat()
        self.export_email()
        self.set_upload()

    def set_upload(self):
        shutil.copy(os.path.join(self.work_dir, 'data.json'), os.path.join(self.output_dir, 'data.json'))
        shutil.copy(self.option('group_table').path, os.path.join(self.output_dir, 'group.txt'))
        shutil.copy(self.option('control_file').path, os.path.join(self.output_dir, 'control.txt'))
        shutil.rmtree(os.path.join(self.output_dir, 'mirna_qc/raw_data'))
        self.pack_and_compress()
        self.add_upload_dir(self.output_dir)
        self.end()

    def pack_and_compress(self):
        self.logger.warn('start packing folders and compressing them, it will take a long time')
        os.chdir(self.output_dir)
        dirs = os.walk(self.output_dir).next()[1]
        for dirname in dirs:
            tarname = '{}.tar.gz'.format(dirname)
            if dirname in ('file_check', 'gzfastq2fastq'):
                continue
            if not os.path.isfile(os.path.join(self.output_dir, tarname)):
                with tarfile.open(tarname, 'w:gz') as tar:
                    tar.add(dirname)
            shutil.rmtree(dirname)
        os.chdir(self.work_dir)
        self.logger.warn('succeed in packing folders and compressing them')

    def end(self):
        super(SmallrnaWorkflow, self).end()

    def export_task_info(self):
        api = self.api.api('whole_transcriptome.task_info')
        api.add_task_info(os.path.join(self.work_dir, 'data.json'))

    def export_genome_info(self):
        api = self.api.api('whole_transcriptome.genome_info')
        file_path = os.path.join(self.db_path, self.genome_doc['gene_stat'])
        species_name = self.option('organism_name')
        ref_anno_version = self.genome_doc['assembly']
        hyperlink = self.genome_doc['ensemble_web']
        api.add_genome_info(file_path, species_name, ref_anno_version, hyperlink)

    def export_qc(self):
        api = self.api.api('whole_transcriptome.qc')
        sample_list = os.path.join(self.option('fastq_dir').path, 'list.txt')
        group_file = self.option('group_table').path
        compare_file = self.option('control_file').path
        qc_stat_before = self.modules['hiseq_reads_stat_raw'].output_dir
        qc_stat_after = self.modules['hiseq_reads_stat_use'].output_dir
        qc_result_dir = os.path.join(self.modules['mirna_qc'].output_dir, 'clean_data')
        api.add_sample_info(sample_list=sample_list, library='small')
        group_id, specimen_names, category_names = api.add_sample_group(group_file=group_file, library='small')
        control_id, compare_names = api.add_group_compare(compare_file=compare_file, library='small', group_id=group_id)
        qc_id = api.add_qc(fq_type='PE', library='small')
        api.add_qc_detail(qc_id, qc_stat_before, qc_stat_after, 'small', qc_result_dir)
        api.add_qc_graph(qc_id, qc_stat_before, 'small', 'before')
        api.add_qc_graph(qc_id, qc_stat_after, 'small', 'after', qc_result_dir)

    def export_mapping(self):
        api = self.api.api('whole_transcriptome.mapping')
        stat_file = self.tools['mapper_and_stat'].output_dir
        distribution = self.tools['mapper_and_stat'].output_dir
        sample_list_file = os.path.join(self.modules['mirna_qc'].output_dir, 'clean_data/list.txt')
        method = 'bowtie'
        params = json.dumps({'task_id': self.task_id, 'submit_location': 'mapping', 'task_type': 2}, sort_keys=True,
                            separators=(',', ':'))
        api.add_mapping_stat(stat_file=stat_file, library='small', method=method, sample_list_file=sample_list_file)
        api.add_chrom_distribution_table(distribution=distribution, params=params, library='small',
                                         sample_list_file=sample_list_file)

    def export_srna(self):
        api = self.api.api('whole_transcriptome.srna')
        mirna_stat = os.path.join(self.modules['srna'].output_dir, 'srna_stat/mirna_stat.xls')
        srna_stat = os.path.join(self.modules['srna'].output_dir, 'srna_stat/srna_stat.xls')
        srna_stat_for_graph = os.path.join(self.modules['srna'].output_dir, 'srna_stat/srna_stat_for_graph.xls')
        params = json.dumps(
            {'task_id': self.task_id, 'submit_location': 'srnastat', 'task_type': 2, 'method': 'mirdeep'},
            sort_keys=True, separators=(',', ':'))
        api.add_mirna_stat(mirna_stat=mirna_stat, task_id=self.task_id, project_sn=self.project_sn, params=params)
        api.add_srna_stat(srna_stat=srna_stat, srna_stat_for_graph=srna_stat_for_graph, task_id=self.task_id,
                          project_sn=self.project_sn, params=params)

    def export_diff_exp_stat(self):
        api = self.api.api('whole_transcriptome.diff_exp_stat')
        map_dict = {'control': self.option('control_file').path, 's_output_dir': self.modules['diff_exp'].output_dir}
        if self.option('diff_method').lower() in ['degseq', 'deseq2', 'edger', 'limma']:
            arg_dict = {
                'program': self.option('diff_method'),
                'fc': self.option('fc'),
                'qvalue': self.option('diff_fdr_ci'),
                'method': self.option('padjust_way')
            }
        if self.option('diff_method').lower() in ['noiseq']:
            arg_dict = {
                'program': self.option('diff_method'),
                'fc': self.option('fc'),
                'prob': self.option('diff_fdr_ci'),
            }
        api.add_diff_exp_stat(map_dict, arg_dict, self.task_id, self.project_sn, 'small')

    def export_email(self):
        if self.modules['diff_exp'].option('email'):
            receiver = ['meng.luo@majorbio.com', 'chunxiang.xue@majordx.com', 'jiameng.li@majorbio.com',
                        'daokuan.zhang@majorbio.com', 'dongmei.fu@majorbio.com', 'rna@majorbio.com']
            a = SendEmail("1274095891@qq.com", "smtp.qq.com", "ocfnjnhnsbqubaej", "1274095891@qq.com",
                          ",".join(receiver), "WARNING - project_sn ({}), task_id ({})".format(
                              self._sheet.project_sn, self._sheet.id), 465)
            a.send_msg('The number of differentially expressed miRNAs is less than 10.')
            a.send_email()


class TestFunction(unittest.TestCase):
    """
    This is test for the workflow. Just run this script to do test.
    """

    def test_ath(self):
        from mbio.workflows.whole_transcriptome.smallrna import SmallrnaWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.smallrna',
            'options': {
                'quality_score_system': 'phred_33',
                'is_duplicate': 'True',
                'fq_type': 'SE',
                'read_length': 50,
                'lib_type': 'other',
                'cut_5': 0,
                'extract_length': 50,
                'cut_3': 0,
                'adapter': 'TGGAATTCTCGGGTGCCAAGG',
                'fastq_dir': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/whole_transcriptome/ath/raw_data',
                'group_table': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/whole_transcriptome/ath/group.txt',
                'control_file': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/whole_transcriptome/ath/control.txt',
                'taxonmy': 'Plant',
                'organism_name': 'Arabidopsis_thaliana',
                'genome_id': 'GM0348',
                'mirbase_category': 'Plant',
                'mirbase_specie': 'ath',
                'mismatch': 0,
            }
        }
        wsheet = Sheet(data=data)
        wf = SmallrnaWorkflow(wsheet)
        wf.sheet.id = 'smallrna'
        wf.sheet.project_sn = 'smallrna'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()

    def test(self):
        from mbio.workflows.whole_transcriptome.smallrna import SmallrnaWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.smallrna',
            'options': {
                'quality_score_system': 'phred_33',
                'is_duplicate': 'True',
                'fq_type': 'SE',
                'read_length': 150,
                'lib_type': 'type I',
                'cut_5': 4,
                'extract_length': 75,
                'cut_3': 4,
                'adapter': 'TGGAATTCTCGGGTGCCAAGG',
                'fastq_dir': '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/whole_transcriptome/demo_human/smallrna/rawdata',
                'group_table': '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/whole_transcriptome/demo_human/group.txt',
                'control_file': '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/whole_transcriptome/demo_human/control.txt',
                'taxonmy': 'Animal',
                'organism_name': 'Homo_sapiens',
                'genome_id': 'GM0259',
                'mirbase_category': 'Animal',
                'mirbase_specie': 'hsa',
                'mismatch': 1,
            }
        }
        wsheet = Sheet(data=data)
        wf = SmallrnaWorkflow(wsheet)
        wf.sheet.id = 'smallrna'
        wf.sheet.project_sn = 'smallrna'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_ath')])
    unittest.TextTestRunner(verbosity=2).run(suite)
