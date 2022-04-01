# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# __last_modified__ = '20190119'

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError, FileError
from bson import ObjectId
from biocluster.config import Config
import os,re
import json
import shutil
import time
import datetime
import gevent
import functools
from mbio.packages.metagbin.common_function import link_dir
from mbio.packages.meta.delete_mongo import DeleteDemoMongo

def time_count(func):  # 统计导表时间
    @functools.wraps(func)
    def wrapper(self, *args, **kw):
        start = time.time()
        func_name = func.__name__
        start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))
        self.main_collection_delete(func_name)
        func(self, *args, **kw)
        end = time.time()
        end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end))
        print('End %s at %s' % (func_name, end_time))
        print("{}函数执行完毕，该阶段导表已进行{}s".format(func_name, end - start))
    return wrapper

def tryforgood(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except:
            return wrapper(*args, **kwargs)
    return wrapper

class MetagbinWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        宏基因组binning参数设置
        """
        self._sheet = wsheet_object
        super(MetagbinWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'in_fastq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq文件夹
            {'name': 'specimen_info', 'type': 'infile', 'format': 'meta_genomic.specimen_info'},  # 样本集信息表
            {'name': 'qc', 'type': 'bool', 'default': False},  # 是否需要质控
            {'name': 'qc_quality', 'type': 'int', 'default': 20},  # 质控质量值标准
            {'name': 'qc_length', 'type': 'int', 'default': 30},  # 质控最短序列长度
            {'name': 'rm_host', 'type': 'bool', 'default': False},  # 是否需要去除宿主
            {'name': 'ref_database', 'type': 'string', 'default': ''},  # 宿主参考序列库中对应的物种名，eg：E.coli ,B.taurus
            {'name': 'ref_undefined', "type": 'infile', 'format': 'sequence.fasta_dir'},
            {'name': 'ref_undefined_name', 'type': 'string', 'default': 'undefined'},  # 自定义参考宿主名称，适应页面参数
            {'name': 'cdhit_identity', 'type': 'float', 'default': 0.95},  # 基因序列聚类相似度
            {'name': 'cdhit_coverage', 'type': 'float', 'default': 0.9},  # 基因序列聚类覆盖度
            {'name': 'min_contig', 'type': 'string', 'default': "1000"},  # 基因序列聚类覆盖度
            {'name': 'sofware_bin', 'type': 'string', 'default': 'metabat'},  # bin运行软件方法
            {"name": "assemble_dir", "type": "infile", "format": "metagbin.assemble_dir"},  # 组装文件fasta文件夹
            {'name': 'clean_fq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq文件夹
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.sequence = self.add_module('sequence.meta_genomic')
        self.qc = self.add_module('meta.qc.qc_and_stat')
        self.rm_host = self.add_module('meta.qc.bwa_remove_host')
        self.ass_seq = self.add_module('metagbin.metagbin_assemble_seq')
        self.bin_tax=self.add_module("annotation.bin_tax")
        self.clean_data = self.add_module("sequence.metagbin_clean_fq")
        self.read_mapping = self.add_module('metagbin.metagbin_mapping')
        self.scaf_abund = self.add_module("metagbin.scaffold_abund")
        self.spring = self.add_module("metagbin.spring_files")
        self.bin = self.add_module('metagbin.bin')
        self.assemble = self.add_module('metagbin.bin_assembly')
        self.predict = self.add_module('metagbin.stat_predict')
        self.anno = self.add_module('metagbin.annotation')
        self.bin_stat = self.add_module("metagbin.metagbin_bin_stat")
        self.list = [self.clean_data, self.ass_seq]
        self.list2 = [self.sequence, self.qc, self.rm_host, self.ass_seq]
        self.list3 = [self.sequence, self.qc, self.ass_seq]
        self.step.add_steps('sequence', "spring", 'ass_seq', 'qc_', 'anno', 'bin', 'scaf_abund', 'read_mapping', 'clean_data', 'bin_tax', 'assemble', 'predict', 'bin_stat',"rm_host",'anno_1')
        self.bin_name = ''
        self.remote_dir = self._sheet.output

        try:
            self.rerun = self._sheet.rerun
        except:
            self.rerun = False
        if self.rerun:
            self.logger.info("该项目重运行中，先删除mongo库中已有数据")
            self.delete_mongo_data()

    @tryforgood
    def delete_mongo_data(self):
        delete = DeleteDemoMongo(self._sheet.id, 'metagbin')
        try:
            delete.run()
            self.logger.info("删库已完成")
        except:
            raise Exception("删除记录失败")

    def check_options(self):
        """
        检查参数设置
        """
        self.logger.info("check_sheet_data...")
        if not self.option('assemble_dir').is_set:
            raise OptionError('请输入组装序列文件夹！')
        if not self.option('in_fastq').is_set and not self.option('clean_fq').is_set:
            raise OptionError('请输入原始序列文件夹或者质控优质序列文件夹！')
        if not self.option('specimen_info').is_set:
            raise OptionError('质控需提供样本集信息表！')
        if self.option('rm_host'):
            if self.option('ref_database') == '' and not self.option('ref_undefined').is_set:
                raise OptionError('已选择去宿主，需输入参考数据库或参考序列')
            if self.option('ref_database') not in ['', 'Custom'] and self.option('ref_undefined').is_set:
                raise OptionError('去宿主不可同时提供参考数据库及参考序列')
        if self.option('min_contig') < 1000:
            raise OptionError('最小Contig长度参数超出范围1000~!')
        if not 0.75 <= self.option("cdhit_identity") <= 1:
            raise OptionError("cdhit identity必须在0.75，1之间")
        if not 0 <= self.option("cdhit_coverage") <= 1:
            raise OptionError("cdhit coverage必须在0,1之间")
        return True

    def check_specimen_names(self):
        error_message = ""
        if self.option('specimen_info').is_set and self.option('in_fastq').is_set:
            list_set = self.read_file_to_set(self.option('in_fastq').prop['path'] + '/list.txt', header=False, list_index=1)
            specimen_info_set = self.read_file_to_set(self.option('specimen_info').prop['path'], header=True)
            if list_set == specimen_info_set:
                self.logger.info("specimen_info sample check pass!!!")
            else:
                self.logger.error("\n specimen_info sample check has error:\n\t in_fastq samples: \t%s, \n\tspecimen_info samples: \t%s" % (list_set, specimen_info_set))
                error_message += "specimen_info, "
        elif self.option('specimen_info').is_set and self.option('clean_fq').is_set:
            list_set = self.read_file_to_set(self.option('clean_fq').prop['path'] + '/list.txt', header=False,list_index=1)
            specimen_info_set = self.read_file_to_set(self.option('specimen_info').prop['path'], header=True)
            if list_set == specimen_info_set:
                self.logger.info("specimen_info sample check pass!!!")
            else:
                self.logger.error("\n specimen_info sample check has error:\n\t clean_fq samples: \t%s, \n\tspecimen_info samples: \t%s" % (list_set, specimen_info_set))
                error_message += "specimen_info, "
        return error_message

    def read_file_to_set(self, file, header=True, list_index=0):
        specimen_set = set()
        try:
            list_index = int(list_index)
        except:
            raise OptionError("list_index is not integer!")
        with open(file, 'r') as rf:
            lines = rf.readlines()
            if header:
                lines.pop(0)
            for line in lines:
                t = line.split('\t')
                if len(t) <= list_index:
                    raise OptionError("%s 只有 %s 列, 不存在第 %s 列", variables=(file, len(t), list_index + 1))
                specimen_set.add(t[list_index])
        return specimen_set

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def set_run(self, opts, module, event, step, start=True):
        module.set_options(opts)
        module.on('start', self.set_step, {'start': step})
        module.on('end', self.set_step, {'end': step})
        module.on('end', self.set_output, event)
        if start:
            module.run()

    def run(self):
        """
        运行 metagbin workflow
        :return:
        """
        task_info = self.api.api('task_info.metagbin_task_info')
        task_info.add_task_info()
        self.logger.info(self.option('qc'))
        if self.option('qc'):
            self.sequence.on('end', self.run_qc)
            self.run_sequence()
            self.run_assemble_seq()
            if self.option('rm_host'):
                self.qc.on('end', self.run_rm_host)
                self.on_rely(self.list2, self.run_mapping)
            else:
                self.on_rely(self.list3, self.run_mapping)
            self.read_mapping.on('end', self.run_scaffold_cov)
            self.scaf_abund.on('end', self.run_bin)
            self.bin.on('end', self.run_bin_stat)
            self.bin_stat.on('end', self.run_genome_ass)
            self.assemble.on('end', self.run_spring_data)
            self.spring.on("end", self.run_gene_pre)
            self.predict.on('end', self.run_anno)
        else:
            self.run_assemble_seq()
            self.run_clean_sequence()
            self.on_rely(self.list,self.run_mapping)
            self.read_mapping.on('end', self.run_scaffold_cov)
            self.scaf_abund.on('end', self.run_bin)
            self.bin.on('end', self.run_bin_stat)
            self.bin_stat.on('end', self.run_genome_ass)
            self.assemble.on('end', self.run_spring_data)
            self.spring.on("end", self.run_gene_pre)
            self.predict.on('end', self.run_anno)
        super(MetagbinWorkflow, self).run()

    def run_sequence(self):
        opts = {
            'fastq_dir': self.option('in_fastq'),
            'fastp':"fastp",
        }
        self.set_run(opts, self.sequence, 'sequence', self.step.sequence)

    def run_qc(self):
        opts = {
            'fastq_dir': self.sequence.output_dir + '/data',
            'insert_size': self.option('specimen_info').prop['path'],
            'fastp': "fastp",
        }
        self.set_run(opts, self.qc, 'qc', self.step.qc_)

    def run_rm_host(self):
        opts = {
            'fq_type': 'PE',
            'ref_database': self.option('ref_database'),
            'ref_undefined': self.option('ref_undefined'),
        }
        if opts['ref_database'] == 'Custom':
            opts['ref_database'] = ""
        if self.option('qc'):
            opts['fastq_dir'] = self.qc.option('after_qc_dir')
        else:
            opts['fastq_dir'] = self.option('in_fastq')
        self.set_run(opts, self.rm_host, 'rm_host', self.step.rm_host)

    def run_clean_sequence(self):
        opts = {
            'fastq_dir': self.option('clean_fq'),
        }
        self.set_run(opts, self.clean_data, 'clean_data', self.step.clean_data)

    def run_assemble_seq(self):
        opts = {
            'assemble_dir': self.option('assemble_dir'),
            'cdhit_identity': self.option('cdhit_identity'),
            'cdhit_coverage': self.option('cdhit_coverage'),
        }
        self.set_run(opts, self.ass_seq, 'ass_seq', self.step.ass_seq)

    def run_mapping(self):
        fastq_dir = ''
        if self.option("qc"):
            if self.option('rm_host'):
                fastq_dir = self.rm_host.output_dir
            else:
                fastq_dir = self.qc.output_dir + '/after_qc_dir'
        else:
            fastq_dir = self.clean_data.output_dir + '/data'
        asse_seq = self.ass_seq.output_dir + '/all.uniContigs.fa'
        opts = {
            'ref_fa': asse_seq,
            'fastq_dir': fastq_dir,
        }
        self.set_run(opts, self.read_mapping, 'read_mapping', self.step.read_mapping)

    def run_scaffold_cov(self):
        opts = {
            "bam_dir":self.read_mapping.option('bam_dir'),
        }
        self.set_run(opts, self.scaf_abund, 'scaf_abund', self.step.scaf_abund)

    def run_bin(self):
        opts = {
            "contig_fa": self.ass_seq.output_dir + '/all.uniContigs.fa',
            "sofware_bin": self.option('sofware_bin'),
            "minContig": self.option("min_contig")
        }
        if re.search(r'metabat', self.option('sofware_bin')):
            opts['metabat_depth'] = self.scaf_abund.option('metabat_depth')
        if re.search(r'maxbin', self.option('sofware_bin')):
            opts['maxbin_depth'] = self.scaf_abund.option('maxbin_depth')
        if re.search(r'concoct', self.option('sofware_bin')):
            opts['specimen_info'] = self.option('specimen_info')
            opts['bam_dir'] = self.read_mapping.option('bam_dir')
        self.set_run(opts, self.bin, 'bin', self.step.bin)

    def run_bin_stat(self):
        opts = {
            'bin_dir': self.bin.output_dir + '/bin',
            'input_genome': self.ass_seq.output_dir + '/all.uniContigs.fa',
            'metabat_depth': self.scaf_abund.option('metabat_depth'),
        }
        self.set_run(opts, self.bin_stat, 'bin_stat', self.step.bin_stat)

    def run_genome_ass(self):
        self.bin_name = self.get_binname(self.bin_stat.output_dir + '/all.bin.summary.xls')
        if self.option('qc'):
            if self.option('rm_host'):
                clean_fq = self.rm_host.option('result_fq_dir')
            else:
                clean_fq = self.qc.option('after_qc_dir')
        else:
            clean_fq = self.clean_data.output_dir + '/data'
        opts = {
            "ref_fa": self.bin.output_dir + '/bin',
            "fastq_dir": clean_fq,
            "bin_id": self.bin_name,
            "qc_stat": self.option('specimen_info').prop['path'],
        }
        self.logger.info(opts)
        self.set_run(opts, self.assemble, 'assemble', self.step.assemble)

    def run_spring_data(self):
        if self.option('qc'):
            if self.option('rm_host'):
                clean_fq = self.rm_host.option('result_fq_dir').prop["path"]
            else:
                clean_fq = self.qc.option('after_qc_dir').prop["path"]
        else:
            clean_fq = self.clean_data.output_dir + '/data'
        opts = {
            "file_dir": clean_fq,
            "type": "compress"
        }
        self.logger.info(opts)
        self.set_run(opts, self.spring, 'spring', self.step.spring)

    def run_gene_pre(self):
        taxon = self.get_taxon(self.bin_stat.output_dir + '/all.bin.summary.xls')
        opts = {
            "genome": self.assemble.option("scaffold"),
            "taxon": taxon,
            "bin_id": self.bin_name,
            "fastq1": self.assemble.option("fastq1"),
            "fastq2": self.assemble.option("fastq2"),
            "fastqs": self.assemble.option("fastqs"),
        }
        self.set_run(opts, self.predict, 'predict', self.step.predict)

    def run_anno(self):
        opts = {
            'gene_seq': self.predict.option('sample_gene_faa'),
            'gene_gff': self.predict.option('sample_gene_gff'),
            'sample': 'G_' + self.bin_name,
        }
        self.anno.set_options(opts)
        self.anno.on('start', self.set_step, {'start': self.step.anno})
        self.anno.on('end', self.set_step, {'end': self.step.anno})
        if self.predict.option('sample_rrna_fnn').is_set:
            self.anno.on("end", self.run_bin_taxon)
            self.anno.on('end', self.set_output, 'anno')
        else:
            self.anno.on("end", self.set_output, 'anno_1')
        self.anno.run()

    def run_bin_taxon(self):
        opts = {
            "blasr_db": "silva132",
            "method": "blasr",
            "query": self.predict.option('sample_rrna_fnn')
        }
        self.bin_tax.set_options(opts)
        self.bin_tax.on('start', self.set_step, {'start': self.step.bin_tax})
        self.bin_tax.on('end', self.set_step, {'end': self.step.bin_tax})
        self.bin_tax.on('end', self.set_output, 'bin_tax')
        self.bin_tax.on('end', self.end)
        self.bin_tax.run()

    def set_output(self, event):
        """
        将各个模块的结果输出至output
        """
        bin_name = 'G_' + self.bin_name
        obj = event['bind_object']
        if event['data'] == 'bin':
            link_dir(obj.output_dir, self.output_dir + '/Binning_result')
        if event['data'] == 'spring':
            link_dir(obj.output_dir, self.output_dir + '/Binning_result/clean_data')
        if event['data'] == 'ass_seq':
            if not os.path.exists(self.output_dir + '/Binning_result'):
                os.mkdir(self.output_dir + '/Binning_result')
            if os.path.exists(self.output_dir + '/Binning_result/all.uniContigs.fa.gz'):
                os.remove(self.output_dir + '/Binning_result/all.uniContigs.fa.gz')
            os.link(obj.output_dir + '/all.uniContigs.fa.gz', self.output_dir + '/Binning_result/all.uniContigs.fa.gz')
        if event['data'] == 'bin_stat':
            files = os.listdir(obj.output_dir)
            for i in files:
                if os.path.exists(self.output_dir + '/Binning_result/' + i):
                    os.remove(self.output_dir + '/Binning_result/' + i)
                os.link(obj.output_dir + '/' + i, self.output_dir + '/Binning_result/' + i)
        if event['data'] == 'assemble':
            link_dir(obj.output_dir, self.output_dir + '/' + bin_name + '/Assembly_summary')
        if event['data'] == 'predict':
            link_dir(obj.output_dir, self.output_dir + '/' + bin_name)
        if event['data'] == 'anno':
            link_dir(obj.output_dir, self.output_dir + '/' + bin_name + '/annotation')
        if event['data'] == 'anno_1':
            link_dir(obj.output_dir, self.output_dir + '/' + bin_name + '/annotation')
            self.end()
        if event['data'] == 'bin_tax':
            if len(os.listdir(obj.output_dir)) >0:
               link_dir(obj.output_dir, self.output_dir + '/' + bin_name + '/Taxon_identify/16S_identify')

    def get_binname(self, file):
        bin = []
        bin_id = ''
        with open(file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                bins = line.rstrip('\n\r\t').split('\t')
                if float(bins[7]) < 10:
                    bin.append(bins[0])
                if len(bin) ==1:
                    bin_id = bin[0]
        return bin_id

    def get_taxon(self, file):
        with open(file, 'r') as f:
            lines = f.readlines()
            bin_name = lines[1].rstrip('\n\r\t').split('\t')[-1]
        return bin_name

    def get_sof(self, file):
        with open(file, 'rb') as f2:
            lines = f2.readlines()
            line = lines[1].rstrip('\r\n\t').split('\t')
            result = line[0]
        return result

    def get_in(self, file):
        with open (file, 'r') as f:
            lines = f.readlines()
            line = lines[1].rstrip('\r\n\t').split('\t')
            insert = line[0]
            max = line[1]
        return insert, max

    def end(self):
        self.run_api()
        self.run_move()
        self.send_files()
        super(MetagbinWorkflow, self).end()

    def run_api(self):
        self.api_bin_stat()
        self.api_assemble()
        self.api_predict()
        self.api_anno()
        self.api_16s_taxon()
        self.update_sg_task()

    def run_move(self):
        bin_name = 'G_' + self.bin_name
        if os.path.exists(self.output_dir + '/' + bin_name + "/Assembly_summary"):
            shutil.rmtree(self.output_dir + '/' + bin_name + "/Assembly_summary")
        if os.path.exists(self.output_dir + '/' + bin_name + "/Assembly/Assembly_summary/Genome_coverage.xls"):
            os.remove(self.output_dir + '/' + bin_name + "/Assembly/Assembly_summary/Genome_coverage.xls")
        if os.path.exists(self.output_dir + '/' + bin_name + "/Assembly/Genome_assessment/Assess_soft.txt"):
            os.remove(self.output_dir + '/' + bin_name + "/Assembly/Genome_assessment/Assess_soft.txt")
        if os.path.exists(self.output_dir + '/' + bin_name + "/Gene_predict/length_distribute.txt"):
            os.remove(self.output_dir + '/' + bin_name + "/Gene_predict/length_distribute.txt")
        if os.path.exists(self.output_dir + '/' + bin_name + "/annotation/KEGG/kegg_pathway_img.tar.gz"):
            os.remove(self.output_dir + '/' + bin_name + "/annotation/KEGG/kegg_pathway_img.tar.gz")

    def api_bin_stat(self):
        min_contig=self.option('min_contig')
        sofware_bin = self.option('sofware_bin')
        len_cut = int(self.option("min_contig"))
        pi_path = self.api.api("metagbin.common_api")
        sofware_bin.replace("metabat","metabat2")
        bin_params = {"soft": sofware_bin, "cuto_ff": min_contig}
        bin_id = pi_path.add_main('bin', name="Bin_origin", params= bin_params, desc ="binning的结果内容!", others={"len_cut":len_cut})
        mongo_key = 'bin_id,len,scf_num,lon_len,n50,aver_scf,comp,cont,strain_heter,domain'
        pi_path.add_main_detail(self.output_dir + '/Binning_result/all.bin.summary.xls', 'bin_detail', bin_id,
                                mongo_key, has_head=True, main_name='bi_id', convert_dic={'comp':float, 'cont':float})
        mongo_key2 = 'bin_id,total_marker,uniq_marker,a,b,c,d,e,f'
        pi_path.add_main_detail(self.output_dir + '/Binning_result/all.marker.xls', 'bin_marker', bin_id,
                                mongo_key2, has_head=True, main_name = 'bi_id')
        mongo_key3 = 'bin_id,taxon'
        pi_path.add_main_detail(self.output_dir + '/Binning_result/summary.anno.xls', 'bin_taxon', bin_id,
                                mongo_key3, has_head=True, main_name='bi_id')
        mongo_key4 = 'sca_id,bin_id,start,end,gc,coverage'
        pi_path.add_main_detail(self.output_dir + '/Binning_result/all_bins.scf.xls', 'bin_scaf', bin_id,
                                mongo_key4, has_head=True, main_name='bi_id')
        mongo_key4 = 'sca_id,bin_id,start,end,strand,seq'
        with open(self.output_dir + '/Binning_result/all_bins.16s.xls') as f:
            lines =f.readlines()
            if len(lines) >1:
                pi_path.add_main_detail(self.output_dir + '/Binning_result/all_bins.16s.xls', 'bin_s16', bin_id,
                                mongo_key4, has_head = True, main_name='bi_id',main_table='bin', update_dic=
                                {'main_id': bin_id, 'seq_path': self.remote_dir + 'Binning_result/clean_data/',
                                 'ref_seq': self.remote_dir + 'Binning_result/all.uniContigs.fa.gz', 'path': self.remote_dir + 'Binning_result/bin/'})

    def api_assemble(self):
        bin_name = 'G_' + self.bin_name
        api_assembly = self.api.api('metagbin.assembly')
        soft = self.get_sof(self.output_dir + '/' + bin_name + "/Assembly_summary/Assembly_soft.txt")
        insert,max = self.get_in(self.output_dir + '/' + bin_name + "/Assembly_summary/config.txt")
        assembly_params = {'genome_id': bin_name, 'soft': soft, 'task_type': 2, 'submit_location': "assemble"}
        assembly_id = api_assembly.add_assembly(bin_name, params=assembly_params, name="Assembly_" + bin_name + "_Origin")
        file_path = self.output_dir + '/' + bin_name + "/Assembly/Assembly_summary/" + bin_name + \
                    "_assembly_summary.xls"
        bin_cover_path1 = self.output_dir + '/' + bin_name + "/Assembly_summary/Bin_coverage.xls"
        assem_cover_path2 = self.output_dir + '/' + bin_name + "/Assembly/Assembly_summary/Genome_coverage.xls"
        genome_path = self.remote_dir + bin_name + "/Assembly/Assembly_summary/" + bin_name + "_scaffold.fa"
        api_assembly.add_assembly_detail(assembly_id, bin_name, insert, max, file_path=file_path, bin_cover_path1=bin_cover_path1,
                                         assem_cover_path2=assem_cover_path2, genome_path=genome_path, stat=self.output_dir + '/' + bin_name + "/Assembly/Assembly_summary/assembly.anno.xls")
        self.logger.info('组装结果统计写入mongo成功')

        api_assembly_assess = self.api.api('metagbin.assembly_assess')
        assess_path = self.output_dir + '/' + bin_name + '/Assembly/Genome_assessment/' + bin_name + '_assess.xls'
        soft = self.get_sof(self.output_dir + '/' + bin_name + '/Assembly/Genome_assessment/Assess_soft.txt')
        assess_params = { 'genome_id': bin_name, 'soft': soft, 'task_type': 2, 'submit_location': "assemble_assess"}
        assess_id = api_assembly_assess.add_assembly_assess(bin_name, params=assess_params)
        api_assembly_assess.add_assembly_assess_detail(assess_id, soft, bin_name, file_path=assess_path)
        self.logger.info('组装评估结果写入mongo成功')

        api_gc_depth = self.api.api("metagbin.assess_gc")
        gc_params = {"genome_id": bin_name, 'soft': 'Bowtie2', 'task_type': 2, 'submit_location': "assemble_gc"}
        gc_path = self.remote_dir + bin_name + "/Assembly/GC_depth/"
        api_gc_depth.add_assembly_gc(bin_name, gc_path, params=gc_params)
        self.logger.info('组装Gc_depth结果写入mongo成功')

    def api_predict(self):
        bin_name = 'G_' + self.bin_name
        api_predict_cds = self.api.api('metagbin.predict_cds')
        cds_params = {"genome_id": bin_name, 'soft': 'Glimmer', 'task_type': 2, 'submit_location': "predict_cds"}
        pre_id = api_predict_cds.add_predict_cds(bin_name, params=cds_params)
        sample_stat = self.output_dir + '/' + bin_name + "/Gene_predict/CDS_predict/" + bin_name + "_CDS_statistics.xls"
        api_predict_cds.add_predict_cds_stat(pre_id, bin_name, sample_stat)
        predict_gff = self.output_dir + '/' + bin_name + "/Gene_predict/CDS_predict/" + bin_name + "_CDS.gff"
        fnn_path = self.remote_dir + bin_name + "/Gene_predict/CDS_predict/" + bin_name + "_CDS.fnn"
        faa_path = self.remote_dir + bin_name + "/Gene_predict/CDS_predict/" + bin_name + "_CDS.faa"
        gff_path = self.remote_dir + bin_name + "/Gene_predict/CDS_predict/" + bin_name + "_CDS.gff"
        api_predict_cds.add_predict_cds_detail(pre_id, bin_name, predict_gff = predict_gff, fnn_path =fnn_path, faa_path = faa_path, gff_path = gff_path)
        distribute = self.output_dir + '/' + bin_name + "/Gene_predict/length_distribute.txt"
        api_predict_cds.add_predict_cds_bar(pre_id, bin_name, distribute)
        self.logger.info('基因预测cds结果写入mongo成功')

        api_predict_repeat = self.api.api('metagbin.predict_repeat')
        repeat_params = {"genome_id": bin_name, 'soft': 'Tandem Repeats Finder', "parameters":"Match=2, Mismatch=7, Delta=7, PM=80, PI=10, Minscore=80, MaxPeriod=500", 'task_type': 2, 'submit_location': "predict_repeat"}
        repeat_id = api_predict_repeat.add_predict_repeat(bin_name, params=repeat_params)
        predict_gff = self.output_dir + '/' + bin_name + "/Gene_predict/repeats/" + bin_name + "_TRF.gff"
        if os.path.exists(predict_gff):
            api_predict_repeat.add_predict_repeat_detail(repeat_id, bin_name, predict_gff)
        self.logger.info('重复序列预测结果写入mongo成功')

        api_predict_rrna = self.api.api('metagbin.predict_rrna')
        rrna_params = {"genome_id": bin_name, 'soft': 'Barrnap', 'task_type': 2, 'submit_location': "predict_rrna"}
        rrna_id = api_predict_rrna.add_predict_rrna(bin_name, params=rrna_params)
        predict_gff = self.output_dir + '/' + bin_name + "/Gene_predict/rRNA/" + bin_name + "_rRNA.gff"
        rna_path = self.remote_dir + bin_name + "/Gene_predict/rRNA/" + bin_name + '-16S_rRNA.fa'
        fnn = self.output_dir + '/' + bin_name + "/Gene_predict/rRNA/" + bin_name + "_rRNA.fnn"
        if os.path.exists(predict_gff):
            api_predict_rrna.add_predict_rrna_detail(rrna_id, bin_name, predict_gff = predict_gff)
            if os.path.exists(self.output_dir + '/' + bin_name + "/Gene_predict/rRNA/" + bin_name + '-16S_rRNA.fa'):
                api_predict_rrna.add_predict_rrna_seq(rrna_id, bin_name, fnn=fnn, predict_gff=rna_path)
            else:
                api_predict_rrna.add_predict_rrna_seq(rrna_id, bin_name, fnn=fnn)
        self.logger.info('rRNA预测结果写入mongo成功')

        api_predict_trna = self.api.api('metagbin.predict_trna')
        trna_params = {"genome_id": bin_name, 'soft': 'tRNAscan-SE', 'task_type': 2, 'submit_location': "predict_trna"}
        trna_id = api_predict_trna.add_predict_trna(bin_name, params=trna_params)
        predict_gff = self.output_dir + '/' + bin_name + "/Gene_predict/tRNA/" + bin_name + "_tRNA.gff"
        if os.path.exists(predict_gff):
            api_predict_trna.add_predict_trna_detail(trna_id, bin_name, predict_gff)
        self.logger.info('tRNA预测结果写入mongo成功')

        api_genome_path = self.api.api("metagbin.common_api")
        soft = self.get_sof(self.output_dir + '/' + bin_name + '/Assembly/Genome_assessment/Assess_soft.txt')
        assess_path = self.output_dir + '/' + bin_name + '/Assembly/Genome_assessment/' + bin_name + '_assess.xls'
        complete = ''
        if soft == 'BUSCO':
            with open(assess_path, 'rb') as f:
                line = f.readlines()
                line = line[1].strip().split('\t')
                complete = line[0]
        else:
            with open(assess_path, 'rb') as f1:
                line = f1.readlines()
                line = line[1].strip().split('\t')
                complete = line[0]
        taxon_table = self.output_dir + '/' + bin_name + '/Assembly/Assembly_summary/assembly.anno.xls'
        with open(taxon_table, 'r') as f:
            lines = f.readlines()
            assess_taxon = lines[1].strip('\r\n').split('\t')[-1]
        genome_params = {"genome_id": bin_name, "complete": complete}
        api_genome_path.add_main("genome", name="Genome_origin", params=genome_params, others={"bin_id": self.bin_name, 'genome_id': bin_name, "taxon": assess_taxon})
        self.logger.info('基因组管理表导入mongo成功')

    def api_anno(self):
        bin_name = 'G_' + self.bin_name
        anno_nr = self.api.api("metagbin.anno_nr")
        nr_params = {"genome_id": bin_name, "nr": "Diamond", 'task_type': 2, 'submit_location': "anno_nr"}
        nr_id = anno_nr.add_anno_nr(params=nr_params)
        anno_nr.add_anno_nr_detail(nr_id, bin_name, self.output_dir + '/' + bin_name + "/annotation/NR/" + bin_name + "_anno_nr.xls")
        anno_cog = self.api.api("metagbin.anno_cog")
        cog_params = {"genome_id": bin_name, "cog": "Diamond", 'task_type': 2, 'submit_location': "anno_cog"}
        cog_id = anno_cog.add_anno_cog(params=cog_params)
        anno_cog.add_anno_cog_detail(cog_id, bin_name,
                                          self.output_dir + '/' + bin_name + "/annotation/COG/" + bin_name + "_cog_anno.xls")
        anno_kegg = self.api.api("metagbin.anno_kegg")
        kegg_params = {"genome_id": bin_name, "kegg": "Diamond", 'task_type': 2, 'submit_location': "anno_kegg"}
        kegg_id = anno_kegg.add_anno_kegg(self.remote_dir, "/annotation/KEGG/", "_kegg_pathway_img",
                                               params = kegg_params)
        anno_kegg.add_anno_kegg_detail(kegg_id, bin_name, self.output_dir + '/' + bin_name + "/annotation/KEGG/" + bin_name + "_kegg_anno.xls")
        anno_kegg.add_anno_kegg_level(kegg_id, bin_name, self.output_dir + '/' + bin_name + "/annotation/KEGG/" + bin_name + "_kegg_level_stat.xls")
        anno_cazy = self.api.api("metagbin.anno_cazy")
        cazy_params = {"genome_id": bin_name, "cazy": "hmmscan", 'task_type': 2, 'submit_location': "annocazy"}
        cazy_id = anno_cazy.add_anno_cazy(params=cazy_params)
        anno_cazy.add_anno_cazy_detail(cazy_id, bin_name, self.output_dir + '/' + bin_name + "/annotation/CAZY/" + bin_name + "_anno_cazy.xls")
        api_card = self.api.api("metagbin.anno_card")
        card_params = {"genome_id": bin_name, "card": "Diamond", 'task_type': 2, 'submit_location': "annocard"}
        api_card.add_anno_card(self.output_dir + '/' + bin_name + "/annotation/CARD", main=True, genome_id=bin_name, params=card_params)
        api_sum = self.api.api("metagbin.anno_summary")
        api_stat = self.api.api("metagbin.anno_stat")
        sum_params = {"genome_id": bin_name, 'task_type': 2, 'submit_location': "annosummmary"}
        stat_params = {"genome_id": bin_name, 'task_type': 2, 'submit_location': "anno_stat"}
        sum_id = api_sum.add_anno_summary(params=sum_params)
        stat_id = api_stat.add_anno_summary(params=stat_params)
        api_sum.add_anno_summary_detail(sum_id, bin_name, self.output_dir + '/' + bin_name + "/annotation/Summary/" + bin_name + "_anno_summary.xls")
        api_stat.add_anno_stat_detail(stat_id, bin_name, self.output_dir + '/' + bin_name + "/annotation/Summary/" + bin_name + "_anno_stat.xls")

    def api_16s_taxon(self):
        bin_name = 'G_' + self.bin_name
        if os.path.exists(self.output_dir + '/' + bin_name + "/Taxon_identify/16S_identify/blasr.xls"):
            pi_path = self.api.api("metagbin.common_api")
            taxonomy_params = {"query": bin_name, "db": "silva128", 'task_type': 2, 'submit_location': "identif_sixteen"}
            taxonomy_params = json.dumps(taxonomy_params, sort_keys=True, separators=(',', ':'))
            s16_id = pi_path.add_main("identif_16s", name="16S_origin", params=taxonomy_params)
            pi_path.add_main_detail(self.output_dir + '/' + bin_name + "/Taxon_identify/16S_identify/blasr.xls", "identif_16s_detail", s16_id, "genome_id,taxon,identify", has_head=True,
                                main_name="s16_id", main_table='identif_16s', update_dic={'main_id': s16_id})
            pi_path.updat_genome_taxon(self._sheet.id, bin_name, self.output_dir + '/' + bin_name + "/Taxon_identify/16S_identify/blasr.xls")

    def update_sg_task(self):
        """
        数据库更新sg_task表中，记录版本信息
        :return:
        """
        api_path = self.api.api("metagbin.common_api")
        task_id = self._sheet.id
        database_type = json.dumps({
            "nr": "v20200604",
            "cazy": "v8",
            "card": "v3.0.9",
            "kegg": "v94.2",
            "eggnog": "v4.5.1"
        },sort_keys=True, separators=(',', ':'))
        api_path.update_mongo('sg_task',{"task_id":task_id}, {"database":database_type})

    def send_files(self):
        bin_name = 'G_' + self.bin_name
        dir_o = self.output_dir
        dir_up = os.path.join(self.work_dir, 'upload_results')
        if os.path.exists(dir_up):
            shutil.rmtree(dir_up)
        os.mkdir(dir_up)
        repaths = []
        regexps = []
        if os.path.exists(os.path.join(dir_o, "Binning_result")):
            link_dir(os.path.join(dir_o, "Binning_result"), os.path.join(dir_up, "Binning_result"))
        if os.path.exists(os.path.join(dir_o, bin_name)):
            link_dir(os.path.join(dir_o, bin_name), os.path.join(dir_up, bin_name))
        repaths += [
            [".", "", "工作流分析结果目录", 0],
            ["Binning_result", "", "Binning的相关结果目录", 0],
            ["Binning_result/all_bins.scf.xls  ", "xls", "所有scaffold的统计信息", 0],
            ["Binning_result/all.marker.xls", "xls", "所有bin的marker基因统计信息", 0],
            ["Binning_result/all.bin.summary.xls", "xls", "每个bin的统计信息", 0],
            ["Binning_result/all_bins.16s.xls  ", "xls", "存在16sRNA的bin的统计信息", 0],
            ["Binning_result/summary.anno.xls", "xls", "每个bin的物种注释信息", 0],
            ["Binning_result/bins.rRNA.fnn  ", "xls", "预测16s的序列文件", 0],
            ["Binning_result/bins.rRNA.gff  ", "xls", "预测出的rRNA的gff文件", 0],
            ["Binning_result/summary.anno.xls", "xls", "每个bin的物种注释信息", 0],
            ["Binning_result/bin", "", "bin的序列文件目录", 0],
            ["%s" %bin_name, "", "%s的结果目录" %bin_name, 0],
            ["%s/Assembly" % bin_name, "", "%s的基因组组装结果目录" % bin_name, 0],
            ["%s/Assembly/Assembly_summary" % bin_name, "", "%s的基因组组装统计结果目录" % bin_name, 0],
            ["%s/Assembly/Genome_assessment" % bin_name, "", "%s的基因组评估结果目录" % bin_name, 0],
            ["%s/Assembly/GC_depth" % bin_name, "", "%s的基因组GC_depth" % bin_name, 0],
            ["%s/Gene_predict" % bin_name, "", "%s的基因组的基因预测目录" % bin_name, 0],
            ["%s/Gene_predict/CDS_predict" % bin_name, "", "%s基因组编码基因预测目录" % bin_name, 0],
            ["%s/Gene_predict/rRNA" % bin_name, "", "%s基因组rRNA预测结果目录" % bin_name, 0],
            ["%s/Gene_predict/tRNA" % bin_name, "", "%s基因组tRNA预测结果目录" % bin_name, 0],
            ["%s/Gene_predict/repeats" % bin_name, "", "%s基因组重复序列预测结果目录" % bin_name, 0],
            ["%s/annotation" % bin_name, "", "%s基因组的功能注释结果目录" % bin_name, 0],
            ["%s/annotation/NR" % bin_name, "", "%s基因组的NR注释结果目录" % bin_name, 0],
            ["%s/annotation/COG" % bin_name, "", "%s基因组的COG注释结果目录" % bin_name, 0],
            ["%s/annotation/KEGG" % bin_name, "", "%s基因组的KEGG注释结果目录" % bin_name, 0],
            ["%s/annotation/CAZY" % bin_name, "", "%s基因组的CAZY注释结果目录" % bin_name, 0],
            ["%s/annotation/CARD" % bin_name, "", "%s基因组的CARD注释结果目录" % bin_name, 0],
            ["%s/annotation/Summary" % bin_name, "", "%s基因组的注释结果汇总目录" % bin_name, 0],
            ["%s/Taxon_identify" % bin_name, "", "%s基因组物种鉴定模块目录" % bin_name, 0],
            ["%s/Taxon_identify/16S_identify" % bin_name, "", "%s基因组的16S物种分类结果目录" % bin_name, 0],
            ["%s/Taxon_identify/16S_identify/blasr_result.xls" % bin_name, "xls", "%s基因组的16S物种分类结果表" % bin_name, 0],
        ]
        regexps += [
            [r"Binning_result/bin/bin[0-9]+\.fa", "", "bin的序列文件", 0],
            [r"%s/Assembly/Assembly_summary/.+_assembly_scaffold_details.xls" % bin_name, "", "%s的基因组组装scaffold统计结果" % bin_name, 0],
            [r"%s/Assembly/Assembly_summary/.+_assembly_contig_details.xls" % bin_name, "", "%s的基因组组装contig统计结果" % bin_name, 0],
            [r"%s/Assembly/Assembly_summary/.+_assembly_summary.xls" % bin_name, "", "%s的基因组组装统计信息" % bin_name, 0],
            [r"%s/Assembly/Assembly_summary/.+_scaffold.fa" % bin_name, "", "%s的基因组组装序列结果文件" % bin_name, 0],
            [r"%s/Assembly/Genome_assessment/.+_assess.xls" % bin_name, "xls", "%s的基因组评估结果文件" % bin_name, 0],
            [r"%s/Assembly/GC_depth/gc_depth.png" % bin_name, "", "%s的基因组GC_depth_1000统计图(PNG)" % bin_name, 0],
            [r"%s/Assembly/GC_depth/gc_depth.svg" % bin_name, "", "%s的基因组GC_depth_1000统计图(SVG)" % bin_name, 0],
            [r"%s/Gene_predict/CDS_predict/.+_CDS.faa" % bin_name, "", "%s基因组编码基因预测蛋白序列" % bin_name, 0],
            [r"%s/Gene_predict/CDS_predict/.+_CDS.fnn" % bin_name, "", "%s基因组编码基因预测核苷酸序列" % bin_name, 0],
            [r"%s/Gene_predict/CDS_predict/.+_CDS.gff" % bin_name, "", "%s基因组编码基因预测统计表(gff格式)" % bin_name, 0],
            [r"%s/Gene_predict/CDS_predict/.+_CDS_statistics.xls" % bin_name, "", "%s基因组编码基因预测统计表" % bin_name, 0],
            [r"%s/Gene_predict/rRNA/.+_rRNA.fnn" % bin_name, "", "%s基因组rRNA预测核苷酸序列文件" % bin_name, 0],
            [r"%s/Gene_predict/rRNA/.+_rRNA.gff" % bin_name, "", "%s基因组rRNA预测统计表(gff格式)" % bin_name, 0],
            [r"%s/Gene_predict/rRNA/.+-16s_rRNA.fnn" % bin_name, "", "%s基因组rRNA预测16S核苷酸序列文件" % bin_name, 0],
            [r"%s/Gene_predict/tRNA/.+_tRNA.gff" % bin_name, "", "%s基因组tRNA预测统计表(gff格式)" % bin_name, 0],
            [r"%s/Gene_predict/tRNA/.+_tRNA.fnn" % bin_name, "", "%s基因组tRNA预测核苷酸序列文件" % bin_name, 0],
            [r"%s/Gene_predict/repeats/.+_TRF.gff" % bin_name, "", "%s基因组重复序列预测统计表(gff格式)" % bin_name, 0],
            [r"%s/annotation/NR/.+_anno_nr.xls" % bin_name, "", "%s基因组的NR注释结果表" % bin_name, 0],
            [r"%s/annotation/COG/.+_cog_summary.xls" % bin_name, "", "%s基因组的COG注释结果分类统计表" % bin_name, 0],
            [r"%s/annotation/COG/.+_cog_anno.xls" % bin_name, "", "%s基因组的COG注释结果表" % bin_name, 0],
            [r"%s/annotation/KEGG/.+_kegg_pathway_img" % bin_name, "", "%s基因组的KEGG通路图文件目录" % bin_name, 0],
            [r"%s/annotation/KEGG/.+_kegg_anno.xls  " % bin_name, "", "%s基因组的KEGG注释结果文件" % bin_name, 0],
            [r"%s/annotation/KEGG/.+_kegg_level_stat.xls" % bin_name, "", "%s基因组的KEGG注释level统计文件" % bin_name, 0],
            [r"%s/annotation/CAZY/.+_anno_cazy.xls" % bin_name, "", "%s基因组的CAZY注释分类结果表" % bin_name, 0],
            [r"%s/annotation/CAZY/.+_cazy_class_stat.xls" % bin_name, "", "%s基因组的CAZY注释class分类结果表" % bin_name, 0],
            [r"%s/annotation/CAZY/.+_cazy_family_stat.xls" % bin_name, "", "%s基因组的CAZY注释family分类结果表" % bin_name, 0],
            [r"%s/annotation/CAZY/.+_cazy_parse_anno.xls" % bin_name, "", "%s基因组的CAZY注释结果表" % bin_name, 0],
            [r"%s/annotation/CARD/.+_card_align.xls" % bin_name, "", "%s基因组的CARD比对结果表" % bin_name, 0],
            [r"%s/annotation/CARD/.+_card_anno.xls" % bin_name, "", "%s基因组的CARD注释结果表" % bin_name, 0],
            [r"%s/annotation/CARD/.+_card_category.xls" % bin_name, "", "%s基因组的CARD注释分类结果表" % bin_name, 0],
            [r"%s/annotation/Summary/.+_anno_summary.xls" % bin_name, "", "%s基因组的注释结果汇总表" % bin_name, 0],
            [r"%s/annotation/Summary/.+_anno_stat" % bin_name, "", "%s基因组的基因注释统计表" % bin_name, 0],
        ]
        sdir = self.add_upload_dir(dir_up)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)