# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.denovo_rna.assemble.trinity_stat import *
import os
import re


class AssembleAgent(Agent):
    """
    无参转录组拼接
    version v1.0.1
    author: qiuping
    last_modify: 2016.08.31
    """
    def __init__(self, parent):
        super(AssembleAgent, self).__init__(parent)
        options = [
            {"name": "fq_type", "type": "string"},  # PE OR SE
            {"name": "fq_l", "type": "infile", "format": "sequence.fastq"},  # PE测序，所有样本fastq左端序列文件
            {"name": "fq_r", "type": "infile", "format": "sequence.fastq"},  # PE测序，所有样本fastq右端序列文件
            {"name": "fq_s", "type": "infile", "format": "sequence.fastq"},  # SE测序，所有样本fastq序列文件
            {"name": "length", "type": "string", "default": "100,200,400,600,800,1000,1200,1500,2000"},  # 统计步长
            {"name": "cpu", "type": "int", "default": 6},  # trinity软件所分配的cpu数量
            {"name": "max_memory", "type": "string", "default": '100G'},  # trinity软件所分配的最大内存，单位为GB
            {"name": "min_contig_length", "type": "int", "default": 200},  # trinity报告出的最短的contig长度。默认为200
            {"name": "kmer_size", "type": "int", "default": 25},
            {"name": "min_kmer_cov", "type": "int", "default": 2},
            {"name": "SS_lib_type", "type": "string", "default": 'none'},  # reads的方向，成对的reads: RF or FR; 不成对的reads: F or R，默认情况下，不设置此参数
            {"name": "gene_fa", "type": "outfile", "format": "sequence.fasta"},
            {"name": "trinity_fa", "type": "outfile", "format": "sequence.fasta"},
            {"name": "gene_full_name", "type": "outfile", "format": "denovo_rna.express.gene_list"}
        ]
        self.add_option(options)
        self.step.add_steps("assemble")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.assemble.start()
        self.step.update()

    def stepfinish(self):
        self.step.assemble.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('fq_type'):
            raise OptionError('必须设置测序类型：PE OR SE')
        if self.option('fq_type') not in ['PE', 'SE']:
            raise OptionError('测序类型不在所给范围内')
        if not self.option("fq_l").is_set and not self.option("fq_r").is_set and not self.option("fq_s").is_set:
            raise OptionError("必须设置PE测序输入文件或者SE测序输入文件")
        if self.option("fq_type") == "PE" and not self.option("fq_r").is_set and not self.option("fq_l").is_set:
            raise OptionError("PE测序时需设置左端序列和右端序列输入文件")
        if self.option("fq_type") == "SE" and not self.option("fq_s").is_set:
            raise OptionError("SE测序时需设置序列输入文件")
        if self.option("SS_lib_type") != 'none' and self.option("SS_lib_type") not in ['F', 'R', 'FR', 'RF']:
            raise OptionError("所设reads方向不在范围值内")
        if self.option('kmer_size') > 32 or self.option('kmer_size') < 1:
            raise OptionError("所设kmer_size不在范围内，请检查")
        if self.option('min_kmer_cov') < 1:
            raise OptionError("所设min_kmer_cov不在范围内，请检查")
        # if self.option("fq_type") == "SE":
        #     self.option('fq_s').check_content()
        # if self.option("fq_type") == "PE":
        #     self.option('fq_r').check_content()
        #     self.option('fq_l').check_content()
        self.set_resource()
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        file_size = 0
        if self.option('fq_type') == 'SE':
            file_size = float(os.path.getsize(self.option('fq_s').prop['path'])) / 1024 / 1024
        else:
            file_size = float(os.path.getsize(self.option('fq_r').prop['path'])) / 1024 / 1024 + float(os.path.getsize(self.option('fq_l').prop['path'])) / 1024 / 1024
        if file_size <= 1024 * 0.5:
            self._cpu = 8
            self._memory = '53G'
        elif file_size <= 1024 * 4 and file_size > 1024 * 0.5:
            self._cpu = 12
            self._memory = '83G'
        elif file_size <= 1024 * 10 and file_size > 1024 * 4:
            self._cpu = 14
            self._memory = '103G'
        elif file_size < 1024 * 20 and file_size > 1024 * 10:
            self._cpu = 18
            self._memory = '128G'
        else:
            self._cpu = 26
            self._memory = '200G'
        self.option('cpu', self._cpu - 2)
        mem = str(int(self._memory.strip('G')) - 3) + 'G'
        self.option('max_memory', mem)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["transcript.iso.txt", "txt", "按照拥有的isoform（可变剪接体）数目，统计转录本的数量分布的文件"],
            ["trinity.fasta.stat.xls", "xls", "trinity.fasta文件统计信息，信息包括（转录本、基因的）序列总数，碱基总数，"
                                              "GC含量，最长（短）转录本长度，平均长度，N50，N90"]
        ])
        result_dir.add_regexp_rules([
            [r"length.distribut.txt$", "txt", "长度分布信息统计文件"]
        ])
        super(AssembleAgent, self).end()


class AssembleTool(Tool):
    def __init__(self, config):
        super(AssembleTool, self).__init__(config)
        self._version = "v1.0.1"
        self.trinity_path = '/bioinfo/rna/trinityrnaseq-2.2.0/'
        self.bowtie = self.config.SOFTWARE_DIR + '/bioinfo/align/bowtie-1.1.2/'
        self.samtools = self.config.SOFTWARE_DIR + '/bioinfo/align/samtools-1.3.1/'
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.java = self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin'
        self.set_environ(PATH=self.bowtie)
        self.set_environ(PATH=self.samtools)
        self.set_environ(PATH=self.java)
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)

    def run(self):
        """
        运行
        :return:
        """
        super(AssembleTool, self).run()
        self.run_trinity()

    def run_trinity(self):
        """
        运行trinity软件，进行拼接组装
        """
        if self.option('fq_type') == 'SE':
            if self.option('SS_lib_type') == 'none':
                cmd = self.trinity_path + 'Trinity --seqType fq --max_memory %s --min_contig_length %s --CPU %s --single %s --no_version_check --KMER_SIZE %s --min_kmer_cov %s --bflyCalculateCPU' % (self.option('max_memory'), self.option('min_contig_length'), self.option('cpu'), self.option('fq_s').prop['path'], self.option('kmer_size'), self.option('min_kmer_cov'))
            else:
                cmd = self.trinity_path + 'Trinity --seqType fq --max_memory %s --min_contig_length %s --CPU %s --single %s --SS_lib_type %s --no_version_check --KMER_SIZE %s --min_kmer_cov %s --bflyCalculateCPU' % (self.option('max_memory'), self.option('min_contig_length'), self.option('cpu'), self.option('fq_s').prop['path'], self.option('SS_lib_type'), self.option('kmer_size'), self.option('min_kmer_cov'))
            if float(self.option('fq_s').get_size()) / 1024 / 1024 / 1024 > 15:
                cmd += ' --normalize_reads'
        else:
            if self.option('SS_lib_type') == 'none':
                cmd = self.trinity_path + 'Trinity --seqType fq --max_memory %s --min_contig_length %s --CPU %s --left %s --right %s --no_version_check --KMER_SIZE %s --min_kmer_cov %s --bflyCalculateCPU' % (self.option('max_memory'), self.option('min_contig_length'), self.option('cpu'), self.option('fq_l').prop['path'], self.option('fq_r').prop['path'], self.option('kmer_size'), self.option('min_kmer_cov'))
            else:
                cmd = self.trinity_path + 'Trinity --seqType fq --max_memory %s --min_contig_length %s --CPU %s --left %s --right %s --SS_lib_type %s --no_version_check --KMER_SIZE %s --min_kmer_cov %s --bflyCalculateCPU' % (self.option('max_memory'), self.option('min_contig_length'), self.option('cpu'), self.option('fq_l').prop['path'], self.option('fq_r').prop['path'], self.option('SS_lib_type'), self.option('kmer_size'), self.option('min_kmer_cov'))
            if float(self.option('fq_l').get_size()) / 1024 / 1024 / 1024 + float(self.option('fq_r').get_size()) / 1024 / 1024 / 1024 > 15:
                cmd += ' --normalize_reads'
        self.logger.info('运行trinity软件，进行组装拼接')
        command = self.add_command("trinity_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("trinity运行完成")
            trinity_info = get_trinity_info(self.work_dir + '/trinity_out_dir/Trinity.fasta')
            result = stat_info(trinity_info, self.work_dir + '/gene.fasta', self.work_dir + '/trinity.fasta.stat.xls', self.work_dir, self.option('length'), self.work_dir + '/gene_full_name.txt')  # 运行trnity_stat.py，对trinity.fatsa文件进行统计
            if result:
                self.logger.info('trinity_stat.py运行成功，统计trinity.fasta信息成功')
                self.set_output()
            else:
                self.logger.info('trinity_stat.py运行出错，统计trinity.fasta信息失败')
        else:
            self.set_error("trinity运行出错!")

    def trinity_cmd_check(self, command, line):
        line = line.strip("\n")
        if re.match(r"succeeded.+?completed\.$", line):
            self.logger.debug(line)
        else:
            pass

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        files = os.listdir(self.work_dir)
        try:
            for f in files:
                if re.search(r'length.distribut.txt$', f):
                    os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
            os.link(self.work_dir + '/trinity.fasta.stat.xls', self.output_dir + '/trinity.fasta.stat.xls')
            self.option('gene_fa', self.work_dir + '/gene.fasta')
            self.option('trinity_fa', self.work_dir + '/trinity_out_dir/Trinity.fasta')
            self.option('gene_full_name', self.work_dir + '/gene_full_name.txt')
            self.logger.info("设置组装拼接分析结果目录成功")
            self.end()
        except Exception as e:
            self.logger.info("设置组装拼接分析结果目录失败{}".format(e))
            self.set_error("设置组装拼接分析结果目录失败{}".format(e))
