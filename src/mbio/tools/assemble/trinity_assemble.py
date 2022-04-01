# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.denovo_rna.assemble.trinity_stat import *
import os
import shutil


class TrinityAssembleAgent(Agent):
    """
    宏转录组 Trinity拼接
    """
    def __init__(self, parent):
        super(TrinityAssembleAgent, self).__init__(parent)
        options = [
            {"name": "fq_type", "type": "string"},  # PE OR SE
            {"name": "fq_l", "type": "infile", "format": "sequence.fastq"},  # PE测序，所有样本fastq左端序列文件
            {"name": "fq_r", "type": "infile", "format": "sequence.fastq"},  # PE测序，所有样本fastq右端序列文件
            {"name": "fq_s", "type": "infile", "format": "sequence.fastq"},  # SE测序，所有样本fastq序列文件
            {"name": "length", "type": "string", "default": "100,200,400,600,800,1000,1200,1500,2000"},  # 统计步长
            {"name": "cpu", "type": "int", "default": 6},  # trinity软件所分配的cpu数量
            {"name": "max_memory", "type": "string", "default": '100G'},  # trinity软件所分配的最大内存，单位为GB
            {"name": "min_contig_length", "type": "int", "default": 200},  # trinity报告出的最短的contig长度。默认为200
            {"name": "kmer_size", "type": "int", "default": 25}, #kmer大小，默认为25
            {"name": "min_kmer_cov", "type": "int", "default": 2}, # 最小kmer的丰度，默认为2
            {"name": "sample", "type": "string"}, ## 样本名称
        ]
        self.add_option(options)
        self.step.add_steps("assemble")

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        if not self.option('fq_type'):
            raise OptionError('必须设置测序类型：PE OR SE', code="31302101")
        if self.option('fq_type') not in ['PE', 'SE']:
            raise OptionError('测序类型不在所给范围内', code="31302102")
        if (not self.option("fq_l").is_set) and (not self.option("fq_r").is_set) and (not self.option("fq_s").is_set):
            raise OptionError("必须设置PE测序输入文件或者SE测序输入文件", code="31302103")
        if (self.option("fq_type") == "PE") and (not self.option("fq_r").is_set) and (not self.option("fq_l").is_set):
            raise OptionError("PE测序时需设置左端序列和右端序列输入文件", code="31302104")
        if (self.option("fq_type") == "SE") and (not self.option("fq_s").is_set):
            raise OptionError("SE测序时需设置序列输入文件", code="31302105")
        if self.option('kmer_size') > 32 or self.option('kmer_size') < 1:
            raise OptionError("所设kmer_size不在范围内，请检查", code="31302106")
        if self.option('min_kmer_cov') < 1:
            raise OptionError("所设min_kmer_cov不在范围内，请检查", code="31302107")
        if not self.option("sample"):
            raise OptionError("没有设置样本名称，请检查", code="31302108")
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
        if file_size <= 1024 * 5:
            self._cpu = 6
            self._memory = '50G'
        elif file_size <= 1024 * 10 and file_size > 1024 * 5:
            self._cpu = 8
            self._memory = '80G'
        elif file_size <= 1024 * 15 and file_size > 1024 * 10:
            self._cpu = 10
            self._memory = '100G'
        elif file_size < 1024 * 20 and file_size > 1024 * 15:
            self._cpu = 12
            self._memory = '128G'
        else:
            self._cpu = 16
            self._memory = '250G'
        self.option('cpu', self._cpu)
        mem = str(int(self._memory.strip('G')) - 3) + 'G'
        self.option('max_memory', mem)

    def end(self):
        """
        结束
        """
        self.logger.info("开始上传结果文件")
        super(TrinityAssembleAgent, self).end()


class TrinityAssembleTool(Tool):
    """
    运行tool
    """
    def __init__(self, config):
        super(TrinityAssembleTool, self).__init__(config)
        self._version = "v1.0.1"
        self.trinity_path = '/bioinfo/rna/trinityrnaseq-2.2.0/'
        self.bowtie = self.config.SOFTWARE_DIR + '/bioinfo/align/bowtie-1.1.2/'
        self.samtools = self.config.SOFTWARE_DIR + '/bioinfo/align/samtools-1.3.1/'
        self.perl_path = '/program/perl-5.24.0/bin/perl'
        self.transter_fasta = self.config.PACKAGE_DIR + '/assemble/scripts/get_id_from_fasta.pl'
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
        super(TrinityAssembleTool, self).run()
        self.run_trinity()
        self.fix_fasta()
        self.run_stat()
        self.set_output()
        self.end()

    def run_trinity(self):
        """
        运行trinity软件，进行拼接组装
        """
        if self.option('fq_type') == 'SE':
            cmd = self.trinity_path + 'Trinity --seqType fq --max_memory {0} --min_contig_length {1} --CPU {2} --single {3} --no_version_check --KMER_SIZE {4} --min_kmer_cov {5} --bflyCalculateCPU --inchworm_cpu {2} --bflyCPU {2}'.format(self.option('max_memory'), self.option('min_contig_length'), self.option('cpu'), self.option('fq_s').prop['path'], self.option('kmer_size'), self.option('min_kmer_cov'))
            if float(self.option('fq_s').get_size()) / 1024 / 1024 / 1024 > 15:
                cmd += ' --normalize_reads'
        else: ## PE数据
            cmd = self.trinity_path + 'Trinity --seqType fq --max_memory {0} --min_contig_length {1} --CPU {2} --left {3} --right {4} --no_version_check --KMER_SIZE {5} --min_kmer_cov {6} --bflyCalculateCPU --inchworm_cpu {2} --bflyCPU {2}'.format(self.option('max_memory'), self.option('min_contig_length'), self.option('cpu'), self.option('fq_l').prop['path'], self.option('fq_r').prop['path'], self.option('kmer_size'), self.option('min_kmer_cov'))
            if float(self.option('fq_l').get_size()) / 1024 / 1024 / 1024 + float(self.option('fq_r').get_size()) / 1024 / 1024 / 1024 > 15:
                cmd += ' --normalize_reads'
        self.logger.info('运行trinity软件，进行组装拼接')
        self.logger.info(cmd)
        command = self.add_command("trinity_cmd", cmd).run()
        self.wait(command)
        if command.return_code in [0]:
            self.logger.info("trinity运行完成")
        else:
            self.set_error("trinity运行出错!", code="31302101")

    def fix_fasta(self):
        """
        对Trinity软件结果进行处理
        :return:
        """
        self.logger.info("开始对fasta序列进行处理")
        fasta_path = os.path.join(self.work_dir, "Trinity.fasta")
        input_fasta = os.path.join(self.work_dir, "trinity_out_dir/Trinity.fasta")
        cmd = '{} {} {} {}'.format(self.perl_path, self.transter_fasta, input_fasta, fasta_path)
        self.logger.info(cmd)
        self.logger.info("开始对fasta序列进行处理文件")
        command = self.add_command("fix_cmd", cmd).run()
        self.wait(command)
        if command.return_code in [0]:
            self.logger.info("fix序列完成！")
        else:
            self.set_error("fix序列失败!", code="31302102")

    def run_stat(self):
        """
        对拼接结果进行统计
        :return:
        """
        try:
            self.logger.info("开始对组装拼接结果进行统计")
            trinity_info = get_trinity_info(self.work_dir + '/Trinity.fasta')
            self.logger.info("统计拼接结果完成")
        except Exception, e:
            self.set_error("统计失败%s", variables=(e), code="31302103")
        try:
            self.logger.info("开始对拼接结果进行步长统计")
            result = stat_info(trinity_info, self.work_dir + '/gene.fasta', self.work_dir + '/trinity.fasta.stat.xls', self.work_dir, self.option('length'), self.work_dir + '/gene_full_name.txt')  # 运行trnity_stat.py，对trinity.fatsa文件进行统计
            self.logger.info("统计步长成功：%s"%result)
        except Exception, e:
            self.set_error('trinity_stat.py运行出错，统计trinity.fasta信息失败: %s', variables=(e), code="31302104")

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        self.logger.info("开始设置结果文件目录")
        if os.path.exists(os.path.join(self.output_dir, self.option("sample"))):
            shutil.rmtree(os.path.join(self.output_dir, self.option("sample")))
            os.mkdir(os.path.join(self.output_dir, self.option("sample")))
        else:
            os.mkdir(os.path.join(self.output_dir, self.option("sample")))

        trinity_fasta = os.path.join(self.output_dir, self.option("sample"), "Trinity.fasta")
        if os.path.exists(trinity_fasta):
            os.remove(trinity_fasta)
        os.link(os.path.join(self.work_dir, "Trinity.fasta"), trinity_fasta)
        trinity_timing = os.path.join(self.output_dir, self.option("sample"), "Trinity.timing")
        if os.path.exists(trinity_timing):
            os.remove(trinity_timing)
        os.link(os.path.join(self.work_dir, "trinity_out_dir/Trinity.timing"), trinity_timing)
        if os.path.exists(os.path.join(self.output_dir, self.option("sample"), "Trinity.fasta.stat.xls")):
            os.remove(os.path.join(self.output_dir, self.option("sample"), "Trinity.fasta.stat.xls"))
        os.link(os.path.join(self.work_dir, 'trinity.fasta.stat.xls'), os.path.join(self.output_dir, self.option("sample"), 'Trinity.fasta.stat.xls'))
        self.logger.info("设置结果文件目录完成")

