# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.denovo_rna.assemble.trinity_stat import *
import os
import re


class TrinityAgent(Agent):
    """
    Trinity拼接
    author: 刘彬旭
    last_modify: 2017.10.23
    """
    def __init__(self, parent):
        super(TrinityAgent, self).__init__(parent)
        options = [
            {"name": "fq_type", "type": "string", "default": "PE"},  # PE OR SE
            {"name": "fq_l", "type": "infile", "format": "sequence.fastq"},  # PE测序，所有样本fastq左端序列文件
            {"name": "fq_r", "type": "infile", "format": "sequence.fastq"},  # PE测序，所有样本fastq右端序列文件
            {"name": "fq_s", "type": "infile", "format": "sequence.fastq"},  # SE测序，所有样本fastq序列文件
            {"name": "cpu", "type": "int", "default": 6},  # trinity软件所分配的cpu数量
            {"name": "max_memory", "type": "string", "default": '100G'},  # trinity软件所分配的最大内存，单位为GB
            {"name": "min_contig_length", "type": "int", "default": 200},  # trinity报告出的最短的contig长度。默认为200
            {"name": "kmer_size", "type": "int", "default": 25},  # kmer 长度 <=32
            {"name": "min_kmer_cov", "type": "int", "default": 2},
            {"name": "jaccard_clip", "type": "bool", "default": False}, #分割高密度基因区间基因
            {"name": "no_normalize_reads", "type": "bool", "default": False}, #不做reads均一化
            {"name": "normalize_max_read_cov", "type": "int", "default": 50}, #reads均一化覆盖倍数
            {"name": "SS_lib_type", "type": "string", "default": 'none'},  # reads的方向，成对的reads: RF or FR; 不成对的reads: F or R，默认情况下，不设置此2参数
            {"name": "gene_fa", "type": "outfile", "format": "denovo_rna_v2.trinity_fasta"},
            {"name": "trinity_fa", "type": "outfile", "format": "denovo_rna_v2.trinity_fasta"},
            {"name": "edge_thr", "type": "float", "default": 0.1}, #butterfly 参数
            {"name": "flow_thr", "type": "float", "default": 0.04}, #butterfly 参数
            {"name": "max_number_of_paths_per_node_init", "type": "int", "default": 20}, #butterfly 参数
            {"name": "min_per_id_same_path", "type": "int", "default": 98}, #butterfly 参数
            {"name": "gene2trans", "type": "outfile", "format": "denovo_rna_v2.common"}, #组装后gene和转录本的对应列表
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
            raise OptionError('必须设置测序类型：PE OR SE', code="32008001")
        if self.option('fq_type') not in ['PE', 'SE']:
            raise OptionError('测序类型不在所给范围内', code="32008002")
        if not self.option("fq_l").is_set and not self.option("fq_r").is_set and not self.option("fq_s").is_set:
            raise OptionError("必须设置PE测序输入文件或者SE测序输入文件", code="32008003")
        if self.option("fq_type") == "PE" and not self.option("fq_r").is_set and not self.option("fq_l").is_set:
            raise OptionError("PE测序时需设置左端序列和右端序列输入文件", code="32008004")
        if self.option("fq_type") == "SE" and not self.option("fq_s").is_set:
            raise OptionError("SE测序时需设置序列输入文件", code="32008005")
        if self.option("SS_lib_type") != 'none' and self.option("SS_lib_type") not in ['F', 'R', 'FR', 'RF']:
            raise OptionError("所设reads方向不在范围值内", code="32008006")
        if self.option('kmer_size') > 32 or self.option('kmer_size') < 1:
            raise OptionError("所设kmer_size不在范围内，请检查", code="32008007")
        if self.option('min_kmer_cov') < 1:
            raise OptionError("所设min_kmer_cov不在范围内，请检查", code="32008008")

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
            file_size = float(os.path.getsize(self.option('fq_r').prop['path'])) / 1024 / 1024 + \
            float(os.path.getsize(self.option('fq_l').prop['path'])) / 1024 / 1024
        if file_size <= 1024 * 0.5:
            self._cpu = 6
            self._memory = '80G'
        elif file_size <= 1024 * 4 and file_size > 1024 * 0.5:
            self._cpu = 6
            self._memory = '80G'
        else:
            self._cpu = 8
            self._memory = '120G'
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
        super(TrinityAgent, self).end()


class TrinityTool(Tool):
    def __init__(self, config):
        super(TrinityTool, self).__init__(config)
        self._version = "v1.0.1"
        self.trinity_path = '/bioinfo/denovo_rna_v2/trinityrnaseq-Trinity-v2.5.0/'
        self.perl = self.config.SOFTWARE_DIR + '/miniconda2/bin'
        self.bowtie = self.config.SOFTWARE_DIR + '/bioinfo/denovo_rna_v2/bowtie2-2.3.3.1-linux-x86_64/'
        self.java = self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin'
        self.set_environ(PATH=self.perl)
        self.set_environ(PATH=self.bowtie)
        self.set_environ(PATH=self.java)


    def run(self):
        """
        运行
        :return:
        """
        super(TrinityTool, self).run()
        self.run_trinity()
        self.set_output()



    def run_trinity(self):
        """
        运行trinity软件，进行拼接组装
        """
        #trinity command
        if self.option('fq_type') == 'SE':
            cmd = '{}Trinity --seqType fq --single {}'.format(
                self.trinity_path, self.option('fq_s').prop['path'])
        else:
            cmd = '{}Trinity --seqType fq --left {} --right {}'.format(
                self.trinity_path, self.option('fq_l').prop['path'], self.option('fq_r').prop['path'])
        if self.option('SS_lib_type') != 'none':
            cmd += ' --SS_lib_type {}'.format(self.option('SS_lib_type'))

        bfly_CPU = self.option('cpu')
        #bflyCalculateCPU = min( int(0.8 * int(self.option('max_memory')[:-1])), bfly_CPU )
        #trinity resource set
        cmd += ' --max_memory {} --CPU {} --bflyCPU {} --bflyCalculateCPU --bflyHeapSpaceMax {} --grid_node_CPU {} --grid_node_max_memory {}'.format(
            self.option('max_memory'), self.option('cpu'), bfly_CPU, '2G', self.option('cpu')-2 , '100G')

        cmd += ' --min_contig_length {} --KMER_SIZE {} --min_kmer_cov {} --no_version_check'.format(
            self.option('min_contig_length'), self.option('kmer_size'), self.option('min_kmer_cov'))

        if self.option('no_normalize_reads') == True:
            cmd += ' --no_normalize_reads'
        else:
            cmd += ' --normalize_max_read_cov {}'.format(self.option('normalize_max_read_cov'))
        if self.option('jaccard_clip') == True:
            cmd += ' --jaccard_clip'

        cmd += ' --bfly_opts " --edge-thr={} --flow-thr={} --max_number_of_paths_per_node_init={} --min_per_id_same_path={} "'.format(
            self.option('edge_thr'), self.option('flow_thr'), self.option('max_number_of_paths_per_node_init'), self.option('min_per_id_same_path'))

        self.logger.info('运行trinity软件，进行组装拼接')
        command = self.add_command("trinity_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("trinity运行完成")
            trinity_info = get_trinity_info(self.work_dir + '/trinity_out_dir/Trinity.fasta')
        else:
            self.set_error("trinity运行出错!", code="32008001")

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
        trinity_fa = self.work_dir + '/trinity_out_dir/Trinity.fasta'
        gene2trans = self.work_dir + '/trinity_out_dir/Trinity.fasta.gene_trans_map'
        if os.path.exists(self.output_dir + '/Trinity.fasta'):
            os.remove(self.output_dir + '/Trinity.fasta')
        if os.path.exists(self.output_dir + '/Trinity.fasta.gene_trans_map'):
            os.remove(self.output_dir + '/Trinity.fasta.gene_trans_map')
        os.link(trinity_fa, os.path.join(self.output_dir, 'Trinity.fasta'))
        os.link(gene2trans, os.path.join(self.output_dir, 'Trinity.fasta.gene_trans_map'))

        self.logger.info("设置结果目录")

        try:
            self.option('trinity_fa', os.path.join(self.output_dir, 'Trinity.fasta'))
            self.option('gene2trans', os.path.join(self.output_dir, 'Trinity.fasta.gene_trans_map'))
            self.logger.info("设置组装拼接分析结果目录成功")
            self.logger.info("序列文件为 {}".format(self.option('trinity_fa').prop['path']))
            self.logger.info("设置组装拼接分析结果目录成功 {}".format(self.option('gene2trans')))
            self.end()
        except Exception as e:
            self.logger.info("设置组装拼接分析结果目录失败{}".format(e))
            self.set_error("设置组装拼接分析结果目录失败%s", variables=(e), code="32008002")
