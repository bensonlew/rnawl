# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.denovo_rna.assemble.trinity_stat import *
import os
import re


class Trinity2Agent(Agent):
    """
    Trinity拼接
    author: 刘彬旭
    last_modify: 2017.10.23
    """
    def __init__(self, parent):
        super(Trinity2Agent, self).__init__(parent)
        options = [
            {"name": "fq_type", "type": "string", "default": "PE"},  # PE OR SE
            {"name": "fq_l", "type": "infile", "format": "sequence.fastq"},  # PE测序，所有样本fastq左端序列文件
            {"name": "fq_r", "type": "infile", "format": "sequence.fastq"},  # PE测序，所有样本fastq右端序列文件
            {"name": "fq_s", "type": "infile", "format": "sequence.fastq"},  # SE测序，所有样本fastq序列文件
            {"name": "cpu", "type": "int", "default": 6},  # trinity软件所分配的cpu数量
            {"name": "max_memory", "type": "string", "default": '100G'},  # trinity软件所分配的最大内存，单位为GB
            {"name": "min_contig_length", "type": "int", "default": 200},  # trinity报告出的最短的contig长度。默认为200
            {"name": "kmer_size", "type": "int", "default": 25},  # kmer 长度 <=32
            {"name": "lines", "type": "int", "default": 1000},  # 分布式命令每个文件行数
            {"name": "min_kmer_cov", "type": "int", "default": 5},
            {"name": "jaccard_clip", "type": "bool", "default": False}, #分割高密度基因区间基因
            {"name": "no_normalize_reads", "type": "bool", "default": False}, #不做reads均一化
            {"name": "normalize_max_read_cov", "type": "int", "default": 50}, #reads均一化覆盖倍数
            {"name": "SS_lib_type", "type": "string", "default": 'none'},  # reads的方向，成对的reads: RF or FR; 不成对的reads: F or R，默认情况下，不设置此2参数
            {"name": "gene_fa", "type": "outfile", "format": "denovo_rna_v2.trinity_fasta"},
            # {"name": "trinity_fa", "type": "outfile", "format": "denovo_rna_v2.trinity_fasta"},
            {"name": "edge_thr", "type": "float", "default": 0.1}, #butterfly 参数
            {"name": "flow_thr", "type": "float", "default": 0.04}, #butterfly 参数
            {"name": "max_number_of_paths_per_node_init", "type": "int", "default": 20}, #butterfly 参数
            {"name": "min_per_id_same_path", "type": "int", "default": 98}, #butterfly 参数
            {"name": "trinity_cmd", "type": "outfile", "format": "denovo_rna_v2.common"},
            {"name": "partion_path", "type": "outfile", "format": "denovo_rna_v2.common_dir"},
            {"name": "trinity_version", "type": "string", "default": "2.8.5"},
        ]
        self.add_option(options)
        self.step.add_steps("assemble")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)
        self._memory_increase_step = 30

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
            raise OptionError('必须设置测序类型：PE OR SE', code = "32006101")
        if self.option('fq_type') not in ['PE', 'SE']:
            raise OptionError('测序类型不在所给范围内', code = "32006102")
        if not self.option("fq_l").is_set and not self.option("fq_r").is_set and not self.option("fq_s").is_set:
            raise OptionError("必须设置PE测序输入文件或者SE测序输入文件", code = "32006103")
        if self.option("fq_type") == "PE" and not self.option("fq_r").is_set and not self.option("fq_l").is_set:
            raise OptionError("PE测序时需设置左端序列和右端序列输入文件", code = "32006104")
        if self.option("fq_type") == "SE" and not self.option("fq_s").is_set:
            raise OptionError("SE测序时需设置序列输入文件", code = "32006105")
        if self.option("SS_lib_type") != 'none' and self.option("SS_lib_type") not in ['F', 'R', 'FR', 'RF']:
            raise OptionError("所设reads方向不在范围值内", code = "32006106")
        if self.option('kmer_size') > 32 or self.option('kmer_size') < 1:
            raise OptionError("所设kmer_size不在范围内，请检查", code = "32006107")
        if self.option('min_kmer_cov') < 1:
            raise OptionError("所设min_kmer_cov不在范围内，请检查", code = "32006108")

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
            self._cpu = 12
            self._memory = '60G'
        elif file_size <= 1024 * 4 and file_size > 1024 * 0.5:
            self._cpu = 16
            self._memory = '80G'
        else:
            self._cpu = 20
            self._memory = '200G'
        self.option('cpu', self._cpu - 2)
        mem = str(int(self._memory.strip('G')) - 10) + 'G'
        self.option('max_memory', mem)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([

        ])
        result_dir.add_regexp_rules([

        ])
        super(Trinity2Agent, self).end()


class Trinity2Tool(Tool):
    def __init__(self, config):
        super(Trinity2Tool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self._version = "v1.0.1"
        if self.option("trinity_version") == "2.5.0":
            self.trinity_path = self.config.SOFTWARE_DIR + '/bioinfo/denovo_rna_v2/trinityrnaseq-Trinity-v2.5.0/'
        elif self.option("trinity_version") == "2.8.5":
            self.trinity_path = self.config.SOFTWARE_DIR + '/bioinfo/denovo_rna_v2/trinityrnaseq-2.8.5/'
        else:
            self.set_error("trinity version error", code="32006111")

        self.perl =  '/program/perl-5.24.0/bin/'
        # self.bowtie = self.config.SOFTWARE_DIR + '/bioinfo/denovo_rna_v2/bowtie2-2.3.3.1-linux-x86_64/'
        self.bowtie = self.config.SOFTWARE_DIR + '/bioinfo/align/bowtie2-2.3.4.3-linux-x86_64/'
        self.java = self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin'
        python_path = self.config.SOFTWARE_DIR + '/program/Python/bin/'
        self.set_environ(PATH=python_path)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + self.perl)
        self.set_environ(PATH=self.bowtie)
        self.set_environ(PATH=self.java)

        if self.option("trinity_version") == "2.8.5":
            self.salmon = self.config.SOFTWARE_DIR + '/bioinfo/rna/salmon-0.14.1/bin/'
            self.jellyfish = self.config.SOFTWARE_DIR + '/bioinfo/denovo_rna_v2/jellyfish-2.3.0/bin/'
            self.set_environ(PATH=self.salmon)
            self.set_environ(PATH=self.jellyfish)
            self.gcc = software_dir + '/gcc/5.1.0/bin'
            self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
            self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)


    def run(self):
        """
        运行
        :return:
        """
        super(Trinity2Tool, self).run()
        self.run_trinity()
        self.set_output()

    def run_trinity(self):
        """
        运行trinity软件，进行拼接组装
        """
        #trinity command
        if self.option('fq_type') == 'SE':
            cmd = '{}perl {}Trinity --seqType fq --single {} --no_distributed_trinity_exec'.format(
                self.perl, self.trinity_path, self.option('fq_s').prop['path'])
        else:
            cmd = '{}perl {}Trinity --seqType fq --left {} --right {} --no_distributed_trinity_exec'.format(
                self.perl, self.trinity_path, self.option('fq_l').prop['path'], self.option('fq_r').prop['path'])
        if self.option('SS_lib_type') != 'none':
            cmd += ' --SS_lib_type {}'.format(self.option('SS_lib_type'))

        bfly_CPU = 2
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

        command = self.add_command("trinity_cmd", cmd, ignore_error=True)
        for run_times in [1, 2, 3]:
            self.logger.info('运行trinity软件第{}次，进行组装拼接'.format(run_times))
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("trinity运行完成")
                break
            elif command.return_code in [2,255]:
                self.logger.info("return code: {}".format(command.return_code))
                self.add_state('memory_limit', 'memory is low!')
            else:
                self.logger.info("第{}次尝试，返回值为{}".format(run_times, command.return_code))
                if run_times == 3:
                    self.set_error("trinity三次运行都没有完成，退出!", code = "32006109")
                else:
                    self.logger.info("trinity运行出错,进行下一次尝试")

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
        #trinity_cmd = self.work_dir + '/trinity_out_dir/recursive_trinity.cmds'
        with open ('trinity_out_dir/recursive_trinity.cmds', 'r') as cmds:
            cmd_lines = cmds.readlines()
        line_num = len(cmd_lines)
        cmd_per_file = 1000
        if line_num > 300000:
            cmd_per_file = 3000
        if line_num > 500000:
            cmd_per_file = 5000
        file_num = int(line_num/cmd_per_file) + 1
        file_objects = []
        for i in range(1, file_num + 1):
            file_object  = open( os.path.join(self.output_dir, 'recursive_trinity.cmds'+str(i)), 'w')
            file_objects.append(file_object)
        for line in  range(0, line_num-1):
            file_index = line % file_num
            line_context = cmd_lines[line]
            line_context = re.subn(r'--max_memory [0-9]*G', '--max_memory 12G', line_context)[0]
            file_objects[file_index].write(line_context)
        for file_object in file_objects:
            file_object.close()

        # split_cmd = "split -d -l {}  trinity_out_dir/recursive_trinity.cmds  {}".format(self.option('lines'), os.path.join(self.output_dir, 'recursive_trinity.cmds'))
        # self.logger.info("分割文件{}".format(split_cmd))
        # os.system(split_cmd)
        self.logger.info("设置结果目录")

        try:
            #self.option('trinity_cmd', os.path.join(self.output_dir, 'recursive_trinity.cmds'))
            #self.logger.info("设置组装，命令结果目录成功 {}".format(self.option('trinity_cmd')))
            self.option('partion_path', os.path.join(self.work_dir, 'trinity_out_dir/read_partitions'))
            self.logger.info("设置组装，结果运行目录成功 {}".format(self.option('partion_path')))
            self.end()
        except Exception as e:
            self.logger.info("设置组装拼接分析结果目录失败{}".format(e))
            self.set_error("设置组装拼接分析结果目录失败%s", variables = (e), code = "32006110")
