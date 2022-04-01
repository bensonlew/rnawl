# -*- coding: utf-8 -*-
# __author__ = 'hongyu.chen'

# import pysam
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import gzip
import os
import subprocess
import re
import datetime


class Fastq2bamAgent(Agent):
    def __init__(self, parent=None):
        super(Fastq2bamAgent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "string"},  # 输入F/M/S的fastq文件的样本名,fastq_gz_dir/WQ235F
            {"name": "is_dcpt", "type": "bool", "default": False},  # 是否是dcpt样品
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},  # hg38.chromosomal_assembly/ref.fa
            {"name": "targets_bedfile", "type": "infile", "format": "paternity_test_V2.rda"},
            {"name": "seq_path", "type": "infile", "format": "sequence.fastq_dir"},  # fastq所在路径
            {"name": "cpu_number", "type": "int", "default": 5},
            {"name": "memory", "type": "int", "default": 10}
        ]
        self.add_option(options)
        self.step.add_steps("Bam2tab")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.Bam2tab.start()
        self.step.update()

    def stepfinish(self):
        self.step.Bam2tab.finish()
        self.step.update()

    def check_options(self):
        if not self.option("fastq"):
            raise OptionError("必须输入fastq文件的样本名")
        # if not self.option("ref_fasta").is_set:
        #     raise OptionError("必须输入参考基因组的fastq文件")
        if not self.option('seq_path'):
            raise OptionError('必须提供fastq文件所在的路径')
        # if not self.option('targets_bedfile'):
        #     raise OptionError('必须提供target_bedfile文件')
        return True

    def set_resource(self):
        self._cpu = self.option("cpu_number")
        self._memory = "{}G".format(self.option("memory"))

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            [r".mem.sort.hit.filter.bam", "bam", "所有位点的信息"],
            [r".qc", "qc", "数据质控文件"]
        ])
        super(Fastq2bamAgent, self).end()


class Fastq2bamTool(Tool):
    def __init__(self, config):
        super(Fastq2bamTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/ruby-2.3.1') # 测试机
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/ruby-2.4.1/bin')  # 正式机
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/seq/bioawk')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/seq/seqtk-master')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/align/bwa-0.7.15')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/seq/samblaster-0.1.24')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/align/samtools-1.4/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/medical/bedtools-2.24.0/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/seq/bcftools-1.4/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/medical/picard-tools-2.2.4')
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/medical/prinseq-lite-0.20.4')
        self.picard_path = self.config.SOFTWARE_DIR + '/bioinfo/medical/picard-tools-2.2.4/picard.jar'
        self.java_path = 'program/sun_jdk1.8.0/bin/java'
        self.cmd_path = "bioinfo/medical/scripts/"
        self.samblaster = "bioinfo/seq/samblaster-0.1.24/samblaster"
        self.samtools = "bioinfo/align/samtools-1.4/bin/samtools"
        self.bedtools = "bioinfo/medical/bedtools-2.24.0/bin/bedtools"

        self.trimfq_left = 0

        self.script_path = "bioinfo/medical/scripts/pt_scripts"
        self.script_absolute_path = self.config.SOFTWARE_DIR + '/' + self.script_path

    def sixn_tag(self):
        """计算self.trim_left的值，为R1、R2添加X6:Z标签"""
        r1_path = self.option('seq_path').prop['path'] + '/' + self.option('fastq') + '_R1.fastq.gz'
        r2_path = self.option('seq_path').prop['path'] + '/' + self.option('fastq') + '_R2.fastq.gz'

        first6 = dict()
        first9 = dict()
        n = 0
        with gzip.open(r1_path, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    line = f.readline()
                    if line[6:12] not in first6:
                        first6[line[6:12]] = 0
                    else:
                        first6[line[6:12]] += 1
                    if line[6:15] not in first9:
                        first9[line[6:15]] = 0
                    else:
                        first9[line[6:15]] += 1
                    n += 1
                if n >= 1000:
                    break
        if max(first6.values()) >= 500:
            if max(first9.values()) >= 500:
                self.trimfq_left = 15
            else:
                self.trimfq_left = 12

        sample_name = self.option('fastq')
        r1_outfile = '{}.with6N_1.fq'.format(sample_name)
        out = open(r1_outfile, 'w')
        with gzip.open(r1_path, 'r') as f1:
            for line in f1:
                if re.match(r'^@(\S+)(\s+)(\S+)', line):
                    name = re.match(r'^@(\S+)(\s+)(\S+)', line).group(1)
                    seq = f1.readline()
                    strand = f1.readline()
                    qual = f1.readline()
                    tag = seq[0:6]
                    out.write('@{} X6:Z:{}\n'.format(name, tag))
                    out.write(seq + strand + qual)
        out.close()

        r2_outfile = '{}.with6N_2.fq'.format(sample_name)
        out1 = open(r2_outfile, 'w')
        with gzip.open(r2_path, 'r') as f2:
            for line in f2:
                if re.match(r'^@(\S+)(\s+)(\S+)', line):
                    name = re.match(r'^@(\S+)(\s+)(\S+)', line).group(1)
                    seq = f2.readline()
                    strand = f2.readline()
                    qual = f2.readline()
                    tag = seq[0:6]
                    out1.write('@{} X6:Z:{}\n'.format(name, tag))
                    out1.write(seq + strand + qual)
        out1.close()

    def fastq2bam(self):
        name = self.option('fastq')
        if not os.path.exists(self.script_absolute_path):
            os.mkdir(self.script_absolute_path)
        now_time = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
        full_path = self.script_absolute_path + '/' + name + '_fastq2bam.sh' + now_time
        path = self.script_path + '/' + name + '_fastq2bam.sh' + now_time
        with open(full_path, 'w') as out:
            header_line = '#!/bin/bash'
            cmd1 = 'seqtk mergepe {}.with6N_1.fq {}.with6N_2.fq > {}.merge.fq'.format(name, name, name)

            cmd2 = 'seqtk trimfq -b {} {}.merge.fq > {}.trim.merge.fq'.format(self.trimfq_left, name, name)

            cmd3 = "bwa mem -t {} -p -C -R \'@RG\tID:{}\tSM:{}\tPL:illumina\tPU:illumina\tLB:illumina\' {} " \
                   "{}.trim.merge.fq > {}.sam".format(
                    self.option('cpu_number')-1, name, name, self.option('ref_fasta').prop['path'], name, name)

            if self.option('is_dcpt'):
                cmd4 = 'cp -f {}.sam {}.markdup.sam'.format(name, name)
            else:
                cmd4 = 'samblaster -i {}.sam -o {}.markdup.sam'.format(name, name)

            cmd5 = 'samtools view -@ {} -Sb {}.markdup.sam -o {}.mem.bam'.format(self.option('cpu_number')-1, name, name)

            cmd6 = 'samtools sort -T {} -@ {} {}.mem.bam -o {}.mem.sort.bam'.format(
                                                                            name, self.option('cpu_number')-1, name, name)

            cmd7 = 'samtools index {}.mem.sort.bam'.format(name)

            cmd8 = 'bedtools intersect -abam {}.mem.sort.bam -b {} > {}.mem.sort.hit.bam'.format(
                name, self.option('targets_bedfile').prop['path'], name)

            out.write('{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(
                header_line, cmd1, cmd2, cmd3, cmd4, cmd5, cmd6, cmd7, cmd8))
        self.run_single_command('chmod +x {}'.format(full_path))
        return path, full_path

    def sam_filter(self):
        """
        __author__ = 'kefei.huang'
        程序是对sam处理，
        1. 删除没有序列信息的
        2. 删除比对到了softcliping的
        3. 删除mapping quality < 30的
        :return:
        """
        # 由于线上与线下的pysam版本不一致，所以上线的时候 要启用 200~201行代码  add by hd 20180112
        import sys
        sys.path = ["/mnt/ilustre/users/sanger/app/program/Python/lib/python2.7/site-packages/test123"] + sys.path
        import pysam
        name = self.option('fastq')
        samfile = pysam.AlignmentFile("{}.mem.sort.hit.bam".format(name), "rb")
        outfile = pysam.AlignmentFile("{}.mem.sort.hit.filter.bam".format(name), "wb", template=samfile)
        for line in samfile:
            if line.alen < 10:
                continue
            if line.mapping_quality < 30:
                continue
            if re.search(r"S", line.cigarstring):
                continue
            outfile.write(line)

        samfile.close()
        outfile.close()

    def run_fastq2bam(self):
        self.sixn_tag()

        cmd = self.fastq2bam()
        self.logger.info('开始运行fastq2bam')
        command = self.add_command('fastq2bam', cmd[0])
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行fastq2bam完成")
        else:
            if command.return_code is None:
                command.rerun()
                self.wait()
                if command.return_code == 0:
                    self.logger.info("重新运行fastq2bam完成")
                elif command.return_code is None:
                    self.logger.info("重新运行fastq2bam返回码仍然为None")
                else:
                    self.logger.info('Run Return Code: {}'.format(command.return_code))
                    raise Exception("fastq2bam运行出错!")
            else:
                self.logger.info('Run Return Code: {}'.format(command.return_code))
                raise Exception("fastq2bam运行出错!")
        os.system('rm {}'.format(cmd[1]))
        self.sam_filter()

    def run_single_command(self, cmd):
        exitcode = subprocess.call(cmd, shell=True)
        if exitcode == 0:
            self.logger.info(" 运行成功：" + cmd)
        else:
            exitcode = subprocess.call(cmd, shell=True)
            if exitcode == 0:
                self.logger.info(" 重运行成功：" + cmd)
            else:
                raise Exception("运行失败：" + cmd)

    def set_output(self):
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        name = self.option('fastq')
        file0 = '/{}.mem.sort.hit.filter.bam'.format(name)
        file1 = '/{}.mem.sort.bam'.format(name)
        file2 = '/{}.mem.sort.hit.bam'.format(name)
        os.link(self.work_dir + file0, self.output_dir + file0)
        os.link(self.work_dir + file1, self.output_dir + file1)
        os.link(self.work_dir + file2, self.output_dir + file2)
        self.logger.info("设置结果目录：成功")
        os.system('rm *.bam *.sam *.fq')                              # 清理work_dir下的大文件

    def run(self):
        super(Fastq2bamTool, self).run()
        self.run_fastq2bam()
        self.set_output()
        self.end()
