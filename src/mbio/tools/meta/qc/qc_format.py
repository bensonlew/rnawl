# -*- coding: utf-8 -*-
# __author__ = 'xuting'
from __future__ import division
import math
import traceback
import os
import subprocess
import re
import gzip
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.files.sequence.fastq_dir import FastqDirFile
from mbio.files.sequence.fasta_dir import FastaDirFile


class QcFormatAgent(Agent):
    """
    author: xuting
    last_modify: 2015.11.06
    接收fastq文件或者fastq文件夹，用于格式化用户的输入和产生一个fasta文件共下游OTU分析使用
    """
    def __init__(self, parent):
        super(QcFormatAgent, self).__init__(parent)
        options = [
            {'name': 'in_fastq', 'type': 'infile', 'format': 'sequence.fastq, sequence.fastq_dir'},  # 输入的fastq文件
            {'name': 'otu_fasta', 'type': 'outfile', 'format': 'sequence.fasta'},  # 输出的合并到一起的fasta，供后续的otu分析用
            {'name': 'renamed_fastq_dir', 'type': 'outfile', 'format': 'sequence.fastq_dir'},  # 按样本名进行重命名或者拆分的fastq文件夹
            {'name': 'fasta_dir', 'type': 'outfile', 'format': 'sequence.fasta_dir'}]  # 由fastq文件夹转化而来的fasta文件
        self.add_option(options)
        self.step.add_steps("fastq_format")
        self.on('start', self.start_qc)
        self.on('end', self.end_qc)

    def start_qc(self):
        self.step.fastq_format.start()
        self.step.update()

    def end_qc(self):
        self.step.fastq_format.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        """
        if not self.option('in_fastq').is_set:
            raise OptionError("必须输入in_fastq参数")
        self.option('in_fastq').get_info()
        if self.get_option_object("in_fastq").format == 'sequence.fastq_dir':
            if not self.option('in_fastq').prop['has_list_file']:
                raise OptionError('fastq文件夹中必须含有一个名为list.txt的文件名--样本名的对应文件')
        if self.get_option_object('in_fastq').format == 'sequence.fastq':
            if not self.option('in_fastq').prop['has_sample_info']:
                raise OptionError("fastq文件中必须在序列名中带有样本名称(以下划线分隔)")

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["fastq_dir", "sequence.fastq_dir", "样本所对应的fastq目录文件夹"],
            ["converted_fastas", "sequence.fasta_dir", "从fastq转换而来的fasta目录文件夹"],
            ["cat_meta.fasta", "sequence.fasta", "所有样本的fasta合并到一起的fasta文件"]
        ])
        result_dir.add_regexp_rules([
            ['\./fastq_dir/.+\.fastq$', "sequence.fastq", "样本对应的fastq文件"],
            ['\./converted_fastas/.+\.fasta$', "sequence.fasta", "由fastq转化而来的fasta文件"]
        ])
        super(QcFormatAgent, self).end()

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 10
        total = 0
        if self.get_option_object("in_fastq").format == 'sequence.fastq_dir':
            for f in self.option("in_fastq").prop["fastq_basename"]:
                f = os.path.join(self.option("in_fastq").prop["path"],f)
                total += os.path.getsize(f)
        if self.get_option_object("in_fastq").format == 'sequence.fastq':
            total = os.path.getsize(self.option("in_fastq").prop["path"])
        total = total / (1024 * 1024 * 1024)
        total = total * 4
        total = math.ceil(total)
        self._memory = '{}G'.format(int(total))


class QcFormatTool(Tool):
    def __init__(self, config):
        super(QcFormatTool, self).__init__(config)
        self.fastq_dir = os.path.join(self.work_dir, "output", "fastq_dir")
        self.fasta_dir = ''
        self.fasta = ''
        self.seq_sample = dict()

    def rename_fastq(self):
        """
        输入是文件夹的时候，根据文件名，样本名的对应文件，重命名文件
        """
        filename_sample = os.path.join(self.option('in_fastq').prop['path'], "list.txt")
        with open(filename_sample, 'r') as r:
            for line in r:
                line = line.rstrip('\n')
                line = re.split('\t', line)
                if re.search(r'\.(fastq|fq)\.gz', line[0]):
                    gz_file = os.path.join(self.option('in_fastq').prop['path'], line[0])
                    target_file = os.path.join(self.fastq_dir, line[1] + ".fastq")
                    try:
                        subprocess.check_call('gunzip -c ' + gz_file + " >> " + target_file, shell=True)
                    except subprocess.CalledProcessError:
                        self.set_error("解压缩文件失败!检查输入是否是正确的gz文件")
                else:
                    file_ = os.path.join(self.option('in_fastq').prop['path'], line[0])
                    target_file = os.path.join(self.fastq_dir, line[1] + ".fastq")
                    with open(file_, "r") as f:
                        txt = f.read()
                    with open(target_file, "a") as a:
                        a.write(txt)

    def seprate_fastq(self):
        """
        无对应文件，拆分fastq
        """
        self.option('in_fastq').get_info()
        if self.option('in_fastq').prop["is_gz"]:
            seq_name2file_name = self._open_gzip()
        else:
            seq_name2file_name = self._open_file()
        handler = dict()
        self.logger.debug("正在准备文件...")
        all_file = list(set(seq_name2file_name.values()))
        """
        for v in seq_name2file_name.itervalues():
            c += 1
            handler[v] = open(v, 'wb')
            if c % 10000 == 0:
                self.logger.debug("正在遍历第{}个key".format(c))
        self.logger.debug(handler.iteritems())
        """
        self.logger.debug(all_file)
        for l in all_file:
            handler[l] = open(l, 'wb')
        count = 0
        with open(self.option('in_fastq').prop['path'], 'r') as f:
            for line in f:
                count += 1
                line = line.rstrip('\r\n')
                name = re.split('\s+', line)[0]
                name = re.sub(r'@', '', name)
                head = name
                line = re.split(r'_', name)
                handler[seq_name2file_name[head]].write("@" + line[-1] + "\n")
                for i in range(1, 4):
                    line = f.next()
                    handler[seq_name2file_name[head]].write(line)
                if count % 10000 == 0:
                    self.logger.info("正在输出第" + str(count) + "条序列")
        for v in seq_name2file_name.itervalues():
            handler[v].close()
        self.logger.info("fastq 文件拆分完毕 ")

    def _open_file(self):
        warninglog = False
        count = 0
        seq_name2file_name = dict()  # 记录序列名应该被输出到哪一个文件里面去
        self.logger.info("正在遍历输入的fastq文件")
        with open(self.option('in_fastq').prop['path'], 'r') as f:
            for line in f:
                count += 1
                line = line.rstrip('\r\n')
                name = re.split('\s+', line)[0]
                name = re.sub(r'@', '', name)
                head = name
                line = re.split(r'_', name)
                if len(line) > 2:
                    warninglog = True
                line.pop(-1)
                filename = "_".join(line)
                filename = os.path.join(self.fastq_dir, filename + ".fastq")
                seq_name2file_name[head] = filename
                f.next()
                f.next()
                f.next()
                if count % 10000 == 0:
                    self.logger.info("正在遍历第" + str(count) + "条序列")
                """
                with open(filename, 'a') as a:
                    a.write("@" + head + "\n")
                    for i in range(1, 4):
                        line = f.next()
                        a.write(line)
                if count % 10000 == 0:
                    self.logger.info("正在输出第" + str(count) + "条序列")
        self.logger.info("fastq 文件拆分完毕 ")
                """
        if warninglog:
            self.logger.warning("fastq文件里包含有两个以上的下划线，程序将取最后一个下划线之前的所有内容作为样本名！")
        return seq_name2file_name

    def _open_gzip(self):
        warninglog = False
        count = 0
        self.logger.info("正在遍历输入的fastq文件")
        seq_name2file_name = dict()  # 记录序列名应该被输出到哪一个文件里面去
        with gzip.open(self.option('in_fastq').prop['path'], 'r') as f:
            for line in f:
                count += 1
                line = line.rstrip('\r\n')
                name = re.split('\s+', line)[0]
                name = re.sub(r'@', '', name)
                line = re.split(r'_', name)
                if len(line) > 2:
                    warninglog = True
                head = line[-1]
                line.pop(-1)
                filename = "_".join(line)
                filename = os.path.join(self.fastq_dir, filename + ".fastq")
                seq_name2file_name[head] = filename
                f.next()
                f.next()
                f.next()
                if count % 10000 == 0:
                    self.logger.info("正在遍历第" + str(count) + "条序列")
                """
                with open(filename, 'a') as a:
                    a.write("@" + head + "\n")
                    for i in range(1, 4):
                        line = f.next()
                        a.write(line)
                if count % 10000 == 0:
                    self.logger.info("正在输出第" + str(count) + "条序列")
        self.logger.info("fastq 文件拆分完毕 ")
        """
        if warninglog:
            self.logger.warning("fastq文件里包含有两个以上的下划线，程序将取最后一个下划线之前的所有内容作为样本名！")
        return seq_name2file_name

    def get_fastq_dir(self):
        """
        生成fastq文件夹
        """
        if not os.path.exists(self.fastq_dir):
            os.mkdir(self.fastq_dir)
        if self.get_option_object("in_fastq").format == 'sequence.fastq_dir':
            self.rename_fastq()
        if self.get_option_object("in_fastq").format == 'sequence.fastq':
            self.seprate_fastq()

    def get_fasta_dir(self):
        """
        生成fasta文件夹
        """
        fq_dir = FastqDirFile()
        fq_dir.set_path(self.fastq_dir)
        fq_dir.get_full_info(os.path.join(self.work_dir, "output"))
        try:
            self.fasta_dir = fq_dir.covert_to_fasta()
            self.logger.info("fasta 文件夹生成完毕")
        except Exception:
            self.logger.error("fastq转化fasta失败！" + traceback.format_exc())
            raise Exception("软件fastq_to_fasta运行出错，请查看输入的fastq是否正确")

    def get_fasta(self):
        """
        生成fasta供otu分析
        """
        fa = FastaDirFile()
        fa.set_path(self.fasta_dir)
        fa.get_full_info(os.path.join(self.work_dir, "output"))
        try:
            self.fasta = fa.cat_fastas_for_meta()
            self.logger.info("后续OTU分析的fasta文件生成完毕")
        except Exception:
            self.set_error("生成fasta失败！")

    def set_output(self):
        self.logger.info("set output begin")
        self.option('otu_fasta').set_path(self.fasta)
        self.option('renamed_fastq_dir').set_path(self.fastq_dir)
        self.option('renamed_fastq_dir').check()
        self.option('fasta_dir').set_path(self.fasta_dir)
        self.option('fasta_dir').check()

    def run(self):
        """
        运行
        """
        super(QcFormatTool, self).run()
        self.get_fastq_dir()
        self.get_fasta_dir()
        self.get_fasta()
        self.set_output()
        self.logger.info("程序完成，即将退出")
        self.end()
