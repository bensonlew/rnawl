# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict
class BacQcModule(Module):
    """
    微生物基因组质控
    author: gaohao
    last_modify: 2018.03.21
    """
    def __init__(self, work_id):
        super(BacQcModule, self).__init__(work_id)
        options = [
            {'name': 'fq1', 'type': "infile", "format": "sequence.fastq"},  # 原始数据
            {'name': 'fq2', 'type': "infile", "format": "sequence.fastq"},  # 原始数据
            {"name": 'rm_single', 'type': 'bool', "default": False},  # 是否舍弃掉single端
            {"name": "phix_tool", "type": "string", "default": "bwa", "choose": ["bwa", "bowtie"]},
            {'name': 'sample_name', "type": "string"},  # 样本名
            {'name': 'insert_size', "type": "string"},  # 插入片段长度
            {'name': 'flag', 'type': "string", "default": "4"},  # 提取没有比对上的reads,此处固定取值4，软件默认0
            {'name': 'readl', "type": "string"},  # 切除序列的阈值
            {'name': 'illuminaclip', 'type': "string", "default": "2:30:10"},  # 2:30:10
            {'name': 'leading', 'type': "string", "default": "3"},  # 切除首端碱基质量小于0的碱基或者N
            {'name': 'tailing', 'type': "string", "default": "3"},  # 切除末端碱基质量小于20的碱基或者N
            {'name': 'sliding_window', 'type': "string", "default": "4:15"},# 例50:20  Windows的size是50个碱基，其平均碱基质量小于20，则切除
            {'name': 'minlen', 'type': "string", "default": "36"},  # 最低reads长度
            {"name": "seqprep_quality", "type": "string", "default": '20'},
            {"name": "seqprep_length", "type": "string", "default": '25'},
            {"name": "adapter_a", "type": "string", "default": "AGATCGGAAGAGCACACGTC"},
            {"name": "adapter_b", "type": "string", "default": "AGATCGGAAGAGCGTCGTGT"},
            {"name": "sickle_quality", "type": "string", "default": '20'},
            {"name": "sickle_length", "type": "string", "default": '20'},
            {"name": "qual_type", "type": "string", "default": 'sanger'},
        ]
        self.step.add_steps('phix_filter', 'trimmomatic', 'seqprep', 'sickle')
        self.phix_filter = self.add_tool('bacgenome.phix_filter')
        self.trimmomatic = self.add_tool('datasplit.trimmomatic')
        self.seqprep = self.add_tool('datasplit.seq_prep')
        self.sickle = self.add_tool('bacgenome.sickle_mg')
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('fq1'):
            raise OptionError('必须输入1端序列', code="21400701")
        if not self.option('fq2'):
            raise OptionError('必须输入2端序列', code="21400702")
        if not self.option('sample_name'):
            raise OptionError('必须输入样本名', code="21400703")
        if not self.option('insert_size'):
            raise OptionError('必须输入插入片段长度', code="21400704")
        return True
    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()
    def phix_filter_run(self):
        opts = {
            "fq1": self.option('fq1'),
            "fq2": self.option('fq2'),
            "phix_tool": self.option("phix_tool"),
            "sample_name": self.option('sample_name'),
            "insert_size": self.option('insert_size'),
            "flag": self.option('flag'),
            'readl': self.option('readl'),
        }
        self.phix_filter.set_options(opts)
        self.phix_filter.on('end', self.set_output, 'phix_filter')
        self.phix_filter.run()
        self.step.phix_filter.finish()
        self.step.trimmomatic.start()
        self.step.update()
    def trimmomatic_run(self):
        self.trimmomatic.set_options({
            "fq1": self.phix_filter.option('out_fq1'),
            "fq2": self.phix_filter.option('out_fq2'),
            "fq_type": "PE",
            "illuminaclip": self.option('illuminaclip'),
            "leading": self.option('leading'),
            "tailing": self.option('tailing'),
            "sliding_window": self.option('sliding_window'),
            "minlen": self.option('minlen'),
            "lib_name": self.option('sample_name'),
        })
        self.trimmomatic.on('end', self.set_output, 'trimmomatic')
        self.trimmomatic.run()
        self.step.trimmomatic.finish()
        self.step.seqprep.start()
        self.step.update()
    def seqprep_run(self):
        opts = {
            "fastq_l": self.trimmomatic.option('out_fq1'),
            "fastq_r": self.trimmomatic.option('out_fq2'),
            "quality": int(self.option('seqprep_quality')),
            "length": int(self.option('seqprep_length')),
            "adapter_a": self.option('adapter_a'),
            "adapter_b": self.option('adapter_b'),
        }
        self.seqprep.set_options(opts)
        self.seqprep.on('end', self.set_output, 'seqprep')
        self.seqprep.run()
        self.step.seqprep.finish()
        self.step.update()
    def sickle_run(self):
        opts = {
            "fastq_l": self.seqprep.option('seqprep_l'),
            "fastq_r": self.seqprep.option('seqprep_r'),
            "quality": self.option('sickle_quality'),
            "length": self.option('sickle_length'),
            "qual_type": self.option('qual_type'),
            "fq_type": "PE",
        }
        self.sickle.set_options(opts)
        self.sickle.on('end', self.set_output, 'sickle')
        self.sickle.run()
        self.step.sickle.finish()
        self.step.update()

    def run(self):
        """
        运行
        :return:
        """
        super(BacQcModule, self).run()
        self.phix_filter.on('end', self.trimmomatic_run)
        self.trimmomatic.on('end', self.seqprep_run)
        self.seqprep.on('end', self.sickle_run)
        self.phix_filter_run()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        newdir = os.path.join(self.work_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        allfiles = os.listdir(dirpath)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.link(oldfiles[i], newdir)
    def set_output(self, event):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if event["data"] == "phix_filter":
            self.linkdir(self.phix_filter.output_dir, self.output_dir + '/phix_filter')
        if event["data"] == "trimmomatic":
            self.linkdir(self.trimmomatic.output_dir, self.output_dir + '/trimmomatic')
            unpair_fq = []
            all_files = os.listdir(self.trimmomatic.work_dir)
            for fq in all_files:
                if re.search(r'\.unpair\.', fq):
                    unpair_fq.append(self.trimmomatic.work_dir + '/' + fq)
            self.logger.info(unpair_fq)
            os.system('cat {} {} > {}'.format(unpair_fq[0], unpair_fq[1], self.trimmomatic.work_dir + '/' + self.option(
                'sample_name') + '.trimmomatic.unpaired.fq'))
        if event["data"] == "seqprep":
            self.linkdir(self.seqprep.output_dir, self.output_dir + '/seqprep')
        if event["data"] == "sickle":
            if not os.path.exists(self.output_dir + '/sickle/'):
                os.mkdir(self.output_dir + '/sickle/')
            if os.path.exists(self.output_dir + '/sickle/' + self.option('sample_name') + '.clean.1.fq'):
                os.remove(self.output_dir + '/sickle/' + self.option('sample_name') + '.clean.1.fq')
            os.link(self.sickle.output_dir + '/sickle_l.fastq',
                      self.output_dir + '/sickle/' + self.option('sample_name') + '.clean.1.fq')
            if os.path.exists(self.output_dir + '/sickle/' + self.option('sample_name') + '.clean.2.fq'):
                os.remove(self.output_dir + '/sickle/' + self.option('sample_name') + '.clean.2.fq')
            os.link(self.sickle.output_dir + '/sickle_r.fastq',
                      self.output_dir + '/sickle/' + self.option('sample_name') + '.clean.2.fq')
            if not self.option("rm_single"):
                os.system('cat {} {} > {}'.format(
                    self.trimmomatic.work_dir + '/' + self.option('sample_name') + '.trimmomatic.unpaired.fq',
                    self.sickle.output_dir + '/' + 'sickle.unpaired.fq',
                    self.output_dir + '/sickle/' + self.option('sample_name') + '.unpaired.fq'))
                if os.path.getsize(self.output_dir + '/sickle/' + self.option('sample_name') + '.unpaired.fq') == 0:  #add 202005
                    os.remove(self.output_dir + '/sickle/' + self.option('sample_name') + '.unpaired.fq')

            self.logger.info("设置结果目录成功")
            self.end()
    def end(self):
        super(BacQcModule, self).end()