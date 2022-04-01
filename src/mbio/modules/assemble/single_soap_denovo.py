# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

import os
from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class SingleSoapDenovoModule(Module):
    """
    宏基因运用SOAPdenovo2进行单个样本单个kmer组装
    author: wangzhaoyue
    last_modify: 2017.06.05
    """

    def __init__(self, work_id):
        super(SingleSoapDenovoModule, self).__init__(work_id)
        options = [
            {"name": "fastq1", "type": "infile", "format": "sequence.fastq"},  # 输入文件,sample.sickle.l.fastq
            {"name": "fastq2", "type": "infile", "format": "sequence.fastq"},  # 输入文件,sample.sickle.r.fastq
            {"name": "fastqs", "type": "infile", "format": "sequence.fastq"},  # 输入文件,sample.sickle.s.fastq
            {"name": "sample_name", "type": "string"},  # 样品名称
            {"name": "mem", "type": "int", "default": 100},  # 拼接内存
            {"name": "max_rd_len", "type": "string"},  # read最大读长
            {"name": "insert_size", "type": "string"},  # 平均插入片段长度
            {"name": "reverse_seq", "type": "string", "default": "0"},  # 配置文件的其他参数
            {"name": "asm_flags", "type": "string", "default": "3"},  # 配置文件的其他参数
            {"name": "rank", "type": "string", "default": "1"},  # 配置文件的其他参数
            {"name": "kmer", "type": "string"},  # k_mer值，例"39"
            {"name": "min_contig", "type": "string", "default": "500"},  # 输入最短contig长度，默认500
            {"name": "scafSeq", "type": "outfile", "format": "sequence.fasta"},  # 输出文件,sample.scafSeq
            {"name": "scaftig", "type": "outfile", "format": "sequence.fasta"},  # 输出文件，scaffold去掉N后的序列
            {"name": "cut_more_scaftig", "type": "outfile", "format": "sequence.fasta"},
            # 输出文件，去掉小于最短contig长度的序列
        ]
        self.add_option(options)
        # self.tools = []
        self.sum_tools = []
        self.step.add_steps("SOAPdenovo2", "GetContig")
        # self.on('start', self.stepstart)
        # self.on('end', self.stepfinish)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('fastq1'):
            raise OptionError('必须输入*l.fastq文件', code="21300401")
        if not self.option('fastq2'):
            raise OptionError('必须输入*r.fastq文件', code="21300402")
        if not self.option('max_rd_len'):
            raise OptionError('必须输入read的最大长度', code="21300403")
        if not self.option('insert_size'):
            raise OptionError('必须输入平均插入片段的长度', code="21300404")
        return True

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def SOAPdenovo2_run(self):
        self.SOAPdenovo2 = self.add_tool('assemble.soap_denovo')
        opts = ({
            "fastq1": self.option('fastq1'),
            "fastq2": self.option('fastq2'),
            "sample_name": self.option('sample_name'),
            "max_rd_len": self.option('max_rd_len'),
            "mem": self.option('mem'),
            "insert_size": self.option('insert_size'),
            "reverse_seq": self.option('reverse_seq'),
            "asm_flags": self.option('asm_flags'),
            "rank": self.option('rank'),
            "kmer": self.option('kmer'),
        })
        if self.option('fastqs'):
            opts['fastqs'] = self.option('fastqs')
        self.SOAPdenovo2.set_options(opts)
        self.SOAPdenovo2.on('end', self.get_contig_run)
        self.SOAPdenovo2.run()
        self.sum_tools.append(self.SOAPdenovo2)
        self.step.SOAPdenovo2.finish()
        self.step.GetContig.start()
        self.step.update()

    def get_contig_run(self):
        self.get_contig = self.add_tool('assemble.get_contig')
        self.get_contig.set_options({
            "scafSeq": self.SOAPdenovo2.option('scafSeq'),
            "min_contig": self.option('min_contig'),
        })
        self.get_contig.on('end', self.set_output)
        self.get_contig.run()
        self.sum_tools.append(self.get_contig)
        self.step.GetContig.finish()
        self.step.update()

    def run(self):
        """
        运行
        :return:
        """
        super(SingleSoapDenovoModule, self).run()
        self.SOAPdenovo2_run()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.work_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
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

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for tool in self.sum_tools:
            self.linkdir(tool.output_dir, self.output_dir)
        self.logger.info("设置SingleSoapDenovo结果目录成功")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(SingleSoapDenovoModule, self).end()
