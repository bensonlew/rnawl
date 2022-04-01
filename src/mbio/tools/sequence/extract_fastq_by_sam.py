# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import pandas as pd


class ExtractFastqBySamAgent(Agent):
    """
    ExtractFastqBySam:根据sam文件/sequence ID 挑选或剔除序列
    version 1.0
    author: zhujuan
    last_modify: 2017.09.15
    """
    def __init__(self, parent):
        super(ExtractFastqBySamAgent, self).__init__(parent)
        options = [
            {"name": "fq_type", "type": "string", "default": "PSE"},  # fq类型，PE、SE、PSE（即PE+SE，单端加双端）
            {"name": "sam", "type": "infile", "format": "align.bwa.sam_dir"},     # sam格式文件,内含对应list文件
            {"name": "extract_type", "type": "string", "default": "unmap"},  # 提取的fq结果类型是 'map'的还是'unmap'的
            {"name": "reasult_dir", "type": "outfile", 'format': "sequence.fastq_dir"}  # 输出文件夹
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        if self.option("fq_type") not in ['PE', 'SE', 'PSE']:
            raise OptionError("请传入fq类型，PE or SE or PSE?", code="34000601")
        if not self.option("sam").is_set:
            raise OptionError("请提供需要提取的文件sam路径", code="34000602")

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ])
        super(ExtractFastqBySamAgent, self).end()

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 1
        self._memory = '10G'


class ExtractFastqBySamTool(Tool):
    def __init__(self, config):
        super(ExtractFastqBySamTool, self).__init__(config)
        self.samples = self.get_sam_list()
        self.perl_path = '/program/perl-5.24.0/bin/perl'
        #self.perl_script = self.config.SOFTWARE_DIR + '/bioinfo/seq/scripts/get_fq_bysam.pl'
        self.perl_script = self.config.PACKAGE_DIR + '/sequence/scripts/get_fq_bysam.pl'

    def fastq_form_sam(self):
        n = 0  # add n to count by guhaidong 20170918
        samples = self.samples
        extract_list = os.path.join(self.output_dir, "list.txt")
        with open(extract_list, "wb") as w:
            for sample in samples:
                n += 1  # added by guhaidong 20170918
                if self.option("fq_type") in ["PE", "PSE"]:
                    sam_pe = os.path.join(self.option("sam").prop["path"], samples[sample]["pe"])
                    sam_pe_cmd = '{} {} {} {} {} {}'.format(self.perl_path, self.perl_script, sam_pe, 'pe',
                                                            self.option("extract_type"),
                                                            os.path.join(self.output_dir, sample))
                    # command = self.add_command('pe_cmd_{}'.format('pe'), sam_pe_cmd).run()
                    command = self.add_command('pe_cmd_{}'.format(n), sam_pe_cmd).run()  # modified by guhaidong 20170918
                    self.wait(command)
                    if command.return_code == 0:
                        self.logger.info("运行{}完成".format(command.name))
                        w.write(sample + '.1.fq\t'+sample+'\tl\n')
                        w.write(sample + '.2.fq\t'+sample+'\tr\n')
                    else:
                        self.set_error("运行%s出错", variables=(command.name), code="34000601")
                        raise Exception("运行{}出错".format(command.name))
                if self.option("fq_type") in ["SE", "PSE"]:
                    sam_se = os.path.join(self.option("sam").prop["path"], samples[sample]["se"])
                    sam_se_cmd = '{} {} {} {} {} {}'.format(self.perl_path, self.perl_script, sam_se, 'se',
                                                            self.option("extract_type"),
                                                            os.path.join(self.output_dir, sample))
                    # command = self.add_command('se_cmd_{}'.format('se'), sam_se_cmd).run()
                    command = self.add_command('se_cmd_{}'.format(n), sam_se_cmd).run()  # modified by guhaidong 20170918
                    self.wait(command)
                    if command.return_code == 0:
                        self.logger.info("运行{}完成".format(command.name))
                        w.write(sample + '.s.fq\t'+sample+'\ts\n')
                    else:
                        self.set_error("运行{}出错", variables=(command.name), code="34000602")
                        raise Exception("运行{}出错".format(command.name))

    def get_sam_list(self):
        list_path = os.path.join(self.option("sam").prop["path"], "list.txt")
        if os.path.exists(list_path):
            self.logger.info(list_path)
        sample = {}
        with open(list_path, "rb") as l:
            for line in l:
                line = line.strip().split()
                if line[1] not in sample:
                    sample[line[1]] = {line[2]: line[0]}
                else:
                    sample[line[1]][line[2]] = line[0]
        return sample

    def set_output(self):
        self.logger.info('开始设置输出结果文件')
        try:
            # self.option('reasult_dir', self.output_dir)
            self.option('reasult_dir').set_path(self.output_dir)  # modified by guhaidong 20170915
            self.logger.info("设置输出结果文件正常")
        except Exception as e:
            raise Exception("设置输出结果文件异常——{}".format(e))

    def fq_stat(self):
        stat_list = os.path.join(self.option("reasult_dir").prop['path'], "rehost_stat.list.txt")
        stat_file = pd.read_table(stat_list, sep='\t', header=None)
        stat_file.columns = ['samples', 'Total_Reads', 'Total_Bases', 'type']
        stat = stat_file.groupby(stat_file['samples']).sum()
        stat.to_csv(os.path.join(self.output_dir, "stat.list.txt"), sep="\t")

    def run(self):
        """
        运行
        """
        super(ExtractFastqBySamTool, self).run()
        if os.path.exists(self.output_dir + '/rehost_stat.list.txt'):
            os.remove(self.output_dir + '/rehost_stat.list.txt')
        self.fastq_form_sam()
        self.set_output()
        self.fq_stat()
        self.end()

