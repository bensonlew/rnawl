# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict


class MicrobialGenomeQcModule(Module):
    """
    Fastp质控
    author: wangzhaoyue
    last_modify: 2017.12.13
    """

    def __init__(self, work_id):
        super(MicrobialGenomeQcModule, self).__init__(work_id)
        options = [
            {"name": "sample_path", "type": "infile", "format": "sequence.file_sample"},
            # fastq路径list.txt文件，第一列路径，第二列样本名，第三列序列类型 l or r
            {'name': "sample_info", "type": "infile", "format": "sequence.barcode_info"},  # 样本信息表
            {'name': 'flag', 'type': "string", "default": "4"},  # 提取没有比对上的reads,此处固定取值4，软件默认0
            {'name': 'readl', "type": "string"},  # 切除序列的阈值
            {'name': 'illuminaclip', 'type': "string", "default": "2:30:10"},  # 2:30:10
            {'name': 'leading', 'type': "string", "default": "3"},  # 切除首端碱基质量小于0的碱基或者N
            {'name': 'tailing', 'type': "string", "default": "3"},  # 切除末端碱基质量小于20的碱基或者N
            {'name': 'sliding_window', 'type': "string", "default": "4:15"},
            # 例50:20  Windows的size是50个碱基，其平均碱基质量小于20，则切除
            {'name': 'minlen', 'type': "string", "default": "36"},  # 最低reads长度
            {"name": "seqprep_quality", "type": "string", "default": '20'},
            {"name": "seqprep_length", "type": "string", "default": '25'},
            {"name": "adapter_a", "type": "string", "default": "AGATCGGAAGAGCACACGTC"},
            {"name": "adapter_b", "type": "string", "default": "AGATCGGAAGAGCGTCGTGT"},
            {"name": "sickle_quality", "type": "string", "default": '20'},
            {"name": "sickle_length", "type": "string", "default": '20'},
            {"name": "qual_type", "type": "string", "default": 'sanger'},

        ]
        self.sample_path = defaultdict(list)
        self.sample_info = {}
        self.modules = []
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('sample_path'):
            raise OptionError('必须输入样本文件夹对应的路径信息')
        row_num = len(open(self.option("sample_path").prop['path'], "r").readline().split())
        if row_num != 3:
            raise OptionError("PE序列list文件应该包括文件名、样本名和左右端说明三列")
        if not self.option('sample_info'):
            raise OptionError('必须输入样本信息文件，包括样本，文库类型、插入片段长度三列')
        return True

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def single_microbial_genome_qc_run(self):
        n = 0
        for sample in self.sample_path:
            self.single_microbial_genome_qc = self.add_module('datasplit.single_microbial_genome_qc')
            self.step.add_steps('single_microbial_genome_qc{}'.format(n))
            opts = {
                "fq1": self.sample_path[sample][0],
                "fq2": self.sample_path[sample][1],
                "sample_name": sample,
                "insert_size": self.sample_info[sample],
                "flag": self.option('flag'),
                "illuminaclip": self.option('illuminaclip'),
                "leading": self.option('leading'),
                "tailing": self.option('tailing'),
                "sliding_window": self.option('sliding_window'),
                "minlen": self.option('minlen'),
                "seqprep_quality": self.option('seqprep_quality'),
                "seqprep_length": self.option('seqprep_length'),
                "adapter_a": self.option('adapter_a'),
                "adapter_b": self.option('adapter_b'),
                "sickle_quality": self.option('sickle_quality'),
                "sickle_length": self.option('sickle_length'),
                "qual_type": self.option('qual_type'),
            }
            if self.option('readl'):
                opts.update({'readl': self.option('readl')})
            self.single_microbial_genome_qc.set_options(opts)
            step = getattr(self.step, 'single_microbial_genome_qc{}'.format(n))
            step.start()
            self.step.update()
            self.single_microbial_genome_qc.on('end', self.finish_update, 'single_microbial_genome_qc{}'.format(n))
            self.modules.append(self.single_microbial_genome_qc)
            n += 1
        self.logger.info(self.modules)
        self.on_rely(self.modules, self.set_output)
        self.step.update()
        for module in self.modules:
            module.run()

    def run(self):
        """
        运行
        :return:
        """
        self.get_info()
        time.sleep(2)
        self.single_microbial_genome_qc_run()
        super(MicrobialGenomeQcModule, self).run()

    def link_has_dir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件一层目录文件夹移动到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.work_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        for file_dir in allfiles:
            dir_path = os.path.join(dirpath, file_dir)
            allfqs = os.listdir(dir_path)
            for fq in allfqs:
                new_path = os.path.join(newdir, file_dir)
                if not os.path.exists(new_path):
                    os.mkdir(new_path)
                oldfile = os.path.join(dir_path, fq)
                newfile = os.path.join(new_path, fq)
                if os.path.exists(newfile):
                    os.remove(newfile)
                size = os.path.getsize(oldfile)
                if size > 0:
                    os.link(oldfile, newfile)
                else:
                    self.logger.info("{}文件为空".format(oldfile))

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for module in self.modules:
            module_name = str(module.work_dir).strip().split("/")[-1]
            self.link_has_dir(module.output_dir, self.output_dir + '/' + module_name)
        self.logger.info("设置结果目录成功")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(MicrobialGenomeQcModule, self).end()

    def get_info(self):
        """
        获得样本对应的路径信息，以及样本的
        :return:
        """
        with open(self.option("sample_info").prop['path'])as f:
            lines = f.readlines()
            for line in lines[1:]:
                tmp = line.strip().split('\t')
                self.sample_info[tmp[0]] = tmp[2]  # 样本的插入片段长度
        with open(self.option('sample_path').prop['path'])as fr:
            for line in fr:
                self.logger.info(self.sample_path)
                tmp = line.strip().split('\t')
                if tmp[1] in self.sample_info.keys():
                    if tmp[1] in self.sample_path.keys():
                        if tmp[2] == 'l':
                            self.sample_path[tmp[1]].insert(0, tmp[0])
                        else:
                            self.sample_path[tmp[1]].append(tmp[0])
                    else:
                        self.sample_path[tmp[1]].append(tmp[0])
                else:
                    raise Exception('需要质控的序列样本{}没有相关的样本信息，请核实！'.format(tmp[1]))
            for key in self.sample_path.keys():
                if len(self.sample_path[key]) > 2:
                    raise Exception('需要质控的序列样本{}有重名，请改样本名或分开质控！'.format(key))
                elif len(self.sample_info[key]) < 2:
                    raise Exception('样本{}对应的R1,R2序列不全,请核实！'.format(key))
