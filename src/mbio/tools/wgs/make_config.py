# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180423

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import datetime
import random
import time
import os
import re


class MakeConfigAgent(Agent):
    """
    数据组装接口--将fastq文件夹中的所有的fastq文件整理成一个fastq 列表，然后在生成soap的运行配置文件
    lasted modified by hongdong@20190311 添加是否将同一个样本的不同的区域合并在一起进行组装的判断， 默认合并
    """
    def __init__(self, parent):
        super(MakeConfigAgent, self).__init__(parent)
        options = [
            {"name": "fastq_dir", "type": "string"},
            {"name": "fastq_list", "type": "string"},  # 增加c参数用于ssr生成配置文件
            {"name": "make_method", "type": 'string', 'default': "combine"}
        ]
        self.add_option(options)
        self.step.add_steps('snpeff')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.snpeff.start()
        self.step.update()

    def step_end(self):
        self.step.snpeff.finish()
        self.step.update()

    def check_options(self):
        if not self.option("fastq_dir") and not self.option("fastq_list"):
            raise OptionError("缺少fastq_dir参数或fastq_list参数", code="34503801")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(MakeConfigAgent, self).end()


class MakeConfigTool(Tool):
    def __init__(self, config):
        super(MakeConfigTool, self).__init__(config)
        self.sample_ids = []
        self.perl_path = "miniconda2/bin/perl"
        self.config_path = self.config.PACKAGE_DIR + "/wgs/denovo.make.config.pl"

    def get_fastq_list(self):
        """
        根据文件夹中的fastq文件整理成一个文件列表   name R1.fastq.gz R2.fastq.gz

        Lands_chr1:10000-40000.R2.fastq.gz
        :return:
        """
        with open(self.output_dir + "/fastq_list.txt", "w") as w:
            for files in os.listdir(self.option("fastq_dir")):
                # m = re.match(r'(.*)\.R2\.fastq\.gz$', files)
                m = re.match(r'(.*)\.fastq\.gz$', files)
                if m:
                    w.write(m.group(1) + '\t' + os.path.join(self.option("fastq_dir"),
                                                             "{}.fastq.gz".format(m.group(1))) + '\n')
                    self.sample_ids.append(m.group(1).split('.')[0])
        self.logger.info(self.sample_ids)

    def set_config_file(self):
        """
        根据fastq_list.txt文件生成sample.config文件
        :return:
        """
        time.sleep(2)
        new_sample_ids = list(set(self.sample_ids))  # 去了重复
        with open(self.output_dir + "/fastq_list.txt", "r") as r:
            data = r.readlines()
            self.logger.info(data)
            if self.option('make_method') == 'combine':
                for sample_id in new_sample_ids:
                    with open("{}/{}.config".format(self.output_dir, sample_id), "w") as w:
                        n = 1
                        for line in data:
                            temp = line.strip().split("\t")
                            if temp[0].split('.')[0] == sample_id:
                                self.logger.info("qqq.{}".format(temp))
                                # w.write("max_rd_len=150" + "\n" + "[LIB]" + "\n" + "avg_ins=400" + "\n" +
                                #         "reverse_seq=0" + "\n" + "asm_flags=3" + "\n" + "rank={}\n".format(n) +
                                #         "q1={}\n".format(temp[1]) + "q2={}\n".format(temp[2]))
                                w.write("max_rd_len=150" + "\n" + "[LIB]" + "\n" + "avg_ins=400" + "\n" +
                                        "reverse_seq=0" + "\n" + "asm_flags=3" + "\n" + "rank={}\n".format(n) +
                                        "q={}\n".format(temp[1]))
                            n += 1
            else:
                for line in data:
                    temp = line.strip().split("\t")
                    nn = os.path.basename(temp[1]).split('.')
                    outfile_name = '_'.join([nn[0], nn[1]])   # Lands_chr1:10000-40000
                    with open("{}/{}.config".format(self.output_dir, outfile_name), 'w') as w:
                        w.write("max_rd_len=150" + "\n" + "[LIB]" + "\n" + "avg_ins=400" + "\n" +
                                "reverse_seq=0" + "\n" + "asm_flags=3" + "\n" + "rank={}\n".format(1) +
                                "q={}\n".format(temp[1]))

    def run_ssr_config(self):
        """
        denovo.make.config.pl
        生成配置文件
        """
        cmd = "{} {} -f {} -o {}".format(self.perl_path, self.config_path, self.option("fastq_list"), self.output_dir)
        command = self.add_command("soapdenovo_127mer_config", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("soapdenovo_127mer_config完成")
        else:
            self.set_error("soapdenovo_127mer_config失败", code="34503801")

    def run(self):
        super(MakeConfigTool, self).run()
        if self.option("fastq_dir"):
            self.get_fastq_list()
            self.set_config_file()
        else:
            self.run_ssr_config()
        self.end()
