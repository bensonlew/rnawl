# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue & guhaidong'

import os, sys
# import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool


class SoapDenovoAgent(Agent):
    """
    宏基因SOAPdenovo2组装
    version: SOAPdenovo2-src-r240
    author: wangzhaoyue & guhaidong
    last_modify: 2017.09.04
    """

    def __init__(self, parent):
        super(SoapDenovoAgent, self).__init__(parent)
        options = [
            {"name": "fastq1", "type": "infile", "format": "sequence.fastq"},  # 输入文件,sample.sickle.l.fastq
            {"name": "fastq2", "type": "infile", "format": "sequence.fastq"},  # 输入文件,sample.sickle.r.fastq
            {"name": "fastqs", "type": "infile", "format": "sequence.fastq"},  # 输入文件,sample.sickle.s.fastq
            {"name": "sample_name", "type": "string"},  # 样品名称
            {"name": "max_rd_len", "type": "string"},  # read最大读长
            {"name": "mem", "type": "int", "default": 100},  # 拼接使用内存
            {"name": "insert_size", "type": "string"},  # 平均插入片段长度
            {"name": "reverse_seq", "type": "string", "default": "0"},  # 配置文件的其他参数
            {"name": "asm_flags", "type": "string", "default": "3"},  # 配置文件的其他参数
            {"name": "rank", "type": "string", "default": "1"},  # 配置文件的其他参数
            {"name": "kmer", "type": "string"},  # k_mer值，例"39"
            {"name": "scafSeq", "type": "outfile", "format": "sequence.fasta"},  # 输出文件,sample.scafSeq
        ]
        self.add_option(options)
        self.step.add_steps("SOAPdenovo")
        self._memory_increase_step = 50  # 每次重运行增加内存50G by guhaidong @ 20180428
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.SOAPdenovo.start()
        self.step.update()

    def stepfinish(self):
        self.step.SOAPdenovo.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('fastq1'):
            raise OptionError('必须输入*l.fastq文件', code="31301701")
        if not self.option('fastq2'):
            raise OptionError('必须输入*r.fastq文件', code="31301702")
        if not self.option('sample_name'):
            raise OptionError('必须输入样品名称', code="31301703")
        if not self.option('max_rd_len'):
            raise OptionError('必须输入read的最大长度', code="31301704")
        if not self.option('insert_size'):
            raise OptionError('必须输入平均插入片段的长度', code="31301705")
        if not self.option('kmer'):
            raise OptionError('必须输入kmer值', code="31301706")
        if self.option('mem') > 250:
            raise OptionError('内存设置不可超过250', code="31301707")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 16
        # if self.option('mem') < 2:
        #     tmp_mem = '2'
        # else:
        #     tmp_mem = int(self.option('mem')) + 50 * self._rerun_time  # 每次因拼接失败而重运行的内存增加50G by GHD @ 20180320
        # self._memory = '%sG' % tmp_mem
        self._memory = '%sG' % self.option('mem')  # 改回 by GHD @ 20180428
        # self._memory = mem_str + "G"
        # self.logger.info('soapdenovo use memory : ' + self._memory)
        # self._memory = "100G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(SoapDenovoAgent, self).end()


class SoapDenovoTool(Tool):
    def __init__(self, config):
        super(SoapDenovoTool, self).__init__(config)
        self.sample_name = self.option('sample_name')
        self._version = "SOAPdenovo2-src-r240"
        self.SOAPdenovo_path = '/bioinfo/metaGenomic/SOAPdenovo2/bin/'

    def run(self):
        """
        运行
        :return:
        """
        super(SoapDenovoTool, self).run()
        self.init_config()
        have_result = self.run_SOAPdenovo2()  # 如果已有拼接结果，则stat为1，如果没有，则stat为0
        self.set_output(have_result)
        self.end()

    def init_config(self):
        # sample_name = os.path.basename(self.option('fastq1').prop['path']).split('.sickle.l.fastq')[0]
        config_file = self.work_dir + "/" + self.sample_name + ".config"
        with open(config_file, "w+") as fw:
            first = "max_rd_len=" + self.option('max_rd_len') + "\n"
            second = "[LIB]" + "\n"
            third = "avg_ins=" + self.option('insert_size') + "\n"
            forth = "reverse_seq=" + self.option('reverse_seq') + "\n"
            fifth = "asm_flags=" + self.option('asm_flags') + "\n"
            sixth = "rank=" + self.option('rank') + "\n"
            q1 = "q1=" + self.option('fastq1').prop['path'] + "\n"
            q2 = "q2=" + self.option('fastq2').prop['path']
            if not self.option('fastqs').is_set:
                fw.write(first + second + third + forth + fifth + sixth + q1 + q2)
            else:
                qs = "\nq=" + self.option('fastqs').prop['path']
                fw.write(first + second + third + forth + fifth + sixth + q1 + q2 + qs)

    def run_SOAPdenovo2(self):
        # sample_name = os.path.basename(self.option('fastq1').prop['path']).split('.sickle.l.fastq')[0]
        # output_dir = self.work_dir + "/" + self.sample_name + "_K" + self.option('kmer')
        # if not os.path.exists(output_dir):
        #    os.makedirs(output_dir)
        cmd = self.SOAPdenovo_path + 'SOAPdenovo2-63mer all -s %s -o %s -K %s -p 16 -d 1 -D 1 -F -u' % \
                                     (self.work_dir + "/" + self.sample_name + ".config",
                                      self.work_dir + "/" + self.sample_name + '.kmer' + self.option('kmer'),
                                      self.option('kmer'))
        if os.path.exists(self.output_dir + "/" + self.sample_name + '.kmer' + self.option('kmer') + '.scafSeq'):
            self.logger.info("%s.kmer%s.scafSeq已存在，跳过拼接" % (self.sample_name, self.option('kmer')))
            result_stat = 1
        else:
            command = self.add_command("soapdenovo", cmd, ignore_error=True)
            command.run()
            self.wait(command)
            exitcode_list = [-6,137]
            if command.return_code == 0:
                self.logger.info("运行cmd完成")
            elif command.return_code in exitcode_list:
                self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
            else:
                self.set_error("运行cmd运行出错!")
            result_stat = 0
        return result_stat

    def set_output(self, status):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        # self.sample_name = os.path.basename(self.option('fastq1').prop['path']).split('.sickle.l.fastq')[0]
        self.option('scafSeq').set_path(
            self.work_dir + "/" + self.sample_name + '.kmer' + self.option('kmer') + '.scafSeq')
        if status == 0:
            self.logger.info("现在复制结果")
            if os.path.exists(self.output_dir + "/" + self.sample_name + '.kmer' + self.option('kmer') + '.scafSeq'):
                os.remove(self.output_dir + "/" + self.sample_name + '.kmer' + self.option('kmer') + '.scafSeq')
            # shutil.copy2(self.work_dir + "/" + sample_name + '_K' + self.option('kmer') + '.scafSeq', self.output_dir +
            #         "/" + sample_name + '_K' + self.option('kmer') + '.scafSeq')
            os.link(self.work_dir + "/" + self.sample_name + '.kmer' + self.option('kmer') + '.scafSeq',
                    self.output_dir +
                    "/" + self.sample_name + '.kmer' + self.option('kmer') + '.scafSeq')
        self.logger.info("设置SOAPdenovo2分析结果目录成功")
