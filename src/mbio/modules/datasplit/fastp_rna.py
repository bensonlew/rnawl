# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

import os
import re
import time
import unittest
from collections import defaultdict

from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class FastpRnaModule(Module):
    """
    FastpRna质控
    author: wangzhaoyue
    last_modify: 2017.12.13
    """

    def __init__(self, work_id):
        super(FastpRnaModule, self).__init__(work_id)
        options = [
            # fastq路径list.txt文件，第一列路径，第二列样本名，第三列序列类型 l or r
            {"name": "sample_path", "type": "infile", "format": "datasplit.list_file"},
            {'name': 'qualified_quality_phred', 'type': "string", "default": "20"},
            # -q,一个碱基合格的质量值,默认表示phred质量> = Q是合格的。
            {'name': 'fq_type', 'type': "string", "default": "PE"},
            # -l,长度过滤参数，比此值短的读取将被丢弃
            {'name': 'length_required', "type": "string", "default": "30"},
            # phred+64 OR phred+33 #add by 20200401 fwy
            {"name": "quality_score_system", "type": "string", "default": "phred+33"},
            # -5,根据前面(5 ')的质量，允许每个读切割，默认是禁用的
            {'name': "cut_by_quality5", "type": "string", "default": "20"},
            # -3,根据后面(3 ')的质量，允许每个读切割，默认是禁用的
            {'name': "cut_by_quality3", "type": "string", "default": "3"},
            # -M,在滑动窗口的基础上的平均质量低于切割质量将被切割，默认是Q20
            {'name': 'cut_mean_quality', "type": "string", "default": "20"},
            # -n,如果reads的碱基数大于该值，那么这个reads就被丢弃了
            {'name': 'n_base_limit', "type": "string", "default": "5"},
            # -z,gzip输出的压缩级别(1 ~ 9). 1是最快的，9是最小的
            {'name': 'compression', "type": "string", "default": "6"},
            # -w,线程数
            {'name': 'thread', "type": "string", "default": "8"},
            # --adapter_sequence,the adapter for read1
            {"name": "adapter_sequence", "type": "string", "default": "AGATCGGAAGAGCACACGTC"},
            # --adapter_sequence_r2,the adapter for read2 (PE data only)
            {"name": "adapter_sequence_r2", "type": "string", "default": "AGATCGGAAGAGCGTCGTGT"},
            {"name": "adapter_sequence_s", "type": "string", "default": "AGATCGGAAGAGCACACGTC"},
            {"name": "sickle_dir", "type": "outfile", "format": "sequence.fastq_dir"},
            {"name": "fq_list", "type": "outfile", "format": "denovo_rna_v2.fastq_list"},
        ]
        self.sample_info = defaultdict(list)
        self.tools = []
        self.end_times = 0
        self.add_option(options)
        self.list_info = self.fastq_info = ''

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('sample_path'):
            raise OptionError('必须输入文库文件夹对应的路径信息')

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def fastp_run(self):
        n = 0
        for sample in sorted(self.sample_info):
            fastp = self.add_tool('ref_rna_v3.fastp_rna')
            self.step.add_steps('fastp{}'.format(n))
            opts = {
                "qualified_quality_phred": self.option('qualified_quality_phred'),
                "length_required": self.option('length_required'),
                "cut_mean_quality": self.option('cut_mean_quality'),
                "quality_score_system": self.option("quality_score_system"),
                "n_base_limit": self.option('n_base_limit'),
                "compression": self.option('compression'),
                "thread": self.option('thread'),
                'fq_type': self.option('fq_type'),
            }
            if self.option('fq_type') == "PE":
                opts.update({'fq1': self.sample_info[sample][0]})
                opts.update({'fq2': self.sample_info[sample][1]})
                opts.update({'adapter_sequence': self.option("adapter_sequence")})
                opts.update({'adapter_sequence_r2': self.option("adapter_sequence_r2")})
            else:
                opts.update({'fqs': self.sample_info[sample][0]})
                opts.update({'adapter_sequence_s': self.option("adapter_sequence_s")})
            if self.option('cut_by_quality5'):
                opts.update({'cut_by_quality5': self.option('cut_by_quality5')})
            if self.option('cut_by_quality3'):
                opts.update({'cut_by_quality3': self.option('cut_by_quality3')})
            fastp.set_options(opts)
            step = getattr(self.step, 'fastp{}'.format(n))
            step.start()
            self.step.update()
            fastp.on('end', self.finish_update, 'fastp{}'.format(n))
            # modified by zengjing 2017.12.28(目的输出文件以样本命名)
            fastp.on('end', self.set_output, 'fastp_{}'.format(sample))
            self.tools.append(fastp)
            n += 1
        self.logger.info(self.tools)
        self.step.update()
        for tool in self.tools:
            tool.run()

    def run(self):
        """
        运行
        :return:
        """
        super(FastpRnaModule, self).run()
        self.get_info()
        time.sleep(2)
        self.fastp_run()

    def set_output(self, event):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        obj = event["bind_object"]
        self.end_times += 1
        m = re.match(r"fastp_(.+)", event["data"])
        if m:
            s = m.group(1)
            self.fastq_info += s + '\t'
        else:
            raise Exception("fastp没匹配到样本名")
        json_dir = os.path.join(self.output_dir, "json_dir")
        fastq_dir = os.path.join(self.output_dir, "fastq")
        if not os.path.exists(json_dir):
            os.mkdir(json_dir)
        if not os.path.exists(fastq_dir):
            os.mkdir(fastq_dir)

        for f in sorted(os.listdir(obj.output_dir)):
            f1 = os.path.join(obj.output_dir, f)
            m = re.match(r".+(.clean.1.fastq)$", f)
            n = re.match(r".+(.clean.2.fastq)$", f)
            k = re.match(r".+(.clean.s.fastq)$", f)
            if m:
                f2 = os.path.join(fastq_dir, s + m.group(1))
                self.list_info += s + m.group(1) + '\t' + s + '\t' + 'l\n'
                self.fastq_info += f2 + '\t'
            elif n:
                f2 = os.path.join(fastq_dir, s + n.group(1))
                self.list_info += s + n.group(1) + '\t' + s + '\t' + 'r\n'
                self.fastq_info += f2 + '\n'
            elif k:
                f2 = os.path.join(fastq_dir, s + k.group(1))
                self.list_info += s + k.group(1) + '\t' + s + '\t' + 's\n'
                self.fastq_info += f2 + '\n'
            else:
                f2 = os.path.join(json_dir, s + ".json")
            if os.path.exists(f2):
                os.remove(f2)
            os.link(f1, f2)
        if len(self.tools) == self.end_times:
            # 重写list.txt和fq_list.txt，避免因重运行导致每次生成的文件中样本顺序不一致而导致的snp分析合并报错
            list_info = ''
            fastq_info = ''
            list_txt = os.path.join(self.output_dir, "fastq", "list.txt")
            fastq_txt = os.path.join(self.output_dir, "fastq", "fq_list.txt")
            if os.path.exists(list_txt):
                os.remove(list_txt)
            if os.path.exists(fastq_txt):
                os.remove(fastq_txt)
            for sample in sorted(self.sample_info):
                if self.option('fq_type') == "PE":
                    list_info += sample + '.clean.1.fastq' + '\t' + sample + '\t' + 'l\n'
                    list_info += sample + '.clean.2.fastq' + '\t' + sample + '\t' + 'r\n'
                    fastq_info += sample + '\t' + os.path.join(self.output_dir, 'fastq',
                                              sample + '.clean.1.fastq') + '\t' + os.path.join(self.output_dir, 'fastq',
                                              sample + '.clean.2.fastq') + '\n'
                else:
                    list_info += sample + '.clean.s.fastq' + '\t' + sample + '\t' + 's\n'
                    fastq_info += sample + '\t' + os.path.join(self.output_dir, 'fastq',
                                                               sample + '.clean.s.fastq') + '\n'
            with open(list_txt, 'w') as list_w, open(fastq_txt, 'w') as fastq_w:
                list_w.write(list_info)
                fastq_w.write(fastq_info)
            self.option("sickle_dir", fastq_dir)
            self.option("fq_list", fastq_txt)
            self.end()

    def end(self):
        super(FastpRnaModule, self).end()

    def get_info(self):
        """
        :return:
        """
        with open(self.option('sample_path').prop['path'])as fr:
            for line in fr:
                tmp = line.strip().split('\t')
                if self.option('fq_type') == "PE":
                    if len(tmp) != 3:
                        raise Exception(
                            '进行fastp的sample_path：{}文件必须是三列，请检查！'.format(self.option('sample_path').path))
                    if tmp[1] in self.sample_info.keys():
                        if tmp[2] == 'l':
                            self.sample_info[tmp[1]].insert(0, tmp[0])
                        else:
                            self.sample_info[tmp[1]].append(tmp[0])
                    else:
                        self.sample_info[tmp[1]].append(tmp[0])
                else:
                    self.sample_info[tmp[1]].append(tmp[0])
            for key in self.sample_info.keys():
                if len(self.sample_info[key]) > 2:
                    raise Exception('需要质控的序列样本{}有重名，请改样本名或分开质控！'.format(key))


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "RnaQc_fastp_33_pe_sad2_" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "datasplit.fastp_rna",
            "instant": False,
            "options": dict(
                    quality_score_system="phred+33",
                    # adapter_sequence_s="AGTTCCCGGAAGAGCACACGTC",
                    adapter_sequence="AGTTCCCGGAAGAGCACACGTC",
                    adapter_sequence_r2="AAGTTCCCGGAAGAGCACACTGT",
                    fq_type="PE",
                    # sample_path="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/lncRNA/test_data_new/new_test_rawdata1/Human_10smaples/rawdata/list.txt"
                    sample_path="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/STR/data/abs.list.txt",
                )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
