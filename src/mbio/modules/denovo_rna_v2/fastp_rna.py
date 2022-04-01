# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict
import shutil
import unittest


class FastpRnaModule(Module):
    """
    FastpRna质控
    author: wangzhaoyue
    last_modify: 2017.12.13
    """

    def __init__(self, work_id):
        super(FastpRnaModule, self).__init__(work_id)
        options = [
            {"name": "sample_path", "type": "infile", "format": "datasplit.list_file"},
            # fastq路径list.txt文件，第一列路径，第二列样本名，第三列序列类型 l or r
            {'name': 'qualified_quality_phred', 'type': "string", "default": "20"},
            {'name': 'fq_type', 'type': "string", "default": "PE"},
            # -q,一个碱基合格的质量值,默认表示phred质量> = Q是合格的。
            {'name': 'length_required', "type": "string", "default": "30"},  # -l,长度过滤参数，比此值短的读取将被丢弃
            {'name': "cut_by_quality5", "type": "string", "default": "20"},  # -5,根据前面(5 ')的质量，允许每个读切割，默认是禁用的
            {'name': "cut_by_quality3", "type": "string", "default": "3"},  # -3,根据后面(3 ')的质量，允许每个读切割，默认是禁用的
            {'name': 'cut_mean_quality', "type": "string", "default": "20"},  # -M,在滑动窗口的基础上的平均质量低于切割质量将被切割，默认是Q20
            {'name': 'n_base_limit', "type": "string", "default": "5"},  # -n,如果reads的碱基数大于该值，那么这个reads就被丢弃了
            {'name': 'compression', "type": "string", "default": "6"},  # -z,gzip输出的压缩级别(1 ~ 9). 1是最快的，9是最小的
            {'name': 'thread', "type": "string", "default": "8"},  # -w,线程数
            {"name": "adapter_sequence", "type": "string", "default": "AGATCGGAAGAGCACACGTC"},  # --adapter_sequence,the adapter for read1
            {"name": "adapter_sequence_r2", "type": "string", "default": "AGATCGGAAGAGCGTCGTGT"},  # --adapter_sequence_r2,the adapter for read2 (PE data only)
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
        # row_num = len(open(self.option("sample_path").prop['path'], "r").readline().split())
        # if row_num != 3:
        #     raise OptionError("PE序列list文件应该包括文件名、样本名和左右端说明三列")
        # return True

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def fastp_run(self):
        n = 0
        for sample in self.sample_info:
            self.fastp = self.add_tool('datasplit.fastp_rna')
            self.step.add_steps('fastp{}'.format(n))
            opts = {
                "qualified_quality_phred": self.option('qualified_quality_phred'),
                "length_required": self.option('length_required'),
                "cut_mean_quality": self.option('cut_mean_quality'),
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

            self.fastp.set_options(opts)
            step = getattr(self.step, 'fastp{}'.format(n))
            step.start()
            self.step.update()
            self.fastp.on('end', self.finish_update, 'fastp{}'.format(n))
            self.fastp.on('end', self.set_output, 'fastp_{}'.format(sample))  # modified by zengjing 2017.12.28(目的输出文件以样本命名)
            self.tools.append(self.fastp)
            n += 1
        self.logger.info(self.tools)
        # self.on_rely(self.tools, self.set_output)
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
            # shutil.rmtree(json_dir)
            os.mkdir(json_dir)
        if not os.path.exists(fastq_dir):
            # shutil.rmtree(fastq_dir)
            os.mkdir(fastq_dir)

        for f in sorted(os.listdir(obj.output_dir)):
            # fastq_info += s + '\t'
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
            list_txt = os.path.join(self.output_dir, "fastq", "list.txt")
            fastq_txt = os.path.join(self.output_dir, "fastq", "fq_list.txt")
            self.samples=[]
            self.sample2p={}
            for i in self.fastq_info.strip().split("\n"):
                i=i.strip().split("\t")
                self.samples.append(i[0])
                sample_path="\t".join(i[1:])
                self.sample2p[i[0]] = sample_path
            with open(list_txt, 'w') as list_w, open(fastq_txt, 'w') as fastq_w:
                list_w.write(self.list_info)
                for i in sorted(self.samples):
                    fastq_w.write(i+"\t"+self.sample2p[i]+"\n")
            self.option("sickle_dir", fastq_dir)
            self.option("fq_list", fastq_txt)
            self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
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
                        raise Exception('进行fastp的sample_path：{}文件必须是三列，请检查！'.format(self.option('sample_path').prop['path']))
                    if tmp[1] in self.sample_info.keys():
                        if tmp[2] == 'l':
                            self.sample_info[tmp[1]].insert(0, tmp[0])
                        else:
                            self.sample_info[tmp[1]].append(tmp[0])
                    else:
                        self.sample_info[tmp[1]].append(tmp[0])
                else:
                    # if len(tmp) != 2:
                    #     raise Exception(
                    #         '进行fastp的sample_path：{}文件必须是两列，请检查！'.format(self.option('sample_path').prop['path']))
                    self.sample_info[tmp[1]].append(tmp[0])

            for key in self.sample_info.keys():
                if len(self.sample_info[key]) > 2:
                    raise Exception('需要质控的序列样本{}有重名，请改样本名或分开质控！'.format(key))
                # elif len(self.sample_info[key]) < 2:
                #     raise Exception('样本{}对应的R1,R2序列不全,请核实！'.format(key))

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "DenovornaQc_fastp_" + str(random.randint(1, 10000))+"yyy",
            "type": "module",
            "name": "denovo_rna_v2.fastp_rna",
            "instant": False,
            "options": dict(
                sample_path="/mnt/ilustre/users/sanger-dev/workspace/20190911/Refrna_i-sanger_203264/remote_input/fastq_dir/s2/abs.list.txt",
                # sample_path="/mnt/ilustre/users/sanger-dev/workspace/20190930/Denovorna_tsg_35690/remote_input/fastq_dir/rawdata/tmp.txt",
                fq_type="SE",
                length_required="30",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()