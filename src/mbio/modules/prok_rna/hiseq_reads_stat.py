#!/usr/bin/env python
# -*- coding: utf-8 -*-
# last modified by shicaiping at 20180507
import os
import shutil
import glob
from Bio import SeqIO
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.files.sequence.file_sample import FileSampleFile
import unittest
from biocluster.config import Config
#import subprocess

class HiseqReadsStatModule(Module):
    """
    对hiseq的PE或者SE测序数据做统计，包括GC含量，总reads数目等
    """
    def __init__(self, work_id):
        super(HiseqReadsStatModule, self).__init__(work_id)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # fastq文件夹
            {"name": "fq_type", "type": "string"},  # PE OR SE
            {"name": "dup", "type": "bool", "default": False},  # 质控前统计为False，质控后统计为True
            {"name": "rfam", "type": "bool", "default": False},  # 质控前统计为False，质控后统计为True
            {"name": "rrna", "type": "string", "default": None},  # 质控前统计为False，质控后统计为True
            {"name": "quality", "type": "int", "default": 33},
            {"name": "strand_specific", "type": "bool", "default": False}  # 链特异性参数
        ]
        self.add_option(options)
        self.samples = {}
        self.tools = []
        self.stat = self.add_tool('prok_rna.fastq_stat')
        self.step.add_steps("dup", "stat", "rfam")
        self.gz2fastq = self.add_tool('prok_rna.gzfastq2fastq')
        self.if_gz = 0
        self.fake_fastq = """@ST-E00575:252:HJ7HCCCXY:7:1101:8674:1713 1:N:0:CGTACG
GNCCCAACTTGCCATCAAGGATATCTATCTCGGCAACCGCTTCGTTAAATGTCTCTTCGTGGTCAGCTTCAATAGCCAATTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTTGAAA
+
A#AF-7FFJFJJAJJJJJJJJJJJJJJ-AJFJJJJJJJJJJJJJ<JJJJAJJ7FJJFFFJFFFAFJJJJJAJJJJFFJJJJJJFJFFJJFJJJJFJJJ<FJJJJFJFJJJJJFJJJJJJJJJJJJJJ7FJJJJJJJJ7AJJJJJJJFJJJJ
@ST-E00575:252:HJ7HCCCXY:7:1101:9810:1713 1:N:0:CGTACG
CNGTAATCGTTTGTGGCGTTAGAAATAAAGCCTCAGCCGCCCCGACGACAGAGCCTTCCTTACAAACTTGCCAAAAATAATAAAGATGATTGAAATTGATGTGCGACATTCGCATGTTGTTATCCCCAGATCGGAAGAGCCACACGTCTGA
+
A#AFFJJJJJJJJJJJJJJJJJFJJJJ<JJJJJJJJJJJJAJJJJJJJJJJJJJJJ<-JJJJJJJJJJJJAAJJJJJJJJJJJJJJJJFFJJJJJJJJFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJA
@ST-E00575:252:HJ7HCCCXY:7:1101:12408:1713 1:N:0:CGTACG
CNCACCAATCATCCTGGACTGGCTCTCAATCTCCATCCTGGAGGTGTCCTTTGTTTCTTCCTGAAACATCCCTTCACTCATCCTAAGCAGTCCCTGAGTCCTTCATCCTGAAGTGGCACCATCCTGATACCGTCCTTTAGATCGGAAGAGC
+
A#AFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJFFJFJJJJFJJJJJFFJJJJJ<AJJJJJJJJJJJJJJJJJJJFJFJJJ<AFFJFJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJF
"""

    def check_options(self):
        """
        检查参数
        """
        if not self.option("fastq_dir").is_set:
            raise OptionError("需要传入fastq文件或者文件夹", code = "25000601")
        # if self.option("fastq_dir"):
        #     self.option("fastq_dir").get_full_info(self.option("fastq_dir").prop["path"])
        if self.option("fastq_dir").is_set:
            list_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
            if not os.path.exists(list_path):
                OptionError("缺少list文件", code = "25000602")
            with open(list_path, 'r') as list_r:
                list_info = list_r.read()
                if u'.gz' in list_info:
                    self.if_gz = 1
                    for line in list_info.split('\n'):
                        line = line.strip().split('\t')
                        if u'.gz' in line[0] and not os.path.exists(self.option("fastq_dir").prop["path"] + "/" + ".".join(line[0].split(".")[:-2]) + ".fastq"):
                            with open(self.option("fastq_dir").prop["path"] + "/" + ".".join(line[0].split(".")[:-2]) + ".fastq" ,'w') as fake:
                                fake.write(self.fake_fastq)
            self.samples = self.get_list()
            row_num = len(open(list_path, "r").readline().split())
            self.logger.info(row_num)
            if self.option('fq_type') == "PE" and row_num != 3:
                raise OptionError("PE序列list文件应该包括文件名、样本名和左右端说明三列", code = "25000603")
            elif self.option('fq_type') == "SE" and row_num != 2:
                raise OptionError("SE序列list文件应该包括文件名、样本名两列", code = "25000604")
            if self.option('fq_type') == "PE":
                for s in sorted(self.samples):
                    if 'l' not in self.samples[s]:
                        self.set_error("list文件中样本{}没有左端测序数据".format(s))
                    elif 'r' not in self.samples[s]:
                        self.set_error("list文件中样本{}没有右端测序数据".format(s))
                    if self.samples[s]["l"].split(".")[-1] in ["gz"]:
                        self.samples[s]["l"] = ".".join(self.samples[s]["l"].split(".")[:-2]) + ".fastq"
                    if self.samples[s]["r"].split(".")[-1] in ["gz"]:
                        self.samples[s]["r"] = ".".join(self.samples[s]["r"].split(".")[:-2]) + ".fastq"
            if self.option('fq_type') == "SE":
                for s in sorted(self.samples):
                    if self.samples[s].split(".")[-1] in ["gz"]:
                        self.samples[s] = ".".join(self.samples[s].split(".")[:-2]) + ".fastq"
            print(self.samples)

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def stat_finish_update(self):
        self.step.stat.finish()
        self.step.update()

    def dup_finish_update(self):
        self.step.dup.finish()
        self.step.update()

    def stat_run(self):
        self.stat.set_options({
            'fastq': self.option('fastq_dir').prop["path"],
            'fq_type': self.option('fq_type'),
            'quality': self.option('quality')
        })
        # self.on_rely(estimators, self.rarefaction_run)
        self.step.stat.start()
        self.stat.on("end", self.stat_finish_update)
        self.stat.on("end", self.stat_finish_update)
        self.stat.run()
        self.tools.append(self.stat)
        self.run_tools()

    def dup_run(self):
        n = 0
        for f in sorted(self.samples):
            options = {}
            if self.option("fq_type") == "PE":
                fq_l = os.path.join(self.option("fastq_dir").prop["path"], self.samples[f]["l"])
                fq_r = os.path.join(self.option("fastq_dir").prop["path"], self.samples[f]["r"])
                options = {
                    'fastq_l': fq_l,
                    'fastq_r': fq_r,
                    'fq_type': self.option('fq_type')
                }
            elif self.option("fq_type") == "SE":
                fq_s = os.path.join(self.option("fastq_dir").prop["path"], self.samples[f])
                options = {
                    'fastq_s': fq_s,
                    'fq_type': self.option('fq_type')
                }
            dup = self.add_tool('sequence.fastq_dup')
            self.step.add_steps('dup_{}'.format(n))
            dup.set_options(options)
            step = getattr(self.step, 'dup_{}'.format(n))
            step.start()
            dup.on("end", self.finish_update, "dup_{}".format(n))
            dup.on("end", self.rename, f)
            dup.run()
            self.tools.append(dup)
            n += 1

    def rfam_run(self):
        n = 0
        for f in sorted(self.samples):
            options = {}
            if self.option("fq_type") == "PE":
                fq_l = os.path.join(self.option("fastq_dir").prop["path"], self.samples[f]["l"])
                # small_fq = os.path.join(self.work_dir, f + ".fq")
                # small_fa = os.path.join(self.work_dir, f + ".fa")
                # self.fastq_convert_fasta(fq_l, small_fa)
                # with open (fq_l, "r") as f, open (small_fq, "w") as w:
                #     line_nu = 0
                #     for line in f:
                #         if line_nu < 400000:
                #             w.writelines(line)
                #             line_nu += 1
                #         else:
                #             break
                # self.fastq_convert_fasta(small_fq, small_fa)
                options = {
                    'query_fq': fq_l,
                    'sample_name': f
                }
            elif self.option("fq_type") == "SE":
                fq_s = os.path.join(self.option("fastq_dir").prop["path"], self.samples[f])
                # small_fq = os.path.join(self.work_dir, f + ".fq")
                # small_fa = os.path.join(self.work_dir, f + ".fa")
                # self.fastq_convert_fasta(fq_l, small_fa)
                # with open (fq_s, "r") as f, open (small_fq, "w") as w:
                #     line_nu = 0
                #     for line in f:
                #         if line_nu < 400000:
                #             w.writelines(line)
                #             line_nu += 1
                #         else:
                #             break
                # self.fastq_convert_fasta(small_fq, small_fa)
                options = {
                    'query_fq': fq_s,
                    'sample_name': f
                }
            rfam = self.add_tool('prok_rna.blast_to_rfam')
            self.step.add_steps('rfam_{}'.format(n))
            rfam.set_options(options)
            step = getattr(self.step, 'rfam_{}'.format(n))
            step.start()
            rfam.on("end", self.finish_update, "rfam_{}".format(n))
            #rfam.on("end", self.rename, f)
            rfam.run()
            self.tools.append(rfam)
            n += 1

    def rrna_run(self):
        n = 0
        for f in sorted(self.samples):
            options = {}
            mapping_tool = self.add_tool('prok_rna.bowtie2')
            if self.option("fq_type") == "PE":
                fq_l = os.path.join(self.option("fastq_dir").prop["path"], self.samples[f]["l"])
                fq_r = os.path.join(self.option("fastq_dir").prop["path"], self.samples[f]["r"])
                options = {
                    "ref_genome": self.option("rrna"),
                    'left_reads': fq_l,
                    'right_reads': fq_r,
                    'sample': f,
                    "seq_method": self.option("fq_type")
                }
            elif self.option("fq_type") == "SE":
                fq_s = os.path.join(self.option("fastq_dir").prop["path"], self.samples[f])
                options = {
                    "ref_genome": self.option("rrna"),
                    'single_end_reads': fq_s,
                    'sample': f,
                    "seq_method": self.option("fq_type")
                }

            self.logger.info("tool {} prok_rna.bowtie2".format(options))
            mapping_tool.set_options(options)
            self.tools.append(mapping_tool)
            # mapping_tool.on("end", self.finish_update, "mapping_{}".format(n))
            mapping_tool.run()
            n += 1

    def fastq_convert_fasta(self, input, output, n_lines=100000):
        """
        将fastq转化成fasta
        """
        sequences = list()
        n_total = 0
        for i, r in enumerate(SeqIO.parse(input, 'fastq')):
            if n_total >= n_lines:
                break
            if i % 100 == 0:
                sequences.append(r)
                n_total += 1
        SeqIO.write(sequences, output, 'fasta')
        # self.fastq_to_fasta_path = os.path.join(Config().SOFTWARE_DIR, "bioinfo/seq/fastx_toolkit_0.0.14/fastq_to_fasta")
        # try:
        #     convert_str = (self.fastq_to_fasta_path + ' -n -i '
        #                    + input + ' -o ' + output)
        #     subprocess.check_call(convert_str, shell=True)
        #     self.is_convert = True
        # except subprocess.CalledProcessError:
        #     self.set_error("fastq转化fasta失败", code = "25000605")

    def draw_run(self):
        n = 1
        if self.option("fq_type") == "PE":
            for f in sorted(self.samples):
                fq_l = os.path.join(self.option("fastq_dir").prop["path"], self.samples[f]["l"])
                fq_r = os.path.join(self.option("fastq_dir").prop["path"], self.samples[f]["r"])
                draw_l = self.add_tool('sequence.draw_fastq_info')
                draw_r = self.add_tool('sequence.draw_fastq_info')
                self.step.add_steps('drawL_{}'.format(n))
                self.step.add_steps('drawR_{}'.format(n))
                draw_l.set_options({
                    "fastq": fq_l
                })
                draw_r.set_options({
                    "fastq": fq_r
                })
                step = getattr(self.step, 'drawL_{}'.format(n))
                step.start()
                step = getattr(self.step, 'drawR_{}'.format(n))
                step.start()
                draw_l.on("end", self.finish_update, "drawL_{}".format(n))
                draw_r.on("end", self.finish_update, "drawR_{}".format(n))
                draw_l.on("end", self.rename, "{}.l".format(f))
                draw_r.on("end", self.rename, "{}.r".format(f))
                draw_l.run()
                draw_r.run()
                self.tools.append(draw_l)
                self.tools.append(draw_r)
                n += 1
        else:
            for f in sorted(self.samples):
                fq_s = os.path.join(self.option("fastq_dir").prop["path"], self.samples[f])
                draw = self.add_tool('sequence.draw_fastq_info')
                self.step.add_steps('draw_{}'.format(n))
                draw.set_options({
                    "fastq": fq_s,
                })
                step = getattr(self.step, 'draw_{}'.format(n))
                step.start()
                draw.on("end", self.finish_update, "draw_{}".format(n))
                draw.on("end", self.rename, f)
                draw.run()
                self.tools.append(draw)
                n += 1

    def ungzfastq(self):
        self.logger.info('需要解压fastq文件，正在解压')
        self.gz2fastq.set_options({'fastq_path': self.option('fastq_dir').prop["path"]})
        self.gz2fastq.on("end", self.stat_run)
        self.gz2fastq.run()
        self.logger.info('解压完成')

    def run_tools(self):
        if self.option("dup") is True:
            self.dup_run()
        if self.option("rfam") is True:
            self.rfam_run()
        if self.option("rrna"):
            self.rrna_run()
        self.draw_run()
        self.on_rely(self.tools, self.set_output)
        self.logger.info('{}'.format(self.events))
        for eve in self.events.values():
            self.logger.info('{}'.format(eve.is_start))

    def run(self):
        super(HiseqReadsStatModule, self).run()
        if self.if_gz:
            self.ungzfastq()
        else:
            self.stat_run()

    def set_output(self):
        self.logger.info("set output")
        for f in glob.glob(r"{}/*".format(self.output_dir)):
            if os.path.isdir(f):
                shutil.rmtree(f)
            else:
                os.remove(f)
        draw_dir = os.path.join(self.output_dir, "qualityStat")
        rfam_dir = os.path.join(self.output_dir, "RfamStat")
        os.mkdir(draw_dir)
        os.mkdir(rfam_dir)
        rrna_dir = os.path.join(self.output_dir, "RrnaStat")
        if self.option("rrna"):
            os.mkdir(rrna_dir)
        dup_out = []
        for tool in self.tools:
            out_files = os.listdir(tool.output_dir)
            for f in out_files:
                f_path = os.path.join(tool.output_dir, f)
                if "qual_stat" in f:
                    target = os.path.join(draw_dir, f)
                    if os.path.exists(target):
                        os.remove(target)
                    os.link(f_path, target)
                elif "rfam" in f:
                    target = os.path.join(rfam_dir, f)
                    if os.path.exists(target):
                        os.remove(target)
                    os.link(f_path, target)
                elif ".stat" in f:
                    target = os.path.join(rrna_dir, f)
                    if os.path.exists(target):
                        os.remove(target)
                    os.link(f_path, target)
                elif "dup" in f:
                    dup_out.append(f_path)
                else:
                    target = self.output_dir + "/{}".format("fastq_stat.xls")
                    if os.path.exists(target):
                        os.remove(target)
                    os.link(f_path, target)
        if self.option("dup") is True:
            with open(self.work_dir + "/dup.xls", "w") as w:
                if self.option("fq_type") == "PE":
                    w.write("sample\tread1Dup\tread2Dup\tPairedDup\n")
                else:
                    w.write("sample\treadDup\n")
                for f in dup_out:
                    if '.fastq_dup.xls' in os.path.basename(f):
                        sample_name = os.path.basename(f).split(".fastq_dup.xls")[0]
                    else:
                        sample_name = os.path.basename(f).split(".")[0]
                    f = open(f, "r")
                    f.readline()
                    w.write("{}\t{}".format(sample_name, f.next()))
            if os.path.exists(self.output_dir + "/dup.xls"):
                os.remove(self.output_dir + "/dup.xls")
            os.link(self.work_dir + "/dup.xls", self.output_dir + "/dup.xls")
        self.logger.info("done")
        self.end()

    def rename(self, event):
        obj = event["bind_object"]
        for f in os.listdir(obj.output_dir):
            old_name = os.path.join(obj.output_dir, f)
            if not event["data"] in f:
                new_name = os.path.join(obj.output_dir, event["data"] + "." + f)
                os.rename(old_name, new_name)

    def get_list(self):
        list_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
        file_sample = FileSampleFile()
        file_sample.set_path(list_path)
        samples = file_sample.get_list()
        self.logger.info(samples)
        return samples

    def end(self):
        self.logger.info('%s' % self.upload_dir)
        if self.upload_dir:
            for i in self.upload_dir:
                self.logger.info('%s' % i._parent)
                self.logger.info('%s' % i.file_list)
        result_dir = self.add_upload_dir(self.output_dir)
        self.logger.info('%s' % self.upload_dir)
        result_dir.add_relpath_rules([
                [".", "", "结果输出目录"],
                ["./qualityStat/", "文件夹", "质量统计文件夹"],
                ["./fastq_stat.xls", "xls", "fastq信息统计表"]
            ])
        if self.option("dup") is True:
            result_dir.add_relpath_rules([
                ["./dup.xls", "xls", "fastq序列重复信息"]
            ])
        super(HiseqReadsStatModule, self).end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "HiseqReadsStat_" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "prok_rna.hiseq_reads_stat",
            "instant": False,
            "options": dict(
                fastq_dir="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/data",
                # fastq_dir="/mnt/ilustre/users/sanger-dev/workspace/20180806/Single_HiseqQc_3546/HiseqQc/output/sickle_dir",
                fq_type="PE",
                dup=True,
                rfam=True
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
