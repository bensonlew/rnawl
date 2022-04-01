# -*- coding: utf-8 -*-
# __author__ = 'fengyitong,fuwenyao,qinjincheng'

import os
import shutil
import unittest

import pandas as pd
from Bio import SeqIO

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.files.sequence.file_sample import FileSampleFile

FAKE_FASTQ = """@ST-E00575:252:HJ7HCCCXY:7:1101:8674:1713 1:N:0:CGTACG
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


class HiseqReadsStatModule(Module):
    """
    对hiseq的PE或者SE测序数据做统计，包括GC含量，总reads数目等
    modified by fengyitong at 20180808 -- 增加了用tool运行解压缩的功能
    modified by fuwenyao at 20190506 -- 新增rfam比对；新增dup数据数据展示
    """

    def __init__(self, work_id):
        super(HiseqReadsStatModule, self).__init__(work_id)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # fastq文件夹
            {"name": "fq_type", "type": "string"},  # PE OR SE
            {"name": "dup", "type": "bool", "default": False},  # PE OR SE
            {"name": "rfam", "type": "bool", "default": False},  # 质控前统计为False，质控后统计为True 20190506 by fwy
            {"name": "quality", "type": "int", "default": 33},
            {"name": "rfam_version", "type": "string", "default": "14.1"},
        ]
        self.add_option(options)
        self.samples = {}
        self.tools = []
        self.gz2fastq = self.add_tool('medical_transcriptome.gzfastq2fastq')
        self.stat = self.add_tool('medical_transcriptome.fastq_stat')
        self.step.add_steps("dup", "stat")
        self.if_gz = 0
        self.fake_fastq = FAKE_FASTQ

    def check_options(self):
        """
        检查参数
        """
        if not self.option("fastq_dir").is_set:
            raise OptionError("需要传入fastq文件或者文件夹", code="23700601")
        files = os.listdir(self.option("fastq_dir").prop['path'])
        if 'list.txt' not in files:
            with open(os.path.join(self.option("fastq_dir").prop['path'], 'list.txt'), 'w') as list_w:
                for _file in sorted(files):
                    if u'fastq' in _file or u'fq' in _file:
                        name = _file.split('.')[0]
                        rl = 'l'
                        if u'.2.' in _file:
                            rl = 'r'
                        list_w.write(_file + '\t' + str(name) + '\t' + rl + '\n')
        if self.option("fastq_dir").is_set:
            list_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
            if not os.path.exists(list_path):
                OptionError("缺少list文件", code="23700602")
            with open(list_path, 'r') as list_r:
                list_info = list_r.read()
                if u'.gz' in list_info:
                    self.if_gz = 1
                    for line in list_info.split('\n'):
                        line = line.strip().split('\t')
                        if u'.gz' in line[0] and not os.path.exists(
                                self.option("fastq_dir").prop["path"] + "/" + ".".join(
                                    line[0].split(".")[:-2]) + ".fastq"):
                            with open(self.option("fastq_dir").prop["path"] + "/" + ".".join(
                                    line[0].split(".")[:-2]) + ".fastq", 'w') as fake:
                                fake.write(self.fake_fastq)
            self.samples = self.get_list()
            row_num = len(open(list_path, "r").readline().split())
            self.logger.info(row_num)
            if self.option('fq_type') == "PE" and row_num != 3:
                raise OptionError("PE序列list文件应该包括文件名、样本名和左右端说明三列", code="23700603")
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
                    if isinstance(self.samples[s], dict):
                        if self.samples[s]['s'].split(".")[-1] in ["gz"]:
                            self.samples[s]['s'] = ".".join(self.samples[s]['s'].split(".")[:-2]) + ".fastq"
                    else:
                        if self.samples[s].split(".")[-1] in ["gz"]:
                            self.samples[s] = ".".join(self.samples[s].split(".")[:-2]) + ".fastq"

    def get_list(self):
        list_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
        file_sample = FileSampleFile()
        file_sample.set_path(list_path)
        samples = file_sample.get_list()
        self.logger.info(samples)
        return samples

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

    def run(self):
        super(HiseqReadsStatModule, self).run()
        if self.if_gz:
            self.ungzfastq()
        else:
            self.stat_run()

    def ungzfastq(self):
        self.logger.info('需要解压fastq文件，正在解压')
        self.gz2fastq.set_options({'fastq_path': self.option('fastq_dir').prop["path"]})
        self.gz2fastq.on("end", self.stat_run)
        self.gz2fastq.run()

    def stat_run(self):
        self.stat.set_options({
            'fastq': self.option('fastq_dir').prop["path"],
            'fq_type': self.option('fq_type'),
            'quality': self.option('quality')
        })
        self.step.stat.start()
        self.stat.on("end", self.stat_finish_update)
        self.stat.run()
        self.tools.append(self.stat)
        self.run_tools()

    def run_tools(self):
        if self.option("dup") is True:
            self.dup_run()
        if self.option("rfam") is True:
            self.rfam_run()
        self.draw_run()
        self.on_rely(self.tools, self.set_output)
        self.logger.info('{}'.format(self.events))
        for eve in self.events.values():
            self.logger.info('{}'.format(eve.is_start))

    def dup_run(self):
        n = 0
        for sample_name in sorted(self.samples):
            options = {}
            if self.option("fq_type") == "PE":
                fq_l = os.path.join(self.option("fastq_dir").prop["path"], self.samples[sample_name]["l"])
                fq_r = os.path.join(self.option("fastq_dir").prop["path"], self.samples[sample_name]["r"])
                options = {'sample_name': sample_name,
                           'fastq_l': fq_l,
                           'fastq_r': fq_r,
                           'fq_type': self.option('fq_type')}
            elif self.option("fq_type") == "SE":
                fq_s = os.path.join(self.option("fastq_dir").prop["path"], self.samples[sample_name]['s'])
                options = {'sample_name': sample_name,
                           'fastq_s': fq_s,
                           'fq_type': self.option('fq_type')}
            dup = self.add_tool('medical_transcriptome.fastq_dup')
            self.step.add_steps('dup_{}'.format(n))
            dup.set_options(options)
            step = getattr(self.step, 'dup_{}'.format(n))
            step.start()
            dup.on("end", self.finish_update, "dup_{}".format(n))
            dup.run()
            self.tools.append(dup)
            n += 1

    def rfam_run(self):
        # fasta_dir = os.path.join(self.work_dir, 'FastaDir')
        # if not os.path.isdir(fasta_dir):
        #     os.mkdir(fasta_dir)
        # if not os.path.isfile(os.path.join(fasta_dir, 'done.txt')):
        #     k = 'l' if self.option('fq_type') == 'PE' else 's'
        #     for sample_name in self.samples:
        #         fq_fp = os.path.join(self.option('fastq_dir').path, self.samples[sample_name][k])
        #         fa_fp = os.path.join(fasta_dir, '{}.fa'.format(sample_name))
        #         self.fastq_to_fasta(fq_fp, fa_fp)
        #     open(os.path.join(fasta_dir, 'done.txt'), 'w').close()
        for n, sample_name in enumerate(self.samples):
            k = 'l' if self.option('fq_type') == 'PE' else 's'
            rfam = self.add_tool('ref_rna_v2.blast_to_rfam')
            self.step.add_steps('rfam_{}'.format(n))
            rfam.set_options({
                'query': os.path.join(self.option('fastq_dir').path,self.samples[sample_name][k]),
                'sample_name': sample_name,
                'rfam_version': self.option('rfam_version')
            })
            step = getattr(self.step, 'rfam_{}'.format(n))
            step.start()
            rfam.on("end", self.finish_update, "rfam_{}".format(n))
            rfam.run()
            self.tools.append(rfam)

    @staticmethod
    def fastq_to_fasta(fq_fp, fa_fp, n_lines=100000):
        sequences = list()
        n_total = 0
        for i, r in enumerate(SeqIO.parse(fq_fp, 'fastq')):
            if n_total >= n_lines:
                break
            if i % 100 == 0:
                sequences.append(r)
                n_total += 1
        SeqIO.write(sequences, fa_fp, 'fasta')

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
                if isinstance(self.samples[f], dict):
                    fq_s = os.path.join(self.option("fastq_dir").prop["path"], self.samples[f]['s'])
                else:
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

    def set_output(self):
        self.logger.info("set output")
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        draw_dir = os.path.join(self.output_dir, "qualityStat")
        os.mkdir(draw_dir)
        dup_out = []
        rfam_detail = {}
        for tool in self.tools:
            out_files = os.listdir(tool.output_dir)
            for f in out_files:
                f_path = os.path.join(tool.output_dir, f)
                if "qual_stat" in f:
                    target = os.path.join(draw_dir, f)
                    if os.path.exists(target):
                        os.remove(target)
                    os.link(f_path, target)
                elif "dup" in f:
                    dup_out.append(f_path)
                elif "rfam.stat.xls" in f:
                    rafm_details = open(f_path, "r")
                    rafm_details.readline()
                    for line in rafm_details:
                        rfam_values = line.split()
                        rfam_detail[rfam_values[0]] = rfam_values[1]
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
                    sample_name = os.path.basename(f).split(".")[0]
                    sample_name_final = sample_name
                    f = open(f, "r")
                    f.readline()
                    w.write("{}\t{}".format(sample_name_final, f.next()))
            if os.path.exists(self.output_dir + "/dup.xls"):
                os.remove(self.output_dir + "/dup.xls")
            os.link(self.work_dir + "/dup.xls", self.output_dir + "/dup.xls")
        stat_results = open(self.output_dir + "/stat_results", "w")
        if self.option("rfam") is True:
            with open(self.output_dir + "/{}".format("fastq_stat.xls"), "r") as fastq_stat_results:
                first_line = fastq_stat_results.readline()
                first_lines = first_line.split()
                stat_results.write(
                    first_lines[0] + "\t" + first_lines[1] + "\t" + first_lines[2] + "\t" + first_lines[10] + "\t" +
                    first_lines[11] + "\t" + first_lines[12] + "\t" + first_lines[13] + "\t" + "rRNA Ratio(%)" + "\n")
                for line in fastq_stat_results.readlines():
                    values = line.split()
                    stat_results.write(
                        values[0] + "\t" + values[1] + "\t" + values[2] + "\t" + values[10] + "\t" + values[11] + "\t" +
                        values[12] + "\t" + values[13] + "\t" + rfam_detail[values[0]] + "\n")
        else:
            with open(self.output_dir + "/{}".format("fastq_stat.xls"), "r") as fastq_stat_results:
                first_line = fastq_stat_results.readline()
                first_lines = first_line.split()
                stat_results.write(
                    first_lines[0] + "\t" + first_lines[1] + "\t" + first_lines[2] + "\t" + first_lines[10] + "\t" +
                    first_lines[11] + "\t" + first_lines[12] + "\t" + first_lines[13] + "\n")
                for line in fastq_stat_results.readlines():
                    values = line.split()
                    stat_results.write(
                        values[0] + "\t" + values[1] + "\t" + values[2] + "\t" + values[10] + "\t" + values[11] + "\t" +
                        values[12] + "\t" + values[13] + "\n")
        stat_results.close()
        if self.option("dup"):
            df1 = pd.read_table(self.output_dir + "/dup.xls", index_col="sample")
            df2 = pd.read_table(self.output_dir + "/stat_results", index_col="#Sample_ID")
            df3 = pd.concat([df2, df1], axis=1, join_axes=[df2.index])
            df3.to_csv(self.output_dir + "/final_results", sep='\t')
        else:
            os.link(self.output_dir + "/stat_results", self.output_dir + "/final_results")
        self.logger.info("done")
        self.end()

    @staticmethod
    def rename(event):
        obj = event["bind_object"]
        for f in os.listdir(obj.output_dir):
            old_name = os.path.join(obj.output_dir, f)
            if not event["data"] in f:
                new_name = os.path.join(obj.output_dir, event["data"] + "." + f)
                os.rename(old_name, new_name)

    def end(self):
        super(HiseqReadsStatModule, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """

    def test_before_qc(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "HiseqReadsStat_" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "medical_transcriptome.hiseq_reads_stat",
            "instant": False,
            "options": dict(
                fastq_dir="/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38312/remote_input/fastq_dir/mRNA",
                fq_type="PE",
                dup=False,
                rfam=False,
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_after_qc(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'hiseq_reads_stat_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'medical_transcriptome.hiseq_reads_stat',
            'instant': False,
            'options': {
                'fastq_dir': '/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38312/FastpRna/output/fastq',
                'fq_type': 'PE',
                'dup': True,
                'rfam': True,
                'quality': 33
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
    # suite = unittest.TestSuite()
    # suite.addTests([TestFunction('test_after_qc')])
    # unittest.TextTestRunner(verbosity=2).run(suite)
