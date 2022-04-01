#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.files.sequence.file_sample import FileSampleFile
import re
import glob
import unittest

class HiseqQcModule(Module):
    """
    hiseq数据指控模块，主要调用seqprep、sickle软件做质量剪切与去接头
    version 1.0
    author: shicaiping
    last_modify: 20180507
    """
    def __init__(self, work_id):
        super(HiseqQcModule, self).__init__(work_id)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # fastq文件夹
            {"name": "fq_type", "type": "string"},  # PE OR SE
            {"name": "quality_score_system", "type": "string", "default": "phred+33"},  # phred+64 OR phred+33
            {"name": "clip_dir", "type": "outfile", "format": "sequence.fastq_dir"},  # SE去接头输出结果文件夹
            {"name": "sickle_dir", "type": "outfile", "format": "sequence.fastq_dir"},  # 质量剪切输出结果文件夹(包括左右段)
            {"name": "sickle_r_dir", "type": "outfile", "format": "sequence.fastq_dir"},  # 质量剪切右端输出结果文件夹
            {"name": "sickle_l_dir", "type": "outfile", "format": "sequence.fastq_dir"},  # 质量剪切左端输出结果文件夹
            {"name": "seqprep_dir", "type": "outfile", "format": "sequence.fastq_dir"},  # PE的去接头输出结果文件
            {"name": "quality_q", "type": "int", "default": 20},  # 质量剪切碱基质量
            {"name": "length_q", "type": "int", "default": 30},  # 质量剪切碱基长度
            {"name": "fq_list", "type": "outfile", "format": "prok_rna.fastq_list"},
        ]
        self.add_option(options)
        self.samples = {}
        self.seqprep = []
        self.clipper = []
        self.fastp = []
        self.sickle = []
        self.end_times = 0
        self.adapt_rate = []

    def check_options(self):
        """
        检查参数
        """
        if not self.option("fastq_dir").is_set:
            raise OptionError("需要传入fastq文件或者文件夹", code = "25000501")
        if self.option("fastq_dir").is_set:
            list_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
            if not os.path.exists(list_path):
                OptionError("缺少list文件", code = "25000502")
            row_num = len(open(list_path, "r").readline().split())
            if self.option('fq_type') == "PE" and row_num != 3:
                raise OptionError("PE序列list文件应该包括文件名、样本名和左右端说明三列", code = "25000503")
            elif self.option('fq_type') == "SE" and row_num != 2:
                raise OptionError("SE序列list文件应该包括文件名", code = "25000504")

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def get_list(self):
        list_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
        file_sample = FileSampleFile()
        file_sample.set_path(list_path)
        samples = file_sample.get_list()
        self.logger.info("sample is {}".format(samples))
        return samples

    def fastp_run(self):
        n = 1
        self.samples = self.get_list()
        for f in sorted(self.samples):
            fq_s = os.path.join(self.option("fastq_dir").prop["path"], self.samples[f])
            fastp = self.add_tool('ref_rna_v2.fastp')
            self.step.add_steps('fastp_{}'.format(n))
            fastp.set_options({
                "in": fq_s,
            })
            step = getattr(self.step, 'fastp_{}'.format(n))
            step.start()
            fastp.on("end", self.finish_update, "fastp_{}".format(n))
            n += 1
            self.fastp.append(fastp)
        if len(self.fastp) == 1:
            self.fastp[0].on(self.clipper_run_phred64)
            self.fastp[0].run()
        else:
            self.on_rely(self.fastp, self.clipper_run_phred64)
            for tool in self.fastp:
                tool.run()

    def fastp_output(self):
        fastp_output = self.work_dir + "/fastp_output"
        if os.path.exists(fastp_output):
                    shutil.rmtree(fastp_output)
        os.mkdir(fastp_output)
        files = glob.glob(r"{}/Fastp*/output/*".format(self.work_dir))
        for file in files:
            link = os.path.join(fastp_output, os.path.basename(file))
            os.link(file, link)

    def clipper_run(self):
        n = 1
        self.samples = self.get_list()
        for f in sorted(self.samples):
            fq_s = os.path.join(self.option("fastq_dir").prop["path"], self.samples[f])
            clipper = self.add_tool('ref_rna_v2.fastx_clipper')
            self.step.add_steps('clipper_{}'.format(n))
            clipper.set_options({
                "fastq_s": fq_s,
            })
            step = getattr(self.step, 'clipper_{}'.format(n))
            step.start()
            clipper.on("end", self.finish_update, "clipper_{}".format(n))
            clipper.on("end", self.adapt, f)
            clipper.on("end", self.sickle_se_run, f)
            n += 1
            self.clipper.append(clipper)
        if len(self.clipper) == 1:
            self.clipper[0].on(self.adapt_write)
            self.clipper[0].run()
        else:
            self.on_rely(self.clipper, self.adapt_write)
            for tool in self.clipper:
                tool.run()

    def clipper_run_phred64(self):
        self.fastp_output()
        n = 1
        self.samples = self.get_list()
        for f in sorted(self.samples):
            fq_s = os.path.join(self.work_dir + "/fastp_output", self.samples[f])
            clipper = self.add_tool('ref_rna_v2.fastx_clipper')
            self.step.add_steps('clipper_{}'.format(n))
            clipper.set_options({
                "fastq_s": fq_s,
            })
            step = getattr(self.step, 'clipper_{}'.format(n))
            step.start()
            clipper.on("end", self.finish_update, "clipper_{}".format(n))
            clipper.on("end", self.adapt, f)
            clipper.on("end", self.sickle_se_run, f)
            n += 1
            self.clipper.append(clipper)
        if len(self.clipper) == 1:
            self.clipper[0].on(self.adapt_write)
            self.clipper[0].run()
        else:
            self.on_rely(self.clipper, self.adapt_write)
            for tool in self.clipper:
                tool.run()

    def seqprep_run(self):
        self.samples = self.get_list()
        n = 1
        for f in sorted(self.samples):
            fq_l = os.path.join(self.option("fastq_dir").prop["path"], self.samples[f]["l"])
            fq_r = os.path.join(self.option("fastq_dir").prop["path"], self.samples[f]["r"])
            seqprep = self.add_tool('prok_rna.seq_prep')
            self.step.add_steps('seqprep_{}'.format(n))
            seqprep.set_options({
                "fastq_l": fq_l,
                "fastq_r": fq_r,
                "quality_score_system": self.option("quality_score_system")
            })
            step = getattr(self.step, 'seqprep_{}'.format(n))
            step.start()
            seqprep.on("end", self.finish_update, "seqprep_{}".format(n))
            seqprep.on("end", self.adapt, f)
            seqprep.on("end", self.sickle_pe_run, f)
            n += 1
            self.seqprep.append(seqprep)
        self.logger.info(self.seqprep)
        if len(self.seqprep) == 1:
            self.seqprep[0].on("end", self.adapt_write)
            self.seqprep[0].run()
        else:
            self.on_rely(self.seqprep, self.adapt_write)
            for tool in self.seqprep:
                tool.run()

    def sickle_se_run(self, event):
        obj = event["bind_object"]
        os.rename(obj.output_dir + "/clip_s.fastq", obj.output_dir + "/" + event["data"] + "_clip_s.fastq")
        clip_s = os.path.join(obj.output_dir, event["data"] + "_clip_s.fastq")
        self.logger.info(clip_s)
        sickle = self.add_tool('sequence.sickle')
        self.step.add_steps('sickle_{}'.format(self.end_times))
        sickle.set_options({
            "fq_type": self.option("fq_type"),
            "fastq_s": clip_s,
            "quality": self.option("quality_q"),
            "length": self.option("length_q")
        })
        step = getattr(self.step, 'sickle_{}'.format(self.end_times))
        step.start()
        sickle.on("end", self.finish_update, 'sickle_{}'.format(self.end_times))
        sickle.on("end", self.set_output, event["data"])
        sickle.run()
        self.sickle.append(sickle)

    def sickle_pe_run(self, event):
        obj = event["bind_object"]
        sickle = self.add_tool('sequence.sickle')
        self.step.add_steps('sickle_{}'.format(self.end_times))
        sickle.set_options({
            "fq_type": self.option("fq_type"),
            "fastq_l": obj.option("seqprep_l"),  # modified by shijin on 20170623，减少阻塞
            "fastq_r": obj.option("seqprep_r"),
            "quality": self.option("quality_q"),
            "length": self.option("length_q")
        })
        step = getattr(self.step, 'sickle_{}'.format(self.end_times))
        step.start()
        sickle.on("end", self.finish_update, 'sickle_{}'.format(self.end_times))
        sickle.on("end", self.set_output, event["data"])
        sickle.run()
        self.sickle.append(sickle)

    def rename(self, event):
        obj = event["bind_object"]
        for f in os.listdir(obj.output_dir):
            old_name = os.path.join(obj.output_dir, f)
            new_name = os.path.join(obj.output_dir, event["data"] + "_" + f)
            os.rename(old_name, new_name)

    def adapt(self, event):
        obj = event["bind_object"]
        adapt_file = obj.work_dir + "/adapter.xls"
        if os.path.exists(adapt_file):
            with open(adapt_file, "r") as f:
                f.readline()
                adapt_rate = f.next().split()[-1]
                self.adapt_rate.append(["{}".format(event["data"]), adapt_rate])

    def adapt_write(self):
        with open(self.output_dir + "/adapter.xls", "w") as w:
            for a in self.adapt_rate:
                w.write("{}\t{}\n".format(a[0], a[1]))

    def set_output(self, event):
        self.logger.info("set output {}".format(event["data"]))
        obj = event["bind_object"]
        if self.end_times < len(self.samples):
            self.end_times += 1
        for f in os.listdir(obj.output_dir):
            if not f.startswith("sickle"):
                continue
            old_name = os.path.join(obj.output_dir, f)
            new_name = os.path.join(obj.output_dir, event["data"] + "_" + f)
            os.rename(old_name, new_name)
        if self.end_times == len(self.samples):
            sickle_dir = os.path.join(self.output_dir, "sickle_dir")
            sickle_r_dir = os.path.join(self.work_dir, "sickle_r_forRSEM")
            sickle_l_dir = os.path.join(self.work_dir, "sickle_l_forRSEM")
            seqprep_dir = os.path.join(self.work_dir, "seqprep_dir")
            clip_dir = os.path.join(self.work_dir, "clip_dir")
            dir_list = [sickle_dir, seqprep_dir, clip_dir, sickle_r_dir, sickle_l_dir]
            for d in dir_list:
                if os.path.exists(d):
                    shutil.rmtree(d)
                os.mkdir(d)
            sickle_out = []
            seqprep_out = []
            clip_out = []
            for sic in self.sickle:
                for f in os.listdir(sic.output_dir):
                    f_path = os.path.join(sic.output_dir, f)
                    sickle_out.append(f_path)
            for seq in self.seqprep:
                for f in os.listdir(seq.output_dir):
                    f_path = os.path.join(seq.output_dir, f)
                    seqprep_out.append(f_path)
            for clip in self.clipper:
                for f in os.listdir(clip.output_dir):
                    f_path = os.path.join(clip.output_dir, f)
                    clip_out.append(f_path)
            self.logger.info(os.path.join(sickle_dir, "list.txt"))
            with open(os.path.join(sickle_dir, "list.txt"), "w") as w:
                for f in sickle_out:
                    f_name = f.split("/")[-1]
                    if "sickle_r.fastq" in f:
                        sample_name = f_name.split("_sickle_r.fastq")[0]
                        w.write("{}\t{}\t{}\n".format(f_name, sample_name, "r"))
                        #os.link(f, os.path.join(sickle_r_dir, f_name))
                    elif "sickle_l.fastq" in f:
                        sample_name = f_name.split("_sickle_l.fastq")[0]
                        w.write("{}\t{}\t{}\n".format(f_name, sample_name, "l"))
                        #os.link(f, os.path.join(sickle_l_dir, f_name))
                    elif "sickle_s.fastq" in f:
                        sample_name = f_name.split("_sickle_s.fastq")[0]
                        w.write("{}\t{}\n".format(f_name, sample_name))
                    else:
                        w.write("\n")
                    target_path = os.path.join(sickle_dir, f_name)
                    if os.path.exists(target_path):
                        os.remove(target_path)
                    os.link(f, target_path)
            self.logger.info(os.path.join(sickle_dir, "fq_list.txt"))
            fq_list = os.path.join(sickle_dir, "fq_list.txt")
            with open(fq_list, "w") as w:
                sample_list = list()
                for f in sorted(os.listdir(sickle_dir)):
                    if (re.match("(.+)_sickle_*", f)):
                        sample_name = re.match("(.+)_sickle_*", f).group(1)
                        if sample_name not in sample_list:
                            sample_list.append(sample_name)
                            if self.option("fq_type") in ["PE"]:
                                sample_l = sample_name + "_sickle_l.fastq"
                                sample_l_path = os.path.join(sickle_dir, sample_l)
                                sample_r = sample_name + "_sickle_r.fastq"
                                sample_r_path = os.path.join(sickle_dir, sample_r)
                                w.write("{}\t{}\t{}\n".format(sample_name, sample_l_path, sample_r_path))
                            else :
                                sample_s = sample_name + "_sickle_s.fastq"
                                sample_s_path = os.path.join(sickle_dir, sample_s)
                                w.write("{}\t{}\n".format(sample_name, sample_s_path))
            self.option("fq_list", fq_list)
            self.option("sickle_dir", sickle_dir)
            self.end()

    def run(self):
        self.logger.info('{}'.format(self.events))
        if self.option("fq_type") in ["PE"]:
            self.seqprep_run()
        else:
            if self.option("quality_score_system").lower() == "phred+33" or self.option("quality_score_system").lower() == "phred 33":
                self.clipper_run()
            else:
                self.fastp_run()
        for eve in self.events.values():
            self.logger.info('{}'.format(eve.is_start))
        super(HiseqQcModule, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        if self.option("fq_type") == "PE":
            result_dir.add_relpath_rules([
                [r".", "", "结果输出目录"],
                [r"./seqprep_dir/", "文件夹", "PE去接头后fastq文件输出目录"],
                [r"./sickle_dir/", "文件夹", "质量剪切后fastq文件输出目录"]
            ])
        else:
            result_dir.add_relpath_rules([
                [r".", "", "结果输出目录"],
                [r"./clip_dir/", "文件夹", "SE去接头后fastq文件输出目录"],
                [r"./sickle_dir/", "文件夹", "质量剪切后fastq文件输出目录"]
            ])
        super(HiseqQcModule, self).end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "HiseqQc_" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "ref_rna_v2.hiseq_qc",
            "instant": False,
            "options": dict(
                fastq_dir="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/data",
                fq_type="PE",
                quality_score_system="phred+33",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
