# -*- coding: utf-8 -*-
# __author__ = 'shijin'

from __future__ import division
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from collections import defaultdict
import os
import shutil

class HiseqOutputAgent(Agent):
    """
    cog功能分类
    """
    # dir_list = [sickle_dir, seqprep_dir, clip_dir, sickle_r_dir, sickle_l_dir]
    def __init__(self, parent):
        super(HiseqOutputAgent, self).__init__(parent)
        options = [
            {"name": "sickle_dir", "type": "string", "default": ""},
            {"name": "seqprep_dir", "type": "string", "default": ""},
            {"name": "clip_dir", "type": "string", "default": ""},
            {"name": "sickle_r_dir", "type": "string", "default": ""},
            {"name": "sickle_l_dir", "type": "string", "default": ""},
            {"name": "fq_type", "type":"string", "default": "PE"},
            {"name": "sickle_out", "type": "string", "default":""},
            {"name": "seqprep_out", "type": "string", "default":""},
            {"name": "clip_out", "type": "string", "default":""},
            {"name": "fq_s_out", "type": "outfile", "format": "sequence.fastq"},  # SE所有样本cat集合
            {"name": "fq_r_out", "type": "outfile", "format": "sequence.fastq"},  # PE所有右端序列样本cat集合
            {"name": "fq_l_out", "type": "outfile", "format": "sequence.fastq"},  # PE所有左端序列样本cat集合
            {"name": "clip_dir_out", "type": "outfile", "format": "sequence.fastq_dir"},  # SE去接头输出结果文件夹
            {"name": "sickle_dir_out", "type": "outfile", "format": "sequence.fastq_dir"},  # 质量剪切输出结果文件夹(包括左右段)
            {"name": "sickle_r_dir_out", "type": "outfile", "format": "sequence.fastq_dir"},  # 质量剪切右端输出结果文件夹
            {"name": "sickle_l_dir_out", "type": "outfile", "format": "sequence.fastq_dir"},  # 质量剪切左端输出结果文件夹
            {"name": "seqprep_dir_out", "type": "outfile", "format": "sequence.fastq_dir"},  # PE的去接头输出结果文件

        ]
        self.add_option(options)
        self.step.add_steps("hiseq_output")
        self.on('start', self.start_hiseq_output)
        self.on("end", self.end_hiseq_output)

    def start_hiseq_output(self):
        self.step.hiseq_output.start()
        self.step.update()

    def end_hiseq_output(self):
        self.step.hiseq_output.finish()
        self.step.update()

    def check_options(self):
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = "4G"


class HiseqOutputTool(Tool):
    def __init__(self, config):
        super(HiseqOutputTool, self).__init__(config)

    def cmd1(self):
        sickle_dir = self.option("sickle_dir").prop["path"]
        sickle_r_dir = self.option("sickle_r_dir").prop["path"]
        sickle_l_dir = self.option("sickle_l_dir").prop["path"]
        seqprep_dir = self.option("seqprep_dir").prop["path"]
        clip_dir = self.option("clip_dir").prop["path"]
        dir_list = [sickle_dir, seqprep_dir, clip_dir, sickle_r_dir, sickle_l_dir]
        # self.logger.info(dir_list)
        for d in dir_list:
            if os.path.exists(d):
                shutil.rmtree(d)
            os.mkdir(d)
        sickle_out = self.option("sickle_out").split(";")
        seqprep_out = self.option("seqprep_out").split(";")
        clip_out = self.option("clip_out").split(";")
        self.logger.info(os.path.join(sickle_dir, "list.txt"))
        with open(os.path.join(sickle_dir, "list.txt"), "w") as w:
            for f in sickle_out:
                f_name = f.split("/")[-1]
                if "sickle_r.fastq" in f:
                    sample_name = f_name.split("_sickle_r.fastq")[0]
                    w.write("{}\t{}\t{}\n".format(f_name, sample_name, "r"))
                    os.link(f, os.path.join(sickle_r_dir, f_name))
                elif "sickle_l.fastq" in f:
                    sample_name = f_name.split("_sickle_l.fastq")[0]
                    w.write("{}\t{}\t{}\n".format(f_name, sample_name, "l"))
                    os.link(f, os.path.join(sickle_l_dir, f_name))
                elif "sickle_s.fastq" in f:
                    sample_name = f_name.split("_sickle_s.fastq")[0]
                    w.write("{}\t{}\n".format(f_name, sample_name))
                else:
                    w.write("\n")
                target_path = os.path.join(sickle_dir, f_name)
                if os.path.exists(target_path):
                    os.remove(target_path)
                os.link(f, target_path)
        self.option("sickle_dir_out", sickle_dir)
        if self.option("fq_type") == "PE":
            shutil.rmtree(clip_dir)
            for f in seqprep_out:
                f_name = f.split("/")[-1]
                target_path = os.path.join(seqprep_dir, f_name)
                os.link(f, target_path)
            self.option("seqprep_dir_out").set_path(seqprep_dir)
            self.option('sickle_r_dir_out', sickle_r_dir)
            self.option('sickle_l_dir_out', sickle_l_dir)
            r_files = os.listdir(self.option('sickle_r_dir_out').prop['path'])
            l_files = os.listdir(self.option('sickle_l_dir_out').prop['path'])
            r_file = ' '.join(r_files)
            l_file = ' '.join(l_files)
            os.system('cd {} && cat {} > {}/left.fq && cd {} && cat {} > {}/right.fq'.format(sickle_l_dir, l_file, self.work_dir, sickle_r_dir, r_file, self.work_dir))
            self.logger.info('cd {} && cat {} > {}/left.fq && cd {} && cat {} > {}/right.fq'.format(sickle_l_dir, l_file, self.work_dir, sickle_r_dir, r_file, self.work_dir))
            self.option('fq_l_out', self.work_dir + '/left.fq')
            self.option('fq_r_out', self.work_dir + '/right.fq')
        elif self.option('fq_type') == 'SE':
            shutil.rmtree(seqprep_dir)
            shutil.rmtree(sickle_r_dir)
            shutil.rmtree(sickle_l_dir)
            for f in clip_out:
                f_name = f.split("/")[-1]
                target_path = os.path.join(clip_dir, f_name)
                os.link(f, target_path)
            self.option("clip_dir_out").set_path(clip_dir)
            files = self.option('sickle_dir_out').prop['fastq_basename']
            s_file = ' '.join(files)
            os.system('cd {} && cat {} > {}/single.fq'.format(sickle_dir, s_file, self.work_dir))
            self.option('fq_s_out', self.work_dir + '/single.fq')
            self.logger.info("done")

    def run(self):
        super(HiseqOutputTool, self).run()
        self.cmd1()
        self.end()