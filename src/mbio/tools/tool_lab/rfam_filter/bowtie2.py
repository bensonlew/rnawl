# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'
# for prok_rna_mapping -----20180808
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import glob
from biocluster.core.exceptions import OptionError
import json
import unittest


class Bowtie2Agent(Agent):
    """
    tophat  
    version 2.0
    last_modify: by shicaiping at 20180508
    """
    def __init__(self, parent):
        super(Bowtie2Agent, self).__init__(parent)
        options = [
            {"name": "ref_genome", "type": "string"},
            {"name": "seq_method", "type": "string"},
            {"name": "single_end_reads", "type": "infile", "format": "sequence.fastq"},
            {"name": "left_reads", "type": "infile", "format": "sequence.fastq"},
            {"name": "right_reads", "type": "infile", "format": "sequence.fastq"},
            {"name": "bam_output", "type": "outfile", "format": "align.bwa.bam"},
            # {"name": "seq_num", "type": "int", "default": 1000},  # 抽取评估序列长度
            {"name": "sample", "type": "string"},
            {"name": "strand_specific", "type": "bool", "default": False}
            ]
        self.add_option(options)
        self.step.add_steps('bowtie2')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        # self.step.bowtie2.start()
        self.step.update()

    def step_end(self):
        # self.step.bowtie2.finish()
        self.step.update()

    def check_options(self):
        self.logger.info("开始运行check步骤")
        if self.option("seq_method") == "PE":
            if self.option("single_end_reads").is_set:
                raise OptionError("上传的是单端测序的序列，请上传双端序列", code = "35002501")
            elif not (self.option("left_reads").is_set and self.option("right_reads").is_set):
                raise OptionError("缺少某端序列", code = "35002502")
            else:
                pass
        else:
            if not self.option("single_end_reads").is_set:
                raise OptionError("请上传单端序列", code = "35002503")
            elif self.option("left_reads").is_set or self.option("right_reads").is_set:
                raise OptionError("只需要有单端的序列", code = "35002504")
            else:
                pass
        return True

    def set_resource(self):
        self._cpu = 10
        self._memory = '20G'

    def end(self):
        super(Bowtie2Agent, self).end()


class Bowtie2Tool(Tool):
    def __init__(self, config):
        super(Bowtie2Tool, self).__init__(config)
        # self.bowtie_path = self.config.SOFTWARE_DIR + "/bioinfo/align/bowtie2-2.2.9/"
        self.bowtie_path = self.config.SOFTWARE_DIR + "/bioinfo/ref_rna_v2/miniconda2/bin/"
        self.samtool_path = self.config.SOFTWARE_DIR + "/program/Python/bin/samtools"
        self.seqtk = self.config.SOFTWARE_DIR + "/bioinfo/seq/seqtk-master/seqtk"
        self.left_reads = ""
        self.right_reads = ""
        self.single_reads = ""

    def run_build_index(self):
        ref_name = os.path.basename(self.option("ref_genome"))
        new_ref_genome_custom_path = os.path.join(self.work_dir, ref_name)
        if os.path.exists(new_ref_genome_custom_path):
            os.remove(new_ref_genome_custom_path)
        os.link(self.option("ref_genome"), new_ref_genome_custom_path)
        ref_genome_custom_path = os.path.split(new_ref_genome_custom_path)[0]
        self.ref_path = os.path.join(ref_genome_custom_path, "ref_index")  # 保证生成的索引都在和基因组同一级的目录下
        cmd = "{}/bowtie2-build {} {}".format(self.bowtie_path, new_ref_genome_custom_path, self.ref_path)
        index_command_obj = self.add_command("build_index", cmd)
        index_command_obj.software_dir = ""
        index_command_obj.run()
        self.wait(index_command_obj)
        if index_command_obj.return_code == 0:
            self.logger.info("索引建立完成")
            self.run_tophat(self.ref_path)
        else:
            self.set_error("索引建立出错", code = "35002505")

    def run_prepare_fastq(self):
        self.ref_path = self.option("ref_genome")
        self.run_tophat(self.ref_path)


    def extract_fastq(self, fq_r, fq_f, n_lines=1000):
        cmd = "{} sample -s11 {} {} > {}".format(self.seqtk, fq_r, str(n_lines), fq_f)
        command = self.add_command("seqtk_extract_fastq", cmd, shell=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("seqtk提取{}条fastq序列运行完成".format(str(n_lines)))
        else:
            self.set_error("seqtk提取{}条fastq序列运行失败".format(str(n_lines)))
        return fq_f

    def extract_fastq_r(self, fq_r, fq_f, n_lines=1000):
        cmd = "{} sample -s11 {} {} > {}".format(self.seqtk, fq_r, str(n_lines), fq_f)
        command = self.add_command("extract_fastq_r", cmd, shell=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("seqtk提取{}条fastq序列运行完成".format(str(n_lines)))
        else:
            self.set_error("seqtk提取{}条fastq序列运行失败".format(str(n_lines)))
        return fq_f

    def run_tophat(self, index_ref):
        self.logger.info(self.option("sample"))
        name = self.option("sample")
        if name.find("_sickle") != -1:
            name = name[:-9]
        if self.option("seq_method") == "PE" and self.option("strand_specific"):
            cmd = '{}/bowtie2 --fr -x {} --un-conc {}_input_rmrRNA.fastq -1 {} -2 {}  -p 10 -S {}.sam > {}.stat\n'.format(self.bowtie_path, self.ref_path, name,self.option("left_reads").prop['path'],  self.option("right_reads").prop['path'], self.work_dir + '/' + name, self.work_dir + '/' + name)
        elif self.option("seq_method") == "PE" and not self.option("strand_specific"):
            cmd = '{}/bowtie2 -x {}  --un-conc {}_input_rmrRNA.fastq -1 {} -2 {}  -p 10 -S {}.sam > {}.stat\n'.format(self.bowtie_path, self.ref_path,name, self.option("left_reads").prop['path'], self.option("right_reads").prop['path'], self.work_dir + '/' + name, self.work_dir + '/' + name)
        else:
            cmd = '{}/bowtie2 -x {} --un-conc {}_input_rmrRNA.fastq {}  -p 10 -S {}.sam > {}.stat\n'.format(self.bowtie_path, self.ref_path, name,self.option("single_end_reads").prop['path'], self.work_dir + '/' + name, self.work_dir + '/' + name)
        cmd += '{} view -bS -F 0x4 {}.sam -o {}.tmp.bam\n'.format(self.samtool_path, self.work_dir + '/' + name, self.work_dir + '/' + name)
        cmd += '{} sort {}.tmp.bam -o {}.bam'.format(self.samtool_path, self.work_dir + '/' + name, self.work_dir + '/' + name)
        with open(self.work_dir + '/bowtie2.bash', 'w') as bow_w:
            bow_w.write(cmd)
        cmd = "bash {}".format(self.work_dir + '/bowtie2.bash')
        tophat_command = self.add_command("bowtie2", cmd)
        tophat_command.software_dir = "/bin"
        self.logger.info("开始运行tophat")
        tophat_command.run()
        self.wait()
        if tophat_command.return_code == 0:
            outfile_bam = os.path.join(self.work_dir, name + ".bam")
            outfile_stat = os.path.join(self.work_dir, "bowtie2.o")
            if os.path.exists(self.output_dir + "/" + name + ".bam"):
                os.remove(self.output_dir + "/" + name + ".bam")
            os.link(outfile_bam, self.output_dir + "/" + name + ".bam")
            os.link(outfile_stat, self.output_dir + "/" + name + ".stat")
        return True

    def set_output(self):
        out_fastqs = glob.glob(self.work_dir+"/*rmrRNA*fastq")
        for fastq in out_fastqs:
            os.link(fastq, os.path.join(self.output_dir,os.path.basename(fastq)))

    def run(self):
        super(Bowtie2Tool, self).run()
        self.run_prepare_fastq()
        # self.run_build_index()
        self.set_output()
        self.end()


class TestFunction(unittest.TestCase):
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "Botie2Rfamfilter" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.rfam_filter.bowtie2",
            "instant": False,
            "options": dict(
                ref_genome ="/mnt/ilustre/users/sanger-dev/app/database/Annotation/other2019/rfam14.1/Rfam_bowtie2_index/Rfam",
                left_reads="/mnt/ilustre/users/sanger-dev/workspace/20210203/Refrna_tsg_249609/remote_input/fastq_dir/rawdata/CL1.1.fq",
                right_reads = "/mnt/ilustre/users/sanger-dev/workspace/20210203/Refrna_tsg_249609/remote_input/fastq_dir/rawdata/CL1.2.fq",
                seq_method = "PE",
                sample ="test"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

