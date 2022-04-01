# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'
# for prok_rna_mapping -----20180808
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import json


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
        python_path = self.config.SOFTWARE_DIR + '/program/Python/bin/'
        self.set_environ(PATH=python_path)
        # self.bowtie_path = self.config.SOFTWARE_DIR + "/bioinfo/align/bowtie2-2.2.9/"
        self.bowtie_path = self.config.SOFTWARE_DIR + "/bioinfo/ref_rna_v2/miniconda2/bin/"
        self.samtool_path = self.config.SOFTWARE_DIR + "/program/Python/bin/samtools"
    
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

    def run_tophat(self, index_ref):
        self.logger.info(self.option("sample"))
        name = self.option("sample")
        if name.find("_sickle") != -1:
            name = name[:-9]
        if self.option("seq_method") == "PE" and self.option("strand_specific"):
            cmd = '{}/bowtie2 --fr -x {} -1 {} -2 {}  -p 10 -S {}.sam > {}.stat\n'.format(self.bowtie_path, self.ref_path, self.option("left_reads").prop['path'], self.option("right_reads").prop['path'], self.work_dir + '/' + name, self.work_dir + '/' + name)
        elif self.option("seq_method") == "PE" and not self.option("strand_specific"):
            cmd = '{}/bowtie2 -x {} -1 {} -2 {}  -p 10 -S {}.sam > {}.stat\n'.format(self.bowtie_path, self.ref_path, self.option("left_reads").prop['path'], self.option("right_reads").prop['path'], self.work_dir + '/' + name, self.work_dir + '/' + name)
        else:
            cmd = '{}/bowtie2 -x {} {}  -p 10 -S {}.sam > {}.stat\n'.format(self.bowtie_path, self.ref_path, self.option("single_end_reads").prop['path'], self.work_dir + '/' + name, self.work_dir + '/' + name)
        # cmd += '{} view -bS -F 0x4 {}.sam -o {}.tmp.bam\n'.format(self.samtool_path, self.work_dir + '/' + name, self.work_dir + '/' + name)
        cmd += '{} view -bS {}.sam -o {}.tmp.bam\n'.format(self.samtool_path, self.work_dir + '/' + name, self.work_dir + '/' + name)   # Obtain unmapped reads
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
    
    def run(self):
        super(Bowtie2Tool, self).run()
        self.run_build_index()
        self.end()
