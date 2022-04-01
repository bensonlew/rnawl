# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
# version 1.0

import os
import re
import glob
import subprocess
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest


class RepeatmaskerAgent(Agent):
    """
    RepeatModeler 构建denovo数据库
    RepeatMasker 查找散在重复序列(interspersed repeats)
    TRF 寻找串联重复序列 (tandem repeat)
    """

    def __init__(self, parent):
        super(RepeatmaskerAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "small_rna.fasta"},  # 组装拼接好的scaffold文件
            {"name": "gff", "type": "outfile", "format": "gene_structure.gff3"},  # 预测生成的GFF格式
            {"name": "repeat", "type": "outfile", "format": "small_rna.fasta"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome")

    def set_resource(self):
        self._cpu = 10
        self._memory = '100G'
        self._memory_increase_step = 150

    def end(self):
        super(RepeatmaskerAgent, self).end()


class RepeatmaskerTool(Tool):
    def __init__(self, config):
        super(RepeatmaskerTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        self.sample_name = (self.genome_fasta).split('/')[-1].rsplit('.')[0]
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/small_rna/"
        self.repeatmasker_path = "/bioinfo/Genomic/Sofware/RepeatMasker/"
        self.cufflinks_path = '/bioinfo/rna/cufflinks-2.2.1/'
        self.bedtools = "/bioinfo/rna/bedtools2-master/bin/"
        self.perl_bin = self.config.SOFTWARE_DIR + "/program/perl/perls/perl-5.24.0/bin"
        self.set_environ(PATH=self.perl_bin)


    def run_repeatmasker(self):
        cmd = "{}RepeatMasker -parallel 50 -engine ncbi -nolow -no_is -norna -dir {} {}".format(
            self.repeatmasker_path, self.work_dir, self.genome_fasta)
        command = self.add_command("repeatmasker", cmd, ignore_error=True).run()
        self.wait(command)
        self.logger.info(command.return_code)
        if command.return_code == 0:
            self.logger.info("Repeatmasker运行完成")
            #self.run_repeat_to_gff()
        elif command.return_code == 255:
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("Repeatmasker运行出错!")

    def run_repeat_to_gff(self):
        file = (self.genome_fasta).split('/')[-1] + ".out"
        os.link(self.work_dir + "/" + file, self.work_dir + "/" + self.sample_name + ".out")
        os.link(self.work_dir + "/" + (self.genome_fasta).split('/')[-1] + ".tbl",
                self.work_dir + "/" + self.sample_name + ".tbl")
        file = self.work_dir + "/" + self.sample_name + ".out"
        cmd = '{} {}repeat_to_gff.pl {}'.format(self.perl_path, self.perl_script, file)
        self.logger.info(cmd)
        command = self.add_command("run_repeat_to_gff", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("生成gff文件运行完成")
            self.run_extract_seq()
        else:
            self.set_error("生成gff文件运行出错!")

    def run_extract_seq(self):
        """
        根据repeatmasker的输出文件提取fasta序列
        """
        gff = self.work_dir + "/" + self.sample_name + ".out.gff"
        index = self.work_dir + "/" + self.sample_name + ".out.bed"
        with open (gff, "r") as f, open (index, "w") as w:
            f.readline()
            for line in f:
                items = line.strip().split("\t")
                name = items[8].split(";")[0].split("=")[1]
                item = [items[0], str(int(items[3])-1), items[4], name, items[5], items[6]]
                w.write("\t".join(item) + "\n")
        cmd = self.bedtools + "bedtools getfasta -fi %s -bed %s -fo %s -name" % (self.genome_fasta, index, self.work_dir + "/" + self.sample_name + ".out.fa")
        command = self.add_command("extract_fasta", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("fasta序列提取完成")
            self.set_output()
        else:
            self.set_error("fasta序列提取出错!")

    def set_output(self):
        gff = self.work_dir + "/" + self.sample_name + ".out.gff"
        repeat = self.work_dir + "/" + self.sample_name + ".out.fa"
        repeat_detail = self.work_dir + "/" + self.sample_name + ".out"
        gff_path = self.output_dir + "/" + self.sample_name + ".SSR.gff"
        repeat_path = self.output_dir + "/" + self.sample_name + ".repeat.fa"
        repeat_detail_path = self.output_dir + "/" + self.sample_name + ".detail.xls"
        if os.path.exists(gff_path):
            os.remove(gff_path)
        os.link(gff, gff_path)
        if os.path.exists(repeat_path):
            os.remove(repeat_path)
        os.link(repeat, repeat_path)
        if os.path.exists(repeat_detail_path):
            os.remove(repeat_detail_path)
        os.link(repeat_detail, repeat_detail_path)
        self.option('gff', gff_path)
        self.option('repeat', repeat_path)

    def run(self):
        super(RepeatmaskerTool, self).run()
        self.run_repeatmasker()
        self.end()

class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "RepeatMasker_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna.srna.repeatmasker",
            "instant": False,
            "options": dict(
                input_genome="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/dna/test2.fa",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
