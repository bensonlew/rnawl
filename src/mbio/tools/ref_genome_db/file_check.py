# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.files.sequence.fastq import FastqFile
from mbio.files.sequence.file_sample import FileSampleFile
import os
import re
import unittest

class FileCheckAgent(Agent):
    """
    用于workflow开始之前对所输入的文件进行详细的内容检测
    """

    def __init__(self, parent):
        super(FileCheckAgent, self).__init__(parent)
        options = [
            {"name": "in_gtf", "type": "infile", "format": "ref_genome_db.gtf"},
            {"name": "gff", "type": "infile", "format": "ref_genome_db.gtf"},
            {"name": "genome", "type": "infile", "format": "ref_genome_db.fasta"},
            {"name": "bed", "type": "outfile", "format": "ref_genome_db.bed"},
            {"name": "out_gtf", "type": "outfile", "format": "ref_genome_db.gtf"},
        ]
        self.add_option(options)
        self.step.add_steps("file_check")
        self.on('start', self.start_file_check)
        self.on('end', self.end_file_check)

    def start_file_check(self):
        self.step.file_check.start()
        self.step.update()

    def end_file_check(self):
        self.step.file_check.finish()
        self.step.update()

    def check_option(self):
        if not self.option("in_gtf").is_set and not self.option("gff").is_set:
            raise OptionError("GTF/GFF文件必须至少输入其中一个")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = "15G"


class FileCheckTool(Tool):
    """
    检查输入文件的格式是否符合要求
    """
    def __init__(self, config):
        super(FileCheckTool, self).__init__(config)
        self.gffread_path = "bioinfo/rna/cufflinks-2.2.1/gffread"

    def gff_to_gtf(self):
        self.logger.info("转换gff文件为gtf文件")
        gff = self.option("gff").prop["path"]
        tmp_gtf = os.path.join(self.work_dir, os.path.basename(gff) + ".tmp.gtf")
        gtf = os.path.join(self.work_dir, os.path.basename(gff) + ".gtf")
        genome = self.option("genome").prop["path"]
        to_gtf_cmd = '%s %s -g %s -T -O -o %s' % (self.gffread_path, gff, genome, tmp_gtf)
        command = self.add_command("gff_to_gtf", to_gtf_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行gff_to_gtf完成")
        else:
            self.set_error("运行gff_to_gtf运行出错!")
            return False
        ## 过滤gtf文件，去掉不含gene_id、transcript_id、链信息的注释行
        with open(tmp_gtf, "r") as f1, open(gtf, "w") as w1:
            for line in f1:
                items = line.strip().split("\t")
                if len(items) != 9:
                    continue
                if not re.search(r'transcript_id', items[8]):
                    continue
                if not re.search(r'gene_id', items[8]):
                    continue
                if not re.search(r'[+|-]', items[6]):
                    continue
                if not line.strip().endswith(";"):
                    w1.write(line.strip() + ";\n")
                else:
                    w1.write(line)
        self.option("out_gtf", gtf)

    def gtf_to_bed(self):
        self.logger.info("转换gtf文件为bed文件")
        if self.option("in_gtf").is_set:
            new_gtf = self.work_dir + "/" + os.path.basename(self.option("in_gtf").prop["path"])
            in_gtf = self.option("in_gtf").prop["path"]
            # if os.path.exists(new_gtf):
            #     os.remove(new_gtf)
            # os.link(self.option("in_gtf").prop["path"], new_gtf)
            ## 过滤gtf文件，去掉不含gene_id、transcript_id、链信息的注释行
            with open(in_gtf, "r") as f1, open(new_gtf, "w") as w1:
                for line in f1:
                    items = line.strip().split("\t")
                    if len(items) != 9:
                        continue
                    if not re.search(r'transcript_id', items[8]):
                        continue
                    if not re.search(r'gene_id', items[8]):
                        continue
                    if not re.search(r'[+|-]', items[6]):
                        continue
                    if not line.strip().endswith(";"):
                        w1.write(line.strip() + ";\n")
                    else:
                        w1.write(line)
            self.option("out_gtf").set_path(new_gtf)
        self.option("out_gtf").to_bed()

    def make_true_bed(self):
        '''
        modified by qinjincheng at 20190402
        create a legal bed by ucsc utilities
        '''
        gtf = self.option('out_gtf').prop['path']
        genepred = '{}.genepred'.format(self.option('out_gtf').prop['path'])
        bed = '{}.bed'.format(self.option('out_gtf').prop['path'])
        gtfToGenePred = 'bioinfo/align/ucsc_tools/gtfToGenePred'
        genePredToBed = 'bioinfo/align/ucsc_tools/genePredToBed'
        cmd1 = '{} {} {}'.format(gtfToGenePred, gtf, genepred)
        self.run_code('gtftogenepred', cmd1)
        cmd2 = '{} {} {}'.format(genePredToBed, genepred, bed)
        self.run_code('genepredtobed', cmd2)

    def run_code(self, cmd_name, cmd, shell=False):
        if shell:
            cmd = self.config.SOFTWARE_DIR + '/' + cmd
        command = self.add_command(cmd_name, cmd, shell=shell)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('succeed in running {}'.format(cmd_name))
        elif command.return_code is None:
            self.logger.warn('fail to run {}, try again'.format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info('succeed in rerunning {}'.format(cmd_name))
            else:
                self.set_error('fail to rerun {}, abord'.format(cmd_name))
        else:
            self.set_error('fail to run {}, abord'.format(cmd_name))

    def filter_bed(self):
        ## modified by shicaiping at 20181225
        ## 新增对bed文件的过滤，基于以下两点：
        # ① 基因单条染色体序列长度超过500000000，无法进行Bamdistribution分析，过滤超过500000000部分的转录本
        # ② 过滤exon存在交叉、包含等关系的转录本
        with open(self.option("out_gtf").prop["path"] + ".bed", "r") as f, open(self.option("out_gtf").prop["path"] + ".filter.bed", "w") as w:
            for line in f:
                if int(line.split("\t")[2]) < 500000000:
                    flag = 1
                    start_list = line.strip().split("\t")[11].split(",")[:-1]
                    length_list = line.strip().split("\t")[10].split(",")[:-1]
                    for index, value in enumerate(length_list[:-1]):
                        if ((int(length_list[index]) + int(start_list[index])) >= int(start_list[index + 1])):
                            flag = 0
                            break
                    if flag == 1:
                        w.write(line)
        self.option("bed").set_path(self.option("out_gtf").prop["path"] + ".filter.bed")

    def set_output(self):
        bed = self.option("out_gtf").prop["path"] + ".filter.bed"
        gtf = self.option("out_gtf").prop["path"]
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(bed))):
            os.remove(os.path.join(self.output_dir, os.path.basename(bed)))
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(gtf))):
            os.remove(os.path.join(self.output_dir, os.path.basename(gtf)))
        os.link(bed, os.path.join(self.output_dir, os.path.basename(bed)))
        os.link(gtf, os.path.join(self.output_dir, os.path.basename(gtf)))

    def run(self):
        super(FileCheckTool, self).run()
        if not self.option("in_gtf").is_set:
            self.gff_to_gtf()
        self.gtf_to_bed()
        self.make_true_bed()
        self.filter_bed()
        self.set_output()
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "FileCheck_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_genome_db.file_check",
            "instant": False,
            "options": dict(
                gff="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/musa_acuminata/GCF_000313855.2_ASM31385v2_genomic.gff",
                genome="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/musa_acuminata/GCF_000313855.2_ASM31385v2_genomic.fna",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()