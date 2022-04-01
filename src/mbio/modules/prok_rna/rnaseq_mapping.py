#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.files.sequence.file_sample import FileSampleFile
import glob
import shutil
import unittest

class RnaseqMappingModule(Module):
    """
    refRNA的比对工具: tophat、hisat
    author: shicaiping
    """
    def __init__(self, work_id):
        super(RnaseqMappingModule, self).__init__(work_id)
        options = [
            {"name": "ref_genome", "type": "string"},  # 参考基因组的路径
            {"name": "seq_method", "type": "string"},  # 双端测序还是单端测序
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # fastq文件夹
            {"name": "bam_output", "type": "outfile", "format": "align.bwa.bam_dir"},  # 输出的bam
            {"name": "bamlist", "type": "outfile", "format": "prok_rna.common"},  # 输出的bamlist
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},
            {"name": "strand_specific", "type": "bool", "default": False}  # 链特异性参数
        ]
        self.add_option(options)
        self.samples = {}
        self.tool_opts = {}
        self.tools = []
        
    def check_options(self):
        """
        检查参数
        """
        if self.option("fastq_dir").is_set:
            self.samples = self.get_list()
            list_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
            row_num = len(open(list_path, "r").readline().split())
            self.logger.info(row_num)
            if self.option('seq_method') == "PE" and row_num != 3:
                raise OptionError("PE序列list文件应该包括文件名、样本名和左右端说明三列", code = "25000901")
            elif self.option('seq_method') == "SE" and row_num != 2:
                raise OptionError("SE序列list文件应该包括文件名、样本名两列", code = "25000902")
        if self.option('seq_method') not in ['PE', 'SE']:
            raise OptionError("请说明序列类型为PE或SE", code = "25000903")
        if not self.option("fastq_dir").is_set and self.option('seq_method') in ["PE"]:
            if self.option("single_end_reads").is_set:
                raise OptionError("您上传的是单端测序的序列，请上传双端序列", code = "25000904")
            elif not (self.option("left_reads").is_set and self.option("right_reads").is_set):
                raise OptionError("您漏了某端序列", code = "25000905")
        if not self.option("fastq_dir").is_set and self.option('seq_method') == "SE":
            if not self.option("single_end_reads").is_set:
                raise OptionError("请上传单端序列", code = "25000906")
            elif self.option("left_reads").is_set or self.option("right_reads").is_set:
                raise OptionError("有单端的序列就够啦", code = "25000907")
        return True
        
    def get_opts(self):
        self.tool_opts = {
            "ref_genome": self.option("ref_genome"),
        }
        if self.option("strand_specific"):
            self.tool_opts.update(
                {
                    "strand_specific": True
                }
            )
        return True
    
    def run(self):
        super(RnaseqMappingModule, self).run()
        self.get_opts()
        self.tool_run()
        
    def tool_run(self):
        if self.option("seq_method") == "PE":
            for f in sorted(self.samples):
                fq_l = os.path.join(self.option('fastq_dir').prop["path"], self.samples[f]["l"])
                fq_r = os.path.join(self.option('fastq_dir').prop["path"], self.samples[f]["r"])
                mapping_tool = self.add_tool('prok_rna.bowtie2')
                self.tool_opts.update({
                    'left_reads': fq_l,
                    'right_reads': fq_r,
                    'sample': f,
                    "seq_method": self.option("seq_method")
                })
                mapping_tool.set_options(self.tool_opts)
                self.tools.append(mapping_tool)

        else:
            for f in sorted(self.samples):
                fq_s = os.path.join(self.option('fastq_dir').prop["path"], self.samples[f])
                mapping_tool = self.add_tool('prok_rna.bowtie2')
                self.tool_opts.update({
                    'single_end_reads': fq_s,
                    'sample': f,
                    "seq_method": self.option("seq_method")
                })
                mapping_tool.set_options(self.tool_opts)
                self.tools.append(mapping_tool)
        self.on_rely(self.tools, self.set_output, 'bowtie2')
        for tool in self.tools:
            tool.run()
            
    def get_list(self):
        list_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
        file_sample = FileSampleFile()
        file_sample.set_path(list_path)
        samples = file_sample.get_list()
        return samples
        
    def end(self):
        super(RnaseqMappingModule, self).end()
        
    def set_output(self, event):
        self.logger.info("set output")
        for f in glob.glob(r"{}/*".format(self.output_dir)):
            if os.path.isdir(f):
                shutil.rmtree(f)
            else:
                os.remove(f)
        if not os.path.exists(self.output_dir + "/bam"):
            os.mkdir(self.output_dir + "/bam")
        if not os.path.exists(self.output_dir + "/stat"):
            os.mkdir(self.output_dir + "/stat")
        new_path_bam = self.output_dir + "/bam"
        new_path_stat = self.output_dir + "/stat"
        for tool in self.tools:
            out_files = os.listdir(tool.output_dir)
            for f in out_files:
                if f.endswith('bam'):
                    f_path = os.path.join(tool.output_dir, f)
                    target = os.path.join(new_path_bam, f)
                    if os.path.exists(target):
                        os.remove(target)
                    os.link(f_path, target)
                if f.endswith('stat'):
                    f_path = os.path.join(tool.output_dir, f)
                    target = os.path.join(new_path_stat, f)
                    if os.path.exists(target):
                        os.remove(target)
                    os.link(f_path, target)
        bam_list = self.output_dir + "/bamlist"
        with open(bam_list, "w") as w:
            for f in sorted(os.listdir(self.output_dir + "/bam")):
                f_path = self.output_dir + "/bam/" + f
                w.write(f_path + "\n")
        self.option("bam_output").set_path(self.output_dir + "/bam")
        self.option("bamlist").set_path(self.output_dir + "/bamlist")
        self.logger.info("done")
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
            "id": "RnaseqMapping_" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "prok_rna.rnaseq_mapping",
            "instant": False,
            "options": dict(
                ref_genome="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/pipline/ref/GCF_000009345.1_ASM934v1_genomic.fna",
                seq_method="PE",
                fastq_dir="/mnt/ilustre/users/sanger-dev/workspace/20180806/Single_HiseqQc_3546/HiseqQc/output/sickle_dir",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()