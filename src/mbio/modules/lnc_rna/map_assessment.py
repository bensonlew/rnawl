#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import shutil
import glob
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import unittest
import pandas as pd
import random

class MapAssessmentModule(Module):
    """
    denovoRNA比对后质量评估:基因覆盖率、比对结果统计、冗余序列分析
    author: qindanhua
    last_modified by shicaiping at 20180509
    """
    def __init__(self, work_id):
        super(MapAssessmentModule, self).__init__(work_id)
        options = [
            {"name": "bed", "type": "infile", "format": "gene_structure.bed"},  # bed格式文件
            {"name": "bam", "type": "infile", "format": "align.bwa.bam,align.bwa.bam_dir"},  # bam格式文件,排序过的
            {"name": "analysis", "type": "string", "default": "saturation,distribution,coverage,chr_stat"},  # 分析类型
            {"name": "quality_satur", "type": "int", "default": 30},  # 测序饱和度分析质量值
            {"name": "quality_dup", "type": "int", "default": 30},  # 冗余率分析质量值
            {"name": "low_bound", "type": "int", "default": 5},  # Sampling starts from this percentile
            {"name": "up_bound", "type": "int", "default": 100},  # Sampling ends at this percentile
            {"name": "step", "type": "int", "default": 5},  # Sampling frequency
            {"name": "rpkm_cutof", "type": "float", "default": 0.01},  # RPKM阈值
            {"name": "min_len", "type": "int", "default": 100}  # Minimum mRNA length (bp).
        ]
        self.add_option(options)
        self.tools = []
        self.files = []

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        """
        analysis = self.option("analysis").split(",")
        for an in analysis:
            if an.lower() in ["saturation", "coverage"]:
                if not self.option("bed").is_set:
                    raise OptionError("请传入bed文件", code = "23700701")
        for an in analysis:
            if an.lower() in ["saturation", "coverage", "distribution", "chr_stat"]:
                self.files = self.get_files()
                if not self.option("bam").is_set:
                    raise OptionError("请传入bam文件", code = "23700702")
            if an.lower() not in ["rrna", "saturation", "coverage", "distribution", "chr_stat"]:
                raise OptionError("所选质量评估分析方法不在范围内", code = "23700703")
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        with open(self.option('bed').path, 'r') as b:
            self.bedlines = len(b.readlines())

    def satur_run(self):
        n = 0
        for f in self.files:
            satur = self.add_tool('denovo_rna.mapping.rpkm_saturation')
            self.step.add_steps('satur{}'.format(n))
            satur.set_options({
                'bam': f,
                "bed": self.option('bed').prop["path"],
                "low_bound": self.option("low_bound"),
                "up_bound": self.option("up_bound"),
                "step": self.option("step"),
                "rpkm_cutof": self.option("rpkm_cutof"),
                "quality": self.option("quality_satur")
                })
            step = getattr(self.step, 'satur{}'.format(n))
            step.start()
            satur.on("end", self.finish_update, 'satur{}'.format(n))
            self.tools.append(satur)
            n += 1

    def coverage_new_run(self):
        bed_file = pd.read_table(self.option('bed').path, header=None, index_col=None, sep='\t')
        lines_number = bed_file.shape[0]
        need_number = lines_number//10
        if lines_number > 10000:
            if need_number < 10000:
                need_number = 10000
        else:
            need_number = lines_number
        split_name = self.work_dir + '/' + os.path.basename(self.option('bed').path) + 'filter'
        need_list = random.sample(range(0, lines_number), need_number)
        need_bed = bed_file.iloc[need_list]
        need_bed.to_csv(split_name, header=None, index=None, sep='\t')
        n = 0
        for f in self.files:
            coverage_new = self.add_tool('denovo_rna.mapping.coverage')
            self.step.add_steps('coverage_{}'.format(n))
            coverage_new.set_options({
                'bam': f,
                'bed': split_name,
                'min_len': self.option('min_len')
            })
            step = getattr(self.step, 'coverage_{}'.format(n))
            step.start()
            coverage_new.on('end', self.finish_update, 'coverage_{}'.format(n))
            self.tools.append(coverage_new)
            n += 1

    def coverage_run(self):
        n = 0
        for f in self.files:
            coverage = self.add_tool('denovo_rna.mapping.coverage')
            self.step.add_steps('coverage_{}'.format(n))
            coverage.set_options({
                'bam': f,
                "bed": self.option('bed').prop["path"],
                "min_len": self.option("min_len")
                })
            step = getattr(self.step, 'coverage_{}'.format(n))
            step.start()
            coverage.on("end", self.finish_update, 'coverage_{}'.format(n))
            self.tools.append(coverage)
            n += 1

    def distribution_run(self):
        n = 0
        for f in self.files:
            distribution = self.add_tool("gene_structure.bam_readsdistribution")
            self.step.add_steps("distribution_{}".format(n))
            distribution.set_options({
                "bam": f,
                "bed": self.option("bed").prop["path"]
            })
            step = getattr(self.step, "distribution_{}".format(n))
            step.start()
            distribution.on("end", self.finish_update, "distribution_{}".format(n))
            self.tools.append(distribution)
            n += 1

    def chr_stat_run(self):
        n = 0
        for f in self.files:
            chr_stat = self.add_tool("gene_structure.chr_distribution")
            self.step.add_steps("chr_distribution_{}".format(n))
            chr_stat.set_options({
                "bam": f
            })
            step = getattr(self.step, "chr_distribution_{}".format(n))
            step.start()
            chr_stat.on("end", self.finish_update, "chr_distribution_{}".format(n))
            self.tools.append(chr_stat)
            n += 1

    def get_files(self):
        files = []
        if self.option("bam").format == "align.bwa.bam":
            files.append(self.option("bam").prop["path"])
        elif self.option("bam").format == "align.bwa.bam_dir":
            for f in glob.glob(r"{}/*.bam".format(self.option("bam").prop["path"])):
                files.append(os.path.join(self.option("bam").prop["path"], f))
        self.logger.info(files)
        return files

    def set_output(self):
        self.logger.info("set output")
        dirs = ["coverage", "saturation", "distribution", "chr_stat"]
        for f in os.listdir(self.output_dir):
            f_path = os.path.join(self.output_dir, f)
            if os.path.exists(f_path):
                if os.path.isdir(f_path):
                    shutil.rmtree(f_path)
                else:
                    os.remove(f_path)
        for d in dirs:
            f_path = os.path.join(self.output_dir, d)
            if not os.path.exists(f_path):
                os.makedirs(f_path)
        for tool in self.tools:
            out = os.listdir(tool.output_dir)
            for f_name in out:
                fp = os.path.join(tool.output_dir, f_name)
                if "eRPKM.xls" in f_name:
                    target = os.path.join(self.output_dir, "saturation", f_name)
                    if os.path.exists(target):
                        os.remove(target)
                    os.link(fp, target)
                if "geneBodyCoverage" in f_name:
                    target = os.path.join(self.output_dir, "coverage", f_name)
                    if os.path.exists(target):
                        os.remove(target)
                    os.link(fp, target)
                if "reads_distribution" in f_name:
                    self.logger.info(fp)
                    target = os.path.join(self.output_dir, "distribution", f_name)
                    if os.path.exists(target):
                        os.remove(target)
                    os.link(fp, target)
                if "chr_stat" in f_name:
                    self.logger.info(fp)
                    target = os.path.join(self.output_dir, "chr_stat", f_name)
                    if os.path.exists(target):
                        os.remove(target)
                    os.link(fp, target)
        self.end()

    def run(self):
        super(MapAssessmentModule, self).run()
        analysis = self.option("analysis").split(",")
        for m in analysis:
            if m.lower() == "saturation":
                self.satur_run()
            if m.lower() == "coverage":
                # self.coverage_run()
                self.coverage_new_run()
            if m.lower() == "distribution":
                self.distribution_run()
            if m.lower() == "chr_stat":
                self.chr_stat_run()

        if len(self.tools) > 1:
            self.on_rely(self.tools, self.set_output)
            for t in self.tools:
                t.run()
        else:
            self.tools[0].on("end", self.set_output)
            self.tools[0].run()

    def end(self):
        super(MapAssessmentModule, self).end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "MapAssessment_lnc_fromtophat" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "lnc_rna.map_assessment",
            "instant": False,
            "options": dict(
                bam="/mnt/ilustre/users/sanger-dev/workspace/20190318/Single_RnaseqMapping_star_test6499/RnaseqMapping/output/bam",
                analysis="saturation,distribution,coverage,chr_stat",
                bed="/mnt/ilustre/users/isanger/workspace/20190222/Single_FileCheck_lnc3075/FileCheck/Homo_sapiens.GRCh38.89.gtf.bed",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
