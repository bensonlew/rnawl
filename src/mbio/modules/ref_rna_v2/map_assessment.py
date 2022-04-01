# -*- coding: utf-8 -*-
# __author__ = 'qindanhua,shicaiping'

import glob
import os
import shutil
import unittest

from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class MapAssessmentModule(Module):
    """
    denovoRNA比对后质量评估:基因覆盖率、比对结果统计、冗余序列分析
    author: qindanhua
    last_modified by shicaiping at 20180509
    """

    def __init__(self, work_id):
        super(MapAssessmentModule, self).__init__(work_id)
        options = [
            {"name": "bed", "type": "infile", "format": "gene_structure.bed"},  # bed file path
            {"name": "bam", "type": "infile", "format": "align.bwa.bam,align.bwa.bam_dir"},  # dir contains sorted bam
            {"name": "analysis", "type": "string", "default": "saturation,distribution,coverage,chr_stat"},  # type
            {"name": "quality_satur", "type": "int", "default": 30},  # minimum mapping quality for an alignment
            {"name": "quality_dup", "type": "int", "default": 30},
            {"name": "low_bound", "type": "int", "default": 5},  # sampling starts from this percentile
            {"name": "up_bound", "type": "int", "default": 100},  # sampling ends at this percentile
            {"name": "step", "type": "int", "default": 5},  # sampling frequency
            {"name": "rpkm_cutof", "type": "float", "default": 0.01},
            {"name": "min_len", "type": "int", "default": 100}  # minimum mRNA length (bp)
        ]
        self.add_option(options)
        self.tools = []
        self.files = []

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def check_options(self):
        analysis = self.option("analysis").split(",")
        for an in analysis:
            if an.lower() in ["saturation", "coverage"]:
                if not self.option("bed").is_set:
                    raise OptionError("请传入bed文件", code="23700701")
        for an in analysis:
            if an.lower() in ["saturation", "coverage", "distribution", "chr_stat"]:
                self.files = self.get_files()
                if not self.option("bam").is_set:
                    raise OptionError("请传入bam文件", code="23700702")
            if an.lower() not in ["saturation", "coverage", "distribution", "chr_stat"]:
                raise OptionError("所选质量评估分析方法不在范围内", code="23700703")

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
                self.coverage_run()
            if m.lower() == "distribution":
                self.distribution_run()
            if m.lower() == "chr_stat":
                self.chr_stat_run()
        self.logger.info('zjx_test')
        self.logger.info(len(self.tools))
        self.logger.info(self.tools)
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
            "id": "MapAssessment_" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "ref_rna_v2.map_assessment",
            "instant": False,
            "options": dict(
                bam="/mnt/ilustre/users/sanger-dev/workspace/20190617/Single_RnaseqMapping_5868/RnaseqMapping/output/bam",
                analysis="saturation,coverage,chr_stat",
                bed="/mnt/ilustre/users/isanger/workspace/20190222/Single_FileCheck_lnc3075/FileCheck/Homo_sapiens.GRCh38.89.gtf.bed",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
