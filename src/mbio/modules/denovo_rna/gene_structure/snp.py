#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class SnpModule(Module):
    """
    samtools:sam格式文件处理工具
    varscan:snp calling 工具
    version 1.0
    author: qindanhua
    last_modify: 2016.07.25
    """
    def __init__(self, work_id):
        super(SnpModule, self).__init__(work_id)
        options = [
            {"name": "bed", "type": "infile", "format": "denovo_rna.gene_structure.bed"},  # bed格式文件
            {"name": "bam", "type": "infile", "format": "align.bwa.bam,align.bwa.bam_dir"},  # bam格式文件,排序过的
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},  # 参考序列
            # {"name": "pileup", "type": "outfile", "format": "gene_structure.pileup"},  # pileup格式文件
        ]
        self.add_option(options)
        self.bam_files = []
        self.end_info = 1
        self.samtools = []
        self.varscan = []

    def check_options(self):
        """
        检查参数
        """
        if not self.option("ref_fasta").is_set:
            raise OptionError("请传入参考序列")
        if not self.option("bam").is_set:
            raise OptionError("bam文件")

    def samtools_finish_update(self, event):
        # obj = event['bind_object']
        # self.logger.info(event)
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def samtools_run(self):
        self.bam_files = []
        if self.option("bam").format == "align.bwa.bam":
            self.bam_files.append(self.option("bam").prop["path"])
        else:
            self.bam_files = glob.glob("{}/*.bam".format(self.option("bam").prop["path"]))
        n = 1
        self.logger.info(len(self.bam_files))
        for f_path in self.bam_files:
            sample_name = f_path.split("/")[-1].split(".")[0]
            self.logger.info(f_path)
            samtools = self.add_tool('denovo_rna.gene_structure.samtools')
            self.step.add_steps('samtools_{}'.format(n))
            samtools.set_options({
                "ref_fasta": self.option('ref_fasta').prop["path"],
                "in_bam": f_path,
                "method": "mpileup"
            })
            step = getattr(self.step, 'samtools_{}'.format(n))
            step.start()
            samtools.on("end", self.samtools_finish_update, 'samtools_{}'.format(n))
            self.logger.info(n)
            samtools.on("end", self.varscan_single, sample_name)
            n += 1
            self.samtools.append(samtools)
        self.logger.info(self.samtools)
        if len(self.samtools) == 1:
            self.samtools[0].run()
        else:
            for tool in self.samtools:
                tool.run()

    def varscan_single(self, event):
        self.logger.info("varscan")
        obj = event["bind_object"]
        pileup_path = obj.option("pileup").prop["path"]
        self.logger.info(self.end_info)
        self.logger.info(event["data"])
        # file_name = event["data"].split("/")[-1].split(".")[0]
        varscan = self.add_tool('denovo_rna.gene_structure.varscan')
        self.step.add_steps('varscan_{}'.format(self.end_info))
        varscan.set_options({
                "pileup": pileup_path,
                "method": "pileup2snp",
                "bed": self.option("bed").prop["path"]
            })
        step = getattr(self.step, 'varscan_{}'.format(self.end_info))
        step.start()
        varscan.on("end", self.samtools_finish_update, 'varscan_{}'.format(self.end_info))
        # varscan.on("end", self.rename, event["data"])
        varscan.on("end", self.set_output, event["data"])
        self.logger.info("varscanrun")
        varscan.run()
        self.varscan.append(varscan)

    # def rename(self, event):
    #     obj = event["bind_object"]
    #     for f in os.listdir(obj.output_dir):
    #         old_name = os.path.join(obj.output_dir, f)
    #         new_name = os.path.join(obj.output_dir, event["data"] + "." + f)
    #         os.rename(old_name, new_name)

    def set_output(self, event):
        self.logger.info("set output")
        obj = event["bind_object"]
        for f in os.listdir(obj.output_dir):
            old_name = os.path.join(obj.output_dir, f)
            new_name = os.path.join(obj.output_dir, event["data"] + "." + f)
            os.rename(old_name, new_name)
        if self.end_info < len(self.bam_files):
            self.logger.info(self.end_info)
            self.end_info += 1
        elif self.end_info == len(self.bam_files):
            for f in os.listdir(self.output_dir):
                f_path = os.path.join(self.output_dir, f)
                if os.path.isdir(f_path):
                    shutil.rmtree(f_path)
                else:
                    os.remove(f_path)
            for var in self.varscan:
                for f in os.listdir(var.output_dir):
                    f_path = os.path.join(var.output_dir, f)
                    target_path = os.path.join(self.output_dir, f)
                    os.link(f_path, target_path)
            self.logger.info("done")
            self.end()

    def run(self):
        """
        运行
        :return:
        """
        # super(SnpModule, self).run()
        self.samtools_run()
        super(SnpModule, self).run()
        # self.on_rely(self.varscan, self.set_output)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [r".", "", "结果输出目录"]
        ])
        if self.option("bed").is_set:
            result_dir.add_regexp_rules([
                [r".*snp_position_stat\.xls", "xls", "样本snp编码位置信息统计表"],
                [r".*snp_type_stat\.xls", "xls", "样本snp类型统计表"],
                [r".*snp\.xls", "xls", "样本snp信息表"]
            ])
        else:
            result_dir.add_regexp_rules([
                [r".*pileup_out\.xls", "xls", "样本snp信息表"]
            ])
        # print self.get_upload_files()
        super(SnpModule, self).end()
