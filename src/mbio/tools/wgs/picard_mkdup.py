# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.08

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class PicardMkdupAgent(Agent):
    """
    软件: picard
    picard的MarkDuplicates方法,进行标记重复，针对WGS文库类型
    lasted modified by hongdong@20190321针对大基因组建的csi索引，使用samtools去进行去重复
    """
    def __init__(self, parent):
        super(PicardMkdupAgent, self).__init__(parent)
        options = [
            {"name": "sort_bam_file", "type": "infile", "format": "align.bwa.bam"},  # bam文件
            {"name": "mkdup_method", "type": "string", "default": "picard"},
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sort_bam_file").is_set:
            raise OptionError("请设置bam文件", code="34504301")
        if self.option('mkdup_method') not in ['picard', 'samtools']:
            raise OptionError('mkdup method must be picard or samtools ')

    def set_resource(self):
        self._cpu = 9
        self._memory = "50G"

    def end(self):
        super(PicardMkdupAgent, self).end()


class PicardMkdupTool(Tool):
    def __init__(self, config):
        super(PicardMkdupTool, self).__init__(config)
        self.java_path = "program/sun_jdk1.8.0/bin/java"
        self.picard_path = self.config.SOFTWARE_DIR + "/bioinfo/WGS/picard.jar"
        self.samtools_path = '/bioinfo/align/samtools-1.10/bin/samtools'

    def run_picard_mkdup(self):
        """
        picard MarkDuplicates
        """
        sample_name = os.path.basename(self.option("sort_bam_file").prop["path"]).split(".")[0]
        cmd = "{} -Xmx30G -jar {} MarkDuplicates TMP_DIR={}".format(self.java_path, self.picard_path,
                                                                    self.work_dir + "/MKDUP")
        cmd += " MAX_FILE_HANDLES=100 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=true"
        cmd += " I={} O={}".format(self.option("sort_bam_file").prop["path"], self.work_dir + "/" +
                                   sample_name + ".mkdup.bam")
        cmd += " M={} CREATE_INDEX=TRUE".format(self.work_dir + "/" + sample_name + ".metric")
        command = self.add_command("picard_mkdup", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("picard MarkDuplicates完成")
        else:
            self.set_error("picard MarkDuplicates失败", code="34504303")
        if not os.path.exists(self.output_dir + "/" + sample_name + ".metric"):
            os.link(self.work_dir + "/" + sample_name + ".metric", self.output_dir + "/" + sample_name + ".metric")
        if not os.path.exists(self.output_dir + "/" + sample_name + ".mkdup.bam"):
            os.link(self.work_dir + "/" + sample_name + ".mkdup.bam",
                    self.output_dir + "/" + sample_name + ".mkdup.bam")
        if not os.path.exists(self.output_dir + "/" + sample_name + ".mkdup.bai"):
            os.link(self.work_dir + "/" + sample_name + ".mkdup.bai",
                    self.output_dir + "/" + sample_name + ".mkdup.bai")

    def run_samtool_mkdup(self):
        """
        samtools sort -n  xxx.bam -o xxx.sort.bam
        samtools fixmate -m xxx.sort.bam xxx.fixmate.bam
        samtools sort  xxx.fixmate.bam -o xxx.positionsort.bam
        samtools markdup -r xxx.positionsort.bam xxx.markdup.bam
        :return:
        """
        sample_name = os.path.basename(self.option("sort_bam_file").prop["path"]).split(".")[0]
        self.samtools_sort(sample_name)
        self.samtools_fixmate(sample_name)
        self.samtools_sort(sample_name, False)
        self.samtools_markdup(sample_name)
        if not os.path.exists(self.output_dir + "/" + sample_name + ".mkdup.bam"):
            os.link(self.work_dir + "/" + sample_name + ".mkdup.bam",
                    self.output_dir + "/" + sample_name + ".mkdup.bam")

    def samtools_sort(self, sample_name, is_sort_by_n=True):
        if is_sort_by_n:
            cmd = "{} sort -n --output-fmt CRAM {} --reference {} -@ 8 -o {}"\
                .format(self.samtools_path, self.option("sort_bam_file").prop["path"],
                        self.option("ref_fa").prop["path"], '{}.sort.bam'.format(sample_name))
            cmd_name = 'samtools_sort_n'
        else:
            cmd = "{} sort --output-fmt CRAM {} --reference {} -@ 8 -o {}"\
                .format(self.samtools_path, '{}.fixmate.bam'.format(sample_name), self.option("ref_fa").prop["path"],
                                            '{}.positionsort.bam'.format(sample_name))
            cmd_name = 'samtools_sort'
        command = self.add_command(cmd_name, cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{}完成".format(cmd_name))
        else:
            self.set_error("{}失败".format(cmd_name))

    def samtools_fixmate(self, sample_name):
        cmd = "{} fixmate -m -@ 8 --output-fmt CRAM --reference {} {} {}"\
            .format(self.samtools_path, self.option("ref_fa").prop["path"], '{}.sort.bam'.format(sample_name),
                    '{}.fixmate.bam'.format(sample_name))
        command = self.add_command("samtools_fixmate", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("samtools_fixmate完成")
        else:
            self.set_error("samtools_fixmate失败")

    def samtools_markdup(self, sample_name):
        cmd = "{} markdup -r --output-fmt CRAM -@ 8 --reference {} {} {}"\
            .format(self.samtools_path, self.option("ref_fa").prop["path"], '{}.positionsort.bam'.format(sample_name),
                    '{}.mkdup.bam'.format(sample_name))
        command = self.add_command("samtools_markdup", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("samtools_markdup完成")
        else:
            self.set_error("samtools_markdup失败")

    def run(self):
        super(PicardMkdupTool, self).run()
        if self.option("mkdup_method") == 'picard':
            self.run_picard_mkdup()
        else:
            self.run_samtool_mkdup()
        self.end()
