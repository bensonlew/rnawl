# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# __last_modify__ = '20191226'


from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil

class SortmernaAgent(Agent):
    """
    Sortmerna基于reads进行丰度统计和功能计算
    """

    def __init__(self, parent):
        super(SortmernaAgent, self).__init__(parent)
        options = [
            {"name": "read1", "type": "infile", "format": "sequence.fastq"},
            {"name": "read2", "type": "infile", "format": "sequence.fastq"},
            {"name": "sample", "type": "string"},
            {"name": "database", "type": 'string', "default": "rfam_5.8s,rfam_5s,arc_16s,arc_23s,bac_16s,bac_23s,euk_18s,euk_28s"},
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("database"):
            raise OptionError("请提供数据库类型！", code="")
        else:
            for file in self.option("database").split(","):
                if file not in ["rfam_5.8s","rfam_5s","arc_16s","arc_23s","bac_16s","bac_23s","euk_18s","euk_28s"]:
                    raise OptionError("请提供正确的数据库类型！{}不正确!".format(file), code="")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 4
        self._memory = '30G'

class SortmernaTool(Tool):
    def __init__(self, config):
        super(SortmernaTool, self).__init__(config)
        self.sormerna = "/bioinfo/metaGenomic/sortmerna-2.1b/sortmerna"
        self.script = "/bioinfo/metaGenomic/sortmerna-2.1b/scripts/"
        self.out =self.work_dir + "/result"
        self.database = ''
        for type in self.option("database").split(","):
            if type in ["rfam_5.8s"]:
                self.database += self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/sortmerna-2.1b/rRNA_databases/rfam-5.8s-database-id98.fasta," + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/sortmerna-2.1b/index/rfam-5.8s-database-id98:"
            elif type in ["rfam_5s"]:
                self.database += self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/sortmerna-2.1b/rRNA_databases/rfam-5s-database-id98.fasta," + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/sortmerna-2.1b/index/rfam-5s-database-id98:"
            elif type in ["arc_16s"]:
                self.database += self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/sortmerna-2.1b/rRNA_databases/silva-arc-16s-id95.fasta," + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/sortmerna-2.1b/index/silva-arc-16s-id95:"
            elif type in ["arc_23s"]:
                self.database += self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/sortmerna-2.1b/rRNA_databases/silva-arc-23s-id98.fasta," + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/sortmerna-2.1b/index/silva-arc-23s-id98:"
            elif type in ["bac_16s"]:
                self.database += self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/sortmerna-2.1b/rRNA_databases/silva-bac-16s-id90.fasta," + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/sortmerna-2.1b/index/silva-bac-16s-id90:"
            elif type in [ "bac_23s"]:
                self.database += self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/sortmerna-2.1b/rRNA_databases/silva-bac-23s-id98.fasta," + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/sortmerna-2.1b/index/silva-bac-23s-id98:"
            elif type in [ "euk_18s"]:
                self.database += self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/sortmerna-2.1b/rRNA_databases/silva-euk-18s-id95.fasta," + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/sortmerna-2.1b/index/silva-euk-18s-id95:"
            elif type in ["euk_28s"]:
                self.database += self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/sortmerna-2.1b/rRNA_databases/silva-euk-28s-id98.fasta," + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/sortmerna-2.1b/index/silva-euk-28s-id98:"


    def merge_fastq(self):
        cmd = "{}merge-paired-reads.sh {} {} {}".format(self.script, self.option("read1").prop["path"], self.option("read2").prop["path"], self.work_dir + "/" + self.option("sample") + ".merge_reads.fq")
        command = self.add_command("merge_fastq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("merge_fastq运行完成！")
        else:
            self.set_error("merge_fastq运行完成运行出错!")

    def run_sortmerna(self):
        """
        description
        :return:
        """
        cmd = "{} -a 8 --ref {} --reads {} --paired_in --other {} --num_alignments 1 --fastx --aligned {} --log".format(self.sormerna, self.database, self.work_dir + "/" + self.option("sample") + ".merge_reads.fq", self.work_dir + "/" + self.option("sample") + ".unalign", self.work_dir + "/" + self.option("sample") + ".align")
        command = self.add_command("run_sortmerna", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_sortmerna运行完成！")
        else:
            self.set_error("run_sortmerna运行完成运行出错!")

    def unmerge_align_fastq(self):
        cmd = "{}unmerge-paired-reads.sh {} {} {}".format(self.script, self.work_dir + "/" + self.option("sample") + ".align.fq", self.output_dir + "/" + self.option("sample") + ".align.1.fq", self.output_dir + "/" + self.option("sample") + ".align.2.fq")
        command = self.add_command("unmerge_align_fastq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("unmerge_align_fastq运行完成！")
        else:
            self.set_error("unmerge_align_fastq运行完成运行出错!")

    def unmerge_unalign_fastq(self):
        cmd = "{}unmerge-paired-reads.sh {} {} {}".format(self.script, self.work_dir + "/" + self.option("sample") + ".unalign.fq", self.output_dir + "/" + self.option("sample") + ".unalign.1.fq", self.output_dir + "/" + self.option("sample") + ".unalign.2.fq")
        command = self.add_command("unmerge_unalign_fastq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("unmerge_unalign_fastq运行完成！")
        else:
            self.set_error("unmerge_unalign_fastq运行完成运行出错!")


    def run(self):
        super(SortmernaTool, self).run()
        self.merge_fastq()
        self.run_sortmerna()
        self.unmerge_align_fastq()
        self.unmerge_unalign_fastq()
        self.end()