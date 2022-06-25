# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.04

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class SamtoolsIndexAgent(Agent):
    """
    软件: samtools
    samtools的index方法
    """
    def __init__(self, parent):
        super(SamtoolsIndexAgent, self).__init__(parent)
        options = [
            {"name": "sort_bam_file", "type": "infile", "format": "align.bwa.bam"},  # bam文件s
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sort_bam_file").is_set:
            raise OptionError("请设置bam文件", code="34505101")

    def set_resource(self):
        self._cpu = 9
        self._memory = "50G"

    def end(self):
        super(SamtoolsIndexAgent, self).end()


class SamtoolsIndexTool(Tool):
    def __init__(self, config):
        super(SamtoolsIndexTool, self).__init__(config)
        # self.samtools_path = "miniconda2/bin/samtools"
        self.samtools_path = "bioinfo/align/samtools-1.10/samtools"   # 换成新的版本samtools
        self.small_bam = True
        self.bamtosrcam = False
        self.bam_name = ""

    def run_samtools_index(self):
        """
        samtools index
        """
        self.bam_name = os.path.basename(self.option("sort_bam_file").prop["path"])
        if os.path.exists(self.work_dir + "/" + self.bam_name):
            os.remove(self.work_dir + "/" + self.bam_name)
        os.link(self.option("sort_bam_file").prop["path"], self.work_dir + "/" + self.bam_name)
        cmd = "{} index {}".format(self.samtools_path, self.work_dir + "/" + self.bam_name)
        command = self.add_command("samtools_index", cmd, ignore_error=True).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("samtools index完成")
        else:
            self.check_large_bam()
            self.logger.info("self.small_bam: {}".format(self.small_bam))
            if not self.small_bam:
                self.logger.info("开始csi建索引")
                # 注释64-77之前这里是通过判断是否建立csi索引，如果建立csi索引失败的时候，再去转为cram接着建索引，现在直接
                # 大基因组就建cram索引
                # cmd1 = "{} index -c {}".format(self.samtools_path, self.work_dir + "/" + self.bam_name)
                # command1 = self.add_command("samtools_index_csi", cmd1, ignore_error=True).run()
                # self.wait()
                # if command1.return_code == 0:
                #     self.logger.info("samtools_index_csi success")
                # else:
                #     self.logger.info('开始检测是否需要转cram')
                #     self.check_cram_bam()
                #     self.logger.info("bamtosrcam:{}".format(self.bamtosrcam))
                #     if self.bamtosrcam:
                #         self.run_bam2cram()
                #     else:
                #         self.set_error("samtools_index_csi failed", code="34505104")
                self.run_bam2cram()
            else:
                self.set_error("samtools index失败", code="34505101")

    def check_large_bam(self):
        """
        解析samtools_index.o
        :return:
        """
        with open(self.work_dir + "/samtools_index.o", "r") as r:
            data = r.readlines()
            for line in data:
                if re.match(r'.*index: failed to create index for.*Numerical result out of range', line):
                    self.small_bam = False
                    self.logger.info("genome is so large，will use csi!")

    def check_cram_bam(self):
        """
        解析samtools_index_csi.o，因为[E::hts_idx_push] NO_COOR reads not in a single block at the end 0 -1报错，
        所以在samtools index -c 需要将bam转为cram格式
        /mnt/lustre/users/sanger/app/miniconda2/bin/samtools view -C -T
        /mnt/lustre/users/sanger/app/database/dna_geneome/Ginkgo_biloba/GIGADB/HIC_V1/2019.06.07/ref.fa
        G83-1.sort.bam -o G83-1.sort.cram
        /mnt/lustre/users/sanger/app/miniconda2/bin/samtools index -c  G83-1.sort.cram
        :return:
        """
        with open(self.work_dir + "/samtools_index_csi.o", "r") as r:
            data = r.readlines()
            for line in data:
                if re.match(r'.*NO_COOR reads not in a single block at the end.*', line):
                    self.bamtosrcam = True
                    self.logger.info("need bam to cram，then make index of csi!")

    def run_bam2cram(self):
        """
        /mnt/lustre/users/sanger/app/miniconda2/bin/samtools view -C -T
        /mnt/lustre/users/sanger/app/database/dna_geneome/Ginkgo_biloba/GIGADB/HIC_V1/2019.06.07/ref.fa
        G83-1.sort.bam -o G83-1.sort.cram
        /mnt/lustre/users/sanger/app/miniconda2/bin/samtools index -c  G83-1.sort.cram
        注意该步骤后得到的G83.sort.bam其实不是bam文件，是转换后的cram文件，主要是为了不影响后续的分析模块进行调用，
        将cram文件命名为bam文件
        :return:
        """
        # 这里要排序下，目的是后面mappping stat的时候会报错,直接排序后输出cram文件，不要将bam转成cram再对cram排序。
        cmd2 = "{} sort -@ 8 --output-fmt CRAM --reference {} {} -o {}".format(self.samtools_path,
                                                                               self.option("ref_fa").prop["path"],
                                                                               self.work_dir + "/" + self.bam_name,
                                                                               self.output_dir + "/" + self.bam_name)
        command2 = self.add_command("samtools_sort_cram", cmd2).run()
        self.wait()
        if command2.return_code == 0:
            self.logger.info("samtools_sort_cram success")
        else:
            self.set_error("samtools_sort_cram failed")

        cmd1 = "{} index -c {}".format(self.samtools_path, self.output_dir + "/" + self.bam_name)
        command1 = self.add_command("samtools_index_csi_2", cmd1).run()
        self.wait()
        if command1.return_code == 0:
            self.logger.info("samtools_index_csi_2 success")
        else:
            self.set_error("samtools_index_csi_2 failed", code="34505104")

    def set_output(self):
        if not self.small_bam:
            return
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        os.link(self.option("sort_bam_file").prop["path"], self.output_dir + "/" + self.bam_name)
        if self.small_bam:
            os.link(self.work_dir + "/" + self.bam_name + ".bai", self.output_dir + "/" + self.bam_name + ".bai")
        else:
            pass
            # os.link(self.work_dir + "/" + self.bam_name + ".csi", self.output_dir + "/" + self.bam_name + ".csi")

    def run(self):
        super(SamtoolsIndexTool, self).run()
        self.run_samtools_index()
        self.set_output()
        self.end()
