# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.08

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import datetime
import random


class SamtoolsDepthAgent(Agent):
    """
    软件: samtools
    samtools的depth方法
    """
    def __init__(self, parent):
        super(SamtoolsDepthAgent, self).__init__(parent)
        options = [
            {"name": "bam_file", "type": "infile", "format": "align.bwa.bam"},  # bam文件
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "large_genome", "type": "string", "default": "false"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("bam_file").is_set:
            raise OptionError("请设置bam文件", code="34504901")

    def set_resource(self):
        self._cpu = 9
        self._memory = "30G"

    def end(self):
        super(SamtoolsDepthAgent, self).end()


class SamtoolsDepthTool(Tool):
    def __init__(self, config):
        super(SamtoolsDepthTool, self).__init__(config)
        self.samtools_sh_path = "bioinfo/WGS/samtools_depth.sh"
        self.samtools_path = self.config.SOFTWARE_DIR + "/bioinfo/align/samtools-1.10/samtools"

    def run_samtools_depth(self):
        """
        samtools depth
        """
        sample_name = os.path.basename(self.option("bam_file").prop["path"]).split(".")[0]
        depth_file = sample_name + ".coverage"
        cmd = "{} {} {} {}".format(self.samtools_sh_path, self.samtools_path, self.option("bam_file").prop["path"],
                                   self.work_dir + "/" + depth_file)
        command = self.add_command("samtools_depth", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("samtools depth完成")
        else:
            self.set_error("samtools depth失败", code="34504901")
        if os.path.exists(self.output_dir + "/" + depth_file):
            os.remove(self.output_dir + "/" + depth_file)
        os.link(self.work_dir + "/" + depth_file, self.output_dir + "/" + depth_file)

    def run_samtools_depth_cram(self):
        """
        run samtools depth cram
        /mnt/lustre/users/sanger/app/miniconda2/bin/samtools depth
        --reference /mnt/lustre/users/sanger/app/database/dna_geneome/Ginkgo_biloba/GIGADB/HIC_V1/2019.06.07/ref.fa
        ./G83.sort.cram  > /mnt/lustre/users/sanger/workspace/20200320/WgsV2_majorbio_245724/MappingStat/
        SamtoolsDepth1/G83.coverage
        :return:
        """
        sample_name = os.path.basename(self.option("bam_file").prop["path"]).split(".")[0]
        depth_file = sample_name + ".coverage"
        cmd1 = "{} depth --reference {} {} > {}".format(self.samtools_path, self.option("ref_fa").prop["path"],
                                                        self.option("bam_file").prop["path"],
                                                        self.work_dir + "/" + depth_file)
        self.logger.info("cmd:{}".format(cmd1))
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + "_" + str(random.randint(1, 10000))
        script_path = self.config.SOFTWARE_DIR + '/bioinfo/WGS/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "samtools_depth_{}.sh".format(now_time)
        with open(file_path, 'w') as w:
            w.write('#!/bin/bash' + "\n")
            w.write(cmd1)
        code = os.system('/bin/chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改%s为可执行文件失败！")
        shell = "/bioinfo/WGS/script_temp/{}".format(os.path.basename(file_path))
        self.logger.info("开始进行run_samtools_depth")
        command1 = self.add_command("run_samtools_depth", shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("run_samtools_depth完成！")
        else:
            self.set_error("run_samtools_depth出错！")
        os.system('rm {}'.format(file_path))
        if os.path.exists(self.output_dir + "/" + depth_file):
            os.remove(self.output_dir + "/" + depth_file)
        os.link(self.work_dir + "/" + depth_file, self.output_dir + "/" + depth_file)

    def run(self):
        super(SamtoolsDepthTool, self).run()
        if self.option("large_genome") == 'false':
            self.run_samtools_depth()
        else:
            self.run_samtools_depth_cram()
        self.end()
