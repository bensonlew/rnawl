## !/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/python
# -*- coding: utf-8 -*-
# __author__ = "hongdongxuan"
#last_modify:20161128

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import os
import re
import shutil


class Family2bamAgent(Agent):
    """
    调用fastq2bam.sh脚本，完成无创亲子鉴定的生信分析流程中将fastq转为bam文件
    version v1.0
    author: hongdongxuan
    modified: moli.zhou
    last_modify: 2016.11.28
    """
    def __init__(self, parent):
        super(Family2bamAgent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "string"},  #输入F/M/S的fastq文件的样本名,fastq_gz_dir/WQ235F
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"}, #hg38.chromosomal_assembly/ref.fa
            {"name": "targets_bedfile","type": "infile","format":"paternity_test.rda"},
            {"name": "seq_path", "type": "infile","format":"sequence.fastq_dir"}, #fastq所在路径
            {"name": "cpu_number", "type": "int", "default": 4}
        ]
        self.add_option(options)
        self.step.add_steps("family2bam_tool")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.family2bam_tool.start()
        self.step.update()

    def stepfinish(self):
        self.step.family2bam_tool.finish()
        self.step.update()


    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("fastq"):
            raise OptionError("必须输入fastq文件的样本名")
        if not self.option("ref_fasta").is_set:
            raise OptionError("必须输入参考基因组的fastq文件")
        if not self.option('seq_path'):
            raise OptionError('必须提供fastq文件所在的路径')
        if not self.option('targets_bedfile'):
            raise OptionError('必须提供target_bedfile文件')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        self._memory = '50G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
                    ])
        result_dir.add_regexp_rules([
            [r".mem.sort.hit.filter.bam", "bam", "所有位点的信息"],
            [r".qc", "qc", "数据质控文件"]

        ])
        super(Family2bamAgent, self).end()


class Family2bamTool(Tool):
    """
    运行fastq2bam.sh sample_id cpu_num ref fastq_dir targets_bedfile
    """
    def __init__(self, config):
        super(Family2bamTool, self).__init__(config)
        self._version = '1.0.1'
        self.cmd_path = "bioinfo/medical/scripts/"
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/ruby-2.3.1') # 测试机
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/ruby-2.4.1/bin')  # 正式机
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/seq/bioawk')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/seq/seqtk-master')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/align/bwa-0.7.15')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/seq/samblaster-0.1.24')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/align/samtools-1.4/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/medical/bedtools-2.24.0/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/seq/bcftools-1.4/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/medical/picard-tools-2.2.4')
        self.picard_path = self.config.SOFTWARE_DIR + '/bioinfo/medical/picard-tools-2.2.4/picard.jar'
        self.java_path = self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin/java'



    def run_Family2bam(self):
        fastq2bam_cmd = "{}fastq2bam.sh {} {} {} {} {} {} {}".format(self.cmd_path, self.option("fastq"), self.option("cpu_number"),
                            self.option("ref_fasta").prop["path"], self.option("seq_path").prop['path'], self.option("targets_bedfile").prop['path']
                            ,self.picard_path,self.java_path)
        # fastq2bam_cmd = "{}fastq2bam.sh {} {} {} {} {}".format(self.cmd_path, self.option("fastq"), self.option("cpu_number"),
        #                     self.option("ref_fasta").prop["path"], self.option("seq_path").prop['path'], self.option("targets_bedfile").prop['path'])

        self.logger.info(fastq2bam_cmd)
        self.logger.info("开始运行转bam文件")
        cmd = self.add_command("fastq2bam_cmd", fastq2bam_cmd).run()
        self.wait(cmd)

        if cmd.return_code == 0:
            self.logger.info("运行转bam文件成功")
        else:
            self.set_error('运行转bam文件出错')
            raise Exception("运行转bam文件出错")

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        file_path = self.option("seq_path").prop['path'] + "/"
        print file_path
        results = os.listdir(file_path)
        for f in results:
            if str(f) == "%s.mem.sort.hit.filter.bam"%(self.option("fastq")) or str(f) == "%s.qc"%(self.option("fastq")):
            #if re.search(r'.*\.filter\.bam$', f) or re.search(r'.*\.qc$', f):
                shutil.move(file_path + f, self.output_dir)
        self.logger.info('设置文件夹路径成功')

    def run(self):
        super(Family2bamTool, self).run()
        self.run_Family2bam()
        self.set_output()
        self.end()
