# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180423

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import datetime
import random
import os
import re


class Bam2fastqAgent(Agent):
    """
    数据组装接口--bam文件按照要求提取成fastq文件
    """
    def __init__(self, parent):
        super(Bam2fastqAgent, self).__init__(parent)
        options = [
            {"name": "bam_file", "type": "string"},     # bam文件
            {"name": "sample_id", "type": "string"},   # 页面传进来的样本id
            {"name": "pos", "type": "string"},   # 页面中的选择区域 chr1,chr1::100,000,chr1:1:,chr1:10000:20000
            {"name": "unmapping", "type": "string", "default": "true"}  # 页面选择添加的时候为true，否则为false
        ]
        self.add_option(options)
        self.step.add_steps('snpeff')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.snpeff.start()
        self.step.update()

    def step_end(self):
        self.step.snpeff.finish()
        self.step.update()
        
    def check_options(self):
        if not self.option("bam_file"):
            raise OptionError("缺少bam_file参数", code="34500401")
        if not self.option("sample_id"):
            raise OptionError("缺少sample_id参数", code="34500402")
        if not self.option("pos"):
            raise OptionError("pos", code="34500403")
        if not self.option("unmapping"):
            raise OptionError("缺少unmapping参数", code="34500404")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '20G'
        
    def end(self):
        super(Bam2fastqAgent, self).end()


class Bam2fastqTool(Tool):
    def __init__(self, config):
        super(Bam2fastqTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/align/samtools-1.7/')
        self.samtools = 'bioinfo/align/samtools-1.7/samtools '
        self.pos = ""
        
    def samtools_view(self):
        """
        samtools view -hb Lands.sort.bam chr1:10000-20000
        samtools fastq -c 9 -ns result_chr1:chr1:chr2_2/01.fastq/Lands.chr1:10000-20000.fastq.gz -
        :return:
        """
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + "_" + str(random.randint(1, 10000))
        script_path = self.config.SOFTWARE_DIR + '/bioinfo/WGS/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "script_{}.sh".format(now_time)
        cmd = "samtools view -hb {} {}".format(self.option("bam_file"), self.pos)
        # cmd += "|samtools fastq -c 9 -n -1 {}_{}.R1.fastq.gz -2 {}_{}.R2.fastq.gz - "\
        #     .format(self.option("sample_id"), self.pos, self.option("sample_id"), self.pos)
        cmd += "|samtools fastq -c 9 -ns {}.{}.fastq.gz - ".format(self.option("sample_id"), self.pos)
        self.logger.info("开始进行samtools_view_fastq")
        self.logger.info(cmd)
        with open(file_path, 'w') as w:
            w.write('#!/bin/bash' + "\n")
            w.write(cmd)
        code = os.system('/bin/chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改%s为可执行文件失败！",variables=(file_path), code="34500401")
        shell = "/bioinfo/WGS/script_temp/{}".format(os.path.basename(file_path))
        command1 = self.add_command("samtools_view", shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("脚本1完成！")
        else:
            self.set_error("脚本1出错！", code="34500402")
        os.system('rm {}'.format(file_path))

    def samtools_view_unmap(self):
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + "_" + str(random.randint(1, 10000))
        script_path = self.config.SOFTWARE_DIR + '/bioinfo/WGS/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "script_{}.sh".format(now_time)
        cmd = "samtools view -hbf 4 {}".format(self.option("bam_file"))
        # cmd += "|samtools fastq -c 9 -n -1 {}_unmap.R1.fastq.gz -2 {}_unmap.R2.fastq.gz - " \
        #     .format(self.option("sample_id"), self.option("sample_id"))
        cmd += "|samtools fastq -c 9 -ns {}.unmapping.fastq.gz - ".format(self.option("sample_id"))
        self.logger.info("开始进行samtools_view_fastq_unmapping")
        self.logger.info(cmd)
        with open(file_path, 'w') as w:
            w.write('#!/bin/bash' + "\n")
            w.write(cmd)
        code = os.system('chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改%s为可执行文件失败！",variables=(file_path), code="34500403")
        shell = "/bioinfo/WGS/script_temp/{}".format(os.path.basename(file_path))
        command1 = self.add_command("samtools_view_unmapping", shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("脚本2完成！")
        else:
            self.set_error("脚本2出错！", code="34500404")
        os.system('rm {}'.format(file_path))

    def set_pos_data(self):
        """
        chr1,chr1::100,000,chr1:1:,chr1:10000:20000
        前端默认传进来的格式chr1,2,3000, 若没有start则是chr1,,3000, 若没有end则是chr1,2,
        chr1:10000-20000
        :return:
        """
        temp = self.option("pos").split(",")
        if temp[0]:
            if temp[1] and temp[2]:
                self.pos = "{}:{}-{}".format(temp[0], temp[1], temp[2])
            elif temp[1] and not temp[2]:
                self.pos = "{}:{}".format(temp[0], temp[1])
            elif not temp[1] and temp[2]:
                self.pos = "{}:{}-{}".format(temp[0], "1", temp[2])
            else:
                self.pos = "{}".format(temp[0])
            self.logger.info(self.pos)
        else:
            self.set_error("%s格式不正确！必须为chr1,2,2000", variables=(self.option("pos")), code="34500405")

    def set_output(self):
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info('开始设置文件夹路径')
        results = os.listdir(self.work_dir)
        for f in results:
            if re.match(r'.*\.fastq\.gz$', f):
                os.link(os.path.join(self.work_dir, f), os.path.join(self.output_dir, f))
        self.logger.info('设置文件夹路径成功')

    def run(self):
        super(Bam2fastqTool, self).run()
        if self.option("pos") == ',,':  # pos为空的时候，添加unmapping数据必须为添加
            if self.option("unmapping") == "true":
                self.samtools_view_unmap()
            else:
                self.set_error("pos参数为空的时候，添加unmapping数据参数必须为添加", code="34500406")
        else:
            if self.option("unmapping") != "true":
                self.set_pos_data()
                self.samtools_view()
            else:
                self.set_pos_data()
                self.samtools_view()
                self.samtools_view_unmap()
        self.set_output()
        self.end()
