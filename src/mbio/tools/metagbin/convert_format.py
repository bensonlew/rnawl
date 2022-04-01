# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'@20190114
import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError

class ConvertFormatAgent(Agent):
    """
    宏基因组binning将sam文件转为fasta格式进行去冗余，将比对结果进行合并，并转为fastq和fasta文件
    """
    def __init__(self, parent):
        super(ConvertFormatAgent, self).__init__(parent)
        options = [
            {"name": "sam", "type": "infile", "format":"align.bwa.sam"}, #输入比对的文件
            {"name": "analysis", "type": "string", "default": ""}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sam"):
            raise OptionError("必须添加sam文件！", code="")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(ConvertFormatAgent, self).end()


class ConvertFormatTool(Tool):
    def __init__(self, config):
        super(ConvertFormatTool, self).__init__(config)
        self.samtools = "bioinfo/align/samtools-1.4/bin/samtools"
        #self.depth = self.work_dir + "/" + 'cov.depth.txt'
        self.sample_name = os.path.basename(self.option('sam').prop['path']).split('.')[0]

    def run_convert(self):
        """
        sam转bam
        :return:
        """
        input_sam = self.option('sam').prop['path']
        self.logger.info('开始生成bam文件')
        if not self.option("analysis") in ['metagbin']:
            out_bam = self.work_dir + "/" +self.sample_name + '_map.bam'
        else:
            out_bam = self.work_dir +"/"+ self.sample_name + '_before_map.bam'
        self.logger.info(out_bam)
        cmd_sam = "{} view -b -S {} -o {}".format(self.samtools, input_sam, out_bam)
        self.logger.info(cmd_sam)
        to_bam = 'sam_to_bam'
        command1 = self.add_command(to_bam, cmd_sam).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("%s运行完成"%to_bam)
        else:
            self.set_error("运行失败")

    def run_pick(self):
        """
        从比对结果bam文件中将比对上的序列提取
        :return:
        """
        self.input_bam = self.work_dir +"/"+ self.sample_name + '_before_map.bam'
        self.logger.info('开始生成bam文件')
        self.out_bam = self.work_dir +"/"+ self.sample_name + '_map.bam'
        self.logger.info(self.out_bam)
        cmd_bam = "{} view -F 4 -S {} -b -o {}".format(self.samtools, self.input_bam, self.out_bam)
        self.logger.info(cmd_bam)
        to_bam = 'to_bam'
        command1 = self.add_command(to_bam, cmd_bam).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("%s运行完成"%to_bam)
        else:
            self.set_error("运行失败")

    def run_sort(self):
        """
        对提取的文件进行排序
        :return:
        """
        self.out_bam = self.work_dir + "/" +self.sample_name + '_map.bam'
        self.sort_bam = self.output_dir + "/" +self.sample_name +'_sorted.bam'
        cmd_sort = '{} sort {} -o {}'.format(self.samtools, self.out_bam, self.sort_bam)
        self.logger.info(cmd_sort)
        to_sort = 'to_sort'
        command2 = self.add_command(to_sort, cmd_sort).run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("排序%s运行成功" %to_sort)
        else:
            self.set_error("运行失败")
        self.logger.info("排序运行结束")

    def run_index(self):
        """
        对bam文件提取index
        :return:
        """
        cmd_index = '{} index {}'.format(self.samtools, self.sort_bam)
        self.logger.info(cmd_index)
        to_index = 'to_index'
        command2 = self.add_command(to_index, cmd_index).run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("index%s运行成功" %to_index)
        else:
            self.set_error("运行失败")
        self.logger.info("index运行结束")

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        if self.option("analysis") in ['metagbin']:
            if os.path.exists(self.output_dir + "/" + self.sample_name + '_sorted.bam'):
                self.logger.info('sort结果文件生成成功')
            else:
                self.set_error('结果文件未生成')
        else:
            if os.path.exists(self.output_dir + "/" + self.sample_name + '_map.bam'):
                os.remove(self.output_dir + "/" + self.sample_name + '_map.bam')
            os.link(self.work_dir + "/" + self.sample_name + '_map.bam',self.output_dir + "/" + self.sample_name + '_map.bam')


    def run(self):
        super(ConvertFormatTool, self).run()
        if self.option("analysis") in ['metagbin']:
            self.run_convert()
            self.run_pick()
            self.run_sort()
            self.run_index()
            self.set_output()
            self.end()
        else:
            self.run_convert()
            self.set_output()
            self.end()
