# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'


import unittest
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class PositionDepthAgent(Agent):
    def __init__(self, parent):
        super(PositionDepthAgent, self).__init__(parent)
        options = [
            {"name": "in_bam", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "in_bed", "type": "infile", "format": "ref_rna_v2.bed"},
        ]
        self.add_option(options)
        self.step.add_steps("position_depth")
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.position_depth.start()
        self.step.update()

    def step_finish(self):
        self.step.position_depth.finish()
        self.step.update()

    def check_options(self):
        if not self.option("in_bam").is_set:
            raise OptionError("必须设置输入BAM文件")
        if not self.option("in_bed").is_set:
            raise OptionError("必须设置输入BED文件")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = "20G"

    def end(self):
        super(PositionDepthAgent, self).end()


class PositionDepthTool(Tool):
    def __init__(self, config):
        super(PositionDepthTool, self).__init__(config)
        self._version = "v1.0"
        self.sample_name = os.path.splitext(os.path.basename(self.option("in_bam").prop["path"]))[0]
        self.program = {
            'samtools': 'bioinfo/align/samtools-1.8/samtools',
        }
        self.file = {
            'outfile':  os.path.join(self.output_dir, self.sample_name + '_depth.xls')
        }

    def run(self):
        super(PositionDepthTool, self).run()
        self.bam_sort()
        self.bam_index()
        self.bam_depth()
        self.set_output()
        self.end()

    def bam_sort(self):
        """

        """
        # 增加MAX_FILE_HANDLES_FOR_READ_ENDS_MAP参数， 原因是投递节点打开文件数目限制, 限制缓存
        self.logger.info(self.sample_name)
        cmd = "{} sort {} -o {} --output-fmt BAM".format(self.program['samtools'], self.option("in_bam").prop["path"],
                                                         'sorted_' + self.sample_name + '.bam')
        self.logger.info("使用samtools对bam文件进行排序")
        command = self.add_command("bam_sort", cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("使用samtools对bam文件进行排序完成!")
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("使用samtools对bam文件进行排序出错！")

    def bam_index(self):
        # 增加MAX_FILE_HANDLES_FOR_READ_ENDS_MAP参数， 原因是投递节点打开文件数目限制, 限制缓存
        self.logger.info(self.sample_name)
        if os.path.isfile(os.path.join(self.work_dir, '{}.bam.bai'.format(self.sample_name))):
            return
        cmd = "{} index {}".format(self.program['samtools'],  'sorted_' + self.sample_name + '.bam')
        self.logger.info("使用samtools对bam文件建立索引")
        command = self.add_command("bam_index", cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("使用samtools对bam文件建立索引完成!")
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("使用samtools对bam文件建立索引出错！")

    def bam_depth(self):
        # 增加MAX_FILE_HANDLES_FOR_READ_ENDS_MAP参数， 原因是投递节点打开文件数目限制, 限制缓存
        self.logger.info(self.sample_name)
        samtools_path = os.path.join(self.config.SOFTWARE_DIR, self.program['samtools'])
        cmd = "{} depth {} -b {} ".format(samtools_path,  'sorted_' + self.sample_name + '.bam',
                                          self.option('in_bed').prop['path'])
        cmd += '> {}'.format(self.file['outfile'])
        self.logger.info("使用samtools统计bam文件位点深度")
        command = self.add_command("bam_depth", cmd, ignore_error=True, shell=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("使用samtools统计bam文件位点深度完成!")
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("使用samtools统计bam文件位点深度出错！")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        pass
