# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'


import unittest
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class SamtoolsDepthAgent(Agent):
    def __init__(self, parent):
        super(SamtoolsDepthAgent, self).__init__(parent)
        options = [
            {"name": "in_sam", "type": "infile", "format": "align.bwa.sam"},  # sam file
            {'name': 'fq_type', 'type': 'string', 'default': None},
            {'name': 'sample_name', 'type': 'string', 'default': None},
            {'name': 'out', 'type': 'string', 'default': None},
        ]
        self.add_option(options)
        self.step.add_steps("samtools")
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.samtools.start()
        self.step.update()

    def step_finish(self):
        self.step.samtools.finish()
        self.step.update()

    def check_options(self):
        if not self.option("in_sam").is_set:
            raise OptionError("必须设置输入sam文件")
        if not self.option("fq_type"):
            raise OptionError("必须设置输入测序类型：PE/SE")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = "20G"

    def end(self):
        super(SamtoolsDepthAgent, self).end()


class SamtoolsDepthTool(Tool):
    def __init__(self, config):
        super(SamtoolsDepthTool, self).__init__(config)
        self._version = "v1.0"
        self.program = {
            'samtools': 'bioinfo/align/samtools-1.8/samtools',
            'minimap2': 'bioinfo/align/minimap2/minimap2-2.17_x64-linux/minimap2',
            'bwa': 'bioinfo/align/bwa-0.7.17/bwa',
        }
        self.file = {
            'bam': os.path.join(self.work_dir, self.option('sample_name') + '.bam'),
            'sorted_bam': os.path.join(self.work_dir, self.option('sample_name') + '_sorted.bam'),
            'outfile':  os.path.join(self.output_dir, self.option('sample_name') + '_map.xls')
        }

    def run(self):
        super(SamtoolsDepthTool, self).run()
        self.sam_view()
        self.bam_sort()
        with open(self.file['outfile'], 'w') as o:
            o.write('Reference\tPosition\tDepth\n')
        self.bam_depth()
        self.set_output()
        self.end()

    def sam_view(self):
        # 增加MAX_FILE_HANDLES_FOR_READ_ENDS_MAP参数， 原因是投递节点打开文件数目限制, 限制缓存
        self.logger.info(self.option('sample_name'))
        cmd = "{} view ".format(self.program['samtools'])
        if self.option('fq_type') == 'PE':
            cmd += '-F 12 '
        else:
            cmd += '-F 4 '
        cmd += '-bS {} -o {} --output-fmt BAM'.format(self.option('in_sam').prop['path'], self.file['bam'])
        self.logger.info("使用samtools转换sam文件")
        command = self.add_command("sam_view", cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("使用samtools转换sam文件完成!")
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("使用samtools转换sam文件出错！")

    def bam_sort(self):
        """

        """
        # 增加MAX_FILE_HANDLES_FOR_READ_ENDS_MAP参数， 原因是投递节点打开文件数目限制, 限制缓存
        self.logger.info(self.option('sample_name'))
        cmd = "{} sort {} -o {} --output-fmt BAM".format(self.program['samtools'], self.file['bam'],
                                                         self.file['sorted_bam'])
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

    def bam_depth(self):
        # 增加MAX_FILE_HANDLES_FOR_READ_ENDS_MAP参数， 原因是投递节点打开文件数目限制, 限制缓存
        self.logger.info(self.option('sample_name'))
        samtools_path = os.path.join(self.config.SOFTWARE_DIR, self.program['samtools'])
        cmd = "{} depth {} ".format(samtools_path,  self.file['sorted_bam'])
        cmd += '>> {}'.format(self.file['outfile'])
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


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "target_depth_samtools_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.target_depth.samtools_depth",
            "instant": False,
            "options": dict(
                sample_name='HY_1',
                in_sam='/mnt/ilustre/users/sanger-dev/workspace/20210511/Single_target_depth_minimap2_7979_7596/Minimap2/output/HY_1_map.sam',
                fq_type='SE',
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)