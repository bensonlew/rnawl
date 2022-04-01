# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import os
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool

class NewblerAgent(Agent):
    """
    进行newbler拼接
    version: v1.0
    author: guhaidong
    last_modify: 2017.09.13
    """
    def __init__(self, parent):
        super(NewblerAgent, self).__init__(parent)
        options = [
            {"name": "contig", "type": "infile", "format": "sequence.fasta"},  # 输入fasta文件
            {"name": "cpu", "type": "int", "default": 5},  # 拼接线程数，默认5
            {"name": "mem", "type": "int", "default": 10},  # 拼接使用内存，默认10
            {"name": "mi", "type": "int", "default": 98},  # 拼接相似度0-100，默认98
            {"name": "ml", "type": "int", "default": 40},  # 拼接比对长度，默认40
            {"name": "all_length","type": "int", "default": 300},  # 拼接结果最小contig长度
            {"name": "large_length","type": "int", "default": 1000},  # 拼接结果认为是长contig的长度
            #{"name": "output", "type": "string"},  # 输出拼接结果状态文件路径
        ]
        self.add_option(options)
        self.step.add_steps("Newbler")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)
        self._memory_increase_step = 50  # 每次重运行增加内存50G by GHD @ 20180502

    def stepstart(self):
        self.step.Newbler.start()
        self.step.update()

    def stepfinish(self):
        self.step.Newbler.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('contig'):
            raise OptionError('必须输入拼接结果文件', code="31301301")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = self.option('cpu')
        self._memory = "{}G".format(self.option('mem'))

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(NewblerAgent, self).end()


class NewblerTool(Tool):
    def __init__(self, config):
        super(NewblerTool, self).__init__(config)
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.newbler_path = '/bioinfo/metaGenomic/454-newbler/apps/mapper/bin/'

    def newbler_run(self):
        """
        进行newbler拼接
        :return:
        """
        if os.path.exists(self.output_dir + '/454ReadStatus.txt'):
            return
        cmd = self.newbler_path + 'runAssembly -o %s -force -cpu %s -mi %s -ml %s -a %s -l %s %s'\
           % (self.work_dir,
           self.option('cpu'),
           self.option('mi'),
           self.option('ml'),
           self.option('all_length'),
           self.option('large_length'),
           self.option('contig').prop['path'])
        self.logger.info("运行newbler拼接")
        command = self.add_command("newbler", cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行newbler完成")
        elif command.return_code == 255:  # newbler拼接的线程问题 @ 20180328
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'thread error')  # 借用memory_limit状态，将tool重新投递到其他节点 @ 20180330
        else:
            self.logger.info("return code: %s" % command.return_code)
            self.set_error("newbler运行出错!", code="31301301")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        out_fa = self.output_dir + '/newbler.contig.fa'
        out_stat = self.output_dir + '/454ReadStatus.txt'
        if os.path.exists(out_fa):
            os.remove(out_fa)
        if os.path.exists(out_stat):
            os.remove(out_stat)
        os.link(self.work_dir + '/454AllContigs.fna', out_fa)
        os.link(self.work_dir + '/454ReadStatus.txt', out_stat)
        #self.option('output', out_stat)
        #self.logger.info(self.option('output'))
        self.logger.info("设置newbler分析结果目录成功")

    def run(self):
        """
        运行
        :return:
        """
        super(NewblerTool, self).run()
        self.newbler_run()
        self.set_output()
        self.end()