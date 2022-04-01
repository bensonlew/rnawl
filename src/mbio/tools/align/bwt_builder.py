# -*- coding: utf-8 -*-
# __author__ == zouxuan

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess
import shutil


class BwtBuilderAgent(Agent):
    """
    create index by 2bwt-builder
    author: zouxuan
    modified at date 20170911
    """

    def __init__(self, parent):
        super(BwtBuilderAgent, self).__init__(parent)
        options = [
            {"name": "fafile", "type": "infile", "format": "sequence.fasta"},  # 非冗余基因集fasta文件
            {"name": "build_dir", "type": "outfile", "format": "align.bwt_index_dir"}  # 输出的索引文件夹
        ]
        self.add_option(options)
        self.step.add_steps('bwtbuilder')
        self._memory_increase_step = 20  # 每次重运行增加内存20G by guhaidong @ 20180427
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.bwtbuilder.start()
        self.step.update()

    def step_end(self):
        self.step.bwtbuilder.finish()
        self.step.update()

    def check_options(self):
        if not self.option("fafile").is_set:
            raise OptionError("必须提供非冗余基因集", code="31100801")

    def set_resource(self):
        self._cpu = 2
        if os.path.getsize(self.option("fafile").prop['path'])/500000000 < 5:
            self._memory = '25G'  # 改回 by GHD @ 20180502
            # self._memory = 25
        else:
            self._memory = str(os.path.getsize(self.option("fafile").prop['path'])/500000000 + 40) + 'G'  # 改回 by GHD @ 20180502
            # self._memory = os.path.getsize(self.option("fafile").prop['path'])/500000000 + 40
        # self._memory = "%sG" % (self._memory + 20 * self._rerun_time)  # 每次因拼接失败而重运行的内存增加10G by GHD @ 20180323
        # self.logger.info('bwt_builder use memory : ' + self._memory)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        super(BwtBuilderAgent, self).end()


class BwtBuilderTool(Tool):
    def __init__(self, config):
        super(BwtBuilderTool, self).__init__(config)
        self._version = "1.0"
        self.bwt_path = "bioinfo/uniGene/soap2.21release/2bwt-builder"
        self.out_dir = os.path.join(self.work_dir, '2bwt_index')
        self.fafile_name = os.path.basename(self.option('fafile').prop['path'])

    def run(self):
        super(BwtBuilderTool, self).run()
        self.makefile()
        self.build_run()
        self.set_output()

    def makefile(self):
        if os.path.exists(self.out_dir):
            pass
        else:
            os.mkdir(self.out_dir)
        shutil.copyfile(self.option('fafile').prop['path'], os.path.join(self.out_dir, self.fafile_name))

    def build_run(self):
        cmd = '%s %s' % (self.bwt_path, os.path.join(self.out_dir, self.fafile_name))
        self.logger.info("运行2bwt-builder生成索引文件")
        self.logger.info(cmd)
        command1 = self.add_command("build", cmd, ignore_error=True)
        command1.run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("2bwt-builder succeed")
        elif command1.return_code in [-7, -9, -11]:  # 加入return_code检测，在sanger超出内存的返回值为-9, add -11 @20190124
            self.logger.info("return code: %s" % command1.return_code)
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180323
        else:
            self.logger.info("return code: %s" % command1.return_code)
            self.set_error("2bwt-builder failed", code="31100801")
            # raise Exception("2bwt-builder failed")

    def set_output(self):
        self.logger.info("set output")
        if os.path.exists(os.path.join(self.output_dir, "2bwt_index")):
            shutil.rmtree(os.path.join(self.output_dir, "2bwt_index"))
        shutil.copytree(self.out_dir, os.path.join(self.output_dir, "2bwt_index"))
        self.option("build_dir").set_path(self.output_dir + '/2bwt_index')
        self.end()
