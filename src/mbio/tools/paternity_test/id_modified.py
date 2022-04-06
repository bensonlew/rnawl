## !/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/python
# -*- coding: utf-8 -*-
# __author__ = "moli.zhou"
#last_modify:20161205

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import os
import re

class IdModifiedAgent(Agent):
    """
    id.py
    将重送样或者二次上机的样本名进行改动
    tab文件中的FC改为F，但存入数据库的tab文件名不变
    """
    def __init__(self, parent):
        super(IdModifiedAgent, self).__init__(parent)
        options = [#输入的参数
            {"name": "sample_id", "type": "string"},  # 输入F/M/S的样本ID
            {"name": "fastq_path", "type": "string"},  # fastq所在路径
        ]
        self.add_option(options)
        self.step.add_steps("id_modified")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.id_modified.start()
        self.step.update()

    def stepfinish(self):
        self.step.id_modified.finish()
        self.step.update()


    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("sample_id"):
            raise OptionError("必须输入样本编号")
        if not self.option("fastq_path"):
            raise OptionError("必须输入fastq的路径")

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        self._memory = '100G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            [".mem.sort.hit.vcf.id.tab", "tab", "修改重送样的样本名"],
        ])
        super(IdModifiedAgent, self).end()


class IdModifiedTool(Tool):
    """

    """
    def __init__(self, config):
        super(IdModifiedTool, self).__init__(config)
        self._version = '1.0.1'

        self.python_path = 'program/Python/bin/'
        self.script_path = Config().SOFTWARE_DIR + '/bioinfo/medical/scripts/'

    def modified_run(self):
        idmodified_cmd = "{}python {}id.py {} {}".\
            format(self.python_path,self.script_path,
                   self.option("fastq_path"), self.option("sample_id"))
        self.logger.info(idmodified_cmd)
        self.logger.info("开始运行")
        cmd = self.add_command("tab2family_cmd", idmodified_cmd).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("运行成功")
        else:
            self.logger.info("运行出错")

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        results = os.listdir(self.work_dir)
        for f in results:
            if re.search(r'.*.id\.tab$', f):
                os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
        self.logger.info('设置文件夹路径成功')

    def run(self):
        super(IdModifiedTool, self).run()
        self.modified_run()
        self.set_output()
        self.end()
