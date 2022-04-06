## !/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/python
# -*- coding: utf-8 -*-
# __author__ = "hongdong"
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError

import os


class MsMatchAgent(Agent):
    """
    该模块用于检测MS是否匹配
    version v1.0
    author: HONGDONG
    last_modify: 20180821
    """
    def __init__(self, parent):
        super(MsMatchAgent, self).__init__(parent)
        options = [
            {"name": "dad_tab", "type": "infile", "format": "paternity_test.tab"},
            {"name": "mom_tab", "type": "infile", "format": "paternity_test.tab"},
            {"name": "preg_tab", "type": "infile", "format": "paternity_test.tab"},
            {"name": "ref_point", "type": "infile", "format": "paternity_test.rda"},
            {"name": "err_min", "type": "int", "default": 5}
        ]
        self.add_option(options)
        self.step.add_steps("Tab2family")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)
        self.queue = "gada"

    def stepstart(self):
        self.step.Tab2family.start()
        self.step.update()

    def stepfinish(self):
        self.step.Tab2family.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if self.option('err_min') < 1:
            raise OptionError("err_min值超出范围，必须大于1")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(MsMatchAgent, self).end()


class MsMatchTool(Tool):
    def __init__(self, config):
        super(MsMatchTool, self).__init__(config)
        self._version = '1.0.1'
        self.R_path = 'program/R-3.3.1/bin/'
        self.script_path = self.config.PACKAGE_DIR + '/medical/paternity_test_v3/'
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')

    def run_tf(self):
        """
        /mnt/ilustre/users/sanger-dev/app/program/R-3.3.1/bin/Rscript
        /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/medical/paternity_test_v3/ms_match.R
        /mnt/ilustre/users/sanger-dev/app/database/human/pt_ref/tab_all/WQ181495-F1.tab
        /mnt/ilustre/users/sanger-dev/app/database/human/pt_ref/tab_all/WQ181495-M-1.tab
        /mnt/ilustre/users/sanger-dev/app/database/human/pt_ref/tab_all/WQ181495-S-1.tab
        2 /mnt/ilustre/users/sanger-dev/app/database/human/pt_ref/targets.bed.rda
        :return:
        """
        tab2family_cmd = "{}Rscript {}ms_match.R {} {} {} {} {}".\
            format(self.R_path, self.script_path, self.option("dad_tab").prop['path'],
                   self.option("mom_tab").prop['path'], self.option("preg_tab").prop['path'],
                   self.option("err_min"), self.option("ref_point").prop['path'])
        self.logger.info(tab2family_cmd)
        self.logger.info("开始运行家系合并")
        cmd = self.add_command("tab2family_cmd", tab2family_cmd).run()
        self.wait(cmd)
        self.logger.info("tab2family_cmd的返回码是{}".format(cmd.return_code))
        if cmd.return_code == 0:
            self.logger.info("运行家系合并成功")
        else:
            self.set_error("运行家系合并出错")

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
            if re.search(r'.*info_show\.txt$', f):
                os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
        self.logger.info('设置文件夹路径成功')

    def run(self):
        super(MsMatchTool, self).run()
        self.run_tf()
        self.set_output()
        self.end()
