## !/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "moli.zhou"
#last_modify:20161121
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import os

class FamilyMergeAgent(Agent):
    """
    tab2family
    将已知的tab文件合并转化为整个家系在一起的表格
    包含family_joined.R
    version v1.0
    author: moli.zhou
    last_modify: 2016.11.21
    """
    def __init__(self, parent):
        super(FamilyMergeAgent, self).__init__(parent)
        options = [#输入的参数
            # {"name": "dad_tab", "type": "infile", "format": "tab"},
            # {"name": "mom_tab", "type": "infile", "format": "tab"},
            # {"name": "preg_tab", "type": "infile", "format": "tab"},
            # {"name": "err_min",  "type": "int", "default": 2},
            # {"name": "tab_merged", "type": "outfile", "format": "Rdata"}

            {"name": "dad_tab", "type": "infile", "format": "paternity_test.tab"},
            {"name": "mom_tab", "type": "infile", "format": "paternity_test.tab"},
            {"name": "preg_tab", "type": "infile", "format": "paternity_test.tab"},
            {"name": "ref_point", "type": "infile","format":"paternity_test.rda"},
            {"name": "err_min", "type": "int", "default": 2},
            {"name": "tab_merged", "type": "infile", "format": "paternity_test.rdata"}
        ]
        self.add_option(options)
        self.step.add_steps("Tab2family")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

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
        # if not self.option('query_amino'):
        #     raise OptionError("必须输入氨基酸序列")
        if self.option('err_min') < 1 :
            raise OptionError("err_min值超出范围，必须大于1")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 6
        self._memory = '50G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            ["family_merge_tab.Rdata", "Rdata", "家系合并后的表格"],
        ])
        super(FamilyMergeAgent, self).end()


class FamilyMergeTool(Tool):
    """

    """
    def __init__(self, config):
        super(FamilyMergeTool, self).__init__(config)
        self._version = '1.0.1'

        self.R_path = 'program/R-3.3.1/bin/'
        self.script_path = self.config.SOFTWARE_DIR + '/bioinfo/medical/scripts/'
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')

    def run_tf(self):
        tab2family_cmd = "{}Rscript {}family_joined.R {} {} {} {} {} {}".\
            format(self.R_path,self.script_path,self.option("dad_tab").prop['path'],
                   self.option("mom_tab").prop['path'],self.option("preg_tab").prop['path'],
                   self.option("err_min"),self.option("ref_point").prop['path'], self.work_dir)
        self.logger.info(tab2family_cmd)
        self.logger.info("开始运行家系合并")
        cmd = self.add_command("tab2family_cmd", tab2family_cmd).run()
        self.wait(cmd)

        self.logger.info("tab2family_cmd的返回码是{}".format(cmd.return_code))
        if cmd.return_code == 0:
            self.logger.info("运行家系合并成功")
        else:
            self.set_error("运行家系合并出错")
            raise Exception("运行家系合并出错")

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
            if re.search(r'.*family_joined_tab\.Rdata$',f):
                os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
            if re.search(r'.*family_joined_tab\.txt$', f):
                os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
        self.logger.info('设置文件夹路径成功')

    def run(self):
        super(FamilyMergeTool, self).run()
        self.run_tf()
        self.set_output()
        self.end()
