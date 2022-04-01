## !/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "hongdong.xuan"
import re
import os
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class DedupV2Agent(Agent):
    """
    新的查重模块
    一个家系对应一个去重tool，输入的是母本tab，胎儿tab，错配位点，参考位点, 父本的id，以及要进行查重的父本列表。
    包括脚本：pt_dup_new.R
    version v1.0
    author: hongdong.xuan
    last_modify: 20171204
    """
    def __init__(self, parent):
        super(DedupV2Agent, self).__init__(parent)
        options = [
            {"name": "mom_tab", "type": "infile", "format": "paternity_test.tab"},
            {"name": "preg_tab", "type": "infile", "format": "paternity_test.tab"},
            {"name": "ref_point", "type": "infile", "format": "paternity_test.rda"},
            {"name": "err_min", "type": "int", "default": 2},
            {"name": "dad_id", "type": "string"},
            {"name": "father_path", "type": "string"},  # 输入父本tab文件的所在路径
            {"name": "dad_list", "type": "string"}  # 要进行分析的父本列表
        ]
        self.add_option(options)
        self.step.add_steps("dedup_analysis")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)
        self.queue = "gada"

    def stepstart(self):
        self.step.dedup_analysis.start()
        self.step.update()

    def stepfinish(self):
        self.step.dedup_analysis.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('ref_point'):
            raise OptionError("必须提供参考位点文件")
        if not self.option('mom_tab'):
            raise OptionError("必须提供母本tab")
        if not self.option('preg_tab'):
            raise OptionError("必须提供胎儿tab")
        if not self.option('father_path'):
            raise OptionError("必须提供查重部分父本tab")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = '40G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        super(DedupV2Agent, self).end()


class DedupV2Tool(Tool):
    """
    查重tool
    # use: pt_dup_new.R WQ123M.tab WQ123S.tab 2
    /mnt/ilustre/users/sanger-dev/sg-users/zhoumoli/pt/targets.bed.rda WQ123F WQ123-F1,WQ123-F2
    """
    def __init__(self, config):
        super(DedupV2Tool, self).__init__(config)
        self._version = '1.0.1'
        self.R_path = 'program/R-3.3.1/bin/'
        self.script_path = self.config.PACKAGE_DIR + '/medical/paternity_test_v3/'
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')

    def get_father_list(self):
        """
        获取dad_id前后200个父本tab文件，暂时不用
        :return:
        """
        dad_list = []
        results = os.listdir(self.option("father_path"))
        m = re.match(r'WQ([0-9]*)-F.*', self.option("dad_id"))
        if m:
            for dad in results:
                if int(m.group(1)) - 200 <= int(re.match(r'WQ([0-9]*)-F.*', str(dad)).group(1)) \
                        <= int(m.group(1)) + 201:
                    dad_list.append(dad)
        else:
            raise Exception("父本{}命名不规范！".format(self.option("dad_id")))
        return dad_list

    def run_tf(self):
        # dad_list = self.get_father_list()
        # self.logger.info(self.option("dad_list"))
        dedup_cmd = "{}Rscript {}pt_dup_new.R {} {} {} {} {} {} {} {}".format(self.R_path, self.script_path,
                                                                              self.option("mom_tab").prop['path'],
                                                                              self.option("preg_tab").prop['path'],
                                                                              self.option("err_min"),
                                                                              self.option("ref_point").prop['path'],
                                                                              "result", self.option("father_path"),
                                                                              self.option("dad_id"),
                                                                              self.option("dad_list"))
        self.logger.info(dedup_cmd)
        self.logger.info("开始进行查重分析")
        cmd = self.add_command("dedup_cmd", dedup_cmd).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("运行查重成功")
        else:
            self.set_error('运行查重出错')
            raise Exception("运行查重出错")

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        results = os.listdir(self.work_dir + '/result')
        for f in results:
            os.link(self.work_dir + '/result/' + f, self.output_dir + '/' + f)
        self.logger.info('设置文件夹路径成功')

    def run(self):
        super(DedupV2Tool, self).run()
        self.run_tf()
        self.set_output()
        self.end()
