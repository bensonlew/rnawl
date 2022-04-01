# -*- coding: utf-8 -*-


import os
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
import shutil



class MetabolomeCleandataAgent(Agent):
    """
    代谢组通路注释
    """

    def __init__(self, parent):
        super(MetabolomeCleandataAgent, self).__init__(parent)
        options = [
            {"name" : "pos_table", "type":"infile","format":"tool_lab.simple"},
            {"name" : "neg_table", "type":"infile","format":"tool_lab.simple"},
            {"name" : "order_file", "type":"infile","format":"tool_lab.simple"},
            {"name": "column", "type": "string", "default":"id"},
            {"name": "qc_filter", "type": "float", "default":0.5},
            {"name": "sample_filter", "type": "float", "default":0.5},
            {"name": "method", "type": "string", "default":"svr"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('pos_table').is_set:
            raise OptionError('必须输入表格')
        if not self.option('order_file').is_set:
            raise OptionError('必须输入顺序表')
        return True

    def set_resource(self):
        """
        设置所需资源，需在子类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = "15G"

    def end(self):
        super(MetabolomeCleandataAgent, self).end()


class MetabolomeCleandataTool(Tool):
    def __init__(self, config):
        super(MetabolomeCleandataTool, self).__init__(config)
        self.python_path =  'program/Python/bin/python'
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/R-3.5.1/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')

    def run(self):
        """
        运行
        :return:
        """
        super(MetabolomeCleandataTool, self).run()
        if self.option('neg_table').is_set:
            cmd =self.python_path + ' ' + self.config.PACKAGE_DIR + "/tool_lab/meta_Norm.py -pos {0}  -neg {1}  -mark {2}  -order {3} -qc_f {4} -smp_f {5} -method {6} -top {7} -rsrc_path {8}".format(
                self.option('pos_table').path, self.option('neg_table').path, self.option('column'), self.option('order_file').path,
                self.option('qc_filter'), self.option('sample_filter'), self.option('method'), 5, self.config.PACKAGE_DIR
            )
        else:
            cmd =self.python_path + ' ' + self.config.PACKAGE_DIR + "/tool_lab/meta_Norm.py -pos {0}  -mark {1}  -order {2} -qc_f {3} -smp_f {4} -method {5} -top {6} -rsrc_path {7}".format(
                self.option('pos_table').path, self.option('column'), self.option('order_file').path,
                self.option('qc_filter'), self.option('sample_filter'), self.option('method'), 5, self.config.PACKAGE_DIR
            )
        self.logger.info(cmd)
        self.logger.info("开始运行")

        command = self.add_command("run_clean", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_clean完成")
        else:
            self.set_error("运行run_clean运行出错!")


        if self.option('neg_table').is_set:
            if os.path.exists(self.output_dir + '/neg'):
                shutil.rmtree(self.output_dir+'/neg')
            os.mkdir(self.output_dir + '/neg')
            f_name = os.path.basename(self.option('neg_table').path)
            os.link(self.work_dir+'/neg/Dealed_%s'%f_name, self.output_dir+'/neg/Dealed_%s'%f_name)
            shutil.copytree(self.work_dir+'/neg/svr_normalization_result', self.output_dir+'/neg/svr_normalization_result')

            if os.path.exists(self.output_dir + '/pos'):
                shutil.rmtree(self.output_dir+'/pos')
            os.mkdir(self.output_dir + '/pos')
            f_name = os.path.basename(self.option('pos_table').path)
            os.link(self.work_dir+'/pos/Dealed_%s'%f_name, self.output_dir+'/pos/Dealed_%s'%f_name)
            shutil.copytree(self.work_dir+'/pos/svr_normalization_result', self.output_dir+'/pos/svr_normalization_result')

        else:
            if os.path.exists(self.output_dir):
                shutil.rmtree(self.output_dir)
            os.mkdir(self.output_dir)
            #os.link(self.work_dir+'/neg/Dealed_pos_measurement.xls', self.output_dir+'/neg/Dealed_pos_measurement.xls')
            f_name = os.path.basename(self.option('pos_table').path)
            os.link(self.work_dir+'/pos/Dealed_%s'%f_name, self.output_dir+'/Dealed_%s'%f_name)
            shutil.copytree(self.work_dir+'/pos/svr_normalization_result', self.output_dir+'/svr_normalization_result')

        self.end()


