## !/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "hongdongxuan"

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import os
import re


class CorrNetworkCalcAgent(Agent):
    """
    调用corr_net_calc.py，计算物种相似性网络的拓扑属性
    version v1.0
    author: hongdongxuan
    last_modify: 2016.12.07
    """
    def __init__(self, parent):
        super(CorrNetworkCalcAgent, self).__init__(parent)
        options = [
            {"name": "corr_table", "type": "string"},  # 物种相似性表
            {"name": "coefficient", "type": "float", "default": 0.08},  # 相关性系数阈值
            {"name": "significance", "type":"float", "default": 0.05}
        ]
        self.add_option(options)
        self.step.add_steps("CorrNetworkCalc")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.CorrNetworkCalc.start()
        self.step.update()

    def stepfinish(self):
        self.step.CorrNetworkCalc.finish()
        self.step.update()


    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("corr_table"):
            raise OptionError("必须输入物种相似性网络文件", code="32704601")
        if self.option('coefficient') <= 0 or self.option('coefficient') >= 1:
            raise OptionError('coefficient值必须是在0~1范围之内的值', code="32704602")
        if self.option('significance') < 0 or self.option('significance') > 1:
            raise OptionError('PARAMETERS ERROR: wrong value of significance (%s), [0,1] expected!' % self.option('significance'))
        return True

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
            [".", "", "结果输出目录"],
            ["corr_network_attributes.txt", "txt", "网络的单值属性表"],
            ["corr_network_by_cut.txt", "txt", "相关系数筛选后网络边文件"],
            ["corr_network_centrality.txt", "txt", "网络节点的中心系数表"],
            ["corr_network_clustering.txt", "txt", "网络节点的聚类系数表"],
            ["corr_network_degree_distribution.txt", "txt", "网络节点的度分布表"],
            ["corr_network_node_degree.txt", "txt", "网络节点的度统计表"]
        ])

        super(CorrNetworkCalcAgent, self).end()


class CorrNetworkCalcTool(Tool):
    """
    将otu转为shared.txt，使用mothur就算物种间相似性
    python corr_net_calc.py -i shared1.0.03.spearman.otu.corr -c 0.08 -o corr_result
    """
    def __init__(self, config):
        super(CorrNetworkCalcTool, self).__init__(config)
        self._version = '1.0.1'
        self.python_path = 'program/Python/bin/python '
        self.script_path = os.path.join(Config().SOFTWARE_DIR, "bioinfo/meta/scripts/corr_net_calc.py")

    def run_CorrNetworkCalc(self):
        one_cmd = self.python_path + self.script_path + " -i %s -c %s -sig %s -o %s " % (self.option('corr_table'), self.option('coefficient'), self.option('significance'), "corr_result")
        self.logger.info(one_cmd)
        self.logger.info("开始运行networkcalc")
        cmd = self.add_command("networkcalc", one_cmd).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("运行networkcalc成功")
        elif cmd.return_code is None:
            self.logger.info("返回码为None，重新运行一次")
            re_cmd = self.add_command("networkcalc", one_cmd).rerun()
            if re_cmd.return_code == 0:
                self.logger.info("重新运行一次成功！")
            else:
                self.logger.info("运行networkcalc出错")
        #     return True
        # else:
        #     self.set_error("运行networkcalc出错")
        #     return False

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        results = os.listdir(self.work_dir + '/' + "corr_result/")
        for f in results:
            os.link(self.work_dir + '/' + "corr_result/" + f, self.output_dir + "/" + f)
        self.logger.info('设置文件夹路径成功')


    def run(self):
        super(CorrNetworkCalcTool, self).run()
        self.run_CorrNetworkCalc()
        self.set_output()
        self.end()
