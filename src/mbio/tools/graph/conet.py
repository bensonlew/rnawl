# -*- coding: utf-8 -*-
# __author__ = 'yuguo'

"""CoNet网络图分析工具"""

from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
import os
# import subprocess


class ConetAgent(Agent):
    """
    Conet
    version v3
    """
    def __init__(self, parent=None):
        """
        """
        super(ConetAgent, self).__init__(parent)
        options = [
            {"name": "data_file", "type": "infile", "format": "meta.otu.otu_table"},  # 输入数据矩阵
            {"name": "feature_file", "type": "infile", "format": "meta.env_table"},  # 输入环境特征文件
            {"name": "stand", "type": "string"},  # 数据标准化方法
            {"name": "method", "type": "string", "default": "correl_spearman"},  # Cooccurrence方法
            {"name": "lower_threshold", "type": "float", "default": -0.7},  # Cooccurrence 阈值，最小值
            {"name": "upper_threshold", "type": "float", "default": 0.7},  # Cooccurrence 阈值，最大值
            {"name": "network_file", "type": "outfile", "format": "graph.gml"},  # 输出网络图文件
            {"name": "randomization", "type": "bool", "default": True},  # 是否进行网络图随机化计算
            {"name": "iterations", "type": "int", "default": 100},  # 随机化迭代次数
            {"name": "resamplemethod", "type": "string", "default": "permute"},  # 重抽样方法
            {"name": "pval_threshold", "type": "float", "default": 0.05}  # 重抽样方法
            ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数设置
        """
        if not self.option("data_file").is_set:
            raise OptionError("必须设置输入数据矩阵")
        if self.option("stand"):
            if self.option("stand") not in ("row_stand", "row_stand_robust", "row_norm", "row_downsample", "row_norm_external", "col_norm", "col_downsample", "col_norm_external", "log2"):
                raise OptionError("数据标准化方法不在范围内")
        if self.option("resamplemethod") not in ("permute", "bootstrap"):
            raise OptionError("重抽样方法只能选择'permute', 'bootstrap'")
        else:
            if self.option("resamplemethod") == "permute":
                self.option("resamplemethod", value="shuffle_rows")

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        self._memory = 6000


class ConetTool(Tool):
    def __init__(self, config):
        super(ConetTool, self).__init__(config)
        self.script_path = self.config.SOFTWARE_DIR+"/meta/CoNet3/lib/CoNet.jar"
        self.java_path = "/sun_jdk1.8.0/bin/"

    def run(self):
        super(ConetTool, self).run()
        self.run_conet()

    def run_conet(self):
        cmd = self.java_path+"/java -Xmx6000m -cp "+self.script_path+" be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser "
        cmd += "--method ensemble --format GML --matrixtype abundance --output network.gml"
        cmd += " --input " + self.option("data_file").prop['path']
        cmd += " --ensemblemethods "+self.option("method")
        cmd += " --ensembleparams "+self.option("method")+"~upperThreshold="+str(self.option("upper_threshold"))+"/"+self.option("method")+"~lowerThreshold="+str(self.option("lower_threshold"))
        if self.option("feature_file").is_set:
            cmd += " --features "+self.option("feature_file")
        if self.option("randomization"):
            cmd += " --resamplemethod "+self.option("resamplemethod")+" --iterations "+str(self.option("iterations"))+" --edgethreshold "+str(self.option("pval_threshold"))
        self.logger.info(u"生成命令: "+cmd)
        conet = self.add_command("conet", cmd)
        self.logger.info("开始运行conet")
        conet.run()
        self.wait(conet)
        if conet.return_code == 0:
            self.logger.info("conet运行完成")
            os.link(self.work_dir+'/network.gml', self.output_dir+'/network.gml')
            self.option('network_file').set_path(self.output_dir+'/network.gml')
            self.end()
        else:
            self.set_error("conet运行出错!")
