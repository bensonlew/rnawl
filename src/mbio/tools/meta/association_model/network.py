# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess
import os
import re


class NetworkAgent(Agent):
    """
    宏基因组分布网络计算
    author: shaohua.yuan
    last_modify:
    """

    def __init__(self, parent):
        super(NetworkAgent, self).__init__(parent)
        options = [
            {"name": "profile_table", "type": "infile", "format": "sequence.profile_table"},
            #{"name": "grouptable", "type": "infile", "format": "meta.otu.group_table"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("profile_table").is_set:
            raise OptionError("必须设置输入丰度文件", code="32701301")
        return True

    def set_resource(self):
        self._cpu = 10
        self._memory = '15G'

    def end(self):
        super(NetworkAgent, self).end()


class NetworkTool(Tool):
    def __init__(self, config):
        super(NetworkTool, self).__init__(config)
        self.python_path = "/miniconda2/bin/python"
        self.cmd_path = self.config.SOFTWARE_DIR + '/bioinfo/meta/scripts/calc_otu_network.py'

    def run(self):
        """
        运行
        :return:
        """
        super(NetworkTool, self).run()
        self.run_network()
        self.set_output()
        self.end()

    def run_network(self):
        table = self.option("profile_table").prop["path"]
        cmd = '{} {} -i {} -o {}'.format(self.python_path, self.cmd_path, table, self.output_dir)
        self.logger.info('开始生成网络并进行计算')
        command = self.add_command("network", cmd).run()
        self.logger.info(cmd)
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行Network完成")
        else:
            self.set_error("运行Network出错", code="32701301")
            raise Exception("运行Network出错")

    def set_output(self):
        self.logger.info("修改文件")
        try:
            with open(self.output_dir + "/real_node_table.txt", "rb") as f1, open(\
                            self.output_dir + "/network_nodes_table.txt", "w") as f2:
                for line in f1:
                    line = line.strip()
                    line1 = line.split("\t")
                    if "user_node" in line1:
                        line = line.replace("user_node", "sample_node")
                    elif "otu_node" in line1:
                        line = line.replace("otu_node", "tax_fun_node")
                    f2.write(line + "\n")
            open(self.output_dir + "/network_nodes_degree.txt", 'w').write(
                re.sub(r'OTU', 'Tax_Fun', open(self.output_dir + '/real_dc_sample_otu_degree.txt').read()))
            open(self.output_dir + "/network_tax_fun_degree.txt", 'w').write(
                re.sub(r'OTU', 'Tax_Fun', open(self.output_dir + '/real_dc_otu_degree.txt').read()))
            open(self.output_dir + "/real_degree.txt", 'w').write(
                re.sub(r'OTU', 'Tax_Fun', open(self.output_dir + '/network_degree.txt').read()))
            os.remove(self.output_dir + "/real_dc_sample_otu_degree.txt")
            os.remove(self.output_dir + "/real_dc_otu_degree.txt")
            os.remove(self.output_dir + "/real_node_table.txt")
            os.remove(self.output_dir + "/network_degree.txt")
        except Exception as e:
            raise Exception("修改文件header失败——{}".format(e))
