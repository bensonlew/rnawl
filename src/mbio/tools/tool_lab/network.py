# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class NetworkAgent(Agent):
    """
    将成对fastq文件转换成fasta格式的文件
    version v1
    author：qiuping
    last_modify:2015.01.06
    """
    def __init__(self, parent):
        super(NetworkAgent, self).__init__(parent)
        options = [
            {"name": "link_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "node_table", "type": "infile", "format": "tool_lab.simple"},
        ]
        self.add_option(options)

    def end(self):
        super(NetworkAgent, self).end()

    def check_options(self):
        """
        检查参数设置
        :return:
        """
        if not self.option("link_table").is_set:
            raise OptionError("必须输入link_table文件")
        if not self.option("node_table").is_set:
            raise OptionError("必须输入node_table文件")

    def set_resource(self):
        """
        设置所需资源
        :return:
        """
        self._cpu = 2
        self._memory = "5G"


class NetworkTool(Tool):
    def __init__(self, config):
        super(NetworkTool, self).__init__(config)

    def run(self):
        """
        运行
        :return:
        """
        super(NetworkTool, self).run()
        self.network()
        self.end()

    def network(self):
        self.logger.info("network")
        with open(self.option("node_table").prop["path"]) as f, \
             open(self.option("link_table").prop["path"]) as f1,\
                open(self.output_dir + "/newlinks", "w") as w,\
                open(self.output_dir + "/newnodes", "w") as w1:

            lines = f.readlines()
            dict_nodes = {}
            list_node = []
            n = 0
            w1.write(lines[0])
            for line in lines[1:]:
                item = line.strip().split("\t")
                if (len(item) > 1) & (item[1] != ""):
                    w1.write(str(item[0]) + "\t" + str(item[1]) + "\n")
                else:
                    w1.write(str(item[0]) + "\t" + "NA" + "\n")
                if item[0] not in dict_nodes.keys():
                    list_node.append(item[0])
                    dict_nodes[item[0]] = n
                    n += 1
                else:
                    raise Exception("node表中有重复节点！")
            rows = f1.readlines()
            w.write(rows[0])
            for row in rows[1:]:
                ele = row.strip().split("\t")
                if (ele[0] in dict_nodes.keys()) & (ele[1] in dict_nodes.keys()):
                    if (len(ele) > 2) & (ele[2] != ""):
                        w.write(str(dict_nodes[ele[0]]) + "\t" + str(dict_nodes[ele[1]]) + "\t" + str(ele[2]) + "\n")
                    else:
                        w.write(str(dict_nodes[ele[0]]) + "\t" + str(dict_nodes[ele[1]]) + "\t" + "NA" + "\n")
                else:
                    raise Exception("link表中的节点在node表中没有写全！")
