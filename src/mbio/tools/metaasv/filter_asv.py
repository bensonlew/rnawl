# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import json
import re
import os
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import calculate_asv


class FilterAsvAgent(Agent):
    """
    根据传入的json，对一张OTU表进行过滤，过滤的条件有5种:
    species_filter: 物种过滤，用于保留或者滤去特定的物种
    sample_filter: 用于滤去在x个样本中相对丰度数小于y的ASVs
    reads_filter: 用于滤去所有序列数小于x的ASV
    func_filter: 去除线粒体或叶绿体的ASVs
    set_list：筛选asv集合；
    """
    def __init__(self, parent):
        super(FilterAsvAgent, self).__init__(parent)
        options = [
            {"name": "in_otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入的ASV文件
            {"name": "filter_json", "type": "string", "default": ""},  # 输入的json文件
            {"name": "group_table", "type": "infile", 'format': "meta.otu.group_table"},#输入Group表
            {"name": "out_otu_table", "type": "outfile", "format": "meta.otu.otu_table"}  # 输出的结果OTU表
        ]
        self.add_option(options)
        self.step.add_steps("filter_otu")

    def start_filter_otu(self):
        self.step.filter_otu.start()
        self.step.update()

    def end_filter_otu(self):
        self.step.filter_otu.end()
        self.step.update()

    def check_options(self):
        """
        参数检测
        """
        if not self.option("in_otu_table").is_set:
            raise OptionError("输入的OTU文件不能为空")
        if self.option("filter_json") == "":
            raise OptionError("输入的筛选条件不能为空")

    def end(self):
        """
        结束
        :return:
        """
        super(FilterAsvAgent, self).end()

    def set_resource(self):
        """
        设置所需的资源
        """
        self._cpu = 2
        self._memory = "10G"


class FilterAsvTool(Tool):
    def __init__(self, config):
        super(FilterAsvTool, self).__init__(config)
        self.keep_list = list()  # 处理物种筛选中的保留的逻辑
        self.python = "miniconda2/bin/python"
        self.python_script = os.path.join(self.config.PACKAGE_DIR, 'metaasv/filter_asv.py')

    def run_get_relative(self):
        """
        运行 函数 calculate_abundance 得到相对丰度数据
        :return:
        """
        self.logger.info("开始得到相对丰度文件！")
        asv_file = os.path.join(self.work_dir, "asv_table.xls")
        if os.path.exists(asv_file):
            os.remove(asv_file)
        os.link(self.option("in_otu_table").prop["path"], asv_file)
        self.logger.info(asv_file)
        try:
            self.asv_filter_input = calculate_asv(asv_file, abundance_method="relative")
        except:
            self.set_error("创建相对丰度表失败！")

    def run_filter(self):
        """
        运行脚本 filter_asv.py 过滤
        :return:
        """
        outfile = os.path.join(self.work_dir, "ASV_abundance.xls")
        if os.path.exists(outfile):
            os.remove(outfile)
        input_file = os.path.join(self.work_dir,"abundance.percents.xls")
        if not os.path.exists(input_file):
            self.set_error("没有生成正确的相对丰度结果文件！")
        self.logger.info("filter_json: {}".format(self.option("filter_json")))
        filter_json = json.dumps(self.option("filter_json"))
        cmd = "{} {} -i {} -j {} -o {} ".format(self.python, self.python_script, input_file,filter_json , outfile)
        if hasattr(self.config, "DBVersion"):
            cmd += ' -d ' + str(self.config.DBVersion)
        self.logger.info(cmd)
        command = self.add_command("filter", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行完成！")
        else:
            self.set_error("运行出错！")

    def pick_absolute(self):
        """
        根据过滤后的相对丰度表，从绝对丰度中挑选出过滤后的绝对丰度表
        :return:
        """
        relative = os.path.join(self.work_dir, "ASV_abundance.xls")
        absolute = os.path.join(self.work_dir, "ASV_absolute.xls")
        asv_list = []
        with open(relative, "r") as f:
            lines = f.readlines()
            for line in lines:
                if re.search(r"ASV ID", line):
                    pass
                else:
                    line = line.strip().split("\t")
                    asv_name = line[0]
                    if asv_name not in asv_list:
                        asv_list.append(asv_name)
        with open(self.option("in_otu_table").prop["path"], 'r') as m, open(absolute, 'w') as w:
            lins = m.readlines()
            for lin in lins:
                if re.search(r"ASV ID", lin):
                    w.write(lin)
                else:
                    lin = lin.strip().split("\t")
                    abso_asv = lin[0]
                    if abso_asv in asv_list:
                        w.write("\t".join(lin) + "\n")

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("正在设置结果文件目录")
        relative = os.path.join(self.work_dir, "ASV_abundance.xls")
        absolute = os.path.join(self.work_dir, "ASV_absolute.xls")
        output_relative = os.path.join(self.output_dir, "ASV_abundance.xls")
        if os.path.exists(output_relative):
            os.remove(output_relative)
        os.link(relative, output_relative)
        output_absolute = os.path.join(self.output_dir, "ASV_absolute.xls")
        if os.path.exists(output_absolute):
            os.remove(output_absolute)
        os.link(absolute, output_absolute)
        self.option("out_otu_table").set_path(output_absolute)
        self.logger.info("完成结果文件的设置")

    def run(self):
        """
        运行
        :return:
        """
        super(FilterAsvTool, self).run()
        self.run_get_relative()
        self.run_filter()
        self.pick_absolute()
        self.set_output()
        self.logger.info("OTU过滤完成")
        self.end()
