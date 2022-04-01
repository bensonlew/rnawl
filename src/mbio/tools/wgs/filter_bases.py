# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180410

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import re
import os


class FilterBasesAgent(Agent):
    """
     用于对vcf结果文件进行简并碱基的过滤，简并碱基的存在会影响到gatk筛选突变位点的计算，原理是所有的REF与ALT中碱基必须都是ATCG
    """
    def __init__(self, parent):
        super(FilterBasesAgent, self).__init__(parent)
        options = [
            {"name": "pop_var_vcf", "type": "string"}
        ]
        self.add_option(options)
        self.step.add_steps('gatk')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.gatk.start()
        self.step.update()

    def step_end(self):
        self.step.gatk.finish()
        self.step.update()
        
    def check_options(self):
        if not self.option("pop_var_vcf"):
            raise OptionError("缺少pop_var_vcf参数", code="34502601")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 5
        self._memory = '15G'
        
    def end(self):
        super(FilterBasesAgent, self).end()


class FilterBasesTool(Tool):
    def __init__(self, config):
        super(FilterBasesTool, self).__init__(config)
        
    def filter_vcf_run(self):
        file_path = self.output_dir + "/pop.variant.vcf"
        # if os.path.exists(file_path):
        #     os.remove(file_path)
        bases = ["A", "T", "C", "G"]
        with open(self.option("pop_var_vcf"), 'r') as r, open(file_path, 'w') as w:
            for line in r:
                if re.match(r'^#', line) or re.match(r'^##', line):
                    w.write(line)
                else:
                    temp_bases = []
                    tip = True
                    temp = line.strip().split("\t")
                    for m in temp[3].strip().split(','):
                        temp_bases.append(m)
                    for n in temp[4].strip().split(','):
                        temp_bases.append(n)
                    for l in temp_bases:
                        for m in l:
                            if m.upper() not in bases:
                                tip = False
                                self.logger.info("{}不在ATGC,该行要过滤掉".format(m.upper()))
                                break
                    if tip:
                        w.write(line)

    def run(self):
        super(FilterBasesTool, self).run()
        self.logger.info("开始进行vcf过滤！")
        self.filter_vcf_run()
        self.logger.info("完成vcf过滤！")
        self.end()
