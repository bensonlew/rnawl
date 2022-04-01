#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__: qingchen.zhang
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError



class RankAbundanceAgent(Agent):
    """
    Metaasv rank_abundance分析
    """
    def __init__(self, parent):
        super(RankAbundanceAgent, self).__init__(parent)
        options = [
            {"name": "otu_table", "type": "infile", "format": "metaasv.otu_table,metaasv.tax_summary_dir"},  # 输入文件
            {"name": "step", "type": "int", "default": 1},  # 取样步长
            {"name": "method", "type": "string", "default": "relative"},  # 计算相对丰度
        ]
        self.add_option(options)
        self.step.add_steps('rank_abundance')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.rank_abundance.start()
        self.step.update()

    def step_end(self):
        self.step.rank_abundance.finish()
        self.step.update()

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option("otu_table").is_set:
            raise OptionError("请选择ASV表")

    def set_resource(self):
        """
            所需资源
        """
        self._cpu = 1
        self._memory = '15G'

    def end(self):
        """
        结束
        """
        super(RankAbundanceAgent, self).end()


class RankAbundanceTool(Tool):
    """
    计算排序曲线图
    """

    def __init__(self, config):
        super(RankAbundanceTool, self).__init__(config)
        self.python = "program/Python/bin/python"
        self.python_script = os.path.join(self.config.PACKAGE_DIR, 'metaasv/rank_abundance.py')
        self.step = self.option("step")

    def get_freq(self, asv_table):
        """
        功能：得到asv表的最小样本数的step
        :param otu_table:
        :return:
        """
        with open(asv_table, "r") as f:
            seq_num_list = []
            sample_num = len(f.readline().strip().split("\t")) - 1
            n = 0
            while n < sample_num:
                seq_num_list.append(0)
                n += 1
            for line in f:
                line = line.strip().split("\t")
                for k, v in enumerate(seq_num_list):
                    seq_num_list[k] += int(float(line[k + 1]))
            min_seq = min(seq_num_list)
            if min_seq < 10000:
                freq = 1
            else:
                freq = int(round(min_seq / 10000.0) * 10)
        return freq

    def run_rank(self):
        """
        运行脚本 rank_abundance.py 计算结果
        :return:
        """
        outfile = os.path.join(self.work_dir, "Rank_abundance.xls")
        if os.path.exists(outfile):
            os.remove(outfile)
        cmd = '{} {} -i {} -m {} -s {} -o {}'.format(self.python, self.python_script, self.option("otu_table").prop['path'], self.option("method"), self.step, outfile)
        self.logger.info(cmd)
        command = self.add_command("rank", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行完成！")
        else:
            self.set_error("运行出错！")

    def set_output(self):
        """
        设置结果文件目录
        """
        self.logger.info("开始设置结果文件目录")
        outfile = os.path.join(self.output_dir, "Rank_abundance.xls")
        if os.path.exists(outfile):
            os.remove(outfile)
        os.link(os.path.join(self.work_dir, "Rank_abundance.xls"), outfile)
        self.logger.info("设置结果文件目录完成")

    def run(self):
        """
        运行
        """
        super(RankAbundanceTool, self).run()
        self.step = self.get_freq(self.option("otu_table").prop['path'])
        self.run_rank()
        self.set_output()
        self.end()
