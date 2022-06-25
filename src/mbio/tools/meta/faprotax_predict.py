# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
# @20200219

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil
import pandas as pd


class FaprotaxPredictAgent(Agent):
    """
    多样性进行FAPROTAX预测
    """
    def __init__(self, parent):
        super(FaprotaxPredictAgent, self).__init__(parent)
        options = [
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table"}, # otu丰度表，这里的otu表最后一列为taxonomy
        ]
        self.add_option(options)
        self.step.add_steps("faprotax_predict")
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.faprotax_predict.start()
        self.step.update()

    def step_end(self):
        self.step.faprotax_predict.finish()
        self.step.update()

    def check_options(self):
        """
        进行参数二次检查
        :return:
        """
        if not self.option("otu_table").is_set:
            raise OptionError("请传入otu丰度文件！")

    def set_resource(self):
        """
        设置资源
        :return:
        """
        self._cpu = 4
        self._memory = '20G'

    def end(self):
        """
        结束啦
        :return:
        """
        super(FaprotaxPredictAgent, self).end()


class FaprotaxPredictTool(Tool):
    def __init__(self, config):
        super(FaprotaxPredictTool, self).__init__(config)
        """
        设置软件、脚本、数据库的路径
        """
        self.set_environ(PYTHONPATH = self.config.SOFTWARE_DIR + "/bioinfo/meta/picrust-1.1.0/")
        self.python_scripts = "/miniconda2/bin/"
        self.faprotax_path = self.config.SOFTWARE_DIR + "/bioinfo/meta/FAPROTAX_1.2.1/" #"/mnt/ilustre/users/sanger-dev/sg-users/zhangqingchen/metagenomic/faprotax/FAPROTAX_1.2.1"

    def faprotax_predict(self):
        """
        开始用faprotax软件进行预测
        :return:
        """
        self.logger.info("开始用faprotax软件进行预测")
        output_dir = os.path.join(self.work_dir, "result")
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        else:
            os.mkdir(output_dir)
        report_path = os.path.join(output_dir, "report.txt")
        function_prediction = os.path.join(output_dir, "function_prediction.txt")
        cmd = "{}/python {}/collapse_table.py -i {} -o {} -g {}/FAPROTAX.txt -v -d 'taxonomy' -c '#' -r {} --column_names_are_in last_comment_line".format(self.python_scripts, self.faprotax_path, self.option("otu_table").prop['path'], function_prediction,self.faprotax_path, report_path)
        self.logger.info(cmd)
        command = self.add_command("faprotax_predict", cmd).run()
        self.wait(command)
        if command.return_code in [0]:
            self.logger.info("FAPROTAX运行完成！")
        else:
            self.set_error("FAPROTAX预测失败！")
        input_dir = os.path.join(self.work_dir, "result")
        all = []
        with open(os.path.join(input_dir, "function_prediction.txt")) as f:
            data = f.readlines()
            for i in data[1:]:
                for x in i.strip().split("\t")[2:]:
                    data.append(str(x))
            if all == ["0"]:
                self.set_error("FAPROTAX预测无结果！")

    def run(self):
        """
        运行
        """
        super(FaprotaxPredictTool, self).run()
        self.logger.info("开始运行预测啦")
        self.faprotax_predict()
        self.set_output()
        self.logger.info("运行tool结束")
        self.end()

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始设置结果文件目录")
        self.logger.info("开始预测bugbase的文件连接")
        input_dir = os.path.join(self.work_dir, "result")
        otu_path = os.path.join(self.output_dir, "function_prediction.txt")
        #if os.path.exists(otu_path):
        #    os.remove(otu_path)
        #os.link(os.path.join(input_dir, "function_prediction.txt"), otu_path)
        with open(input_dir+"/function_prediction.txt") as f,open(otu_path,"w") as t:
            data = f.readlines()
            t.write(data[0].strip().split("\t")[0] + "\t" + "\t".join(data[0].strip().split("\t")[2:]) + "\n")
            for i in data[1:]:
                all = 0
                for x in i.strip().split("\t")[2]:
                    all += int(x)
                if all != 0:
                    t.write(i.strip().split("\t")[0] + "\t" + "\t".join(i.strip().split("\t")[2:]) + "\n")
        report_path = os.path.join(self.output_dir, "report.txt")
        if os.path.exists(report_path):
            os.remove(report_path)
        os.link(os.path.join(input_dir, "report.txt"), report_path)
        self.logger.info("完成预测bugbase的kegg文件连接")