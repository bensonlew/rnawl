# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
import os, re, subprocess, shutil, time
import pandas as pd
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool


class TrainAgent(Agent):
    """
    用户自己想构建疾病数学分类模型的工具，用户的数据和我们的数据整合在一起，然后构建数学模型
    author: guhaidong
    last_modified:
    """

    def __init__(self, parent):
        super(TrainAgent, self).__init__(parent)
        options = [
            {"name": "train_file", "type": "infile", "format": "sequence.profile_table"},
            {"name": "disease", "type": "string"},  # 疾病类型
            {"name": "model_type", "type": "string"}  # 建什么模型
        ]
        self.add_option(options)
        self.step.add_steps("train")
        # self.on('start', self.stepstart)
        # self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.train.start()
        self.step.update()

    def stepfinish(self):
        self.step.train.finish()
        self.step.update()

    def check_option(self):
        if not self.option('train_file').is_set:
            raise OptionError('必须输入模型构建文件')
        if not self.option('disease').is_set:
            raise OptionError('必须输入疾病类型')
        if self.option("disease") not in ["crc", "obesity", "t2d"]:
            raise OptionError("疾病不在范围内")
        if not self.option('model_type').is_set:
            raise OptionError('必须输入模型方法')
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "模型训练结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            [".+/.+model", "", "模型结果"]
        ])
        super(TrainAgent, self).end()


class TrainTool(Tool):
    def __init__(self, config):
        super(TrainTool, self).__init__(config)
        self.out_file_name = self.option('disease') + '_' + self.option("model_type") + '_model'  # 最终的模型结果
        self.our_prof = self.config.SOFTWARE_DIR + '/bioinfo/model/data/' + self.option("disease") + "_reads_percent.xls"  # 我们的丰度表
        self.our_group = self.config.SOFTWARE_DIR + '/bioinfo/model/data/' + self.option("disease") + "_group" # 我们的分组表
        self.script = self.config.PACKAGE_DIR + '/hmdb/scripts/user_train.py'
        self.python_path = '/program/Python/bin/python'
        self.train_prof = 'percent_file.txt'
        self.train_group = 'group.txt'

    def pre_model(self):
        # 丰度、分类分离、分类数必须为2
        # 丰度、分类合并，物种名称按新上传的数据为准
        # 丰度求百分数
        self.logger.info(self.our_prof)
        self.logger.info(self.our_group)
        self.our_data = pd.read_table(self.our_prof, index_col=0)
        self.our_group_data = pd.read_table(self.our_group, index_col=0, names=["Healthy"])
        self.usr_data = pd.read_table(self.option("train_file").path, index_col=0)
        group_label = self.our_group_data["Healthy"].unique().tolist()
        if "Healthy" in self.usr_data.index:
            usr_group = self.usr_data.loc["Healthy",].T
            usr_group = pd.DataFrame(usr_group)
            concat_group = pd.concat([self.our_group_data, usr_group])
            self.logger.info(usr_group)
            usr_group_label = usr_group["Healthy"].unique().tolist()
            if len(concat_group["Healthy"].unique()) != 2:
                self.set_error("错误的分类标签值：%s, 应写为：%s" % (usr_group_label.remove("Healthy")[0], group_label.remove("Healthy")[0]))
            tmp_usr_data = self.usr_data.drop(["Healthy"])
            self.usr_data = tmp_usr_data.apply(pd.to_numeric, errors='ignore')
        else:
            self.set_error("模型数据应含有名称为Healthy的分类信息行")
        new_data = pd.merge(self.usr_data, self.our_data, left_index=True, right_index=True, how="left").fillna(0)
        percent_data = new_data.apply(lambda x: x/x.sum())
        percent_data.to_csv(self.train_prof, sep="\t", index=True)
        concat_group.to_csv(self.train_group, sep="\t", header=0, index=True)

    def run_model(self):
        self.logger.info("运行模型预测")
        cmd = "%s %s -m %s -i %s  -model %s -disease %s -g %s" % (
        self.python_path, self.script, self.option('model_type'), self.train_prof, self.out_file_name,
        self.option('disease'), self.train_group)
        command = self.add_command("model", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("模型预测完成")
        else:
            self.set_error("模型预测出错")

    def set_output(self):
        self.logger.info("设置结果目录")
        outfile = os.path.join(self.output_dir, self.out_file_name)
        os.link(self.out_file_name, outfile)
        self.logger.info("设置结果目录成功")

    def run(self):
        super(TrainTool, self).run()
        self.pre_model()
        self.run_model()
        self.set_output()
        self.end()
