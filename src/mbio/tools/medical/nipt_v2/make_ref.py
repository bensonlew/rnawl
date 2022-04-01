## !/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "hongdongxuan"

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import os
import re


class MakeRefAgent(Agent):
    """
    调用make_ref.R脚本，构建出参考组
    example: ~/app/program/R-3.3.1/bin/Rscript make_ref.R "WS17080003,WS17080006,WS17080010,WS17080015" 1 /mnt/ilustre/
    users/sanger-dev/app/database/human/hg38_nipt/bed_file /mnt/ilustre/users/sanger-dev/app/database/human
    /hg38_nipt/ref_cor
    version v1.0
    author: hongdongxuan
    last_modify: 20170821
    """
    def __init__(self, parent):
        super(MakeRefAgent, self).__init__(parent)
        options = [
            {"name": "sample_txt", "type": "infile", "format": "nipt.bed"},  # 训练组中样本名称
            {"name": "bed_dir", "type": "string"},  # 训练组样本的bed文件
            {"name": "ref_cor", "type": "string"},  # 训练组构建好了存放的路径
            {"name": "ref_group", "type": "int", "default": 1}  # 选择要构建的参考组的名字
        ]
        self.add_option(options)
        self.step.add_steps("bed_analysis")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.bed_analysis.start()
        self.step.update()

    def stepfinish(self):
        self.step.bed_analysis.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("sample_txt").is_set:
            raise OptionError("必须输入sample_txt文件！")
        if not self.option("bed_dir"):
            raise OptionError("必须输入bed_dir！")
        if not self.option("ref_cor"):
            raise OptionError("必须输入ref_cor值！")
        if not self.option("ref_group"):
            raise OptionError("必须输入ref_group值！")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 4
        self._memory = '20G'

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "结果输出目录"],
        # ])
        # result_dir.add_regexp_rules([
        #     ["z.xls", "xls", "z值分析结果"],
        #     ["zz.xls", "xls", "统计zz值"]
        # ])
        super(MakeRefAgent, self).end()


class MakeRefTool(Tool):
    """
    nipt bed文件分析tool
    """
    def __init__(self, config):
        super(MakeRefTool, self).__init__(config)
        self._version = '1.0.1'
        self.r_path = 'program/R-3.3.1/bin/Rscript'
        self.script_path = os.path.join(Config().SOFTWARE_DIR, "bioinfo/medical/scripts/")
        self.ref_cor_cn = Config().SOFTWARE_DIR + "/database/human/hg38" \
                                                  "_nipt/ref_cor/" + str(self.option("ref_group")) + ".ref.cor.Rdata"
        self.ref_cor_cn_new = Config().SOFTWARE_DIR + "/database/human/hg38_nipt/re" \
                                                      "f_cor/" + str(self.option("ref_group")) + ".ref.cor.Rdata_bak"
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.ref_cor_cn_dir = Config().SOFTWARE_DIR + "/database/human/hg38_nipt/ref_cor"

    def run_make_ref(self):
        sample_list = ''
        with open(self.option('sample_txt').prop['path'], 'r') as r:
            data = r.readlines()
            for m in data:
                sample_list += str(m.strip())
        self.logger.info("sample_data:{}".format(sample_list))
        self.file_check(sample_list)
        if os.path.exists(self.ref_cor_cn):
            os.rename(self.ref_cor_cn, self.ref_cor_cn_new)
            self.logger.info("参考组重命名成功！")
        self.logger.info("开始进行构建参考组")
        one_cmd = self.r_path + " %smake_ref.R %s %s %s %s" % (self.script_path, sample_list,
                                                               self.option("ref_group"), self.option("bed_dir"),
                                                               self.option("ref_cor"))
        # self.logger.info(one_cmd)
        self.logger.info("开始运行one_cmd")
        cmd = self.add_command("one_cmd", one_cmd).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("运行one_cmd成功")
        else:
            raise Exception("输入文件的格式不正确，请确保所有的样本名是用“，”隔开")

    def set_output(self):
        """
        将结果文件link到output文件夹下面,暂时这一步没有用
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        results = os.listdir(self.work_dir)
        for m in results:
            if str(m) == str(self.option("ref_group")) + ".ref.cor.Rdata":
                os.link(self.work_dir + '/' + m, self.output_dir + "/" + m)
        self.logger.info('设置文件夹路径成功')

    def run(self):
        super(MakeRefTool, self).run()
        self.run_make_ref()
        # self.set_output()
        self.end()

    def file_check(self, sample_list):
        """
        用于检查样本在bed_dir中是否存在
        :param sample_list:
        :return:
        """
        self.logger.info("开始检查样本在bed_dir中是否存在")
        sample = sample_list.strip().split(",")
        for m in sample:
            if not os.path.exists(os.path.join(self.option("bed_dir"), m + ".bed.2")):
                self.logger.error("bed_dir中样本 {} 不存在！".format(m))
                raise Exception("bed_dir中样本 {} 不存在！".format(m))
        self.logger.info("检查完成，所有的样本都在bed_dir中！")