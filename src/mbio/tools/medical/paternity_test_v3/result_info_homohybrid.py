# !/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/python
# -*- coding: utf-8 -*-
# __author__ = "hongdong"
# last_modify:20161123

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class ResultInfoHomohybridAgent(Agent):
    """
    亲子鉴定的结果输出
    包含家系图，存入报告的图、胎儿浓度等,这里是胎儿位点是纯合与杂合的结果图片生成脚本。生成snp分型图，家系位点匹配情况图以及测试位点
    包含脚本：plot_update.R、convert2png.sh
    如果结果中有样本有问题（如测序深度过低）的，就不生成结果图片，后续判断为异常家系
    version v1.0
    author: hongdong
    last_modify: 2016.11.21
    last modified by hongdong@20180821
    """
    def __init__(self, parent):
        super(ResultInfoHomohybridAgent, self).__init__(parent)
        options = [
            {"name": "tab_merged", "type": "infile", "format": "paternity_test.rdata"},  # format:Rdata
            {"name": "is_free_combine", "type": "bool", "default": False}
        ]
        self.add_option(options)
        self.step.add_steps("result_info")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)
        self.queue = "gada"

    def stepstart(self):
        self.step.result_info.start()
        self.step.update()

    def stepfinish(self):
        self.step.result_info.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('tab_merged'):
            raise OptionError("必须提供合并之后的家系表")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 3
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            ["family.png", "png", "家系图"],
        ])
        super(ResultInfoHomohybridAgent, self).end()


class ResultInfoHomohybridTool(Tool):
    """
    蛋白质互作组预测tool
    """
    def __init__(self, config):
        super(ResultInfoHomohybridTool, self).__init__(config)
        self._version = '1.0.1'
        self.R_path = 'program/R-3.3.1/bin/'
        self.script_path = self.config.PACKAGE_DIR + '/medical/paternity_test_v3/'
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/ImageMagick/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/bin")
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/medical/svg2png')
        self.report_point = self.config.SOFTWARE_DIR + "/database/human/pt_ref"  # 报告背景点的文件所在路径

    def run_tf(self):
        data_name = self.option("tab_merged").prop['path'].split('/')[-1]
        dad = data_name.split('_')[0]
        mom = data_name.split('_')[1]
        preg = data_name.split('_')[2]
        if dad != 'NA' and mom != 'NA'and preg != 'NA':
            is_free_combine = "true" if self.option("is_free_combine") else "false"
            plot_cmd = "{}Rscript {}plot_update.R {} {} {}".\
                format(self.R_path, self.script_path, self.option("tab_merged").prop['path'], self.report_point,
                       is_free_combine)
            self.logger.info(plot_cmd)
            self.logger.info("开始运行结果信息图的绘制")
            cmd = self.add_command("plot_cmd", plot_cmd).run()
            self.wait(cmd)
            if cmd.return_code == 0:
                self.logger.info("运行绘制结果图成功")
            else:
                self.set_error("运行绘制结果图出错")
        else:
            self.logger.info("家系中有样本质控不合格")

        path_fig1 = None
        path_family_match = None

        result = os.listdir(self.work_dir)
        for file in result:
            if re.search('.*NA.*', file):
                break
            fig1 = re.search("(.*fig1_update)\.svg", file)
            family_match = re.search("(.*family_match_update)\.svg", file)
            if fig1:
                fig1_name = fig1.group(1)
                path_fig1 = os.path.join(self.work_dir, fig1_name)
            elif family_match:
                family_match_name = family_match.group(1)
                path_family_match = os.path.join(self.work_dir, family_match_name)
        if path_fig1 and path_family_match:
            convert_cmd = "bioinfo/medical/scripts/convert2png_update.sh {} {}".format(path_fig1, path_family_match)
            self.logger.info(convert_cmd)
            self.logger.info("开始运行结果图的转化")
            cmd = self.add_command("convert_cmd", convert_cmd).run()
            self.wait(cmd)
            if cmd.return_code == 0:
                self.logger.info("运行转化结果图成功")
            else:
                self.set_error("运行转化结果图出错")

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        results = os.listdir(self.work_dir)
        for f in results:
            if re.search(r'.*.png$', f):
                os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
            if re.search(r'.*info_show\.txt$', f):
                os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
            if re.search(r'.*test_pos_update\.txt$', f):
                os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
            if re.search(r'.*num_pos_map\.txt$', f):
                os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
        self.logger.info('设置文件夹路径成功')

    def run(self):
        super(ResultInfoHomohybridTool, self).run()
        self.run_tf()
        self.set_output()
        self.end()
