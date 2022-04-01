# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
import subprocess
import re
from mbio.packages.taxon.mask_taxon import mask_taxon  # add by guhaidong 20171025


class MgPlotEnterotypingAgent(Agent):
    """
    meta之样本菌群分型分析画图数据生成
    version v1
    author: zhouxuan
    last_modify: 2017.11.10 by zouxuan
    """

    def __init__(self, parent):
        super(MgPlotEnterotypingAgent, self).__init__(parent)
        options = [
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入的文件
            {"name": "g", "type": "string", "default": ""},  # 分型组别
            {"name": "s", "type": "string", "default": ""},  # 分组名称（暂定取第一个名称，散点样本不管）
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},  # 输入样本的分组信息
            {"name": "dis_matrix", "type": "infile",
             "format": "meta.beta_diversity.distance_matrix"}   #输入距离矩阵 add by zouxuan 20171120
        ]
        self.add_option(options)
        self.step.add_steps("plotenterotyping")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.plotenterotyping.start()
        self.step.update()

    def stepfinish(self):
        self.step.plotenterotyping.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('otu_table'):
            raise OptionError('必须输入正确的丰度表', code="32702601")
        if not self.option('g'):
            raise OptionError('必须输入正确的分型分组信息', code="32702602")
        if not self.option('s'):
            raise OptionError('输入信息量过小，不能进行正确的分型，请选择较多的样本或者较高的分类级别', code="32702603")
        if not self.option('group'):
            raise OptionError('必须输入正确的样本分组信息', code="32702604")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(MgPlotEnterotypingAgent, self).end()


class MgPlotEnterotypingTool(Tool):
    def __init__(self, config):
        super(MgPlotEnterotypingTool, self).__init__(config)
        self._version = "v1"
        self.plotenterotyping_path = self.config.PACKAGE_DIR + '/beta_diversity'  # 2017.11.9 by zouxuan
        self.perl_path = 'program/perl-5.24.0/bin/perl '
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'  # 2016.12.26 by zhouxuan 3
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.R_path = 'program/R-3.3.1/bin/Rscript'

    def run(self):
        """
        运行
        :return:
        """
        super(MgPlotEnterotypingTool, self).run()
        self.run_plotenterotyping()
        # self.set_output()
        self.end()

    def run_plotenterotyping(self):
        """
        运行perl脚本，进行样本菌群分型分析画图数据的生成
        """
        if self.option("dis_matrix").is_set:
            dist = self.option("dis_matrix").prop['path']
        else:
            dist = ""
        otu_table = self.option("otu_table").prop['path']
        self.name_to_name = mask_taxon(otu_table, self.work_dir + "/tmp_mask_otu.xls")  # add by guhaidong 20171025
        otu_table = self.work_dir + '/tmp_mask_otu.xls'  # add by guhaidong 20171025
        all_num = self.option('g')
        print(all_num)
        a = all_num.strip().split(",")
        num = a[-1]
        all_species = self.option('s')
        b = all_species.strip().split(",")
        species_name = []
        for i in a:
            if i == num:
                ppp = "\\" + "\"" + b[int(i) - 1] + "\\" + "\""
                species_name.append(ppp)
            # species_name = species_name + "\\" + "\"" + b[int(i)-1] + "\\" + "\""
            else:
                ppp = "\\" + "\"" + b[int(i) - 1] + "\\" + "\"" + ","
                species_name.append(ppp)
            # species_name = species_name + "\\" + "\"" + b[int(i)-1] + "\\" + "\"" + ","
            # print(species)
        species_name = ''.join(species_name)

        cmd = self.perl_path + self.plotenterotyping_path + (
            '/plot_enterotyping.pl -i %s -g %s -s %s -o %s -group %s -l T -d %s' % (
           otu_table, self.option('g'), species_name, self.output_dir,
            self.option('group').prop['path'],dist))
        self.logger.info('运行perl脚本，进行样本菌群分型分析画图数据的生成')
        command = self.add_command("plotenterotyping_cmd", cmd).run()
        self.wait(command)
        # if os.path.exists(self.work_dir +"/OUT/final_fusion_genes/" + self.option('sample_name') + "/" + self.option('sample_name') + ".final.Fusion.specific.for.genes"):
        #     if os.path.exists(self.work_dir + "/" + self.option('sample_name') + "/final_fusion_genes/" + self.option('sample_name') + "/" + self.option('sample_name') + ".final.Fusion.specific.for.trans"):
        #         if os.path.exists(self.work_dir + "/" + self.option('sample_name') + "/final_fusion_genes/" + self.option('sample_name') + "/analysis/figures/fusions_expression"):
        #             self.logger.info("soapfuse运行成功!")
        # else:
        #     self.set_error("soapfuse运行出错!")
        if command.return_code == 0:
            self.logger.info("r脚本生成完成!")
        else:
            self.set_error("r脚本生成出错!", code="32702601")
        cmd_1 = '%s %s' % (self.R_path, self.work_dir + "/plot-Enterotyping.r")
        self.logger.info(cmd_1)
        command1 = self.add_command('cmd_1', cmd_1)
        command1.run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("菌群分型分析画图数据生成成功")
        else:
            self.set_error("菌群分型分析画图数据生成失败", code="32702602")
            raise Exception("菌群分型分析画图数据生成失败")
        # def set_output(self):
        # 	"""
        # 	将结果文件复制到output文件夹下面
        # 	:return:
        # 	"""
        # 	self.logger.info("设置结果目录")
        # 	# try:
        # 	# 	shutil.copytree(self.work_dir, self.output_dir)
        # 		# shutil.copy2(self.work_dir + "/" + self.option('sample_name') + "/final_fusion_genes/" + self.option('sample_name') + "/" + self.option('sample_name') + ".final.Fusion.specific.for.trans",self.output_dir + "/final.Fusion.specific.for.trans")
        # 		# shutil.copytree(self.work_dir + "/" + self.option('sample_name') + "/final_fusion_genes/" + self.option('sample_name') + "/analysis/figures/fusions_expression" , self.output_dir + "/fusions_expression")
        # 	self.option('result_dir').set_path(self.work_dir + "/result")
        # 	self.logger.info("设置样本菌群分型分析数据目录成功")

        # except Exception as e:
        # 	self.logger.info("设置样本菌群分型分析数据结果目录失败{}".format(e))
        # 	self.set_error("设置样本菌群分型分析数据结果目录失败{}".format(e))
