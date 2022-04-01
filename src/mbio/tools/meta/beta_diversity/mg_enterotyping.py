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


class MgEnterotypingAgent(Agent):
    """
    meta之样本菌群分型分析
    version v1
    author: zhouxuan
    last_modify: 2017.11.09
    # last_modified by zouxuan
    """

    def __init__(self, parent):
        super(MgEnterotypingAgent, self).__init__(parent)
        options = [
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入的文件
            # {"name": "result_dir", "type": "outfile", "format": "meta.beta_diversity.result_dir"}  # 输出文件夹
            {"name": "dis_matrix", "type": "infile",
             "format": "meta.beta_diversity.distance_matrix"}  # 输入距离矩阵 add by zouxuan 20171120
        ]
        self.add_option(options)
        self.step.add_steps("enterotyping")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.enterotyping.start()
        self.step.update()

    def stepfinish(self):
        self.step.enterotyping.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('otu_table'):
            raise OptionError('必须输入正确的丰度表', code="32702501")
        # if self.option('otu_table').prop['sample_num'] < 18:
        #     raise OptionError('样品数必须大于18')
        # if self.option('otu_table').prop['otu_num'] < 10:
        #     raise OptionError('物种/功能数必须大于10')
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
        super(MgEnterotypingAgent, self).end()


class MgEnterotypingTool(Tool):
    def __init__(self, config):
        super(MgEnterotypingTool, self).__init__(config)
        self._version = "v1"
        self.enterotyping_path = self.config.PACKAGE_DIR + '/beta_diversity'  # 2017.11.9 by zouxuan
        self.perl_path = 'program/perl-5.24.0/bin/perl '
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'  # 2016.12.26 by zhouxuan 3
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.R_path = 'program/R-3.3.1/bin/Rscript'
        # self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')  # 2016.12.26 by zhouxuan

    def run(self):
        """
        运行
        :return:
        """
        super(MgEnterotypingTool, self).run()
        # self.logger.info(self.option('otu_table').prop['sample_num'])
        self.run_enterotyping()
        self.set_output()
        self.end()

    def run_enterotyping(self):
        """
        运行perl脚本，进行样本菌群分型分析
        """
        if self.option("dis_matrix").is_set:
            dist = self.option("dis_matrix").prop['path']
        else:
            dist = ""
        otu_table = self.option("otu_table").prop['path']
        self.name_to_name = mask_taxon(otu_table, self.work_dir + "/tmp_mask_otu.xls")  # add by guhaidong 20171025
        otu_table = self.work_dir + '/tmp_mask_otu.xls'  # add by guhaidong 20171025
        cmd = self.perl_path + self.enterotyping_path + (
            '/enterotyping.pl -i %s -o %s -d %s' % (otu_table, self.output_dir, dist))
        self.logger.info('运行perl脚本，进行样本菌群分型分析')
        command = self.add_command("enterotyping_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("enterotyping生成r文件完成!")
        else:
            self.set_error("enterotyping生成r文件出错!", code="32702501")
        cmd_1 = '%s %s' % (self.R_path, self.work_dir + "/Enterotyping.r")
        self.logger.info(cmd_1)
        command1 = self.add_command('cmd_1', cmd_1)
        command1.run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("R程序计算enterotyping成功")
        else:
            self.set_error("R程序计算enterotyping失败", code="32702502")
            raise Exception("R程序计算enterotyping失败")

    def set_output(self):  # 20170412
        """
        结果已经生成在out_put文件夹下面，判断是否正确生成数据
        :return:
        """
        self.logger.info("判断结果文件是否完整")
        file_name = os.listdir(self.output_dir)
        if len(file_name) == 2:
            self.logger.info("样本量或物种分类过低，不能进行菌群分型分析")
            raise OptionError('样本量或物种分类过低，不能进行菌群分型分析，请选择较多的样本或者较高的分类水平', code="32702502")
        else:
            self.logger.info("样本菌群分型分析成功")
        file_name.remove("ch.txt")
        file_name.remove("cluster.txt")
        for name in file_name:
            self.add_taxon(os.path.join(self.output_dir, name))

    def dashrepl(self, matchobj):
        """
        add func by guhaidong 20171025
        """
        return self.name_to_name[matchobj.groups()[0]]

    def add_taxon(self, taxon_result):
        """
        add func by zouxuan 201711117
        description: 将旧注释的名称，根据词典替换成新注释名称
        """
        old_result = taxon_result + "bak"
        os.rename(os.path.join(self.output_dir, taxon_result), os.path.join(self.output_dir, old_result))
        with open(old_result, "r") as f, open(taxon_result, "w") as w:
            # w.write(old_result)
            for i in f.readlines():
                # line = i.strip()
                new_line = re.sub(r"(name\d+)", self.dashrepl, i)
                w.write(new_line)
        os.remove(os.path.join(self.output_dir, old_result))
