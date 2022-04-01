# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os
import subprocess
import re
import shutil
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.taxon.mask_taxon import mask_taxon
from mbio.packages.metaasv.common_function import link_dir


class PermanovaAgent(Agent):
    """
    metaasv Permanova分析
    与原来不一样的地方在于增加了多因素分析结果
    """
    def __init__(self, parent):
        super(PermanovaAgent, self).__init__(parent)
        options = [
            {"name": "abu_table", "type": "infile",
             "format": "meta.otu.otu_table, meta.otu.tax_summary_dir, toolapps.table"},
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "permutations", "type": "int", "default": 999},
            {"name": "dis_method", "type": "string", "default": "bray"},
            {"name": "binary", "type": "string", "default": "false"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('abu_table').is_set:
            raise OptionError('必须提供丰度表')
        if not self.option('envtable').is_set:
            raise OptionError('必须提供环境因子表')

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(PermanovaAgent, self).end()


class PermanovaTool(Tool):
    def __init__(self, config):
        super(PermanovaTool, self).__init__(config)
        self._version = '1.0'
        self.perl_path = "program/perl-5.24.0/bin/perl"
        self.script1_path = self.config.PACKAGE_DIR + "/metaasv/permanova_multiple.pl"
        self.permanova = self.config.PACKAGE_DIR + "/metaasv/permanova.pl"
        self.R_path = 'program/R-3.3.1/bin/Rscript'

    def run_permanova(self):
        """
        运行脚本 permanova.pl
        :return:
        """
        dir_name = os.path.join(self.work_dir, "permanova")
        if os.path.exists(dir_name):
            shutil.rmtree(dir_name)
        os.mkdir(dir_name)
        otu_table = self.option("abu_table").prop['path']
        self.name_to_name = mask_taxon(otu_table, self.work_dir + "/tmp_mask_otu.xls")
        otu_table = self.work_dir + '/tmp_mask_otu.xls'
        cmd = self.perl_path
        cmd += ' %s -i %s -env %s -o %s  -pe %s -d_method %s -binary %s -multi %s' % (
            self.script1_path, otu_table, os.path.join(self.work_dir, "env_table.xls"),
            self.work_dir + '/permanova/One-factor_PERMANOVA.xls', self.option("permutations"), self.option("dis_method"),
            self.option("binary"), self.work_dir + '/permanova/Multi-factor_PERMANOVA.xls')
        self.logger.info(cmd)
        command = self.add_command('cmd', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("permanova的r文件生成成功")
        else:
            self.set_error("permanova的r文件生成失败")
            raise Exception("permanova的r文件生成失败")
        cmd_1 = '%s %s' % (self.R_path, self.work_dir + "/cmd.r")
        self.logger.info(cmd_1)
        command1 = self.add_command('cmd_1', cmd_1)
        command1.run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("R程序计算permanova成功")
        else:
            self.set_error("R程序计算permanova失败")
            raise Exception("R程序计算permanova失败")

    def dashrepl(self, matchobj):
        """
        替换名称
        """
        return self.name_to_name[matchobj.groups()[0]]

    def add_taxon(self, old_result, taxon_result):
        """
        description: 将旧注释的名称，根据词典替换成新注释名称
        """
        with open(old_result, "r") as f, open(taxon_result, "w") as w:
            # w.write(old_result)
            for i in f.readlines():
                # line = i.strip()
                new_line = re.sub(r"(name\d+)", self.dashrepl, i)
                w.write(new_line)

    def replace_env_name(self):
        """
        替换环境因子名称，避免因为环境因子顺序导致计算结果的差异
        :return:
        """
        self.env_name = {}
        env_table = os.path.join(self.work_dir, "env_table.xls")
        data = pd.read_table(self.option("envtable").prop['path'], sep="\t", header=0)
        data_list = list(data.columns)
        new_data_list = data_list[1:]
        new_data_list.sort()
        new_line = []
        for env_name in new_data_list:
            env_index = new_data_list.index(env_name)
            new_name = "name" + str(env_index + 1)
            if env_name not in self.env_name:
                self.env_name[new_name] = env_name
            if new_name not in new_line:
                new_line.append(new_name)
        all_list = [data_list[0]] + new_data_list
        new_all_list = [data_list[0]] + new_line
        data = data[all_list]
        data.columns = new_all_list
        self.logger.info("correlation_relationship: {}".format(self.env_name))
        data.to_csv(env_table, index=0, sep="\t")


        # with open(self.option("envtable").prop['path'], "r") as f, open(env_table, "w") as w:
        #     line = f.readline()
        #     line = line.strip().split("\t")
        #
        #     for env_name in line[1:]:
        #         env_index = line.index(env_name)
        #         new_name = "name" + str(env_index)
        #         if env_name not in self.env_name:
        #             self.env_name[new_name] = env_name
        #         if new_name not in new_line:
        #             new_line.append(new_name)
        #     w.write(line[0] + "\t" + "\t".join(new_line) + "\n")
        #     for lin in f:
        #         w.write(lin)

    def replace_name(self, infile1, infile2,type):
        """
        对最后的多因素计算结果添加名称
        原因是：最后的多因素结果可能会有部分的因子计算不出结果，导致结果错位
        解决方案是：根据统计结果，对多因素结果进行添加名称
        :param infile1: 统计结果表
        :param infile2: 多因素计算结果表
        :return:
        """
        if type in ["multiple"]:
            tmp_path = os.path.join(self.work_dir, "tmp.xls")
            with open(infile1, "r") as f, open(infile2, "r") as r, open(tmp_path, 'w') as w:
                lines1 = f.readlines()
                lines2 = r.readlines()
                w.write(lines2[0])
                for line in lines1[1:]:
                    line1 = line.strip().split("\t")
                    for lin in lines2[1:]:
                        line2 = lin.strip().split("\t")
                        if line1[2] == line2[1]:
                            w.write(lin)
                            break
            if os.path.exists(infile2):
                os.remove(infile2)
            if os.path.exists(tmp_path):
                os.rename(tmp_path, infile2)
        out_path = os.path.join(self.work_dir, "last_result.xls")
        factor_name_list = []
        with open(infile1, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                try:
                    factor = eval(line[0])
                except:
                    factor = line[0]
                if factor in self.env_name:
                    new_factor = self.env_name[factor]
                    if new_factor not in factor_name_list:
                        factor_name_list.append(new_factor)
        if type in ["multiple"]:
            factor_name_list.append("Residuals")
            factor_name_list.append("Total")
        with open(infile2, "r") as r, open(out_path, 'w') as w:
            lins = r.readlines()
            header = lins[0].strip().split('\t')
            if type in ["multiple"]:
                new_header = ["Name"] + header
            else:
                new_header = header
            header_string = "\t".join(new_header)
            w.write(header_string+"\n")
            n = 0
            for lin in lins[1:]:
                lin = lin.strip().split("\t")
                if type in ["multiple"]:
                    new_line = [factor_name_list[n]] + lin
                else:
                    new_line = [factor_name_list[n]] + lin[1:]
                line_string = "\t".join(new_line)
                w.write(line_string+"\n")
                n += 1
        if os.path.exists(infile2):
            os.remove(infile2)
        if os.path.exists(out_path):
            os.rename(out_path, infile2)

    def run(self):
        """
        运行
        :return:
        """
        super(PermanovaTool, self).run()
        self.replace_env_name()
        self.run_permanova()
        self.replace_name(os.path.join(self.work_dir, "test_result.xls"), os.path.join(self.work_dir, "permanova","Multi-factor_PERMANOVA.xls"), "multiple")
        self.replace_name(os.path.join(self.work_dir, "permanova","One-factor_PERMANOVA.xls"), os.path.join(self.work_dir, "permanova","One-factor_PERMANOVA.xls"), "single")

        self.set_output()

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        link_dir(os.path.join(self.work_dir,"permanova"), self.output_dir)
        self.end()
