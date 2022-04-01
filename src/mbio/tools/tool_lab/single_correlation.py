# -*- coding: utf-8 -*-
# author: xueqinwen
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess
import os
import unittest
import re
import pandas as pd


class SingleCorrelationAgent(Agent):
    """
    单因子相关系分析
    Xue Qinwen
    20210316
    """
    def __init__(self, parent):
        super(SingleCorrelationAgent, self).__init__(parent)
        options = [
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "p_value", "type": "float", "default": 0.05},
            {"name": "cor_value","type":"float","default": 0.6},
            {"name": "method", "type": "string", "default": "pearsonr"},
        ]
        self.add_option(options)
        self.step.add_steps('pearsons_correlation')

    def check_options(self):
        if not self.option("otutable").is_set:
            raise OptionError('必须提供otu表', code="34100902")
        if self.option("method") not in ["pearsonr", "spearmanr", "kendalltau"]:  # add "kendalltau"(kendall) by zhujuan
            raise OptionError('不支持该相关系数方法', code="34100903")
        self.option('otutable').get_info()
        if self.option('envtable').is_set:
            pass
        else:
            raise OptionError('请选择环境因子表', code="34100905")
        if self.option('envtable').is_set:
            env_table = self.option('envtable').prop['path']
            sample_e = []
            with open(env_table, 'r') as e:
                e.next()
                for line in e:
                    line = line.strip().split('\t')
                    # if line[0] != '#SampleID':
                    sample_e.append(line[0])
                    for i in range(1, len(line)):
                        if float(line[i]) or line[i] == '0' or float(line[i]) == 0.0:  #modify by zhujuan 20171115
                            continue
                        else:
                            raise OptionError('环境因子表中存在分类型环境因子', code="34100907")

            otu_path = self.option("otutable").prop['path']
            with open(otu_path, 'r') as o:
                line = o.readline()
                line = line.strip().split('\t')
                sample_o = line[1:]
                self.logger.info(sample_e)
                self.logger.info(sample_o)
            for i in sample_o:
                if i in str(sample_e):
                    continue
                else:
                    sample_name = i
                    self.logger.info(i)
                    raise OptionError('OTU表中的样本%s和环境因子表中的样本不一致，请剔除OTU中非法样本！', variables=(sample_name), code="34100908")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 5
        self._memory = '5G'

    def end(self):
        super(SingleCorrelationAgent, self).end()


class SingleCorrelationTool(Tool):
    def __init__(self, config):
        super(SingleCorrelationTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.cmd_path = '{}/program/Python/bin/python {}/statistical/pearsonsCorrelation.py'\
            .format(self.config.SOFTWARE_DIR, self.config.PACKAGE_DIR)
        self.env_table = self.option("envtable").prop['path']
        self.real_otu = self.option("otutable").prop['path']

    def get_new_env(self):
        """
        根据挑选后的样品生成新的envtable  # add by zhujuan  20171106 解决挑完样品后环境因子数据都为零的时bug
        """
        self.logger.info("获得整理后的环境因子表")
        new_path = self.work_dir + '/temp_env_table.xls'
        sample_o = []
        otu_path = self.option("otutable").prop['path']
        otu_table = pd.DataFrame(pd.read_table(otu_path, sep='\t', index_col=0))
        env_table = self.option('envtable').prop['path']
        env_abund = pd.DataFrame(pd.read_table(env_table, sep='\t', index_col=0))
        name = env_abund.index.name
        env_abund.index = [str(i) for i in env_abund.index]
        new_env_abund = env_abund.ix[otu_table.columns]
        new_env_abund.index.name = name
        a = (new_env_abund != 0).any(axis=0)
        env_empt = []
        for i in range(len(new_env_abund.columns)):
            if a[i]:
                pass
            else:
                sample_name = new_env_abund.columns[i]
                env_empt.append(sample_name)
        if env_empt:
            self.set_error('环境因子：%s在所选样品中均为0，请剔除该环境因子!', variables=(env_empt), code="34100901")
        new_env_abund.to_csv(new_path, sep="\t", encoding="utf-8")
        return new_path
  

    def run(self):
        """
        运行
        """
        super(SingleCorrelationTool, self).run()
        self.env_table = self.get_new_env()
        self.run_pearsonsCorrelation()
        self.set_output()
        self.end()

    def run_pearsonsCorrelation(self):
        """
        run pearsonsCorrelation.py
        """
        self.logger.info("运行脚本开始对otu表与环境因子进行计算相关性")
        self.check_env_table(self.env_table)
        self.logger.info("对环境因子数据检查完成!")
        cmd = self.cmd_path
        cmd += " %s %s %s %s %s" % (self.real_otu, self.env_table, "./{}_correlation_at_otu_level.xls".format(self.option("method"))
                                    , "./{}_pvalue_at_otu_level.xls".format(self.option("method")),
                                    self.option("method"))
        self.logger.info('运行pearsonsCorrelation.py计算correlation')
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('Pearsons Correlation 计算成功')
        except subprocess.CalledProcessError:
            self.logger.info('Pearsons Correlation 计算失败')
            self.set_error('pearsonsCorrelation.py 计算失败')
        self.logger.info('运行pearsonsCorrelation.py计算correlation完成')

    def check_env_table(self,table):
        """
        对环境因子的每一列进行检查，如果相同则报错；如果不相同则正常进行分析
        add by qingchen.zhang @20200305
        """
        self.logger.info("对环境因子数据进行检查")
        data = pd.read_table(table, sep="\t", header=0)
        env_name_list = list(data.columns)
        for env_name in env_name_list:
            env_factor_list = set(list(data[env_name]))
            if len(env_factor_list) == 1:
                self.set_error("所选样本的环境因子值(%s)相等，无法计算，建议重新选择环境因子!", variables=(env_name), code="34100905")
                break
            else:
                continue

    def set_output(self):
        self.logger.info("开始链接和检查结果文件")
        newpath = self.output_dir + "/{}_correlation.xls".format(self.option("method"))
        if os.path.exists(newpath):
            os.remove(newpath)
        os.link(self.work_dir + "/{}_correlation_at_otu_level.xls".format(self.option("method")),
                self.output_dir + "/{}_correlation.xls".format(self.option("method")))
        newpath2 = self.output_dir + "/{}_pvalue.xls".format(self.option("method"))
        if os.path.exists(newpath2):
            os.remove(newpath2)
        os.link(self.work_dir + "/{}_pvalue_at_otu_level.xls".format(self.option("method")),
                self.output_dir + "/{}_pvalue.xls".format(self.option("method")))
        result_dict = {}
        with open(self.work_dir + "/{}_pvalue_at_otu_level.xls".format(self.option("method")),"r") as old_p:
            old_p.readline()
            while 1:
                line1 = old_p.readline()
                # line2 = old_co.readline
                if not line1:
                    break
                # if not line2:
                #     break
                fd1 = line1.rstrip().split('\t')
                # fd2 = line2.rstrip().split('\t')
                result_dict[fd1[0]] = {}
                result_dict[fd1[0]]["p"] = fd1[1]
        with open(self.work_dir + "/{}_correlation_at_otu_level.xls".format(self.option("method")),"r") as old_co:
            old_co.readline()
            while 1:
                # line1 = old_p.readline()
                line2 = old_co.readline()
                if not line2:
                    break
                # fd1 = line1.rstrip().split('\t')
                fd2 = line2.rstrip().split('\t')
                result_dict[fd2[0]]["co"] = fd2[1]
        with open(self.output_dir + "/plot_pvalue.xls", "w") as new_p,open(self.output_dir + "/plot_correlation.xls", "w") as new_co:
            for sn in result_dict.keys():
                if float(result_dict[sn]["p"])< self.option("p_value") and abs(float(result_dict[sn]["co"]))>= self.option("cor_value"):
                    new_co.write("{}\t{}\n".format(sn,result_dict[sn]["co"]))
                    new_p.write("{}\t{}\n".format(sn,result_dict[sn]["p"]))

class TestFunctionf(unittest.TestCase):
    '''
    测试脚本
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "single_co_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.single_correlation",
            "options": dict(

                # bam_dir="/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/de_tools/bam_dir",
                otutable="/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/single_corrrelation/otu_taxon.txt",
                envtable="/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/single_corrrelation/env_table.txt",
                # c_bam_dir="/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/de_tools/c_bam",
                
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
