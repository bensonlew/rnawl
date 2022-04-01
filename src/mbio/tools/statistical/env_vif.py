# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last_modify: 2017.10.17

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
from mbio.packages.statistical.env_vif import env_vif
from mbio.packages.taxon.mask_taxon import mask_taxon,mask_env  # add by zhujuan 20180122
from itertools import islice
import re
import glob
import os

class EnvVifAgent(Agent):
    """
    计算vif方差膨胀因子的工具
    """
    def __init__(self, parent):
        super(EnvVifAgent, self).__init__(parent)
        options = [
            {"name": "abundtable", "type": "infile", "format": "meta.otu.otu_table,meta.otu.tax_summary_dir"},
            # 物种/功能/基因丰度表格
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},  # 环境因子表
            {"name": "viflim", "type": "int", "default": 10},  # 膨胀因子的筛选阈值[10-20],
            {"name": "method", "type": "string", "default": ""},  # rda|cca，默认根据 DCA result（DCA1>=3.5,CCA;DCA1<3.5,RDA）
        ]
        self.add_option(options)

    def gettable(self):
        """
        根据level返回进行计算的丰度表
        :return:
        """
        if self.option('abundtable').format == "meta.otu.tax_summary_dir":
            return self.option('abundtable').get_table(self.option('level'))
        else:
            return self.option('abundtable').prop['path']

    def check_options(self):
        """
        检查参数
        """
        if not self.option("abundtable").is_set:
            raise OptionError("请传入丰度文件！", code="34100401")
        else:
            self.option("abundtable").get_info()
            if self.option('abundtable').prop['sample_num'] < 3:
                raise OptionError('丰度表的样本数目少于3，不可进行VIF方差膨胀因子分析', code="34100402")
        if not self.option("envtable").is_set:
            raise OptionError("请传入环境因子文件！", code="34100403")
        else:
            self.option('envtable').get_info()
            with open(self.option('envtable').prop['path']) as f:
                for line in islice(f, 1, None):
                   lines = line.strip('\n').split('\t')
                   for i in range(1, len(lines)):
                      if re.match("(-{0,1})[0-9.]+$", lines[i]):
                          pass
                      elif re.match("^[a-zA-Z]+[0-9]*$", lines[i]):
                          pass
                      else:
                          raise OptionError('环境因子数据中存在非法字符,如含有空格、下划线等,请检查环境因子文件!', code="34100404")
        samplelist = open(self.gettable()).readline().strip().split('\t')[1:]
        envlist = open(self.option('envtable').prop['path']).readline().strip().split('\t')[1:]
        #for sp in self.option('envtable').prop['sample']:
        #    if sp not in samplelist:
        #         raise OptionError('提供的环境因子表中有不在样品表中存在的样品：%s' % sp)
        common_samples = set(samplelist) & set(self.option('envtable').prop['sample'])
        if len(common_samples) < 3:
            raise OptionError("环境因子表和丰度表的共有样本数必须大于等于3个：%s", variables=(len(common_samples)), code="34100405")
        if len(common_samples) <= len(envlist):
           raise OptionError("环境因子个数必须小于样品个数，请剔除多余的环境因子!", code="34100406")

        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        super(EnvVifAgent, self).end()


class EnvVifTool(Tool):
    def __init__(self, config):
        super(EnvVifTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR+"/gcc/5.1.0/lib64:$LD_LIBRARY_PATH")
        self.r_path = "/program/R-3.3.1/bin/"

    def create_otu_and_env_common(self, T1, T2, new_T1, new_T2):
        T1 = pd.read_table(T1, sep='\t', dtype=str)
        T2 = pd.read_table(T2, sep='\t', dtype=str)
        T1_names = list(T1.columns[1:])
        T2_names = list(T2.iloc[0:, 0])
        T1_T2 = set(T1_names) - set(T2_names)
        T2_T1 = set(T2_names) - set(T1_names)
        T1T2 = set(T2_names) & set(T1_names)
        if len(T1T2) < 3:
            return False
        [T1_names.remove(value) for value in T1_T2]
        T1.to_csv(new_T1, sep="\t", columns=[T1.columns[0]] + T1_names, index=False)
        indexs = [T2_names.index(one) for one in T2_T1]
        T2 = T2.drop(indexs)
        T2.to_csv(new_T2, sep="\t", index=False)
        return True

    def env_vif_r(self):
        old_abundtable = self.option("abundtable").prop['path']
        mask_taxon(old_abundtable, self.work_dir + "/tmp_mask_table.xls")
        new_abundtable = self.work_dir + "/tmp_mask_table.xls"
        env_abund = pd.DataFrame(pd.read_table(self.option("envtable").prop['path'], sep='\t', index_col=0))
        env_abund.index.name = "SampleID"
        env_path = self.work_dir + "/new_nev_table.xls"
        env_abund.to_csv(env_path, sep="\t", encoding="utf-8")
        self.abund_table = self.work_dir + '/new_table.xls'
        self.env_table = self.work_dir + '/new_env.xls'
        if not self.create_otu_and_env_common(new_abundtable, env_path, self.abund_table, self.env_table):
            self.set_error('环境因子表中的样本与丰度表中的样本共有数量少于2个', code="34100401")
        self.env_to_name = mask_env(self.env_table, self.work_dir + '/tmp_mask_env.xls')
        env_vif(self.abund_table, self.work_dir + '/tmp_mask_env.xls', self.option("viflim"), self.option('method'), self.work_dir)
        cmd = self.r_path + "Rscript run_env_vif.r"
        self.logger.info("开始运行VIF方差膨胀因子分析")
        command = self.add_command("env_vif", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("vif方差膨胀因子运行完成")
            vif_file = glob.glob(self.work_dir + "/*_vif.txt")
            for i in vif_file:
                file_name = i.split("/")[-1]
                self.add_taxon(i, self.output_dir + '/' + file_name, "env")
            os.link(self.work_dir + "/DCA.txt",self.output_dir + "/DCA.txt")
        else:
            self.set_error("vif方差膨胀因子运行出错!", code="34100402")

    def dashrepl_env(self, matchobj):
        """
        add by zhujuan 20180122
        将环境因子名称替换成name[0-9]*进行运算，后将结果中的环境因子名改回真实名，避免有些+等特殊字符不能运算
        """
        return self.env_to_name[matchobj.groups()[0]]

    def add_taxon(self, old_result, taxon_result, type):
        """
        add func by guhaidong 20171025, last modify by 20180110
        description: 将旧注释的名称，根据词典替换成新注释名称
        """
        dashrepl = ""
        if type == "tax":
            dashrepl = self.dashrepl
        elif type == "env":
            dashrepl = self.dashrepl_env
        with open(old_result, "r") as f, open(taxon_result, "w") as w:
            for i in f.readlines():
                new_line = re.sub(r"(name\d+)", dashrepl, i)
                w.write(new_line)

    def run(self):
        super(EnvVifTool, self).run()
        self.env_vif_r()
        self.end()
