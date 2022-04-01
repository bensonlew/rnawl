# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last_modify: 2018.04.19

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import os


class MicropitaAgent(Agent):
    """
    Micropita挑选样品的工具
    """

    def __init__(self, parent):
        super(MicropitaAgent, self).__init__(parent)
        options = [
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table,meta.otu.tax_summary_dir"},
            {"name": "level", "type": "string", "default": "otu"},
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "target_feature", "type": "string", "default": ""},  # 筛选某些特定物种，物种间以“；”连接
            {"name": "target_metrics", "type": "string", "default": "rank"},  # target_feature的度量方法，[rank/abundance]
            {"name": "distance_method", "type": "string", "default": "braycurtis"},  # 距离算法,进化距离待确定
            {"name": "diversity_index", "type": "string", "default": "simpson"},  # 多样性指数
            {"name": "filter_nu", "type": "int", "default": 10},  # 挑选的样品量
            {"name": "phy_newick", "type": "infile", "format": "meta.beta_diversity.newick_tree"},
            # 当distance_method 为unifrac时需要用到
        ]
        self.add_option(options)

    def gettable(self):
        """
        根据level返回进行计算的丰度表
        :return:
        """
        if self.option('otutable').format == "meta.otu.tax_summary_dir":
            return self.option('otutable').get_table(self.option('level'))
        else:
            return self.option('otutable').prop['path']

    def check_options(self):
        """
        检查参数
        """
        if self.option("level") not in ['otu', 'domain', 'kindom', 'phylum', 'class', 'order',
                                        'family', 'genus', 'species']:
            raise OptionError("请选择正确的分类水平", code="32704501")
        if not self.option("otutable").is_set:
            raise OptionError("请传入丰度文件！", code="32704502")
        else:
            self.option("otutable").get_info()
            if self.option('otutable').prop['sample_num'] <= self.option('filter_nu'):
                raise OptionError('丰度表的样本数目少于需要挑选的样品量，不可进行Micropita分析', code="32704503")
        if self.option("group_file").is_set:
            group = self.option("group_file").get_group_spname()
            for i in group.keys():
                if len(group[i]) <= self.option('filter_nu'):
                    raise OptionError("组内样品数小于需要挑选的样品量，不可进行Micropita分析！", code="32704504")
                    return
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(MicropitaAgent, self).end()


class MicropitaTool(Tool):
    def __init__(self, config):
        super(MicropitaTool, self).__init__(config)
        self.micropita = self.config.PACKAGE_DIR + "/meta/scripts/"
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/meta/micropita/")
        self.r_path = "/program/R-3.3.1/bin/"
        self.perl_path = "/program/perl-5.24.0/bin/perl"

    def run_micropita(self):
        with open(self.option("otutable").prop['path'], "r") as f, open(self.work_dir + "/input_abund.xls", "w") as g:
            head = f.readline().split("\t", 1)
            head[0] = "OTUID"
            g.write(head[0] + "\t" + head[1])
            for line in f:
                g.write(line)
        cmd = '{} {}microPITA.pl -a {} -b {}  -f  {} -label OTUID  -o  {} -n {}'.format(self.perl_path, self.micropita,
                                                                                        self.option("diversity_index"),
                                                                                        self.option("distance_method"),
                                                                                        self.work_dir + "/input_abund.xls",
                                                                                        self.work_dir,
                                                                                        self.option('filter_nu'))
        if self.option("group_file").is_set:  # 软件需要的group文件与常规的不一样，需要转置处理
            group_file = pd.DataFrame(pd.read_table(self.option("group_file").prop['path'], sep='\t')).T
            group_file.index = ["#sample","group"]
            group_file.to_csv("./m_group.xls", sep='\t', header=None)
            cmd += " -g " + self.work_dir + "/m_group.xls"
        if self.option("target_feature"):
            feature = self.option("target_feature").split(";")
            with open(self.work_dir + "/target.xls", "w") as g:
                for i in feature:
                    g.write(i + "\n")
            cmd += " -target " + self.work_dir + "/target.xls" + " -r " + self.option("target_metrics")
        command = self.add_command("micropita", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("microPITA运行完成")
            cmd1 = self.r_path + "Rscript pita.cmd.r"
            command = self.add_command("pita", cmd1).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("microPITA运行完成")
                self.micropita_tidy()
            else:
                self.set_error("microPITA运行出错!", code="32704501")
        else:
            self.set_error("microPITA运行出错!", code="32704502")

    def micropita_tidy(self):
        tidy = {}
        specimens = {}
        if self.option("group_file").is_set:
            with open(self.work_dir + "/m_group.xls", "r") as g:
                specimen = g.readline().strip().split("\t")[1:]
                group = g.readline().strip().split("\t")[1:]
                for i in range(len(specimen)):
                    specimens[specimen[i]] = group[i]
        diverse = self.work_dir + "/diverse.xls"
        if os.path.exists(diverse):
            with open(diverse, "r") as f:
                sample = f.readline().strip().split("\t")[1:]
                if len(sample) > 0:
                    tidy["Maximum diversity"] = sample
        extreme = self.work_dir + "/extreme.xls"
        if os.path.exists(extreme):
            with open(extreme, "r") as f:
                sample = f.readline().strip().split("\t")[1:]
                if len(sample) > 0:
                    tidy["Most dissimilar"] = sample
        representative = self.work_dir + "/representative.xls"
        if os.path.exists(representative):
            with open(representative, "r") as f:
                sample = f.readline().strip().split("\t")[1:]
                if len(sample) > 0:
                    tidy["Most representative"] = sample
        feature = self.work_dir + "/feature.xls"
        if os.path.exists(feature):
            with open(feature, "r") as f:
                sample = f.readline().strip().split("\t")[1:]
                if len(sample) > 0:
                    tidy["Targeted feature"] = sample
        discriminant = self.work_dir + "/discriminant.xls"
        if os.path.exists(discriminant):
            with open(discriminant, "r") as f:
                sample = f.readline().strip().split("\t")[1:]
                for i in sample:
                    if ("Discriminant " + specimens[i]) in tidy.keys():
                        tidy["Discriminant " + specimens[i]].append(i)
                    else:
                        tidy["Discriminant " + specimens[i]] = []
                        tidy["Discriminant " + specimens[i]].append(i)
        distinct = self.work_dir + "/distinct.xls"
        if os.path.exists(distinct):
            with open(distinct, "r") as f:
                sample = f.readline().strip().split("\t")[1:]
                for i in sample:
                    if ("Distinct " + specimens[i]) in tidy.keys():
                        tidy["Distinct " + specimens[i]].append(i)
                    else:
                        tidy["Distinct " + specimens[i]] = []
                        tidy["Distinct " + specimens[i]].append(i)
        tidy_xls = pd.DataFrame(tidy).T
        tidy_xls.index.name = "Method"
        sample_col = []
        for i in range(self.option('filter_nu')):
            sample_col.append("Sample" + str(i + 1))
        tidy_xls.columns = sample_col
        tidy_xls.to_csv(self.output_dir + "/microPITA_select.xls", sep="\t")
        os.link(self.work_dir + "/pcoa.xls", self.output_dir + "/pcoa.xls")

    def run(self):
        super(MicropitaTool, self).run()
        self.run_micropita()
        self.end()
