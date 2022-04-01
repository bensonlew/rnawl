# -*- coding: utf-8 -*-
# author: qindanhua
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.statistical.correlation import corr_heatmap
import subprocess
import os
import re
import pandas as pd


class PearsonsCorrelationAgent(Agent):
    """
    pearsonsCorrelation:用于生成环境因子和otu/taxon之间的correlation
    version: 0.1
    author: wangbixuan
    last_modified: 20160930 by qindanhua
    """
    def __init__(self, parent):
        super(PearsonsCorrelationAgent, self).__init__(parent)
        options = [
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "level", "type": "string", "default": "otu"},
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "envlabs", "type": "string", "default": ""},
            {"name": "method", "type": "string", "default": "pearsonr"},
            {"name": "env_cluster", "type": "string", "default": ""},  # 默认不聚类 by zhujuan 2017.12.19
            {"name": "species_cluster", "type": "string", "default": ""},  # 默认不聚类 by zhujuan 2017.12.19
            {"name": "cor_table", "type": "outfile", "format": "meta.otu.otu_table"},
            {"name": "pvalue_table", "type": "outfile", "format": "meta.otu.otu_table"},
            {"name": "top_species", "type": "int", "default": 0},
            {"name": "project", "type": "string"},
            {"name": "meta_group_name", "type": "string"},
        ]
        self.add_option(options)
        self.step.add_steps('pearsons_correlation')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.pearsons_correlation.start()
        self.step.update()

    def step_end(self):
        self.step.pearsons_correlation.finish()
        self.step.update()

    def gettable(self):
        """
        根据level返回进行计算的矩阵
        """
        if self.option("otutable").format == "meta.otu.tax_summary_dir":
            return self.option("otutable").get_table(self.option('level'))
        else:
            return self.option('otutable').prop['path']

    def check_options(self):
        if self.option("level") not in ['otu', 'domain', 'kindom', 'phylum', 'class', 'order',
                                        'family', 'genus', 'species']:
            raise OptionError("请选择正确的分类水平", code="34100901")
        if not self.option("otutable").is_set:
            raise OptionError('必须提供otu表', code="34100902")
        if self.option("method") not in ["pearsonr", "spearmanr", "kendalltau"]:  # add "kendalltau"(kendall) by zhujuan
            raise OptionError('不支持该相关系数方法', code="34100903")
        self.option('otutable').get_info()
        if self.option('envtable').is_set:
            self.option('envtable').get_info()
            if self.option('envlabs'):
                labs = self.option('envlabs').split(',')
                for lab in labs:
                    if lab not in self.option('envtable').prop['group_scheme']:
                        raise OptionError('该envlabs中的因子不存在于环境因子表：%s', variables=(lab), code="34100904")
            else:
                pass
        else:
            raise OptionError('请选择环境因子表', code="34100905")
        if self.option("top_species") != 0:
            if self.option("top_species") < 2:
                raise OptionError('至少选择两个物种', code="34100906")
        if self.option('envtable').is_set:  # add by zhouxuan 20170720
            env_table = self.option('envtable').prop['path']
            sample_e = []
            with open(env_table, 'r') as e:
                e.next()
                for line in e:
                    line = line.strip('\n').split('\t')
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
                line = line.strip('\n').split('\t')
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
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "PearsonsCorrelation计算结果输出目录"],
            ["./pearsons_correlation_at_'%s'_level.xls" % self.option('level'), "xls", "PearsonsCorrelation矩阵"],
            ["./pearsons_pvalue_at_'%s'_level.xls" % self.option('level'), "xls", "PearsonsCorrelationPvalues"]
        ])
        super(PearsonsCorrelationAgent, self).end()


class PearsonsCorrelationTool(Tool):
    def __init__(self, config):
        super(PearsonsCorrelationTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.perl_path = self.config.SOFTWARE_DIR + "/program/perl/perls/perl-5.24.0/bin/"
        self.r_path = self.config.SOFTWARE_DIR + '/program/R-3.3.1/bin/Rscript'
        self.hcluster_script_path = self.config.SOFTWARE_DIR + "/bioinfo/statistical/scripts/"
        self.Rscript_path = self.config.SOFTWARE_DIR + "/program/R-3.3.1/bin/"
        self.cmd_path = '{}/program/Python/bin/python {}/statistical/pearsonsCorrelation.py'\
            .format(self.config.SOFTWARE_DIR, self.config.PACKAGE_DIR)
        """
        # self.cmd_path=os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/statistical/scripts/pearsonsCorrelation.py')
        #self.env_table = self.get_new_env()
        #self.real_otu = self.get_otu_table()
        """
        self.env_table = ""
        self.real_otu = ""
        self.name_to_name = {}
        self.env_name = {}

    def get_otu_table(self):
        """
        根据level返回进行计算的otu表路径
        :return:
        """
        self.logger.info("获得新的otu表")
        if self.option('otutable').format == "meta.otu.tax_summary_dir":
            otu_path = self.option('otutable').get_table(self.option('level'))
        else:
            otu_path = self.option('otutable').prop['path']
        if self.option("top_species") != 0:
            self.get_top_specise(otu_path)
            otu_path = self.work_dir + "/sorted_otu_file.xls"
        return otu_path

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
        if self.option('project') == "meta":
            env_abund = pd.DataFrame(pd.read_table(env_table, sep='\t', converters = {'#SampleID': str}))
            env_abund.set_index('#SampleID', inplace=True)
        else:
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
            # raise Exception("环境因子：%s的数据都为零，需要剔除该环境因子!" % env_empt)
            self.set_error('环境因子：%s在所选样品中均为0，请剔除该环境因子!', variables=(env_empt), code="34100901")
        new_env_abund.to_csv(new_path, sep="\t", encoding="utf-8")
        return new_path
        #if self.option('envlabs'):
        #    new_path = self.work_dir + 'temp_env_table.xls'
        #    self.option('envtable').sub_group(new_path, self.option('envlabs').split(','))
        #    return new_path
        #else:
        #   return self.option('envtable').prop['path']


    def run(self):
        """
        运行
        """
        super(PearsonsCorrelationTool, self).run()
        self.env_table = self.get_new_env()
        self.real_otu = self.get_otu_table()
        self.run_pearsonsCorrelation()
        if self.option("species_cluster") != "":
            if self.option("env_cluster") != "":
                self.run_heatmap(self.option("env_cluster"),self.option("species_cluster"))
            else:
                self.run_heatmap("average",self.option("species_cluster"))
        elif self.option("env_cluster") != "":
            self.run_heatmap(self.option("env_cluster"),"average")
        else:
            pass
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
        cmd += " %s %s %s %s %s" % (self.real_otu, self.env_table, "./pearsons_correlation_at_'%s'_level.xls" %
                                    self.option('level'), "./pearsons_pvalue_at_'%s'_level.xls" % self.option('level'),
                                    self.option("method"))
        self.logger.info('运行pearsonsCorrelation.py计算correlation')
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('Pearsons Correlation 计算成功')
        except subprocess.CalledProcessError:
            self.logger.info('Pearsons Correlation 计算失败')
            self.set_error('pearsonsCorrelation.py 计算失败', code="34100902")
        self.logger.info('运行pearsonsCorrelation.py计算correlation完成')
        # self.set_output()
        # self.end()

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

    def get_name(self, table):
        self.logger.info("开始对计算的相关性表进行替换名称")
        with open(table, "r") as f, open(self.work_dir + "/tem.collection.xls", "w") as w, \
                open("name_to_name.xls", "w") as nf, open("env_name.xls", "w") as ew:
            first_line = f.readline()
            # w.write(first_line)
            col_names = first_line.strip().split("\t")
            w.write(col_names[0])
            e = 1
            for c in col_names[1:]:
                env_new_name = "colnew"+str(e)
                w.write("\t" + env_new_name)
                self.env_name[env_new_name] = c
                e += 1
            n = 1
            w.write("\n")
            self.logger.info(self.env_name)
            ew.write(str(self.env_name))
            for line in f:
                line = line.split("\t")
                name = line[0]
                new_name = "name"+str(n)
                nf.write(new_name + "\t" + name + "\n")
                self.name_to_name[new_name] = name
                n += 1
                new_line = new_name+"\t"+"\t".join(line[1:])
                w.write(new_line)
        return n

    def run_heatmap(self,env_cluster,species_cluster):
        self.logger.info("开始进行计算heatmap图")
        line_num = self.get_name(self.work_dir + "/pearsons_correlation_at_%s_level.xls" % self.option('level'))
        if line_num < 3:  # modify by zhujuan
            # raise Exception('相关系数矩阵行数/物种数小于2，请尝试切换水平重新运行') #modified by hongdongxuan 20170406
            self.set_error('相关系数矩阵行数/物种数小于2，请尝试切换水平重新运行', code="34100903")
        corr_heatmap(self.work_dir + "/tem.collection.xls", "env_tree.tre", "species_tree.tre",
                     env_cluster, species_cluster)
        cmd = self.r_path + " run_corr_heatmap.r"
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('heatmap计算成功')
        except subprocess.CalledProcessError:
            self.logger.info('heatmap计算失败')
            self.set_error('heatmap计算失败', code="34100904")
        self.logger.info('生成树文件成功')

    def dashrepl(self, matchobj):
        return self.name_to_name[matchobj.groups()[0]]

    def dashrepl_env(self, matchobj):
        return self.env_name[matchobj.groups()[0]]

    def set_output(self):
        self.logger.info("开始链接和检查结果文件")
        if self.option('meta_group_name'):
            if not os.path.exists(self.output_dir + "/" + self.option('meta_group_name')):
                os.makedirs(self.output_dir + "/" + self.option('meta_group_name'))
            newpath = self.output_dir+"/"+self.option('meta_group_name')+"/pearsons_correlation.xls"
        else:
            newpath = self.output_dir + "/pearsons_correlation.xls"
        if os.path.exists(newpath):
            os.remove(newpath)
        os.link(self.work_dir + "/pearsons_correlation_at_%s_level.xls" % self.option('level'),newpath)
        self.option('cor_table', newpath)

        if self.option('meta_group_name'):
            if not os.path.exists(self.output_dir + "/" + self.option('meta_group_name')):
                os.makedirs(self.output_dir + "/" + self.option('meta_group_name'))
            newpath2 = self.output_dir+"/"+self.option('meta_group_name')+"/pearsons_pvalue.xls"
        else:
            newpath2 = self.output_dir + "/pearsons_pvalue.xls"
        if os.path.exists(newpath2):
            os.remove(newpath2)
        os.link(self.work_dir + "/pearsons_pvalue_at_%s_level.xls" % self.option('level'),newpath2)
        self.option('pvalue_table', newpath2)

        species_tree_path = self.work_dir + "/species_tree.tre"
        if os.path.exists(species_tree_path):
            with open(species_tree_path, "r") as f, open(self.work_dir + "/final_species_tree.tre", "w") as w:
                species_tree = f.readline().strip()
                new_species_tree = re.sub(r"(name\d+)", self.dashrepl, species_tree)
                w.write(new_species_tree)

        env_tree_path = self.work_dir + "/env_tree.tre"
        if os.path.exists(env_tree_path):
            with open(env_tree_path, "r") as f, open(self.work_dir + "/final_env_tree.tre", "w") as w:
                env_tree = f.readline().strip()
                new_species_tree = re.sub(r"(colnew\d+)", self.dashrepl_env, env_tree)
                w.write(new_species_tree)

        if self.option('meta_group_name'):
            if os.path.exists(species_tree_path):
                os.link(self.work_dir + "/species_tree.tre",self.output_dir+"/"+self.option('meta_group_name')+"/species_tree.tre")
            if os.path.exists(env_tree_path):
                os.link(self.work_dir + "/final_env_tree.tre",self.output_dir + "/" + self.option('meta_group_name') + "/final_env_tree.tre")
            if os.path.exists(self.work_dir + "/env_name.xls"):
                os.link(self.work_dir + "/env_name.xls",self.output_dir+"/"+self.option('meta_group_name')+"/env_name.xls")
            if os.path.exists(self.work_dir + "/name_to_name.xls"):
                os.link(self.work_dir + "/name_to_name.xls",self.output_dir + "/" + self.option('meta_group_name') + "/name_to_name.xls")
            if os.path.exists(self.work_dir + "/sorted_otu_file.xls"):
                os.link(self.work_dir + "/sorted_otu_file.xls", self.output_dir + "/" + self.option('meta_group_name') + "/sorted_otu_file.xls")


    def get_top_specise(self, otu_file):  # add by zhujuan 解决丰度表中有样品数据都为零的时bug
        abund = pd.DataFrame(pd.read_table(otu_file, sep='\t', index_col=0))
        abund['Col_sum'] = abund.apply(lambda x: x.sum(), axis=1)
        number = int(self.option("top_species"))
        abund_table = abund.sort_values(by=['Col_sum'], ascending=0).head(number)
        del abund_table['Col_sum']
        """
        a = (abund_table > 0).any(axis=0)
        sample_empt = []
        for i in range(len(abund_table.columns)):
            if a[i]:
                pass
            else:
                sample_name = abund_table.columns[i]
                sample_empt.append(sample_name)
        if sample_empt:
                self.set_error("样品：%s的数据都为零，需要剔除该样品或增大总丰度前N的值!" % sample_empt)
        abund_table = abund_table[abund_table.columns[list(a)]]  # 去除数据都为零的样品   需确定信息传到网页，能否不停止运行
        """
        abund_table = abund_table.ix[list((abund_table > 0).any(axis=1))]  # 去除数据都为零的物种/功能/基因
        abund_table_path = self.work_dir + "/sorted_otu_file.xls"
        abund_table.to_csv(abund_table_path, sep="\t", encoding="utf-8")
        """
        specise_sum = {}
        species_list = []
        sorted_spe_list = []
        total = 0
        with open(otu_file, "r") as f:
            f.readline()
            for line in f:
                line = line.strip("\n").split("\t")
                row_sum = sum(map(float, line[1:]))
                specise_sum[line[0]] = row_sum
                species_list.append(line[0])
                total += row_sum
        sorted_species = sorted(specise_sum.items(), key=lambda item: item[1], reverse=True)
        for s in sorted_species[:int(self.option("top_species"))]:
            sorted_spe_list.append(s[0])
        # print sorted_spe_list
        with open(otu_file, "r") as f, open(self.work_dir + "/sorted_otu_file.xls", "w") as w:
            w.write(f.readline())
            for line in f:
                line = line.strip("\n").split("\t")
                if line[0] in sorted_spe_list:
                    w.write("\t".join(line) + "\n")
                    sorted_spe_list.remove(line[0])
        """
