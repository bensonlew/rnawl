# -*- coding: utf-8 -*-
# author: qindanhua
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.statistical.correlation import corr_heatmap
import subprocess
import os
import re


class CorHeatmapAgent(Agent):
    """
    CorHeatmapAgent:用于生成之间的correlation
    version: 0.1
    author: wangbixuan
    last_modified: 20171025 by gaohao
    """
    def __init__(self, parent):
        super(CorHeatmapAgent, self).__init__(parent)
        options = [
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "level", "type": "string", "default": "otu"},
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "method", "type": "string", "default": "pearsonr"},
            {"name": "env_cluster", "type": "string", "default": "average"},
            {"name": "species_cluster", "type": "string", "default": "average"},
            {"name": "cor_table", "type": "outfile", "format": "meta.otu.otu_table"},
            {"name": "pvalue_table", "type": "outfile", "format": "meta.otu.otu_table"},
            {"name": "top_species", "type": "int", "default": 0},
            {"name": "abu_trans", "type": "string", "default": "column"},
            {"name": "env_trans", "type": "string", "default": "column"}
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
            raise OptionError("请选择正确的分类水平")
        if not self.option("otutable").is_set:
            raise OptionError('必须提供otu表')
        if self.option("method") not in ["pearsonr", "spearmanr", "kendalltau"]:  # add "kendalltau"(kendall) by zhujuan
            raise OptionError('不支持该相关系数方法')
        self.option('otutable').get_info()
        if self.option('envtable').is_set:
            self.option('envtable').get_info()
            if self.option('envlabs'):
                labs = self.option('envlabs').split(',')
                for lab in labs:
                    if lab not in self.option('envtable').prop['group_scheme']:
                        raise OptionError('该envlabs中的因子不存在于环境因子表：%s' % lab)
            else:
                pass
        else:
            raise OptionError('请选择环境因子表')
        if self.option("top_species") != 0:
            if self.option("top_species") < 2:
                raise OptionError('至少选择两个物种')
        if self.option('envtable').is_set:  # add by zhouxuan 20170720
            env_table = self.option('envtable').prop['path']
            sample_e = []
            with open(env_table, 'r') as e:
                for line in e:
                    line = line.strip('\n').split('\t')
                    if line[0] != '#SampleID':
                        sample_e.append(line[0])
                        for i in range(1, len(line)):
                            if float(line[i]) or line[i] == '0':
                                continue
                            else:
                                raise OptionError('环境因子表中存在分类型环境因子')
            otu_path = self.option("otutable").prop['path']
            with open(otu_path, 'r') as o:
                line = o.readline()
                line = line.strip('\n').split('\t')
                sample_o = line[1:]
            for i in sample_o:
                if i in sample_e:
                    continue
                else:
                    raise OptionError('OTU表中的样本和环境因子表中的样本不一致，请剔除OTU中非法样本！')


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
        super(CorHeatmapAgent, self).end()


class CorHeatmapTool(Tool):
    def __init__(self, config):
        super(CorHeatmapTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.perl_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/"
        self.r_path = self.config.SOFTWARE_DIR + '/program/R-3.3.1/bin/Rscript'
        self.hcluster_script_path = self.config.SOFTWARE_DIR + "/bioinfo/statistical/scripts/"
        self.Rscript_path = self.config.SOFTWARE_DIR + "/program/R-3.3.1/bin/"
        self.cmd_path = '{}/miniconda2/bin/python {}/statistical/pearsonsCorrelation.py'\
            .format(self.config.SOFTWARE_DIR, self.config.PACKAGE_DIR)
        # self.cmd_path=os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/statistical/scripts/pearsonsCorrelation.py')
        self.env_table = self.get_new_env()
        self.real_otu = self.get_otu_table()
        self.name_to_name = {}
        self.env_name = {}

    def get_otu_table(self):
        """
        根据level返回进行计算的otu表路径
        :return:
        """
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
        根据envlabs生成新的envtable
        """
        if self.option('envlabs'):
            new_path = self.work_dir + 'temp_env_table.xls'
            self.option('envtable').sub_group(new_path, self.option('envlabs').split(','))
            return new_path
        else:
            return self.option('envtable').prop['path']

    def run(self):
        """
        运行
        """
        super(CorHeatmapTool, self).run()
        self.run_pearsonsCorrelation()
        self.run_heatmap()
        # self.plot_hcluster()
        self.set_output()
        self.end()

    def run_pearsonsCorrelation(self):
        """
        run pearsonsCorrelation.py
        """
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
            self.set_error('pearsonsCorrelation.py 计算失败')
        self.logger.info('运行pearsonsCorrelation.py计算correlation完成')
        # self.set_output()
        # self.end()

    def get_name(self, table):
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

    def run_heatmap(self):
        line_num = self.get_name(self.work_dir + "/pearsons_correlation_at_%s_level.xls" % self.option('level'))
        if line_num < 2:
            # raise Exception('相关系数矩阵行数/物种数小于2，请尝试切换水平重新运行') #modified by hongdongxuan 20170406
            self.set_error('相关系数矩阵行数/物种数小于2，请尝试切换水平重新运行')
        corr_heatmap(self.work_dir + "/tem.collection.xls", "env_tree.tre", "species_tree.tre",
                     self.option("env_cluster"), self.option("species_cluster"))
        cmd = self.r_path + " run_corr_heatmap.r"
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('heatmap计算成功')
        except subprocess.CalledProcessError:
            self.logger.info('heatmap计算失败')
            self.set_error('heatmap计算失败')
        self.logger.info('生成树文件成功')

    def dashrepl(self, matchobj):
        return self.name_to_name[matchobj.groups()[0]]

    def dashrepl_env(self, matchobj):
        return self.env_name[matchobj.groups()[0]]

    def set_output(self):
        newpath = self.output_dir + "/pearsons_correlation.xls"
        if os.path.exists(newpath):
            os.remove(newpath)
        os.link(self.work_dir + "/pearsons_correlation_at_%s_level.xls" % self.option('level'),
                self.output_dir + "/pearsons_correlation.xls")
        self.option('cor_table', newpath)

        newpath2 = self.output_dir + "/pearsons_pvalue.xls"
        if os.path.exists(newpath2):
            os.remove(newpath2)
        os.link(self.work_dir + "/pearsons_pvalue_at_%s_level.xls" % self.option('level'),
                self.output_dir + "/pearsons_pvalue.xls")
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

    def get_top_specise(self, otu_file):
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
