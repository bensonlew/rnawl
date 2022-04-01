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
    last_modified: 20171025 by gaohao
    """
    def __init__(self, parent):
        super(CorHeatmapAgent, self).__init__(parent)
        options = [
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "method", "type": "string", "default": "pearsonr"},
            {"name": "env_cluster", "type": "string", "default": ""},
            {"name": "species_cluster", "type": "string", "default": ""},
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

    def check_options(self):
        if not self.option("otutable").is_set:
            raise OptionError('必须提供数据表表')
        if self.option("method") not in ["pearson", "spearman", "kendalltau"]:
            raise OptionError('不支持该相关系数方法')
        if self.option("top_species") != 0:
            if self.option("top_species") < 2:
                raise OptionError('至少选择两个物种')
        if not self.option('envtable').is_set:
            raise OptionError('必须提供因子文件')
        
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
            ["./pearsons_correlation.xls" , "xls", "PearsonsCorrelation矩阵"],
            ["./pearsons_pvalue.xls", "xls", "PearsonsCorrelationPvalues"]
        ])
        super(CorHeatmapAgent, self).end()


class CorHeatmapTool(Tool):
    def __init__(self, config):
        super(CorHeatmapTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.perl_path = self.config.SOFTWARE_DIR + "/program/perl/perls/perl-5.24.0/bin/"
        self.r_path = self.config.SOFTWARE_DIR + '/program/R-3.3.1/bin/Rscript'
        self.hcluster_script_path = self.config.SOFTWARE_DIR + "/bioinfo/statistical/scripts/"
        self.Rscript_path = self.config.SOFTWARE_DIR + "/program/R-3.3.1/bin/"
        self.cmd_path = '{}/program/Python/bin/python {}/statistical/pearsonsCorrelation.py'\
            .format(self.config.SOFTWARE_DIR, self.config.PACKAGE_DIR)
        # self.cmd_path=os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/statistical/scripts/pearsonsCorrelation.py')
        self.real_otu = self.get_otu_table()
        self.env_table = self.get_env_table()
        self.name_to_name = {}
        self.env_name = {}

    def get_otu_table(self):
        """
        根据level返回进行计算的otu表路径
        :return:
        """
        if self.option('abu_trans') in ['column']:
            otu_path = self.option('otutable').prop['path']
        else:
            otu_path = self.work_dir + "/otu_t_table.xls"
            self.t_table(self.option('otutable').prop['path'],otu_path )
        if self.option("top_species") != 0:
            self.get_top_specise(otu_path)
            otu_path = self.work_dir + "/sorted_otu_file.xls"
        return otu_path

    def get_env_table(self):
        """
        根据level返回进行计算的otu表路径
        :return:
        """
        if self.option('env_trans') in ['column']:
            env_path = self.option('envtable').prop['path']
        else:
            env_path = self.work_dir + "/new_env_table.xls"
            self.t_table(self.option('envtable').prop['path'], env_path)
        sample_e = []
        with open(env_path, 'r') as e:
            for line in e.readlines()[1:]:  # by zhaozhigang 20200918
                line = line.strip('\n\r').split('\t')
                if line[0] != '#SampleID':
                    sample_e.append(line[0])
                    for i in range(1, len(line)):
                        if float(line[i]) or line[i] == '0':
                            continue
                        else:
                            raise OptionError('环境因子表中存在分类型环境因子')
        with open(self.real_otu, 'r') as o:
            line = o.readline()
            line = line.strip('\n\r').split('\t')
            sample_o = line[1:]
        sample_e.sort()
        sample_o.sort()
        for i in sample_o:
            if i in sample_e:
                continue
            else:
                raise OptionError('数据表中的样本和环境因子表中的样本不一致')
        return env_path

    def run(self):
        """
        运行
        """
        super(CorHeatmapTool, self).run()
        self.run_pearsonsCorrelation()
        self.run_heatmap()
        self.set_output()
        self.end()

    def run_pearsonsCorrelation(self):
        """
        run pearsonsCorrelation.py
        """
        if self.option("method") in ['pearson']:
            method = 'pearsonr'
        elif self.option("method") in ['spearman']:
            method = 'spearmanr'
        else:
            method = self.option("method")
        cmd = self.cmd_path
        cmd += " %s %s %s %s %s" % (self.real_otu, self.env_table, "./pearsons_correlation.xls"
                                    , "./pearsons_pvalue.xls" ,method)
        self.logger.info('运行pearsonsCorrelation.py计算correlation')
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('Pearsons Correlation 计算成功')
        except subprocess.CalledProcessError:
            self.logger.info('Pearsons Correlation 计算失败')
            self.set_error('pearsonsCorrelation.py 计算失败')
        self.logger.info('运行pearsonsCorrelation.py计算correlation完成')


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
        line_num = self.get_name(self.work_dir + "/pearsons_correlation.xls")
        if line_num < 2:
            # raise Exception('相关系数矩阵行数/物种数小于2，请尝试切换水平重新运行') #modified by hongdongxuan 20170406
            self.set_error('相关系数矩阵行数/物种数小于2，请尝试切换水平重新运行')
        if self.option("species_cluster") == "" and self.option("env_cluster") == "": # modified by qingchen.zhang @20200219
            pass
        else:
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
        os.link(self.work_dir + "/pearsons_correlation.xls",
                self.output_dir + "/pearsons_correlation.xls")

        newpath2 = self.output_dir + "/pearsons_pvalue.xls"
        if os.path.exists(newpath2):
            os.remove(newpath2)
        os.link(self.work_dir + "/pearsons_pvalue.xls",
                self.output_dir + "/pearsons_pvalue.xls")

        species_tree_path = self.work_dir + "/species_tree.tre"
        if os.path.exists(species_tree_path):
            with open(species_tree_path, "r") as f, open(self.work_dir + "/final_species_tree.tre", "w") as w:
                species_tree = f.readline().strip()
                new_species_tree = re.sub(r"(name\d+)", self.dashrepl, species_tree)
                w.write(new_species_tree)
            os.link(self.work_dir + "/final_species_tree.tre",self.output_dir + "/final_species_tree.tre")
        env_tree_path = self.work_dir + "/env_tree.tre"
        if os.path.exists(env_tree_path):
            with open(env_tree_path, "r") as f, open(self.work_dir + "/final_env_tree.tre", "w") as w:
                env_tree = f.readline().strip()
                new_species_tree = re.sub(r"(colnew\d+)", self.dashrepl_env, env_tree)
                w.write(new_species_tree)
            os.link(self.work_dir + "/final_env_tree.tre", self.output_dir + "/final_env_tree.tre")

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

    def t_table(self, table_file, new_table):  # 表格转置
        """
    	转换颠倒表格内容
    	"""
        with open(table_file) as f, open(new_table, 'w') as w:
            table_list = [i.rstrip().split('\t') for i in f.readlines()]
            table_list = map(lambda *a: '\t'.join(a) + '\n', *table_list)
            w.writelines(table_list)
