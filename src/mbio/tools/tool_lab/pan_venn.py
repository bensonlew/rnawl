# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
# version 1.0
# last_modify: 20119.05.09

import os
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.graph.venn_table import venn_graph_mg
from mbio.packages.graph.venn_table import venn_graph
from mbio.packages.bac_comp_genome.common_function import link_dir

class PanVennAgent(Agent):
    """
    细菌比较基因组venn图
    """
    def __init__(self, parent):
        super(PanVennAgent, self).__init__(parent)
        options = [
            {"name": "cluster", "type": "infile", "format": "sequence.profile_table"},  # 输入的样本丰度计算的表格
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 输入文件夹
            {"name": "version", "type":"string", "default": ''}, #用于区分新老插件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("cluster"):
            raise OptionError("必须输入丰度表!")
        if not self.option("group_table").is_set:
            raise OptionError("必须输入文件分组表")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(PanVennAgent, self).end()

class PanVennTool(Tool):
    def __init__(self, config):
        super(PanVennTool, self).__init__(config)
        self.R_path = '/program/R-3.3.1/bin/'
        # self.R_path2 = self.config.SOFTWARE_DIR + '/program/R-3.3.1/bin/'  # 循环投递时需要全路径
        # self.venn_path = self.config.SOFTWARE_DIR + '/bioinfo/plot/scripts/'
        self.venn_path = self.config.PACKAGE_DIR + '/graph/scripts/venn_table.py'  # 将scripts的脚本移到packages中使用，凡是涉及的都已改动 add by zhujuan 20171113
        self.venn_path_mg = self.config.PACKAGE_DIR + '/graph/scripts/venn_table_mg.py'
        self.python_path = self.config.SOFTWARE_DIR + '/program/Python/bin/python'
        self.python_path2 = '/program/Python/bin/python'  # 用于小工具
        self.petal_path = self.config.PACKAGE_DIR + '/tool_lab/pan_venn.py'
        self.new_venn_path = self.config.PACKAGE_DIR + '/bac_comp_genome/pan_venn.py'
        self.out = self.work_dir + '/venn_table.xls'
        if os.path.exists(self.out):
            os.remove(self.out)

    def get_group(self):
        """
        获取分组数或样本数
        :return:
        """
        group_table = self.option("group_table").prop['path']
        group_list = []
        sample_list = []
        with open(group_table, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                group_name = line[1]
                sample_name = line[0]
                if group_name not in group_list:
                    group_list.append(group_name)
                if sample_name not in sample_list:
                    sample_list.append(sample_name)
        return (group_list, sample_list)

    def run_venn(self):
        """
        当分组小于等于6的时候画venn图
        """
        otu_table = self.get_cluster_enrich()
        group_file = self.option("group_table").prop['path']
        get_cmd_list = []
        cmd_list = []

        if len(self.option("group_table").prop['group_scheme']) == 1:   # 判断分组方案的个数
            if os.path.exists(self.work_dir + '/group_table'):
                os.remove(self.work_dir + '/group_table')
            os.link(group_file, self.work_dir + '/group_table')  # venn_table的结果与分组文件的目录一致，所以需要将分组文件放在工作目录下
            venn_cmd = '%s %s -i %s -g %s -o cmd.r' % (
            self.python_path, self.venn_path_mg, otu_table, self.work_dir + '/group_table')
            self.logger.info(venn_cmd)
            os.system(venn_cmd)
            self.logger.info('运行venn_cmd')
            cmd = self.R_path + 'Rscript cmd.r'
            self.logger.info("开始运行venn_table")
            self.logger.info(cmd)
            command = self.add_command("venn_table", cmd)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("运行venn_table完成")
            else:
                self.set_error("运行venn_table运行出错!")
            # 统计各组所有otu/物种名 add by qindanhua
            venn_graph_mg(otu_table, group_file, "venn_graph.xls")
            self.remove_other_label(self.work_dir + '/venn_table.xls', self.work_dir + '/venn_graph.xls')
            # self.add_taxon(self.work_dir + '/venn_table2.xls', self.work_dir + '/venn_table_new.xls')  # add_taxon换回原名称
            self.out = self.work_dir + '/venn_table2.xls'
        else:  # 小工具专用，用于批量生成多个分组方案对应的结果
            for i in range(len(self.option("group_table").prop['group_scheme'])):
                select_group = []
                sample_dir = self.work_dir + '/' + self.option("group_table").prop['group_scheme'][i]
                os.mkdir(sample_dir)
                select_group.append(self.option("group_table").prop['group_scheme'][i])
                self.option('group_table').sub_group(sample_dir + '/venn_group_' + str(i + 1), select_group)
                self.option("otu_table").get_table_of_main_table(otu_table, sample_dir + '/input_' + str(i + 1),
                                                                 group_file)
                venn_cmd = '%spython %s -i %s -g %s -o %scmd_%s.r' % (
                self.python_path2, self.venn_path, sample_dir + '/input_' + str(i + 1),
                sample_dir + '/venn_group_' + str(i + 1), sample_dir + '/', i + 1)
                get_cmd_list.append(venn_cmd)  # 存放所有生成cmd.r的命令
                cmd_list.append(self.R_path + 'Rscript {}cmd_{}.r'.format(sample_dir + '/', i + 1))  # 存放所有运行cmd.r的命令
            #  循环投递，批量生成cmd.r文件，结果及日志存放在对应分组方案的文件夹下
            n = 0
            for i in get_cmd_list:
                command = self.add_command("get_cmd_{}".format(n), i).run()
                self.wait(command)
                if command.return_code == 0:
                    self.logger.info("运行{}完成".format(command.name))
                    n += 1
                else:
                    self.set_error("运行%s运行出错!", variables=(command.name))

            # 循环投递，批量运行cmd.r文件，结果及日志存放在对应分组方案的文件夹下
            n = 0
            for i in cmd_list:
                command = self.add_command("cmd_{}".format(n), i).run()
                self.wait(command)
                if command.return_code == 0:
                    self.logger.info("运行{}完成".format(command.name))
                else:
                    self.set_error("运行%s运行出错!", variables=(command.name), code="32301203")
                # 统计各组所有otu/物种名
                venn_graph(otu_table, self.work_dir + '/' + self.option("group_table").prop['group_scheme'][
                    n] + '/venn_group_' + str(n + 1),
                           self.work_dir + '/' + self.option("group_table").prop['group_scheme'][n] + "/venn_graph.xls")
                n += 1

    def run_petal(self):
        """
        当分组大于等于6或者当分组等于1且样本数大于等于6的时候画花瓣图
        :return:
        """
        otu_table = self.option("cluster").prop['path']
        cmd = "{} {} -i {} -r {} -out {}".format(self.python_path2, self.petal_path, otu_table,
                                               self.option('group_table').prop['path'], self.out)
        command = self.add_command("run_input", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行venn完成！")
        else:
            self.set_error("运行venn出错!")

    def get_cluster_enrich(self):
        """
        根据输入的cluster表格整理成丰度表
        :return:
        """
        venn_input = os.path.join(self.work_dir, "pan_venn.xls")
        cluster = self.option("cluster").prop["path"]
        with open(cluster, 'r') as f, open(venn_input, 'w') as w:
            lines = f.readlines()
            head = lines[0].strip().split("\t")
            venn_header = head[0] + "\t" + "\t".join(head[3:])
            w.write(venn_header + "\n")
            for line in lines[1:]:
                line = line.strip().split("\t")
                cluster_name = line[0]
                sample_cluster_num_list = []
                for i in range(3, len(head[3:])+3):
                    if line[i] != '-':
                        sample_gene_list = line[i].split(',')
                        sample_gene_num = len(sample_gene_list)
                    else:
                        sample_gene_num = 0
                    sample_cluster_num_list.append(str(sample_gene_num))
                w.write("{}\t{}\n".format(cluster_name, '\t'.join(sample_cluster_num_list)))
        return venn_input


    def run_new_venn(self):
        """
        针对插件升级重新导表的数据准备的结果
        :return:
        """
        otu_table = self.option("cluster").prop['path']
        cmd = "{} {} -i {} -r {} -out {}".format(self.python_path, self.new_venn_path, otu_table,
                                               self.option('group_table').prop['path'], self.out)
        command = self.add_command("run_venn", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行venn完成！")
        else:
            self.set_error("运行venn出错!")

    def remove_other_label(self,infile, graph_infile):
        """
        去掉之前共有的部分
        """
        ori_venn = infile
        group_ele = graph_infile
        group_map = {}
        with open(group_ele, 'r') as f:
            f.readline()
            for line in f:
                spline = line.strip().split('\t')
                group_map[spline[0]] = set(spline[1].split(';'))
        groups = set(group_map.keys())

        new_venn = self.work_dir + '/venn_table2.xls'
        with open(new_venn,'w') as fw, open(ori_venn, 'r') as fr:
            for line in fr:
                spline = line.strip().split('\t')
                if '&' in spline[0]:
                    cg = [i.strip() for i in spline[0].split('&')]
                    cg = set([i.split(" only")[0] for i in cg])
                    other_g = groups - cg
                    if len(spline) == 3:
                        cur_ele = set(spline[2].split(';'))
                    else:
                        cur_ele = []
                    for og in other_g:
                        cur_ele = set(cur_ele) - set(group_map[og])
                    # s_cur_ele = sorted(cur_ele,key=lambda x: int(x.split('_')[1]))
                    s_cur_ele = sorted(cur_ele)
                    fw.write(spline[0]+'\t'+str(len(s_cur_ele))+'\t'+';'.join(s_cur_ele)+'\n')
                else:
                    fw.write(line)

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("正在生成结果文件目录")
        group_list, sample_list = self.get_group()
        if len(group_list) >= 6:
            if os.path.exists(self.output_dir + "/new_venn_graph.xls"):
                os.remove(self.output_dir + "/new_venn_graph.xls")
            if os.path.exists(self.out):
                os.link(self.out, self.output_dir + "/new_venn_graph.xls")
        else:
            if len(group_list) == 1 and group_list[0] in ['all', 'All', 'ALL'] and len(sample_list) >= 6:
                if os.path.exists(self.output_dir + "/new_venn_graph.xls"):
                    os.remove(self.output_dir + "/new_venn_graph.xls")
                if os.path.exists(self.out):
                    os.link(self.out, self.output_dir + "/new_venn_graph.xls")
            else:
                if os.path.exists(self.output_dir + "/new_venn_graph.xls"):
                    os.remove(self.output_dir + "/new_venn_graph.xls")
                if os.path.exists(self.out):
                    os.link(self.work_dir + "/new_venn_graph.xls", self.output_dir + "/new_venn_graph.xls")

    def run(self):
        super(PanVennTool, self).run()
        if self.option("version") == '':
            self.run_new_venn()
        else:
            group_list, sample_list = self.get_group()
            if len(group_list) >= 6:
                self.run_petal()
            else:
                if len(group_list) == 1 and group_list[0] in ['all', 'All', 'ALL'] and len(sample_list) >= 6:
                    self.run_petal()
                else:
                    self.run_venn()
        self.set_output()
        self.end()