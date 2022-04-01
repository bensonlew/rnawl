# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os
import re
from biocluster.workflow import Workflow
from mbio.packages.metaasv.common_function import link_file


class PlotTreeWorkflow(Workflow):
    """
    metaasv 系统进化树分析
    筛选原则，筛选注释的分类学水平丰度为第一的代表序列构建进化树
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PlotTreeWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "level", "type": 'int', "default": 9},##分类水平
            {"name": "topN", "type": 'int', "default": 0},## 选择丰度前50的物种
            {"name": "asv_id", "type": 'string', "default": ''},
            {"name": "main_id", "type": 'string', "default": ''},
            {"name": "params", "type": 'string', "default": ''},
            {"name": "group_id", "type": 'string'},
            {"name": "task_id", "type": 'string'},
            {"name": "update_info", "type": 'string'},
            {"name": "group_detail", "type": 'string', "default": ""},
            {"name": "color_level_id", "type": 'int', "default": 0},##to_file用
            {"name": "sample_group", "type": "infile", "format": "meta.otu.group_table"},##group表
            {"name": "method", "type": "string", "default": "NJ"},  #NJ, MP, ML计算进化树的方法
            {"name": "run_tree", "type": "string" , "default": "part"},  # 是否挑选进化树
            {"name": "otu_seq", "type": "infile", "format": "meta.fasta" },##输入的otu代表序列

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.tree_type_map = {
            'NJ' : 'phylo',
            'ML' : 'phylo_ml',
            'MP' : 'phylo_mp'
        }
        self.tree_type = self.tree_type_map[self.option('method')]

    def check(self):
        """
        参数检查
        :return:
        """
        if self.option("level") >= self.option("color_level_id"):
            self.set_error("颜色设置id必须大于选择的分类水平")
        if not  self.option('otu_seq').is_set:
            self.set_error('run_tree 必须提供otu_seq 文件')

    def run(self):
        """
        运行
        :return:
        """
        self.run_tree()
        super(PlotTreeWorkflow, self).run()

    def run_tree(self):
        """
        从头构建进化树
        :return:
        """
        self.phly_tree_tool = self.add_tool('graph.phy_tree')
        method_map = {
            'ML' : 'Maximum_Likehood(ML)',
            'MP' : 'Maximum_Parsimony(MP)',
            'NJ' : 'Neighbor_Joining(NJ)'
        }
        self.change_name() ##完成改名
        self.logger.info(self.fasta)
        method = method_map[self.option('method')]
        opts =  {
            'method' : method,
            'fasta' : self.fasta,
            'align_method' : 'mafft',
            'sequence_type':'no_coding'
        }
        if self.option('method') == 'ML':
            opts['tree_software'] = 'iqtree'

        if self.option('method') == 'MP':
            opts['bootstrap'] = 100
        else:
            opts['bootstrap'] = 500

        self.phly_tree_tool.set_options(opts)
        if self.option('run_tree') == 'all':
            self.phly_tree_tool.on('end',self.all_tree_api)
        else:
            self.phly_tree_tool.on('end',self.set_db)

        self.phly_tree_tool.run()

    def change_name(self):
        """
        目的是将fasta文件的名称进行替换，主要是减少名称中含有的特殊字符导致的错误
        解决mega运行得到的进化树少了一个下划线导致名称报错的问题
        :return:
        """
        self.fasta = self.work_dir + '/all.fasta'
        self.name_map = {} ##名称与改写后的名称的对应关系
        seq_id = 1
        self.ref_names = [] ##所有改写前的名称list
        otu_path = self.option('otu_seq').prop['path']

        fwr = open(self.fasta, 'w')
        with open(otu_path, 'r') as f:
            for line in f:
                if line[0] == '>':
                    line = line.strip().split(">")
                    name = line[1]
                    new_name = 'seq'+str(seq_id)+"n"
                    fwr.write('>'+new_name+"\n")
                    self.ref_names.append(new_name)
                    self.name_map[new_name] = name
                    seq_id += 1
                else:
                    fwr.write(line)
        fwr.close()

    def change_back_name(self, output_tree, new_output_tree):
        """
        针对软件的结果进行替换为原来的名称
        :return:
        """
        with open(output_tree, 'r') as f, open(new_output_tree, 'w') as w:
            line = f.readline()
            new_line = ''
            for sp in self.ref_names:
                if re.search(r"%s"%sp, line):
                    if sp in self.name_map.keys():
                        new_sp = self.name_map[sp]
                        if new_line == "":
                            new_line = line
                        else:
                            new_line = new_line
                        new_line = new_line.replace('%s'%sp, '%s'%new_sp, 1)
            w.write(new_line)

    def set_db(self):
        """
        连接结果文件、替换名称和导入MongoDB数据
        """
        output_otu = self.output_dir + '/species_table.xls'
        output_tree = self.output_dir + '/phylo_tree.tre'
        output_species = self.output_dir + '/species_group.xls'
        link_file(self.phly_tree_tool.output_dir + '/phylo_tree.nwk', output_tree)
        link_file(self.work_dir + '/species_group.xls', output_species)
        link_file(self.work_dir + '/species_table.xls', output_otu)
        new_output_tree = self.output_dir + '/phylo_tree.tre_new'
        self.change_back_name(output_tree, new_output_tree)
        os.remove(output_tree)
        os.rename(new_output_tree, output_tree)
        api_tree = self.api.api("metaasv.phylo_tree")
        api_tree.add_phylo_tree_info(self.option('main_id'), seq_num=self.option('topN'))  #seq_num 作用是给前端确定展示的窗口
        self.end()

    def end(self):
        """
        上传结果文件夹
        :return:
        """
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "进化分析结果目录", 0, ""],
            ["species_table.xls", "txt", "物种样本统计表", 0, ""],
            ["phylo_tree.tre", "tree", "进化树", 0, ""],
            ["species_group.xls", "txt", "物种在高层级的分类表", 0, ""]
        ])
        super(PlotTreeWorkflow, self).end()