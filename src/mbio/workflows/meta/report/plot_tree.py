# -*- coding: utf-8 -*-
# __author__ = 'shenghe'

""""""

import os
import re
from biocluster.workflow import Workflow
# from bson import ObjectId
from mbio.packages.beta_diversity.filter_newick import get_level_newicktree
# import datetime
# import json
from collections import defaultdict
import time
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class PlotTreeWorkflow(Workflow):
    """
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PlotTreeWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "level", "type": 'int', "default": 9},
            {"name": "topN", "type": 'int', "default": 0},
            {"name": "otu_id", "type": 'string', "default": ''},
            {"name": "main_id", "type": 'string', "default": ''},
            {"name": "params", "type": 'string', "default": ''},
            {"name": "group_id", "type": 'string'},
            {"name": "update_info", "type": 'string'},
            {"name": "group_detail", "type": 'string', "default": ""},
            {"name": "color_level_id", "type": 'int', "default": 0},
            {"name": "sample_group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "method", "type": "string", "default": "NJ"},  #NJ, MP, ML
            {"name": "task_id", "type": "string"},
            {"name": "run_tree", "type": "string" , "default": "part"},  # add v4 201910
            {"name": "otu_seq", "type": "infile", "format": "meta.fasta" },
            {"name": "ori_otu_id", "type": 'string', "default": ''}

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False

    def check(self):
        if self.option("level") >= self.option("color_level_id"):
            self.set_error("颜色设置id必须大于选择的分类水平", code="12703501")
        if not self.option("otu_id"):
            self.set_error("必须提供OTU的主表id", code="12703502")

        if self.option('run_tree') in ['all','part']:
            if not  self.option('otu_seq').is_set:
                self.set_error('run_tree 必须提供otu_seq 文件', code="12703504")


    def run(self):
        self.tree_type_map = {
            'NJ' : 'phylo',
            'ML' : 'phylo_ml',
            'MP' : 'phylo_mp'
        }
        self.tree_type = self.tree_type_map[self.option('method')]

        if self.option('run_tree') in ['all','part'] :
            self.run_tree()
            super(PlotTreeWorkflow, self).run()
        else:
            self.start_listener()
            self.fire("start")
            self.run_pick_pip()


    def run_tree(self):  #v4 add
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


    def all_tree_api(self):  #v4 add
        api_tree = self.api.phylo_tree
        self.tree_file = self.phly_tree_tool.output_dir + '/phylo_tree.nwk'
        self.logger.info('开始导入文件：%s' %self.tree_file)
        api_tree.add_phylo_tree_for_report(self.option('ori_otu_id'),self.option("task_id"), self.tree_type ,self.tree_file)
        self.logger.info('运行完all tree。开始运行pick tree')
        self.run_pick_pip()

    def run_pick_pip(self):
        self.species = []
        otu_format = self.work_dir + '/format_otu_table.xls'
        species_format = self.work_dir + '/species_group.xls'
        tree_file = self.work_dir + '/format.tre'
        group = self.sample_in_group()
        if self.option("group_id") not in ['all','All','ALL',None]:
            self.get_newicktree(tree_file,group)
        else:
            self.get_newicktree(tree_file)
        if self.option('color_level_id'):
            if self.option("group_id") not in ['all', 'All', 'ALL', None]:
                self.format_group_otu_table(otu_format, species_format)
            else:
                self.format_otu_table(otu_format, species_format)
        else:
            if self.option("group_id") not in ['all', 'All', 'ALL', None]:
                self.format_group_otu_table(otu_format)
            else:
                self.format_otu_table(otu_format)

        self.set_db()

    def sample_in_group(self):                 #增加函数，提取分组中的样品名称为列表，传递给get_netwicktree
        with open(self.option("sample_group").path) as g:
            g.readline()
            group = {}
            for i in g:
                split_i = i.strip().split('\t')
                group[split_i[0]] = split_i[1]
            sample_list = list(group.keys())
        return sample_list

    def format_group_otu_table(self, out_otu_file, out_species_group_file=None):
        """
        """
        species_dict = defaultdict(list)
        species_index = self.option("color_level_id") - 1
        with open(self.option('otu_table').path) as f, open(self.option("sample_group").path) as g, open(out_otu_file, 'w') as w:
            g.readline()
            group = {}  # 样本分组
            for i in g:
                split_i = i.strip().split('\t')
                group[split_i[0]] = split_i[1]
            group_names = list(set(group.values()))
            group_index = dict(zip(group_names, [[] for i in range(len(group_names))]))  # 样本index
            group_value = dict(zip(group_names, [[] for i in range(len(group_names))]))  # 包含样本值列表
            all_sample = f.readline().rstrip().split('\t')[1:]
            for m, n in enumerate(all_sample):
                group_index[group[n]].append(m + 1)
            w.write('ID\t' + '\t'.join(group_names) + '\n')
            for i in f:
                line_split = re.split('\t', i.strip())
                name_split = line_split[0].split(';')
                new_name = name_split[-1].strip().replace(':', '-')  # 协同树去除名称中有冒号的样本名
                if new_name not in self.leaves:
                    continue
                # self.species.append(new_name)
                if out_species_group_file:
                    species_dict[name_split[species_index]].append(new_name)
                for key, indexs in group_index.iteritems():
                    for index in indexs:
                        group_value[key].append(int(line_split[index]))
                new_line = new_name
                for group_name in group_names:
                    new_line += '\t' + str(sum(group_value[group_name]))
                    group_value[group_name] = []
                w.write(new_line + '\n')
            if out_species_group_file:
                with open(out_species_group_file, 'w') as w:
                    w.write('#species/OTU\tGROUP\n')
                    for m, n in species_dict.iteritems():
                        m = m.strip()
                        for i in n:
                            w.write(i + '\t' + m + '\n')


    def change_name(self):
        """
        目的是将fasta文件的名称进行替换，主要是减少名称中含有的特殊字符导致的错误
        解决mega运行得到的进化树少了一个下划线导致名称报错的问题
        add by qingchen.zhang
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
        """
        output_otu = self.output_dir + '/species_table.xls'
        output_tree = self.output_dir + '/phylo_tree.tre'
        output_species = self.output_dir + '/species_group.xls'
        if os.path.exists(output_otu):
            os.remove(output_otu)
        if os.path.exists(output_tree):
            os.remove(output_tree)
        if os.path.exists(output_species):
            os.remove(output_species)

        if self.option('run_tree') == 'part':
            os.link(self.work_dir + '/species_table.xls', self.output_dir + '/species_table.xls')
            os.link(self.phly_tree_tool.output_dir + '/phylo_tree.nwk', output_tree)
        else:
            os.link(self.work_dir + '/format_otu_table.xls', self.output_dir + '/species_table.xls')
            ##os.link(self.work_dir + '/format.tre', self.output_dir + '/phylo_tree.tre')
            os.link(self.work_dir + '/temp_filter_newick.tre.bootstrap', self.output_dir + '/phylo_tree.tre')  #guanqing.zou 20180420
        if os.path.exists(self.work_dir + '/species_group.xls'):
            os.link(self.work_dir + '/species_group.xls', self.output_dir + '/species_group.xls')
        # if self.phly_tree_tool.option("tree_software") in ["megacc"]:
        new_output_tree = self.output_dir + '/phylo_tree.tre_new'
        self.change_back_name(output_tree, new_output_tree)
        os.system("rm %s" %(output_tree))
        os.system("mv %s %s" %(new_output_tree, output_tree))
        api_tree = self.api.phylo_tree
        api_tree.add_phylo_tree_info(self.option('main_id'), seq_num=self.option('topN'))  #seq_num 作用是给前端确定展示的窗口
        #self.end()
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "sg_phylo_tree")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "phylo_tree",
                "interaction": 1,
                "main_table": "sg_phylo_tree",
            })
            self.figsave.run()
        else:
            self.end()


    def format_otu_table(self, out_otu_file, out_species_group_file=None):
        """
        配合进化树，拆分otu表名称
        """
        species_dict = defaultdict(list)
        species_index = self.option("color_level_id") - 1
        with open(self.option("otu_table").path) as f, open(out_otu_file, 'w') as w:
            w.write(f.readline())
            for i in f:
                line_split = re.split('\t', i, maxsplit=1)
                name_split = line_split[0].split(';')
                new_name = name_split[-1].strip().replace(':', '-')
                if new_name not in self.leaves:
                    continue
                #if sum([int(i) for i in line_split[1].strip().split('\t')]) == 0:  # hesheng 20161115 去除所有样本为0的情况，目前to_file没有相关功能，暂时添加
                    # raise Exception('存在全部物种/OTU代表序列数量都为0的情况')
                    #continue
                self.species.append(new_name)
                if out_species_group_file:
                    species_dict[name_split[species_index]].append(new_name)
                w.write(new_name + '\t' + line_split[1])
        if out_species_group_file:
            with open(out_species_group_file, 'w') as w:
                w.write('#species/OTU\tGROUP\n')
                for m, n in species_dict.iteritems():
                    m = m.strip()
                    for i in n:
                        w.write(i + '\t' + m + '\n')

    def get_newicktree(self, output_file,group=None):   #增加group参数

        tree = get_level_newicktree(self.option("otu_id"),
                                    group = group,
                                    level=self.option('level'),
                                    tempdir=self.work_dir,
                                    topN=self.option('topN'),
                                    tree_type= self.tree_type)
        if '(' not in tree:
            self.set_error('进化树水平选择过高，导致没有树枝， 请选择较低的水平', code="12703503")
        
        
        def simple_name(name):
            name = name.group()
            return name.split(';')[-1].strip().replace(':', '-').strip('\'')  # replace用于去掉名称中带有冒号

        format_tree = re.sub(r'\'(.+?)\'', simple_name, tree)
        format_tree = re.sub(r'(\[)', '--temp_replace_left--', format_tree)  # 中括号在phylo中的读取会被特别识别，出现错误，后续对中括号进行暂时替换处理
        format_tree = re.sub(r'(\])', '--temp_replace_right--', format_tree)
        from Bio import Phylo

        open(output_file + '.temp', 'w').write(format_tree)
        newick_tree = Phylo.read(output_file + '.temp', 'newick')
        leaves = newick_tree.get_terminals()
        self.leaves = []
        for i in leaves:
            i.name = i.name.replace('--temp_replace_left--', '[')
            i.name = i.name.replace('--temp_replace_right--', ']')
            self.leaves.append(i.name)
        Phylo.write(newick_tree, output_file + '.temp2', 'newick')
        temp_tree = open(output_file + '.temp2').read()

        def replace_fun(matched):
            return ''
        temp_tree = re.sub(r'\'', replace_fun, temp_tree)
        temp_tree = re.sub(r'\"', replace_fun, temp_tree)
        with open(output_file, 'w') as w:
            w.write(temp_tree)


    def end(self):
        if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "进化分析结果目录", 0, "110204"],
            ["species_table.xls", "txt", "物种样本统计表", 0, "110207"],
            ["phylo_tree.tre", "tree", "进化树", 0, "110206"],
            ["species_group.xls", "txt", "物种在高层级的分类表", 0, "110205"],
            ["进化树图.pdf", "pdf", "进化树图", 0, ""],
            ["环形进化树图.pdf", "pdf", "环形进化树图", 0, ""]
        ])

        super(PlotTreeWorkflow, self).end()

        # self.set_end()
        # self.fire('end')
        # self._upload_result()
        # self._import_report_data()
        # self.step.finish()
        # self.step.update()
        # self.logger.info("运行结束!")
        # self._save_report_data()
        # # self._update("end")



if __name__ == '__main__':
    a = get_level_newicktree('57fd7c2b17b2bf377d2d6dae', level=8)
    print a

    def simple(name):
        name = name.group()
        return name.split(';')[-1].strip().replace(':', '-')  # replace用于去掉名称中带有冒号
    b = re.sub(r'\'(.+?)\'', simple, a)
    print b
    #wsheet = Sheet(data=data)
    #a =  PlotTreeWorkflow(wsheet)
    #a.run()
