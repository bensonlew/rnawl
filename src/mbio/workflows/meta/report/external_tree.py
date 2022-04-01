# -*- coding: utf-8 -*-


import os
import re
from biocluster.workflow import Workflow
import copy
from bson import ObjectId
#from mbio.packages.beta_diversity.filter_newick import get_level_newicktree
# import datetime
# import json
from collections import defaultdict
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class ExternalTreeWorkflow(Workflow):
    """
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExternalTreeWorkflow, self).__init__(wsheet_object)
        options = [
            #{"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "level", "type": 'int', "default": 9},
            #{"name": "topN", "type": 'int', "default": 0},
            {"name": "otu_id", "type": 'string', "default": ''},
            {"name": "main_id", "type": 'string', "default": ''},
            {"name": "params", "type": 'string', "default": ''},
            {"name": "group_id", "type": 'string'},
            {"name": "update_info", "type": 'string'},
            {"name": "group_detail", "type": 'string', "default": ""},
            #{"name": "color_level_id", "type": 'int', "default": 0},
            #{"name": "sample_group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "method" , "type" : "string", "default": "ML"},
            {"name": "upload_spe" , "type" : "infile", "format": "meta.fasta"},
            {"name": "select_spe_fasta" , "type" : "infile", "format": "meta.fasta"}

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.phly_tree_tool = self.add_tool('graph.phy_tree')


    def check(self):
        # if self.option("level") >= self.option("color_level_id"):
        #     self.set_error("颜色设置id必须大于选择的分类水平", code="12703501")
        if not self.option("otu_id"):
            self.set_error("必须提供OTU的主表id", code="12704901")

    def run_tree(self):
        method_map = {
            'ML' : 'Maximum_Likehood(ML)',
            'MP' : 'Maximum_Parsimony(MP)',
            'NJ' : 'Neighbor_Joining(NJ)'
        }

        self.name_map = self.change_name()
        if self.option('method') == 'MP':
            if len(self.name_map) > 500 : #序列条数限制
                self.set_error('The number of sequences cannot be greater than 500 when use MP', code="12704904")
        elif self.option('method') in ['ML','NJ']:
            if len(self.name_map) > 500 :
                self.set_error('The number of sequences cannot be greater than 500 when use ML or NJ', code="12704905")
        # os.system('cat {} {} > {}'.format(self.option('select_spe_fasta').path, self.new_upload_fasta, self.fasta))
        # self.logger.info('select_spe_fasta:' + self.option('select_spe_fasta'))
        # self.logger.info('upload_spe: ' + self.option('upload_spe'))

        method = method_map[self.option('method')]
        opts =  {
            'method' : method,
            'fasta' : self.fasta,
            'align_method':'mafft',
            'sequence_type':'no_coding',
        }
        if self.option('method') == 'ML':
            opts['tree_software'] = 'iqtree'

        if self.option('method') == 'MP':
            opts['bootstrap'] = 100
        else:
            opts['bootstrap'] = 500

        self.phly_tree_tool.set_options(opts)
        self.phly_tree_tool.on('end',self.set_db)
        self.phly_tree_tool.run()

    def change_name(self):
        self.fasta = self.work_dir + '/all.fasta'
        name_map = {}
        seq_id = 1
        ref_names = []

        space = re.compile('\s+')

        fwr = open(self.fasta, 'w')
        with open(self.option('select_spe_fasta').path) as f:
            for line in f:
                if line[0] == '>':
                    line = line.strip()
                    name = space.split(line)[0][1:]
                    ref_names.append(name)
                    new_name = 'seq'+str(seq_id)+"n"
                    fwr.write('>'+new_name+"\n")
                    name_map[new_name] = name
                    seq_id += 1
                else:
                    fwr.write(line)

        self.name_change_log = self.work_dir + '/upload_spe_name_change.log'
        with open(self.option('upload_spe').path) as f2 , open(self.name_change_log, 'w') as fw:
            fw.write("#ori_name\tnew_name\n")
            for line in f2:
                if line[0] == '>':
                    line = line.strip()
                    spline = space.split(line,1)
                    name = spline[0][1:]
                    tmp_id = 0
                    ori_name = name
                    while name in ref_names:
                        tmp_id +=1
                        name = ori_name + '_' + str(tmp_id)
                    if name != ori_name:
                        fw.write(ori_name+"\t"+name+"\n")

                    new_name = 'seq'+str(seq_id)+"n"
                    fwr.write(">"+new_name+'\n')
                    name_map[new_name] = name
                    seq_id += 1

                else:
                    fwr.write(line)
        fwr.close()
        return name_map




    def change_tree_name(self,intree,new_file):
        def change(x, mp):
            c_mp = copy.deepcopy(mp)
            for k in mp:
                if k not in x:
                    c_mp.pop(k)
                    continue
                spx = x.split(k)
                for i in range(len(spx)):
                    c_str = spx[i]
                    spx[i] = change(c_str,c_mp)
                new_name = c_mp[k]
                new_x = new_name.join(spx)
                return new_x

            return x

        with open(intree) as fr:
            line = fr.read()
        line = change(line,self.name_map)

        with open(new_file,'w') as fw:
            fw.write(line)

    def run(self):
        self.run_tree()
        super(ExternalTreeWorkflow, self).run()


    def set_db(self):
        api_tree = self.api.phylo_tree
        self.tree_file = self.phly_tree_tool.output_dir + '/phylo_tree.nwk'
        self.new_tree = self.output_dir + '/phylo_tree.nwk'
        self.change_tree_name(self.tree_file,self.new_tree)
        self.logger.info('开始导入文件：%s' %self.new_tree)
        api_tree.add_phylo_tree_info_2(self.option('main_id'),self.new_tree, seq_num=len(self.name_map.keys()))
        os.rename(self.output_dir + '/phylo_tree.nwk', self.output_dir + '/phylo_tree.tre')
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
                "submit_loc": "phylo_tree_extend",
                "interaction": 1,
                "main_table": "sg_phylo_tree",
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "进化分析结果目录", 0, "110204"],
            ["phylo_tree.tre", "tree", "进化树", 0, ""],
            ["环形进化树图.pdf", "pdf", "环形进化树图", 0, ""],
            ["进化树图.pdf", "pdf", "进化树图", 0, ""],
        ])
        super(ExternalTreeWorkflow, self).end()

