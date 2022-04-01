# -*- coding: utf-8 -*-


import re
from biocluster.workflow import Workflow
import copy
from mbio.packages.metaasv.common_function import link_dir


class ExternalTreeWorkflow(Workflow):
    """
    metaasv 个性化系统发生进化树
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExternalTreeWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "level", "type": 'int', "default": 9},
            {"name": "asv_id", "type": 'string', "default": ''},
            {"name": "main_id", "type": 'string', "default": ''},
            {"name": "group_id", "type": 'string'},
            {"name": "update_info", "type": 'string'},
            {"name": "group_detail", "type": 'string', "default": ""},
            {"name": "method" , "type" : "string", "default": "ML"},#进化树构建方法
            {"name": "upload_spe" , "type" : "infile", "format": "meta.fasta"},##外源参考物种
            {"name": "select_spe_fasta" , "type" : "infile", "format": "meta.fasta"}##选择的物种

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.phly_tree_tool = self.add_tool('graph.phy_tree')

    def check(self):
        if not self.option("asv_id"):
            self.set_error("必须提供OTU的主表id")

    def run_tree(self):
        """
        运行 phy_tree进行进化树分析
        :return:
        """
        method_map = {
            'ML' : 'Maximum_Likehood(ML)',
            'MP' : 'Maximum_Parsimony(MP)',
            'NJ' : 'Neighbor_Joining(NJ)'
        }

        self.name_map = self.change_name()
        if self.option('method') == 'MP':
            if len(self.name_map) > 500 : #序列条数限制
                self.set_error('The number of sequences cannot be greater than 500 when use MP')
        elif self.option('method') in ['ML','NJ']:
            if len(self.name_map) > 500 :
                self.set_error('The number of sequences cannot be greater than 500 when use ML or NJ')
        self.method = method_map[self.option('method')]
        opts =  {
            'method' : self.method,
            'fasta' : self.fasta,
            'align_method':'mafft',
            'sequence_type':'no_coding',
        }
        if self.option('method') == 'ML':
            opts['tree_software'] = 'iqtree'
            self.software = "IQ-TREE"
        else:
            self.software = "Mega"

        if self.option('method') == 'MP':
            opts['bootstrap'] = 100
        else:
            opts['bootstrap'] = 500
        self.phly_tree_tool.set_options(opts)
        self.phly_tree_tool.on('end',self.set_db)
        self.phly_tree_tool.run()

    def change_name(self):
        """
        根据序列文件将物种名称进行改换，后面生成结果文件后再重新替换回来
        :return:
        """
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
        """
        替换名称
        :param intree:
        :param new_file:
        :return:
        """
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
        """
        运行
        :return:
        """
        self.run_tree()
        super(ExternalTreeWorkflow, self).run()

    def set_db(self):
        """
        改换名称和导入MongoDB
        :return:
        """
        api_tree = self.api.api("metaasv.phylo_tree")
        self.tree_file = self.phly_tree_tool.output_dir + '/phylo_tree.nwk'
        self.new_tree = self.output_dir + '/phylo_tree.nwk'
        self.change_tree_name(self.tree_file,self.new_tree)
        self.logger.info('开始导入文件：%s' %self.new_tree)
        api_tree.add_phylo_tree_info_2(self.option('main_id'),self.new_tree, seq_num=len(self.name_map.keys()), software=self.software)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.phly_tree_tool.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "进化分析结果目录", 0, ""],
            ["phylo_tree.tre", "tree", "进化树", 0, ""],
        ])

        super(ExternalTreeWorkflow, self).end()

