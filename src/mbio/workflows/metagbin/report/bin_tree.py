# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# __modify__ = '2019.04.11'

from biocluster.workflow import Workflow
import os,re
import json
from biocluster.file import download

class BinTreeWorkflow(Workflow):
    """
    宏基因组进化树
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BinTreeWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "bin_list", "type": "string", "default": ""},  # bin序列序列
            {"name": "genome_list", "type": "string", "default": ""},  # 基因组文件
            {"name": "method", "type": "string"},  # 进化树分析类型
            {"name": "tree_type", "type": "string", "default": ""},  # 进化树的方法
            {"name": "ref", "type": "string", "default": ""},  # 上传的序列文件
            {"name": "num", "type": "int"},  # 样品总的数量
            {"name": "bins", "type": "string", "default": ""},  # bin的名称
            {"name": "genomes", "type": "string", "default": ""},  # 基因组名称
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tree_16s = self.add_tool("graph.phy_tree")
        self.tree_hgene = self.add_tool("graph.phy_tree")
        self.max_tree = self.add_tool("metagbin.plot_tree")
        self.ref_names = []

    def run(self):
        self.get_fasta()
        if self.option("num") <= 100:
            if self.option("method") in ['16s']:
                self.tree_16s.on("end",self.set_db)
                self.run_16s_tree()
            else:
                self.tree_hgene.on("end", self.set_db)
                self.run_coregene_tree()
        else:
            self.max_tree.on("end", self.set_db)
            self.run_max_tree()
        super(BinTreeWorkflow, self).run()

    def run_16s_tree(self):
        tree_type = ''
        if self.option("tree_type") in ["NJ"]:
            tree_type = "Neighbor_Joining(NJ)"
        elif self.option("tree_type") in ["ML"]:
            tree_type = "Maximum_Likehood(ML)"
        elif self.option("tree_type") in ["MP"]:
            tree_type = "Maximum_Parsimony(MP)"
        opts = {
            "fasta": self.work_dir + "/all.16s.fa",
            "method": tree_type,
        }
        self.tree_16s.set_options(opts)
        self.tree_16s.on("end", self.set_output)
        self.tree_16s.run()

    def run_coregene_tree(self):
        tree_type = ''
        if self.option("tree_type") in ["NJ"]:
            tree_type = "Neighbor_Joining(NJ)"
        elif self.option("tree_type") in ["ML"]:
            tree_type = "Maximum_Likehood(ML)"
        elif self.option("tree_type") in ["MP"]:
            tree_type = "Maximum_Parsimony(MP)"
        opts = {
            "fasta": self.work_dir + "/all.coregene.fa",
            "method": tree_type,
            "sequence_type": 'amino_acid',
        }
        self.tree_hgene.set_options(opts)
        self.tree_hgene.on("end", self.set_output)
        self.tree_hgene.run()

    def run_max_tree(self):
        seq = ''
        method = ''
        if self.option("method") in ['16s']:
            seq = self.work_dir + '/all.16s.fa'
            method = "nuc"
        elif self.option("method") in ['corgene']:
            seq = self.work_dir + '/all.coregene.fa'
            method = "pro"
        opts = {
            "seq_fa": seq,
            "method": method,
        }
        self.max_tree.set_options(opts)
        self.max_tree.on("end", self.set_output)
        self.max_tree.run()

    def get_fasta(self):
        if self.option("method") in ['16s']:
            with open(self.work_dir + "/all.16s.fa", 'w') as f:
                if self.option("bin_list") != "":
                    samples = eval(self.option("bin_list"))
                    for sample, seq in samples.items():
                        f.write(">{}\n{}\n".format(sample, seq))
                if self.option("genome_list") != "":
                    genoe = eval(self.option("genome_list"))
                    for genome_id, file in genoe.items():
                        file = file.replace('\\',"")
                        download(file, self.work_dir + '/seq/' + genome_id + ".16s.fa")
                        with open(self.work_dir + '/seq/' + genome_id + ".16s.fa", 'r') as g:
                            lines =g.readlines()
                            seq = lines[1].strip("\n\r")
                            f.write(">{}\n{}\n".format(genome_id, seq))
                if self.option("ref") != "":
                    path = self.option("ref").replace('\\', "")
                    download(path, self.work_dir + '/ref.fa')
                    with open(self.work_dir + '/ref.fa', 'r') as g:
                        lines = g.readlines()
                        for line in lines:
                            line = line.strip('\r\n')
                            f.write("{}\n".format(line))
                            if re.search(r'^>', line):
                                m = re.search('^>(.*)', line)
                                self.ref_names.append(m.group(1))
        elif self.option("method") in ['corgene']:
            with open(self.work_dir + "/all.coregene.fa", 'w') as f:
                samples = eval(self.option("bin_list"))
                for sample, seq in samples.items():
                    f.write(">{}\n{}\n".format(sample, seq))

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        self.set_output()
        if self.option("bins") != "":
            names = self.option("bins").split(",")
            for name in names:
                self.ref_names.append(name)
        if self.option("genomes") != "":
            names = self.option("genomes").split(",")
            for name in names:
                self.ref_names.append(name)
        str = ",".join(self.ref_names)
        main_id = self.option('main_id')
        self.api_path = self.api.api('metagbin.bin_tree')
        self.api_path.add_tree_detail(main_id, self.output_dir + '/phylo_tree.nwk', str)
        self.end()

    def set_output(self):
        if self.option("num") <= 100 and self.option("method") in ['16s']:
            if os.path.exists(self.output_dir + '/phylo_tree.nwk'):
                os.remove(self.output_dir + '/phylo_tree.nwk')
            os.link(self.tree_16s.work_dir + '/phylo_tree.nwk',self.output_dir + '/phylo_tree.nwk')
        if self.option("num") <= 100 and self.option("method") in ['corgene']:
            if os.path.exists(self.output_dir + '/phylo_tree.nwk'):
                os.remove(self.output_dir + '/phylo_tree.nwk')
            os.link(self.tree_hgene.work_dir + '/phylo_tree.nwk',self.output_dir + '/phylo_tree.nwk')
        if self.option("num") > 100:
            if os.path.exists(self.output_dir + '/phylo_tree.nwk'):
                os.remove(self.output_dir + '/phylo_tree.nwk')
            os.link(self.tree_hgene.output_dir + '/phylo_tree.nwk', self.output_dir + '/phylo_tree.nwk')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "进化树结果目录"],
            ["phylo_tree.nwk", "", "进化树结果nwk文件"],
        ])
        super(BinTreeWorkflow, self).end()