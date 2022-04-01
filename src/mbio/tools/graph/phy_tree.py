# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
# import copy
import collections
import pandas as pd
import subprocess
import math
import shutil


class PhyTreeAgent(Agent):
    """
    version 1.0
    author zouxuan
    last_modified:20180123
    """

    def __init__(self, parent):
        super(PhyTreeAgent, self).__init__(parent)
        options = [
            {"name": "abundance_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "otu_group", "type": "infile", "format": "toolapps.group_table"},
            {"name": "sample_group", "type": "infile", "format": "toolapps.group_table"},
            {"name": "sequence_type", "type": "string", "default": "coding"},
            {"name": "method", "type": "string", "default": "Neighbor_Joining(NJ)"},
            {"name": "align", "type": "string", "default": "true"},
            {"name": "newicktree", "type": "outfile", "format": "meta.beta_diversity.newick_tree"},
            {"name": "abundance_table_out", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "otu_group_out", "type": "infile", "format": "toolapps.group_table"},
            {"name": "align_method", "type": "string", "default": "muscle"} , #muscle ,mafft
            {"name": "tree_software", "type": "string", "default": "megacc"} , #megacc, iqtree
            {"name": "bootstrap" , "type": "int","default":1000}
        ]
        self.add_option(options)

    def check_options(self):
        """
        重写参数检查
        """
        # if not self.option('abundance_table').is_set:
        #     raise OptionError('必须提供输入文件:物种丰富文件')
        # else:
        #     otulist = [line.split('\t')[0] for line in open(self.option('abundance_table').prop['path'])][
        #               1:]  # 获取所有OTU/物种名
        #     sample = open(self.option('abundance_table').prop['path']).readline().rstrip().split('\t')[1:]
        if not self.option('fasta').is_set:
            raise OptionError('必须提供输入序列文件', code="32300601")
        self.option('fasta').get_info()
        if int(self.option('fasta').prop['seq_number']) < 2:
            raise OptionError('序列条数必须大于等于2', code="32300602")
        #if int(self.option('fasta').prop['seq_number']) > 300:
        #    raise OptionError('%s序列条数必须小于等于300', variables=(self.option('fasta').prop['seq_number'],), code="32300603")
        if self.option('method') not in ["Neighbor_Joining(NJ)", "Maximum_Likehood(ML)", "Maximum_Parsimony(MP)"]:
            raise OptionError('计算方法只能选择Neighbor_Joining(NJ)、Maximum_Likehood(ML)、Maximum_Parsimony(MP)', code="32300604")
        if self.option('align') not in ["true", "false"]:
            raise OptionError('是否多重比对只能传入true或者false', code="32300605")
        if self.option('sequence_type') not in ["coding", "no_coding", "amino_acid"]:
            raise OptionError('序列类型只能传入coding、no_coding、amino_acid', code="32300606")
        seq_name = self.option('fasta').get_all_seq_name()
        if self.option('otu_group').is_set:
            self.option('otu_group').get_info()
            otu_seq = self.option('otu_group').prop['sample_name']
            for seq in otu_seq:
                if seq not in seq_name:
                    raise OptionError('序列分组文件中的序列:%s不存在于序列文件中', variables=(seq), code="32300607")
        # if self.option('abundance_table').is_set:
        #     abu_sample = self.option('abundance_table').get_sample_info
        #     abu_seq = []
        #     with open(self.option('abundance_table').prop['path'], 'r') as f:
        #         heads = f.readline().rstrip().split('\t')
        #         while 1:
        #             line = f.readline().rstrip()
        #             abu_seq.append(line.split('\t')[0])
        #     for seq in abu_seq:
        #         if seq not in seq_name:
        #             raise OptionError('丰度文件中的序列:{}不存在于序列文件中'.format(seq))
        #     if self.option('sample_group').is_set:
        #         self.option('sample_group').get_info()
        #         sample_group = self.option('sample_group').prop['sample_name']
        #         for sample in sample_group:
        #             if sample not in abu_sample:
        #                 raise OptionError('样品分组文件中的样品:{}不存在于丰度文件中'.format(sample))

    def set_resource(self):
        """
        设置所需资源
        """

        if self.option("align") in ['true']:
            self._cpu = 10
            total = os.path.getsize(self.option("fasta").prop["path"])
            total = int(math.ceil(total / (1024 * 1024 * 1024)))
            total = 30 + int(total * 45)
            #total = '40G'
            self._memory = "{}G".format(total)
        else:
            self._cpu = 2
            self._memory = '20G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "系统进化分析结果目录"],
            ["phylo_tree.nwk", "nwk", "进化树文件"],
            ["abundance_table.xls", "xls", "序列丰度文件"]
        ])
        super(PhyTreeAgent, self).end()


class PhyTreeTool(Tool):
    def __init__(self, config):
        super(PhyTreeTool, self).__init__(config)
        self.muscle = '/bioinfo/phylogenetic/muscle-3.8.31-release/muscle'
        self.mega_cc = '/bioinfo/phylogenetic/mega_10.1.7/megacc'
        self.mega_path = self.config.SOFTWARE_DIR + '/bioinfo/phylogenetic/mega_10.1.7'
        self.mafft_path = self.config.SOFTWARE_DIR+'/bioinfo/align/mafft-7.299-with-extensions/bin/'
        self.colors_num = 21  # 目前支持的颜色数量
        self.max_species_group_name = 0
        self.iqtree = '/bioinfo/metaGenomic/iqtree-1.6.8-Linux/bin/iqtree'
        if self.option('method') == "Neighbor_Joining(NJ)":
            method = "NJ"
        if self.option('method') == "Maximum_Likehood(ML)":
            method = "ML"
        if self.option('method') == "Maximum_Parsimony(MP)":
            method = "MP"
        if self.option('sequence_type') == 'coding':
            self.mao_file = self.mega_path + '/infer_' + method + '_coding.mao'
        if self.option('sequence_type') == 'no_coding':
            self.mao_file = self.mega_path + '/infer_' + method + '_nucleotide.mao'
        if self.option('sequence_type') == 'amino_acid':
            self.mao_file = self.mega_path + '/infer_' + method + '_protein.mao'

    def run(self):
        """
        运行
        """
        super(PhyTreeTool, self).run()
        self.run_plot_tree()
        if self.option('abundance_table').is_set:
            self.get_abun()
        self.end()

    def get_abun(self):
        table = pd.DataFrame(pd.read_table(self.option("abundance_table").prop["path"], sep='\t', index_col=0))
        table_name = table.index.name
        # table = table.ix[list((table > 0).any(axis=1))]  # 去除都为0的物种/功能/基因
        if self.option("sample_group").is_set:
            group = pd.DataFrame(pd.read_table(self.option("sample_group").prop["path"], sep='\t', index_col=0))
            group["sample"] = group.index
            group_sample = group.join(table.T, on="sample").groupby(group.columns[0]).mean()  # 求均值
            abund = group_sample.T
            abund.index.name = table_name
        else:
            abund = table
        abund['Col_sum'] = abund.apply(lambda x: x.sum(), axis=1)
        abund_table = abund.sort_values(by=['Col_sum'], ascending=0)
        del abund_table["Col_sum"]
        abund_table = abund_table.ix[list((abund_table > 0).any(axis=1))]  # 去除都为0的物种/功能/基因
        if len(abund_table) < 1:
            self.set_error('在所选参数下数据为空，请重新设置水平或分组方案参数!', code="32300601")
            self.set_error('在所选参数下数据为空，请重新设置水平或分组方案参数!', code="32300602")
        self.abund_table_path = os.path.join(self.output_dir, "abundance_table.xls")
        abund_table.to_csv(self.abund_table_path, sep="\t", encoding="utf-8")

    def run_plot_tree(self):
        fasta_fiel = self.option('fasta').prop['path']
        cmd= "sed -i 's/\s.*//g' "+fasta_fiel
        try :
            subprocess.check_output(cmd,shell=True)
            self.logger.info("fasta文件生成成功")
        except subprocess.CalledProcessError:
            self.set_error("fasta文件生成失败", code="32300603")
            self.set_error("fasta文件生成失败", code="32300604")
        if self.option('align') == "true":
            align = self.work_dir + '/align.fasta'
            if self.option("align_method") in ["mafft"]:  #add 20191015
                cmd = "{}mafft --thread 10 {} > {}".format(self.mafft_path, fasta_fiel, align)
                self.logger.info("开始运行{}软件，进行比对".format(self.option("align_method")))
                command = subprocess.Popen(cmd, shell=True)
                command.communicate()
                if command.returncode == 0:
                    self.logger.info("完成比对！")
                else:
                    self.set_error("mafft运行出错！", code="32300610")

            else:
                cmd1 = '%s -in %s -out %s' % (self.muscle, fasta_fiel, align)
                command1 = self.add_command('cmd_1', cmd1)
                command1.run()
                self.wait(command1)
                if command1.return_code == 0:
                    self.logger.info("align succeed")
                else:
                    self.set_error("align failed", code="32300605")
                    self.set_error("align failed", code="32300606")
            fasta_fiel = align

        else:
            pass
        out = self.work_dir + '/phylo_tree'
        if self.option("tree_software") == 'megacc':
            new_mao = self.work_dir+'/config.mao'
            shutil.copy(self.mao_file, new_mao)
            ret = os.system("sed -i 's/No. of Bootstrap Replications.*/No. of Bootstrap Replications        = %s/' %s "%(self.option('bootstrap'),new_mao))
            if ret == 0:
                self.logger.info('修改bootstrap成功')
            else:
                self.logger.warning('修改bootstrap失败，用原配置参数运行')
            cmd2 = '%s -a %s -d %s -o %s' % (self.mega_cc, new_mao, fasta_fiel, out)
        else:
            cmd2 = '%s -s %s -m MFP -pre phylo_tree  -bb %s' %(self.iqtree, fasta_fiel, 1000)
        self.logger.info(cmd2)
        command2 = self.add_command('cmd_2', cmd2)
        command2.run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("phylo_tree succeed")
        else:
            self.set_error("phylo_tree failed", code="32300607")
        tree_file = self.output_dir + "/phylo_tree.nwk"
        if os.path.exists(tree_file):
            os.remove(tree_file)

        if self.option("tree_software") == 'megacc':
            if not os.path.exists(self.work_dir + "/phylo_tree.nwk"):
                self.set_error("phylo_tree failed,该序列无法构建进化树", code="32300609")
            os.link(self.work_dir + "/phylo_tree.nwk",tree_file)
        else:
            if not os.path.exists(self.work_dir + "/phylo_tree.treefile"):
                self.set_error("phylo_tree failed,该序列无法构建进化树", code="32300609")
            os.link(self.work_dir + "/phylo_tree.treefile",tree_file)
