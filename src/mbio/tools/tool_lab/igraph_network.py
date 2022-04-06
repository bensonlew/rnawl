# -*- coding: utf-8 -*-
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import unittest
import re
import math
import numpy as np
from mbio.packages.whole_transcriptome.utils import runcmd
__author__ = 'gdq'

# normal_dict =["#E64B35B2" ,"#4DBBD5B2", "#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2", "#B09C85B2"]
normal_dict = ["#374E55FF", "#DF8F44FF" ,"#00A1D5FF" ,"#B24745FF" ,"#79AF97FF", "#6A6599FF", "#80796BFF"]
large_dict = ["#7FC97F","#BEAED4", "#FDC086" ,"#FFFF99" ,"#386CB0" ,"#F0027F", "#BF5B17", "#666666", "#1B9E77" ,"#D95F02" ,"#7570B3" ,
         "#E7298A", "#66A61E" ,"#E6AB02", "#A6761D", "#666666" ,"#A6CEE3" ,"#1F78B4", "#B2DF8A" ,"#33A02C", "#FB9A99" ,"#E31A1C",
         "#FDBF6F", "#FF7F00" ,"#CAB2D6" ,"#6A3D9A", "#FFFF99" ,"#B15928" ,"#FBB4AE" ,"#B3CDE3", "#CCEBC5" ,"#DECBE4" ,"#FED9A6",
         "#FFFFCC" ,"#E5D8BD", "#FDDAEC" ,"#F2F2F2" ,"#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4" ,"#E6F5C9" ,"#FFF2AE", "#F1E2CC",
         "#CCCCCC", "#E41A1C" ,"#377EB8" ,"#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33" ,"#A65628" ,"#F781BF" ,"#999999", "#66C2A5",
         "#FC8D62", "#8DA0CB" ,"#E78AC3" ,"#A6D854" ,"#FFD92F" ,"#E5C494", "#B3B3B3" ,"#8DD3C7" ,"#FFFFB3" ,"#BEBADA", "#FB8072",
         "#80B1D3" ,"#FDB462" ,"#B3DE69" ,"#FCCDE5" ,"#D9D9D9" ,"#BC80BD" ,"#CCEBC5", "#FFED6F"]

class IgraphNetworkAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(IgraphNetworkAgent, self).__init__(parent)
        options = [
            # 输入表格 分三列,第一列为id,第二列为log2fc,第三列为pajust
            dict(name="project_type", type='string', default="custom"),
            dict(name="nodes_file", type="infile", format="ref_rna_v2.common"),
            dict(name="edges_file", type="infile", format="ref_rna_v2.common"),
            dict(name="layout", type="string", default="layout_in_circle"),
            dict(name="line_style", type="string", default="solid"),
            dict(name="node_style", type='string', default="circle"),
            dict(name="string", type='string', default="circle"),
            dict(name="line_width_method", type='string', default="no_change"),#["no_change","add","minus","multiply","divide","log"]
            dict(name="line_width_scale", type='float', default= 0),
            dict(name="node_size_method", type='string', default="no_change"),#["no_change","add","minus","multiply","divide","log"]
            dict(name="node_size_scale", type='float', default=0),
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(IgraphNetworkAgent, self).end()


class IgraphNetworkTool(Tool):
    """
    Differential analysis based on edgeR, DEGseq, DEseq2.
    """
    def __init__(self, config):
        super(IgraphNetworkTool, self).__init__(config)
        # self.rscript = '/bioinfo/tool_lab/MEGENA/miniconda3/bin/Rscript'
        self.rscript = '/bioinfo/rna/miniconda2/bin/Rscript'
        self.r_igraph = self.config.PACKAGE_DIR + "/tool_lab/igraph_network/igraph_network.py"
        self.color_dict = ""
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'miniconda2/bin/python'
        # self.gcc = software_dir + '/gcc/5.1.0/bin'
        # self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        # self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        # if "sanger-dev" in self.config.SOFTWARE_DIR:
        #     self.r_path = software_dir + "/program/R-3.3.3/bin:$PATH"
        # else:
        #     self.r_path = software_dir + "/program/R-3.3.1_gcc5.1/bin:$PATH"
        # self.r_path = software_dir + "/bioinfo/tool_lab/MEGENA/miniconda3/bin:$PATH"
        self.r_path = software_dir + "/bioinfo/rna/miniconda2/bin:$PATH"
        self.set_environ(PATH=self.r_path)
        # # self.r_path = software_dir + "/program/R-3.3.1_gcc5.1/bin:$PATH"
        # self._r_home = software_dir + "/program/R-3.3.1_gcc5.1/lib64/R/"
        # self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1_gcc5.1/lib64/R/lib:$LD_LIBRARY_PATH"
        # self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
        self.node_dict = {'blank':0,'solid':1,'dashed':2,'dotted':3,'dotdash':4,'longdash':5,'twodash':6}
        self.nodes_file = ""
        self.edges_file = ""


    def igraph(self):
        cmd = '{} {} '.format(self.python_path, self.r_igraph)
        cmd += '-{} {} '.format("nodes_file", self.nodes_file)
        cmd += '-{} {} '.format("edges_file",  self.edges_file)
        cmd += '-{} {} '.format("layout", self.option("layout"))
        cmd += '-{} {} '.format("node_style", self.option("node_style"))
        cmd += '-{} {} '.format("line_style",  self.node_dict[self.option("line_style")])
        print cmd
        cmd_name = 'igraph'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} 运行成功".format(cmd_name))
            else:
                self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704306")
        else:
            self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704307")

    def prepare_plot_file(self):
        if self.option("project_type") == "custom" :
            node_df = pd.read_table(self.option("nodes_file").prop["path"])
            edge_df = pd.read_table(self.option("edges_file").prop["path"])
            node_df.columns = ["id", "name", "group", "size"]
            edge_df.columns = ["from", "to", "type", "weight"]
            node_df_f = self.scale_info_deal(node_df,"size",self.option("node_size_method"),self.option("node_size_scale"))
            edge_df_df = self.scale_info_deal(edge_df,"size",self.option("line_width_method"),self.option("line_width_scale"))
            #准备点的配色
            node_group_num = node_df_f.drop_duplicates("group").shape[0]
            if node_group_num <=7:
                self.color_dict = normal_dict
            elif node_group_num <= 73:
                self.color_dict = large_dict
            else:
                self.set_error("分组过多,无法绘图")
            target_color = self.color_dict[:node_group_num]
            group_list = node_df_f.drop_duplicates("group")["group"]
            color_map ={k:v for k,v in zip(group_list,target_color)}
            node_df_f["plot_color"] = node_df_f["group"].map(color_map)
            node_df_f = node_df_f.sort_values("group")
            t_node_df_f = node_df_f.set_index("id").to_dict("index")
            edge_df_df["plot_color"] = edge_df_df["from"].apply(lambda x: t_node_df_f[x]["plot_color"])
            node_df_f.to_csv(os.path.join(self.work_dir,"nodes.txt"), sep="\t", index=False)
            edge_df_df.to_csv(os.path.join(self.work_dir, "edges.txt"), sep="\t", index=False)
            self.modify(os.path.join(self.work_dir,"nodes.txt"),os.path.join(self.work_dir,"nodes_final.txt"))
            self.modify(os.path.join(self.work_dir, "edges.txt"), os.path.join(self.work_dir, "edges_final.txt"))
            self.nodes_file = os.path.join(self.work_dir,"nodes_final.txt")
            self.edges_file = os.path.join(self.work_dir, "edges_final.txt")

    def modify(self,raw_file,target_file):
            with open(raw_file, "r") as r, open(target_file, "w") as w:
                header = r.readline()
                w.write(header)
                for line in r.readlines():
                    line = line.strip().split("\t")
                    line[-1] = "\"" + line[-1] + "\""
                    line = "\t".join(line)
                    w.write(line + "\n")



    def scale_info_deal(self,df,column,method,scale):
        if method == "no_change":
            df =df
        elif method == "add":
            df[column] = df[column].apply(lambda x: x + scale)
        elif method == "minus":
            df[column] = df[column].apply(lambda x: x - scale)
        elif method == "multiply":
            df[column] = df[column].apply(lambda x: x * scale)
        elif method == "divide":
            df[column] = df[column].apply(lambda x: float(x/scale))
        elif method == "log":
            df[column] = df[column].apply(lambda x :math.log(x,scale))
        else:
            self.set_error(u"不存在{}方法".format(scale))
        return  df



    def set_output(self):
        for i in ["igraph.pdf","igraph.png",'igraph.svg']:
            if os.path.exists(os.path.join(self.output_dir,i)):
                os.remove(os.path.join(self.output_dir,i))
            os.link(os.path.join(self.work_dir,i),os.path.join(self.output_dir,i))


    def run(self):
        super(IgraphNetworkTool, self).run()
        self.prepare_plot_file()
        self.igraph()
        self.set_output()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    test_dir = '/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/tools/denovo_rna_v2/test_files'

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "IgraphNetwork" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.diff_ma_new",
            "instant": False,
            "options": dict(
                raw_file='s3://refrnav2/files/m_188/188_5ffbaeead3e00/mbs6_v5o6eq0967dj319dmpsmrq/workflow_results/07DiffExpress_G/HFL_vs_HGL.edger.annot.xls',
                # exp='/mnt/ilustre/users/sanger-dev/workspace/20200113/Refrna_tsg_36819/Quant/output/gene.tpm.matrix',
                # method="edgeR",
                project_type = "ref_rna_v2",
                pvalue=0.05,
                fc= 2,
                x_axis_name="log10(TPM)",
                y_axis_name= "log2(FC)",
                title_name="MA Plot",
                color = "ref_blue_grey"
            )
        }

        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()


