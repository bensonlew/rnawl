# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# last_modify:20210817

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import unittest
import gevent.subprocess as subprocess
import json
import glob
import unittest
import shutil
import pandas as pd
from collections import OrderedDict
from biocluster.config import Config
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.medical_transcriptome.copy_file import CopyFile
import numpy as np


class DiffGenesetAllPiplineModule(Module):
    """
    该Module用于基因融合分析，默认使用方法
    """
    def __init__(self, work_id):
        super(DiffGenesetAllPiplineModule, self).__init__(work_id)
        options = [
            # {"name": "geneset_names", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "diff_path", "type": "string", "default": None},
            {"name": "annot_result", "type": "string", "default": None},
            {"name": "diff_method", "type": "string"},
            {"name": "transcript_exp_file", "type": "string"},
            {"name": "group", "type": "string"},
            {"name": "level", "type": "string", "default": "G"},
            {"name": "kegg_version", "type": "string", "default": None},
            {"name": "go_version", "type": "string", "default": '20200628'},
            # {"name": "reactome_version", "type": "string", "default": None},
            {"name": "species", "type": "string", "default": "Homo_sapiens"}
        ]
        self.add_option(options)
        self.inter_dirs = []
        self.genest_infos = OrderedDict()
        self.file_prepare = self.add_module("prok_rna.workflow_diffgt.diff_geneset_all")
        self.cluster = self.add_tool("prok_rna.exp_cluster")
        self.ipath = self.add_tool("prok_rna.geneset.ipath")
        self.kegg_level_path = ""
        self.geneset_list = ""
        self.all_list = ""
        self.select_geneset = ""
        self.geneset_analysis = self.add_module("prok_rna.workflow_diffgt.diff_geneset_analysis")
        self.file_dict = {}


    def check_options(self):
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run_file_prepare(self):
        opts = {
            "annot_result": self.option("annot_result"),
            "diff_path" : self.option("diff_path"),
            "diff_method": self.option("diff_method"),
            "transcript_exp_file": self.option("transcript_exp_file"),
            "kegg_version": self.option("kegg_version"),
            'go_version' :self.option("go_version"),
            "level" : self.option("level"),
            "species" :self.option("species")
        }
        self.file_prepare.set_options(opts)
        self.file_prepare.on("end",self.run_genesets_analysis)
        self.file_prepare.on("end", self.run_exp_cluster)
        self.file_prepare.on("end", self.run_ipath)
        self.file_prepare.run()

    def run_genesets_analysis(self):
        # self.check_genesets_infos()
        file_json_path = os.path.join(self.file_prepare.output_dir,"prepare_json")
        if not os.path.exists(file_json_path):
            raise Exception("未找到文件准备汇总json文件")
        with open(file_json_path,"r") as j:
            file_dict = json.load(j)
        self.kegg_level_path = file_dict["common_file"]["common_annot_file"]["kegg_level_table"]
        gensets_infos = file_dict["genesets"]
        geneset_info = sorted(gensets_infos.keys())[0]
        opts = self.opts_prepare(file_dict, gensets_infos[geneset_info])

        self.geneset_analysis.set_options(opts)
        # self.geneset_analysis.on('end', self.run_exp_cluster)
        self.geneset_analysis.run()


    def opts_prepare(self,file_dict,genesets):
        analysis_names = ["go","kegg"]
        opts={
                # "task_id":self.option("task_id"),
                "annot_result": self.option("annot_result"),
                "level" :self.option("level"),
                "species": self.option("species"),
                "go_list":file_dict["common_file"]["common_annot_file"]["go_list"],
                "kegg_table":file_dict["common_file"]["common_annot_file"]["kegg_table"],
                "kegg_table2": file_dict["common_file"]["common_annot_file"]["kegg_level_table"],
                "all_list": file_dict["common_file"]["common_annot_file"]["all_list"],
                "add_info" : file_dict["common_file"]["common_annot_file"]["add_info"],
                "kegg_version" : file_dict["common_file"]["common_annot_file"]["kegg_version"],
                "geneset_name": genesets["geneset_name"],
                "genes_num": genesets["gene_num"],
                "geneset_path":genesets["geneset_path"],
                "regulate": genesets["regulate"],
                "gene_list":genesets["file_path"]["go_enrich"]["gene_list_path"],
                "gene_multi" : genesets["file_path"]["kegg_class"]["multi_gene_list_path"],
                "go_class" :genesets["file_path"]["go_class"]["go_class_path"],
                "go_class_detail": genesets["file_path"]["go_class"]["go_class_detail_path"],
                "cog_class": genesets["file_path"]["cog_class"]["cog_class_path"]
            }
        return opts

    def run_exp_cluster(self):
        self.select_geneset, selecet_exp_file = self.get_cluster_data()
        options = dict(
            exp=selecet_exp_file,
            group=self.option('group'),
        )
        self.cluster.set_options(options)
        # self.cluster.on('end', self.set_output)
        self.cluster.run()

    def run_ipath(self):
        file_json_path = os.path.join(self.file_prepare.output_dir, "prepare_json")
        if not os.path.exists(file_json_path):
            raise Exception("未找到文件准备汇总json文件")
        with open(file_json_path, "r") as j:
            file_dict = json.load(j)
        gensets_infos = file_dict["genesets"]
        geneset_info = sorted(gensets_infos.keys())[0]
        genesets = gensets_infos[geneset_info]
        options = dict(
            geneset_kegg=genesets["file_path"]["kegg_class"]["multi_gene_list_path"],
            kegg_table=file_dict["common_file"]["common_annot_file"]["kegg_table"],
        )
        self.ipath.set_options(options)
        # self.cluster.on('end', self.set_output)
        self.ipath.run()

    def check_genesets_infos(self):
        # 准备聚类分析表达量文件
        file_json_path = os.path.join(self.file_prepare.output_dir, "prepare_json")
        if not os.path.exists(file_json_path):
            raise Exception("未找到文件准备汇总json文件")
        with open(file_json_path, "r") as j:
            file_dict = json.load(j)
        gensets_infos = file_dict["genesets"]
        geneset_num = dict()
        for geneset in gensets_infos:
            gene_num = gensets_infos[geneset]["gene_num"]
            geneset_num[geneset] = gene_num

        # select_geneset = min(geneset_num, key=lambda x:geneset_num[x])

        #临时测试修改
        # select_geneset = min(geneset_num, key=lambda x:geneset_num[x])
        #从这开始需要保存
        def min_geneset_select(geneset_num):
            max_geneset_num = max([geneset_num[x] for x in geneset_num])
            if max_geneset_num <= 50:
                return False
            else:
                select_geneset = min(geneset_num, key=lambda x: geneset_num[x] if geneset_num[x] >= 50 else 10000000)
                return select_geneset

        select_geneset = min_geneset_select(geneset_num)
        if not select_geneset:
            self.set_error("本次分析没有任何对比组超过50个差异基因,请核查分组或更换差异分析软件")
        # 到这结束

    def get_cluster_data(self):
        # 准备聚类分析表达量文件
        file_json_path = os.path.join(self.file_prepare.output_dir, "prepare_json")
        if not os.path.exists(file_json_path):
            raise Exception("未找到文件准备汇总json文件")
        with open(file_json_path, "r") as j:
            file_dict = json.load(j)
        gensets_infos = file_dict["genesets"]
        with open(self.option("group")) as r:
            r.readline()
            samples = []
            for line in r.readlines():
                samples.append(line.strip().split("\t")[0])

        geneset_path = gensets_infos[self.select_geneset]["geneset_path"]
        select_ids = pd.read_table(geneset_path)["seq_id"].tolist()
        exp_df = pd.read_table(self.option("transcript_exp_file"),index_col=0)
        exp_df = exp_df[samples]
        final_df = exp_df.loc[select_ids,]
        selecet_exp_file = os.path.join(self.work_dir,"cluset_exp")
        final_df.to_csv(selecet_exp_file,sep="\t")
        return self.select_geneset,selecet_exp_file


    def set_output(self):
        if not self.select_geneset:
            with open(os.path.join(self.output_dir,"results_info"),"w") as f:
                result_infos = {"has_result":False}
                json.dump(result_infos,f)
        else:
            analysis_json = os.path.join(self.geneset_analysis.output_dir,"analysis_json")
            with open(analysis_json,"r") as f:
                 analysis_dict = json.load(f)
            geneset_name = analysis_dict["geneset_name"]
            genset_result_dir = os.path.join(self.output_dir,geneset_name)
            if os.path.exists(genset_result_dir):
                shutil.rmtree(genset_result_dir)
            #add by fwy 20210125
            CopyFile().linkdir(self.geneset_analysis.output_dir, genset_result_dir)
            # os.system('ln -s {} {}'.format(geneset_analysis.output_dir, genset_result_dir))
            # shutil.copytree(geneset_analysis.output_dir,genset_result_dir)
            if os.path.exists(os.path.join(self.output_dir,"cluster")):
                cluster_dir = os.path.join(self.output_dir,"cluster")
                shutil.rmtree(cluster_dir)
            cluster_dir = os.path.join(self.output_dir,"cluster",self.select_geneset)
            shutil.copytree(self.cluster.output_dir, cluster_dir)
            if os.path.exists(os.path.join(self.output_dir,"ipath")):
                cluster_dir = os.path.join(self.output_dir,"ipath")
                shutil.rmtree(cluster_dir)
            cluster_dir = os.path.join(self.output_dir,"ipath",self.select_geneset)
            shutil.copytree(self.ipath.output_dir, cluster_dir)

        self.end()

    def get_genesets_info(self):
        diff_dir = self.option("diff_path")
        all_detail = pd.read_table(os.path.join(diff_dir, 'all_detail.txt'), header=0, sep='\t')
        all_detail['type'] = np.where(all_detail['seq_id'].str.contains("sRNA"), 'sRNA', 'mRNA')
        all_detail.fillna('', inplace=True)
        diff_dict_list = all_detail.to_dict('records')
        all_diff_list = list()
        for i in all_detail.groupby('compare'):
            compare = i[0]
            diff_pd = i[1]
            ctrl, test = compare.split('|')
            name = ctrl + '_vs_' + test + '_mRNA'
            sig_seqs = diff_pd[(~diff_pd['seq_id'].str.contains("sRNA")) & (diff_pd['significant'] == 'yes')][
                'seq_id'].tolist()
            sig_regulate = list(diff_pd['regulate'][diff_pd['significant'] == 'yes'])
            all_diff_list += sig_seqs
        all_diff_list = list(set(all_diff_list))
        geneset_name = 'All_Diff_mRNA'
        length = len(all_diff_list)
        self.genest_infos[geneset_name] = {
            "geneset_name": geneset_name,
            "compare": "All",
            "gene_num": length,
            "regulate": "all",
            "compare_path": os.path.join(diff_dir, 'all_detail.txt')
        }

    def get_selected_genset_info(self):
        self.get_genesets_info()
        geneset_name = self.genest_infos.keys()[0]
        if self.genest_infos[geneset_name]["gene_num"]<50:
            return False
        else:
            return geneset_name


    def run(self):
        super(DiffGenesetAllPiplineModule, self).run()
        self.logger.info("开始运行差异基因集数据挖掘数据准备")
        self.logger.info("首先解析基因集列表文件")
        self.select_geneset = self.get_selected_genset_info()
        if self.select_geneset:
            self.run_file_prepare()
            relys = [ self.file_prepare,self.geneset_analysis,self.cluster,self.ipath]
            self.on_rely(relys,self.set_output)
        else:
            self.set_output()

    def end(self):
        super(DiffGenesetAllPiplineModule, self).end()



class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna'
        data = {
            "id": "geneset_prepare" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "medical_transcriptome.workflow_diffgt.diff_geneset_all_pipline",
            "instant": False,
            "options": dict(
                diff_path= "/mnt/ilustre/users/sanger-dev/workspace/20200923/MedicalTranscriptome_medical_transcriptome_workflow_20200923_153802009/DiffexpBatch/output",
                annot_result ="/mnt/ilustre/users/sanger-dev/workspace/20200923/MedicalTranscriptome_medical_transcriptome_workflow_20200923_153802009/AnnotMerge__1/output",
                diff_method="DESeq2",
                kegg_version ="202007",
                reactome_version="202007",
                level = "G",
                gene_count_file ="/mnt/ilustre/users/sanger-dev/workspace/20200923/MedicalTranscriptome_medical_transcriptome_workflow_20200923_153802009/Quant/output/gene.count.matrix",
                gene_exp_file = "/mnt/ilustre/users/sanger-dev/workspace/20200923/MedicalTranscriptome_medical_transcriptome_workflow_20200923_153802009/Quant/output/ref.gene.tpm.matrix",
                group = "/mnt/ilustre/users/sanger-dev/workspace/20200907/Refrna_tsg_218672/remote_input/group_table/example_group_1528169151.txt",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()