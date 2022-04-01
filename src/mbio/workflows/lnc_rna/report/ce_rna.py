# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modifiy = modified 2019.03.18

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from biocluster.wsheet import Sheet
import os
import re
import shutil
import json
import glob
import pandas as pd
from biocluster.config import Config
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from biocluster.core.function import filter_error_info, link, CJsonEncoder


class CeRnaWorkflow(Workflow):
    """
    ceRNA分析工作流
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        # self.rpc = False
        super(CeRnaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "smallrna", "type": "infile", "format": "lnc_rna.fasta"},  # 输入文件
            {"name": "mrna", "type": "infile", "format": "lnc_rna.fasta"},  # 输入文件
            {"name": "utr3", "type": "outfile", "format": "lnc_rna.fasta"},  # 输入文件
            {"name": "lncrna", "type": "infile", "format": "lnc_rna.fasta"},  # 输入文件
            {"name": "known_lnc_fa", "type": "infile", "format": "lnc_rna.fasta"},  # 输入文件
            {"name": "novol_lnc_fa", "type": "infile", "format": "lnc_rna.fasta"},  # 输入文件
            {"name": "genome", "type": "infile", "format": "lnc_rna.fasta"},  # 输入文件
            {"name": "mrna_gtf", "type": "infile", "format": "lnc_rna.gtf"},  # 输入文件
            {"name": "known_lnc_gtf", "type": "infile", "format": "lnc_rna.gtf"},  # 输入文件
            {"name": "novol_lnc_gtf", "type": "infile", "format": "lnc_rna.gtf"},  # 输入文件
            {"name": "lncrna_gtf", "type": "infile", "format": "lnc_rna.gtf"},  # 输入文件
            {'name': 'smallrna_exp', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'cerna_exp', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'annotation', 'type': 'infile', 'format': 'lnc_rna.common'},
            {"name": "use_samples", "type": "string", "default": None},
            {"name": "group_dict", "type": "string", "default": None},
            {"name": "group_id", "type": "string", "default": None},
            {"name": "exp_level", "type": "string", "default": None},

            {'name': 'corr_pvalue_cutoff', 'type': 'float', 'default': 0.05},
            {'name': 'corr_pvalue_type', 'type': 'string', 'default': "pvalue"},

            {'name': 'corr_cutoff', 'type': 'float', 'default': 0.0},
            {'name': 'corr_method', 'type': 'string', 'default': "spearman"},
            {'name': 'corr_padjust_method', 'type': 'string', 'default': "fdr_bh"},
            {'name': 'micorr_pvalue_cutoff', 'type': 'float', 'default': 1},
            {'name': 'micorr_pvalue_type', 'type': 'string', 'default': "pvalue"},
            {'name': 'micorr_cutoff', 'type': 'float', 'default': 0},
            {'name': 'micorr_method', 'type': 'string', 'default': "spearman"},
            {'name': 'micorr_padjust_method', 'type': 'string', 'default': "fdr_bh"},

            {"name": "geneset_id", "type": "string", "default": None},
            {"name": "lncset_id", "type": "string", "default": None},
            {"name": "type", "type": "string", "default": "G"},
            {"name": "outtable", "type": "outfile", "format": "lnc_rna.common"},
            {"name": "task_id", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "trans2gene", "type": "string"},
            {"name": "taxon", "type": "string", "default": "Animal"},

            {"name": "miranda_score", "type": "string", "default": "160"},
            {"name": "miranda_energy", "type": "string", "default": "-20"},
            {"name": "miranda_strict", "type": "string", "default": "on"},
            {"name": "ps_robot_score", "type": "string", "default": "2.5"},

        ]
        self.add_option(options)

        self.set_options(self._sheet.options())
        self.mrna_target_predict = self.add_module("lnc_rna.target_predict")
        self.lncrna_target_predict = self.add_module("lnc_rna.target_predict")
        self.extract_utr = self.add_tool("lnc_rna.cerna.extract_utr3")
        self.target2ce = self.add_tool("lnc_rna.cerna.mitarget2ce")
        self.mirna_family = self.add_tool("lnc_rna.cerna.smallrna_family_analyse")
        self.small_mrna_corr = self.add_tool("lnc_rna.cerna.expcorr_mitarget")
        self.small_lncrna_corr = self.add_tool("lnc_rna.cerna.expcorr_mitarget")
        self.ce_corr = self.add_tool("lnc_rna.cerna.expcorr_cerna")
        self.tool_merge = self.add_tool("lnc_rna.cerna.merge_ce")

        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/07 Advanced_Analysis/01 ceRNA_Analysis')
        self.inter_dirs = []

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(CeRnaWorkflow, self).send_log(data)

    def get_utr(self):
        '''
        获取utr序列
        '''
        options = {
            "ref": self.option("genome").prop['path'],
            "gtf": self.option("mrna_gtf").prop['path'],
            "choose_list": self.work_dir + "/choose_m.list"
        }
        self.extract_utr.set_options(options)
        self.extract_utr.run()

    def run_small_family(self):
        '''
        小RNA家族
        '''
        options = {
            "matfa": self.option("smallrna").prop['path'],
        }
        self.mirna_family.set_options(options)
        self.mirna_family.run()

    def get_lncrna_list(self):
        with open(self.option("geneset_id"), "r") as geneset_f:
            self.mrna_list = [line.strip() for line in geneset_f.readlines()]

    def get_mrna_list(self):
        with open(self.option("lncset_id"), "r") as geneset_f:
            self.lncrna_list = [line.strip() for line in geneset_f.readlines()]

    def filter_exp(self):
        if self.option("cerna_exp").is_set:
            samples = self.option("use_samples").split("|")
            exp = pd.read_table(self.option("cerna_exp").prop['path'], header=0, index_col=0)
            exp_list = list(exp.index)
            gene_filter = list(set(self.mrna_list + self.lncrna_list).intersection(set(exp_list)))
            exp_filter = exp.loc[gene_filter, samples]
            exp_filter.to_csv(self.work_dir + "/exp_matrix_filter.xls", sep="\t")
        else:
            pass

        if self.option("smallrna_exp").is_set:
            samples = self.option("use_samples").split("|")
            exp = pd.read_table(self.option("smallrna_exp").prop['path'], header=0, index_col=0)
            exp_filter = exp.loc[:, samples]
            exp_filter.to_csv(self.work_dir + "/exp_s_matrix_filter.xls", sep="\t")
        else:
            pass

    def convert_genelist2tranlist(self, gene_list, lncrna_list):
        '''
        根据基因列表提取转录本列表
        '''
        lnc_tran2gene = self.option("known_lnc_gtf").get_txpt_gene_dic()
        lnc_tran2gene.update(self.option("novol_lnc_gtf").get_txpt_gene_dic())
        lnc_trans_list = []
        for t, g in lnc_tran2gene.items():
            if g in lncrna_list:
                lnc_trans_list.append(t)
        m_tran2gene = self.option("mrna_gtf").get_txpt_gene_dic()
        m_trans_list = []
        for t, g in m_tran2gene.items():
            if g in gene_list:
                m_trans_list.append(t)
        return m_trans_list, lnc_trans_list

    def prepare_file(self):
        self.get_mrna_list()
        self.get_lncrna_list()
        self.filter_exp()
        if self.option("type") == "G":
            self.mrna_list, self.lncrna_list = self.convert_genelist2tranlist(self.mrna_list, self.lncrna_list)
        with open(self.work_dir + "/choose_m.list", 'w') as f:
            f.write("\n".join(self.mrna_list))
        if self.option("lncrna").is_set:
            self.option("lncrna").choose_seq_by_list(self.lncrna_list, self.work_dir + "/choose_lnc.fa")
        else:
            if self.option("known_lnc_fa").is_set:
                self.option("known_lnc_fa").choose_seq_by_list(self.lncrna_list, self.work_dir + "/choose_lnc.fa")
            if self.option("novol_lnc_fa").is_set:
                self.option("novol_lnc_fa").choose_seq_by_list(self.lncrna_list, self.work_dir + "/choose_lnc.fa",
                                                               mode='a')

        self.option("mrna").choose_seq_by_list(self.mrna_list, self.work_dir + "/choose_m.fa")

    def run_mrna_target(self):
        '''
        mirna 靶基因(mrna)分析
        '''
        if self.option("taxon").lower() == "animal":
            target = self.extract_utr.output_dir + "/utr3.fa"
        else:
            target = self.work_dir + "/choose_m.fa"
        method = {
            "Animal": "miranda",
            "Plant": "psrobot"
        }
        options = {
            "known": self.option("smallrna").prop['path'],
            "ref": target,
            "method": method[self.option("taxon")],
            "anno_detail": self.option("annotation").prop['path'],
            "type": self.option("taxon").lower(),
            "miranda_score": self.option("miranda_score"),
            "miranda_energy": self.option("miranda_energy"),
            "miranda_strict": self.option("miranda_strict"),
            "ps_robot_score": self.option("ps_robot_score"),
            "min_support": 1
        }
        self.mrna_target_predict.set_options(options)
        self.mrna_target_predict.run()

    def run_lncrna_target(self):
        '''
        mirna 靶基因(lnc)分析
        '''
        method = {
            "Animal": "miranda",
            "Plant": "psrobot"
        }
        options = {
            "known": self.option("smallrna").prop['path'],
            "ref": self.work_dir + "/choose_lnc.fa",
            "method": method[self.option("taxon")],
            "anno_detail": self.option("annotation").prop['path'],
            "type": self.option("taxon").lower(),
            "miranda_score": self.option("miranda_score"),
            "miranda_energy": self.option("miranda_energy"),
            "miranda_strict": self.option("miranda_strict"),
            "ps_robot_score": self.option("ps_robot_score"),
            "min_support": 1
        }
        self.lncrna_target_predict.set_options(options)
        self.lncrna_target_predict.run()

    def run_target_ce(self):
        '''
        根据靶基因判断是否为ceRNA
        '''
        options = {
            "mrna_target": self.mrna_target_predict.output_dir + "/known_target.xls",
            "lncrna_target": self.lncrna_target_predict.output_dir + "/known_target.xls",
            "mirna_family": self.mirna_family.output_dir + "/known_miR_family.xls",
            "type": self.option("type")
        }
        self.target2ce.set_options(options)
        self.target2ce.run()

    def run_ce_corr(self):
        '''
        cerna 相关性分析
        '''
        options = {
            "exp_cerna": self.option("cerna_exp").prop['path'],
            "ce_pair": self.target2ce.output_dir + "/ce.xls",
            "corr_way": self.option("corr_method"),
            "corr_cutoff": self.option("corr_cutoff"),
            "output": 'ceRNA_corr.xls',
            "type": self.option("type")
        }
        if self.option("corr_pvalue_type") == "padjust":
            options.update({"qvalue_cutoff": self.option("corr_pvalue_cutoff")})
            options.update({"padjust_way": self.option("corr_padjust_method")})
        else:
            options.update({"pvalue_cutoff": self.option("corr_pvalue_cutoff")})

        if self.option("smallrna_exp").is_set:
            options.update({
                "exp_mirna": self.option("smallrna_exp").prop['path'],
                "corr_mirna2mrna": self.small_mrna_corr.output_dir + "/miRNA_mRNA_corr.xls",
                "corr_mirna2lncrna": self.small_lncrna_corr.output_dir + "/miRNA_lncRNA_corr.xls",

            })
        else:
            options.update({
                "corr_mirna2mrna": self.mrna_target_predict.output_dir + "/known_target.xls",
                "corr_mirna2lncrna": self.lncrna_target_predict.output_dir + "/known_target.xls",

            })
        self.ce_corr.set_options(options)
        self.ce_corr.run()

    def run_mrna_corr(self):
        '''
        smallrna mrna 相关性分析
        '''
        options = dict(
            exp=self.work_dir + "/exp_s_matrix_filter.xls",
            exp_target=self.work_dir + "/exp_matrix_filter.xls",
            rna_target=self.mrna_target_predict.output_dir + "/known_target.xls",
            corr_way=self.option("micorr_method"),
            corr_cutoff=float(self.option("micorr_cutoff")),
            type=self.option("type"),
            output='miRNA_mRNA_corr.xls',
        )
        if self.option("micorr_pvalue_type") == "pvalue":
            options.update({
                'pvalue_cutoff': self.option("micorr_pvalue_cutoff")
            })
        else:
            options.update({
                'qvalue_cutoff': self.option("micorr_pvalue_cutoff")
            })
        self.small_mrna_corr.set_options(options)
        self.small_mrna_corr.run()

    def run_lncrna_corr(self):
        '''
        smallrna lncrna 相关性分析
        '''
        options = dict(
            exp=self.work_dir + "/exp_s_matrix_filter.xls",
            exp_target=self.work_dir + "/exp_matrix_filter.xls",
            rna_target=self.lncrna_target_predict.output_dir + "/known_target.xls",
            corr_way=self.option("micorr_method"),
            corr_cutoff=float(self.option("micorr_cutoff")),
            type=self.option("type"),
            output='miRNA_lncRNA_corr.xls',
        )

        if self.option("micorr_pvalue_type") == "pvalue":
            options.update({
                'pvalue_cutoff': self.option("micorr_pvalue_cutoff")
            })
        else:
            options.update({
                'qvalue_cutoff': self.option("micorr_pvalue_cutoff")
            })
        self.small_lncrna_corr.set_options(options)
        self.small_lncrna_corr.run()

    def run_merge(self):
        options = {
            "nodes": "ceRNA_nodes.xls",
            "edges": "ceRNA_edges.xls",
            "ce_corr": self.ce_corr.output_dir + "/ceRNA_corr.xls"
        }
        if self.option("smallrna_exp").is_set:
            options.update({
                "mi2mrna": self.small_mrna_corr.output_dir + "/miRNA_mRNA_corr.xls",
                "mi2lncrna": self.small_lncrna_corr.output_dir + "/miRNA_lncRNA_corr.xls"
            })
        else:
            options.update({
                "mi2mrna": self.mrna_target_predict.output_dir + "/known_target.xls",
                "mi2lncrna": self.lncrna_target_predict.output_dir + "/known_target.xls"
            })

        self.tool_merge.set_options(options)
        self.tool_merge.run()

    def run(self):
        # 设置运行逻辑
        self.get_run_log()
        self.prepare_file()

        self.on_rely([self.mirna_family, self.mrna_target_predict, self.lncrna_target_predict], self.run_target_ce)
        if self.option("taxon") == "Animal":
            self.extract_utr.on('end', self.run_mrna_target)
        if self.option("smallrna_exp").is_set:
            self.on_rely([self.small_mrna_corr, self.small_lncrna_corr, self.target2ce], self.run_ce_corr)
            self.mrna_target_predict.on('end', self.run_mrna_corr)
            self.lncrna_target_predict.on('end', self.run_lncrna_corr)
            self.on_rely([self.ce_corr, self.small_mrna_corr, self.small_lncrna_corr], self.run_merge)
        else:
            self.target2ce.on('end', self.run_ce_corr)
            self.on_rely([self.ce_corr, self.mrna_target_predict, self.lncrna_target_predict], self.run_merge)

        self.tool_merge.on('end', self.set_output)

        if self.option("taxon") == "Animal":
            self.get_utr()
        else:
            self.run_mrna_target()
        self.run_lncrna_target()
        self.run_small_family()
        super(CeRnaWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("lnc_rna", table="sg_cerna", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_output(self):
        '''
        self.logger.info("结果路径为{}".format(self.target_predict.output_dir))
        #output_dir = self.annotation.output_dir
        self.logger.info("结果路径为{}".format(self.target_predict.output_dir))
        '''
        files = [self.tool_merge.output_dir + "/ceRNA_nodes.xls",
                 self.tool_merge.output_dir + "/ceRNA_edges.xls",
                 self.ce_corr.output_dir + "/ceRNA_corr.xls"]
        if self.option("smallrna_exp").is_set:
            files += [
                self.small_mrna_corr.output_dir + "/miRNA_mRNA_corr.xls",
                self.small_lncrna_corr.output_dir + "/miRNA_lncRNA_corr.xls"
            ]
            for source in files:
                basename = os.path.basename(source)
                link_name = os.path.join(self.output_dir, basename)
                if os.path.exists(link_name):
                    os.remove(link_name)
                os.link(source, link_name)
                self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        else:
            for source in files:
                basename = os.path.basename(source)
                link_name = os.path.join(self.output_dir, basename)
                if os.path.exists(link_name):
                    os.remove(link_name)
                os.link(source, link_name)
                self.logger.info('succeed in linking {} to {}'.format(source, link_name))
            if os.path.exists(os.path.join(self.output_dir, "miRNA_mRNA_target.xls")):
                os.remove(os.path.join(self.output_dir, "miRNA_mRNA_target.xls"))
            os.link(self.mrna_target_predict.output_dir + "/known_target.xls",
                    os.path.join(self.output_dir, "miRNA_mRNA_target.xls"))
            if os.path.exists(os.path.join(self.output_dir, "miRNA_lncRNA_target.xls")):
                os.remove(os.path.join(self.output_dir, "miRNA_lncRNA_target.xls"))
            os.link(self.lncrna_target_predict.output_dir + "/known_target.xls",
                    os.path.join(self.output_dir, "miRNA_lncRNA_target.xls"))
        self.set_db()

    def set_db(self):
        cerna_api = self.api.api("lnc_rna.cerna")
        if self.option("smallrna_exp").is_set:
            mi2mrna = self.small_mrna_corr.output_dir + "/miRNA_mRNA_corr.xls"
            mi2lncrna = self.small_lncrna_corr.output_dir + "/miRNA_lncRNA_corr.xls"
        else:
            mi2mrna = os.path.join(self.output_dir, "miRNA_mRNA_target.xls")
            mi2lncrna = os.path.join(self.output_dir, "miRNA_lncRNA_target.xls")

        cerna_api.import_net_webroot(main_id=self.option("main_id"),
                                     result_dir=self.output_dir,
                                     mi2mrna=mi2mrna,
                                     mi2lncrna=mi2lncrna,
                                     ce_corr=self.ce_corr.output_dir + "/ceRNA_corr.xls",
                                     nodes=self.tool_merge.output_dir + "/ceRNA_nodes.xls",
                                     edges=self.tool_merge.output_dir + "/ceRNA_edges.xls")
        self.end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if re.match(r'tsanger:', workflow_output):
            workflow_output = workflow_output.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        elif re.match(r'sanger:', workflow_output):
            workflow_output = workflow_output.replace('sanger:', '/mnt/ilustre/data/')
        elif re.match(r'^\w+://\S+/.+$', workflow_output):
            pass
        else:
            self.set_error("json output wrong")
        self.workflow_output = workflow_output
        return workflow_output

    def end(self):
        target_dir = self.get_workflow_output_dir()
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.rename(self.output_dir + "/ceRNA_corr.xls", self.output_dir + "/ceRNA_ceRNA_corr.xls")
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["07 Advanced_Analysis", "", "高级分析结果目录", 0],
            ["07 Advanced_Analysis/01 ceRNA_Analysis", "", "ceRNA关联分析", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "ceRNA关联分析文件", 0],
            ["./*_nodes.xls", "xls", "调控网络点信息表", 0],
            ["./*_edges.xls", "xls", "调控网络边信息列表", 0],
            ["./miRNA_mRNA_target.xls", "xls", "miRNA-mRNA靶基因预测表", 0],
            ["./miRNA_lncRNA_target.xls", "xls", "miRNA-lncRNA靶基因预测表", 0],
            ["./miRNA_mRNA_corr.xls", "xls", "miRNA-mRNA相关性分析表", 0],
            ["./miRNA_lncRNA_corr.xls", "xls", "miRNA-lncRNA相关性分析表", 0],
            ["./ceRNA_ceRNA_corr.xls", "xls", "ceRNA-ceRNA对应关系表", 0],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        db = Config().get_mongo_client(mtype="lnc_rna")[Config().get_mongo_dbname("lnc_rna")]
        col1 = db["sg_cerna"]
        col1.update({"_id": ObjectId(self.option("main_id"))}, {"$set": {"result_dir": target_dir}}, upsert=True)
        super(CeRnaWorkflow, self).end()
