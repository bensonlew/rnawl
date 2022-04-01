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
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class CeRnaWorkflow(Workflow):
    """
    ceRNA分析工作流
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        #self.rpc = False
        super(CeRnaWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'mirna_exp', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'mi_fa', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'cerna_exp', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'annotation', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'target', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {"name": "use_samples", "type": "string", "default": None},
            {"name": "group_dict", "type": "string", "default": None},
            {"name": "group_id", "type": "string", "default": None},
            {"name": "exp_level", "type": "string", "default": None},

            {'name': 'cecorr_pvalue_cutoff', 'type': 'float', 'default': 0.05},
            {'name': 'cecorr_pvalue_type', 'type': 'string', 'default': "pvalue"},
            {'name': 'cecorr_cutoff', 'type': 'float', 'default': 0.0},
            {'name': 'cecorr_method', 'type': 'string', 'default': "spearman"},
            {'name': 'cecorr_pvalue_method', 'type': 'string', 'default': "fdr_bh"},
            {'name': 'micorr_pvalue_cutoff', 'type': 'float', 'default': 1},
            {'name': 'micorr_pvalue_type', 'type': 'string', 'default': "pvalue"},
            {'name': 'micorr_cutoff', 'type': 'float', 'default': 0},
            {'name': 'micorr_method', 'type': 'string', 'default': "spearman"},
            {'name': 'micorr_pvalue_method', 'type': 'string', 'default': "fdr_bh"},
            {'name': 'mi_min_num', 'type': 'string', 'default': "3"},
            {'name': 'retain_type', 'type': 'string', 'default': "all"},

            {"name": "miset_id", "type": "string", "default":None},
            {"name": "level", "type": "string", "default":"T"},
            {"name": "type", "type": "string", "default": "T"},
            {"name": "group_dict", "type": "string", "default": None},
            {"name": "group_id", "type": "string", "default": None},
            {"name": "outtable", "type": "outfile", "format": "whole_transcriptome.common"},
            {"name": "task_id", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "trans2gene", "type": "string"},
            {"name": "taxon", "type": "string", "default": "Animal"}
        ]
        self.add_option(options)

        print "sheet opt is {}".format(self._sheet.options())
        self.set_options(self._sheet.options())
        # self.mrna_target_predict = self.add_module("whole_transcriptome.target_predict")
        # self.lncrna_target_predict = self.add_module("lnc_rna.target_predict")
        # self.extract_utr = self.add_tool("lnc_rna.cerna.extract_utr3")
        self.target2ce = self.add_tool("whole_transcriptome.cerna.mitarget2ce")
        self.mirna_family = self.add_tool("whole_transcriptome.cerna.smallrna_family_analyse")
        self.small_mrna_corr = self.add_tool("whole_transcriptome.cerna.expcorr_mitarget")
        # self.small_lncrna_corr = self.add_tool("whole_transcriptome.cerna.expcorr_mitarget")
        self.ce_corr = self.add_tool("whole_transcriptome.cerna.expcorr_cerna")
        self.tool_merge = self.add_tool("whole_transcriptome.cerna.merge_ce")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/07 Advanced_Analysis/01 ceRNA_analysis')
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
            "ref" : self.option("genome").prop['path'],
            "gtf" : self.option("mrna_gtf").prop['path'],
            "choose_list": self.work_dir + "/choose_m.list"
        }
        self.extract_utr.set_options(options)
        self.extract_utr.run()

    def run_small_family(self):
        '''
        小RNA家族
        '''
        options = {
            "matfa": self.option("mi_fa").prop['path'],
        }
        self.mirna_family.set_options(options)
        self.mirna_family.run()


    def run_target_ce(self):
        '''
        根据靶基因判断是否为ceRNA
        '''
        options = {
            "mrna_target": self.small_mrna_corr.output_dir + "/miRNA_mRNA_corr.xls",
            "mirna_family": self.mirna_family.output_dir  + "/known_miR_family.xls",
            "mi_min_num": self.option("mi_min_num"),
            "type": self.option("type"),
            "retain_type": self.option("retain_type")
        }
        self.target2ce.set_options(options)
        self.target2ce.run()

    def run_ce_corr(self):
        '''
        cerna 相关性分析
        '''
        options = {
            "exp_cerna" : self.option("cerna_exp").prop['path'],
            "ce_pair" : self.target2ce.output_dir + "/ce.xls",
            "corr_way" : self.option("cecorr_method"),
            "corr_cutoff" : self.option("cecorr_cutoff"),
            "output" : 'ceRNA_corr.xls',
            "type" : self.option("type")
        }
        if self.option("cecorr_pvalue_type") == "padjust":
            options.update({"qvalue_cutoff": self.option("cecorr_pvalue_cutoff")})
            options.update({"padjust_way": self.option("cecorr_pvalue_method")})
        else:
            options.update({"pvalue_cutoff": self.option("cecorr_pvalue_cutoff")})

        if self.option("mirna_exp").is_set:
            options.update({
                "exp_mirna" : self.option("mirna_exp").prop['path'],
                "corr_mirna2mrna" : self.small_mrna_corr.output_dir + "/miRNA_mRNA_corr.xls",

            })
        else:
            options.update({
                "corr_mirna2mrna" : self.mrna_target_predict.output_dir  + "/known_target.xls",

            })
        self.ce_corr.set_options(options)
        self.ce_corr.run()

    def run_mrna_corr(self):
        '''
        smallrna mrna 相关性分析
        '''
        options = dict(
            exp = self.option("mirna_exp").prop['path'],
            exp_target = self.option("cerna_exp").prop['path'],
            rna_target = self.option("target").prop['path'],
            corr_way = self.option("micorr_method"),
            corr_cutoff = float(self.option("micorr_cutoff")),
            type = self.option("type"),
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


    def run_merge(self):
        options = {
            "nodes" : "ceRNA_nodes.xls",
            "edges" : "ceRNA_edges.xls",
            "ce_corr" : self.ce_corr.output_dir + "/ceRNA_corr.xls"
        }
        if self.option("mirna_exp").is_set:
            options.update({
                "mi2mrna" : self.small_mrna_corr.output_dir + "/miRNA_mRNA_corr.xls",
            })
        else:
            options.update({
                "mi2mrna" :  self.mrna_target_predict.output_dir  + "/known_target.xls",
            })

        self.tool_merge.set_options(options)
        self.tool_merge.run()


    def run(self):
        # 设置运行逻辑
        # self.prepare_file()
        self.on_rely([self.mirna_family, self.small_mrna_corr], self.run_target_ce)
        self.target2ce.on("end", self.run_ce_corr)
        self.ce_corr.on("end", self.run_merge)
        self.tool_merge.on('end', self.set_output)
        self.get_run_log()
        self.run_mrna_corr()
        self.run_small_family()
        super(CeRnaWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("whole_transcriptome", table="cerna", main_id=self.option('main_id'),
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
        if self.option("mirna_exp").is_set:
            # files += [
            #     self.small_mrna_corr.output_dir + "/miRNA_mRNA_corr.xls",
            # ]
            for source in files:
                basename = os.path.basename(source)
                link_name = os.path.join(self.output_dir, basename)
                if os.path.exists(link_name):
                    os.remove(link_name)
                os.link(source, link_name)
                self.logger.info('succeed in linking {} to {}'.format(source, link_name))
            source_path = os.path.join(self.small_mrna_corr.output_dir, "miRNA_mRNA_corr.xls")
            link_path = os.path.join(self.output_dir, "miRNA_tar_RNA_corr.xls")
            if os.path.exists(link_path):
                os.remove(link_path)
            os.link(source_path, link_path)
            self.logger.info('succeed in linking {} to {}'.format(source_path, link_path))
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
            os.link(self.mrna_target_predict.output_dir  + "/known_target.xls", os.path.join(self.output_dir, "miRNA_mRNA_target.xls"))
            if os.path.exists(os.path.join(self.output_dir, "miRNA_lncRNA_target.xls")):
                os.remove(os.path.join(self.output_dir, "miRNA_lncRNA_target.xls"))
            os.link(self.lncrna_target_predict.output_dir  + "/known_target.xls", os.path.join(self.output_dir, "miRNA_lncRNA_target.xls"))
        self.set_db()



    def set_db(self):
        cerna_api = self.api.api("whole_transcriptome.cerna")
        if self.option("mirna_exp").is_set:
            mi2mrna = self.small_mrna_corr.output_dir + "/miRNA_mRNA_corr.xls"
        else:
            mi2mrna = os.path.join(self.output_dir, "miRNA_mRNA_target.xls")


        cerna_api.import_net_webroot(main_id=self.option("main_id"),
                                     result_dir=self.output_dir,
                                     mi2mrna=mi2mrna,
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
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        target_dir = self.get_workflow_output_dir()
        result_dir = self.add_upload_dir(self.output_dir)

        self.inter_dirs = [
            ["07 Advanced_Analysis", "", "高级分析结果目录",0],
            ["07 Advanced_Analysis/01 ceRNA_analysis", "", "ceRNA关联分析", 0]
        ]

        result_dir.add_relpath_rules([
            [".", "", "mirna关联分析结果目录", 0],
            ["./*_nodes.xls", "xls", "调控网络点信息表", 0],
            ["./*_edges.xls", "xls", "调控网络边信息列表", 0],
            ["./miRNA_mRNA_target.xls", "xls", "miRNA-mRNA靶基因预测表", 0],
            ["./miRNA_lncRNA_target.xls", "xls", "miRNA-lncRNA靶基因预测表", 0],
            ["./miRNA_mRNA_corr.xls", "xls", "miRNA-mRNA相关性分析表", 0],
            ["./miRNA_lncRNA_corr.xls", "xls", "miRNA-lncRNA相关性分析表", 0],
            ["./miRNA_tar_RNA_corr.xls", "xls", "miRNA-靶RNA对应关系表", 0],
            ["./ceRNA_corr.xls", "xls", "ceRNA-ceRNA对应关系表", 0],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        db = Config().get_mongo_client(mtype="whole_transcriptome")[Config().get_mongo_dbname("whole_transcriptome")]
        col1 = db["cerna"]
        col1.update({"_id": ObjectId(self.option("main_id"))}, {"$set": {"result_dir": target_dir}}, upsert=True)
        super(CeRnaWorkflow, self).end()
