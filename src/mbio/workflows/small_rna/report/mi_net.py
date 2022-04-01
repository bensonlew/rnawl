# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
import re
import os
from bson.objectid import ObjectId
from biocluster.config import Config
import pandas as pd
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class MiNetWorkflow(Workflow):
    """
    表达量相关性
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MiNetWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'exp_matrix', 'type': 'infile', 'format': 'small_rna.common'},
            {'name': 'exp_target', 'type': 'infile', 'format': 'small_rna.common'},
            {'name': 'known_target', 'type': 'infile', 'format': 'small_rna.common'},
            {'name': 'novol_target', 'type': 'infile', 'format': 'small_rna.common'},
            {'name': 'all_target', 'type': 'infile', 'format': 'small_rna.common'},
            {"name": "ref", "type": "infile", "format": "small_rna.fasta"},
            {"name": "use_samples", "type": "string", "default": None},
            {"name": "known_pre", "type": "infile", "format": "small_rna.common"},
            {"name": "novol_pre", "type": "infile", "format": "small_rna.common"},
            {"name": "mirnaset_id", "type": "string", "default":None},
            {"name": "geneset_id", "type": "string", "default":None},
            {'name': 'species', 'type': 'string', 'default': 'vertebrates'},

            {'name': 'tf_sub_db', 'type': 'string', 'default': "All"},
            {'name': 'tss_up', 'type': 'int', 'default': 5000},
            {'name': 'tss_down', 'type': 'int', 'default': 1000},
            {"name": "thresh", "type": "float", "default": 0.0001},
            {"name": "pro_dis", "type": "int", "default": 500},
            {"name": "mood_pvalue", "type": "float", "default": 0.0001},
            {"name": "known_species", "type": "string", "default": None},
            {'name': 'gene_link', 'type': 'infile', 'format': 'small_rna.common'},
            {'name': 'gene_detail', 'type': 'infile', 'format': 'small_rna.common'},
            {"name": "db", "type": "string", "default": None},
            {"name": "genome", "type": "string", "default": None},
            {'name': 'output', 'type': 'string', 'default': None},
            {'name': 'pvalue_cutoff', 'type': 'float', 'default': 1},
            {'name': 'qvalue_cutoff', 'type': 'float', 'default': 1},
            {'name': 'corr_cutoff', 'type': 'float', 'default': 0.0},
            {'name': 'corr_way', 'type': 'string', 'default': "spearman"},
            {'name': 'padjust_way', 'type': 'string', 'default': "fdr_bh"},
            {'name': 'sig_type', 'type': 'int', 'default': 1},
            {'name': 'group_dict', 'type': 'string', 'default': None},
            {'name': 'group_id', 'type': 'string', 'default': None},
            {'name': 'main_id', 'type': 'string', 'default': None},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/05 Advanced_Analysis/01 (TF-)miRNA-mRNA_Analysis')
        self.inter_dirs = []
        self.tool_corr = self.add_tool("small_rna.expcorr_mitarget")
        self.tool_tfbinding = self.add_tool("small_rna.tf_binding")
        self.tool_merge = self.add_tool("small_rna.merge_tf_corr")
        self.smallrna_list = list()
        self.target_list = list()

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
        super(MiNetWorkflow, self).send_log(data)

    def run(self):
        self.get_run_log()
        if self.option("exp_target").is_set and self.option("db"):
            self.on_rely([self.tool_tfbinding, self.tool_corr], self.run_merge)
        elif self.option("exp_target").is_set:
            self.tool_corr.on("end", self.run_merge)
        elif self.option("db"):
            self.tool_tfbinding.on("end", self.run_merge)
        else:
            pass

        self.tool_merge.on('end', self.set_db)
        self.get_smallrna_list()
        self.get_target_list()
        self.filter_target()
        self.filter_pre()
        self.get_non_cor_target()
        if self.option("exp_target").is_set:
            self.run_corr()
        else:
            pass
        if self.option("db"):
            self.run_tf()
        else:
            pass

        if (not self.option("exp_target").is_set) and (not self.option("db")):
            self.run_merge()

        super(MiNetWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("small_rna", table="sg_minet", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def get_smallrna_list(self):
        if self.option("mirnaset_id") in [None, ""]:
            self.smallrna_list = None
        else:
            with open(self.option("mirnaset_id"), "r") as mirnaset_f:
                self.smallrna_list = [line.strip() for line in mirnaset_f.readlines()]
    def get_target_list(self):
        if self.option("geneset_id") in [None, ""]:
            self.target_list = None
        else:
            with open(self.option("geneset_id"), "r") as geneset_f:
                self.target_list = [line.strip() for line in geneset_f.readlines()]

    def get_non_cor_target(self):
        with open(self.work_dir + "/target_filter.xls", "r") as target_filter_f, open(self.work_dir + "/target_filter_noncor.xls", "w") as target_filter_w:
            header = target_filter_f.readline()
            headers = header.split("\t")
            if headers[-1].startswith("mirtarbase"):
                target_filter_w.write("\t".join(headers[:5] + [headers[-1]]))
            else:
                target_filter_w.write("\t".join(headers[:5] + [headers[-2]]) + "\n")
            for line in target_filter_f:
                cols = line.split("\t")
                if headers[-1].startswith("mirtarbase"):
                    target_filter_w.write("\t".join(cols[:5] + [cols[-1]]))
                else:
                    target_filter_w.write("\t".join(cols[:5] + [cols[-2]]) + "\n")


    def filter_pre(self):
        with open(self.option("known_pre").prop["path"], "r") as known_f,  open(self.work_dir + "/known_prefilter.xls", "w") as known_w:
            header = known_f.readline()
            known_w.write(header)
            for line in known_f.readlines():
                cols = line.strip().split("\t")
                if self.smallrna_list:
                    if cols[0] in self.smallrna_list:
                        known_w.write(line)
                else:
                    known_w.write(line)
        with open(self.option("novol_pre").prop["path"], "r") as known_f,  open(self.work_dir + "/novol_prefilter.xls", "w") as known_w:
            header = known_f.readline()
            known_w.write(header)
            for line in known_f.readlines():
                cols = line.strip().split("\t")
                if self.smallrna_list:
                    if cols[0] in self.smallrna_list:
                        known_w.write(line)
                else:
                    known_w.write(line)


    def filter_target(self):
        if self.option("exp_target").is_set:
            samples = self.option("use_samples").split("|")
            self.logger.info("exp_target is {}".format(self.option("exp_target").prop['path']))
            self.logger.info("target list is {}".format(self.target_list))
            exp = pd.read_table(self.option("exp_target").prop['path'], header=0, index_col=0)
            exp_list = list(exp.index)
            if self.target_list:
                gene_filter = list(set(self.target_list).intersection(set(exp_list)))
            else:
                gene_filter = exp_list
            exp_filter = exp.loc[gene_filter, samples]
            exp_filter.to_csv(self.work_dir + "/exp_matrix_filter.xls", sep="\t")
        else:
            pass

        if self.option("all_target").is_set:
            with open(self.option("all_target").prop['path'], "r") as target_f, open(self.work_dir + "/target_filter.xls", "w") as target_filter_w:
                header = target_f.readline()
                target_filter_w.write(header)
                pair = set()
                for line in target_f.readlines():
                    cols = line.strip().split("\t")
                    # 历史项目上传靶基因结果中可能没有gene_id, 用转录本代替
                    if cols[2] == "_":
                        cols[2] = cols[1]

                    if self.smallrna_list and cols[0] not in self.smallrna_list:
                        pass
                    elif self.target_list and cols[2] not in self.target_list:
                        pass
                    else:
                        if cols[0] + "__" + cols[2] in pair:
                            pass
                        else:
                            target_filter_w.write(line)
                            pair.add(cols[0] + "__" + cols[2])
        else:
            with open(self.option("known_target").prop['path'], "r") as known_target_f, open(self.option("novol_target").prop['path'], "r") as novol_target_f, open(self.work_dir + "/target_filter.xls", "w") as target_filter_w:
                    header = known_target_f.readline()
                    target_filter_w.write(header)
                    pair = set()
                    for line in known_target_f.readlines():
                        cols = line.strip().split("\t")
                        if self.smallrna_list and cols[0] not in self.smallrna_list:
                            pass
                        elif self.target_list and cols[2] not in self.target_list:
                            pass
                        else:
                            if cols[0] + "__" + cols[2] in pair:
                                pass
                            else:
                                target_filter_w.write(line)
                                pair.add(cols[0] + "__" + cols[2])
                    _ = novol_target_f.readline()
                    for line in novol_target_f.readlines():
                        cols = line.strip().split("\t")
                        if self.smallrna_list and cols[0] not in self.smallrna_list:
                            pass
                        elif self.target_list and cols[2] not in self.target_list:
                            pass
                        else:
                            if cols[0] + "__" + cols[2] in pair:
                                pass
                            else:
                                target_filter_w.write(line)
                                pair.add(cols[0] + "__" + cols[2])

    def run_corr(self):
        options = dict(
            exp=self.option("exp_matrix").prop['path'],
            exp_target = self.work_dir + "/exp_matrix_filter.xls",
            rna_target = self.work_dir + "/target_filter.xls",
            # corr_way = self.option("corr_way"),
            corr_cutoff = self.option("corr_cutoff"),
            pvalue_cutoff = self.option("pvalue_cutoff"),
            qvalue_cutoff = self.option("qvalue_cutoff"),
            output='corr.xls',
        )
        '''
            known_target = self.option("known_target").prop['path'],
            novol_target = self.option("novol_target").prop['path'],
        '''

        self.tool_corr.set_options(options)
        self.tool_corr.run()

    def run_tf(self):
        genome = Config().SOFTWARE_DIR + "/database/Genome_DB_finish/" + self.option("genome")
        options = dict(
            ref = genome,
            known_pre = self.work_dir + "/known_prefilter.xls",
            # self.option("known_pre").prop["path"],
            novol_pre = self.work_dir + "/novol_prefilter.xls",
            # self.option("novol_pre").prop["path"],
            species = self.option("species"),
            gene_link = self.option("gene_link"),
            gene_detail = self.option("gene_detail"),
            db = self.option("db"),
            mood_pvalue = self.option("mood_pvalue"),
            thresh = self.option("mood_pvalue"),
        )
        if self.option("tf_sub_db"):
            options.update({
                "sub_species": self.option("tf_sub_db")
            })
        self.tool_tfbinding.set_options(options)
        self.tool_tfbinding.run()

    def run_merge(self):
        options = dict(
            nodes="nodes.xls",
            edges="edges.xls",
        )

        if self.option("exp_target").is_set:
            options.update({
                "corr" : self.tool_corr.output_dir + "/corr.xls"
            })
        else:
            options.update({
                "corr" : self.work_dir + "/target_filter_noncor.xls"
            })
        if self.option("db"):
            options.update({
                "tfbinding" : self.tool_tfbinding.output_dir + "/binding.xls"
            })

        self.tool_merge.set_options(options)
        self.tool_merge.run()

    def target_rename(self, target):
        with open(target, "rb") as target_f, open(target + ".rename", "w") as target_fo:
            head_line = target_f.readline()
            if len(head_line.strip("\n").split("\t")) == 9:
                target_fo.write("miRNA name\ttarget gene id\ttarget gene name\ttarget gene description\tevidence_1(prediction_result)\tevidence_2 (miRTarBase_result)\tcorr\tp value\tq value\n")
            else:
                target_fo.write("miRNA name\ttarget gene id\ttarget gene name\ttarget gene description\tevidence_1(prediction_result)\tevidence_2 (miRTarBase_result)\t\n")

            for line in target_f:
                cols = line.strip("\n").split("\t")

                if len(cols) == 9:
                    target_fo.write("\t".join([cols[0],
                                               cols[2],
                                               cols[3],
                                               cols[4],
                                               "yes",
                                               cols[5],
                                               cols[6],
                                               cols[7],
                                               cols[8]
                    ]) + "\n")
                else:
                    target_fo.write("\t".join([cols[0],
                                               cols[2],
                                               cols[3],
                                               cols[4],
                                               "yes",
                                               cols[5]
                    ]) + "\n")
            return target + ".rename"



    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        minet_api = self.api.api("small_rna.minet")
        tf = self.tool_tfbinding.output_dir + "/binding.xls"
        if self.option("exp_target").is_set:
            target = self.tool_corr.output_dir + "/corr.xls"
        else:
            target = self.work_dir + "/target_filter_noncor.xls"
        nodes = self.tool_merge.output_dir + "/nodes.xls"
        edges = self.tool_merge.output_dir + "/edges.xls"
        minet_api.import_net_webroot(self.option("main_id"), self.output_dir, tf=tf, target=target, nodes=nodes, edges=edges)

        target = self.target_rename(target)


        base = os.path.basename(tf)
        if os.path.exists(self.output_dir + "/TF_miRNA_detail.xls"):
            os.remove(self.output_dir + "/TF_miRNA_detail.xls")
        if os.path.exists(tf):
            os.link(tf, self.output_dir + "/TF_miRNA_detail.xls")

        if os.path.exists(self.output_dir + "/miRNA_mRNA_detail.xls"):
            os.remove(self.output_dir + "/miRNA_mRNA_detail.xls")
        if os.path.exists(target):
            os.link(target, self.output_dir + "/miRNA_mRNA_detail.xls")

        if os.path.exists(self.output_dir + "/All_nodes.xls"):
            os.remove(self.output_dir + "/All_nodes.xls")
        if os.path.exists(nodes):
            os.link(nodes, self.output_dir + "/All_nodes.xls")

        if os.path.exists(self.output_dir + "/All_edges.xls"):
            os.remove(self.output_dir + "/All_edges.xls")
        if os.path.exists(edges):
            os.link(edges, self.output_dir + "/All_edges.xls")

        '''
        allfiles = [tf, target, nodes, edges]
        for i in allfiles:
            base = os.path.basename(i)
            if os.path.exists(self.output_dir + "/" + base):
                os.remove(self.output_dir + "/" + base)
            if os.path.exists(i):
                os.link(i, self.output_dir + "/" + base)
        '''
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
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["05 Advanced_Analysis", "", "高级分析结果目录",0],
            ["05 Advanced_Analysis/01 (TF-)miRNA-mRNA_Analysis", "", "(TF-)miRNA-mRNA关联分析", 0]
        ]

        result_dir.add_relpath_rules([
            [".", "", "(TF-)miRNA-mRNA关联分析文件", 0],
            ["./All_nodes.xls", "xls", "调控网络点信息表", 0],
            ["./All_edges.xls", "xls", "调控网络边信息列表", 0],
            ["./miRNA_mRNA_detail.xls", "xls", "miRNA-mRNA对应关系表", 0],
            ["./TF_miRNA_detail.xls", "xls", "miRNA-靶基因对应关系表 ", 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        db = Config().get_mongo_client(mtype="small_rna")[Config().get_mongo_dbname("small_rna")]
        col1 = db["sh_minet"]
        col1.update({"_id": ObjectId(self.option("main_id"))}, {"$set": {"result_dir": target_dir}}, upsert=True)
        super(MiNetWorkflow, self).end()
