# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
import re
import os
from bson.objectid import ObjectId
from biocluster.config import Config
import pandas as pd
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class MiNetWorkflow(Workflow):
    """
    表达量相关性
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MiNetWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'exp_matrix', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'exp_target', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'target', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {"name": "ref", "type": "infile", "format": "whole_transcriptome.fasta"},
            {"name": "task_id", "type": "string", "default": None},
            {"name": "use_samples", "type": "string", "default": None},
            {"name": "known_pre", "type": "infile", "format": "whole_transcriptome.common"},
            {"name": "novol_pre", "type": "infile", "format": "whole_transcriptome.common"},
            {"name": "mirnaset_id", "type": "string", "default":None},
            {"name": "geneset_id", "type": "string", "default":None},
            {'name': 'species', 'type': 'string', 'default': 'vertebrates'},
            {'name': 'tf_sub_db', 'type': 'string', 'default': "all"},
            {'name': 'tss_up', 'type': 'int', 'default': 5000},
            {'name': 'tss_down', 'type': 'int', 'default': 1000},
            {"name": "thresh", "type": "float", "default": 0.0001},
            {"name": "pro_dis", "type": "int", "default": 500},
            {"name": "mood_pvalue", "type": "float", "default": 0.0001},
            {"name": "known_species", "type": "string", "default": None},
            {'name': 'gene_link', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'gene_detail', 'type': 'infile', 'format': 'whole_transcriptome.common'},
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
        self.tool_corr = self.add_tool("whole_transcriptome.expcorr_mitarget")
        self.tool_tfbinding = self.add_tool("small_rna.tf_binding")
        self.tool_merge = self.add_tool("small_rna.merge_tf_corr")
        self.smallrna_list = list()
        self.target_list = list()

    def run(self):
        if self.option("exp_target").is_set and self.option("db"):
            self.on_rely([self.tool_tfbinding, self.tool_corr], self.run_merge)
        elif self.option("exp_target").is_set:
            self.tool_corr.on("end", self.run_merge)
        elif self.option("db"):
            self.tool_tfbinding.on("end", self.run_merge)
        else:
            pass

        self.tool_merge.on('end', self.set_db)
        self.get_run_log()
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
        get_run_log = GetRunLog("whole_transcriptome", table="minet", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def get_smallrna_list(self):
        with open(self.option("mirnaset_id"), "r") as mirnaset_f:
            self.smallrna_list = [line.strip() for line in mirnaset_f.readlines()]

    def get_target_list(self):
        with open(self.option("geneset_id"), "r") as geneset_f:
            self.target_list = [line.strip() for line in geneset_f.readlines()]

    def get_non_cor_target(self):
        with open(self.work_dir + "/target_filter.xls", "r") as target_filter_f, open(self.work_dir + "/target_filter_noncor.xls", "w") as target_filter_w:
            for line in target_filter_f:
                cols = line.split("\t")
                target_filter_w.write("\t".join(cols[:5] + [cols[-1]]))

    def filter_pre(self):
        with open(self.option("known_pre").prop["path"], "r") as known_f,  open(self.work_dir + "/known_prefilter.xls", "w") as known_w:
            header = known_f.readline()
            known_w.write(header)
            for line in known_f.readlines():
                cols = line.strip().split("\t")
                if cols[0] in self.smallrna_list:
                    known_w.write(line)
        with open(self.option("novol_pre").prop["path"], "r") as known_f,  open(self.work_dir + "/novol_prefilter.xls", "w") as known_w:
            header = known_f.readline()
            known_w.write(header)
            for line in known_f.readlines():
                cols = line.strip().split("\t")
                if cols[0] in self.smallrna_list:
                    known_w.write(line)


    def filter_target(self):
        if self.option("exp_target").is_set:
            samples = self.option("use_samples").split("|")
            self.logger.info("exp_target is {}".format(self.option("exp_target").prop['path']))
            self.logger.info("target list is {}".format(self.target_list))
            exp = pd.read_table(self.option("exp_target").prop['path'], header=0, index_col=0)
            exp_list = list(exp.index)
            gene_filter = list(set(self.target_list).intersection(set(exp_list)))
            exp_filter = exp.loc[gene_filter, samples]
            exp_filter.to_csv(self.work_dir + "/exp_matrix_filter.xls", sep="\t")
        else:
            pass

        with open(self.option("target").prop['path'], "r") as known_target_f, open(self.work_dir + "/target_filter.xls", "w") as target_filter_w:
            header = known_target_f.readline()
            target_filter_w.write("small_rna\ttarget\tgene\tname\tcategory\tmirtarbase\n")
            pair = set()
            # self.logger.info(" small {} target {} ".format(self.smallrna_list, self.target_list))
            for line in known_target_f.readlines():
                cols = line.strip().split("\t")
                # self.logger.info(" choose {} {}".format( cols[3], cols[4]))
                if cols[3] in self.smallrna_list and cols[4] in self.target_list:
                    if cols[3] + "__" + cols[4] in pair:
                        pass
                    else:
                        target_filter_w.write("\t".join([cols[3], cols[4], cols[1], cols[2], cols[0], "no" + "\n"]))
                        pair.add(cols[3] + "__" + cols[4])

    def run_corr(self):
        options = dict(
            exp=self.option("exp_matrix").prop['path'],
            exp_target = self.work_dir + "/exp_matrix_filter.xls",
            rna_target = self.work_dir + "/target_filter.xls",
            # corr_way = self.option("corr_way"),
            corr_cutoff = self.option("corr_cutoff"),
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
            gene_link = self.option("gene_link").prop['path'],
            gene_detail = self.option("gene_detail").prop['path'],
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


    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        minet_api = self.api.api("whole_transcriptome.minet")
        tf = self.tool_tfbinding.output_dir + "/binding.xls"
        if self.option("exp_target").is_set:
            target = self.tool_corr.output_dir + "/corr.xls"
        else:
            target = self.work_dir + "/target_filter_noncor.xls"
        nodes = self.tool_merge.output_dir + "/nodes.xls"
        edges = self.tool_merge.output_dir + "/edges.xls"
        minet_api.import_net_webroot(self.option("main_id"), self.output_dir, tf=tf, target=target, nodes=nodes, edges=edges)


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

        result_dir.add_relpath_rules([
            [".", "", "mirna关联分析结果目录", 0],
            ["./All_nodes.xls", "xls", "调控网络点信息表", 0],
            ["./All_edges.xls", "xls", "调控网络边信息列表", 0],
            ["./miRNA_mRNA_detail.xls", "xls", "miRNA-mRNA对应关系表", 0],
            ["./TF_miRNA_detail.xls", "xls", "TF-miRNA对应关系表", 0],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        db = Config().get_mongo_client(mtype="whole_transcriptome")[Config().get_mongo_dbname("whole_transcriptome")]
        col1 = db["sh_minet"]
        col1.update({"_id": ObjectId(self.option("main_id"))}, {"$set": {"result_dir": target_dir}}, upsert=True)
        super(MiNetWorkflow, self).end()
