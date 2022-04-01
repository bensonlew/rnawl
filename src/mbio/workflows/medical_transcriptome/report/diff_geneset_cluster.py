# -*- coding: utf-8 -*-
# __author__ = "chenyanyan, 2016.10.12"
# last_modify by fwy 2020120204

from biocluster.workflow import Workflow
import os,re,glob
import pandas as pd
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.medical_transcriptome.chart.chart_geneset import ChartGeneset
import shutil
from collections import OrderedDict


class DiffGenesetClusterWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffGenesetClusterWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name="cluster_main_id", type="string"),
            dict(name="n_clusters", type='int'),
            dict(name="use_group", type="string"),
            dict(name="group", type="string"),
            dict(name="scm", type="string"),
            dict(name="scd", type="string"),
            dict(name="sct", type="string"),
            dict(name="gct", type="string"),
            dict(name="gcm", type="string"),
            dict(name="gcd", type="string"),
            dict(name="gene_detail", type="string"),
            dict(name="group_id", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("medical_transcriptome.exp_cluster")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/01 Diff_Express/02 DiffExpress_ Cluster_Analysis')
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
        super(DiffGenesetClusterWorkflow, self).send_log(data)

    def run(self):
        self.tool.on("end", self.set_db)
        self.get_run_log()
        self.run_tool()
        super(DiffGenesetClusterWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("medical_transcriptome", table="sg_diff_geneset_cluster", main_id=self.option('cluster_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("medical_transcriptome.diff_geneset")
        # add result info
        all_exp.add_geneset_cluster(self.tool.output_dir, self.option("gene_detail"), main_id=self.option('cluster_main_id'), )
        self.end()

    def get_seq_id2name(self):
        seq_id2name = dict()
        with open(self.option("gene_detail"), 'rb') as f:
            for linen in f.readlines()[1:]:
                line = linen.strip("\n")
                if len(line.split("\t")) > 3:
                    if line.split("\t")[2].strip() in ["-", "_", ""]:
                        seq_id2name[line.split("\t")[0]] = line.split("\t")[0]
                    else:
                        seq_id2name[line.split("\t")[0]] = line.split("\t")[2].strip()
                else:
                    if line.split("\t")[1].strip() in ["-", "_", ""]:
                        seq_id2name[line.split("\t")[0]] = line.split("\t")[0]
                    else:
                        seq_id2name[line.split("\t")[0]] = line.split("\t")[1]
        return seq_id2name

    def chart(self):
        seq_id2name = self.get_seq_id2name()
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"
        group_dict = json.loads(self.option("group_dict"), object_pairs_hook=OrderedDict)
        cluster_exp = self.tool.output_dir + "/expression_matrix.xls"
        cluster_tree = self.tool.output_dir + "/seq.cluster_tree.txt"
        sample_tree = self.tool.output_dir + "/sample.cluster_tree.txt"
        subcluster_list = glob.glob(self.tool.output_dir + "/*subcluster_*.xls")
        chart.chart_geneset_cluster(cluster_exp, cluster_tree, sample_tree, subcluster_list, group_dict=group_dict,seq_id2name = seq_id2name)
        chart.to_pdf()

        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            print "copy", p, self.tool.output_dir + "/" + os.path.basename(p)
            shutil.copyfile(p, self.tool.output_dir + "/" + os.path.basename(p))

    def end(self):
        self.chart()
        subcluster = glob.glob(self.tool.output_dir + "/seq.subcluster_*.xls")
        for file in subcluster:
            sub = os.path.basename(file).strip().split("_")[1]
            os.rename(file, os.path.join(os.path.dirname(file), "subcluster_" + sub + ".xls"))
        seq_annot_pd = pd.read_table(self.work_dir + '/seq_annot.xls', header=0, index_col=0)
        seq_matrix_pd = pd.read_table(self.tool.output_dir + '/expression_matrix.xls', header=0, index_col=0)
        seq_matrix_add_annot_pd = seq_matrix_pd.join(seq_annot_pd)
        df1 = seq_matrix_add_annot_pd[["gene_name", "gene_desc"]]
        df2 = seq_matrix_add_annot_pd.drop(columns=["gene_name", "gene_desc"])
        seq_matrix_add_annot_pd = pd.concat([df1, df2], axis=1)
        seq_matrix_add_annot_file = os.path.join(self.tool.output_dir, 'expression_matrix_annot.xls')
        header = ['seq_id']
        header.extend(seq_matrix_add_annot_pd.columns.tolist())
        with open(seq_matrix_add_annot_file, "w") as w:
            w.write("\t".join(header) + "\n")
        seq_matrix_add_annot_pd.to_csv(seq_matrix_add_annot_file, header=False, index=True, sep='\t', mode='a')
        os.rename(os.path.join(self.tool.output_dir, 'expression_matrix_annot.xls'), os.path.join(self.tool.output_dir, 'expression_matrix.xls'))

        df = pd.read_table(os.path.join(self.tool.output_dir, 'expression_matrix.xls'))
        df = df.drop_duplicates()
        df.to_csv(os.path.join(self.tool.output_dir, 'expression_matrix.xls'), sep='\t', index=False)

        if os.path.exists(os.path.join(self.tool.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.tool.output_dir)
        self.inter_dirs = [
            ["01 Diff_Express", "", "差异基因数据挖掘结果目录",0],
            ["01 Diff_Express/02 DiffExpress_ Cluster_Analysis", "", "差异基因聚类分析", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "聚类分析文件",0,"211530"],
            ["./seq.cluster_tree.txt", "txt", "基因/转录本聚类树文件",0,"211531"],
            ["./seq.kmeans_cluster.txt", "txt", "基因/转录本聚类树文件", 0],
            ["./sample.cluster_tree.txt", "txt", "样本聚类树文件",0,"211532"],
            ["./expression_matrix.xls", "xls", "聚类热图分析表",0,"211533"],
            ["./heatmap.pdf", 'pdf', "聚类热图",0],
            ["./cluster.*.pdf", 'pdf', "子聚类趋势图", 0],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        result_dir.add_regexp_rules([
            [r'.*subcluster_.*\.xls', 'xls', '子聚类分析表',0,"211534"],
            [r".*corr\.pdf", 'pdf', "聚类热图", 0]
        ])
        super(DiffGenesetClusterWorkflow, self).end()

    def run_tool(self):
        options = dict(
            exp=self.option('exp_matrix'),
            group=self.option('group'),
            n_clusters=int(self.option('n_clusters')),
            sct=self.option('sct'),
            gct=self.option('gct'),
            scm=self.option('scm'),
            gcm=self.option('gcm'),
            scd=self.option('scd'),
            gcd=self.option('gcd'),
            use_group=self.option('use_group'),
        )
        self.tool.set_options(options)
        self.tool.run()
