# -*- coding: utf-8 -*-
# __author__ = "chenyanyan, 2016.10.12"
# last_modify by khl 20170504

from biocluster.workflow import Workflow
import os,re,glob
import pandas as pd
import json
import shutil
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.ref_rna_v2.chart_geneset import ChartGeneset
from mbio.packages.project_demo.interaction_rerun.interaction_delete import InteractionDelete
from collections import OrderedDict

class GenesetClusterWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetClusterWorkflow, self).__init__(wsheet_object)
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
        self.tool = self.add_tool("ref_rna_v2.exp_cluster")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/03 GeneSet/01 Cluster_Analysis')
        self.inter_dirs = []
        if self._sheet.rerun:
            self.logger.info("该交互为重运行项目,先删除mongo库内容")
            interactiondelete =  InteractionDelete(bind_object=self,project_type="ref_rna_v2",main_id=self.option('cluster_main_id'))
            interactiondelete.delete_interactions_records()

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
        super(GenesetClusterWorkflow, self).send_log(data)

    def run(self):
        self.tool.on("end", self.set_db)
        self.get_run_log()
        self.run_tool()
        super(GenesetClusterWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("ref_rna_v2.all_exp")
        # add result info
        all_exp.add_geneset_cluster(self.tool.output_dir, self.option("gene_detail"), main_id=self.option('cluster_main_id'), )
        self.end()

    def get_run_log(self):
        get_run_log = GetRunLog("ref_rna_v2", table="sg_geneset_cluster", main_id=self.option('cluster_main_id'), dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def chart(self):
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"
        group_dict = json.loads(self.option("group_dict"),object_pairs_hook=OrderedDict)
        cluster_exp = self.tool.output_dir + "/expression_matrix.xls"
        cluster_tree = self.tool.output_dir + "/seq.cluster_tree.txt"
        sample_tree = self.tool.output_dir + "/sample.cluster_tree.txt"
        subcluster_list = glob.glob(self.tool.output_dir + "/*subcluster_*.xls")
        chart.chart_geneset_cluster(cluster_exp, cluster_tree, sample_tree, subcluster_list, group_dict=group_dict)
        chart.to_pdf()

        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            print "copy", p, self.tool.output_dir + "/" + os.path.basename(p)
            # shutil.copyfile(p, self.tool.output_dir + "/" + os.path.basename(p))
            #兼容交互重运行,让文件更改等操作均避免再tool结果文件中进行
            shutil.copyfile(p, self.output_dir + "/" + os.path.basename(p))



    def end(self):
        self.chart()

        subcluster = glob.glob(self.tool.output_dir + "/seq.subcluster_*.xls")
        for file in subcluster:
            sub = os.path.basename(file).strip().split("_")[1]
            if os.path.exists(os.path.join(self.output_dir, "subcluster_" + sub + ".xls")):
                os.remove(os.path.join(self.output_dir, "subcluster_" + sub + ".xls"))
            os.link(file, os.path.join(self.output_dir, "subcluster_" + sub + ".xls"))
        seq_annot_pd = pd.read_table(self.work_dir + '/seq_annot.xls', header=0, index_col=0)
        seq_matrix_pd = pd.read_table(self.tool.output_dir + '/expression_matrix.xls', header=0, index_col=0)
        seq_matrix_add_annot_pd = seq_matrix_pd.join(seq_annot_pd)
        df1 = seq_matrix_add_annot_pd[["gene_name", "gene_desc"]]
        df2 = seq_matrix_add_annot_pd.drop(columns=["gene_name", "gene_desc"])
        seq_matrix_add_annot_pd = pd.concat([df1, df2], axis=1)
        seq_matrix_add_annot_file = os.path.join(self.output_dir, 'expression_matrix.xls')
        header = ['seq_id']
        header.extend(seq_matrix_add_annot_pd.columns.tolist())
        with open(seq_matrix_add_annot_file, "w") as w:
            w.write("\t".join(header) + "\n")
        seq_matrix_add_annot_pd.to_csv(seq_matrix_add_annot_file, header=False, index=True, sep='\t', mode='a')

        df = pd.read_table(os.path.join(self.output_dir, 'expression_matrix.xls'))
        df = df.drop_duplicates()
        df.to_csv(os.path.join(self.output_dir, 'expression_matrix.xls'), sep='\t', index=False)
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))

        #兼容代码以及为了重运行
        if os.path.exists(os.path.join(self.tool.output_dir,"sample.cluster_tree.txt")):
            if os.path.exists(os.path.join(self.output_dir,"sample.cluster_tree.txt")):
                os.remove(os.path.join(self.output_dir,"sample.cluster_tree.txt"))
            os.link(os.path.join(self.tool.output_dir,"sample.cluster_tree.txt"),
                     os.path.join(self.output_dir,"sample.cluster_tree.txt"))
        if os.path.exists(os.path.join(self.tool.output_dir,"seq.cluster_tree.txt")):
            if os.path.exists(os.path.join(self.output_dir,"seq.cluster_tree.txt")):
                os.remove(os.path.join(self.output_dir,"seq.cluster_tree.txt"))
            os.link(os.path.join(self.tool.output_dir,"seq.cluster_tree.txt"),
                     os.path.join(self.output_dir,"seq.cluster_tree.txt"))

        if os.path.exists(os.path.join(self.tool.output_dir,"seq.kmeans_cluster.txt")):
            if os.path.exists(os.path.join(self.output_dir,"seq.kmeans_cluster.txt")):
                os.remove(os.path.join(self.output_dir,"seq.kmeans_cluster.txt"))
            os.link(os.path.join(self.tool.output_dir,"seq.kmeans_cluster.txt"),
                     os.path.join(self.output_dir,"seq.kmeans_cluster.txt"))

        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["03 GeneSet", "", "基因集分析结果目录",0],
            ["03 GeneSet/01 Cluster_Analysis", "", "聚类分析", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "聚类分析文件",0,"211530"],
            ["./seq.cluster_tree.txt", "txt", "基因/转录本聚类树文件",0,"211531"],
            ["./seq.kmeans_cluster.txt", "txt", "基因/转录本聚类树文件", 0],
            ["./sample.cluster_tree.txt", "txt", "样本聚类树文件",0,"211532"],
            ["./expression_matrix.xls", "xls", "聚类热图分析表",0,"211533"],
            ["./heatmap.pdf", 'pdf', "聚类热图",0],
            ["./*heat_corr.pdf", 'pdf', "聚类热图", 0],
            ["./*.line.pdf", 'pdf', "子聚类折线图",0],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        result_dir.add_regexp_rules([
            [r'.*subcluster_.*\.xls', 'xls', '子聚类分析表',0,"211534"],
        ])
        super(GenesetClusterWorkflow, self).end()

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
