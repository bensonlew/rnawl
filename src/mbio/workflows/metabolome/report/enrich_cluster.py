
# -*- coding: utf-8 -*-


from biocluster.workflow import Workflow
from biocluster.config import Config
import glob
import os
from bson.objectid import ObjectId


class EnrichClusterWorkflow(Workflow):
    """
    基因集富集聚类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        # self.rpc = False
        super(EnrichClusterWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "anno_overview", "type": "infile", "format": "metabolome.overview_table"},
            {"name": "ko_overview", "type": "infile", "format": "sequence.profile_table"},  # 代谢集总览表ko表 by ghd @20191015
            {"name": "metabset", "type": "infile", "format": "metabolome.mul_metabset"},
            {"name": "correct", "type": "string", "default": "BH"},  #bonferroni,holm,BY,BH
            {"name": "bg", "type": "string", "default": "project"},  #project,species,kegg
             # 背景，project:本项目鉴定到的代谢物合集; species:本物种全部代谢物合集; kegg:KEGG数据库全部代谢物合集
            {"name": "species", "type": "string", "default": "all"},
            {"name": "select", "type": "string", "default": "pvalue"},  # 挑选那列做聚类分析 pvalue, qvaule,in_pop,in_study
            {"name": "pathway_method", "type": "string", "default": "hierarchy"}, #hierarchy, kmeans
            {"name": "pathway_dist", "type": "string", "default": "euclidean"},
            {"name": "pathway_ctype", "type": "string", "default": "complete"},
            {"name": "n_cluster", "type": "int", "default": 10},
            {"name": "metabset_method", "type": "string", "default": "hierarchy"},
            {"name": "metabset_dist", "type": "string", "default": "euclidean"},
            {"name": "metabset_ctype", "type": "string", "default": "complete"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "scale", "type": "string", "default": "T"},  ## 是否标准化聚类数据，默认做标准化 T ，不做： F
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
            {"name": "name", "type": "string"},

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.enrich_cluster_module = self.add_module("metabolome.enrich_cluster")
        self.output_dir = self.enrich_cluster_module.output_dir

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        options = {
            "anno_overview": self.option('anno_overview'),
            "ko_overview": self.option("ko_overview"),
            "metabset": self.option("metabset"),
            "correct": self.option("correct"),
            "bg": self.option("bg"),
            "species":  self.option("species"),
            "sct" : self.option("metabset_method"),
            "scd" : self.option("metabset_dist"),
            "scm" : self.option("metabset_ctype"),
            "mct" : self.option("pathway_method"),
            "mcd" : self.option("pathway_dist"),
            "mcm" : self.option("pathway_ctype"),
            "n_cluster" : self.option("n_cluster"),
            "select" : self.option('select'),
            "scale" : self.option("scale"),
            "task_id": "_".join(self._sheet.id.split("_")[0:2])
        }
        self.enrich_cluster_module.set_options(options)
        self.enrich_cluster_module.on('end', self.set_db)
        self.enrich_cluster_module.run()
        super(EnrichClusterWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """

        api_name = self.api.api('metabolome.enrich_cluster')

        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
             main_id = ObjectId(main_id)

        metabset_tree = None
        pathway_tree = None
        if os.path.exists(self.output_dir + "/metabset.cluster_tree.xls"):
            metabset_tree = os.path.join(self.output_dir,"metabset.cluster_tree.xls")
        if os.path.exists(self.output_dir + "/pathway.cluster_tree.xls"):
            pathway_tree = os.path.join(self.output_dir,"pathway.cluster_tree.xls")


        expression_file = os.path.join(self.output_dir,"cluster_exp.xls")
        api_name.add_enrich_cluster(metabset_tree=metabset_tree, pathway_tree=pathway_tree, main_id=main_id)
        pvlaue_file =  os.path.join(self.output_dir,"enrich_pvalue.xls")
        metab_ids_file = os.path.join(self.enrich_cluster_module.work_dir, 'metab_ids.xls')
        api_name.add_enrich_cluster_detail(main_id, expression_file, pvlaue_file, metab_ids_file)
        # self.option("save_pdf", 0)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "metabsetkeggheatmap",
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        if self.option("save_pdf"):
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        relpath_rules =[
            [".", "", "代谢物聚类结果文件夹", 0],
            ["metabset.cluster_tree.xls", "xls", "代谢集聚类树文件", 0],
            ["pathway.cluster_tree.xls", "xls", "通路聚类树文件", 0],
            ["pathway.kmeans_cluster.xls", "xls", "通路kmeans分类文件", 0],
            ["cluster_exp.xls", "xls", "矩阵值", 0]
        ]
        regexps = [
            [r"pathway.subcluster_.*\.xls", "xls", "通路各子类结果表", 0],
        ]
        result_dir.add_relpath_rules(relpath_rules)
        result_dir.add_regexp_rules(regexps)

        super(EnrichClusterWorkflow, self).end()



if __name__ == "__main__":
    from biocluster.wsheet import Sheet
    data = {
        'name': 'test_enrich_cluster',
        'id' : "tsg_36965",
        'type' : "workflow",
        "options" : {
            "anno_overview" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Anno/AnnoOverview/output/anno.xls2",
            "ko_overview" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Anno/AnnoOverview/output/ko.xls",
            "metabset" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/venn_input.xls",
            "main_table_id" : "5e58c2ba17b2bf54005e3145",
            "select" : "in_pop"
        }

    }
    data = Sheet(data=data)
    wf = EnrichClusterWorkflow(data)
    wf.run()