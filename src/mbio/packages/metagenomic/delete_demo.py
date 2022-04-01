# -*- coding: utf-8 -*-
from biocluster.config import Config
from pymongo import MongoClient
from biocluster.api.database.base import Base


class DeleteDemo(Base):
    def __init__(self):
        super(DeleteDemo, self).__init__()
        self._project_type = 'metagenomic'

    def remove_collection(self, task_id, main_coll, detail_coll=[], change_option=''):
        """
        删除主表和对应的详细表
        task_id:task_id，detail_coll:主表对应的详细表, change_option:详细表中主表字段
        """
        results = self.db[main_coll].find({"task_id": task_id})
        if results:
            for result in results:
                main_id = result["_id"]
                for coll in detail_coll:
                    try:
                        items = self.db[coll].remove({change_option: main_id})
                        print "成功删除task_id为 {} 的细节表 {}".format(task_id, coll)
                    except:
                        print "删除task_id为 {} 的细节表 {} 失败".format(task_id, coll)
                self.db[main_coll].remove({"_id": result["_id"]})
                print "成功删除task_id为 {} 的主表 {}".format(task_id, main_coll)
        else:
            print "没有找到task_id为 {} 的主表 {}".format(task_id, main_coll)

    def remove(self, task_id):
        self.remove_collection(task_id=task_id, main_coll="sg_task", detail_coll=[], change_option="")
        self.remove_collection(task_id=task_id, main_coll="data_stat", detail_coll=["data_stat_detail"], change_option="data_stat_id")
        self.remove_collection(task_id=task_id, main_coll="specimen_group", detail_coll=[], change_option="")
        self.remove_collection(task_id=task_id, main_coll="specimen_graphic", detail_coll=[], change_option="")
        self.remove_collection(task_id=task_id, main_coll="env", detail_coll=["env_detail"], change_option="env_id")
        self.remove_collection(task_id=task_id, main_coll="assemble_stat", detail_coll=["assemble_stat_bar", "assemble_stat_detail"], change_option="assem_id")
        self.remove_collection(task_id=task_id, main_coll="predict_gene", detail_coll=["predict_gene_bar", "predict_gene_detail"], change_option="predict_gene_id")
        self.remove_collection(task_id=task_id, main_coll="predict_gene_total", detail_coll=[], change_option="")
        self.remove_collection(task_id=task_id, main_coll="geneset", detail_coll=["geneset_bar", "geneset_readsn", "geneset_readsr"], change_option="geneset_id")
        anno_list = ["ardb", "card", "cazy", "cog", "kegg", "nr", "vfdb", "overview"]
        anno_detail = {
            "ardb": ["ardb_arg", "ardb_class", "ardb_type"],
            "card": ["card_aro", "card_class"],
            "cazy": ["cazy_class", "cazy_family"],
            "cog": ["cog_category", "cog_function", "cog_nog"],
            "kegg": ["kegg_enzyme", "kegg_gene", "kegg_module", "kegg_orthology", "kegg_pathway"],
            "nr": ["nr_detail"],
            "vfdb": ["vfdb_pie", "vfdb_vfs"],
            "overview": [],
        }
        for anno in anno_list:
            collection = "anno_" + anno
            position = "" if anno == "overview" else anno + "_id"
            self.remove_collection(task_id=task_id, main_coll=collection, detail_coll=anno_detail[anno], change_option=position)
        self.remove_collection(task_id=task_id, main_coll="hcluster_tree", detail_coll=[], change_option="")
        self.remove_collection(task_id=task_id, main_coll="anosim", detail_coll=["anosim_detail"], change_option="anosim_id")
        self.remove_collection(task_id=task_id, main_coll="beta_diversity", detail_coll=["beta_diversity_detail"], change_option="beta_diversity_id")
        self.remove_collection(task_id=task_id, main_coll="composition", detail_coll=["composition_detail"], change_option="composition_id")
        self.remove_collection(task_id=task_id, main_coll="enterotype", detail_coll=['enterotype_detail','enterotype_detail_cluster'], change_option="enterotype_id")
        self.remove_collection(task_id=task_id, main_coll="env_vif", detail_coll=["env_vif_detail"], change_option="vif_id")
        self.remove_collection(task_id=task_id, main_coll="heatmap_cor", detail_coll=["heatmap_cor_detail"], change_option="heatmap_cor_id")
        self.remove_collection(task_id=task_id, main_coll="lefse", detail_coll=["lefse_detail"], change_option="species_lefse_id")
        self.remove_collection(task_id=task_id, main_coll="mantel_test", detail_coll=['mantel_test_detail'], change_option="mantel_id")
        self.remove_collection(task_id=task_id, main_coll="metastat", detail_coll=['metastat_detail', 'metastat_plot'], change_option="metastat_id")
        self.remove_collection(task_id=task_id, main_coll="network", detail_coll=['network_degree', 'network_link', 'network_node'], change_option="network_id")
        self.remove_collection(task_id=task_id, main_coll="network_cor", detail_coll=['network_cor_degree', 'network_cor_link', 'network_cor_node'], change_option="network_cor_id")
        self.remove_collection(task_id=task_id, main_coll="permanova", detail_coll=["permanova_detail"], change_option="permanova_id")
        self.remove_collection(task_id=task_id, main_coll="regression", detail_coll=["regression_curve", "regression_line"], change_option="regression_id")
        self.remove_collection(task_id=task_id, main_coll="specimen_distance", detail_coll=["specimen_distance_detail"], change_option="specimen_distance_id")
        self.remove_collection(task_id=task_id, main_coll="contribute", detail_coll=["contribute_detail"], change_option="contribute_id")

    def find_task_id(self, task_id):
        # results = self.db["sg_task"].find({"task_id": {"$regex": task_id + "_.*_.*"}})
        results = self.db["sg_task"].find({"task_id": task_id})
        if results or results.count() != 0:
            for result in results:
                target_task_id = result["task_id"]
                print "开始删除"
                self.remove(target_task_id)
        else:
            print "没有找到以task_id为{}备份的demo,请检查!".format(task_id)


if __name__ == "__main__":
    test = DeleteDemo()
    test.find_task_id(task_id="tsg_28198")
