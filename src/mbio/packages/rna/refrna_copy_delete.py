# -*- coding: utf-8 -*-
from biocluster.config import Config
from biocluster.api.database.base import Base
from pymongo import MongoClient


class RefrnaCopyDelete(Base):
    def __init__(self):
        super(RefrnaCopyDelete, self).__init__()
        self._project_type = 'ref_rna'

    def remove_collection(self, task_id, main_coll, detatil_coll=[], change_option=''):
        """
        删除主表和对应的详细表
        task_id:task_id，detatil_coll:主表对应的详细表, change_option:详细表中主表字段
        """
        results = self.db[main_coll].find({"task_id": task_id})
        if results:
            for result in results:
                main_id = result["_id"]
                for coll in detatil_coll:
                    items = self.db[coll].remove({change_option: main_id})
                    print "成功删除task_id为{}的细节表{}".format(task_id, coll)
                self.db[main_coll].remove({"_id": result["_id"]})
                print "成功删除task_id为{}的主表{}".format(task_id, main_coll)
        else:
            print "没有找到task_id为{}的主表{}".format(task_id, main_coll)

    def remove(self, task_id):
        self.remove_collection(task_id=task_id, main_coll="sg_task", detatil_coll=[], change_option="")
        self.remove_collection(task_id=task_id, main_coll="sg_specimen", detatil_coll=[], change_option="")
        self.remove_collection(task_id=task_id, main_coll="sg_specimen_graphic", detatil_coll=[], change_option="")
        self.remove_collection(task_id=task_id, main_coll="sg_specimen_group", detatil_coll=[], change_option="")
        self.remove_collection(task_id=task_id, main_coll="sg_specimen_group_compare", detatil_coll=[], change_option="")
        self.remove_collection(task_id=task_id, main_coll="sg_specimen_info", detatil_coll=[], change_option="")
        self.remove_collection(task_id=task_id, main_coll="sg_specimen_mapping", detatil_coll=[], change_option="")
        self.remove_collection(task_id=task_id, main_coll="sg_annotation_stat", detatil_coll=["sg_annotation_stat_detail"], change_option="stat_id")
        self.remove_collection(task_id=task_id, main_coll="sg_express", detatil_coll=["sg_express_detail", "sg_express_gragh", "sg_express_box"], change_option="express_id")
        self.remove_collection(task_id=task_id, main_coll="sg_geneset", detatil_coll=["sg_geneset_detail"], change_option="geneset_id")
        self.remove_collection(task_id=task_id, main_coll="sg_software_para", detatil_coll=[], change_option="")
        self.remove_collection(task_id=task_id, main_coll="sg_annotation_blast", detatil_coll=["sg_annotation_blast_detail"], change_option="blast_id")
        self.remove_collection(task_id=task_id, main_coll="sg_annotation_nr", detatil_coll=["sg_annotation_nr_pie"], change_option="nr_id")
        self.remove_collection(task_id=task_id, main_coll="sg_annotation_swissprot", detatil_coll=["sg_annotation_swissprot_pie"], change_option="swissprot_id")
        self.remove_collection(task_id=task_id, main_coll="sg_annotation_pfam", detatil_coll=["sg_annotation_pfam_bar", "sg_annotation_pfam_detail"], change_option="pfam_id")
        self.remove_collection(task_id=task_id, main_coll="sg_annotation_cog", detatil_coll=["sg_annotation_cog_detail"], change_option="cog_id")
        self.remove_collection(task_id=task_id, main_coll="sg_annotation_go", detatil_coll=["sg_annotation_go_detail", "sg_annotation_go_graph", "sg_annotation_go_level", "sg_annotation_go_list"], change_option="go_id")
        self.remove_collection(task_id=task_id, main_coll="sg_annotation_kegg", detatil_coll=["sg_annotation_kegg_categories", "sg_annotation_kegg_level", "sg_annotation_kegg_table"], change_option="kegg_id")
        self.remove_collection(task_id=task_id, main_coll="sg_annotation_query", detatil_coll=["sg_annotation_query_detail"], change_option="query_id")
        self.remove_collection(task_id=task_id, main_coll="sg_assessment_chrom_distribution", detatil_coll=["sg_assessment_chrom_distribution_detail"], change_option="chrom_distribution_id")
        self.remove_collection(task_id=task_id, main_coll="sg_assessment_coverage", detatil_coll=["sg_assessment_coverage_detail"], change_option="coverage_id")
        self.remove_collection(task_id=task_id, main_coll="sg_assessment_distribution", detatil_coll=["sg_assessment_distribution_detail"], change_option="distribution_id")
        self.remove_collection(task_id=task_id, main_coll="sg_assessment_duplicate", detatil_coll=["sg_assessment_duplicate_detail"], change_option="dup_id")
        self.remove_collection(task_id=task_id, main_coll="sg_assessment_saturation", detatil_coll=["sg_assessment_saturation_curve"], change_option="saturation_id")
        self.remove_collection(task_id=task_id, main_coll="sg_express_diff", detatil_coll=["sg_express_diff_detail", "sg_express_diff_summary"], change_option="express_diff_id")
        self.remove_collection(task_id=task_id, main_coll="sg_express_correlation", detatil_coll=["sg_express_correlation_detail"], change_option="correlation_id")
        self.remove_collection(task_id=task_id, main_coll="sg_express_pca", detatil_coll=["sg_express_pca_rotation"], change_option="pca_id")
        self.remove_collection(task_id=task_id, main_coll="sg_express_venn", detatil_coll=["sg_express_venn_detail", "sg_express_venn_graph"], change_option="venn_id")
        self.remove_collection(task_id=task_id, main_coll="sg_express_class_code", detatil_coll=["sg_express_class_code_detail"], change_option="class_code_id")
        self.remove_collection(task_id=task_id, main_coll="sg_geneset_venn", detatil_coll=["sg_geneset_venn_detail", "sg_geneset_venn_graph"], change_option="venn_id")
        self.remove_collection(task_id=task_id, main_coll="sg_geneset_cluster", detatil_coll=["sg_geneset_cluster_detail"], change_option="cluster_id")
        self.remove_collection(task_id=task_id, main_coll="sg_geneset_cog_class", detatil_coll=["sg_geneset_cog_class_detail"], change_option="geneset_cog_id")
        self.remove_collection(task_id=task_id, main_coll="sg_geneset_go_class", detatil_coll=["sg_geneset_go_class_detail"], change_option="go_regulate_id")
        self.remove_collection(task_id=task_id, main_coll="sg_geneset_go_enrich", detatil_coll=["sg_geneset_go_enrich_detail"], change_option="go_enrich_id")
        self.remove_collection(task_id=task_id, main_coll="sg_geneset_kegg_class", detatil_coll=["sg_geneset_kegg_class_detail", "sg_geneset_kegg_class_pathway"], change_option="kegg_id")
        self.remove_collection(task_id=task_id, main_coll="sg_geneset_kegg_enrich", detatil_coll=["sg_geneset_kegg_enrich_detail"], change_option="kegg_enrich_id")
        self.remove_collection(task_id=task_id, main_coll="sg_ppinetwork", detatil_coll=["sg_ppinetwork_centrality_node", "sg_ppinetwork_distribution_node", "sg_ppinetwork_node_table", "sg_ppinetwork_structure_attributes", "sg_ppinetwork_structure_link", "sg_ppinetwork_structure_node"], change_option="ppi_id")
        self.remove_collection(task_id=task_id, main_coll="sg_snp", detatil_coll=["sg_snp_detail", "sg_snp_stat"], change_option="snp_id")
        self.remove_collection(task_id=task_id, main_coll="sg_species_information", detatil_coll=["sg_species_information_detail"], change_option="species_id")
        self.remove_collection(task_id=task_id, main_coll="sg_splicing_rmats", detatil_coll=["sg_splicing_rmats_detail", "sg_splicing_rmats_graph", "sg_splicing_rmats_psi", "sg_splicing_rmats_stats"], change_option="splicing_id")
        self.remove_collection(task_id=task_id, main_coll="sg_transcripts", detatil_coll=["sg_transcripts_seq_type", "sg_transcripts_step", "sg_transcripts_relations"], change_option="transcripts_id")

    def find_task_id(self, task_id):
        results = self.db["sg_task"].find({"task_id": {"$regex": task_id + "_.*_.*"}})
        if results:
            for result in results:
                target_task_id = result["task_id"]
                self.remove(target_task_id)
        else:
            print "没有找到以task_id为{}备份的demo,请检查!".format(task_id)


if __name__ == "__main__":
    test = RefrnaCopyDelete()
    test.find_task_id(task_id="demo_zj")
