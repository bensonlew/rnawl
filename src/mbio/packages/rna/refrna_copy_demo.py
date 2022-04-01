# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modifiy: 2017.06.05
import json
import os
import gevent
import datetime
import time
import sys
from bson import ObjectId
from gevent import Greenlet
#from gevent.monkey import patch_all
from biocluster.config import Config
from biocluster.api.database.base import Base
from mainapp.libs.param_pack import group_detail_sort


class RefrnaCopyMongo(Base):
    def __init__(self, old_task_id, new_task_id, new_project_sn, new_member_id, new_bam_path=None, new_ref_gtf=None, db='sanger_ref_rna'):
        # self.db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        super(RefrnaCopyMongo, self).__init__()
        self._project_type = 'ref_rna'
        self._old_task_id = old_task_id
        self._new_task_id = new_task_id
        self._new_project_sn = new_project_sn
        self._new_member_id = new_member_id
        self._new_bam_path = new_bam_path
        self._new_ref_gtf = new_ref_gtf
        self.specimen_id_dict = {}
        self.specimen_group_id_dict = {}
        self.specimen_compare_id_dict = {}
        self.stat_id_dict = {}
        self.express_id_dict = {}
        self.geneset_id_dict = {}
        self.nr_id_dict = {}
        self.all_greenlets = []
        self._exchange_dict = {  # 根据特定字段名称，进行特定的ID新旧替换
            'specimen_id': self.specimen_id_dict,
            'specimen_group_id': self.specimen_group_id_dict,
            'compare_id': self.specimen_compare_id_dict,
            'stat_id': self.stat_id_dict,
            'express_id': self.express_id_dict,
            'geneset_id': self.geneset_id_dict
            }

    def run(self):
        """
        运行执行复制特定ID数据的操作，如果有新的分析请参照下面的写法添加代码，不同分析表结构不同，所有需要手动添加。
        """
        #patch_all()
        self.copy_member_id()
        self.copy_sg_specimen()  # specimen_id
        self.copy_collection_with_change('sg_specimen_graphic', change_positions=['specimen_id'])
        self.copy_sg_specimen_group()  # specimen_group_id
        self.copy_sg_specimen_group_compare()  # compare_id
        self.stat_id_dict = self.copy_sg_annotation_stat()  # stat_id
        self.copy_main_details('sg_annotation_stat_detail', 'stat_id', self.stat_id_dict, join=False)
        self.express_id_dict = self.copy_sg_express()  # express_id
        self.geneset_id_dict = self.copy_sg_geneset()  # geneset_id
        self.copy_main_details("sg_geneset_detail", "geneset_id", self.geneset_id_dict, join=False)
        self.copy_collection_with_change('sg_specimen_mapping')
        self.copy_collection_with_change('sg_specimen_info')
        self.copy_collection_with_change('sg_software_para')
        greenlet = Greenlet(self.annotation_blast)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.annotation_nr)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.annotation_swissprot)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.annotation_pfam)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.annotation_cog)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.annotation_go)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.annotation_kegg)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.annotation_query)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.assessment_chrom_distribution)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.assessment_coverage)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.assessment_distribution)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.assessment_duplicate)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.assessment_saturation)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.express)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.express_class_code)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.express_corr)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.express_diff)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.express_pca)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.express_venn)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.geneset_cluster)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.geneset_cog_class)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.geneset_go_class)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.geneset_go_enrich)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.geneset_kegg_class)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.geneset_kegg_enrich)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.geneset_venn)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.ppinetwork)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.snp)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.species_information)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.splicing_rmats)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        greenlet = Greenlet(self.transcripts)
        greenlet.start()
        self.all_greenlets.append(greenlet)
        gevent.joinall(self.all_greenlets)
        gevent.joinall(self.all_greenlets)
        import socket
        reload(socket)

    def annotation_blast(self):
        annotation_blast_dict = self.copy_collection_with_change('sg_annotation_blast', change_positions=['stat_id'], update_sg_status=False)
        self.copy_main_details('sg_annotation_blast_detail', 'blast_id', annotation_blast_dict, join=False)

    def annotation_nr(self):
        annotation_nr_dict = self.copy_collection_with_change('sg_annotation_nr', change_positions=['stat_id'], update_sg_status=False)
        self.copy_main_details('sg_annotation_nr_pie', 'nr_id', annotation_nr_dict, join=False)

    def annotation_swissprot(self):
        annotation_swissprot_dict = self.copy_collection_with_change('sg_annotation_swissprot', change_positions=['stat_id'], update_sg_status=False)
        self.copy_main_details('sg_annotation_swissprot_pie', 'swissprot_id', annotation_swissprot_dict, join=False)

    def annotation_pfam(self):
        annotation_pfam_dict = self.copy_collection_with_change('sg_annotation_pfam', change_positions=[], update_sg_status=False)
        self.copy_main_details('sg_annotation_pfam_bar', 'pfam_id', annotation_pfam_dict, join=False)
        self.copy_main_details('sg_annotation_pfam_detail', 'pfam_id', annotation_pfam_dict, join=False)

    def annotation_cog(self):
        annotation_cog_dict = self.copy_collection_with_change('sg_annotation_cog', change_positions=[], update_sg_status=False)
        self.copy_main_details('sg_annotation_cog_detail', 'cog_id', annotation_cog_dict, join=False)
        # self.copy_main_details('sg_annotation_cog_table', 'cog_id', annotation_cog_dict, join=False)

    def annotation_go(self):
        annotation_go_dict = self.copy_collection_with_change('sg_annotation_go', change_positions=[], update_sg_status=False)
        self.copy_main_details('sg_annotation_go_detail', 'go_id', annotation_go_dict, join=False)
        self.copy_main_details('sg_annotation_go_graph', 'go_id', annotation_go_dict, join=False)
        self.copy_main_details('sg_annotation_go_level', 'go_id', annotation_go_dict, join=False)
        self.copy_main_details('sg_annotation_go_list', 'go_id', annotation_go_dict, join=False)

    def annotation_kegg(self):
        annotation_kegg_dict = self.copy_collection_with_change('sg_annotation_kegg', change_positions=[], update_sg_status=False)
        self.copy_main_details('sg_annotation_kegg_categories', 'kegg_id', annotation_kegg_dict, join=False)
        self.copy_main_details('sg_annotation_kegg_level', 'kegg_id', annotation_kegg_dict, join=False)
        self.copy_main_details('sg_annotation_kegg_table', 'kegg_id', annotation_kegg_dict, join=False)

    def annotation_query(self):
        annotation_query_dict = self.copy_collection_with_change('sg_annotation_query', change_positions=[], update_sg_status=False)
        # self.copy_main_details('sg_annotation_query_detail', 'query_id', annotation_query_dict, join=False)

    def assessment_chrom_distribution(self):
        assessment_chrom_distribution_dict = self.copy_collection_with_change('sg_assessment_chrom_distribution', change_positions=[], update_sg_status=False)
        self.copy_main_details('sg_assessment_chrom_distribution_detail', 'chrom_distribution_id', assessment_chrom_distribution_dict, join=False)
        # self.copy_main_details('sg_assessment_chrom_distribution_circos', 'chrom_distribution_id', assessment_chrom_distribution_dict, join=False)

    def assessment_coverage(self):
        assessment_coverage_dict = self.copy_collection_with_change('sg_assessment_coverage', change_positions=[], update_sg_status=False)
        self.copy_main_details('sg_assessment_coverage_detail', 'coverage_id', assessment_coverage_dict, join=False)

    def assessment_distribution(self):
        assessment_distribution_dict = self.copy_collection_with_change('sg_assessment_distribution', change_positions=[], update_sg_status=False)
        self.copy_main_details('sg_assessment_distribution_detail', 'distribution_id', assessment_distribution_dict, join=False)

    def assessment_duplicate(self):
        assessment_duplicate_dict = self.copy_collection_with_change('sg_assessment_duplicate', change_positions=[], update_sg_status=False)
        self.copy_main_details('sg_assessment_duplicate_detail', 'dup_id', assessment_duplicate_dict, join=False)

    def assessment_saturation(self):
        assessment_saturation_dict = self.copy_collection_with_change('sg_assessment_saturation', change_positions=[], update_sg_status=False)
        self.copy_main_details('sg_assessment_saturation_curve', 'saturation_id', assessment_saturation_dict, join=False)

    def express(self):
        # self.express_id_dict = {"596da5cfa4e1af72a2c853b4": "5982d445a4e1af3e78860d94", "596da83ea4e1af74636a93f1": "5982d445a4e1af3e78860d95", "596da9fca4e1af762e41a251": "5982d445a4e1af3e78860d96", "596dad0ba4e1af7855c0ef83": "5982d445a4e1af3e78860d97"}
        self.copy_main_details("sg_express_detail", 'express_id', self.express_id_dict, join=False)
        self.copy_main_details("sg_express_gragh", "express_id", self.express_id_dict, join=False)
        self.copy_main_details("sg_express_box", "express_id", self.express_id_dict, join=False)

    def express_diff(self):
        express_diff_dict = self.copy_collection_with_change("sg_express_diff", change_positions=["express_id"], update_sg_status=False)
        self.copy_main_details("sg_express_diff_detail", "express_diff_id", express_diff_dict, join=False)
        self.copy_main_details("sg_express_diff_summary", "express_diff_id", express_diff_dict, join=False)

    def express_corr(self):
        express_corr_dict = self.copy_collection_with_change("sg_express_correlation", change_positions=[], update_sg_status=False)
        self.copy_main_details("sg_express_correlation_detail", "correlation_id", express_corr_dict, join=False)

    def express_pca(self):
        express_pca_dict = self.copy_collection_with_change("sg_express_pca", change_positions=[], update_sg_status=False)
        self.copy_main_details("sg_express_pca_rotation", "pca_id", express_pca_dict, join=False)

    def express_venn(self):
        express_venn_dict = self.copy_collection_with_change("sg_express_venn", change_positions=[], update_sg_status=False)
        self.copy_main_details("sg_express_venn_detail", "venn_id", express_venn_dict, join=False)
        self.copy_main_details("sg_express_venn_graph", "venn_id", express_venn_dict, join=False)

    def express_class_code(self):
        class_code_dict = self.copy_collection_with_change("sg_express_class_code", change_positions=[], update_sg_status=False)
        # self.copy_main_details("sg_express_class_code_detail", "class_code_id", class_code_dict, join=False)

    def geneset_venn(self):
        geneset_venn_dict = self.copy_collection_with_change("sg_geneset_venn", change_positions=[],update_sg_status=False)
        self.copy_main_details("sg_geneset_venn_detail", "venn_id", geneset_venn_dict, join=False)
        self.copy_main_details("sg_geneset_venn_graph", "venn_id", geneset_venn_dict, join=False)

    def geneset_cluster(self):
        geneset_cluster_dict = self.copy_collection_with_change("sg_geneset_cluster", change_positions=["express_id"], update_sg_status=False)
        self.copy_main_details("sg_geneset_cluster_detail", "cluster_id", geneset_cluster_dict, join=False)

    def geneset_cog_class(self):
        geneset_cog_class_dict = self.copy_collection_with_change("sg_geneset_cog_class", change_positions=[], update_sg_status=False)
        self.copy_main_details("sg_geneset_cog_class_detail", "geneset_cog_id", geneset_cog_class_dict, join=False)

    def geneset_go_class(self):
        geneset_go_class_dict = self.copy_collection_with_change("sg_geneset_go_class", change_positions=[], update_sg_status=False)
        self.copy_main_details("sg_geneset_go_class_detail", "go_regulate_id", geneset_go_class_dict, join=False)

    def geneset_go_enrich(self):
        geneset_go_enrich_dict = self.copy_collection_with_change("sg_geneset_go_enrich", change_positions=[], update_sg_status=False)
        self.copy_main_details("sg_geneset_go_enrich_detail", "go_enrich_id", geneset_go_enrich_dict, join=False)

    def geneset_kegg_class(self):
        geneset_kegg_class_dict = self.copy_collection_with_change("sg_geneset_kegg_class", change_positions=[], update_sg_status=False)
        self.copy_main_details("sg_geneset_kegg_class_detail", "kegg_id", geneset_kegg_class_dict, join=False)
        self.copy_main_details("sg_geneset_kegg_class_pathway", "kegg_id", geneset_kegg_class_dict, join=False)

    def geneset_kegg_enrich(self):
        geneset_kegg_enrich_dict = self.copy_collection_with_change("sg_geneset_kegg_enrich", change_positions=[], update_sg_status=False)
        self.copy_main_details("sg_geneset_kegg_enrich_detail", "kegg_enrich_id", geneset_kegg_enrich_dict, join=False)

    def ppinetwork(self):
        ppinetwork_dict = self.copy_collection_with_change("sg_ppinetwork", change_positions=["geneset_id"], update_sg_status=False)
        self.copy_main_details("sg_ppinetwork_centrality_node", "ppi_id", ppinetwork_dict, join=False)
        self.copy_main_details("sg_ppinetwork_distribution_node", "ppi_id", ppinetwork_dict, join=False)
        self.copy_main_details("sg_ppinetwork_node_table", "ppi_id", ppinetwork_dict, join=False)
        self.copy_main_details("sg_ppinetwork_structure_attributes", "ppi_id", ppinetwork_dict, join=False)
        self.copy_main_details("sg_ppinetwork_structure_link", "ppi_id", ppinetwork_dict, join=False)
        self.copy_main_details("sg_ppinetwork_structure_node", "ppi_id", ppinetwork_dict, join=False)

    def snp(self):
        snp_dict = self.copy_collection_with_change('sg_snp', change_positions=[], update_sg_status=False)
        # self.copy_main_details('sg_snp_detail', 'snp_id', snp_dict, join=False)
        # # self.copy_main_details('sg_snp_graphic', 'snp_id', snp_dict, join=False)
        # self.copy_main_details('sg_snp_stat', 'snp_id', snp_dict, join=False)

    def species_information(self):
        species_information_dict = self.copy_collection_with_change('sg_species_information', change_positions=[], update_sg_status=False)
        self.copy_main_details('sg_species_information_detail', 'species_id', species_information_dict, join=False)

    def splicing_rmats(self):
        splicing_rmats_dict = self.copy_collection_with_change('sg_splicing_rmats', change_positions=[], update_sg_status=False)
        self.update_sg_splicing_rmats()
        self.copy_main_details('sg_splicing_rmats_detail', 'splicing_id', splicing_rmats_dict, join=False)
        self.copy_main_details('sg_splicing_rmats_graph', 'splicing_id', splicing_rmats_dict, join=False)
        self.copy_main_details('sg_splicing_rmats_psi', 'splicing_id', splicing_rmats_dict, join=False)
        # self.copy_main_details('sg_splicing_rmats_model', 'splicing_id', splicing_rmats_dict, join=False)
        self.copy_main_details('sg_splicing_rmats_stats', 'splicing_id', splicing_rmats_dict, join=False)

    def transcripts(self):
        transcripts_dict = self.copy_collection_with_change('sg_transcripts', change_positions=[], update_sg_status=False)
        self.copy_main_details('sg_transcripts_seq_type', 'transcripts_id', transcripts_dict, join=False)
        self.copy_main_details('sg_transcripts_step', 'transcripts_id', transcripts_dict, join=False)
        self.copy_main_details('sg_transcripts_relations', 'transcripts_id', transcripts_dict, join=False)

    def copy_member_id(self):
        """
        复制sg_task的数据
        """
        coll = self.db["sg_task"]
        find = coll.find_one({'task_id': self._old_task_id})
        if not find:
            raise Exception('运行错误：找不到demo任务相关信息')
        find['task_id'] = self._new_task_id
        find['member_id'] = self._new_member_id
        find.pop('_id')
        find['project_sn'] = self._new_project_sn
        find['is_demo'] = 1
        try:
            find['demo_id'] = self._old_task_id
        except:
            find.pop('demo_id')
            find['demo_id'] = self._old_task_id
        self.db["sg_task"].insert_one(find)

    def copy_sg_specimen(self):
        """
        复制样本表
        """
        finds = self.db["sg_specimen"].find({"task_id": self._old_task_id})
        news = []
        old_specimen_ids = []
        for i in finds:
            old_specimen_ids.append(str(i.pop('_id')))
            old_bam_path = os.path.basename(i["bam_path"])
            i['task_id'] = self._new_task_id
            i['project_sn'] = self._new_project_sn
            if self._new_bam_path:
                i['bam_path'] = os.path.join(self._new_bam_path, old_bam_path)
            news.append(i)
        if news:
            result = self.db["sg_specimen"].insert_many(news)
        else:
            raise Exception('没有任何样本信息，请核对任务结果是否完整')
        self.specimen_id_dict = dict(zip(old_specimen_ids, [str(one) for one in result.inserted_ids]))
        self._exchange_dict['specimen_id'] = self.specimen_id_dict
        return self.specimen_id_dict

    def copy_sg_specimen_group(self):
        """
        复制sg_specimen_group表
        """
        finds = self.db["sg_specimen_group"].find({"task_id": self._old_task_id})
        news = []
        old_group_ids = []
        for i in finds:
            i['task_id'] = self._new_task_id
            i['project_sn'] = self._new_project_sn
            old_group_ids.append(str(i.pop('_id')))
            for one in i['specimen_names']:
                for sp in one.copy():
                    one[self.specimen_id_dict[sp]] = one[sp]
                    one.pop(sp)
            news.append(i)
        if news:
            result = self.db.sg_specimen_group.insert_many(news)
            self.specimen_group_id_dict = dict(zip(old_group_ids, [str(one) for one in result.inserted_ids]))
        else:
            self.specimen_group_id_dict = {}
        self.specimen_group_id_dict['all'] = 'all'  # 特殊ID
        self.specimen_group_id_dict[None] = None  # 特殊ID
        self.specimen_group_id_dict[''] = None  # 特殊ID
        self._exchange_dict['specimen_group_id'] = self.specimen_group_id_dict
        return self.specimen_group_id_dict

    def copy_sg_specimen_group_compare(self):
        """
        复制sg_specimen_group_compare
        """
        finds = self.db["sg_specimen_group_compare"].find({"task_id": self._old_task_id})
        news = []
        old_compare_ids = []
        for i in finds:
            old_compare_ids.append(str(i.pop('_id')))
            i['task_id'] = self._new_task_id
            i['project_sn'] = self._new_project_sn
            if "specimen_group_id" in i:
                i["specimen_group_id"] = self.specimen_group_id_dict[i["specimen_group_id"]]
            news.append(i)
        if news:
            result = self.db["sg_specimen_group_compare"].insert_many(news)
        else:
            raise Exception('没有任何样本分组对照信息，请核对任务结果是否完整')
        self.specimen_compare_id_dict = dict(zip(old_compare_ids, [str(one) for one in result.inserted_ids]))
        self._exchange_dict['compare_id'] = self.specimen_compare_id_dict
        return self.specimen_compare_id_dict

    def copy_sg_annotation_stat(self):
        """
        复制sg_annotation_stat
        """
        finds = self.db["sg_annotation_stat"].find({"task_id": self._old_task_id})
        news = []
        old_stat_ids = []
        for i in finds:
            i['task_id'] = self._new_task_id
            i['project_sn'] = self._new_project_sn
            old_stat_ids.append(str(i.pop('_id')))
            news.append(i)
        if news:
            result = self.db["sg_annotation_stat"].insert_many(news)
        else:
            raise Exception('没有注释统计信息，请核对任务结果是否完整')
        self.stat_id_dict = dict(zip(old_stat_ids, [str(one) for one in result.inserted_ids]))
        self._exchange_dict['stat_id'] = self.stat_id_dict
        return self.stat_id_dict

    def copy_sg_express(self):
        """
        复制sg_express
        """
        finds = self.db["sg_express"].find({"task_id": self._old_task_id})
        news = []
        old_express_ids = []
        for i in finds:
            i['task_id'] = self._new_task_id
            i['project_sn'] = self._new_project_sn
            old_express_ids.append(str(i.pop('_id')))
            news.append(i)
            if 'params' in i:
                i['params'] = self.params_exchange(i['params'])
        if news:
            result = self.db["sg_express"].insert_many(news)
        else:
            raise Exception('没有任何表达信息，请核对任务结果是否完整')
        self.express_id_dict = dict(zip(old_express_ids, [str(one) for one in result.inserted_ids]))
        self._exchange_dict['express_id'] = self.express_id_dict
        return self.express_id_dict

    def copy_sg_geneset(self):
        """
        复制sg_geneset表
        """
        finds = self.db["sg_geneset"].find({"task_id": self._old_task_id})
        news = []
        old_geneset_ids = []
        for i in finds:
            old_geneset_ids.append(str(i.pop('_id')))
            i['task_id'] = self._new_task_id
            i['project_sn'] = self._new_project_sn
            try:
                i['group_id'] = self.specimen_group_id_dict[str(i['group_id'])]
            except:
                print "该基因集没有分组"
                sys.stdout.flush()
            news.append(i)
        if news:
            result = self.db.sg_geneset.insert_many(news)
        else:
            raise Exception('不存在基因集,请检查数据完整性')
        self.geneset_id_dict = dict(zip(old_geneset_ids, [str(one) for one in result.inserted_ids]))
        self._exchange_dict['geneset_id'] = self.geneset_id_dict
        return self.geneset_id_dict

    def _copy_main_details(self, collection, main_field, change_dict, others_position=[]):
        """
        公共模块，一般用于更新detail表，根据提供的主表id字段名，和主表新旧ID字典，进行查找，再复制替换，others_position用于更新主表ID之外其他需要更新的ID
        params collection: detail表名称
        params main_field: 主表字段名称
        params change_dict: 主表新旧替换字典，一般来源于 copy_collection_with_change 的返回字典
        params others_position: detail表中除了主表还需要更新的字段，
            只能是 specimen_id,group_id
        """
        time_start = datetime.datetime.now()
        coll = self.db[collection]
        for old, new in change_dict.items():
            finds = coll.find({main_field: ObjectId(old)})
            news = []
            for i in finds:
                i.pop('_id')
                i[main_field] = ObjectId(new)
                for position in others_position:
                    i[position] = self.exchange_ObjectId(position, i[position])
                news.append(i)
            if news:
                coll.insert_many(news)
            else:
                print 'WARNING: 主表:{}没有detail表信息，请注意数据合理性,collection:{}'.format(old, collection)
                sys.stdout.flush()
        time_end = datetime.datetime.now()
        run_time = (time_end - time_start).seconds
        print "{}复制运行时间: {}s".format(collection, run_time)
        sys.stdout.flush()

    def copy_main_details(self, collection, main_field, change_dict, others_position=[], join=True):
        greenlet = Greenlet(self._copy_main_details, collection, main_field, change_dict, others_position)
        greenlet.start()
        if join is True:
            greenlet.join()
            return greenlet.value
        self.all_greenlets.append(greenlet)
        return greenlet

    def copy_collection_with_change(self, collection, change_positions=[], update_sg_status=False, targetcoll=None):
        """
        公共模块，一般用于导入主表数据，依靠task_id进行查询，修改change_positions提供的字段，相应修改ID为新的，同时更新params中的数据ID
        params collection: 主表名称
        params change_positions: 需要替换的ID,可用为specimen_id,group_id...
        params update_sg_status: 更新sg_status表
        params targetcoll: 更新到特定集合， 默认与collection参数相同
        """
        coll = self.db[collection]
        if targetcoll:
            targetcoll = self.db[targetcoll]
        else:
            targetcoll = self.db[collection]
        finds = coll.find({'task_id': self._old_task_id})
        news = []
        olds = []
        for i in finds:
            i['task_id'] = self._new_task_id
            if 'project_sn' in i:
                i['project_sn'] = self._new_project_sn
            olds.append(str(i.pop('_id')))
            for position in change_positions:
                if position in i:
                    i[position] = self.exchange_ObjectId(position, i[position])
            if 'params' in i:
                i['params'] = self.params_exchange(i['params'])
            news.append(i)
        if news:
            result = targetcoll.insert_many(news)
            if update_sg_status:
                self.insert_new_status(collection, news, result.inserted_ids)
            return dict(zip(olds, [str(one) for one in result.inserted_ids]))
        else:
            return {}

    def exchange_ObjectId(self, key, thisObjectId):
        """
        用于替换id，key是该ID的字段名，thisObjectId是旧的ID(ObjectId类型)
        """
        if isinstance(thisObjectId, ObjectId):
            return ObjectId(self._exchange_dict[key][str(thisObjectId)])
        else:
            return self._exchange_dict[key][thisObjectId]  # 不是ObjectId时直接返回也是字符串

    def params_exchange(self, params_str):
        """
        专门用于params的数据ID替换
        """
        try:
            params = json.loads(params_str)
        except Exception:
            print("WRANNING：非json格式的params：{}".format(params_str))
            sys.stdout.flush()
            # return params_str
            return "null"
        if not params:
            return "null"
        if 'group_detail' in params:
            for one_group in params['group_detail']:
                params['group_detail'][one_group] = [self.specimen_id_dict[one_sp] for one_sp in params['group_detail'][one_group]]
            params['group_detail'] = group_detail_sort(params['group_detail'])
            if 'second_group_detail' in params:
                if params['second_group_detail']:
                    for one_group in params['second_group_detail']:
                        params['second_group_detail'][one_group] = [self.specimen_id_dict[one_sp] for one_sp in params['second_group_detail'][one_group]]
                    params['second_group_detail'] = group_detail_sort(params['second_group_detail'])
                    if 'second_group_id' in params:
                        params['second_group_id'] = self.group_id_dict[params['second_group_id']]
        if 'group_id' in params:
            params['group_id'] = self.specimen_group_id_dict[params['group_id']]
        if 'express_id' in params:
            params['express_id'] = self.express_id_dict[params['express_id']]
        if 'control_id' in params:
            params['control_id'] = self.specimen_compare_id_dict[params['control_id']]
        if 'geneset_id' in params:
            try:
                geneset_ids = []
                for one_id in params['geneset_id'].split(','):
                    geneset_ids.append(self.geneset_id_dict[one_id])
                params['geneset_id'] = ','.join(geneset_ids)
            except:
                params['geneset_id'] = self.geneset_id_dict
        return json.dumps(params, sort_keys=True, separators=(',', ':'))

    def update_sg_splicing_rmats(self):
            """
            更新sg_splicing_rmats的ref_gtf和params里的splicing_id
            """
            find = self.db["sg_splicing_rmats"].find_one({"task_id": self._new_task_id})  # 同一个task_id里params里的splcing_id皆为第一个主表的_id
            if 'params' in find:
                params_str = find['params']
                try:
                    params = json.loads(params_str)
                    if 'splicing_id' in params:
                        splicing_id = str(find['_id'])
                except Exception:
                    raise Exception("sg_splicing_rmats表的params非json格式")
            finds = self.db["sg_splicing_rmats"].find({"task_id": self._new_task_id})
            for i in finds:
                if self._new_ref_gtf:
                    self.db["sg_splicing_rmats"].update_one({"_id": i["_id"]}, {"$set": {"ref_gtf": self._new_ref_gtf}})
                if 'params' in i:
                    params_str = i['params']
                    try:
                        params = json.loads(params_str)
                        if 'splicing_id' in params:
                            params['splicing_id'] = splicing_id
                    except Exception:
                        raise Exception("sg_splicing_rmats表的params非json格式")
                    new_params = json.dumps(params, sort_keys=True, separators=(',', ':'))
                    self.db["sg_splicing_rmats"].update_one({"_id": i["_id"]}, {"$set": {"params": new_params}})

    def insert_new_status(self, collection, main_docs, ids):
        """
        导入mongo表sg_status数据信息
        """
        coll = self.db["sg_status"]
        news = []
        for index, doc in enumerate(main_docs):
            try:
                submit_location = json.loads(doc['params'])['submit_location']
            except Exception:
                print("WARNING: params参数没有submit_location字段, Doc:{}".format(doc))
                sys.stdout.flush()
                submit_location = None
            status = {
                "status": doc['status'],
                "table_id": ids[index],
                "time": doc['created_ts'],
                "task_id": self._new_task_id,
                "params": doc['params'],
                "table_name": doc['name'],
                "submit_location": submit_location,
                "type_name": collection,
                "is_new": "new",
                "desc": doc['desc'] if 'desc' in doc else None
            }
            news.append(status.copy())
        if news:
            coll.insert_many(news)


if __name__ == '__main__':
    new_bam_path = "/mnt/ilustre/users/sanger/test/bam/"
    new_ref_gtf = "/mnt/ilustre/users/sanger/test/Mus_musculus.GRCm38.87.gff3.gtf"
    # time.sleep(21600)  # 指定休眠时间
    start_time =  time.time()
    num = 1
    for i in range(num):
        task_id = "refrna_demo_mouse_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        copy_task = RefrnaCopyMongo('sanger_21455', task_id, 'refrna_demo', 'refrna_demo', new_bam_path, new_ref_gtf)
        copy_task.run()
    end_time = time.time()
    print "total time: {}".format(end_time - start_time)
