# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
import os
import re
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import types
import gridfs
import json
import unittest
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
from mbio.packages.denovo_rna_v2.nr_stat import nr_stat
from mbio.api.database.denovo_rna_v2.api_base import ApiBase


class DenovoAnnotation(ApiBase):
    def __init__(self, bind_object):
        super(DenovoAnnotation, self).__init__(bind_object)
        self.result_dir = ''
        self.result_file = {}
        self.trans_gene = {}
        self.trans_isgene = {}
        self.version = "v1"
        self.task_id = self.bind_object.sheet.id
        self.kegg_json = Config().SOFTWARE_DIR + "/database/KEGG/br08901.json"
        self.anno_type = 'origin'
        #self._db_name = Config().MONGODB + '_ref_rna'

    def set_result_dir(self, annotation_mudule_dir):
        '''
        根据注释模块结果导入结果路径
        '''
        self.result_dir = annotation_mudule_dir
        self.bind_object.logger.info("导入**** {}".format(annotation_mudule_dir))
        self.result_file['new_stat_path'] = os.path.join(annotation_mudule_dir, "anno_stat/all_annotation_statistics.xls")
        self.result_file['new_venn_path'] = os.path.join(annotation_mudule_dir, "anno_stat/venn")
        blast_path = annotation_mudule_dir + "/anno_stat/blast"
        # for db in ["nr", "swissprot", "string", "kegg"]:
        for db in ["nr", "swissprot"]:
            db_gene_name = db + "gene_blast_path"
            db_trans_name = db + "trans_blast_path"
            self.result_file[db_gene_name] = blast_path + '/gene_' + db + '.xls'
            self.result_file[db_trans_name] = blast_path + '/' + db + '.xls'

        nr_path = annotation_mudule_dir + "/anno_stat/blast_nr_statistics"
        self.result_file['evalue_nr_path'] = nr_path + "/nr_evalue.xls"
        self.result_file['similar_nr_path'] = nr_path + "/nr_similar.xls"
        self.result_file['gene_evalue_nr_path'] = nr_path + "/gene_nr_evalue.xls"
        self.result_file['gene_similar_nr_path'] = nr_path + "/gene_nr_similar.xls"
        self.result_file['gene_species_nr_path'] = nr_path + "/gene_nr_species.xls"

        swissprot_path = annotation_mudule_dir + "/anno_stat/blast_swissprot_statistics"
        self.result_file['evalue_swissprot_path'] = swissprot_path + "/swissprot_evalue.xls"
        self.result_file['similar_swissprot_path'] = swissprot_path + "/swissprot_similar.xls"
        self.result_file['gene_evalue_swissprot_path'] = swissprot_path + "/gene_swissprot_evalue.xls"
        self.result_file['gene_similar_swissprot_path'] = swissprot_path + "/gene_swissprot_similar.xls"

        self.result_file['gene_pfam_path'] = annotation_mudule_dir + "/anno_stat/pfam_stat/gene_pfam_domain"
        self.result_file['pfam_path'] = annotation_mudule_dir + "/blast_xml/pfam_domain"
        self.result_file['query_path'] = annotation_mudule_dir + "/anno_stat/trans_anno_detail.xls"
        self.result_file['gene_query_path'] = annotation_mudule_dir + "/anno_stat/gene_anno_detail.xls"

        self.result_file['n_sum_path'] = annotation_mudule_dir + "/cog/cog_summary.xls"
        self.result_file['n_gene_sum_path'] = annotation_mudule_dir + "/anno_stat/cog_stat/gene_cog_summary.xls"

        go_path = annotation_mudule_dir + "/go"
        gene_go_path = annotation_mudule_dir + "/anno_stat/go_stat"

        self.result_file['stat_level2'] = go_path + "/go12level_statistics.xls"
        self.result_file['stat_level3'] = go_path + "/go123level_statistics.xls"
        self.result_file['stat_level4'] = go_path + "/go1234level_statistics.xls"
        self.result_file['gene_stat_level2'] = gene_go_path + "/gene_go12level_statistics.xls"
        self.result_file['gene_stat_level3'] = gene_go_path + "/gene_go123level_statistics.xls"
        self.result_file['gene_stat_level4'] = gene_go_path + "/gene_go1234level_statistics.xls"
        self.result_file['gos_path'] = go_path + "/query_gos.list"
        self.result_file['gene_gos_path'] = gene_go_path + "/gene_gos.list"

        self.result_file['n_stat_level2'] = annotation_mudule_dir + "/go/go12level_statistics.xls"
        self.result_file['n_stat_level3'] = annotation_mudule_dir + "/go/go123level_statistics.xls"
        self.result_file['n_stat_level4'] = annotation_mudule_dir + "/go/go1234level_statistics.xls"

        self.result_file['n_gene_stat_level2'] = annotation_mudule_dir + "/anno_stat/go_stat/gene_go12level_statistics.xls"
        self.result_file['n_gene_stat_level3'] = annotation_mudule_dir + "/anno_stat/go_stat/gene_go123level_statistics.xls"
        self.result_file['n_gene_stat_level4'] = annotation_mudule_dir + "/anno_stat/go_stat/gene_go1234level_statistics.xls"

        kegg_path = annotation_mudule_dir + "/kegg"
        gene_kegg_path = annotation_mudule_dir + "/anno_stat/kegg_stat"

        self.result_file['layer_path'] = kegg_path + "/kegg_layer.xls"
        self.result_file['pathway_path'] = kegg_path + "/pathway_table.xls"
        self.result_file['png_path'] = kegg_path + "/pathways"
        self.result_file['table_path'] = kegg_path + "/kegg_table.xls"
        self.result_file['gene_layer_path'] = gene_kegg_path + "/gene_kegg_layer.xls"
        self.result_file['gene_pathway_path'] = gene_kegg_path + "/gene_pathway_table.xls"
        self.result_file['gene_png_path'] = gene_kegg_path + "/gene_pathway"
        self.result_file['gene_table_path'] = gene_kegg_path + "/gene_kegg_table.xls"


        for key, value in self.result_file.items():
            if os.path.exists(value):
                pass
            else:
                self.bind_object.set_error('%s对应的结果文件%s 不存在，请检查', variables=(key, value), code="52000501")
        self.bind_object.logger.info("数据路径正确，文件完整 {}")


    def set_result_dir_v2(self, annotation_mudule_dir):
        '''
        根据注释模块结果导入结果路径
        '''
        self.result_dir = annotation_mudule_dir
        self.bind_object.logger.info("导入**** {}".format(annotation_mudule_dir))
        self.result_file['new_stat_path'] = os.path.join(annotation_mudule_dir, "all_stat.xls")
        self.result_file['new_venn_path'] = os.path.join(annotation_mudule_dir)

        for db in ["nr", "swissprot"]:
            db_gene_name = db + "gene_blast_path"
            db_trans_name = db + "trans_blast_path"
            self.result_file[db_gene_name] = annotation_mudule_dir + "/{}/{}_blast_gene.xls".format(db, db)
            self.result_file[db_trans_name] = annotation_mudule_dir + "/{}/{}_blast_tran.xls".format(db, db)

        nr_path = annotation_mudule_dir + "/nr"
        self.result_file['evalue_nr_path'] = nr_path + "/trans_nr_evalue.xls"
        self.result_file['similar_nr_path'] = nr_path + "/trans_nr_similar.xls"
        self.result_file['species_nr_path'] = nr_path + "/tran_nr_species_stat.xls"
        self.result_file['gene_evalue_nr_path'] = nr_path + "/gene_nr_evalue.xls"
        self.result_file['gene_similar_nr_path'] = nr_path + "/gene_nr_similar.xls"
        self.result_file['gene_species_nr_path'] = nr_path + "/gene_nr_species_stat.xls"
        self.result_file['gene_species_nr_detail'] = nr_path + "/nr.G.tsvspecies.txt"
        self.result_file['species_nr_detail'] = nr_path + "/nr.T.tsvspecies.txt"

        swissprot_path = annotation_mudule_dir + "/swissprot"
        self.result_file['evalue_swissprot_path'] = swissprot_path + "/trans_swissprot_evalue.xls"
        self.result_file['similar_swissprot_path'] = swissprot_path + "/trans_swissprot_similar.xls"
        self.result_file['gene_evalue_swissprot_path'] = swissprot_path + "/gene_swissprot_evalue.xls"
        self.result_file['gene_similar_swissprot_path'] = swissprot_path + "/gene_swissprot_similar.xls"

        self.result_file['gene_pfam_path'] = annotation_mudule_dir + "/pfam/pfam_domain_gene.xls"
        self.result_file['pfam_path'] = annotation_mudule_dir + "/pfam/pfam_domain_tran.xls"
        self.result_file['query_path'] = annotation_mudule_dir + "/all_annot.xls"
        self.result_file['gene_query_path'] = annotation_mudule_dir + "/all_annot.xls"

        self.result_file['n_sum_path'] = annotation_mudule_dir + "/cog/summary.T.tsv"
        self.result_file['n_gene_sum_path'] = annotation_mudule_dir + "/cog/summary.G.tsv"

        go_path = annotation_mudule_dir + "/go"
        gene_go_path = annotation_mudule_dir + "/anno_stat/go_stat"

        self.result_file['stat_level2'] = go_path + "/go_lev2_tran.stat.xls"
        self.result_file['stat_level3'] = go_path + "/go_lev3_tran.stat.xls"
        self.result_file['stat_level4'] = go_path + "/go_lev4_tran.stat.xls"
        self.result_file['gene_stat_level2'] = go_path + "/go_lev2_gene.stat.xls"
        self.result_file['gene_stat_level3'] = go_path + "/go_lev3_gene.stat.xls"
        self.result_file['gene_stat_level4'] = go_path + "/go_lev4_gene.stat.xls"
        self.result_file['gos_path'] = go_path + "/go_list_tran.xls"
        self.result_file['gene_gos_path'] = go_path + "/go_list_gene.xls"

        '''
        self.result_file['n_stat_level2'] = annotation_mudule_dir + "/go/go12level_statistics.xls"
        self.result_file['n_stat_level3'] = annotation_mudule_dir + "/go/go123level_statistics.xls"
        self.result_file['n_stat_level4'] = annotation_mudule_dir + "/go/go1234level_statistics.xls"

        self.result_file['n_gene_stat_level2'] = annotation_mudule_dir + "/anno_stat/go_stat/gene_go12level_statistics.xls"
        self.result_file['n_gene_stat_level3'] = annotation_mudule_dir + "/anno_stat/go_stat/gene_go123level_statistics.xls"
        self.result_file['n_gene_stat_level4'] = annotation_mudule_dir + "/anno_stat/go_stat/gene_go1234level_statistics.xls"
        '''
        kegg_path = annotation_mudule_dir + "/kegg"

        self.result_file['layer_path'] = kegg_path + "/kegg_layer_tran.xls"
        self.result_file['pathway_path'] = kegg_path + "/kegg_pathway_tran.xls"
        self.result_file['png_path'] = kegg_path + "/kegg_pathway_tran_dir"
        self.result_file['table_path'] = kegg_path + "/kegg_gene_tran.xls"
        self.result_file['gene_layer_path'] = kegg_path + "/kegg_layer_gene.xls"
        self.result_file['gene_pathway_path'] = kegg_path + "/kegg_pathway_gene.xls"
        self.result_file['gene_png_path'] = kegg_path + "/kegg_pathway_gene_dir"
        self.result_file['gene_table_path'] = kegg_path + "/kegg_gene_gene.xls"


        for key, value in self.result_file.items():
            if os.path.exists(value):
                pass
            else:
                self.bind_object.set_error('%s对应的结果文件%s 不存在，请检查', variables=(key, value), code="52000502")
        self.bind_object.logger.info("数据路径正确，文件完整 {}")



    def check_id(self, object_id):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('assemble_id必须为ObjectId对象或其对应的字符串！', code="52000503")
        return object_id



    def get_trans2gene(self, trans2gene):
        trans2gene = str(trans2gene)
        if os.path.exists(trans2gene):
            pass
        else:
            self.bind_object.set_error('转录本基因对应的结果文件%s不存在，请检查', variables=(trans2gene), code="52000504")
        self.bind_object.logger.info("读入基因转录本对应关系文件 {}".format(trans2gene))
        with open(trans2gene, 'rb') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split("\t")
                self.trans_gene[line[0]] = line[1]
                if line[2] == "yes":
                    self.trans_isgene[line[0]] = True
                else:
                    self.trans_isgene[line[0]] = False
        self.bind_object.logger.info("读入基因转录本对应关系文件结束")


    def run_nr_stat(self, trans_taxonfile, gene_taxonfile, trans_exp, gene_exp):
        samples = nr_stat().nr_stat_species2(
            tran_taxonfile = trans_taxonfile,
            gene_list = "",
            exp = trans_exp,
            outpath = trans_taxonfile + '.exp.xls'
        )

        samples = nr_stat().nr_stat_species2(
            tran_taxonfile = gene_taxonfile,
            gene_list = "",
            exp = gene_exp,
            outpath = gene_taxonfile + ".exp.xls"
        )
        return samples


    def run(self, result_dir, trans2gene, params_dict, taxon='Animals', version="v1", exp_level="T",  gene_exp = None, trans_exp = None):
        """
        new_anno_path: 新序列注释的结果文件夹
        pfam_path:转录本的pfam_domain
        merge_tran_output: 转录本的merge_annot tool输出结果路径
        merge_gene_output: 基因的merge_annot tool"
        """
        self.bind_object.logger.info("开始到表情数据路径为 {}".format(result_dir))
        if version == "v2":
            self.set_result_dir_v2(result_dir)
        else:
            self.set_result_dir(result_dir)
        self.get_trans2gene(trans2gene)

        self.version = version

        task_id = self.task_id
        params_dict.update({"exclude_taxon": ""})

        params = json.dumps(params_dict, sort_keys=True, separators=(',', ':'))
        self.remove_table_by_main_record(main_table='sg_annotation_stat', task_id=task_id, type=self.anno_type, detail_table=['sg_annotation_stat_detail'], detail_table_key='stat_id')
        stat_id = self.add_annotation_stat(name=None, params=params, seq_type="new", database="nr,swissprot,pfam,cog,go,kegg", taxon=taxon, result_dir=self.result_dir)
        self.add_annotation_stat_detail(stat_id=stat_id, stat_path=self.result_file['new_stat_path'], venn_path=self.result_file['new_venn_path'], exp_level=exp_level, gene_exp = gene_exp, trans_exp = trans_exp)
        self.update_db_record('sg_annotation_stat', stat_id, status="end", main_id=stat_id)

        # self.update_db_record('sg_assembly', obj_id, status="end", main_id=obj_id)
        # blast_id = self.add_annotation_blast(name=None, params=params, stat_id=stat_id, result_dir=self.result_dir)
        # for db in ["nr", "swissprot"]:
            #self.add_annotation_blast_detail(blast_id=blast_id, seq_type="new", anno_type="T", database=db, blast_path=self.result_file[db + "trans_blast_path"])
            #self.add_annotation_blast_detail(blast_id=blast_id, seq_type="new", anno_type="G", database=db, blast_path=self.result_file[db + "gene_blast_path"])

        params_select_nr = dict([(k,params_dict.get(k,None)) for k in ('nr_evalue', 'nr_similarity', 'nr_identity')])
        params_select_nr = json.dumps(params_select_nr, sort_keys=True, separators=(',', ':'))
        self.remove_table_by_main_record(main_table='sg_annotation_nr', task_id=task_id, type=self.anno_type, detail_table=['sg_annotation_nr_detail', 'sg_annotation_nr_pie'], detail_table_key='nr_id')

        samples = self.run_nr_stat(self.result_file["species_nr_detail"], self.result_file["gene_species_nr_detail"], trans_exp, gene_exp)
        nr_id = self.add_annotation_nr(name=None, params=params_select_nr, stat_id=stat_id, result_dir=self.result_dir)
        self.add_annotation_blast_nr_detail(blast_id=nr_id, seq_type="new", anno_type="T", database='nr', blast_path=self.result_file["nr" + "trans_blast_path"])
        self.add_annotation_nr_pie(nr_id=nr_id, species_path=self.result_file['species_nr_path'] ,  evalue_path=self.result_file['evalue_nr_path'], similar_path=self.result_file['similar_nr_path'], seq_type="new", anno_type="T")
        self.add_annotation_nr_pie(nr_id=nr_id, species_path=self.result_file['gene_species_nr_path'],  evalue_path=self.result_file['gene_evalue_nr_path'], similar_path=self.result_file['gene_similar_nr_path'], seq_type="new", anno_type="G")
        self.add_annotation_nr_species_detail(nr_id=nr_id, gene_species_count=self.result_file["gene_species_nr_detail"] + ".exp.xls", trans_species_count=self.result_file["species_nr_detail"] + ".exp.xls", seq_type="new", anno_type="G")
        self.update_db_record('sg_annotation_nr', nr_id, status="end", samples=samples,  main_id=nr_id)

        self.remove_table_by_main_record(main_table='sg_annotation_swissprot', task_id=task_id, type=self.anno_type,  detail_table=['sg_annotation_swissprot_detail', 'sg_annotation_swissprot_pie'], detail_table_key='swissprot_id')
        params_select_swissprot = dict([(k,params_dict.get(k,None)) for k in ('swissprot_evalue', 'swissprot_similarity', 'swissprot_identity')])
        params_select_swissprot = json.dumps(params_select_swissprot, sort_keys=True, separators=(',', ':'))
        swissprot_id = self.add_annotation_swissprot(name=None, params=params_select_swissprot, stat_id=stat_id, result_dir=self.result_dir)
        self.add_annotation_blast_swissprot_detail(blast_id=swissprot_id, seq_type="new", anno_type="T", database='swissprot', blast_path=self.result_file["swissprot" + "trans_blast_path"])
        self.add_annotation_swissprot_pie(swissprot_id=swissprot_id, evalue_path=self.result_file['evalue_swissprot_path'], similar_path=self.result_file['similar_swissprot_path'], seq_type="new", anno_type="T")
        self.add_annotation_swissprot_pie(swissprot_id=swissprot_id, evalue_path=self.result_file['gene_evalue_swissprot_path'], similar_path=self.result_file['gene_similar_swissprot_path'], seq_type="new", anno_type="G")
        self.update_db_record('sg_annotation_swissprot', swissprot_id, status="end", main_id=swissprot_id)

        self.remove_table_by_main_record(main_table='sg_annotation_pfam', task_id=task_id, type=self.anno_type, detail_table=['sg_annotation_pfam_detail', 'sg_annotation_pfam_bar'], detail_table_key='pfam_id')
        params_select_pfam = dict([('pfam_evalue',params_dict['pfam_evalue'])])
        params_select_pfam = json.dumps(params_select_pfam, sort_keys=True, separators=(',', ':'))
        pfam_id = self.add_annotation_pfam(name=None, params=params_select_pfam, stat_id=stat_id, result_dir=self.result_dir)
        self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=self.result_file['pfam_path'], seq_type="new", anno_type="T")
        #self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=self.result_file['gene_pfam_path'], seq_type="new", anno_type="G")
        self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=self.result_file['pfam_path'], seq_type="new", anno_type="T")
        self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=self.result_file['gene_pfam_path'], seq_type="new", anno_type="G")
        self.update_db_record('sg_annotation_pfam', pfam_id, status="end", main_id=pfam_id)
        self.remove_table_by_main_record(main_table='sg_annotation_query', task_id=task_id, type=self.anno_type, detail_table=['sg_annotation_query_detail'], detail_table_key='query_id')
        query_id = self.add_annotation_query(name=None, params=params, stat_id=stat_id, result_dir=self.result_dir)
        self.add_annotation_query_denovo_detail(query_id=query_id, query_path=self.result_file['query_path'], anno_type="T")
        self.update_db_record('sg_annotation_query',query_id, status="end", main_id=query_id)
        #self.add_annotation_gene_query_denovo_detail(query_id=query_id, query_path=self.result_file['gene_query_path'], anno_type="G")
        self.remove_table_by_main_record(main_table='sg_annotation_cog', task_id=task_id, type=self.anno_type, detail_table=['sg_annotation_cog_detail'], detail_table_key='cog_id')
        params_select_cog = dict([(k,params_dict.get(k,None)) for k in ('cog_evalue', 'cog_similarity', 'cog_identity')])
        params_select_cog = json.dumps(params_select_cog, sort_keys=True, separators=(',', ':'))

        cog_id = self.add_annotation_cog(name=None, params=params_select_cog, result_dir=self.result_dir)
        self.add_annotation_cog_detail(cog_id=cog_id, cog_path=self.result_file['n_sum_path'], seq_type="new", anno_type="T")
        self.add_annotation_cog_detail(cog_id=cog_id, cog_path=self.result_file['n_gene_sum_path'], seq_type="new", anno_type="G")
        self.update_db_record('sg_annotation_cog', cog_id, status="end", main_id=cog_id)
        self.remove_table_by_main_record(main_table='sg_annotation_go', task_id=task_id, type=self.anno_type, detail_table=['sg_annotation_go_detail', 'sg_annotation_go_graph','sg_annotation_go_level', 'sg_annotation_go_list'], detail_table_key='go_id')
        go_id = self.add_annotation_go(name=None, params=params_select_nr, result_dir=self.result_dir)
        seq_type = "new"
        self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, level_path=self.result_file['stat_level2'])
        self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, level_path=self.result_file['gene_stat_level2'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, go_path=self.result_file['stat_level2'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=3, go_path=self.result_file['stat_level3'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=4, go_path=self.result_file['stat_level4'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, go_path=self.result_file['gene_stat_level2'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, go_path=self.result_file['gene_stat_level3'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, go_path=self.result_file['gene_stat_level4'])
        self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, go_path=self.result_file['stat_level2'])
        self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=3, go_path=self.result_file['stat_level3'])
        self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=4, go_path=self.result_file['stat_level4'])
        self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, go_path=self.result_file['gene_stat_level2'])
        self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, go_path=self.result_file['gene_stat_level3'])
        self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, go_path=self.result_file['gene_stat_level4'])
        self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="T", gos_path=self.result_file['gos_path'])
        self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="G", gos_path=self.result_file['gene_gos_path'])
        self.update_db_record('sg_annotation_go', go_id, status="end", main_id=go_id)

        self.remove_table_by_main_record(main_table='sg_annotation_kegg', task_id=task_id, type=self.anno_type, detail_table=['sg_annotation_kegg_categories', 'sg_annotation_kegg_level','sg_annotation_kegg_table'], detail_table_key='kegg_id')
        params_select_kegg = dict([(k,params_dict.get(k,None)) for k in ('kegg_evalue', 'kegg_similarity', 'kegg_identity')])
        params_select_kegg = json.dumps(params_select_kegg, sort_keys=True, separators=(',', ':'))
        kegg_id = self.add_annotation_kegg(name=None, params=params_select_kegg, result_dir=self.result_dir)
        self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", categories_path=self.result_file['layer_path'])
        self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", categories_path=self.result_file['gene_layer_path'])
        self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", level_path=self.result_file['pathway_path'], png_dir=self.result_file['png_path'])
        self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", level_path=self.result_file['gene_pathway_path'], png_dir=self.result_file['gene_png_path'])
        self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", table_path=self.result_file['table_path'])
        self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", table_path=self.result_file['gene_table_path'])
        self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type="new", anno_type="G", level_path=self.result_file['gene_pathway_path'], png_dir=self.result_file['gene_png_path'])
        self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type="new", anno_type="T", level_path=self.result_file['pathway_path'], png_dir=self.result_file['gene_png_path'])
        self.update_db_record('sg_annotation_kegg', kegg_id, status="end", main_id=kegg_id)


    def run_webroot(self, result_dir, trans2gene, params_dict, task_id, stat_id, last_id, taxonomy,  version="v1", exp_level="T", gene_exp = None, trans_exp = None):
        """
        用于注释重运行导表
        result_dir: 新序列注释的结果文件夹
        trans2gene: 转录本基因对应关系文件
        params_dict: 参数
        stat_id: 统计结果主表ID
        last_id: 上次重运行ID
        taxonomy: 物种分类
        """
        self.bind_object.logger.info("开始导表webroot数据路径为 {}".format(result_dir))
        if version == "v2":
            self.set_result_dir_v2(result_dir)
        else:
            self.set_result_dir(result_dir)

        self.version = version
        # self.set_result_dir(result_dir)
        self.get_trans2gene(trans2gene)
        self.bind_object.logger.info("开始导表task_id为 {}".format(self.task_id))

        stat_id = ObjectId(stat_id)
        if last_id:
            last_id = ObjectId(last_id)
        else:
            pass
        self.task_id = task_id
        #task_id = self.task_id
        self.bind_object.logger.info("开始导表task_id为 {}".format(task_id))
        # task_id = "denovo_rna_v2"
        params = json.dumps(params_dict, sort_keys=True, separators=(',', ':'))

        self.add_annotation_stat_detail(stat_id=stat_id, stat_path=self.result_file['new_stat_path'], venn_path=self.result_file['new_venn_path'], exp_level="T", gene_exp=gene_exp, trans_exp=trans_exp)
        self.update_db_record('sg_annotation_stat', stat_id, status="end", main_id=stat_id)
        if last_id:
            self.bind_object.logger.info("删除表格为 {}".format(last_id))
            self.remove_table_by_main_record(main_table='sg_annotation_stat', _id=last_id, detail_table=['sg_annotation_stat_detail'], detail_table_key='stat_id')
            self.bind_object.logger.info("删除表格成功 {}".format(last_id))

        #for db in ["nr", "swissprot"]:
            #self.add_annotation_blast_detail(blast_id=blast_id, seq_type="new", anno_type="T", database=db, blast_path=self.result_file[db + "trans_blast_path"])
            #self.add_annotation_blast_detail(blast_id=blast_id, seq_type="new", anno_type="G", database=db, blast_path=self.result_file[db + "gene_blast_path"])

        params_select_nr = dict([(k,params_dict.get(k,None)) for k in ('nr_evalue', 'nr_similarity', 'nr_identity')])
        params_select_nr = json.dumps(params_select_nr, sort_keys=True, separators=(',', ':'))

        nr_old_id = self.get_table_by_main_record(main_table='sg_annotation_nr', task_id=task_id, type=self.anno_type)
        self.bind_object.logger.info("查找表格{}".format(nr_old_id))

        samples = self.run_nr_stat(self.result_file["species_nr_detail"], self.result_file["gene_species_nr_detail"], trans_exp, gene_exp)
        nr_id = self.add_annotation_nr(name=None, params=params_select_nr, stat_id=stat_id, result_dir=self.result_dir)

        self.add_annotation_blast_nr_detail(blast_id=nr_id, seq_type="new", anno_type="T", database='nr', blast_path=self.result_file["nr" + "trans_blast_path"])
        self.bind_object.logger.info("插入新表 {}".format(nr_id))
        self.add_annotation_nr_pie(nr_id=nr_id, species_path=self.result_file['species_nr_path'] ,  evalue_path=self.result_file['evalue_nr_path'], similar_path=self.result_file['similar_nr_path'], seq_type="new", anno_type="T")
        self.add_annotation_nr_pie(nr_id=nr_id, species_path=self.result_file['gene_species_nr_path'],  evalue_path=self.result_file['gene_evalue_nr_path'], similar_path=self.result_file['gene_similar_nr_path'], seq_type="new", anno_type="G")
        self.add_annotation_nr_species_detail(nr_id=nr_id, gene_species_count=self.result_file["gene_species_nr_detail"] + ".exp.xls", trans_species_count=self.result_file["species_nr_detail"] + ".exp.xls", seq_type="new", anno_type="G")
        self.update_db_record('sg_annotation_nr', nr_id, status="end", samples=samples, main_id=nr_id)

        if nr_old_id:
            self.bind_object.logger.info( "删除旧表 {}".format(nr_old_id['_id']))
            self.remove_table_by_main_record(main_table='sg_annotation_nr', _id=nr_old_id['_id'], detail_table=['sg_annotation_nr_detail', 'sg_annotation_nr_pie'], detail_table_key='nr_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        swissprot_old_id = self.get_table_by_main_record(main_table='sg_annotation_swissprot', task_id=task_id, type=self.anno_type)
        self.bind_object.logger.info("查找表格成{}".format(swissprot_old_id))
        params_select_swissprot = dict([(k,params_dict.get(k,None)) for k in ('swissprot_evalue', 'swissprot_similarity', 'swissprot_identity')])
        params_select_swissprot = json.dumps(params_select_swissprot, sort_keys=True, separators=(',', ':'))
        swissprot_id = self.add_annotation_swissprot(name=None, params=params_select_swissprot, stat_id=stat_id, result_dir=self.result_dir)
        self.add_annotation_blast_swissprot_detail(blast_id=swissprot_id, seq_type="new", anno_type="T", database='swissprot', blast_path=self.result_file["swissprot" + "trans_blast_path"])
        self.add_annotation_swissprot_pie(swissprot_id=swissprot_id, evalue_path=self.result_file['evalue_swissprot_path'], similar_path=self.result_file['similar_swissprot_path'], seq_type="new", anno_type="T")
        self.add_annotation_swissprot_pie(swissprot_id=swissprot_id, evalue_path=self.result_file['gene_evalue_swissprot_path'], similar_path=self.result_file['gene_similar_swissprot_path'], seq_type="new", anno_type="G")
        self.update_db_record('sg_annotation_swissprot', swissprot_id, status="end", main_id=swissprot_id)
        if swissprot_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_swissprot', _id=swissprot_old_id['_id'],  detail_table=['sg_annotation_swissprot_detail', 'sg_annotation_swissprot_pie'], detail_table_key='swissprot_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")
        pfam_old_id = self.get_table_by_main_record(main_table='sg_annotation_pfam', task_id=task_id, type=self.anno_type)
        params_select_pfam = dict([('pfam_evalue',params_dict['pfam_evalue'])])
        params_select_pfam = json.dumps(params_select_pfam, sort_keys=True, separators=(',', ':'))
        pfam_id = self.add_annotation_pfam(name=None, params=params_select_pfam, stat_id=stat_id, result_dir=self.result_dir)
        self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=self.result_file['pfam_path'], seq_type="new", anno_type="T")

        #self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=self.result_file['gene_pfam_path'], seq_type="new", anno_type="G")
        self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=self.result_file['pfam_path'], seq_type="new", anno_type="T")
        self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=self.result_file['gene_pfam_path'], seq_type="new", anno_type="G")
        self.update_db_record('sg_annotation_pfam', pfam_id, status="end", main_id=pfam_id)
        if pfam_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_pfam', _id=pfam_old_id['_id'], detail_table=['sg_annotation_pfam_detail', 'sg_annotation_pfam_bar'], detail_table_key='pfam_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        query_old_id = self.get_table_by_main_record(main_table='sg_annotation_query', task_id=task_id, type=self.anno_type)
        self.remove_table_by_main_record(main_table='sg_annotation_query', task_id=task_id, type=self.anno_type, detail_table=['sg_annotation_query_detail'], detail_table_key='query_id')
        query_id = self.add_annotation_query(name=None, params=params, stat_id=stat_id, result_dir=self.result_dir)
        self.add_annotation_query_denovo_detail(query_id=query_id, query_path=self.result_file['query_path'], anno_type="T")
        self.update_db_record('sg_annotation_query',query_id, status="end", main_id=query_id)
        if query_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_query', _id=query_old_id['_id'], detail_table=['sg_annotation_query_detail'], detail_table_key='query_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        cog_old_id = self.get_table_by_main_record(main_table='sg_annotation_cog', task_id=task_id, type=self.anno_type)
        params_select_cog = dict([(k,params_dict.get(k,None)) for k in ('cog_evalue', 'cog_similarity', 'cog_identity')])
        params_select_cog = json.dumps(params_select_cog, sort_keys=True, separators=(',', ':'))
        cog_id = self.add_annotation_cog(name=None, params=params_select_cog, result_dir=self.result_dir)
        self.add_annotation_cog_detail(cog_id=cog_id, cog_path=self.result_file['n_sum_path'], seq_type="new", anno_type="T")
        self.add_annotation_cog_detail(cog_id=cog_id, cog_path=self.result_file['n_gene_sum_path'], seq_type="new", anno_type="G")
        self.update_db_record('sg_annotation_cog', cog_id, status="end", main_id=cog_id)
        if cog_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_cog', _id=cog_old_id['_id'], detail_table=['sg_annotation_cog_detail'], detail_table_key='cog_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        go_old_id = self.get_table_by_main_record(main_table='sg_annotation_go', task_id=task_id, type=self.anno_type)

        go_id = self.add_annotation_go(name=None, params=params_select_nr, result_dir=self.result_dir)
        seq_type = "new"
        self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, level_path=self.result_file['stat_level2'])
        self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, level_path=self.result_file['gene_stat_level2'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, go_path=self.result_file['stat_level2'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=3, go_path=self.result_file['stat_level3'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=4, go_path=self.result_file['stat_level4'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, go_path=self.result_file['gene_stat_level2'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, go_path=self.result_file['gene_stat_level3'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, go_path=self.result_file['gene_stat_level4'])
        self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, go_path=self.result_file['stat_level2'])
        self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=3, go_path=self.result_file['stat_level3'])
        self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=4, go_path=self.result_file['stat_level4'])
        self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, go_path=self.result_file['gene_stat_level2'])
        self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, go_path=self.result_file['gene_stat_level3'])
        self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, go_path=self.result_file['gene_stat_level4'])
        self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="T", gos_path=self.result_file['gos_path'])
        self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="G", gos_path=self.result_file['gene_gos_path'])
        self.update_db_record('sg_annotation_go', go_id, status="end", main_id=go_id)
        if go_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_go', _id=go_old_id['_id'], detail_table=['sg_annotation_go_detail', 'sg_annotation_go_graph','sg_annotation_go_level', 'sg_annotation_go_list'], detail_table_key='go_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        kegg_old_id = self.get_table_by_main_record(main_table='sg_annotation_kegg', task_id=task_id, type=self.anno_type)

        params_select_kegg = dict([(k,params_dict.get(k,None)) for k in ('kegg_evalue', 'kegg_similarity', 'kegg_identity')])
        params_select_kegg = json.dumps(params_select_kegg, sort_keys=True, separators=(',', ':'))
        kegg_id = self.add_annotation_kegg(name=None, params=params_select_kegg, result_dir=self.result_dir)
        self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", categories_path=self.result_file['layer_path'])
        self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", categories_path=self.result_file['gene_layer_path'])
        self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", level_path=self.result_file['pathway_path'], png_dir=self.result_file['png_path'])
        self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", level_path=self.result_file['gene_pathway_path'], png_dir=self.result_file['gene_png_path'])
        self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", table_path=self.result_file['table_path'])
        self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", table_path=self.result_file['gene_table_path'])
        self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type="new", anno_type="G", level_path=self.result_file['gene_pathway_path'], png_dir=self.result_file['gene_png_path'])
        self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type="new", anno_type="T", level_path=self.result_file['pathway_path'], png_dir=self.result_file['gene_png_path'])
        self.update_db_record('sg_annotation_kegg', kegg_id, status="end", main_id=kegg_id)
        if kegg_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_kegg', _id=kegg_old_id['_id'], detail_table=['sg_annotation_kegg_categories', 'sg_annotation_kegg_level','sg_annotation_kegg_table'], detail_table_key='kegg_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

    def add_annotation(self, name=None, params=None, ref_anno_path=None, new_anno_path=None, pfam_path=None, merge_tran_output=None, merge_gene_output=None):
        """
        ref_anno_path: 已知序列注释的结果文件夹
        new_anno_path: 新序列注释的结果文件夹
        pfam_path:转录本的pfam_domain
        merge_tran_output: 转录本的merge_annot tool输出结果路径
        merge_gene_output: 基因的merge_annot tool输出结果路径
        """
        new_stat_path = new_anno_path + "/anno_stat/all_annotation_statistics.xls"
        new_venn_path = new_anno_path + "/anno_stat/venn"
        stat_id = self.add_annotation_stat(name=None, params=params, seq_type="new", database="nr,swissprot,pfam,cog,go,kegg")
        self.add_annotation_stat_detail(stat_id=stat_id, stat_path=new_stat_path, venn_path=new_venn_path)
        blast_id = self.add_annotation_blast(name=None, params=params, stat_id=stat_id)
        blast_path = new_anno_path + "/anno_stat/blast"
        if os.path.exists(blast_path):
            for db in ["nr", "swissprot"]:
                trans_blast_path = blast_path + "/" + db + '.xls'
                gene_blast_path = blast_path + '/gene_' + db + '.xls'
                self.add_annotation_blast_detail(blast_id=blast_id, seq_type="new", anno_type="transcript", database=db, blast_path=trans_blast_path)
                self.add_annotation_blast_detail(blast_id=blast_id, seq_type="new", anno_type="gene", database=db, blast_path=gene_blast_path)
        else:
            self.bind_object.set_error("没有blast的结果文件夹", code="52000505")
        nr_id = self.add_annotation_nr(name=None, params=params, stat_id=stat_id)
        nr_path = new_anno_path + "/anno_stat/blast_nr_statistics"
        if os.path.exists(nr_path):
            evalue_path = nr_path + "/nr_evalue.xls"
            similar_path = nr_path + "/nr_similar.xls"
            gene_evalue_path = nr_path + "/gene_nr_evalue.xls"
            gene_similar_path = nr_path + "/gene_nr_similar.xls"
            self.add_annotation_nr_pie(nr_id=nr_id, evalue_path=evalue_path, similar_path=similar_path, seq_type="new", anno_type="transcript")
            self.add_annotation_nr_pie(nr_id=nr_id, evalue_path=gene_evalue_path, similar_path=gene_similar_path,  seq_type="new", anno_type="gene")
        else:
            self.bind_object.set_error("新序列NR注释结果文件不存在", code="52000506")
        swissprot_id = self.add_annotation_swissprot(name=None, params=params, stat_id=stat_id)
        swissprot_path = new_anno_path + "/anno_stat/blast_swissprot_statistics"
        if os.path.exists(swissprot_path):
            evalue_path = swissprot_path + "/swissprot_evalue.xls"
            similar_path = swissprot_path + "/swissprot_similar.xls"
            gene_evalue_path = swissprot_path + "/gene_swissprot_evalue.xls"
            gene_similar_path = swissprot_path + "/gene_swissprot_similar.xls"
            self.add_annotation_swissprot_pie(swissprot_id=swissprot_id, evalue_path=evalue_path, similar_path=similar_path, seq_type="new", anno_type="transcript")
            self.add_annotation_swissprot_pie(swissprot_id=swissprot_id, evalue_path=gene_evalue_path, similar_path=gene_similar_path, seq_type="new", anno_type="gene")
        else:
            self.bind_object.set_error("新序列Swiss-Prot注释结果文件不存在", code="52000507")
        pfam_id = self.add_annotation_pfam(name=None, params=params, stat_id=stat_id)
        gene_pfam_path = new_anno_path + "/anno_stat/pfam_stat/gene_pfam_domain"
        if os.path.exists(pfam_path) and os.path.exists(gene_pfam_path):
            self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=pfam_path, seq_type="new", anno_type="transcript")
            self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=gene_pfam_path, seq_type="new", anno_type="gene")
            self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=pfam_path, seq_type="new", anno_type="transcript")
            self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=gene_pfam_path, seq_type="new", anno_type="gene")
        else:
            self.bind_object.set_error("pfam注释结果文件不存在", code="52000508")
        ref_stat_path = ref_anno_path + "/anno_stat/all_annotation_statistics.xls"
        ref_venn_path = ref_anno_path + "/anno_stat/venn"
        if os.path.exists(ref_stat_path) and os.path.exists(ref_venn_path):
            stat_id = self.add_annotation_stat(name=None, params=params, seq_type="ref" , database="cog,go,kegg")
            self.add_annotation_stat_detail(stat_id=stat_id, stat_path=ref_stat_path, venn_path=ref_venn_path)
        else:
            self.bind_object.set_error("已知序列注释统计文件和venn图文件夹不存在", code="52000509")
        query_id = self.add_annotation_query(name=None, params=params, stat_id=stat_id)
        query_path = ref_anno_path + "/anno_stat/trans_anno_detail.xls"
        gene_query_path = ref_anno_path + "/anno_stat/gene_anno_detail.xls"
        self.add_annotation_query_detail(query_id=query_id, query_path=query_path, anno_type="transcript")
        self.add_annotation_gene_query_detail(query_id=query_id, query_path=gene_query_path, anno_type="gene")
        query_path = new_anno_path + "/anno_stat/trans_anno_detail.xls"
        gene_query_path = new_anno_path + "/anno_stat/gene_anno_detail.xls"
        self.add_annotation_query_detail(query_id=query_id, query_path=query_path, anno_type="transcript")
        self.add_annotation_gene_query_detail(query_id=query_id, query_path=gene_query_path, anno_type="gene")
        cog_id = self.add_annotation_cog(name=name, params=params)
        r_sum_path = ref_anno_path + "/cog/cog_summary.xls"
        self.add_annotation_cog_detail(cog_id=cog_id, cog_path=r_sum_path, seq_type="ref", anno_type="transcript")
        r_gene_sum_path = ref_anno_path + "/anno_stat/cog_stat/gene_cog_summary.xls"
        self.add_annotation_cog_detail(cog_id=cog_id, cog_path=r_gene_sum_path, seq_type="ref", anno_type="gene")
        n_sum_path = new_anno_path + "/cog/cog_summary.xls"
        self.add_annotation_cog_detail(cog_id=cog_id, cog_path=n_sum_path, seq_type="new", anno_type="transcript")
        n_gene_sum_path = new_anno_path + "/anno_stat/cog_stat/gene_cog_summary.xls"
        self.add_annotation_cog_detail(cog_id=cog_id, cog_path=n_gene_sum_path, seq_type="new", anno_type="gene")
        self.add_annotation_cog_detail_all(cog_id=cog_id, r_cog_path=r_sum_path, n_cog_path=n_sum_path, seq_type="all", anno_type="transcript")
        self.add_annotation_cog_detail_all(cog_id=cog_id, r_cog_path=r_gene_sum_path, n_cog_path=n_gene_sum_path, seq_type="all", anno_type="gene")

        def add_go(go_id, go_path, gene_go_path, seq_type):
            if os.path.exists(go_path) and os.path.exists(gene_go_path):
                stat_level2 = go_path + "/go12level_statistics.xls"
                stat_level3 = go_path + "/go123level_statistics.xls"
                stat_level4 = go_path + "/go1234level_statistics.xls"
                gene_stat_level2 = gene_go_path + "/gene_go12level_statistics.xls"
                gene_stat_level3 = gene_go_path + "/gene_go123level_statistics.xls"
                gene_stat_level4 = gene_go_path + "/gene_go1234level_statistics.xls"
                gos_path = go_path + "/query_gos.list"
                gene_gos_path = gene_go_path + "/gene_gos.list"
                self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=2, level_path=stat_level2)
                self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="gene", level=2, level_path=gene_stat_level2)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=2, go_path=stat_level2)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=3, go_path=stat_level3)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=4, go_path=stat_level4)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="gene", level=2, go_path=gene_stat_level2)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="gene", level=3, go_path=gene_stat_level3)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="gene", level=4, go_path=gene_stat_level4)
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=2, go_path=stat_level2)
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=3, go_path=stat_level3)
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=4, go_path=stat_level4)
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="gene", level=2, go_path=gene_stat_level2)
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="gene", level=3, go_path=gene_stat_level3)
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="gene", level=4, go_path=gene_stat_level4)
                self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="transcript", gos_path=gos_path)
                self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="gene", gos_path=gene_gos_path)
            else:
                self.bind_object.set_error("GO注释的结果文件不存在", code="52000510")
        go_id = self.add_annotation_go(name=name, params=params)
        go_path = ref_anno_path + "/go"
        gene_go_path = ref_anno_path + "/anno_stat/go_stat"
        add_go(go_id=go_id, go_path=go_path, gene_go_path=gene_go_path, seq_type="ref")
        go_path = new_anno_path + "/go"
        gene_go_path = new_anno_path + "/anno_stat/go_stat"
        add_go(go_id=go_id, go_path=go_path, gene_go_path=gene_go_path, seq_type="new")
        r_stat_level2 = ref_anno_path + "/go/go12level_statistics.xls"
        r_stat_level3 = ref_anno_path + "/go/go123level_statistics.xls"
        r_stat_level4 = ref_anno_path + "/go/go1234level_statistics.xls"
        n_stat_level2 = new_anno_path + "/go/go12level_statistics.xls"
        n_stat_level3 = new_anno_path + "/go/go123level_statistics.xls"
        n_stat_level4 = new_anno_path + "/go/go1234level_statistics.xls"
        r_gene_stat_level2 = ref_anno_path + "/anno_stat/go_stat/gene_go12level_statistics.xls"
        r_gene_stat_level3 = ref_anno_path + "/anno_stat/go_stat/gene_go123level_statistics.xls"
        r_gene_stat_level4 = ref_anno_path + "/anno_stat/go_stat/gene_go1234level_statistics.xls"
        n_gene_stat_level2 = new_anno_path + "/anno_stat/go_stat/gene_go12level_statistics.xls"
        n_gene_stat_level3 = new_anno_path + "/anno_stat/go_stat/gene_go123level_statistics.xls"
        n_gene_stat_level4 = new_anno_path + "/anno_stat/go_stat/gene_go1234level_statistics.xls"
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="transcript", level=2, r_go_path=r_stat_level2, n_go_path=n_stat_level2)
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="transcript", level=3, r_go_path=r_stat_level3, n_go_path=n_stat_level3)
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="transcript", level=4, r_go_path=r_stat_level4, n_go_path=n_stat_level4)
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="gene", level=2, r_go_path=r_gene_stat_level2, n_go_path=n_gene_stat_level2)
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="gene", level=3, r_go_path=r_gene_stat_level3, n_go_path=n_gene_stat_level3)
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="gene", level=4, r_go_path=r_gene_stat_level4, n_go_path=n_gene_stat_level4)

        def add_kegg(kegg_id, kegg_path, gene_kegg_path, seq_type):
            if os.path.exists(kegg_path) and os.path.exists(gene_kegg_path):
                layer_path = kegg_path + "/kegg_layer.xls"
                pathway_path = kegg_path + "/pathway_table.xls"
                png_path = kegg_path + "/pathways"
                table_path = kegg_path + "/kegg_table.xls"
                gene_layer_path = gene_kegg_path + "/gene_kegg_layer.xls"
                gene_pathway_path = gene_kegg_path + "/gene_pathway_table.xls"
                gene_png_path = gene_kegg_path + "/gene_pathway"
                gene_table_path = gene_kegg_path + "/gene_kegg_table.xls"
                self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="transcript", categories_path=layer_path)
                self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="gene", categories_path=gene_layer_path)
                self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="transcript", level_path=pathway_path, png_dir=png_path)
                self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="gene", level_path=gene_pathway_path, png_dir=gene_png_path)
                self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="transcript", table_path=table_path)
                self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="gene", table_path=gene_table_path)
            else:
                self.bind_object.set_error("KEGG注释文件不存在", code="52000511")
        kegg_id = self.add_annotation_kegg(name=None, params=params)
        kegg_path = ref_anno_path + "/kegg"
        gene_kegg_path = ref_anno_path + "/anno_stat/kegg_stat"
        add_kegg(kegg_id=kegg_id, kegg_path=kegg_path, gene_kegg_path=gene_kegg_path, seq_type="ref")
        kegg_path = new_anno_path + "/kegg"
        gene_kegg_path = new_anno_path + "/anno_stat/kegg_stat"
        add_kegg(kegg_id=kegg_id, kegg_path=kegg_path, gene_kegg_path=gene_kegg_path, seq_type="new")
        r_cate_path = ref_anno_path + "/kegg/kegg_layer.xls"
        n_cate_path = new_anno_path + "/kegg/kegg_layer.xls"
        r_gene_cate_path = ref_anno_path + "/anno_stat/kegg_stat/gene_kegg_layer.xls"
        n_gene_cate_path = new_anno_path + "/anno_stat/kegg_stat/gene_kegg_layer.xls"
        self.add_annotation_kegg_categories_all(kegg_id=kegg_id, seq_type="all", anno_type="transcript", r_cate_path=r_cate_path, n_cate_path=n_cate_path)
        self.add_annotation_kegg_categories_all(kegg_id=kegg_id, seq_type="all", anno_type="gene", r_cate_path=r_gene_cate_path, n_cate_path=n_gene_cate_path)
        pathway_path = merge_tran_output + "/pathway_table.xls"
        png_path = merge_tran_output + "/all_pathways"
        gene_pathway_path = merge_gene_output + "/pathway_table.xls"
        gene_png_path = merge_gene_output + "/all_pathways"
        self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type="all", anno_type="transcript", level_path=pathway_path, png_dir=png_path)
        self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type="all", anno_type="gene", level_path=gene_pathway_path, png_dir=gene_png_path)


    @report_check
    def add_annotation_stat(self, name=None, params=None, seq_type=None, database=None, result_dir=None, taxon='Animals' ):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationStat_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'version': self.version,
            'result_dir': result_dir,
            "taxonomy": taxon,
            'status': 'start',
            'desc': '注释统计主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            # 'seq_type': seq_type,
            'database': database
        }
        collection = self.db['sg_annotation_stat']
        stat_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add ref_annotation_stat!")
        return stat_id

    @report_check
    def add_annotation_stat_detail(self, stat_id, stat_path, venn_path, exp_level, gene_exp=None, trans_exp=None):
        """
        database: 进行统计的数据库
        stat_path: all_annotation_statistics.xls
        venn_path: venn图目录
        """
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="52000512")
        if not os.path.exists(stat_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(stat_path), code="52000513")
        if not os.path.exists(venn_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(venn_path), code="52000514")


        if self.version == "v2":
            exp_g_set = set()
            with open(gene_exp, 'r') as gene_exp_f:
                gene_exp_f.readline()
                for line in gene_exp_f:
                    exps = line.strip().split("\t")[1:]
                    if sum(map(float, exps)) != 0:
                        exp_g_set.add(line.strip().split("\t")[0])
                    else:
                        pass
            exp_t_set = set()
            with open(trans_exp, 'r') as trans_exp_f:
                trans_exp_f.readline()
                for line in trans_exp_f:
                    exps = line.strip().split("\t")[1:]
                    if sum(map(float, exps)) != 0:
                        exp_t_set.add(line.strip().split("\t")[0])
                    else:
                        pass

        gene_dict = dict()
        trans_dict = dict()

        data_list = []
        database_venn = {
            'NR': 'nr',
            'Swiss-Prot': 'swissprot',
            'Swiss-prot': 'swissprot',
            'Pfam': 'pfam',
            'KEGG': 'kegg',
            'GO': 'go',
            'COG': 'cog',
        }
        tail = "gene"
        database_venn = {
            'NR': 'nr/nr_venn_{}.txt'.format(tail),
            'Swiss-Prot': 'swissprot/swissprot_venn_{}.txt'.format(tail),
            'Swiss-prot': 'swissprot/swissprot_venn_{}.txt'.format(tail),
            'Pfam': 'pfam/pfam_venn_{}.txt'.format(tail),
            'KEGG': 'kegg/kegg_venn_{}.txt'.format(tail),
            'GO': 'go/go_venn_{}.txt'.format(tail),
            'COG': 'cog/cog_venn_{}.txt'.format(tail),
        }
        tail = "tran"
        database_venn_tran = {
            'NR': 'nr/nr_venn_{}.txt'.format(tail),
            'Swiss-Prot': 'swissprot/swissprot_venn_{}.txt'.format(tail),
            'Swiss-prot': 'swissprot/swissprot_venn_{}.txt'.format(tail),
            'Pfam': 'pfam/pfam_venn_{}.txt'.format(tail),
            'KEGG': 'kegg/kegg_venn_{}.txt'.format(tail),
            'GO': 'go/go_venn_{}.txt'.format(tail),
            'COG': 'cog/cog_venn_{}.txt'.format(tail),
        }

        with open(stat_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                if exp_level.lower() == "gene":
                    data = [
                        ('stat_id', stat_id),
                        ('type', line[0]),
                        ('gene', int(line[2])),
                        ('gene_percent', round(float(line[4]), 4)),
                        ('anno_type', 'all')
                    ]
                else:
                    data = [
                        ('stat_id', stat_id),
                        ('type', line[0]),
                        ('transcript', int(line[1])),
                        ('gene', int(line[2])),
                        ('transcript_percent', round(float(line[3]), 4)),
                        ('gene_percent', round(float(line[4]), 4)),
                        ('anno_type', 'all')
                    ]

                venn_list, gene_venn_list = None, None
                database = ["nr", "swissprot", "pfam", "kegg", "go", "string", "cog"]
                # print database_venn
                # print line[0]
                if database_venn.has_key(line[0]):
                    db = line[0]
                    venn = venn_path + "/" + database_venn_tran[line[0]]
                    gene_venn = venn_path + "/" + database_venn[line[0]]
                    # print venn
                    # print gene_venn
                    if os.path.exists(venn) and os.path.exists(gene_venn):
                        with open(venn, "rb") as f:
                            venn_list = f.readline().strip('\n')
                            for line in f:
                                venn_list += ',{}'.format(line.strip('\n'))
                            trans_dict[db] = set(venn_list.split(","))
                        with open(gene_venn, "rb") as f:
                            gene_venn_list = f.readline().strip('\n')
                            for line in f:
                                gene_venn_list += ',{}'.format(line.strip('\n'))
                            gene_dict[db] = set(gene_venn_list.split(","))
                        data.append(("gene_list", gene_venn_list))
                        if len(gene_venn_list) + len(venn_list) < 15000000:
                            data.append(("transcript_list", venn_list))
                    else:
                        self.bind_object.set_error("%s对应的venn.txt文件不存在", variables=(line[0]), code="52000515")
                data = SON(data)
                data_list.append(data)

        if self.version == "v2":
            anno_gene_list = set()
            anno_trans_list = set()
            # print trans_dict
            for db in ['GO', 'KEGG', 'COG', 'NR', 'Swiss-Prot', 'Pfam']:
                trans_set = exp_t_set & set(trans_dict[db])
                gene_set = exp_g_set & set(gene_dict[db])
                data = [
                    ('stat_id', stat_id),
                    ('type', db),
                    ('transcript', len(trans_set)),
                    ('gene', len(gene_set)),
                    ('transcript_percent', round(float(len(trans_set))/float(len(exp_t_set)), 4)),
                    ('gene_percent', round(float(len(gene_set))/float(len(exp_g_set)), 4)),
                    ('anno_type', 'exp'),
                    ('gene_list', ",".join(list(gene_set))),
                    ('transcript_list', ",".join(list(trans_set)))
                ]
                anno_gene_list |= gene_set
                anno_trans_list |= trans_set
                data = SON(data)
                data_list.append(data)
            data = [
                ('stat_id', stat_id),
                ('type', 'Total_anno'),
                ('transcript', len(anno_trans_list)),
                ('gene', len(anno_gene_list)),
                ('transcript_percent', round(float(len(anno_trans_list))/float(len(exp_t_set)), 4)),
                ('gene_percent', round(float(len(anno_gene_list))/float(len(exp_g_set)), 4)),
                ('anno_type', 'exp')
            ]
            data = SON(data)
            data_list.append(data)
            data = [
                ('stat_id', stat_id),
                ('type', 'Total'),
                ('transcript', len(exp_t_set)),
                ('gene', len(exp_g_set)),
                ('transcript_percent', 1.0),
                ('gene_percent', 1.0),
                ('anno_type', 'exp')
            ]
            data = SON(data)
            data_list.append(data)
            with open(stat_path.split(".xls")[0] + "_detail.xls", 'w') as f_stat_o:
                data_dict = dict()
                for data in data_list:
                    db_type = data['type']
                    if db_type in data_dict:
                        pass
                    else:
                        data_dict[db_type] = dict()
                    if data['anno_type'] == 'exp':
                        data_dict[db_type]['exp'] = "{}({})\t{}({})".format(data['gene'], data['gene_percent'], data['transcript'], data['transcript_percent'])
                    else:
                        data_dict[db_type]['all'] = "{}({})\t{}({})".format(data['gene'], data['gene_percent'], data['transcript'], data['transcript_percent'])
                f_stat_o.write("\tExp_Unigene number(percent)\tExp_Transcript number(percent)\tAll_Unigene number(percent)\tAll_Transcript number(percent)")
                for db_type in ['GO', 'KEGG', 'COG', 'NR', 'Swiss-Prot', 'Pfam', 'Total_anno', 'Total']:
                    if db_type in data_dict:
                        f_stat_o.write("\t".join([db_type, data_dict[db_type]['exp'], data_dict[db_type]['all']]) + "\n")


        try:
            collection = self.db['sg_annotation_stat_detail']
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入注释统计信息失败:%s" , variables=( stat_path), code="52000516")
        else:
            self.bind_object.logger.info("导入注释统计信息成功：%s" % (stat_path))


        '''
        try:
            collection = self.db['sg_annotation_stat_detail']
            if len(data_list) > 5000:
                for i in range(0, len(data_list), 3000):
                    tmp_list = data_list[i: i+3000]
                    collection.insert(tmp_list)
            else:
                collection.insert(data_list)
        except Exception, e:
            raise Exception("导入注释统计信息失败:%s" % stat_path)
        else:
            self.bind_object.logger.info("导入注释统计信息成功：%s" % (stat_path))
        '''

    @report_check
    def add_stat_detail(self, old_stat_id, stat_id, nr_evalue, gene_nr_evalue, sw_evalue, gene_sw_evalue):
        """
        注释重运行时注释统计导表sg_annotation_stat_detail
        """
        if not isinstance(old_stat_id, ObjectId):
            if isinstance(old_stat_id, types.StringTypes):
                old_stat_id = ObjectId(old_stat_id)
            else:
                self.bind_object.set_error('old_stat_id必须为ObjectId对象或其对应的字符串！', code="52000517")
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="52000518")
        collection = self.db["sg_annotation_stat_detail"]
        results = collection.find({"stat_id": old_stat_id})
        data_list, data = [], []
        for result in results:
            db = result["type"]
            if db == "total":
                total_tran = result["transcript"]
                total_gene = result["gene"]
            # 增加数据库类型，区分有注释的基因和有分类的基因
            # if db in ["pfam", "total_anno", "total_anno_nsp","total_class", "total"]:
            if db in ["pfam", "total_anno",  "total"]:
                data = [
                    ('stat_id', stat_id),
                    ('type', result["type"]),
                    ('transcript', result["transcript"]),
                    ('gene', result["gene"]),
                    ('transcript_percent', result["transcript_percent"]),
                    ('gene_percent', result["gene_percent"]),
                    ('gene_list', result["gene_list"]),
                    ('transcript_list', result["transcript_list"])
                ]
                data = SON(data)
                data_list.append(data)
        nr_ids = self.stat(stat_path=nr_evalue)
        gene_nr_ids = self.stat(stat_path=gene_nr_evalue)
        data = [
            ('stat_id', stat_id),
            ('type', "nr"),
            ('transcript', len(nr_ids)),
            ('gene', len(gene_nr_ids)),
            ('transcript_percent', round(float(len(nr_ids))/total_tran, 4)),
            ('gene_percent', round(float(len(gene_nr_ids))/total_gene, 4)),
            ('gene_list', ",".join(gene_nr_ids)),
            ('transcript_list', ",".join(nr_ids))
        ]
        data = SON(data)
        data_list.append(data)
        sw_ids = self.stat(stat_path=sw_evalue)
        gene_sw_ids = self.stat(stat_path=gene_sw_evalue)
        data = [
            ('stat_id', stat_id),
            ('type', "swissprot"),
            ('transcript', len(sw_ids)),
            ('gene', len(gene_sw_ids)),
            ('transcript_percent', round(float(len(sw_ids))/total_tran, 4)),
            ('gene_percent', round(float(len(gene_sw_ids))/total_gene, 4)),
            ('gene_list', ",".join(gene_sw_ids)),
            ('transcript_list', ",".join(sw_ids))
        ]
        data = SON(data)
        data_list.append(data)
        try:
            collection = self.db['sg_annotation_stat_detail']
            if len(data_list) > 5000:
                for i in range(0, len(data_list), 3000):
                    tmp_list = data_list[i: i+3000]
                    collection.insert_many(tmp_list)
            else:
                    collection.insert_many(data_list)
        except:
            self.bind_object.set_error("导入注释统计信息出错", code="52000519")
        else:
            self.bind_object.logger.info("导入注释统计信息成功")

    def stat(self, stat_path):
        with open(stat_path, "rb") as f:
            id_list = []
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                q_id = line[5]
                id_list.append(q_id)
        id_list = list(set(id_list))
        return id_list

    @report_check
    def add_annotation_blast(self, name=None, params=None, stat_id=None, result_dir=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="52000520")
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationBlast_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': params,
            'version': self.version,
            'result_dir': result_dir,
            'status': 'start',
            'desc': 'blast最佳比对结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'stat_id': stat_id
        }
        collection = self.db['sg_annotation_blast']
        blast_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add ref_annotation_blast!")
        return blast_id

    @report_check
    def add_annotation_blast_detail(self, blast_id, seq_type, anno_type, database, blast_path):
        if not isinstance(blast_id, ObjectId):
            if isinstance(blast_id, types.StringTypes):
                blast_id = ObjectId(blast_id)
            else:
                self.bind_object.set_error('blast_id必须为ObjectId对象或其对应的字符串！', code="52000521")
        if not os.path.exists(blast_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(blast_path), code="52000522")
        data_list = []
        with open(blast_path, 'r') as f:
            lines = f.readlines()
            flag = None
            for line in lines[1:]:
                line = line.strip().split('\t')
                query_name = line[5]
                hit_name = line[10]
                if flag == query_name:
                    pass
                else:
                    flag = query_name
                    data = {
                        'blast_id': blast_id,
                        'seq_type': seq_type,
                        'anno_type': anno_type,
                        'database': database,
                        'score': float(line[0]),
                        'e_value': float(line[1]),
                        'hsp_len': int(line[2]),
                        'identity_rate': round(float(line[3]), 4),
                        'similarity_rate': round(float(line[4]), 4),
                        'query_id': line[5],
                        'q_len': int(line[6]),
                        'q_begin': line[7],
                        'q_end': line[8],
                        'q_frame': line[9],
                        'hit_name': line[10],
                        'hit_len': int(line[11]),
                        'hsp_begin': line[12],
                        'hsp_end': line[13],
                        'hsp_frame': line[14],
                        'description': line[15]
                    }
                    if self._project_type == 'denovo_rna_v2' and anno_type == 'T':
                        data['gene_id'] = self.trans_gene[line[5]]
                        data['is_gene'] = self.trans_isgene[line[5]]
                    collection = self.db['sg_annotation_blast_detail']
                    collection.insert_one(data).inserted_id
        self.bind_object.logger.info("导入blast信息：%s成功!" % (blast_path))

    @report_check
    def add_annotation_blast_nr_detail(self, blast_id, seq_type, anno_type, database, blast_path):
        if not isinstance(blast_id, ObjectId):
            if isinstance(blast_id, types.StringTypes):
                blast_id = ObjectId(blast_id)
            else:
                self.bind_object.set_error('blast_id必须为ObjectId对象或其对应的字符串！', code="52000523")
        if not os.path.exists(blast_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(blast_path), code="52000524")
        data_list = []
        with open(blast_path, 'r') as f:
            lines = f.readlines()
            flag = None
            for line in lines[1:]:
                line = line.strip().split('\t')
                query_name = line[5]
                hit_name = line[10]
                if flag == query_name:
                    pass
                else:
                    flag = query_name
                    data = {
                        'nr_id': blast_id,
                        # 'seq_type': seq_type,
                        # 'anno_type': anno_type,
                        'database': database,
                        'score': float(line[0]),
                        'e_value': float(line[1]),
                        'hsp_len': int(line[2]),
                        'identity_rate': round(float(line[3]), 4),
                        'similarity_rate': round(float(line[4]), 4),
                        'transcript_id': line[5],
                        'q_len': int(line[6]),
                        # 'q_begin': line[7],
                        # 'q_end': line[8],
                        # 'q_frame': line[9],
                        'hit_name': line[10],
                        # 'hit_len': int(line[11]),
                        # 'hsp_begin': line[12],
                        # 'hsp_end': line[13],
                        # 'hsp_frame': line[14],
                        'description': line[15]
                    }
                    if self._project_type == 'denovo_rna_v2' and anno_type == 'T':
                        data['gene_id'] = self.trans_gene[line[5]]
                        data['is_gene'] = self.trans_isgene[line[5]]
                    collection = self.db['sg_annotation_nr_detail']
                    collection.insert_one(data).inserted_id
        self.bind_object.logger.info("导入blast信息：%s成功!" % (blast_path))


    @report_check
    def add_annotation_nr_species_detail(self, nr_id, gene_species_count, trans_species_count, seq_type="new", anno_type="G"):
        if not isinstance(nr_id, ObjectId):
            if isinstance(nr_id, types.StringTypes):
                nr_id = ObjectId(nr_id)
            else:
                self.bind_object.set_error('nr_id必须为ObjectId对象或其对应的字符串！', code="52000525")
        if not os.path.exists(gene_species_count):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(gene_species_count), code="52000526")
        if not os.path.exists(trans_species_count):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(trans_species_count), code="52000527")
        data_list = []
        with open(gene_species_count, 'r') as gene_spe_f:
            header = gene_spe_f.readline().strip().split("\t")
            for line in gene_spe_f:
                cols = line.strip().split("\t")
                data = zip(header, cols)
                data.extend([("anno_type", "G"), ("nr_id", nr_id)])
                data = SON(data)
                data_list.append(data)
        with open(trans_species_count, 'r') as trans_spe_f:
            header = trans_spe_f.readline().strip().split("\t")
            for line in trans_spe_f:
                cols = line.strip().split("\t")
                data = zip(header, cols)
                data.extend([("anno_type", "T"), ("nr_id", nr_id)])

                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db["sg_annotation_nr_species"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入%s信息出错:%s" % ("sg_annotation_nr_detail", e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % ("sg_annotation_nr_detail"))


    @report_check
    def add_annotation_blast_swissprot_detail(self, blast_id, seq_type, anno_type, database, blast_path):
        if not isinstance(blast_id, ObjectId):
            if isinstance(blast_id, types.StringTypes):
                blast_id = ObjectId(blast_id)
            else:
                self.bind_object.set_error('blast_id必须为ObjectId对象或其对应的字符串！', code="52000528")
        if not os.path.exists(blast_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(blast_path), code="52000529")
        data_list = []
        with open(blast_path, 'r') as f:
            lines = f.readlines()
            flag = None
            for line in lines[1:]:
                line = line.strip().split('\t')
                query_name = line[5]
                hit_name = line[10]
                if flag == query_name:
                    pass
                else:
                    flag = query_name
                    data = {
                        'swissprot_id': blast_id,
                        # 'seq_type': seq_type,
                        # 'anno_type': anno_type,
                        'database': database,
                        'score': float(line[0]),
                        'e_value': float(line[1]),
                        'hsp_len': int(line[2]),
                        'identity_rate': round(float(line[3]), 4),
                        'similarity_rate': round(float(line[4]), 4),
                        'transcript_id': line[5],
                        'q_len': int(line[6]),
                        # 'q_begin': line[7],
                        # 'q_end': line[8],
                        # 'q_frame': line[9],
                        'hit_name': line[10],
                        # 'hit_len': int(line[11]),
                        # 'hsp_begin': line[12],
                        # 'hsp_end': line[13],
                        # 'hsp_frame': line[14],
                        'description': line[15]
                    }
                    if self._project_type == 'denovo_rna_v2' and anno_type == 'T':
                        data['gene_id'] = self.trans_gene[line[5]]
                        data['is_gene'] = self.trans_isgene[line[5]]
                    collection = self.db['sg_annotation_swissprot_detail']
                    collection.insert_one(data).inserted_id
        self.bind_object.logger.info("导入blast信息：%s成功!" % (blast_path))

    @report_check
    def add_annotation_nr(self, name=None, params=None, stat_id=None, result_dir=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="52000530")
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationNr_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'version': self.version,
            'result_dir': result_dir,
            'status': 'start',
            'desc': 'nr注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'stat_id': stat_id
        }
        collection = self.db['sg_annotation_nr']
        nr_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add sg_annotation_nr!")
        return nr_id

    @report_check
    def add_annotation_nr_pie(self, nr_id, species_path, evalue_path, similar_path, seq_type, anno_type):
        if not isinstance(nr_id, ObjectId):
            if isinstance(nr_id, types.StringTypes):
                nr_id = ObjectId(nr_id)
            else:
                self.bind_object.set_error('nr_id必须为ObjectId对象或其对应的字符串！', code="52000531")
        if not os.path.exists(evalue_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(evalue_path), code="52000532")
        if not os.path.exists(similar_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(similar_path), code="52000533")
        evalue, evalue_list,species_list_a, species_a, similar, similar_list = [], [], [], [], [], []
        data_list = []
        with open(evalue_path, "r") as f1, open(similar_path, "r") as f2:
            lines1 = f1.readlines()
            lines2 = f2.readlines()
            gene_num = sum([int(line.strip().split('\t')[1]) for line in lines1[1:]])
            for line1 in lines1[1:]:
                line1 = line1.strip().split('\t')
                value = {"key": line1[0], "value": int(line1[1]), "pct": round(float(line1[1])/gene_num * 100, 4)}
                try:
                    value_list = {"key": line1[0], "value": line1[2]}
                except:
                    value_list = {"key": line1[0], "value": None}
                evalue.append(value)
                evalue_list.append(value_list)
            for line2 in lines2[1:]:
                line2 = line2.strip().split('\t')
                similarity = {"key": line2[0], "value": int(line2[1]), "pct": round(float(line2[1])/gene_num * 100, 4)}
                try:
                    similarity_list = {"key": line2[0], "value": line2[2]}
                except:
                    similarity_list = {"key": line2[0], "value": None}
                similar.append(similarity)
                similar_list.append(similarity_list)

        if self.version == "v2":
            data_list_species = []
            with open(species_path, "r") as f3:
                lines3 = f3.readlines()
                for line3 in lines3[1:15]:
                    line3 = line3.strip().split("\t")
                    species = {"key": line3[0], "value": int(line3[1])}
                    try:
                        species_list = {"key": line3[0], "value": line3[1]}
                    except:
                        species_list = {"key": line3[0], "value": None}
                    species_a.append(species)
                    species_list_a.append(species_list)
                    data_species = [
                        ('nr_id', nr_id),
                        ('anno_type', anno_type),
                        ('species', line3[0]),
                        ('num', line3[1]),
                        ('percent', line3[2]),
                        ('class', line3[3])
                    ]
                    data_species = SON(data_species)
                    data_list_species.append(data_species)
                other = 0
                other_pct = 0
                species_other = []
                species_other_list = []
                for line3 in lines3[15:]:
                    line3 = line3.strip().split("\t")
                    try:
                        species_list = {"key": line3[0], "value": line3[1]}
                    except:
                        species_list = {"key": line3[0], "value": None}
                    other += int(line3[1])
                    other_pct += float(line3[2])
                    species_other.append(species)
                    species_other_list.append(species_list)
                species = {"key": 'other', "value": other}
                species_list = {"key": 'other', "value": None}
                species_a.append(species)
                species_list_a.append(species_list)
                data_species = SON([
                    ('nr_id', nr_id),
                    ('anno_type', anno_type),
                    ('species', "other"),
                    ('num', other),
                    ('percent', other_pct),
                    ('class', "other")
                ])
                data_list_species.append(data_species)

            data = [
                ('nr_id', nr_id),
                #('seq_type', seq_type),
                ('anno_type', anno_type),
                ('e_value', evalue),
                ('similar', similar),
                ('evalue_list', evalue_list),
                ('similar_list', similar_list),
                ('species', species_a),
                ('species_list', species_list_a),
                ('species_other', species_other),
                ('species_other_list', species_other_list),
            ]
        else:
            with open(species_path, "r") as f3:
                lines3 = f3.readlines()
                for line3 in lines3[1:15]:
                    line3 = line3.strip().split("\t")
                    if anno_type == "T":
                        species = {"key": line3[0], "value": int(line3[1])}
                        try:
                            species_list = {"key": line3[0], "value": line3[5]}
                        except:
                            species_list = {"key": line3[0], "value": None}
                    elif anno_type == "G":
                        species = {"key": line3[0], "value": int(line3[2])}
                        try:
                            species_list = {"key": line3[0], "value": line3[6]}
                        except:
                            species_list = {"key": line3[0], "value": None}
                    species_a.append(species)
                    species_list_a.append(species_list)
                other = 0
                species_other = []
                species_other_list = []
                for line3 in lines3[16:]:
                    line3 = line3.strip().split("\t")
                    if anno_type == "T":
                        species = {"key": line3[0], "value": int(line3[1])}
                        try:
                            species_list = {"key": line3[0], "value": line3[5]}
                        except:
                            species_list = {"key": line3[0], "value": None}
                        other += int(line3[1])
                    elif anno_type == "G":
                        try:
                            species_list = {"key": line3[0], "value": line3[6]}
                        except:
                            species_list = {"key": line3[0], "value": None}
                        other += int(line3[2])
                    species_other.append(species)
                    species_other_list.append(species_list)
                species = {"key": 'other', "value": other}
                species_list = {"key": 'other', "value": None}
                species_a.append(species)
                species_list_a.append(species_list)

            data = [
                ('nr_id', nr_id),
                #('seq_type', seq_type),
                ('anno_type', anno_type),
                ('e_value', evalue),
                ('similar', similar),
                ('evalue_list', evalue_list),
                ('similar_list', similar_list),
                ('species', species_a),
                ('species_list', species_list_a),
                ('species_other', species_other),
                ('species_other_list', species_other_list),
            ]

        data = SON(data)
        data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_nr_pie']
                if len(data_list) > 5000:
                    for i in range(0, len(data_list), 3000):
                        tmp_list = data_list[i: i+3000]
                        collection.insert_many(tmp_list)
                else:
                    collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入nr库注释作图信息evalue,similar：%s、%s出错!" , variables=(evalue_path, similar_path), code="52000534")
            else:
                self.bind_object.logger.info("导入nr库注释作图信息evalue,similar：%s、%s成功!" % (evalue_path, similar_path))


        '''
        if data_list_species:
            try:
                collection = self.db['sg_annotation_nr_species']
                if len(data_list_species) > 5000:
                    for i in range(0, len(data_list), 3000):
                        tmp_list = data_list_species[i: i+3000]
                        collection.insert_many(tmp_list)
                else:
                    collection.insert_many(data_list_species)
            except Exception, e:
                raise Exception("导入nr库注释作图信息evalue,similar：%s、%s出错!" % (species_path, e))
            else:
                self.bind_object.logger.info("导入nr库注释作图信息evalue,similar：%s、%s成功!" % (species_path, species_path))
                data_list.append(data)
        '''


    @report_check
    def add_annotation_swissprot(self, name=None, params=None, stat_id=None, result_dir=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="52000535")
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationSwissprot_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'version': self.version,
            'result_dir': result_dir,
            'status': 'start',
            'desc': 'swissprot注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'stat_id': stat_id
        }
        collection = self.db['sg_annotation_swissprot']
        swissprot_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add sg_annotation_swissprot!")
        return swissprot_id

    @report_check
    def add_annotation_swissprot_pie(self, swissprot_id, evalue_path, similar_path, seq_type, anno_type):
        """
        """
        if not isinstance(swissprot_id, ObjectId):
            if isinstance(swissprot_id, types.StringTypes):
                swissprot_id = ObjectId(swissprot_id)
            else:
                self.bind_object.set_error('swissprot_id必须为ObjectId对象或其对应的字符串！', code="52000536")
        if not os.path.exists(evalue_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(evalue_path), code="52000537")
        if not os.path.exists(similar_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(similar_path), code="52000538")
        evalue, evalue_list, similar, similar_list = [], [], [], []
        data_list = []
        with open(evalue_path, "r") as f1, open(similar_path, "r") as f2:
            lines1 = f1.readlines()
            lines2 = f2.readlines()
            gene_num = sum([int(line.strip().split('\t')[1]) for line in lines1[1:]])
            for line1 in lines1[1:]:
                line1 = line1.strip().split('\t')
                value = {"key": line1[0], "value": int(line1[1]), "pct": round(float(line1[1])/gene_num * 100, 4)}

                try:
                    value_list = {"key": line1[0], "value": line1[2]}
                except:
                    value_list = {"key": line1[0], "value": None}
                evalue.append(value)
                evalue_list.append(value_list)
            for line2 in lines2[1:]:
                line2 = line2.strip().split('\t')
                similarity = {"key": line2[0], "value": int(line2[1]), "pct": round(float(line2[1])/gene_num * 100, 4)}
                try:
                    similarity_list = {"key": line2[0], "value": line2[2]}
                except:
                    similarity_list = {"key": line2[0], "value": None}
                similar.append(similarity)
                similar_list.append(similarity_list)
        data = [
            ('swissprot_id', swissprot_id),
            # ('seq_type', seq_type),
            ('anno_type', anno_type),
            ('e_value', evalue),
            ('similar', similar),
            ('evalue_list', evalue_list),
            ('similar_list', similar_list),
        ]
        data = SON(data)
        data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_swissprot_pie']
                if len(data_list) > 5000:
                    for i in range(0, len(data_list), 3000):
                        tmp_list = data_list[i: i+3000]
                        collection.insert_many(tmp_list)
                else:
                    collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入swissprot库注释作图信息evalue,similar：%s、%s出错!" , variables=(evalue_path, similar_path), code="52000539")
            else:
                self.bind_object.logger.info("导入swissprot库注释作图信息evalue,similar：%s、%s成功!" % (evalue_path, similar_path))

    @report_check
    def add_annotation_pfam(self, name=None, params=None, stat_id=None, result_dir=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="52000540")
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationPfam_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'version': self.version,
            'result_dir': result_dir,
            'status': 'start',
            'desc': 'pfam注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'stat_id': stat_id
        }
        collection = self.db['sg_annotation_pfam']
        pfam_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add sg_annotation_pfam!")
        return pfam_id

    def add_annotation_pfam_detail(self, pfam_id, pfam_path, seq_type, anno_type):
        """
        pfam_path: pfam_domain
        """
        if not isinstance(pfam_id, ObjectId):
            if isinstance(pfam_id, types.StringTypes):
                pfam_id = ObjectId(pfam_id)
            else:
                self.bind_object.set_error('pfam_id必须为ObjectId对象或其对应的字符串！', code="52000541")
        if not os.path.exists(pfam_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(pfam_path), code="52000542")
        data_list = []
        with open(pfam_path, "r") as f:
            lines = f.readlines()
            last_seq_id = ''
            last_pfam_id = ''
            for line in lines[1:]:
                line = line.strip().split("\t")
                if line[2] != last_seq_id or line[3] != last_pfam_id:
                    data = [
                        ('pfam_id', pfam_id),
                        # ('seq_type', seq_type),
                        # ('anno_type', anno_type),
                        ('transcript_id', line[0]),
                        ('pfam', line[2]),
                        ('domain', line[3]),
                        ('description', line[4]),
                        ('protein_id', line[1]),
                        ('e_value', float(line[9])),
                        ('length', int(line[6])-int(line[5])),
                        ('protein_start', int(line[5])),
                        ('protein_end', int(line[6])),
                        ('pfam_start', int(line[7])),
                        ('pfam_end', int(line[8])),
                    ]
                    if self._project_type == 'denovo_rna_v2' and anno_type == 'T':
                        data.append(('gene_id', self.trans_gene[line[0]]))
                        data.append(('is_gene', self.trans_isgene[line[0]]))
                    data = SON(data)
                    data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_pfam_detail']
                if len(data_list) > 5000:
                    for i in range(0, len(data_list), 3000):
                        tmp_list = data_list[i: i+3000]
                        collection.insert_many(tmp_list)
                else:
                    collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入pfam注释信息:%s失败！" , variables=( pfam_path), code="52000543")
            else:
                self.bind_object.logger.info("导入pfam注释信息:%s成功" % pfam_path)

    @report_check
    def add_annotation_pfam_bar(self, pfam_id, pfam_path, seq_type, anno_type):
        pfam = []
        domain = {}
        with open(pfam_path, "rb") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                if line[3] not in pfam:
                    pfam.append(line[3])
                    domain[line[3]] = 1
                else:
                    domain[line[3]] += 1
        if not isinstance(pfam_id, ObjectId):
            if isinstance(pfam_id, types.StringTypes):
                pfam_id = ObjectId(pfam_id)
            else:
                self.bind_object.set_error('pfam_id必须为ObjectId对象或其对应的字符串！', code="52000544")
        if not os.path.exists(pfam_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(pfam_path), code="52000545")
        data_list = []
        for i in domain:
            data = [
                ('pfam_id', pfam_id),
                # ('seq_type', seq_type),
                ('anno_type', anno_type),
                ('domain', i),
                ('num', domain[i])
            ]
            data = SON(data)
            data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_pfam_bar']
                if len(data_list) > 5000:
                    for i in range(0, len(data_list), 3000):
                        tmp_list = data_list[i: i+3000]
                        collection.insert_many(tmp_list)
                else:
                    collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入pfam注释信息:%s失败！" , variables=( pfam_path), code="52000546")
            else:
                self.bind_object.logger.info("导入pfam注释信息:%s成功！" % pfam_path)

    @report_check
    def add_annotation_cog(self, name=None, params=None, result_dir=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationCog_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'version': self.version,
            'result_dir': result_dir,
            'status': 'start',
            'desc': 'cog注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        }
        collection = self.db['sg_annotation_cog']
        cog_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add ref_annotation_cog!")
        return cog_id

    @report_check
    def add_annotation_cog_detail(self, cog_id, cog_path, seq_type, anno_type):
        '''
        cog_path: cog_summary.xls
        seq_type: ref/new
        anno_type: transcript/gene
        '''
        if not isinstance(cog_id, ObjectId):
            if isinstance(cog_id, types.StringTypes):
                cog_id = ObjectId(cog_id)
            else:
                self.bind_object.set_error('cog_id必须为ObjectId对象或其对应的字符串！', code="52000547")
        if not os.path.exists(cog_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(cog_path), code="52000548")
        data_list = list()
        with open(cog_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [
                    ('cog_id', cog_id),
                    # ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('type', line[0]),
                    ('function_categories', "[" + line[2] + "]" + " " + line[1]),
                    ('cog', int(line[3])),
                ]
                '''
                    ('function_categories', line[1]),
                    ('cog', int(line[2])),
                    ('nog', int(line[3]))
                '''

                try:
                    data.append(('cog_list', line[4]))
                except:
                    data.append(('cog_list', None))
                '''
                try:
                    data.append(('nog_list', line[5]))
                except:
                    data.append(('nog_list', None))
                '''
                data = SON(data)
                data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_cog_detail']
                if len(data_list) > 5000:
                    for i in range(0, len(data_list), 3000):
                        tmp_list = data_list[i: i+3000]
                        collection.insert_many(tmp_list)
                else:
                    collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入cog注释信息：%s出错!" , variables=(cog_path), code="52000549")
            else:
                self.bind_object.logger.info("导入cog注释信息：%s成功!" % (cog_path))

    @report_check
    def add_annotation_cog_detail_all(self, cog_id, r_cog_path, n_cog_path, seq_type, anno_type):
        '''
        r_cog_path: cog_summary.xls(ref)
        n_cog_path: cog_summary.xls(new)
        '''
        if not isinstance(cog_id, ObjectId):
            if isinstance(cog_id, types.StringTypes):
                cog_id = ObjectId(cog_id)
            else:
                self.bind_object.set_error('cog_id必须为ObjectId对象或其对应的字符串！', code="52000550")
        if not os.path.exists(r_cog_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(r_cog_path), code="52000551")
        if not os.path.exists(n_cog_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(n_cog_path), code="52000552")
        first = ['INFORMATION STORAGE AND PROCESSING', 'CELLULAR PROCESSES AND SIGNALING', 'METABOLISM', 'POORLY CHARACTERIZED']
        func_type = {
            'INFORMATION STORAGE AND PROCESSING': sorted(['J', 'A', 'K', 'L', 'B']),
            'CELLULAR PROCESSES AND SIGNALING': sorted(['D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'W', 'U', 'O']),
            'METABOLISM': sorted(['C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q']),
            'POORLY CHARACTERIZED': sorted(['R', 'S']),
        }
        func_decs = {
            'J': 'Translation, ribosomal structure and biogenesis',
            'A': 'RNA processing and modification', 'K': 'Transcription',
            'L': 'Replication, recombination and repair',
            'B': 'Chromatin structure and dynamics',
            'D': 'Cell cycle control, cell division, chromosome partitioning',
            'Y': 'Nuclear structure', 'V': 'Defense mechanisms', 'T': 'Signal transduction mechanisms',
            'M': 'Cell wall/membrane/envelope biogenesis',
            'N': 'Cell motility', 'Z': 'Cytoskeleton', 'W': 'Extracellular structures',
            'U': 'Intracellular trafficking, secretion, and vesicular transport',
            'O': 'Posttranslational modification, protein turnover, chaperones',
            'C': 'Energy production and conversion', 'G': 'Carbohydrate transport and metabolism',
            'E': 'Amino acid transport and metabolism', 'F': 'Nucleotide transport and metabolism',
            'H': 'Coenzyme transport and metabolism', 'I': 'Lipid transport and metabolism',
            'P': 'Inorganic ion transport and metabolism',
            'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
            'R': 'General function prediction only', 'S': 'Function unknown'
        }
        data_list = list()
        cate = {}
        funlist = {'COG': {}, 'NOG': {}}
        cog_fun, nog_fun = {}, {}
        with open(r_cog_path, 'r') as f, open(n_cog_path, 'r') as n:
            lines = f.readlines()
            items = n.readlines()
            for line in lines[2:]:
                line = line.strip().split('\t')
                m = re.match(r"\[(.+)\].+$", line[1])
                if m:
                    fun1 = m.group(1)
                    try:
                        funlist['COG'][fun1] = line[4].split(";")
                    except:
                        funlist['COG'][fun1] = []
                    try:
                        funlist['NOG'][fun1] = line[5].split(";")
                    except:
                        funlist['NOG'][fun1] = []
            for item in items[2:]:
                item = item.strip().split('\t')
                m = re.match(r"\[(.+)\].+$", item[1])
                if m:
                    fun1 = m.group(1)
                    if fun1 in funlist['COG']:
                        try:
                            cog_ids = item[4].split(";")
                            if cog_ids:
                                for cog in cog_ids:
                                    if cog not in funlist['COG'][fun1]:
                                        funlist['COG'][fun1].append(cog)
                        except:
                            pass
                    else:
                        try:
                            funlist['COG'][fun1] = item[4].split(";")
                        except:
                            funlist['COG'][fun1] = []
                    if fun1 in funlist['NOG']:
                        try:
                            nog_ids = item[5].split(";")
                            if nog_ids:
                                for nog in nog_ids:
                                    if nog not in funlist['NOG'][fun1]:
                                        funlist['NOG'][fun1].append(nog)
                        except:
                            pass
                    else:
                        try:
                            funlist['NOG'][fun1] = item[5].split(";")
                        except:
                            funlist['NOG'][fun1] = []
        for thekey in first:
            for g in func_type[thekey]:
                detail = func_decs[g]
                category = '[' + g + ']' + ' ' + detail
                try:
                    cog_list = list(set(funlist['COG'][g]))
                except:
                    cog_list = []
                try:
                    nog_list = list(set(funlist['NOG'][g]))
                except:
                    nog_list = []
                data = [
                    ('cog_id', cog_id),
                    # ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('type', thekey),
                    ('function_categories', category),
                    ('cog', len([x for x in cog_list if x ])),
                    ('nog', len([x for x in nog_list if x])),
                    ('cog_list', ';'.join(cog_list)),
                    ('nog_list', ';'.join(nog_list)),
                ]
                data = SON(data)
                if len(cog_list) + len(nog_list) != 0:
                    data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_cog_detail']
                if len(data_list) > 5000:
                    for i in range(0, len(data_list), 3000):
                        tmp_list = data_list[i: i+3000]
                        collection.insert_many(tmp_list)
                else:
                    collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入cog注释all出错：%s， %s" , variables=(r_cog_path, n_cog_path), code="52000553")
            else:
                self.bind_object.logger.info("导入cog注释all成功：%s, %s" % (r_cog_path, n_cog_path))

    @report_check
    def add_annotation_cog_table(self, cog_id, table_path, seq_type, anno_type):
        '''
        table_path:cog_table.xls
        '''
        if not isinstance(cog_id, ObjectId):
            if isinstance(cog_id, types.StringTypes):
                cog_id = ObjectId(cog_id)
            else:
                self.bind_object.set_error('cog_id必须为ObjectId对象或其对应的字符串！', code="52000554")
        if not os.path.exists(table_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(table_path), code="52000555")
        data_list = list()
        with open(table_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [
                    ('cog_id', cog_id),
                    # ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('query_name', line[0]),
                    ('query_length', line[1]),
                    ('hsp_start', line[2]),
                    ('hsp_end', line[3]),
                    ('hsp_strand', line[4]),
                    ('hit_name', line[5]),
                    ('hit_description', line[6]),
                    ('hit_length', line[7]),
                    ('hit_start', line[8]),
                    ('hit_end', line[9]),
                    ('group', line[10]),
                    ('group_description', line[11]),
                    ('group_categories', line[12]),
                    ('region_start', line[13]),
                    ('region_end', line[14]),
                    ('region_coverage', line[15]),
                    ('region_identities', line[16]),
                    ('region_positives', line[17])
                ]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_annotation_cog_table']
            if len(data_list) > 5000:
                for i in range(0, len(data_list), 3000):
                    tmp_list = data_list[i: i+3000]
                    collection.insert_many(tmp_list)
            else:
                collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入cog注释table信息：%s出错!" , variables=(table_path), code="52000556")
        else:
            self.bind_object.logger.info("导入cog注释table信息：%s成功!" % (table_path))

    @report_check
    def add_annotation_go(self, name=None, params=None, result_dir=None):
        """
        go注释导表函数
        """
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationGo_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'version': self.version,
            'result_dir': result_dir,
            'status': 'start',
            'desc': 'go注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        }
        collection = self.db['sg_annotation_go']
        go_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add ref_annotation_go!")
        return go_id

    def add_annotation_go_detail(self, go_id, seq_type, anno_type, level, go_path):
        """
        go_path: go1234level_statistics.xls/go123level_statistics.xls/go12level_statistics.xls
        """
        if not isinstance(go_id, ObjectId):
            if isinstance(go_id, types.StringTypes):
                go_id = ObjectId(go_id)
            else:
                self.bind_object.set_error('go_id必须为ObjectId对象或其对应的字符串！', code="52000557")
        if not os.path.exists(go_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(go_path), code="52000558")
        data_list = list()
        with open(go_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [
                    ('go_id', go_id),
                    # ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('level', level),
                    ('goterm', line[0]),
                    ('goterm_2', line[1]),
                    ('goid_2', line[2]),
                    ('seq_number', int(line[-3])),
                    ('percent', round(float(line[-2]), 4)),
                    #('seq_list', line[-1])
                ]
                if level == 2:
                    data.append(('seq_list', line[-1]))
                if level >= 3:
                    data.append(('goterm_3', line[3]))
                    data.append(('goid_3', line[4]))
                if level == 4:
                    data.append(('goterm_4', line[5]))
                    data.append(('goid_4', line[6]))
                data = SON(data)
                data_list.append(data)
            if data_list:
                try:
                    collection = self.db['sg_annotation_go_detail']
                    if len(data_list) > 5000:
                        for i in range(0, len(data_list), 3000):
                            tmp_list = data_list[i: i+3000]
                            collection.insert_many(tmp_list)
                    else:
                        collection.insert_many(data_list)
                except Exception, e:
                    self.bind_object.set_error("导入go注释信息：%s出错!" , variables=(go_path), code="52000559")
                else:
                    self.bind_object.logger.info("导入go注释信息：%s成功!" % (go_path))

    @report_check
    def add_annotation_go_graph(self, go_id, seq_type, anno_type, level, go_path):
        """
        go_path: go1234level_statistics.xls/go123level_statistics.xls/go12level_statistics.xls
        """
        if not isinstance(go_id, ObjectId):
            if isinstance(go_id, types.StringTypes):
                go_id = ObjectId(go_id)
            else:
                self.bind_object.set_error('go_id必须为ObjectId对象或其对应的字符串！', code="52000560")
        if not os.path.exists(go_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(go_path), code="52000561")
        data_list = list()
        with open(go_path, 'r') as f:
            lines = f.readlines()
            term = {}
            term_list = []
            for i in range(1, len(lines)):
                line = lines[i].strip().split('\t')
                if level == 2:
                    term_type = line[0]
                    go_term = line[1]
                    if go_term not in term:
                        term[go_term] = []
                        term[go_term].append(i)
                    else:
                        term[go_term].append(i)
                if level == 3:
                    term_type = line[0]
                    go_term = line[3]
                    if go_term not in term:
                        term[go_term] = []
                        term[go_term].append(i)
                    else:
                        term[go_term].append(i)
                if level == 4:
                    term_type = line[0]
                    go_term = line[5]
                    if go_term not in term:
                        term[go_term] = []
                        term[go_term].append(i)
                    else:
                        term[go_term].append(i)
        with open(go_path, 'r') as f:
            lines = f.readlines()
            for item in term:
                seq_list = []
                for j in term[item]:
                    line = lines[j].strip().split("\t")
                    term_type = line[0]
                    for seq in line[-1].split(";"):
                        if seq not in seq_list:
                            seq_list.append(seq)
                data = [
                    ('go_id', go_id),
                    # ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('level', level),
                    ('term_type', term_type),
                    ('go_term', item),
                    ('seq_number', len(seq_list)),
                    ('percent', line[-2]),
                    #('seq_list', seq_list)
                ]
                data = SON(data)
                data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_go_graph']
                if len(data_list) > 5000:
                    for i in range(0, len(data_list), 3000):
                        tmp_list = data_list[i: i+3000]
                        collection.insert_many(tmp_list)
                else:
                    collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入go注释画图信息出错：%s" , variables=(go_path), code="52000562")
            else:
                self.bind_object.logger.info("导入go注释画图信息成功：%s" % (go_path))

    @report_check
    def add_annotation_go_level(self, go_id, seq_type, anno_type, level, level_path):
        """
        level_path: go2level.xls
        """
        if not isinstance(go_id, ObjectId):
            if isinstance(go_id, types.StringTypes):
                go_id = ObjectId(go_id)
            else:
                self.bind_object.set_error('go_id必须为ObjectId对象或其对应的字符串！', code="52000563")
        if not os.path.exists(level_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(level_path), code="52000564")
        data_list = list()
        with open(level_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [
                    ('go_id', go_id),
                    # ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('level', level),
                    ('term_type', line[1]),
                    ('parent_name', line[0]),
                    ('num', int(line[3])),
                    ('percent', round(float(line[4]), 4)),
                    ('go', line[2]),
                    #('seq_list', line[5]),
                ]
                data = SON(data)
                data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_go_level']
                if len(data_list) > 5000:
                    for i in range(0, len(data_list), 3000):
                        tmp_list = data_list[i: i+3000]
                        collection.insert_many(tmp_list)
                else:
                    collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入go注释第二层级信息：%s出错!" , variables=(level_path), code="52000565")
            else:
                self.bind_object.logger.info("导入go注释第二层级信息：%s成功!" % (level_path))

    @report_check
    def add_annotation_go_list(self, go_id, seq_type, anno_type, gos_path):
        """
        gos_path: go.list
        """
        if not isinstance(go_id, ObjectId):
            if isinstance(go_id, types.StringTypes):
                go_id = ObjectId(go_id)
            else:
                self.bind_object.set_error('go_id须为ObjectId对象或其他对应的字符串！', code="52000566")
        if not os.path.exists(gos_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(gos_path), code="52000567")
        data_list = []
        with open(gos_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split('\t')
                data = [
                    ('go_id', go_id),
                    # ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('gene_id', line[0]),
                    ('gos_list', line[1]),
                ]
                data = SON(data)
                data_list.append(data)
            if data_list:
                try:
                    collection = self.db['sg_annotation_go_list']
                    if len(data_list) > 5000:
                        for i in range(0, len(data_list), 3000):
                            tmp_list = data_list[i: i+3000]
                            collection.insert_many(tmp_list)
                    else:
                        collection.insert_many(data_list)
                except Exception, e:
                    self.bind_object.set_error("导入gos_list注释信息：%s出错:%s" , variables=(gos_path), code="52000568")
                else:
                    self.bind_object.logger.info("导入gos_list注释信息：%s成功!" % (gos_path))

    def add_annotation_go_all(self, go_id, seq_type, anno_type, level, r_go_path, n_go_path):
        """
        r_go_path: go1234level_statistics.xls/go123level_statistics.xls/go12level_statistics.xls(ref)
        n_go_path: go1234level_statistics.xls/go123level_statistics.xls/go12level_statistics.xls(new)
        """
        if not isinstance(go_id, ObjectId):
            if isinstance(go_id, types.StringTypes):
                go_id = ObjectId(go_id)
            else:
                self.bind_object.set_error('go_id必须为ObjectId对象或其对应的字符串！', code="52000569")
        if not os.path.exists(r_go_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(r_go_path), code="52000570")
        if not os.path.exists(n_go_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(n_go_path), code="52000571")
        data_list1, data_list2, query_ids = list(), list(), list()
        funlist, termlist = {}, {}
        with open(r_go_path, 'r') as f, open(n_go_path, 'r') as n:
            lines = f.readlines()
            items = n.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                fun = line[0] + "|||" + line[1] + "|||" + line[2]
                term = line[0] + "|||" + line[1]
                if level == 3:
                    fun += "|||" + line[3] + "|||" + line[4]
                    term = line[0] + "|||" + line[3]
                if level == 4:
                    fun += "|||" + line[3] + "|||" + line[4]
                    fun += "|||" + line[5] + "|||" + line[6]
                    term = line[0] + "|||" + line[5]
                funlist[fun] = line[-1].split(";")
                if term not in termlist:
                    termlist[term] = line[-1].split(";")
                else:
                    for q in line[-1].split(";"):
                        if q not in termlist[term]:
                            termlist[term].append(q)
                for q in line[-1].split(";"):
                    if q not in query_ids:
                        query_ids.append(q)
            for item in items[1:]:
                item = item.strip().split('\t')
                fun = item[0] + "|||" + item[1] + "|||" + item[2]
                term = item[0] + "|||" + item[1]
                if level == 3:
                    fun += "|||" + item[3] + "|||" + item[4]
                    term = item[0] + "|||" + item[3]
                if level == 4:
                    fun += "|||" + item[3] + "|||" + item[4]
                    fun += "|||" + item[5] + "|||" + item[6]
                    term = item[0] + "|||" + item[5]
                if term not in termlist:
                    termlist[term] = item[-1].split(";")
                else:
                    for q in item[-1].split(";"):
                        if q not in termlist[term]:
                            termlist[term].append(q)
                if fun not in funlist:
                    funlist[fun] = item[-1].split(";")
                    for q in item[-1].split(";"):
                        if q not in query_ids:
                            query_ids.append(q)
                else:
                    for q in item[-1].split(";"):
                        if q not in query_ids:
                            query_ids.append(q)
                        if q not in funlist[fun]:
                            funlist[fun].append(q)
        for term in termlist:
            terms = term.split("|||")
            seq_list = termlist[term]
            percent = float(len(seq_list)) / len(query_ids)
            data = [
                ('go_id', go_id),
                # ('seq_type', seq_type),
                ('anno_type', anno_type),
                ('level', level),
                ('term_type', terms[0]),
                ('go_term', terms[1]),
                ('seq_number', len(seq_list)),
                ('percent', round(percent, 4))
                #('seq_list', ";".join(seq_list))
            ]
            data = SON(data)
            data_list1.append(data)
        if data_list1:
            try:
                collection = self.db['sg_annotation_go_graph']
                if len(data_list1) > 5000:
                    for i in range(0, len(data_list1), 3000):
                        tmp_list = data_list1[i: i+3000]
                        collection.insert_many(tmp_list)
                else:
                    collection.insert_many(data_list1)
            except Exception, e:
                self.bind_object.set_error("导入go注释画图all信息出错：%s, %s" % (r_go_path, n_go_path))
                print "导入go注释画图all出错：%s、%s" % (r_go_path, n_go_path)
            else:
                self.bind_object.logger.info("导入go注释画图all信息成功：%s, %s" % (r_go_path, n_go_path))
        for fun in funlist:
            terms = fun.split("|||")
            data = [
                ('go_id', go_id),
                # ('seq_type', seq_type),
                ('anno_type', anno_type),
                ('level', level),
                ('goterm', terms[0]),
                ('goterm_2', terms[1]),
                ('goid_2', terms[2])
            ]
            if level >= 3:
                data.append(('goterm_3', terms[3]))
                data.append(('goid_3', terms[4]))
            if level == 4:
                data.append(('goterm_4', terms[5]))
                data.append(('goid_4', terms[6]))
            seq_list = funlist[fun]
            percent = float(len(funlist[fun])) / len(query_ids)
            data.append(('seq_number', len(seq_list)))
            data.append(('percent', round(percent, 4)))
            data.append(('seq_list', ";".join(seq_list)))
            data = SON(data)
            data_list2.append(data)
        if data_list2:
            try:
                collection = self.db['sg_annotation_go_detail']
                if len(data_list2) > 5000:
                    for i in range(0, len(data_list2), 3000):
                        tmp_list = data_list2[i: i+3000]
                        collection.insert_many(tmp_list)
                else:
                    collection.insert_many(data_list2)
            except Exception, e:
                self.bind_object.set_error("导入go注释all出错：%s、%s" , variables=(r_go_path, n_go_path), code="52000572")
            else:
                self.bind_object.logger.info("导入go注释all信息成功：%s, %s" % (r_go_path, n_go_path))

    @report_check
    def add_annotation_kegg(self, name=None, params=None, result_dir=None):
        """
        kegg注释导表函数
        """
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationKegg_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'version': self.version,
            'result_dir': result_dir,
            'status': 'start',
            'desc': 'kegg注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'categories': ["M","CP","EIP","GIP","OS","HD","DD"]
        }
        collection = self.db['sg_annotation_kegg']
        kegg_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add ref_annotation_kegg!")
        return kegg_id

    @report_check
    def add_annotation_kegg_categories(self, kegg_id, seq_type, anno_type, categories_path):
        """
        categories_path:kegg_layer.xls
        """
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                self.bind_object.set_error('kegg_id必须为ObjectId对象或其对应的字符串！', code="52000573")
        if not os.path.exists(categories_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(categories_path), code="52000574")
        data_list = list()
        first_type = []
        with open(categories_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split('\t')
                type_abr = ''.join([x[0] for x in line[0].split(' ')])
                first_type.append(type_abr)
                data = [
                    ('kegg_id', kegg_id),
                    # ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('first_category', line[0]),
                    ('second_category', line[1]),
                    ('num', int(line[2])),
                    ('seq_list', line[3]),
                ]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_annotation_kegg_categories']
            if len(data_list) > 5000:
                for i in range(0, len(data_list), 3000):
                    tmp_list = data_list[i: i+3000]
                    collection.insert_many(tmp_list)
            else:
                collection.insert_many(data_list)
            self.update_db_record('sg_annotation_kegg', kegg_id, categories=list(set(first_type)))
        except Exception, e:
            self.bind_object.set_error("导入kegg注释分类信息：%s出错!" , variables=(categories_path), code="52000575")
        else:
            self.bind_object.logger.info("导入kegg注释分类信息：%s 成功!" % categories_path)

    @report_check
    def add_annotation_kegg_level(self, kegg_id, seq_type, anno_type, level_path, png_dir):
        """
        level_path: pathway_table.xls
        """
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                self.bind_object.set_error('kegg_id必须为ObjectId对象或其对应的字符串！', code="52000576")
        if not os.path.exists(level_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(level_path), code="52000577")
        if not os.path.exists(png_dir):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(png_dir), code="52000578")
        data_list = []
        with open(level_path, 'rb') as r:
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                fs = gridfs.GridFS(self.db)
                pid = re.sub('path:', '', line[0])
                # pdfid = fs.put(open(png_dir + '/' + pid + '.pdf', 'rb'))
                # graph_png_id = fs.put(open(png_dir + '/' + pid + '.png', 'rb'))
                insert_data = {
                    'kegg_id': kegg_id,
                    # 'seq_type': seq_type,
                    'anno_type': anno_type,
                    'pathway_id': line[0],
                    'first_category': line[1],
                    'second_category': line[2],
                    'pathway_definition': line[3],
                    'number_of_seqs': int(line[4]),
                    'seq_list': line[5],
                    # 'graph_id': pdfid,
                    # 'graph_png_id': graph_png_id,
                    'hyperlink': line[-1]
                }
                data_list.append(insert_data)
        if data_list:
            try:
                collection = self.db['sg_annotation_kegg_level']
                if len(data_list) > 5000:
                    for i in range(0, len(data_list), 3000):
                        tmp_list = data_list[i: i+3000]
                        collection.insert_many(tmp_list)
                else:
                    collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入kegg注释层级信息：%s、%s出错!" , variables=(level_path, png_dir), code="52000579")
            else:
                self.bind_object.logger.info("导入kegg注释层级信息：%s、%s 成功!" % (level_path, png_dir))

    @report_check
    def add_annotation_kegg_table(self, kegg_id, seq_type, anno_type, table_path):
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                self.bind_object.set_error('kegg_id必须为ObjectId对象或其对应的字符串！', code="52000580")
        if not os.path.exists(table_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(table_path), code="52000581")
        with open(table_path, 'rb') as r:
            data_list = []
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                insert_data = {
                    'kegg_id': kegg_id,
                    # 'seq_type': seq_type,
                    'anno_type': anno_type,
                    'transcript_id': line[0],
                    'ko_id': line[1],
                    'ko_name': line[2],
                    'hyperlink': line[3],
                    'paths': line[4],
                }
                data_list.append(insert_data)
        try:
            collection = self.db['sg_annotation_kegg_table']
            if len(data_list) > 5000:
                for i in range(0, len(data_list), 3000):
                    tmp_list = data_list[i: i+3000]
                    collection.insert_many(tmp_list)
            else:
                collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入kegg注释table信息：%s出错!" , variables=(table_path), code="52000582")
        else:
            self.bind_object.logger.info("导入kegg注释table信息：%s成功!" % (table_path))

    @report_check
    def add_annotation_kegg_categories_all(self, kegg_id, seq_type, anno_type, r_cate_path, n_cate_path):
        """
        r_cate_path:kegg_layer.xls(ref)
        n_cate_path:kegg_layer.xls(new)
        """
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                self.bind_object.set_error('kegg_id必须为ObjectId对象或其对应的字符串！', code="52000583")
        if not os.path.exists(r_cate_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(r_cate_path), code="52000584")
        if not os.path.exists(n_cate_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(n_cate_path), code="52000585")
        data_list = list()
        cate = {}
        with open(r_cate_path, 'r') as f, open(n_cate_path, "r") as n:
            lines = f.readlines()
            items = n.readlines()
            for line in lines:
                line = line.strip().split('\t')
                if line[0] == "Metabolism" and line[1] == "Global and overview maps":
                    pass
                else:
                    cate_d = line[0] + "|||" + line[1]
                    cate[cate_d] = line[3].split(";")
            for item in items:
                item = item.strip().split('\t')
                if item[0] == "Metabolism" and item[1] == "Global and overview maps":
                    pass
                else:
                    cate_d = item[0] + "|||" + item[1]
                    if cate_d not in cate:
                        cate[cate_d] = item[3].split(";")
                    else:
                        ids = item[3].split(";")
                        for q in ids:
                            if q not in cate[cate_d]:
                                cate[cate_d].append(q)
        for f in ["Metabolism", "Genetic Information Processing", "Environmental Information Processing", "Cellular Processes", "Organismal Systems", "Human Diseases", "Drug Development"]:
            for c in cate:
                ca = c.split("|||")
                first_category = ca[0]
                second_category = ca[1]
                if f == first_category:
                    num = len(cate[c])
                    seq_list = ";".join(cate[c])
                    data = [
                        ('kegg_id', kegg_id),
                        # ('seq_type', seq_type),
                        ('anno_type', anno_type),
                        ('first_category', first_category),
                        ('second_category', second_category),
                        ('num', num),
                        ('seq_list', seq_list)
                    ]
                    data = SON(data)
                    data_list.append(data)
        try:
            collection = self.db['sg_annotation_kegg_categories']
            if len(data_list) > 5000:
                for i in range(0, len(data_list), 3000):
                    tmp_list = data_list[i: i+3000]
                    collection.insert_many(tmp_list)
            else:
                collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入kegg注释分类alls出错：%s, %s" , variables=(r_cate_path, n_cate_path), code="52000586")
        else:
            self.bind_object.logger.info("导入kegg注释分类all成功：%s, %s" % (r_cate_path, n_cate_path))

    def get_pic(self, path, kos_path, png_path):
        """
        画通路图
        """
        fs = gridfs.GridFS(self.mongodb)
        pid = re.sub("map", "ko", path)
        with open("pathway.kgml", "w+") as k, open("pathway.png", "w+") as p:
            result = self.png_coll.find_one({"pathway_id": pid})
            if result:
                kgml_id = result['pathway_ko_kgml']
                png_id = result['pathway_map_png']
                k.write(fs.get(kgml_id).read())
                p.write(fs.get(png_id).read())
        cmd = "{} {} {} {} {} {} {}".format(self.r_path, self.map_path, path, kos_path, png_path, "pathway.kgml", "pathway.png")
        try:
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError:
            print "{}画图出错".format(path)
            os.system("cp {} {}".format("pathway.png", png_path))

    @report_check
    def add_annotation_kegg_level_all(self, kegg_id, seq_type, anno_type, r_level_path, n_level_path):
        """
        r_level_path: pathway_table.xls(ref)
        n_level_path: pathway_table.xls(new)
        """
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                self.bind_object.set_error('kegg_id必须为ObjectId对象或其对应的字符串！', code="52000587")
        if not os.path.exists(r_level_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(r_level_path), code="52000588")
        if not os.path.exists(n_level_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(n_level_path), code="52000589")
        data_list = []
        path_def = {}
        r_path_list = {}
        n_path_list = {}
        fs = gridfs.GridFS(self.mongodb)
        with open(r_level_path, "rb") as r, open(n_level_path, "rb") as n:
            lines = r.readlines()
            items = n.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                if line[1] == "Metabolism" and line[2] == "Global and overview maps":
                    pass
                else:
                    path = line[0] + "|||" + line[1] + "|||" + line[2] + "|||" + line[3]
                    sp = 'http://www.genome.jp/dbget-bin/show_pathway?' + line[0]
                    k_cols = line[7].split(sp)
                    k_ids = k_cols[1].split("%09yellow")
                    k_list = []
                    for k in k_ids:
                        if k.startswith('/'):
                            k_id = k.split('/')[1]
                            k_list.append(k_id)
                    seqlist = line[5].split(";")
                    path_def[line[0]] = path
                    r_path_list[line[0]] = []
                    r_path_list[line[0]].append(seqlist)
                    r_path_list[line[0]].append(k_list)
            for item in items[1:]:
                item = item.strip().split("\t")
                if item[1] == "Metabolism" and item[2] == "Global and overview maps":
                    pass
                else:
                    path = item[0] + "|||" + item[1] + "|||" + item[2] + "|||" + item[3]
                    sp = 'http://www.genome.jp/dbget-bin/show_pathway?' + item[0]
                    k_cols = item[7].split(sp)
                    k_ids = k_cols[1].split("%09green")
                    k_list = []
                    for k in k_ids:
                        if k.startswith('/'):
                            k_id = k.split('/')[1]
                            k_list.append(k_id)
                    seqlist = item[5].split(";")
                    n_path_list[item[0]] = []
                    n_path_list[item[0]].append(seqlist)
                    n_path_list[item[0]].append(k_list)
                    if item[0] not in path_def:
                        path_def[item[0]] = path
        for map_id in path_def:
            link = []
            r_kos, n_kos, b_kos = [], [], []
            ref, new, both = [], [], []
            paths = path_def[map_id].split("|||")
            first_category = paths[1]
            second_category = paths[2]
            pathway_definition = paths[3]
            try:
                seq_list = r_path_list[map_id][0]
                r_ko = r_path_list[map_id][1]
            except:
                seq_list = []
                r_ko = []
            try:
                for q in n_path_list[map_id][0]:
                    seq_list.append(q)
                n_ko = n_path_list[map_id][1]
            except:
                n_ko = []
            seq_list = list(set(seq_list))
            for r_c in r_ko:
                r = r_c.split("%09tomato")[0]
                if r_c not in n_ko:
                    r_id = r + '%09' + 'yellow'
                    link.append(r_id)
                    r_kos.append(r)
                else:
                    b_id = r + '%09' + 'tomato'
                    link.append(b_id)
                    b_kos.append(r)
            for n_c in n_ko:
                n = n_c.split("%09tomato")[0]
                if n_c not in r_ko:
                    n_id = n + '%09' + 'green'
                    link.append(n_id)
                    n_kos.append(n)
                else:
                    b_id = n + '%09' + 'tomato'
                    link.append(b_id)
                    b_kos.append(n)
            link = list(set(link))
            b_kos = list(set(b_kos))
            link = 'http://www.genome.jp/kegg-bin/show_pathway?' + map_id + '/' + '/'.join(link)
            png_path = os.getcwd() + '/' + map_id + ".png"
            pdf_path = os.getcwd() + '/' + map_id + ".pdf"
            kos_path = os.path.join(os.getcwd(), "KOs.txt")
            with open(kos_path, "w") as w:
                w.write("#KO\tbg\tfg\n")
                for k in n_kos:
                    w.write(k + "\t" + "#00CD00" + "\t" + "NA" + "\n")
                for k in r_kos:
                    w.write(k + "\t" + "#FFFF00" + "\t" + "NA" + "\n")
                for k in b_kos:
                    w.write(k + "\t" + "#FFFF00,#00CD00" + "\t" + "NA" + "\n")
            self.get_pic(map_id, kos_path, png_path)
            cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' + png_path + ' ' + pdf_path
            try:
                subprocess.check_output(cmd, shell=True)
            except subprocess.CalledProcessError:
                print '图片格式pdf转png出错'
            pdfid = fs.put(open(pdf_path, 'rb'))
            graph_png_id = fs.put(open(png_path, 'rb'))
            insert_data = {
                'kegg_id': kegg_id,
                # 'seq_type': seq_type,
                'anno_type': anno_type,
                'pathway_id': map_id,
                'first_category': first_category,
                'second_category': second_category,
                'pathway_definition': pathway_definition,
                'number_of_seqs': len(seq_list),
                'seq_list': ";".join(seq_list),
                'graph_id': pdfid,
                "hyperlink": link,
                'graph_png_id': graph_png_id
            }
            data_list.append(insert_data)
            # os.remove(pdf)
        try:
            collection = self.db['sg_annotation_kegg_level']
            if len(data_list) > 5000:
                for i in range(0, len(data_list), 3000):
                    tmp_list = data_list[i: i+3000]
                    collection.insert_many(tmp_list)
            else:
                collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入kegg注释层级all信息出错：%s、%s" , variables=(r_level_path, n_level_path), code="52000590")
        else:
            self.bind_object.logger.info("导入kegg注释层级all信息成功：%s、%s" % (level_path, png_dir))

    @report_check
    def add_annotation_kegg_pic(self, kegg_id, seq_type, anno_type, level_path, png_dir):
        """
        level_path: pathway_table.xls
        """
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                self.bind_object.set_error('kegg_id必须为ObjectId对象或其对应的字符串！', code="52000591")
        if not os.path.exists(level_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(level_path), code="52000592")
        if not os.path.exists(png_dir):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(png_dir), code="52000593")
        data_list = []
        with open(level_path, 'rb') as r:
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                pid = re.sub('path:', '', line[0])
                if os.path.exists(png_dir + '/' + line[0] + '.html.mark'):
                    with open(png_dir + '/' + line[0] + '.html.mark', 'r') as mark_f:
                        for line_mark in mark_f.readlines():
                            # print len(line_mark.strip("\n").split("\t"))
                            if len(line_mark.strip("\n").split("\t")) == 8:
                                [png, shape, bg_color, fg_color, coords, title, kos, href] = line_mark.strip("\n").split("\t")
                                title = title.replace("\\n", "\n")
                            else:
                                continue

                            insert_data = {
                                'kegg_id': kegg_id,
                                'seq_type': seq_type,
                                'anno_type': anno_type,
                                'pathway_id': line[0],
                                'shape': shape,
                                'bg_colors': bg_color,
                                'fg_colors': fg_color,
                                'coords': coords,
                                'href': href,
                                'kos': kos,
                                'title': title
                            }

                            if bg_color != "" and len(bg_color.split(",")) > 0:
                                insert_data.update({'bg_type': len(bg_color.split(","))})
                            if fg_color != "" and len(fg_color.split(",")) > 0:
                                insert_data.update({'fg_type': len(fg_color.split(","))})
                            data_list.append(insert_data)
                else:
                    self.bind_object.logger.info("kegg 图片{} 不存在html标记!".format(line[0]))

        if data_list:
            try:
                collection = self.db['sg_annotation_kegg_pic']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入kegg注释图片信息：%s、%s出错!" , variables=(level_path, png_dir), code="52000594")
            else:
                self.bind_object.logger.info("导入kegg注释图片信息：%s、%s 成功!" % (level_path, png_dir))


    @report_check
    def add_annotation_query(self, name=None, params=None, stat_id=None, result_dir=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        if not isinstance(stat_id, ObjectId):
            print stat_id
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="52000595")
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationQuery' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'version': self.version,
            'result_dir': result_dir,
            'status': 'start',
            'desc': '注释查询主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'stat_id': stat_id
        }
        collection = self.db['sg_annotation_query']
        query_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add ref_annotation_query!")
        return query_id

    def get_kegg_map2class(self):
        '''
        返回kegg mapid 与 classI classII 对应关系字典
        key: map00010
        value: (classi_name, classii_name)
        '''
        map2class = dict()
        with open(self.kegg_json, "rb") as f:
            root = json.load(f)
        classI = root['children']

        for class1 in classI:
            class1_name = class1['name']
            for class2 in class1['children']:
                class2_name = class2['name']
                paths = class2['children']
                for path in paths:
                    path_name = "map" + str(path['name']).split(" ")[0]
                    map2class[path_name] = (class1_name, class2_name)
        return map2class


    @report_check
    def add_annotation_query_denovo_detail(self, query_id, query_path, anno_type):
        if not isinstance(query_id, ObjectId):
            if isinstance(query_id, types.StringTypes):
                query_id = ObjectId(query_id)
            else:
                self.bind_object.set_error('query_id必须为ObjectId对象或其对应的字符串！', code="52000596")
        if not os.path.exists(query_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(query_path), code="52000597")
        data_list = []
        map2class = self.get_kegg_map2class()
        with open(query_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip("\n").split('\t')
                data = [
                    ('query_id', query_id),
                    ('anno_type', anno_type),
                    ('transcript_id', line[0]),
                    ('gene_id', line[1]),
                    ('is_gene', self.trans_isgene[line[0]]),
                ]
                try:
                    data.append(('length', line[3]))
                except:
                    data.append(('length', ""))
                try:
                    data.append(('cog', line[5]))
                    data.append(('cog_description', line[6]))
                except:
                    data.append(('cog', ""))
                    data.append(('cog_description', ""))
                '''
                try:
                    data.append(('nog', line[4]))
                    data.append(('nog_description', line[6]))
                except:
                    data.append(('nog', ""))
                    data.append(('nog_description', ""))
                '''
                try:
                    data.append(('ko_id', line[7]))
                except:
                    data.append(('ko_id', ""))
                try:
                    data.append(('ko_name', line[8]))
                except:
                    data.append(('ko_name', ""))

                try:
                    data.append(('pathways', line[9]))
                    paths = [pathway.split("(")[0] for pathway in line[9].split("; ")]
                    if line[9] != "":
                        pathway_class1 = [map2class[path][0] for path in paths if path in map2class]
                        pathway_class2 = [map2class[path][1] for path in paths if path in map2class]
                    data.append(('pathways_class1', ";".join(pathway_class1)))
                    data.append(('pathways_class2', ";".join(pathway_class2)))
                except:
                    data.append(('pathways', ""))
                try:
                    data.append(('pfam', line[10]))
                except:
                    data.append(('pfam', ""))
                try:
                    data.append(('go', line[11]))
                except:
                    data.append(('go', ""))
                try:
                    data.append(('nr', line[12]))
                except:
                    data.append(('nr', ""))
                try:
                    data.append(('swissprot', line[13]))
                except:
                    data.append(('swissprot', ""))
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_annotation_query_detail']
            if len(data_list) > 5000:
                for i in range(0, len(data_list), 3000):
                    tmp_list = data_list[i: i+3000]
                    collection.insert_many(tmp_list)
            else:
                collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入转录本注释统计信息：%s出错!" , variables=(query_path), code="52000598")
        else:
            self.bind_object.logger.info("导入转录本注释统计信息：%s成功!" % (query_path))


    @report_check
    def add_annotation_query_detail(self, query_id, query_path, anno_type):
        if not isinstance(query_id, ObjectId):
            if isinstance(query_id, types.StringTypes):
                query_id = ObjectId(query_id)
            else:
                self.bind_object.set_error('query_id必须为ObjectId对象或其对应的字符串！', code="52000599")
        if not os.path.exists(query_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(query_path), code="520005100")
        data_list = []
        with open(query_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [
                    ('query_id', query_id),
                    ('anno_type', anno_type),
                    ('transcript_id', line[0]),
                    ('gene_id', line[1]),
                ]
                try:
                    data.append(('gene_name', line[2]))
                except:
                    data.append(('gene_name', None))
                try:
                    data.append(('length', line[3]))
                except:
                    data.append(('length', None))
                try:
                    data.append(('cog', line[4]))
                    data.append(('cog_description', line[6]))
                except:
                    data.append(('cog', None))
                    data.append(('cog_description', None))
                try:
                    data.append(('nog', line[5]))
                    data.append(('nog_description', line[7]))
                except:
                    data.append(('nog', None))
                    data.append(('nog_description', None))
                try:
                    data.append(('ko_id', line[8]))
                except:
                    data.append(('ko_id', None))
                try:
                    data.append(('ko_name', line[9]))
                except:
                    data.append(('ko_name', None))

                try:
                    data.append(('pathways', line[10]))
                except:
                    data.append(('pathways', None))
                try:
                    data.append(('pfam', line[11]))
                except:
                    data.append(('pfam', None))
                try:
                    data.append(('go', line[12]))
                except:
                    data.append(('go', None))
                try:
                    data.append(('nr', line[13]))
                except:
                    data.append(('nr', None))
                try:
                    data.append(('swissprot', line[14]))
                except:
                    data.append(('swissprot', None))
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_annotation_query_detail']
            if len(data_list) > 5000:
                for i in range(0, len(data_list), 3000):
                    tmp_list = data_list[i: i+3000]
                    collection.insert_many(tmp_list)
            else:
                collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入转录本注释统计信息：%s出错!" , variables=(query_path), code="520005101")
        else:
            self.bind_object.logger.info("导入转录本注释统计信息：%s成功!" % (query_path))

    @report_check
    def add_annotation_gene_query_detail(self, query_id, query_path, anno_type):
        if not isinstance(query_id, ObjectId):
            if isinstance(query_id, types.StringTypes):
                query_id = ObjectId(query_id)
            else:
                self.bind_object.set_error('query_id必须为ObjectId对象或其对应的字符串！', code="520005102")
        if not os.path.exists(query_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(query_path), code="520005103")
        data_list = []
        with open(query_path, 'r') as f:
            lines = f.readlines()
            for j in range(1, len(lines)):
                line = lines[j].strip().split('\t')
                data = [
                    ('query_id', query_id),
                    ('gene_id', line[0]),
                    ('anno_type', anno_type),
                ]
                try:
                    data.append(('gene_name', line[1]))
                except:
                    data.append(('gene_name', None))
                try:
                    data.append(('cog', line[2]))
                    data.append(('cog_description', line[4]))
                except:
                    data.append(('cog', None))
                    data.append(('cog_description', None))
                try:
                    data.append(('nog', line[3]))
                    data.append(('nog_description', line[5]))
                except:
                    data.append(('nog', None))
                    data.append(('nog_description', None))
                try:
                    data.append(('ko_id', line[6]))
                except:
                    data.append(('ko_id', None))
                try:
                    data.append(('ko_name', line[7]))
                except:
                    data.append(('ko_name', None))
                try:
                    data.append(('pathways', line[8]))
                except:
                    data.append(('pathways', None))
                try:
                    data.append(('pfam', line[9]))
                except:
                    data.append(('pfam', None))
                try:
                    data.append(('go', line[10]))
                except:
                    data.append(('go', None))
                try:
                    data.append(('nr', line[11]))
                except:
                    data.append(('nr', None))
                try:
                    data.append(('swissprot', line[12]))
                except:
                    data.append(('swissprot', None))
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_annotation_query_detail']
            if len(data_list) > 5000:
                for i in range(0, len(data_list), 3000):
                    tmp_list = data_list[i: i+3000]
                    collection.insert_many(tmp_list)
            else:
                collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入基因注释查询信息：%s出错!" , variables=(query_path), code="520005104")
        else:
            self.bind_object.logger.info("导入基因注释统计信息：%s成功!" % (query_path))

    @report_check
    def add_annotation_gene_query_denovo_detail(self, query_id, query_path, anno_type):
        if not isinstance(query_id, ObjectId):
            if isinstance(query_id, types.StringTypes):
                query_id = ObjectId(query_id)
            else:
                self.bind_object.set_error('query_id必须为ObjectId对象或其对应的字符串！', code="520005105")
        if not os.path.exists(query_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(query_path), code="520005106")
        data_list = []
        with open(query_path, 'r') as f:
            lines = f.readlines()
            for j in range(1, len(lines)):
                line = lines[j].strip().split('\t')
                data = [
                    ('query_id', query_id),
                    ('gene_id', line[0]),
                    ('anno_type', anno_type),
                ]
                try:
                    data.append(('cog', line[1]))
                    data.append(('cog_description', line[3]))
                except:
                    data.append(('cog', None))
                    data.append(('cog_description', None))
                try:
                    data.append(('nog', line[2]))
                    data.append(('nog_description', line[4]))
                except:
                    data.append(('nog', None))
                    data.append(('nog_description', None))
                try:
                    data.append(('ko_id', line[5]))
                except:
                    data.append(('ko_id', None))
                try:
                    data.append(('ko_name', line[6]))
                except:
                    data.append(('ko_name', None))
                try:
                    data.append(('pathways', line[7]))
                except:
                    data.append(('pathways', None))
                try:
                    data.append(('pfam', line[8]))
                except:
                    data.append(('pfam', None))
                try:
                    data.append(('go', line[9]))
                except:
                    data.append(('go', None))
                try:
                    data.append(('nr', line[10]))
                except:
                    data.append(('nr', None))
                try:
                    data.append(('swissprot', line[11]))
                except:
                    data.append(('swissprot', None))
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_annotation_query_detail']
            if len(data_list) > 5000:
                for i in range(0, len(data_list), 3000):
                    tmp_list = data_list[i: i+3000]
                    collection.insert_many(tmp_list)
            else:
                collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入基因注释查询信息：%s出错!" , variables=(query_path), code="520005107")
        else:
            self.bind_object.logger.info("导入基因注释统计信息：%s成功!" % (query_path))

class TestFunction(unittest.TestCase):
    """
    测试导表函数

    """
    def test_mongo(test):
        from mbio.workflows.denovo_rna_v2.denovo_test_api import DenovoTestApiWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            "id": "denovo_rna_v2_upgrade2",
            #+ str(random.randint(1,10000)),
            #"id": "denovo_rna_v2",
            "project_sn": "denovo_rna_v2_upgrade2",
            #+ str(random.randint(1,10000)),
            "type": "workflow",
            "name": "denovo_rna_v2.denovo_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = DenovoTestApiWorkflow(wsheet)

        test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_denovo_rna_v2/annot_out'
        trans2gene = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_denovo_rna_v2/annot_out/all_tran2gene.txt'
        gene_exp = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_denovo_rna_v2/annot_out/gene_exp'
        trans_exp = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_denovo_rna_v2/annot_out/trans_exp'

        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("denovo_rna_v2.denovo_annotation")
        params = {
            "submit_location": "annotationstat",
            "task_id": "denovo_rna_v2",
            "task_type": 2,
            "nr_evalue": 1e-3,
            "nr_similarity": 0,
            "nr_identity": 0,
            "swissprot_evalue":1e-3,
            "swissprot_similarity": 0,
            "swissprot_identity": 0,
            "cog_evalue": 1e-3,
            "cog_similarity": 0,
            "cog_identity": 0,
            "kegg_evalue": 1e-3,
            "kegg_similarity": 0,
            "kegg_identity": 0,
            "pfam_evalue": 1e-3,
        }
        wf.test_api.anno_type = 'origin'
        wf.test_api.run(test_dir, trans2gene,  params, taxon="Plant", version="v2", exp_level="T", gene_exp=gene_exp, trans_exp=trans_exp)

if __name__ == '__main__':
    unittest.main()
