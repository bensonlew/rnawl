# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import glob
import os
import re
import shutil

import numpy as np
import pandas as pd
from biocluster.config import Config
from bson.objectid import ObjectId
from mbio.packages.medical_transcriptome.copy_file import CopyFile
from mbio.packages.medical_transcriptome.gene_info_supple import GeneInfoSupple


# self.database = Config().get_mongo_client(mtype='medical_transcriptome')[Config().get_mongo_dbname('medical_transcriptome')]

class Upload(object):
    def __init__(self,project_type="medical_transcriptome",task_id=None,level =None):
        # self.database = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        self.task_id = task_id
        self.level = level

    def set_diff_pipline(self,idir=None,odir=None,species=None,level=None):
        if os.path.isdir(odir):
            shutil.rmtree(odir)
        else:
            os.mkdir(odir)
        add_name = GeneInfoSupple(task_id=self.task_id, level=self.level)
        # 03DiffExpress_ Cluster_Analysis
        # 获取基因集信息
        os.makedirs(os.path.join(odir,  '03DiffExpress_Geneset_Annotion'))
        diff_genesets = []
        for gt in os.listdir(os.path.join(idir)):
            if gt != "cluster":
                diff_genesets.append(gt)
        # GO_Annotion
        os.makedirs(os.path.join(odir,  '03DiffExpress_Geneset_Annotion', '01GO'))
        for geneset in diff_genesets:
            os.makedirs(os.path.join(odir,  '03DiffExpress_Geneset_Annotion',
                                     '01GO', geneset))
            # os.link(os.path.join(idir, geneset, 'diff_go_class', 'go_class_table.xls'),
            #         os.path.join(odir,  '03DiffExpress_Geneset_Annotion',
            #                      '01GO', geneset, "go_class_stat.xls"))
            #add by fwy 20210415 新增一列gene_name
            raw_path = os.path.join(idir, geneset, 'diff_go_class', 'go_class_table.xls')
            file_path = os.path.join(odir,  '03DiffExpress_Geneset_Annotion','01GO', geneset, "go_class_stat.xls")
            split = ";"
            add_columns = [geneset + " list"]
            # add_name.add_gene_name(file_path,split,"go",add_columns,)
            try:
                add_name.add_gene_name(raw_path, split, "go", add_columns, file_path)
            except:
                os.link(raw_path, file_path)


        # 02KEGG_Annotion
        os.makedirs(os.path.join(odir, '03DiffExpress_Geneset_Annotion', '02KEGG'))
        for geneset in diff_genesets:
            os.makedirs(
                os.path.join(odir,  '03DiffExpress_Geneset_Annotion', '02KEGG', geneset))
            # os.link(os.path.join(idir, geneset, 'diff_kegg_class', 'kegg_stat.xls'),
            #         os.path.join(odir,  '03DiffExpress_Geneset_Annotion',
            #                      '02KEGG', geneset, "kegg_class_stat.xls"))

            # add by fwy 20210415 新增一列gene_name
            raw_path = os.path.join(idir, geneset, 'diff_kegg_class', 'kegg_stat.xls')
            file_path = os.path.join(odir,  '03DiffExpress_Geneset_Annotion','02KEGG', geneset, "kegg_class_stat.xls")
            split = ";"
            kc = pd.read_table(raw_path, sep="\t", index_col=False)
            add_columns = [x for x in kc.columns if x.endswith("genes")]
            try:
                add_name.add_gene_name(raw_path, split, "kegg", add_columns, file_path)
            except:
                os.link(raw_path, file_path)
            os.link(os.path.join(idir, geneset, 'diff_kegg_class', 'pathways.tar.gz'),
                            os.path.join(odir, '03DiffExpress_Geneset_Annotion',
                                         '02KEGG', geneset,
                                         "kegg_class_pathways.tar.gz"))
            # shutil.copytree(os.path.join(idir, geneset, 'diff_kegg_class', 'pathways'),
            #                 os.path.join(odir,  '03DiffExpress_Geneset_Annotion',
            #                              '02KEGG', geneset,
            #                              "pathways"))

        if level == "G":
            # 03 Reactome _Annotion
            os.makedirs(os.path.join(odir,  '03DiffExpress_Geneset_Annotion', '03Reactome'))
            for geneset in diff_genesets:
                os.makedirs(os.path.join(odir,  '03DiffExpress_Geneset_Annotion',
                                         '03Reactome', geneset))
                # os.link(os.path.join(idir, geneset, 'diff_reactome_class',
                #                      'reactome_class.xls'),
                #         os.path.join(odir, '03DiffExpress_Geneset_Annotion',
                #                      '03Reactome', geneset, "Reactome_class_stat.xls"))
                # add by fwy 20210415 新增一列gene_name
                raw_path = os.path.join(idir, geneset, 'diff_reactome_class','reactome_class.xls')
                file_path = os.path.join(odir, '03DiffExpress_Geneset_Annotion', '03Reactome', geneset, "Reactome_class_stat.xls")
                split = ";"
                rc = pd.read_table(raw_path, sep="\t", index_col=False)
                add_columns = [x for x in rc.columns if x.endswith("genes")]
                try:
                    add_name.add_gene_name(raw_path, split, "reactome", add_columns, file_path)
                except:
                    os.link(raw_path, file_path)

                try:
                    raw_path = os.path.join(idir, geneset, 'diff_reactome_class','reactome_path.xls')
                    file_path = os.path.join(odir, '03DiffExpress_Geneset_Annotion', '03Reactome', geneset, "pathway_class_stat.xls")
                    os.link(raw_path, file_path)
                except Exception as e:
                    print(e)

                os.link(
                    os.path.join(idir, geneset, 'diff_reactome_class', 'svg.tar.gz'),
                    os.path.join(odir, '03DiffExpress_Geneset_Annotion',
                                 '03Reactome', geneset,
                                 "Reactome_pathways.tar.gz"))

                # shutil.copytree(
                #     os.path.join(idir, geneset, 'diff_reactome_class', 'svg'),
                #     os.path.join(odir, '03DiffExpress_Geneset_Annotion',
                #                  '03Reactome', geneset,
                #                  "Reactome_pathways"))
        # 04DO_Annotion
        if species == "Homo_sapiens" and level == "G":
            os.makedirs(os.path.join(odir, '03DiffExpress_Geneset_Annotion',
                                     '04DO'))
            for geneset in diff_genesets:
                os.makedirs(os.path.join(odir,  '03DiffExpress_Geneset_Annotion',
                                         '04DO', geneset))
                # os.link(os.path.join(idir, geneset, 'diff_do_class',
                #                      'do_level2_class.xls'),
                #         os.path.join(odir,  '03DiffExpress_Geneset_Annotion',
                #                      '04DO', geneset,
                #                      "DO_class_stat.xls"))
                # add by fwy 20210415 新增一列gene_name
                raw_path = os.path.join(idir, geneset, 'diff_do_class','do_level2_class.xls')
                file_path = os.path.join(odir,  '03DiffExpress_Geneset_Annotion','04DO', geneset,"DO_class_stat.xls")
                split = ";"
                dc = pd.read_table(raw_path, sep="\t", index_col=False)
                add_columns = [x for x in dc.columns if x.endswith("seqs")]
                try:
                    add_name.add_gene_name(raw_path, split, "do", add_columns, file_path)
                except:
                    os.link(raw_path, file_path)



        os.makedirs(os.path.join(odir,  '04DiffExpress_Geneset_Enrich'))
        # GO_Enrich
        os.makedirs(
            os.path.join(odir,  '04DiffExpress_Geneset_Enrich', '01GO'))
        for geneset in diff_genesets:
            os.makedirs(os.path.join(odir,  '04DiffExpress_Geneset_Enrich',
                                     '01GO', geneset))
            # os.link(os.path.join(idir, geneset, 'diff_go_enrich',
            #                      'go_enrich_geneset_list_gene.xls'),
            #         os.path.join(odir,  '04DiffExpress_Geneset_Enrich',
            #                      '01GO', geneset, "go_enrich_stat.xls"))

            # add by fwy 20210415 新增一列gene_name
            raw_path = os.path.join(idir, geneset, 'diff_go_enrich', 'go_enrich_geneset_list_gene.xls')
            file_path = os.path.join(odir,  '04DiffExpress_Geneset_Enrich','01GO', geneset, "go_enrich_stat.xls")
            split = ";"
            add_columns = ["seq_list"]
            try:
                add_name.add_gene_name(raw_path, split, "go", add_columns, file_path)
            except:
                os.link(raw_path, file_path)

        # 02KEGG_Enrich
        os.makedirs(
            os.path.join(odir, '04DiffExpress_Geneset_Enrich', '02KEGG'))
        for geneset in diff_genesets:
            os.makedirs(os.path.join(odir,  '04DiffExpress_Geneset_Enrich',
                                     '02KEGG', geneset))
            # os.link(os.path.join(idir, geneset, 'diff_kegg_enrich', 'enrich',
            #                      '{}_gene.list.DE.list.check.kegg_enrichment.xls'.format(geneset)),
            #         os.path.join(odir,  '04DiffExpress_Geneset_Enrich',
            #                      '02KEGG', geneset, "kegg_enrich_stat.xls"))

            # add by fwy 20210415 新增一列gene_name
            raw_path =os.path.join(idir, geneset, 'diff_kegg_enrich', 'enrich','{}_gene.list.DE.list.check.kegg_enrichment.xls'.format(geneset))
            file_path = os.path.join(odir,  '04DiffExpress_Geneset_Enrich','02KEGG', geneset, "kegg_enrich_stat.xls")
            split = "|"
            add_columns = ["Genes"]
            try:
                add_name.add_gene_name(raw_path, split, "kegg_enrich", add_columns, file_path)
            except:
                os.link(raw_path, file_path)

            # shutil.copytree(
            #     os.path.join(idir, geneset, 'diff_kegg_enrich', 'class', 'pathways'),
            #     os.path.join(odir,  '04DiffExpress_Geneset_Enrich',
            #                  '02KEGG', geneset,
            #                  "kegg_enrich_pathways"))
            os.link(
                os.path.join(idir, geneset, 'diff_kegg_enrich', 'class', 'pathways.tar.gz'),
                os.path.join(odir, '04DiffExpress_Geneset_Enrich',
                             '02KEGG', geneset,
                             "kegg_enrich_pathways.tar.gz"))
        # 03 Reactome_Enrich
        os.makedirs(
            os.path.join(odir, '04DiffExpress_Geneset_Enrich', '03Reactome'))
        for geneset in diff_genesets:
            os.makedirs(os.path.join(odir,  '04DiffExpress_Geneset_Enrich',
                                     '03Reactome', geneset))
            # os.link(os.path.join(idir, geneset, 'diff_reactome_enrich',
            #                      '{}_gene.list.reactome_enrichment.xls'.format(geneset)),
            #         os.path.join(odir,  '04DiffExpress_Geneset_Enrich',
            #                      '03Reactome', geneset, "Reactome_enrich_stat.xls"))

            # add by fwy 20210415 新增一列gene_name
            raw_path = os.path.join(idir, geneset, 'diff_reactome_enrich','{}_gene.list.reactome_enrichment.xls'.format(geneset))
            file_path =  os.path.join(odir,  '04DiffExpress_Geneset_Enrich','03Reactome', geneset, "Reactome_enrich_stat.xls")
            split = "|"
            add_columns = ["Genes"]
            try:
                add_name.add_gene_name(raw_path, split, "reactome", add_columns, file_path)
            except:
                if os.path.exists(raw_path):
                    os.link(raw_path, file_path)


            # shutil.copytree(
            #     os.path.join(idir, geneset, 'diff_reactome_enrich', 'svg'),
            #     os.path.join(odir,  '04DiffExpress_Geneset_Enrich',
            #                  '03Reactome', geneset,
            #                  "Reactome_enrich_pathways"))
            os.link(
                os.path.join(idir, geneset, 'diff_reactome_enrich', 'svg.tar.gz'),
                os.path.join(odir, '04DiffExpress_Geneset_Enrich',
                             '03Reactome', geneset,
                             "Reactome_enrich_pathways.tar.gz"))

        if species == "Homo_sapiens" and level == "G" :
            # 04 DO_Enrich
            os.makedirs(os.path.join(odir,  '04DiffExpress_Geneset_Enrich', '04DO'))
            for geneset in diff_genesets:
                os.makedirs(os.path.join(odir,  '04DiffExpress_Geneset_Enrich',
                                         '04DO', geneset))
                # os.link(
                #     os.path.join(idir, geneset, 'diff_do_enrich', 'do_enrichment.xls'),
                #     os.path.join(odir,  '04DiffExpress_Geneset_Enrich',
                #                  '04DO', geneset, "do_enrich_stat.xls"))

                # add by fwy 20210415 新增一列gene_name
                raw_path =  os.path.join(idir, geneset, 'diff_do_enrich', 'do_enrichment.xls')
                file_path = os.path.join(odir,  '04DiffExpress_Geneset_Enrich','04DO', geneset, "do_enrich_stat.xls")
                split = "|"
                add_columns = ["Genes"]
                try:
                    add_name.add_gene_name(raw_path, split, "do", add_columns, file_path)
                except:
                    os.link(raw_path, file_path)

    def set_somatic(self,result_dir,tmp_dir,odir):
        #统计部分文件统计
        if os.path.isdir(odir):
            shutil.rmtree(odir)
        os.makedirs(odir)
        # else:
        #     os.makedirs(odir)
        os.makedirs(os.path.join(odir,"SNV_vcf"))
        # CopyFile().linkdir(tmp_dir,odir)
        if os.path.exists(os.path.join(tmp_dir,"indel_anno.xls")):
            os.link(os.path.join(tmp_dir,"indel_anno.xls"),os.path.join(odir,"indel_anno.xls"))
        if os.path.exists(os.path.join(tmp_dir,"snp_anno.xls")):
            os.link(os.path.join(tmp_dir,"snp_anno.xls"),os.path.join(odir,"snp_anno.xls"))
        if os.path.exists(os.path.join(tmp_dir,"indel_position_distribution.xls")):
            indel_p = pd.read_table(os.path.join(tmp_dir,"indel_position_distribution.xls"))
            indel_pf = indel_p.set_index("Region")
            indel_pf.to_csv(os.path.join(odir,"indel_position_distribution.xls"),sep="\t")
        if os.path.exists(os.path.join(tmp_dir,"snp_position_distribution.xls")):
            snp_p = pd.read_table(os.path.join(tmp_dir,"snp_position_distribution.xls"))
            snp_pf = snp_p.set_index("Region")
            snp_pf.to_csv(os.path.join(odir,"snv_position_distribution.xls"),sep="\t")
        if os.path.exists(os.path.join(tmp_dir,"snp_depth_statistics.xls")):
            snp_d = pd.read_table(os.path.join(tmp_dir,"snp_depth_statistics.xls"))
            snp_df = snp_d.set_index("Depth")
            snp_df.to_csv(os.path.join(odir,"snv_depth_statistics.xls"),sep="\t")
        if os.path.exists(os.path.join(tmp_dir,"snp_transition_tranversion_statistics.xls")):
            snp_tt = pd.read_table(os.path.join(tmp_dir,"snp_transition_tranversion_statistics.xls"))
            snp_ttf = snp_tt.set_index("type")
            snp_ttf.to_csv(os.path.join(odir,"snv_transition_tranversion_statistics.xls"),sep="\t")
        if os.path.exists(os.path.join(result_dir,"somatic_predict", "predict_detail")):
            os.link(os.path.join(result_dir,"somatic_predict", "predict_detail"), os.path.join(odir, "Columns_of_dbNSFP_variant.xls"))
        if os.path.exists(os.path.join(result_dir,"somatic_predict", "predict_stat")):
            os.link(os.path.join(result_dir,"somatic_predict", "predict_stat"), os.path.join(odir, "dbNSFP_variant_stat.xls"))
        if os.path.exists(os.path.join(result_dir, "somatic_predict", "predict.png")):
            os.link(os.path.join(result_dir, "somatic_predict", "predict.png"),
                    os.path.join(odir, "predict.png"))
        if os.path.exists(os.path.join(result_dir, "somatic_predict", "predict.pdf")):
            os.link(os.path.join(result_dir, "somatic_predict", "predict.pdf"),
                    os.path.join(odir, "predict.pdf"))
        if os.path.exists(os.path.join(result_dir, "somatic_predict", "predict.svg")):
            os.link(os.path.join(result_dir, "somatic_predict", "predict.svg"),
                    os.path.join(odir, "predict.svg"))
        if os.path.exists(os.path.join(tmp_dir, "tnhaplotyper")):
            vcf_list = os.listdir(os.path.join(tmp_dir, "tnhaplotyper"))
            for vcf in vcf_list:
                if not vcf.endswith("idx"):
                    os.link(os.path.join(tmp_dir, "tnhaplotyper",vcf),os.path.join(odir,"SNV_vcf",vcf))

    def set_gene_fusion(self,idir,odir):
        if os.path.isdir(odir):
            shutil.rmtree(odir)
        os.makedirs(odir)
        # else:
        #     os.makedirs(odir)
        os.link(os.path.join(idir, 'fusion_stat.txt'), os.path.join(odir, 'fusion_stat.txt'))
        if os.path.isdir(os.path.join(odir,"Star_fusion")):
            shutil.rmtree(os.path.join(odir,"Star_fusion"))
        else:
            os.makedirs(os.path.join(odir,"Star_fusion"))
        sample_list = os.listdir(os.path.join(idir,"star_fusion"))
        for sample in sample_list:
            os.makedirs(os.path.join(odir,"Star_fusion",sample))
            os.link(os.path.join(idir,"star_fusion",sample,"star-fusion.fusion_predictions.tsv"),
                    os.path.join(odir,"Star_fusion",sample,"star-fusion.fusion_predictions.tsv"))
            os.link(os.path.join(idir, "star_fusion", sample, "star-fusion.fusion_predictions.abridged.tsv"),
                    os.path.join(odir, "Star_fusion", sample, "star-fusion.fusion_predictions.abridged.tsv"))
            try:
                os.makedirs(os.path.join(odir, "Star_fusion", sample,"fusion_inspector"))
                os.link(os.path.join(idir, "star_fusion", sample, "fusion_inspector","output","finspector.fa"),
                        os.path.join(odir, "Star_fusion", sample, "fusion_inspector","finspector.fa"))
                os.link(os.path.join(idir, "star_fusion", sample, "fusion_inspector","output", "finspector.bed"),
                        os.path.join(odir, "Star_fusion", sample, "fusion_inspector","finspector.bed"))
                os.link(os.path.join(idir, "star_fusion", sample, "fusion_inspector","output","finspector.junction_reads.bam"),
                        os.path.join(odir, "Star_fusion", sample, "fusion_inspector","finspector.junction_reads.bam"))
                os.link(os.path.join(idir, "star_fusion", sample, "fusion_inspector","output","finspector.spanning_reads.bam"),
                        os.path.join(odir, "Star_fusion", sample, "fusion_inspector","finspector.spanning_reads.bam"))
            except:
                pass


    def set_assembly(self,idir, odir):
        if os.path.isdir(odir):
            shutil.rmtree(odir)
        os.makedirs(odir)
        # else:
        #     os.mkdir(odir)
        # transcript length distribution table

        #医学版取消 by fwy 20201215
        # for file in glob.glob(os.path.join(idir, 'Statistics/trans_count_stat_*.txt')):
        #     step = re.search(r'trans_count_stat_(\d+?).txt', file).group(1)
        #     df0 = pd.read_table(file)
        #     df0.index.name = 'Length'
        #     df0 = df0.rename({'all_transcripts': 'Number'}, axis=1)
        #     df0.to_csv(os.path.join(odir, 'length_distribution.{}.xls'.format(step)), sep='\t')
        # new transcript type statistics table
        df1 = pd.DataFrame([
            ('=', 'Complete match of intron chain'),
            ('c', 'Contained'),
            ('e',
             'Single exon transfrag overlapping a reference exon and at least 10 bp of a reference intron,indiction a possible pre-mRNA fragment'),
            ('i', 'Transfrag falling entirely within a reference intron'),
            (
                'j',
                'Potentially novel isoform (fragment):at least one splice junction is shared with a reference transcript'),
            ('o', 'Generic exonic overlap with a referfence transcript'),
            ('p', 'Possible polymerase run-on fragment(within 2Kbases of a reference transcript'),
            ('s',
             'An intron of the transfrag overlaps a reference intron on the opposite strand(likely due to read mapping errors'),
            ('u', 'Unknown,intergenic transcript'),
            ('x', 'Exonic overlap with reference on the opposite strand'),
            ('r',
             'Repeat. Currently determined by looking at the soft-masked reference sequence and applied to transcripts where at least 50% of the bases are lower case'),
            ('.', 'Tracking file only,indicates multiple classifications')], columns=['Class code', 'Description'])
        df2 = pd.read_table(os.path.join(idir, 'Statistics/code_num.txt'), names=['Class code', 'Number'], usecols=[0, 2])
        df3 = pd.merge(df1, df2, how='outer')
        df3 = df3.fillna(0)
        df3['Number'] = df3['Number'].astype(np.int64)
        df3.to_csv(os.path.join(odir, 'classcode_statistics.xls'), sep='\t', index=False)
        # sequence files
        sdir = os.path.join(odir, 'Sequence')
        if os.path.isdir(sdir):
            shutil.rmtree(sdir)
        os.mkdir(sdir)
        os.link(os.path.join(idir, 'NewTranscripts/all_transcripts.fa'), os.path.join(sdir, 'all_transcripts.fa'))
        os.link(os.path.join(idir, 'NewTranscripts/new_transcripts.fa'), os.path.join(sdir, 'new_transcripts.fa'))
        os.link(os.path.join(idir, 'NewTranscripts/ref_and_new.gtf'), os.path.join(sdir, 'all_transcripts.gtf'))
        os.link(os.path.join(idir, 'NewTranscripts/new_transcripts.gtf'), os.path.join(sdir, 'new_transcripts.gtf'))
        os.link(os.path.join(idir, 'NewTranscripts/trans2gene'), os.path.join(sdir, 'trans2gene.txt'))


    def set_rmats(self,idir, odir, task_id, cmp_list):
        if os.path.isdir(odir):
            shutil.rmtree(odir)
        os.mkdir(odir)
        # alternative splicing event statistics table
        # os.link(os.path.join(idir, 'sample.event.count.JC.txt'), os.path.join(odir, 'event_type.JC.xls'))
        # os.link(os.path.join(idir, 'sample.event.count.JCEC.txt'), os.path.join(odir, 'event_type.JCEC.xls'))
        # alternative splicing event detail table
        # gid2des_dct,gid2name_dct = self.get_gid2des_dct(task_id)
        # for document in self.database['sg_splicing_rmats'].find({'task_id': task_id}):
        #     splicing_id = document['main_id']
        #     s1, s2 = document['compare_plan'].split('|')
        for s1, s2 in cmp_list:
            output_dir = os.path.join(odir, '{}_vs_{}_AS'.format(s1, s2))
            self.export_rmats_detail(gid2des_dct,gid2name_dct, s2 = s2, s1=s1, output_dir=output_dir)
        # alternative splicing intra-group event statistics table and pattern statistics table
        # for document in self.database['sg_splicing_rmats_stats'].find({'task_id': task_id}):
        #     stat_id = document['main_id']
        #     print 'stat_id -> ({})'.format(stat_id)
        #     print 'group -> ({})'.format(document['group'])
        #     for k, v in document['group'].items():
        #         if v == 's1':
        #             s1 = k
        #         if v == 's2':
        #             s2 = k
        #     else:
        #         print 's1 -> ({})'.format(s1)
        #         print 's2 -> ({})'.format(s2)
            output_dir = os.path.join(odir, '{}_vs_{}_AS'.format(s1, s2))
            self.export_rmats_diff_stats(output_dir)
            self.export_rmats_psi(output_dir)


    def get_gid2des_dct(self,task_id):
        query_id = self.database['sg_exp'].find_one({'task_id': task_id})['main_id']
        lst = list({(document['gene_id'], document['description']) for document in
                    self.database['sg_exp_detail'].find({'exp_id': query_id})})
        lst2 = list({(document['gene_id'], document['gene_name']) for document in
                    self.database['sg_exp_detail'].find({'exp_id': query_id})})
        return dict(lst),dict(lst2)



    def export_rmats_detail(self,gid2des_dct,gid2name_dct, s2, s1, output_dir):
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        dct = {k: {'JC': dict(), 'JCEC': dict()} for k in ('SE', 'MXE', 'A3SS', 'A5SS', 'RI')}

        def get_columns(event_type, junction_type):
            lst = ['AS ID', 'Gene ID', 'Gene name', 'Gene description', 'Novel AS', 'Chr', 'Strand',
                   'Diff significant', 'IncLevelDiff ({}/{},ΔPSI)'.format(s2, s1)]
            lst.extend({
                           'JC': ['P Value JunctionCountOnly', 'FDR JunctionCountOnly'],
                           'JCEC': ['P Value ReadsOnTargetAndJunctionCounts', 'FDR ReadsOnTargetAndJunctionCounts']
                       }[junction_type])
            lst.extend({
                           'SE': ['InclusionTranscripts', 'SkippingTranscripts', 'ExonStart', 'ExonEnd', 'UpstreamES',
                                  'UpstreamEE', 'DownstreamES', 'DownstreamEE'],
                           'MXE': ['1stExonTranscripts', '2ndExonTranscripts', '1stExonStart', '1stExonEnd', '2ndExonStart',
                                   '2ndExonEnd'],
                           'A3SS': ['LongExonTranscripts', 'ShortExonTranscripts', 'LongExonStart', 'LongExonEnd',
                                    'ShortES', 'ShortEE', 'FlankingES', 'FlankingEE'],
                           'A5SS': ['LongExonTranscripts', 'ShortExonTranscripts', 'LongExonStart', 'LongExonEnd',
                                    'ShortES', 'ShortEE', 'FlankingES', 'FlankingEE'],
                           'RI': ['RetainTranscripts', 'AbandonTranscripts', 'RiExonStart', 'RiExonEnd']
                       }[event_type])
            lst.extend({
                           'JC': ['IJC {}'.format(s2), 'SJC {}'.format(s2), 'IJC {}'.format(s1), 'SJC {}'.format(s1)],
                           'JCEC': ['IC {}'.format(s2), 'SC {}'.format(s2), 'IC {}'.format(s1), 'SC {}'.format(s1)]
                       }[junction_type])
            lst.extend([
                'IncFormLen', 'SkipFormLen', 'IncLevel1 (PSI {})'.format(s2), 'IncLevel2 (PSI {})'.format(s1),
                'Average IncLevel1 (Average PSI {})'.format(s2), 'Average IncLevel2 (Average PSI {})'.format(s1),
                'Increase Inclusion {}'.format(s2), 'Increase Exclusion {}'.format(s2),
                'Increase Inclusion {}'.format(s1), 'Increase Exclusion {}'.format(s1)
            ])
            return lst

        if not os.path.isdir(os.path.join(output_dir, 'JC')):
            os.mkdir(os.path.join(output_dir, 'JC'))
        if not os.path.isdir(os.path.join(output_dir, 'JCEC')):
            os.mkdir(os.path.join(output_dir, 'JCEC'))
        for et, jt2dcts in dct.items():
            print 'event_type -> ({})'.format(et)
            for jt, dcts in jt2dcts.items():
                print 'junction_type -> ({})'.format(jt)
                df = pd.DataFrame(dcts.values())
                df = df.reindex(get_columns(et, jt), axis=1)
                df.to_csv(os.path.join(output_dir, '{}/{}.detail.xls'.format(jt, et)), sep='\t', index=False)
        else:
            print 'succeed in exporting files to {}'.format(output_dir)


    def export_rmats_diff_stats(self,output_dir):
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        # document = self.database['sg_splicing_rmats_diff_stats'].find_one({'stat_id': ObjectId(stat_id)})
        # index = ['JunctionCountOnly(JC)', 'ReadsOnTargetAndJunctionCounts(JCEC)', 'JC&JCEC', 'JC|JCEC']
        # df = pd.DataFrame(document['diff_stats'], index=index).T
        # df = df.rename({i: i.upper() for i in df.index})
        # df = df.reindex(['SE', 'MXE', 'A3SS', 'A5SS', 'RI', 'TOTAL'])
        # df.index.name = 'AS type'
        # df.to_csv(os.path.join(output_dir, 'diff_event_stats.xls'), sep='\t')
        # print 'succeed in exporting {}'.format(os.path.join(output_dir, 'diff_event_stats.txt'))


    def export_rmats_psi(self,stat_id, output_dir):
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        document = self.database['sg_splicing_rmats_psi'].find_one({'stat_id': ObjectId(stat_id)})
        index = ['Exclusion (Increase exclusion in case，ΔPSI<0)',
                 'Inclusion (Increase inclusion in case, ΔPSI>0)',
                 'Total events']
        df = pd.DataFrame(document['s1_jc'], index=index).T
        df = df.rename({i: i.upper() for i in df.index})
        df = df.reindex(['SE', 'MXE', 'A3SS', 'A5SS', 'RI', 'TOTAL'])
        df.index.name = 'AS type'
        jc_ofile = os.path.join(output_dir, 'diff_pattern_stats.JC.xls')
        df.to_csv(jc_ofile, sep='\t')
        print 'succeed in exporting {}'.format(jc_ofile)
        df = pd.DataFrame(document['s1_all'], index=index).T
        df = df.rename({i: i.upper() for i in df.index})
        df = df.reindex(['SE', 'MXE', 'A3SS', 'A5SS', 'RI', 'TOTAL'])
        df.index.name = 'AS type'
        jcec_ofile = os.path.join(output_dir, 'diff_pattern_stats.JCEC.xls')
        df.to_csv(jcec_ofile, sep='\t')
        print 'succeed in exporting {}'.format(jcec_ofile)
