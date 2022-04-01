# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import re
import os
import sys

class LncFileDes(object):
    def __init__(self):
        self.Background = [
            [".", "", "项目背景结果目录", 0],
            ['genome_info.xls', 'xls', '基因组注释信息表', 0],
            ['sample_info.xls', 'xls', '样本信息表', 0],
            ['software_info.xls', 'xls', '软件信息表', 0],
        ]

        self.Basic_Analysis = [
            [".", "", "基础分析结果目录", 0],
            ["01 QC", "", "测序数据质控", 0],
            ['02 Align', '', '序列比对分析',0],
            ['03 Assemble', '', '转录本组装', 0],
            ["01 QC/cleandata", "", "质控结果目录", 1, "211029"],
            ["01 QC/rawdata_statistics.xls", "", "原始数据统计表", 0, "211027"],
            ["01 QC/cleandata_statistics.xls", "", "质控数据统计表", 0, "211028"],
            ["02 Align/AlignBam", "", "比对结果bam目录", 1, "211031"],
            ["02 Align/AlignStat", "", "比对结果统计目录", 0, "211032"],
            ["02 Align/QualityAssessment", "", "比对结果整体评估目录", 0, "211033"],
            ["02 Align/QualityAssessment/", '', '比对结果整体评估', 0],
            [r"02 Align/AlignBam/.*\.bam", "", "样本比对bam文件", 1, "211002"],
            [r"02 Align/AlignStat/.*_align_stat\.txt", "txt", "样本比对统计结果表", 0, "211003"],
            [r"02 Align/QualityAssessment/.*\.chr_distribution\.xls", "xls", "Reads在不同染色体的分布统计表", 0, "211004"],
            [r"02 Align/QualityAssessment/.*\.geneBodyCoverage\.txt", "txt", "基因覆盖度分布结果文件", 0, "211005"],
            [r"02 Align/QualityAssessment/.*\.region_distribution\.xls", "xls", "Reads在不同区域的分布统计表", 0, "211006"],
            [r"02 Align/QualityAssessment/.*\.saturation\.xls", "xls", "测序饱和度分析结果文件", 0, "211007"],
            ["03 Assemble/ref_and_new.gtf", "", "组装结果GTF文件", 0, "211038"],
            ["03 Assemble/new_transcripts.fa", "", "新预测转录本序列", 0, "211039"],
            ["03 Assemble/new_transcripts.gtf", "", "新预测转录本GTF文件", 0, "211040"],
            ["03 Assemble/all_transcripts.fa", "", "所有转录本序列文件", 0, "211041"],
            ['03 Assemble/all_cds.fa', 'fasta', 'CDS序列文件', 0],
            ['03 Assemble/all_pep.fa', 'fasta', '蛋白序列文件', 0],
            ['03 Assemble/all_id.xls', 'xls', '基因转录本蛋白ID对应关系文件', 0],
            ["03 Assemble/trans2gene.txt", "", "转录本与基因的对应关系文件", 0, "211042"]
        ]

        self.Annotation = [
            [".", "", "功能注释结果目录", 0],
            ["allannot_class", "", "已知基因与新预测合并注释结果目录", 0],
            ["allannot_class/all_annot.xls", "", "注释结果详情表", 0],
            ["allannot_class/all_tran2gene.txt", "", "转录本与基因对应关系列表", 0],
            ["allannot_class/all_stat.xls", "", "注释结果统计表", 0],
            ["allannot_class/cog", "", "cog注释结果目录", 0],
            ["allannot_class/cog/", "", "cog注释结果文件", 0],
            ["allannot_class/cog/cog_list_tran.xls", "", "cog注释详情表", 0],
            ["allannot_class/cog/cog_venn_gene.txt", "", "存在cog注释的基因列表", 0],
            ["allannot_class/cog/cog_venn_tran.txt", "", "存在cog注释的转录本列表", 0],
            ["allannot_class/go", "", "go注释结果目录", 0],
            ["allannot_class/go/", "", "go注释结果文件", 0],
            ["allannot_class/go/go_lev2_gene.stat.xls", "", "基因go二级分类统计表", 0],
            ["allannot_class/go/go_lev2_tran.stat.xls", "", "转录本go二级分类统计表", 0],
            ["allannot_class/go/go_lev3_gene.stat.xls", "", "基因go三级分类统计表", 0],
            ["allannot_class/go/go_lev3_tran.stat.xls", "", "转录本go三级分类统计表", 0],
            ["allannot_class/go/go_lev4_gene.stat.xls", "", "基因go四级分类统计表", 0],
            ["allannot_class/go/go_lev4_tran.stat.xls", "", "转录本go四级分类统计表", 0],
            ["allannot_class/go/go_list_gene.xls", "", "基因go注释列表", 0],
            ["allannot_class/go/go_list_tran.xls", "", "转录本go注释列表", 0],
            ["allannot_class/go/go_venn_gene.txt", "", "存在go注释的基因列表", 0],
            ["allannot_class/go/go_venn_tran.txt", "", "存在go注释的转录本列表", 0],
            ["allannot_class/kegg", "", "kegg注释结果目录", 0],
            ["allannot_class/kegg/", "", "kegg注释结果文件", 0],
            ["allannot_class/kegg/kegg_gene_gene.xls", "", "基因kegg注释详情表", 0],
            ["allannot_class/kegg/kegg_gene_tran.xls", "", "转录本kegg注释详情表", 0],
            ["allannot_class/kegg/kegg_pathway_gene_dir", "", "基因kegg注释通路图文件", 0],
            ["allannot_class/kegg/kegg_pathway_gene.xls", "", "基因kegg pathway注释详情表", 0],
            ["allannot_class/kegg/kegg_pathway_tran_dir", "", "转录本kegg通路图文件", 0],
            ["allannot_class/kegg/kegg_pathway_tran.xls", "", "转录本kegg pathway注释详情表", 0],
            ["allannot_class/kegg/kegg_venn_gene.txt", "", "存在kegg注释的基因列表", 0],
            ["allannot_class/kegg/kegg_venn_tran.txt", "", "存在kegg注释的转录本列表", 0],
            [r'allannot_class/kegg/kegg_pathway_(gene|tran)_dir/.*\.png', "png", "kegg注释通路图png", 0],
            [r'allannot_class/kegg/kegg_pathway_(gene|tran)_dir/.*\.html', "html", "kegg注释通路图html", 0],
            ["allannot_class/nr", "", "nr注释结果目录", 0],
            ["allannot_class/nr/", "", "nr注释结果文件", 0],
            ["allannot_class/nr/nr_blast_gene.xls", "", "基因nr注释详情表", 0],
            ["allannot_class/nr/nr_blast_tran.xls", "", "转录本nr注释详情表", 0],
            ["allannot_class/nr/nr_venn_gene.txt", "", "存在nr注释的基因列表", 0],
            ["allannot_class/nr/nr_venn_tran.txt", "", "存在nr注释的转录本列表", 0],
            ["allannot_class/pfam", "", "pfam注释结果目录", 0],
            ["allannot_class/pfam/", "", "pfam注释结果文件", 0],
            ["allannot_class/pfam/pfam_domain_gene.xls", "", "基因pfam注释详情表", 0],
            ["allannot_class/pfam/pfam_domain_tran.xls", "", "转录本pfam注释详情表", 0],
            ["allannot_class/pfam/pfam_venn_gene.txt", "", "存在pfam注释的基因列表", 0],
            ["allannot_class/pfam/pfam_venn_tran.txt", "", "存在pfam注释的转录本列表", 0],
            ["allannot_class/swissprot", "", "swissprot注释结果目录", 0],
            ["allannot_class/swissprot/", "", "swissprot注释结果文件", 0],
            ["allannot_class/swissprot/swissprot_blast_gene.xls", "", "基因swissprot注释详情表", 0],
            ["allannot_class/swissprot/swissprot_blast_tran.xls", "", "转录本swissprot注释详情表", 0],
            ["allannot_class/swissprot/swissprot_venn_gene.txt", "", "存在swissprot注释的基因列表", 0],
            ["allannot_class/swissprot/swissprot_venn_tran.txt", "", "存在swissprot注释的基因列表", 0],
            ["newannot_class", "", "新预测转录本注释结果目录", 0],
            ["newannot_class/all_annot.xls", "", "注释结果详情表", 0],
            ["newannot_class/all_tran2gene.txt", "", "转录本与基因对应关系列表", 0],
            ["newannot_class/all_stat.xls", "", "注释结果统计表", 0],
            ["newannot_class/cog", "", "cog注释结果目录", 0],
            ["newannot_class/cog/", "", "cog注释结果文件", 0],
            ["newannot_class/cog/cog_list_tran.xls", "", "cog注释详情表", 0],
            ["newannot_class/cog/cog_venn_gene.txt", "", "存在cog注释的基因列表", 0],
            ["newannot_class/cog/cog_venn_tran.txt", "", "存在cog注释的转录本列表", 0],
            ["newannot_class/go", "", "go注释结果目录", 0],
            ["newannot_class/go/", "", "go注释结果文件", 0],
            ["newannot_class/go/go_lev2_gene.stat.xls", "", "基因go二级分类统计表", 0],
            ["newannot_class/go/go_lev2_tran.stat.xls", "", "转录本go二级分类统计表", 0],
            ["newannot_class/go/go_lev3_gene.stat.xls", "", "基因go三级分类统计表", 0],
            ["newannot_class/go/go_lev3_tran.stat.xls", "", "转录本go三级分类统计表", 0],
            ["newannot_class/go/go_lev4_gene.stat.xls", "", "基因go四级分类统计表", 0],
            ["newannot_class/go/go_lev4_tran.stat.xls", "", "转录本go四级分类统计表", 0],
            ["newannot_class/go/go_list_gene.xls", "", "基因go注释列表", 0],
            ["newannot_class/go/go_list_tran.xls", "", "转录本go注释列表", 0],
            ["newannot_class/go/go_venn_gene.txt", "", "存在go注释的基因列表", 0],
            ["newannot_class/go/go_venn_tran.txt", "", "存在go注释的转录本列表", 0],
            ["newannot_class/kegg", "", "kegg注释结果目录", 0],
            ["newannot_class/kegg/", "", "kegg注释结果文件", 0],
            ["newannot_class/kegg/kegg_gene_gene.xls", "", "基因kegg注释详情表", 0],
            ["newannot_class/kegg/kegg_gene_tran.xls", "", "转录本kegg注释详情表", 0],
            ["newannot_class/kegg/kegg_pathway_gene_dir", "", "基因kegg注释通路图文件", 0],
            ["newannot_class/kegg/kegg_pathway_gene.xls", "", "基因kegg pathway注释详情表", 0],
            ["newannot_class/kegg/kegg_pathway_tran_dir", "", "转录本kegg注释通路图文件", 0],
            ["newannot_class/kegg/kegg_pathway_tran.xls", "", "转录本kegg pathway注释详情表", 0],
            ["newannot_class/kegg/kegg_venn_gene.txt", "", "存在kegg注释的基因列表", 0],
            ["newannot_class/kegg/kegg_venn_tran.txt", "", "存在kegg注释的转录本列表", 0],
            [r'newannot_class/kegg/kegg_pathway_(gene|tran)_dir/.*\.png', "png", "kegg注释通路图png", 0],
            [r'newannot_class/kegg/kegg_pathway_(gene|tran)_dir/.*\.html', "html", "kegg注释通路图html", 0],
            ["newannot_class/nr", "", "nr注释结果目录", 0],
            ["newannot_class/nr/", "", "nr注释结果文件", 0],
            ["newannot_class/nr/nr_blast_gene.xls", "", "基因nr注释详情表", 0],
            ["newannot_class/nr/nr_blast_tran.xls", "", "转录本nr注释详情表", 0],
            ["newannot_class/nr/nr_venn_gene.txt", "", "存在nr注释的基因列表", 0],
            ["newannot_class/nr/nr_venn_tran.txt", "", "存在nr注释的转录本列表", 0],
            ["newannot_class/pfam", "", "pfam注释结果目录", 0],
            ["newannot_class/pfam/", "", "pfam注释结果文件", 0],
            ["newannot_class/pfam/pfam_domain_gene.xls", "", "基因pfam注释详情表", 0],
            ["newannot_class/pfam/pfam_domain_tran.xls", "", "转录本pfam注释详情表", 0],
            ["newannot_class/pfam/pfam_venn_gene.txt", "", "存在pfam注释的基因列表", 0],
            ["newannot_class/pfam/pfam_venn_tran.txt", "", "存在pfam注释的转录本列表", 0],
            ["newannot_class/swissprot", "", "swissprot注释结果目录", 0],
            ["newannot_class/swissprot/", "", "swissprot注释结果文件", 0],
            ["newannot_class/swissprot/swissprot_blast_gene.xls", "", "基因swissprot注释详情表", 0],
            ["newannot_class/swissprot/swissprot_blast_tran.xls", "", "转录本swissprot注释详情表", 0],
            ["newannot_class/swissprot/swissprot_venn_gene.txt", "", "存在swissprot注释的基因列表", 0],
            ["newannot_class/swissprot/swissprot_venn_tran.txt", "", "存在swissprot注释的基因列表", 0],
            ["refannot_class", "", "已知基因/转录本注释结果目录", 0],
            ["refannot_class/all_annot.xls", "", "注释结果详情表", 0],
            ["refannot_class/all_tran2gene.txt", "", "转录本与基因对应关系列表", 0],
            ["refannot_class/all_stat.xls", "", "注释结果统计表", 0],
            ["refannot_class/cog", "", "cog注释结果目录", 0],
            ["refannot_class/cog/", "", "cog注释结果文件", 0],
            ["refannot_class/cog/cog_list_tran.xls", "", "cog注释详情表", 0],
            ["refannot_class/cog/cog_venn_gene.txt", "", "存在cog注释的基因列表", 0],
            ["refannot_class/cog/cog_venn_tran.txt", "", "存在cog注释的转录本列表", 0],
            ["refannot_class/go", "", "go注释结果目录", 0],
            ["refannot_class/go/", "", "go注释结果文件", 0],
            ["refannot_class/go/go_lev2_gene.stat.xls", "", "基因go二级分类统计表", 0],
            ["refannot_class/go/go_lev2_tran.stat.xls", "", "转录本go二级分类统计表", 0],
            ["refannot_class/go/go_lev3_gene.stat.xls", "", "基因go三级分类统计表", 0],
            ["refannot_class/go/go_lev3_tran.stat.xls", "", "转录本go三级分类统计表", 0],
            ["refannot_class/go/go_lev4_gene.stat.xls", "", "基因go四级分类统计表", 0],
            ["refannot_class/go/go_lev4_tran.stat.xls", "", "转录本go四级分类统计表", 0],
            ["refannot_class/go/go_list_gene.xls", "", "基因go注释列表", 0],
            ["refannot_class/go/go_list_tran.xls", "", "转录本go注释列表", 0],
            ["refannot_class/go/go_venn_gene.txt", "", "存在go注释的基因列表", 0],
            ["refannot_class/go/go_venn_tran.txt", "", "存在go注释的转录本列表", 0],
            ["refannot_class/kegg", "", "kegg注释结果目录", 0],
            ["refannot_class/kegg/", "", "kegg注释结果文件", 0],
            ["refannot_class/kegg/kegg_gene_gene.xls", "", "基因kegg注释详情表", 0],
            ["refannot_class/kegg/kegg_gene_tran.xls", "", "转录本kegg注释详情表", 0],
            ["refannot_class/kegg/kegg_pathway_gene_dir", "", "基因kegg注释通路图文件", 0],
            ["refannot_class/kegg/kegg_pathway_gene.xls", "", "基因kegg pathway注释详情表", 0],
            ["refannot_class/kegg/kegg_pathway_tran_dir", "", "转录本kegg通路图文件", 0],
            ["refannot_class/kegg/kegg_pathway_tran.xls", "", "转录本kegg pathway注释详情表", 0],
            ["refannot_class/kegg/kegg_venn_gene.txt", "", "存在kegg注释的基因列表", 0],
            ["refannot_class/kegg/kegg_venn_tran.txt", "", "存在kegg注释的转录本列表", 0],
            [r'refannot_class/kegg/kegg_pathway_(gene|tran)_dir/.*\.png', "png", "kegg注释通路图png", 0],
            [r'refannot_class/kegg/kegg_pathway_(gene|tran)_dir/.*\.html', "html", "kegg注释通路图html", 0],
            ["refannot_class/nr", "", "nr注释结果目录", 0],
            ["refannot_class/nr/", "", "nr注释结果文件", 0],
            ["refannot_class/nr/nr_blast_gene.xls", "", "基因nr注释详情表", 0],
            ["refannot_class/nr/nr_blast_tran.xls", "", "转录本nr注释详情表", 0],
            ["refannot_class/nr/nr_venn_gene.txt", "", "存在nr注释的基因列表", 0],
            ["refannot_class/nr/nr_venn_tran.txt", "", "存在nr注释的转录本列表", 0],
            ["refannot_class/pfam", "", "pfam注释结果目录", 0],
            ["refannot_class/pfam/", "", "pfam注释结果文件", 0],
            ["refannot_class/pfam/pfam_domain_gene.xls", "", "基因pfam注释详情表", 0],
            ["refannot_class/pfam/pfam_domain_tran.xls", "", "转录本pfam注释详情表", 0],
            ["refannot_class/pfam/pfam_venn_gene.txt", "", "存在pfam注释的基因列表", 0],
            ["refannot_class/pfam/pfam_venn_tran.txt", "", "存在pfam注释的转录本列表", 0],
            ["refannot_class/swissprot", "", "swissprot注释结果目录", 0],
            ["refannot_class/swissprot/", "", "swissprot注释结果文件", 0],
            ["refannot_class/swissprot/swissprot_blast_gene.xls", "", "基因swissprot注释详情表", 0],
            ["refannot_class/swissprot/swissprot_blast_tran.xls", "", "转录本swissprot注释详情表", 0],
            ["refannot_class/swissprot/swissprot_venn_gene.txt", "", "存在swissprot注释的基因列表", 0],
            ["refannot_class/swissprot/swissprot_venn_tran.txt", "", "存在swissprot注释的基因列表", 0],
            ["newannot_mapdb", "", "新转录本与数据库比对结果目录", 0],
            ["newannot_mapdb/eggnog", "", "新转录本eggnog数据库比对结果目录", 0],
            ["newannot_mapdb/go", "", "新转录本GO数据库比对结果目录", 0],
            ["newannot_mapdb/kegg", "", "新转录本kegg数据库比对结果目录", 0],
            ["newannot_mapdb/nr", "", "新转录本nr数据库比对结果目录", 0],
            ["newannot_mapdb/swissprot", "", "新转录本swissprot数据库比对结果目录", 0],

        ]

        self.LncRNA_Analysis = [
            ['.', '', 'lncRNA分析结果目录', 0],
            ['01 Known_LncRNA', '', '已知lncRNA鉴定', 0],
            ['02 Novel_LncRNA', '', '新lncRNA预测', 0],
            ['03 LncRNA_stat', '', 'lncRNA统计', 0],
            ['01 Known_LncRNA/', '', '已知lncRNA注释信息', 0],
            ['01 Known_LncRNA/known_lncrna.gtf', 'gtf', '已知lncRNA的GTF文件', 0],
            ['01 Known_LncRNA/known_lncrna.fa', 'fasta', '已知lncRNA序列', 0],
            ['01 Known_LncRNA/known_lncrna_detail.xls', 'xls', '已知lncRNA详情表', 0],
            ['02 Novel_LncRNA/', '', '新lncRNA注释信息', 0],
            ['02 Novel_LncRNA/novel_lncrna.gtf', 'gtf', '新lncRNA的GTF文件', 0],
            ['02 Novel_LncRNA/novel_lncrna.fa', 'fasta', '新lncRNA序列', 0],
            ['02 Novel_LncRNA/novel_lncrna_predict_stat.xls', 'xls', '新lncRNA预测统计表', 0],
            ['02 Novel_LncRNA/novel_lncrna_predict_detail.xls', 'xls', '新lncRNA预测详情表', 0],
            ['03 LncRNA_stat/', '', 'lncRNA统计', 0],
            ['03 LncRNA_stat/lncrna_stat_in_category.xls', 'xls', 'lncRNA分类统计表', 0],
            ['03 LncRNA_stat/lncrna_stat_in_sample.xls', 'xls', 'lncRNA统计表', 0],
        ]

        self.Express = [
            [".", "", "表达量分析结果目录", 0,  "211386"],
            ["01 ExpAnnalysis", "", "表达定量结果目录", 0,  "211387"],
            ["02 ExpCorr", "", "样本间相关性分析", 0, "211398"],
            ["03 ExpPCA", "", "样本间PCA分析", 0, "211400"],
            ["01 ExpAnnalysis/gene.count.matrix.xls", "", "gene count表达定量结果表", 0,  "211388"],
            ["01 ExpAnnalysis/gene.tpm.matrix.xls", "", "gene tpm表达定量结果表", 0,  "211389"],
            ["01 ExpAnnalysis/gene.fpkm.matrix.xls", "", "gene fpkm表达定量结果表", 0,  "211390"],
            ["01 ExpAnnalysis/transcript.count.matrix.xls", "", "transcript count表达定量结果表", 0,  "211391"],
            ["01 ExpAnnalysis/transcript.tpm.matrix.xls", "", "transcript tpm表达定量结果表", 0,  "211392"],
            ["01 ExpAnnalysis/transcript.fpkm.matrix.xls", "", "transcript fpkm表达定量结果表", 0,  "211393"],
            ["01 ExpAnnalysis/gene.tpm.matrix.annot.xls", "", "gene tpm表达定量注释结果表", 0,  "211394"],
            ["01 ExpAnnalysis/gene.fpkm.matrix.annot.xls", "", "gene fpkm表达定量注释结果表", 0,  "211395"],
            ["01 ExpAnnalysis/transcript.tpm.matrix.annot.xls", "", "transcript tpm表达定量注释结果表", 0,  "211396"],
            ["01 ExpAnnalysis/transcript.fpkm.matrix.annot.xls", "", "transcript fpkm表达定量注释结果表", 0,  "211397"],
            [r'02 ExpCorr/(G|T)_sample_correlation.xls', "xls", "样本间相关系数表", 0],
            [r'03 ExpPCA/(G|T)_PCA.xls', "", "样本间PCA分析结果表", 0,  "211401"],
            [r'03 ExpPCA/(G|T)_Explained_variance_ratio.xls', "", "样本间PCA主成分解释表", 0,  "211402"],
        ]

        self.Diff_Express = [
            [".", "", "表达量差异分析结果目录", 0],
            [r".*_vs_.*\.sizeFactor.xls", "xls", "差异分析均一化因子表", 0, "211008"],
            [r"G_.*\.DE\.list", "xls", "差异基因列表", 0, "211009"],
            [r"T_.*\.DE\.list", "xls", "差异转录本列表", 0, "211009"],
            [r'(G|T)_diff_summary_.*\.xls', 'xls', '表达量差异统计表', 0],
            [r'.*_vs_.*\..*\.xls', 'xls', '表达量差异详情表', 0],
            [r'.*_vs_.*\..*\.annot.xls', 'xls', '表达量差异注释表', 0],
            [r'(G|T)_total_diff_stat.*\.xls', 'xls', '差异表达基因详情总表', 0],
        ]
        # [r".*_vs_.*\.xls", "xls", "差异表达基因详情总表", 0, "211008"],
        # [r".*summary\.xls", "xls", "差异统计结果表", 0, "211010"],
        # [r".*_vs_.*\.*annot\.xls", "xls", "差异统计注释结果表", 0, "211011"],
        # [r".*_vs_.*\.*annot\.xls", "xls", "差异统计注释结果表", 0, "211012"],

        self.LncRNA_Target = [
            [".", "", "lncRNA靶基因预测结果目录", 0],
            ['cis_annot.xls', 'xls', 'cis作用靶基因预测详情表', 0],
            ['trans_annot.xls', 'xls', 'trans作用靶基因预测详情表', 0],
        ]

        self.LncRNA_Structure = [
            ['.', '', 'lncRNA结构分析结果目录', 0],
            ['01 LncRNA_Family', '', 'lncRNA家族分析', 0],
            ['02 miRNA_Pre', '', 'miRNA前体预测', 0],
            ['01 LncRNA_Family/lncRNA_family.xls', 'xls', 'lncRNA家族信息表', 0],
            ['02 miRNA_Pre/miRNA_precursor.xls', 'xls', 'miRNA前体预测详情表', 0],
        ]

        self.AS = [
            [".", "", "可变剪切分析结果目录", 0],
            ['splicing_stats.xls', 'xls', '组内差异可变剪切事件统计表', 0],
            [r".*/all_events_detail_big_table\.txt", "", "结果详情表", 0, "211013"],
            [r".*/event_stats\.file\.txt", "", "差异可变剪切事件统计表", 0, "211014"],
            [r".*/event_type\.file\.txt", "", "可变剪切事件类型统计表", 0, "211015"],
            [r".*/psi_stats\.file\.txt", "", "差异可变剪切模式变化统计表", 0, "211016"],
            [r".*/fromGTF\.(RI|A3SS|A5SS|SE|MXE)\.alter_id\.txt", 'txt', '全部可变剪接事件表', 0, "211017"],
            [r".*/fromGTF\.novelEvents\.(RI|A3SS|A5SS|SE|MXE)\.alter_id\.txt", 'txt', '新发现可变剪接事件表', 0, "211018"],
            [r".*/(RI|A3SS|A5SS|SE|MXE)\.MATS\.JCEC\.alter_id\.psi_info\.txt", 'txt', '差异事件详情表（JCEC）', 0, "211019"],
            [r".*/(RI|A3SS|A5SS|SE|MXE)\.MATS\.JC\.alter_id\.psi_info\.txt", 'txt', '差异事件详情表（JC）', 0, "211020"],
            [r".*/(RI|A3SS|A5SS|SE|MXE)\.MATS\.JCEC\.alter_id\.txt", 'txt', '事件详情表（JCEC）', 0, "211021"],
            [r".*/(RI|A3SS|A5SS|SE|MXE)\.MATS\.JC\.alter_id\.txt", 'txt', '事件详情表（JC）', 0, "211022"],
            [r'.*/.*_splicing_JC_stats\.xls', 'xls', '组内差异可变剪切模式变化统计表(JC)', 0],
            [r'.*/.*_splicing_JCEC_stats\.xls', 'xls', '组内差异可变剪切模式变化统计表(JCEC)', 0],
        ]

        self.SNP_InDel = [
            [".", "", "SNP/InDel分析结果目录", 0],
            ["snp_anno.xls", 'xls', 'SNP分析结果注释详情表', 0],
            ["indel_anno.xls", 'xls', 'InDel分析结果注释详情表', 0],
            ["snp_annotation_statistics.xls", "", "snp分析结果注释详情表格", 0,  "211405"],
            ["snp_transition_tranversion_statistics.xls", "", "SNP类型统计结果表格", 0,  "211406"],
            ["snp_freq_statistics.xls", "", "SNP频率统计结果表格", 0,  "211407"],
            ["snp_depth_statistics.xls", "", "SNP深度统计结果表格", 0,  "211408"],
            ["snp_position_distribution.xls", "", "SNP不同区域布析结果表格", 0,  "211409"],
            ["indel_position_distribution.xls", "", "InDel不同区域布析结果表格", 0,  "211410"],
        ]

        self.Other = [
            ['.', '', '其他中间文件', 0],
            ["Filtered_result", "", "表达量过滤信息", 0],
            ["Sequence_database", "", "序列文件数据库", 0],
        ]

        # self.LncRNA_target_cis = [
        #     [".", "", "lncRNA cis靶基因预测结果", 0],
        # ]
        # self.LncRNA_target_trans = [
        #     [".", "", "lncRNA trans靶基因预测结果", 0],
        # ]


