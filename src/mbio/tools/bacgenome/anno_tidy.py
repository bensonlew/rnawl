# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# version 1.0
# last_modify: 2018.03.15

import os
import re
import shutil
import pandas as pd
from collections import Counter
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bacgenome.anno_table_tidy import anno_table_tidy, get_nr_des



class AnnoTidyAgent(Agent):
    """
    对NR/Ref/Swissprot/Pfam的结果进行整理：
    1. 比对的m6结果文件整理，提取需要的列以及按不同数据库修改表头；
    2. 给数据注释表格加上“Location”、“NR- Description”的信息；
    3. 可能生成一些整理的文件
    """

    def __init__(self, parent):
        super(AnnoTidyAgent, self).__init__(parent)
        options = [
            {"name": "gene_gff", "type": "infile", "format": "gene_structure.gff3"},  # 基因预测结果
            {"name": "anno_nr", "type": "infile", "format": "align.blast.blast_table"},  # NR注释结果
            {"name": "anno_go", "type": "infile", "format": "annotation.go.level2"},  # GO注释结果
            {"name": "anno_swissprot", "type": "infile", "format": "align.blast.blast_table"},  # Swiss-prot注释结果
            {"name": "anno_ref", "type": "infile", "format": "align.blast.blast_table"},  # 参考基因组注释结果
            {"name": "anno_pfam", "type": "infile", "format": "meta_genomic.hmmscan_table"},  # Pfam注释结果
            {"name": "anno_cog", "type": "infile", "format": "sequence.profile_table"},  # COG注释结果
            {"name": "kegg_xml", "type": "infile", "format": "align.blast.blast_xml"},  # KEGGblast结果，用于生成pathway通路图
            {"name": "anno_kegg", "type": "infile", "format": "sequence.profile_table"},  # KEGG注释结果
            {"name": "anno_cazy", "type": "infile", "format": "meta_genomic.hmmscan_table"},  # CAZy注释结果
            {"name": "anno_antismash", "type": "infile", "format": "sequence.profile_table"},  # Antismash注释结果
            {"name": "summary", "type": "string", 'default': "F"},  # 是否作数据库汇总表，默认不做
            {"name": "sample", "type": "string", 'default': ""},
            {"name": "tidy_nr", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_go", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_swissprot", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_ref", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_pfam", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_cog", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_kegg", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_cazy", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_antismash", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_summary", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "pathway_img", "type": "string"},  # pathway通路图的dir文件路径,里面可能含有多个文件夹对应多个染色体或质粒
            {"name": "analysis", "type": "string", "default":""}  #complete  uncomplete


        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("gene_gff").is_set:
            raise OptionError("必须设置参数gene_gff，提供基因位置信息！", code="31403001")
        if not self.option("anno_nr").is_set and not self.option("anno_go").is_set and not self.option(
                "anno_swissprot").is_set and not self.option("anno_ref").is_set and not self.option(
            "anno_pfam").is_set and not self.option("anno_cog").is_set and not self.option(
            "anno_kegg").is_set and not self.option("anno_cazy").is_set and not self.option(
            "anno_antismash").is_set:
            raise OptionError("请至少提供一种注释结果文件！", code="31403002")
        if self.option("anno_cog").is_set or self.option("anno_kegg").is_set or self.option(
                "anno_go").is_set or self.option("anno_pfam").is_set or self.option("anno_cazy").is_set or self.option(
            "anno_antismash").is_set:
            if not self.option("anno_nr").is_set:
                raise OptionError("提供该注释结果的同时，必须提供NR的注释结果！", code="31403003")

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(AnnoTidyAgent, self).end()


class AnnoTidyTool(Tool):
    def __init__(self, config):
        super(AnnoTidyTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/Python/bin', LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/program/Python/lib')
        self.kegg_img = self.config.PACKAGE_DIR + '/annotation/mg_annotation/kegg_pathway_img_v94.py'

    def run_anno_tidy(self):
        gene = self.option("gene_gff").prop["path"]
        if self.option("anno_nr").is_set:
            anno_table_tidy(self.option("anno_nr").prop["path"], "nr", gene, self.output_dir + "/anno_nr.xls")
            self.file_split(self.output_dir + "/anno_nr.xls","NR")
        if self.option("anno_swissprot").is_set:
            anno_table_tidy(self.option("anno_swissprot").prop["path"], "swissprot", gene,
                            self.output_dir + "/anno_swissprot.xls")
            self.file_split(self.output_dir + "/anno_swissprot.xls","Swissprot")
        if self.option("anno_ref").is_set:
            anno_table_tidy(self.option("anno_ref").prop["path"], "nr", gene, self.output_dir + "/anno_ref.xls")
            self.logger.info(self.output_dir + "/anno_ref.xls")
            self.option("tidy_ref", self.output_dir + "/anno_ref.xls")

    def run_get_nr_des(self):
        nr = self.output_dir + "/anno_nr.xls"
        if self.option("anno_cog").is_set:
            get_nr_des(self.option("anno_cog").prop["path"], nr, "anno_cog.xls")
            cog = pd.read_table(self.work_dir + "/anno_cog.xls", sep='\t', header=0)
            cog.rename(columns={"NOG": "COG ID", "Function": "COG Type", "NOG_description": "COG Description",
                                "Fun_description": "Type Description"}, inplace=True)
            cog_summar = cog.ix[:, ["Gene ID","COG Type","Type Description", "Category"]]
            cog_summar1 = cog_summar.drop("COG Type", axis=1).join(cog_summar["COG Type"].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename("Type"))
            cog_summar2 = cog_summar.drop("Category", axis=1).join(cog_summar["Category"].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename("Category"))
            cog_summar3 = cog_summar.drop("Type Description", axis=1).join(cog_summar["Type Description"].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename("Type Description"))
            tmp = pd.DataFrame({"Category":list(cog_summar2["Category"]),"Type":list(cog_summar1["Type"]),"Type Description":list(cog_summar3["Type Description"])})
            type_nu = (Counter(cog_summar1["Type"]))
            cog_tmp = tmp.drop_duplicates(subset=["Type"], keep='first')
            cog_tmp.index = cog_tmp["Type"]
            nu_list = []
            type_list = []
            for i in list(cog_tmp["Type"]):
                type = "[" + i + "] " + cog_tmp.loc[i,"Type Description"]
                nu = type_nu[i]
                type_list.append(type)
                nu_list.append(nu)
            cog_tmp = cog_tmp[["Category"]].copy()
            cog_tmp["Type"] = type_list
            cog_tmp["Genes"] = nu_list
            cog_tmp.sort_values('Category', inplace=True)
            #cog_summary = cog_tmp.ix[:, ["Category","Type","Gene No."]]
            if not os.path.exists(self.output_dir + "/COG"):
                os.makedirs(self.output_dir + "/COG")
            cog_tmp.to_csv(self.output_dir + "/COG/" + self.option("sample") + "_cog_summary.xls", sep='\t', header=True, index=False)
            cog.to_csv(self.output_dir + "/COG/" + self.option("sample") + "_cog_anno.xls", sep='\t', header=True, index=False)
            self.option("tidy_cog", self.output_dir + "/COG/" + self.option("sample") + "_cog_anno.xls")
        if self.option("anno_kegg").is_set:
            if not os.path.exists(self.output_dir + "/KEGG"):
                os.makedirs(self.output_dir + "/KEGG")
            get_nr_des(self.option("anno_kegg").prop["path"], nr, self.output_dir + "/KEGG/" + self.option("sample") + "_kegg_anno.xls")
            kegg = pd.read_table(self.output_dir + "/KEGG/" + self.option("sample") + "_kegg_anno.xls", sep='\t', header=0)
            kegg = kegg.ix[:, ["Location", "Pathway", "KO"]]
            pathway = kegg.drop("Pathway", axis=1).join(
                kegg["Pathway"].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename("Pathway"))
            self.option("tidy_kegg", self.output_dir + "/KEGG/" + self.option("sample") + "_kegg_anno.xls")
            #self.option("pathway_img", self.output_dir + "/KEGG")
            pathway_ko_dir = self.output_dir + "/KEGG/" + self.option("sample") + "_kegg_pathway_img"
            pathway_ko = self.work_dir + "/scaffold_pathway_img"
            with open(pathway_ko, 'w') as tab:
                tab.write("pathway\tKO\n")
                for k in list(set(pathway["Pathway"])):
                    if k == "-":
                        pass
                    else:
                        ko_list = ";".join(list(set(pathway[pathway["Pathway"] == k]["KO"])))
                        tab.write(k + "\t" + ko_list + "\n")
            self.run_pathway_img(self.option("kegg_xml").prop["path"], self.work_dir, pathway_ko,"pathway", "KO", "pathway_img")
            if os.path.exists(pathway_ko_dir):
                shutil.rmtree(pathway_ko_dir)
            os.rename(self.work_dir + "/pathway_img", pathway_ko_dir)
            """
            if re.match("^\Scaffold.*", list(set(pathway['Location']))[0]):
                pathway_ko = self.work_dir + "/scaffold_pathway_img"
                with open(pathway_ko, 'w') as tab:
                    tab.write("pathway\tKO\n")
                    for k in list(set(pathway["Pathway"])):
                        if k == "-":
                            pass
                        else:
                            ko_list = ";".join(list(set(pathway[pathway["Pathway"] == k]["KO"])))
                            tab.write(k + "\t" + ko_list + "\n")
                self.run_pathway_img(self.option("kegg_xml").prop["path"], self.work_dir, pathway_ko,
                                     "pathway", "KO", "pathway_img")
                os.rename(self.work_dir + "/pathway_img", self.output_dir + "/KEGG/scaffold_pathway_img")
            else:
                for loc in list(set(pathway['Location'])):
                    pathway_ko = self.work_dir + "/" + loc
                    with open(pathway_ko, 'w') as tab:
                        tab.write("pathway\tKO\n")
                        pa = pathway[pathway["Location"] == loc]
                        for k in list(set(pa["Pathway"])):
                            if k == "-":
                                pass
                            else:
                                ko_list = ";".join(pa[pa["Pathway"] == k]["KO"])
                                tab.write(k + "\t" + ko_list + "\n")
                    self.run_pathway_img(self.option("kegg_xml").prop["path"], self.work_dir, pathway_ko,
                                         "pathway", "KO", loc.lower())
                    os.rename(self.work_dir + "/pathway_img", self.output_dir + "/KEGG/" + loc + "_pathway_img")
            self.logger.info(self.output_dir + "/anno_kegg.xls")  # 细菌单个质粒和染色体的pathway_img
            """

        if self.option("anno_go").is_set:
            go = pd.read_table(self.option("anno_go").prop["path"], sep='\t', header=0)
            go = go.ix[:,["GO ID", "Seq List","GO Term"]]   ##GO ID (Lev2)改成GO ID
            go.columns = ["GO", "Seq List","GO Term"]
            go = go.drop("Seq List", axis=1).join(go["Seq List"].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename("Gene ID"))
            gene_go = []
            go_desc = []
            for i in sorted(list(set(go["Gene ID"]))):
                gene_go.append(";".join(list(go[go["Gene ID"]==i]["GO"])))
                go_desc.append(";;".join(list(go[go["Gene ID"]==i]["GO Term"])))
            go_anno = pd.DataFrame(columns = ["Seq_id", "GO","GO Term"])
            go_anno["Seq_id"] = sorted(list(set(go["Gene ID"])))
            go_anno["GO"] = gene_go
            go_anno["GO Term"] = go_desc
            go_anno.to_csv("go_tmp.xls", sep="\t", index=None)
            if not os.path.exists(self.output_dir + "/GO"):
                os.makedirs(self.output_dir + "/GO")
            get_nr_des(self.work_dir + "/go_tmp.xls", nr, self.output_dir + "/GO/" + self.option("sample") + "_go_anno.xls")
            if os.path.exists(self.output_dir + "/GO/" + self.option("sample") + "_go_statistics.xls"):    ##20190510  _go_statistics.xls
                os.remove(self.output_dir + "/GO/" + self.option("sample") + "_go_statistics.xls")
            os.link(self.option("anno_go").prop["path"], self.output_dir + "/GO/" + self.option("sample") + "_go_statistics.xls")
            self.option("tidy_go", self.output_dir + "/GO/" + self.option("sample") + "_go_anno.xls")
        if self.option("anno_pfam").is_set:
            get_nr_des(self.option("anno_pfam").prop["path"], nr, "anno_pfam.xls")
            pfam = pd.read_table(self.work_dir + "/anno_pfam.xls", sep='\t', header=0, index_col=0)
            del pfam["Protein_id"]
            pfam = pfam.rename(columns={'ProteinStart': 'Start', 'Protein': 'End'})
            pfam.to_csv(self.output_dir + "/anno_pfam.xls", sep="\t")
            self.file_split(self.output_dir + "/anno_pfam.xls","Pfam")
        if self.option("anno_cazy").is_set:
            if not os.path.exists(self.output_dir + "/CAZy"):
                os.makedirs(self.output_dir + "/CAZy")
            get_nr_des(self.option("anno_cazy").prop["path"], nr, self.output_dir + "/CAZy/" + self.option("sample") + "_anno_cazy.xls")
            self.option("tidy_cazy", self.output_dir + "/CAZy/" + self.option("sample") + "_anno_cazy.xls")

    def run_anno_antismash(self):
        antismash = pd.read_table(self.option("anno_antismash").prop["path"], sep='\t', header=0)
        antismash = antismash.ix[:, ["Cluster ID", "Genes"]]
        antis = antismash.drop("Genes", axis=1).join(
            antismash["Genes"].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename("Gene ID"))
        nr = pd.read_table(self.output_dir + "/anno_nr.xls", sep='\t', header=0)
        nr = nr.ix[:, ["Gene ID", "Location", "Hit-Description"]]
        nr.columns = ["Gene ID", "Location", "Gene Description"]
        nr = nr.drop_duplicates(subset=["Gene ID", "Location"], keep='first')
        table = antis.merge(nr, on="Gene ID", how="left")
        pfam = pd.read_table(self.output_dir + "/anno_pfam.xls", sep='\t', header=0)
        ann_gene = pfam[["Gene ID", "Location"]].drop_duplicates(subset=["Gene ID", "Location"], keep='first')
        gene_li = list(set(ann_gene['Gene ID']))
        pfam_id = []
        for ge in gene_li:
            pfam_g = ";".join(pfam[pfam["Gene ID"] == ge]["Pfam_id"])
            pfam_id.append(pfam_g)
        ann_gene["Pfam_id"] = pfam_id
        ann_gene.index = ann_gene["Gene ID"]
        table = table.merge(ann_gene, on="Gene ID", how="left")
        table = table.fillna("-")
        table.index = table["Gene ID"]
        del table["Gene ID"]
        table.to_csv(self.output_dir + "/anno_antismash.xls", sep='\t')
        self.logger.info(self.output_dir + "/anno_antismash.xls")
        self.option("tidy_antismash", self.output_dir + "/anno_antismash.xls")

    def run_summary(self):
        gene = pd.read_table(self.option("gene_gff").prop["path"], sep='\t', header=0)
        gene = gene.ix[:, ["Gene ID", "Sequence id", "Strand", "Start", "End", "Gene Length(bp)"]]
        gene["Location"] = gene["Sequence id"].str.split("_", expand=True)[0]
        del gene["Sequence id"]
        nr = pd.read_table(self.output_dir + "/anno_nr.xls", sep='\t', header=0)
        nr = nr.drop_duplicates(subset=["Gene ID"], keep='first')
        nr = nr.ix[:, ["Gene ID", "Hit-Description"]]
        nr.columns = ["Gene ID", "NR Description"]
        swiss = pd.read_table(self.output_dir + "/anno_swissprot.xls", sep='\t', header=0)
        swiss = swiss.drop_duplicates(subset=["Gene ID"], keep='first')
        swiss = swiss.ix[:, ["Gene ID", "Hit-Description"]]
        swiss.columns = ["Gene ID", "Swiss-Prot Description"]
        pfam = pd.read_table(self.output_dir + "/anno_pfam.xls", sep='\t', header=0)
        ann_gene = pfam.drop_duplicates(subset=["Gene ID", "Location"], keep='first')
        gene_li = list(set(ann_gene['Gene ID']))
        pfam_id = []
        pfam_domain = []
        pfam_domain_desc = []
        for ge in gene_li:
            pfam_g = ";".join(pfam[pfam["Gene ID"] == ge]["Pfam_id"])
            pfam_d_g =  ";".join(pfam[pfam["Gene ID"] == ge]["Domain"])
            pfam_d_d_g = ";".join(pfam[pfam["Gene ID"] == ge]["DomainDescription"])
            pfam_id.append(pfam_g)
            pfam_domain.append(pfam_d_g)
            pfam_domain_desc.append(pfam_d_d_g)
        ann_gene["Pfam_id"] = pfam_id
        ann_gene['Domain'] = pfam_domain
        ann_gene["DomainDescription"] = pfam_domain_desc
        ann_gene['Gene ID'] = gene_li
        ann_gene.index = ann_gene["Gene ID"]
        Pfam = ann_gene.ix[:, ["Gene ID", "Pfam_id","Domain", "DomainDescription"]]

        Pfam.columns = ["Gene ID", "Pfam_id","Pfam Domain","Domain Description"]
        cog = pd.read_table(self.output_dir + "/COG/" + self.option("sample") + "_cog_anno.xls", sep='\t', header=0)
        cog = cog.ix[:, ["Gene ID", "COG ID","COG Description", "COG Type"]]
        kegg = pd.read_table(self.output_dir + "/KEGG/" + self.option("sample") + "_kegg_anno.xls", sep='\t', header=0)
        kegg = kegg.ix[:, ["Gene ID", "Gene", "KO","Definition","Pathway"]]
        kegg.columns = ["Gene ID", "Gene Name", "KO ID","KO Description","Pathway"]
        go = pd.read_table(self.output_dir + "/GO/" + self.option("sample") + "_go_anno.xls", sep='\t', header=0)
        go = go.ix[:, ["Gene ID", "GO", "GO Term"]]  ##20190510 add GO desc
        go.columns = ["Gene ID", "GO ID","GO Description"]  ##20190510
        a = gene.merge(nr, on='Gene ID', how='left')
        a = a.merge(swiss, on='Gene ID', how='left')
        a = a.merge(cog, on='Gene ID', how='left')  #cog type在第十一列
        a = a.merge(Pfam, on='Gene ID', how='left')
        a = a.merge(kegg, on='Gene ID', how='left')
        a = a.merge(go, on='Gene ID', how='left')
        a = a.fillna("-")
        if not os.path.exists(self.output_dir + "/Summary"):
            os.mkdir(self.output_dir + "/Summary")
        a.to_csv(self.output_dir + "/Summary/" + self.option("sample") + "_anno_summary.xls", sep='\t', header=True, index=False)
        #self.logger.info(self.output_dir + "/anno_summary.xls")
        self.option("tidy_summary", self.output_dir + "/Summary/" + self.option("sample") + "_anno_summary.xls")

    def run_pathway_img(self, kegg_xml, output_dir, pathway_ko, pathway, ko, loc):
        cmd = "{}/program/Python/bin/python {} -i {} -o {} -p {} -ko {} -KO {} -png_file True -bac True".format("", self.kegg_img, kegg_xml,
                                                                                       output_dir, pathway_ko, pathway,
                                                                                       ko)
        cmd += '  -html {}/database/Annotation/all/KEGG/version_202007_meta/html/'.format(self.config.SOFTWARE_DIR)   #-png_file True
        self.logger.info(cmd)
        command = self.add_command(loc, cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("output_kegg_pathway_img succeed")
        else:
            self.set_error("output kegg pathway img failed", code="31403001")

    def file_split(self,anno,type):
        f1 = pd.read_table(anno, sep='\t', header=0)
        path = self.output_dir + "/" + type
        if not os.path.exists(path):
            os.makedirs(path)
        if re.match("^[Ss]caffold[0-9]+$", list(f1["Location"])[1]) or self.option("analysis") in ["uncomplete"] :
            f1.to_csv(path + "/" + self.option("sample") + "_anno_" + type.lower() + ".xls", sep="\t",index=None)
            self.option("tidy_" + type.lower(), path + "/" + self.option("sample") + "_anno_" + type.lower() + ".xls")
        else:
            f1.to_csv(path + "/" + self.option("sample") + "_whole_genome_anno_" + type.lower() + ".xls", sep="\t",index=None)
            self.option("tidy_" + type.lower(),path + "/" + self.option("sample") + "_whole_genome_anno_" + type.lower() + ".xls")
            for loc in list(set(f1['Location'])):
                ta = f1[f1["Location"] == loc]
                location = loc
                ta.to_csv(path + "/" + self.option("sample") + "_" + location + "_" + type.lower() + ".xls",sep="\t", index=None)

    def run(self):
        super(AnnoTidyTool, self).run()
        self.run_anno_tidy()
        self.run_get_nr_des()
        if self.option("anno_antismash").is_set:
            self.run_anno_antismash()
        if self.option("summary") != "F":
            self.run_summary()
        self.end()
