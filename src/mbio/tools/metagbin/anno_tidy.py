# -*- coding: utf-8 -*-
# __author__ = 'gaohao'


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
    注释统计总览表
    """

    def __init__(self, parent):
        super(AnnoTidyAgent, self).__init__(parent)
        options = [
            {"name": "gene_gff", "type": "infile", "format": "gene_structure.gff3"},  # 基因预测结果
            {"name": "anno_nr", "type": "infile", "format": "align.blast.blast_table"},  # NR注释结果
            {"name": "anno_cog", "type": "infile", "format": "sequence.profile_table"},  # COG注释结果
            {"name": "anno_kegg", "type": "infile", "format": "sequence.profile_table"},  # KEGG注释结果
            {"name": "anno_cazy", "type": "infile", "format": "meta_genomic.hmmscan_table"},  # CAZy注释结果
            {"name": "kegg_xml", "type": "infile", "format": "align.blast.blast_xml"},  # KEGGblast结果，用于生成pathway通路图
            {"name": "anno_card", "type": "infile", "format": "sequence.profile_table"},  # CAZy注释结果
            {"name": "summary", "type": "string", 'default': "F"},  # 是否作数据库汇总表，默认不做
            {"name": "sample", "type": "string", 'default': ""},
            {"name": "tidy_nr", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_card", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_cog", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_kegg", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_cazy", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_summary", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "version", "type": "string", 'default': "new"}, ## 兼容新老版本
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("gene_gff").is_set:
            raise OptionError("必须设置参数gene_gff，提供基因位置信息！")
        if not self.option("anno_nr").is_set and not self.option(
            "anno_card").is_set and not self.option("anno_cog").is_set and not self.option(
            "anno_kegg").is_set and not self.option("anno_cazy").is_set:
            raise OptionError("请至少提供一种注释结果文件！")
        if self.option("anno_cog").is_set or self.option("anno_kegg").is_set  or self.option("anno_card").is_set or self.option("anno_cazy").is_set :
            if not self.option("anno_nr").is_set:
                raise OptionError("提供该注释结果的同时，必须提供NR的注释结果！")

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(AnnoTidyAgent, self).end()


class AnnoTidyTool(Tool):
    def __init__(self, config):
        super(AnnoTidyTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/Python/bin', LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/program/Python/lib')
        if self.option("version") in ["new"]:
            self.kegg_img = self.config.PACKAGE_DIR + '/annotation/mg_annotation/kegg_pathway_img_v94.py'
        else:
            self.kegg_img = self.config.PACKAGE_DIR + '/annotation/mg_annotation/kegg_pathway_img.py'

    def run_anno_tidy(self):
        gene = self.option("gene_gff").prop["path"]
        if self.option("anno_nr").is_set:
            anno_table_tidy(self.option("anno_nr").prop["path"], "nr", gene, self.output_dir + "/anno_nr.xls")
            self.file_split(self.output_dir + "/anno_nr.xls","NR")

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
        if self.option("anno_cazy").is_set:
            if not os.path.exists(self.output_dir + "/CAZy"):
                os.makedirs(self.output_dir + "/CAZy")
            get_nr_des(self.option("anno_cazy").prop["path"], nr, self.output_dir + "/CAZy/" + self.option("sample") + "_anno_cazy.xls")
            self.option("tidy_cazy", self.output_dir + "/CAZy/" + self.option("sample") + "_anno_cazy.xls")


    def run_summary(self):
        gene = pd.read_table(self.option("gene_gff").prop["path"], sep='\t', header=0)
        gene = gene.ix[:, ["Gene ID", "Sequence id", "Strand", "Start", "End", "Gene Length(bp)"]]
        gene["Location"] = gene["Sequence id"].str.split("_", expand=True)[0]
        del gene["Sequence id"]
        nr = pd.read_table(self.output_dir + "/anno_nr.xls", sep='\t', header=0)
        nr = nr.drop_duplicates(subset=["Gene ID"], keep='first')
        nr_num =nr[nr['Hit'] != "-"]['Gene ID'].count()
        nr = nr.ix[:, ["Gene ID", "Hit-Description"]]
        nr.columns = ["Gene ID", "NR Description"]
        card_num = 0
        if self.option("anno_card").is_set:
            card = pd.read_table(self.option("anno_card").prop["path"], sep='\t', header=0)
            card_num=card['Gene ID'].count()
            if self.option("version") in ['new']:
                card = card.ix[:, ["Gene ID", "ARO_name","ARO_description","Drug_class", 'Resistance_mechanism']]
            else:
                card = card.ix[:, ["Gene ID", "ARO_name","ARO_description","ARO_category"]]
        cog = pd.read_table(self.output_dir + "/COG/" + self.option("sample") + "_cog_anno.xls", sep='\t', header=0)
        cog_num = cog['Gene ID'].count()
        cog = cog.ix[:, ["Gene ID", "COG ID", "COG Type"]]
        kegg = pd.read_table(self.output_dir + "/KEGG/" + self.option("sample") + "_kegg_anno.xls", sep='\t', header=0)
        kegg_num = kegg['Gene ID'].count()
        kegg = kegg.ix[:, ["Gene ID", "Gene", "KO","Definition"]]
        kegg.columns = ["Gene ID", "Gene Name", "KO ID","KO Description"]
        cazy = pd.read_table(self.output_dir + "/CAZy/" + self.option("sample") + "_anno_cazy.xls", sep='\t', header=0)
        cazy_num = cazy['Gene ID'].count()
        cazy = cazy.ix[:, ["Gene ID", "Family","Class","Class_description"]]
        cazy.columns = ["Gene ID", "Family","Class","Class description"]
        a = gene.merge(nr, on='Gene ID', how='left')
        a = a.merge(cog, on='Gene ID', how='left')
        a = a.merge(kegg, on='Gene ID', how='left')
        a = a.merge(cazy, on='Gene ID', how='left')
        if self.option("anno_card").is_set:
            a = a.merge(card, on='Gene ID', how='left')
        a = a.fillna("-")
        if not os.path.exists(self.output_dir + "/Summary"):
            os.mkdir(self.output_dir + "/Summary")
        a.to_csv(self.output_dir + "/Summary/" + self.option("sample") + "_anno_summary.xls", sep='\t', header=True, index=False)
        with open (self.output_dir + "/Summary/" + self.option("sample") + "_anno_stat.xls",'w') as f:
            f.write("Genome ID" + '\t' + "NR" + '\t' + "COG" + '\t' + "KEGG" + '\t' + "CAZY" + '\t' + "CARD" + '\n')
            f.write(self.option("sample") + '\t' + str(nr_num) + '\t' + str(cog_num) + '\t' + str(kegg_num) + '\t' + str(cazy_num) + '\t' + str(card_num) + '\n')
        self.option("tidy_summary", self.output_dir + "/Summary/" + self.option("sample") + "_anno_summary.xls")

    def run_pathway_img(self, kegg_xml, output_dir, pathway_ko, pathway, ko, loc):
        cmd = "{}/miniconda2/bin/python {} -i {} -o {} -p {} -ko {} -KO {}".format("", self.kegg_img, kegg_xml,
                                                                                       output_dir, pathway_ko, pathway,ko)
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
        if re.match("^[Ss]caffold[0-9]+$", list(f1["Location"])[1]):
            f1.to_csv(path + "/" + self.option("sample") + "_anno_" + type.lower() + ".xls", sep="\t",index=None)
            self.option("tidy_" + type.lower(), path + "/" + self.option("sample") + "_anno_" + type.lower() + ".xls")

    def run(self):
        super(AnnoTidyTool, self).run()
        self.run_anno_tidy()
        self.run_get_nr_des()
        if self.option("summary") != "F":
            self.run_summary()
        self.end()