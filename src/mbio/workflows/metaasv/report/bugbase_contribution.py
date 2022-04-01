# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
import glob
import os
from bson import ObjectId
import re
import types
import pandas as pd
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir,link_file


class BugbaseContributionWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BugbaseContributionWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "bugbase_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "tax_file", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "species_level", "type": "string"},
            {"name": "bugbase_id", "type": "string"},
            {"name": "asv_id", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "top", "type": "int", "default":10},
            {"name": "method", "type": "string", "default":"average"},
            {"name": "normalized_table", "type": "infile", "format": "meta.otu.otu_table"},
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        #self.bugbase_contribution = self.add_tool('meta.bugbase_contribution')
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.contribute_tool = self.add_tool("meta.association_model.contribute")  # 贡献度分析各样本计算
        self.tax_level={"3":"Phylum","4":"Class","5":"Order","6":"Family","7":"Genus"}

    def run_sort_samples(self):
        self.sort_samples.set_options({
            "in_otu_table": self.option("normalized_table"),
            "group_table": self.option("group_table"),
            "method": self.option("method"),
            "get_only_percent": "T"
        })
        self.sort_samples.on("end", self.run_contribute)
        self.sort_samples.run()

    def change_info(self):
        with open(self.sort_samples.output_dir + "/taxa.percents.table.xls","r") as f, open(self.work_dir+"/otu.xls","w") as t:
            data = f.readlines()
            t.write("GeneID" + "\t" + "\t".join(data[0].rstrip().split("\t")[1:]) + "\t" + "Total" +"\n")
            for i in data[1:]:
                tmp = 0
                for x in i.strip().split("\t")[1:]:
                    tmp += float(x)
                t.write(i.strip()+ "\t" +str(tmp) + "\n")

    def run_contribute(self):
        self.change_info()
        self.logger.info(self.option("top"))
        #sel_gene_profile = self.sort_samples.option("out_otu_table")
        options = {
            "taxon_file": self.option("tax_file").prop["path"],
            "function_file": self.option("bugbase_table").prop["path"],
            "gene_profile": self.work_dir+"/otu.xls",
            "tax_level": self.tax_level[self.option("species_level")],
            "fun_level": "category",
            "top_tax": 250,
            "top_fun": 10
        }
        self.contribute_tool.set_options(options)
        # self.select_tools.append(anno_select)
        self.contribute_tool.on('end', self.get_top)
        self.contribute_tool.run()

    def run(self):

        #self.sort_tax_samples.on("end",self.run_bugbase_contribution)
        self.run_sort_samples()
        super(BugbaseContributionWorkflow, self).run()

    def get_top(self):
        fun_tax_file = os.path.join(self.contribute_tool.output_dir, "Function_taxon_abundance.xls")
        new_fun_tax_file = os.path.join(self.work_dir, "Function_taxon_abundance.xls")
        all_function = []
        with open(fun_tax_file) as f, open(new_fun_tax_file,"w") as t:
            data = f.readlines()
            t.write(data[0])
            for i in data[1:]:
                if i.strip().split("\t")[0] not in all_function:
                    all_function.append(i.strip().split("\t")[0])
            for fun in all_function:
                all_dict = []
                for x in data[1:]:
                    if x.strip().split("\t")[0] == fun:
                        all_dict.append({"content": x,"value": float(x.strip().split("\t")[-1])})
                all_dict_sort = sorted(all_dict, key=lambda all_dict: all_dict["value"], reverse=True)
                if int(self.option("top")) >= len(all_dict_sort):
                    for xx in all_dict_sort:
                        t.write(xx["content"])
                else:
                    for xx in all_dict_sort[0:int(self.option("top"))]:
                        t.write(xx["content"])
                    t.write(fun+ "\t" + "others")
                    total = 0
                    for num in range(len(data[0].strip().split("\t")[2:-1])):
                        order = int(num) + 2
                        temp = 0
                        for xxx in all_dict_sort[(int(self.option("top"))+1):]:
                            temp += float(xxx["content"].strip().split("\t")[order])
                        total += temp
                        t.write("\t" + str(temp))
                    t.write("\t" + str(total) + "\n")
        self.set_db()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        if os.path.exists(os.path.join(self.output_dir, "Phenotype_taxon_abundance.xls")):
            os.remove(os.path.join(self.output_dir, "Phenotype_taxon_abundance.xls"))
        os.link(os.path.join(self.work_dir, "Function_taxon_abundance.xls"), os.path.join(self.output_dir, "Phenotype_taxon_abundance.xls"))
        api_contribute = self.api.api("metaasv.bugbase_contribution")
        tax_species_old = os.path.join(self.contribute_tool.output_dir , "Tax_{}_abu.xls".format(self.tax_level[self.option("species_level")]))
        fun_tax_file = os.path.join(self.output_dir , "Phenotype_taxon_abundance.xls")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                main_id = None
        self.logger.info(main_id)
        #api_contribute.add_contribute(main_id, contribute_dir, "core", update_main=False)
        api_contribute.add_contribute_detail(main_id, fun_tax_file,tax_species_old)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "BugBase贡献度分析结果目录", 0,],
            ["Phenotype_taxon_abundance.xls", "xls", "物种-表型贡献度结果表", 0,]
        ])
        super(BugbaseContributionWorkflow, self).end()
