# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
import glob
import os
from bson import ObjectId
import re
from biocluster.core.exceptions import OptionError


class BugbaseWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BugbaseWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_table", "type": "infile","format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "group_table", "type": "infile","format": "meta.otu.group_table"},
            {"name": "otu_fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "asv_id", "type": 'string'},
            {"name": "group_detail", "type": "string"},
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.bugbase = self.add_tool('metaasv.bugbase_predict')
        self.sort_tax_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False


    def run_tax_sort_samples(self):
        abund_table = self.option("otu_table")
        self.sort_tax_samples.set_options({
            "in_otu_table": abund_table,
            "group_table": self.option('group_table')
        })
        self.sort_tax_samples.run()


    def run_bugbase(self):
        mapping_file = self.work_dir + "/mapping_file"
        with open(self.option("group_table").prop["path"],"r") as v,open(mapping_file,"w") as t:
            data = v.readlines()
            all_group = []
            for i in data[1:]:
                if i.strip().split("\t")[1] in all_group:
                    pass
                else:
                    all_group.append(i.strip().split("\t")[1])
            t.write("#SampleID" + "\t" + "BODY_SITE" + "\n")
            if len(all_group) > 1:
                for i in data[1:]:
                    t.write(i)
            else:
                t.write(data[1].strip().split("\t")[0] + "\t" + data[1].strip().split("\t")[1]+"tmp" + "\n")
                for i in data[2:]:
                    t.write(i)
        tax_abund_table =  self.sort_tax_samples.option("out_otu_table").prop['path']
        self.bugbase.set_options({
            "otu_table": tax_abund_table,
            "otu_fasta" : self.option("otu_fasta"),
            "mapping_file": mapping_file,
            })
        self.bugbase.on("end", self.set_db)
        self.bugbase.run()

    def run(self):

        self.sort_tax_samples.on("end",self.run_bugbase)
        self.run_tax_sort_samples()
        super(BugbaseWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.logger.info("正在写入mongo数据库")
        self.api_bugbase = self.api.api("metaasv.bugbase")
        prediction_path = self.bugbase.output_dir+"/bugbase_predictions.txt"
        threshold_path = self.bugbase.output_dir+"/contributing_asv.txt"
        normalized_path = self.bugbase.output_dir+"/16s_normalized_asv.txt"
        self.api_bugbase.add_detail(prediction_file=prediction_path,threshold_file=threshold_path,normalized_file=normalized_path,table_id=self.option("main_id"))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.bugbase.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "BugBase表型预测结果目录", 0],
            ["./16s_normalized_asv.txt", "txt", "Bugbase表型预测标准化ASV丰度表", 0],
            ["./bugbase_predictions.txt", "txt", "Bugbase表型预测结果表", 0],
            ["./bugbase_thresholds.txt", "txt", "Bugbase表型预测阈值表", 0],
            ["./bugbase_variances.txt", "txt", "Bugbase表型预测不同阈值表型相对丰度表", 0],
            ["./contributing_asv.txt", "txt", "Bugbase表型预测贡献度结果表", 0],
        ])
        super(BugbaseWorkflow, self).end()