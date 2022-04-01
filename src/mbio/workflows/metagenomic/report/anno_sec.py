# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2018/11/20'

from biocluster.workflow import Workflow
import os


class AnnoSecWorkflow(Workflow):
    """
    分泌蛋白预测
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AnnoSecWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "params", "type": "string"},
            # {"name": "nr_method", "type": "string"},  # best de_unclassied  lca  导表时用
            {"name": "gene_sec_anno", "type": "infile", "format": "ref_rna_v2.common_dir"},
            {"name": "gene_profile", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_nr_anno", "type": "infile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.sec_anno = self.add_tool("annotation.sec_anno_stat")

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.sec_anno.on("end", self.set_db)
        self.run_sec()
        super(AnnoSecWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_path = self.api.api('metagenomic.mg_anno_sec')
        self.logger.info("开始进行导表")
        file_list = ["all_fisher.txt", "bac_fisher.txt", "euk_fisher.txt", "gramneg_fisher.txt", "grampos_fisher.txt",
                     "all_summary.txt", "bac_summary.txt", "euk_summary.txt", "gramneg_summary.txt", "grampos_summary.txt"]
        model_map = {
            "all": "bac_fun",
            "bac": "bac",
            "euk": "fun",
            "grampos": "bac_pos",
            "gramneg": "bac_neg"
        }
        for file_name in file_list:
            str1,str2 = file_name.split("_")
            model = model_map[str1]
            if "fisher" in str2:
                api_path.add_anno_sec_tax(self.option("main_id"), self.sec_anno.output_dir + "/" + file_name, model=model)
            else:
                api_path.add_anno_sec_stat(self.option("main_id"), self.sec_anno.output_dir + "/" + file_name, model=model)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.sec_anno.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "分泌蛋白预测结果目录", 0, "120270"]
        ])
        result_dir.add_regexp_rules([
            [r".*fisher\.txt", "txt", "分泌蛋白物种注释结果", 0, "120271"],
            [r".*summary\.txt", "txt", "分泌蛋白个数统计", 0, "120272"]
        ])
        super(AnnoSecWorkflow, self).end()

    def run_sec(self):
        if not os.path.isdir(os.path.join(self.work_dir, "predict_dir")):
            os.mkdir(os.path.join(self.work_dir, "predict_dir"))
        for file in os.listdir(self.option("gene_sec_anno").path):
            if not file.startswith("signalp"):
                continue
            old_path = os.path.join(self.option('gene_sec_anno').path, file)
            new_path = os.path.join(self.work_dir, "predict_dir", file)
            if os.path.isfile(new_path):
                os.remove(new_path)
            os.link(old_path, new_path)
        self.sec_anno.set_options({
            "predict_dir": os.path.join(self.work_dir, "predict_dir"),
            "reads_profile_table": self.option("gene_profile"),
            "gene_nr_anno": self.option("gene_nr_anno")
        })
        self.sec_anno.run()
