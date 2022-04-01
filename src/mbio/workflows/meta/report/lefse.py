# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

"""lefse分析"""

from biocluster.workflow import Workflow
import os
import pandas as pd
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class LefseWorkflow(Workflow):
    """
    报告中调用lefse分析时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(LefseWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", 'format': "meta.otu.otu_table"},
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_detail", "type": "string"},
            {"name": "second_group_detail", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "lefse_type", "type": "string"},
            {"name": "lda_filter", "type": "float", "default": 2.0},
            {"name": "strict", "type": "int", "default": 0},
            {"name": "group_name", "type": "string"},
            {"name": "lefse_id", "type": "string"},
            {"name": "start_level", "type": "int", "default": 3},
            {"name": "end_level", "type": "int", "default": 7},
            {"name": "normalization", "type": "int", "default": 1},  # 是否标准化 1，0
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.lefse = self.add_tool("statistical.lefse")
        self.logger.info(self.option("group_name"))

    # def get_percent_abund(self):
    #     in_file = self.option('otu_file').path
    #     data = pd.read_table(in_file,sep='\t',header=0)
    #     samples = list(data.columns)
    #     samples.remove('OTU ID')
    #     samples.remove('taxonomy')
    #
    #     for s in samples:
    #         tmp_sum = float(data[s].sum())
    #         data[s] = data[s].apply(lambda x: x/tmp_sum)
    #
    #     out_file = in_file + '_percent'
    #     data.to_csv(out_file, sep='\t',index=False)
    #     self.percent_file = out_file



    def run_lefse(self):
        options = {
            "lefse_type":self.option("lefse_type"),
            "lefse_input": self.option("otu_file"),
            #"lefse_input" : self.percent_file,
            "lefse_group": self.option("group_file"),
            "lda_filter": self.option("lda_filter"),
            "strict": self.option("strict"),
            "lefse_gname": self.option("group_name"),
            "start_level": self.option("start_level"),
            "end_level": self.option("end_level"),
            "percent_abund": "true"
        }
        if not self.option('normalization'):
            options['percent_abund'] = 'false'

        self.lefse.set_options(options)
        self.lefse.on("end", self.set_db)
        self.output_dir = self.lefse.output_dir
        self.lefse.run()

    def end(self):
        if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "LEfSe差异分析结果目录", 0, "110146"],
            # ["./lefse_LDA.cladogram.png", "png", "LEfSe分析cladogram结果图片"],
            # ["./lefse_LDA.png", "png", "LEfSe分析LDA图片"],
            ["./lefse_LDA.xls", "xls", "LEfSe分析lda数据表", 0, "110147"],
            ["./LEfSe多级物种层级树图.pdf", "pdf", "LEfSe多级物种层级树图", 0, ""],
            ["./LDA判别结果图.pdf", "pdf", "LDA判别结果图", 0, ""]
        ])
        super(LefseWorkflow, self).end()

    def set_db(self):
        """
        保存两组比较分析的结果表保存到mongo数据库中
        """
        api_lefse = self.api.stat_test
        lefse_path = self.output_dir + '/lefse_LDA.xls'
        # lda_png_path = self.output_dir + '/lefse_LDA.png'
        # lda_cladogram_path = self.output_dir + '/lefse_LDA.cladogram.png'
        if not os.path.isfile(lefse_path):
            self.logger.error("找不到报告文件:{}".format(lefse_path))
            self.set_error("找不到报告文件", code="12702101")
        # if not os.path.isfile(lda_png_path):
        #     raise Exception("找不到报告文件:{}".format(lda_png_path))
        # if not os.path.isfile(lda_cladogram_path):
        #     raise Exception("找不到报告文件:{}".format(lda_cladogram_path))
        api_lefse.add_species_difference_lefse_detail(file_path=lefse_path, table_id=self.option("lefse_id"))
        # api_lefse.update_species_difference_lefse(lda_png_path, lda_cladogram_path, self.option("lefse_id"))
        #self.end()
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("lefse_id"), "sg_species_difference_lefse")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("lefse_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "species_lefse_analyse",
                "interaction": 1,
                "main_table": "sg_species_difference_lefse",
            })
            self.figsave.run()
        else:
            self.end()

    def run(self):
        #self.get_percent_abund()
        self.run_lefse()
        super(LefseWorkflow, self).run()
