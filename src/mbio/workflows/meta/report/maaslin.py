# -*- coding: utf-8 -*-


from biocluster.workflow import Workflow
import glob
import os
from bson import ObjectId
import re
from biocluster.core.exceptions import OptionError
import pandas as pd
from mbio.packages.meta.common_function import envname_restore
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class MaaslinWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MaaslinWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_table", "type": "infile","format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "level", "type": "int", "default": 9},
            {"name": "group_table", "type": "infile","format": "meta.otu.group_table"},
            {"name": "envtable", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "env_labs", "type": "string", "default": ""},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "otu_id", "type": 'string'},
            {"name": "group_detail", "type": "string"},
            {"name": "env_id","type":"string"}
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.maaslin = self.add_tool('meta.maaslin')

        self.sort_tax_samples = self.add_tool("meta.otu.sort_samples_mg")


    def run_tax_sort_samples(self):
        abund_table = self.option("otu_table").path
        data = pd.read_table(abund_table,sep='\t',header=0)
        pat = re.compile(';\s*')
        data['OTU ID'] = data['OTU ID'].apply(lambda x: pat.split(str(x))[-1])
        data.to_csv(abund_table+"_new", sep='\t',index=False)

        self.sort_tax_samples.set_options({
            "in_otu_table": abund_table+"_new",
            "group_table": self.option('group_table'),
        })
        self.sort_tax_samples.run()



    def run_maaslin(self):

        env_table = self.option("envtable")
        tax_abund_table =  self.sort_tax_samples.option("out_otu_table").prop['path']
        self.maaslin.set_options({
            "taxon_table": tax_abund_table,
            "env_table": env_table,
            })
        self.maaslin.on("end", self.set_db)
        self.maaslin.run()


    def run(self):

        self.sort_tax_samples.on("end",self.run_maaslin)
        self.run_tax_sort_samples()
        super(MaaslinWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.logger.info("正在写入mongo数据库")
        self.api_maaslin = self.api.api("maaslin")
        data_path = self.maaslin.output_dir+"/Maaslin.txt"
        site_path = self.maaslin.output_dir+"/Maaslin.xls"
        line_path = self.maaslin.output_dir + "/Maaslin.message.xls"

        self.api_maaslin.add_site_detail(file_path=site_path,table_id=self.option("main_id"))
        self.api_maaslin.add_data_detail(file_path=data_path,table_id=self.option("main_id"))
        self.api_maaslin.add_line_detail(file_path=line_path,table_id=self.option("main_id"))
        self.output_dir = self.maaslin.output_dir
        #self.end()
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "sg_maaslin")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "maaslin",
                "interaction": 1,
                "main_table": "sg_maaslin",
            })
            self.figsave.run()
        else:
            self.end()

    @envname_restore
    def end(self):
        if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "MaAsLin分析结果文件",0,"110251"],
            ["./Maaslin.xls", "xls", "散点图数据",0,"110252"],
            ["./Maaslin.txt", "txt", "表格数据",0,"110253"],
            ["./MaAsLin分析结果图.pdf", "txt", "MaAsLin分析结果图", 0, ""]
        ])
        super(MaaslinWorkflow, self).end()
