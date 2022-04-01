# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

from biocluster.workflow import Workflow
import os
import pandas as pd

class MetabsetKeggcWorkflow(Workflow):
    """
    化合物注释
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MetabsetKeggcWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "anno_overview", "type": "infile", "format": "metabolome.overview_table"},
            {"name": "metabset", "type": "infile", "format": "metabolome.metabset"},
            {"name": "update_info", "type": "string"},
            {"name": "name", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "database","type":"string","default":"CBR"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        #self.keggc = self.add_tool("metabolome.metabset.keggc")
        self.keggc = self.add_tool("metabolome.annotation.anno_keggc")
        #self.output_dir = self.keggc.output_dir

    ## 根据metabset 对总览表进行筛选
    def select_overview(self):
        overview = pd.read_table(self.option('anno_overview').path,sep='\t')
        metabs = pd.read_table(self.option('metabset').path,sep='\t',header=-1)
        metabs_list = metabs[0].tolist()
        new_names = {}
        if 'metab' in overview.columns:
            new_names['metab'] = 'Metabolite'
        if 'compound_id' in overview.columns:
            new_names['compound_id'] = 'KEGG Compound ID'
        overview.rename(columns=new_names,inplace=True)
        sub_overview = overview[overview['metab_id'].apply(lambda x: x in metabs_list)]
        self.new_overview = self.work_dir+'/sub_overview.xls'
        sub_overview.to_csv(self.new_overview,sep='\t',index=False)


    def run(self):
        # options = {
        #     "anno_overview": self.option('anno_overview'),
        #     "metabset": self.option("metabset")
        # }
        self.select_overview()
        options = {
            'metab_table': self.new_overview,
            'database_name':self.option('database'),
            "task_id": "_".join(self._sheet.id.split("_")[0:2])
        }
        self.keggc.set_options(options)
        self.keggc.on('end', self.set_db)
        self.keggc.run()
        super(MetabsetKeggcWorkflow, self).run()

    def set_db(self):
        """
        导表程序
        """
        anno_api = self.api.api('metabolome.metabset_keggc')
        anno_api.add_metabsetc_detail(self.option('main_table_id'), self.keggc.option('stat_out').path)
        # self.option("save_pdf", 0)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "metabsetkeggc",
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        if self.option("save_pdf"):
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        out_file = self.output_dir + '/level.xls'
        if os.path.exists(out_file):
            os.remove(out_file)
        os.link(self.keggc.option('stat_out').path,out_file)  #stat.xls 改名成了level.xls ,和老项目统一
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "代谢集KEGG化合物分类结果", 0, "150035"],
            ["level.xls", "xls", "代谢集KEGG化合物分类层级表", 0, "150036"]
        ])
        super(MetabsetKeggcWorkflow, self).end()
