# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

from biocluster.workflow import Workflow
import os


class MetabsetKeggpWorkflow(Workflow):
    """
    化合物注释
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MetabsetKeggpWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "anno_overview", "type": "infile", "format": "sequence.profile_table"},
            {"name": "metabset", "type": "infile", "format": "metabolome.mul_metabset"},
            {"name": "trans_file", "type": "infile", "format": "sequence.profile_table"}, #  转化id使用
            # {"name": "set_name", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "name", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.keggp = self.add_tool("metabolome.metabset.keggp")
        self.output_dir = self.keggp.output_dir

    def run(self):
        options = {
            "anno_overview": self.option('anno_overview'),
            "metabset": self.option("metabset"),
            "trans_file": self.option("trans_file"),
            "task_id": "_".join(self._sheet.id.split("_")[0:2])
        }
        self.keggp.set_options(options)
        self.keggp.on('end', self.set_db)
        self.keggp.run()
        super(MetabsetKeggpWorkflow, self).run()

    def set_db(self):
        """
        导表程序
        """
        anno_api = self.api.api('metabolome.metabset_keggp')
        name_list = self.option('metabset').prop['set_name']
        bset_name = None
        if len(name_list) == 2:
            bset_name = name_list[1]
        anno_api.add_metabsetp_level(self.option('main_table_id'), self.keggp.option('level_out').path,bset_name=bset_name)
        stat_path = self.keggp.output_dir + '/stat.xls'

        anno_api.update_metabsetp(self.option('main_table_id'), name_list)
        if len(name_list) == 1:
            anno_api.add_metabsetp_stat(self.option('main_table_id'), name_list[0], stat_path)
        else:
            anno_api.add_metabsetp_stat(self.option("main_table_id"), name_list[0], self.keggp.output_dir + '/stat1.xls')
            anno_api.add_metabsetp_stat(self.option("main_table_id"), name_list[1], self.keggp.output_dir + '/stat2.xls')
        # self.option("save_pdf", 0)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "metabsetkeggp",
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        if self.option("save_pdf"):
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "代谢集KEGG功能通路结果", 0, "150037"],
            ["level.xls", "xls", "代谢集KEGG通路分类层级表", 0, "150038"],
            ["stat.xls", "xls", "代谢集KEGG通路统计表", 0, "150039"],
            ["pathway_img", "", "代谢集KEGG通路图", 0, "150041"]
        ])
        super(MetabsetKeggpWorkflow, self).end()
