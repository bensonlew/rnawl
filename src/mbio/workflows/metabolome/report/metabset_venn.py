# -*- coding: utf-8 -*-
# __author__ = 'liulinmeng'

from biocluster.workflow import Workflow
import os


class MetabsetVennWorkflow(Workflow):
    """
    代谢集Venn分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MetabsetVennWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metabset", "type": "infile", "format": "metabolome.mul_metabset"},
            {"name": "update_info", "type": "string"},
            {"name": "name", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
            {"name": "name", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.venn = self.add_tool("metabolome.metabset.venn")
        self.output_dir = self.venn.output_dir

    def run(self):
        options = {
            "list_file": self.option("metabset")
        }
        self.venn.set_options(options)
        self.venn.on('end', self.set_db)
        self.venn.run()
        super(MetabsetVennWorkflow, self).run()

    def set_db(self):
        """
        导表程序
        """
        venn_api = self.api.api('metabolome.metabset_venn')
        venn_table_result = self.venn.output_dir + "/venn_table.xls"
        venn_graph_result = self.venn.output_dir + "/metabset_detail.xls"
        venn_api.add_venn_detail(venn_table_result, self.option('main_table_id'))
        venn_api.add_venn_graph(venn_graph_result, self.option('main_table_id'))
        # self.option("save_pdf", 0)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "metabsetvenn",
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
            [".", "", "代谢集venn结果目录", 0, "150024"],
            ["venn_table.xls", "xls", "venn分析结果表",0, "150025"],
            ["metabset_detail.xls", "xls", "各代谢集包含的代谢物列表",0,"150026"]
        ])
        super(MetabsetVennWorkflow, self).end()
