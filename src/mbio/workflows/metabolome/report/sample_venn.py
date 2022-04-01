# -*- coding: utf-8 -*-


from biocluster.workflow import Workflow
import os
import json
import pandas as pd


class SampleVennWorkflow(Workflow):
    """
    代谢集Venn分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SampleVennWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metab_table", "type": "infile", "format": "metabolome.express"},
            {"name": "group_detail", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "threshold", "type": "float", "default":50},
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
            {"name": "name", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.venn = self.add_module("metabolome.sample_venn")
        self.output_dir = self.venn.output_dir


    def run(self):
        options = {
            "metab_table": self.option('metab_table'),
            "group_detail" : self.option("group_detail"),
            "threshold" : self.option("threshold")
        }
        self.venn.set_options(options)
        self.venn.on('end', self.set_db)
        self.venn.run()
        super(SampleVennWorkflow, self).run()

    def set_db(self):
        """
        导表程序
        """
        venn_api = self.api.api('metabolome.sample_venn')
        venn_table_result = self.venn.output_dir + "/venn_table.xls"

        group_num = len(json.loads(self.option('group_detail')).keys())
        venn_api.add_venn_detail(venn_table_result, self.option('main_table_id'),group_num, self.option("metab_desc").path)
        ###venn_api.add_venn_graph(venn_graph_result, self.option('main_table_id'))
        # self.option("save_pdf", 0)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "expvenn",
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
            [".", "", "代谢集venn结果目录", 0],
            ["venn_table.xls", "xls", "venn分析结果表",0],
           # ["metabset_detail.xls", "xls", "各代谢集包含的代谢物列表"]
        ])
        super(SampleVennWorkflow, self).end()

if __name__ == '__main__':
    from biocluster.wsheet import Sheet

    data = {
        'name': 'test_sample_venn',
        'id': 'tsg_36964',
        'type': 'workflow',
        'options': {
            "metab_table" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Preprocess/output/pos/metab_abund.txt",
            "group_detail" : '{"a":["NG_D1_A","NG_D1_B"],"b":["NG_D1_C","NG_D1_D"]}',
            "main_table_id" : "5e4ce5e817b2bf4b326f4b20"
        }
    }

    wsheet = Sheet(data=data)

    wf = SampleVennWorkflow(wsheet)
    wf.run()