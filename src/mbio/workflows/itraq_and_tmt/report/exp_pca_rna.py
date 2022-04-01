# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time


class ExpPcaWorkflow(Workflow):
    """
    表达量分布
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpPcaWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name="pca_main_id", type='string'),
            dict(name="type", type='string'),
            dict(name="group",type="infile",format="itraq_and_tmt.common"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("itraq_and_tmt.exp_pca")
        group_dict = json.loads(self.option("group_dict"))
        self.ellipse = None
        i=0
        for sample_infos in group_dict.values():
           if len(sample_infos) < 3:
               i += 1
           else:
               continue
        if i == 0:
            self.ellipse = self.add_tool("graph.ellipse")

    def run(self):
        if not self.ellipse is None:
            self.tool.on("end", self.run_ellipse)
            self.ellipse.on("end", self.set_db)
            self.run_tool()
        else:
            self.tool.on("end", self.set_db)
            self.run_tool()
        super(ExpPcaWorkflow, self).run()

    def run_ellipse(self):
        options = {}
        if self.option("group"):
            options['group_table'] = self.option("group").prop['path']
        options['analysis'] = 'pca'
        options['pc_table'] = self.tool.output_dir + "/PCA.xls"
        self.ellipse.set_options(options)
        self.ellipse.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("itraq_and_tmt.all_exp")
        # add result info
        all_exp.add_exp_pca2(self.tool.work_dir, main_id=self.option('pca_main_id'), )
        if not self.ellipse is None:
            all_exp.insert_ellipse_table(self.ellipse.work_dir + '/ellipse_out.xls', str(self.option('pca_main_id')))
        else:
            pass
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "PCA分析结果目录", 0, "211068"],
        ])
        super(ExpPcaWorkflow, self).end()

    def run_tool(self):
        options = dict(
            exp=self.option('exp_matrix'),
        )
        self.tool.set_options(options)
        self.tool.run()
