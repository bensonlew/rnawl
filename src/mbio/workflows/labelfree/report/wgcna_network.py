# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import re
import glob
import os


class WgcnaNetworkWorkflow(Workflow):
    """
    666
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(WgcnaNetworkWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='module', type='string'),
            dict(name="main_id", type='string'),
            dict(name="threshold", type='string'),
            dict(name="top", type='string'),
            dict(name="step3output", type="infile", format='labelfree.common_dir'),
            dict(name="step2output", type="infile", format='labelfree.common_dir'),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("labelfree.wgcna.wgcna_network")

    def run(self):
        self.tool.on("end", self.set_db)
        self.run_tool()
        super(WgcnaNetworkWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        dump_tool = self.api.api("labelfree.wgcna")
        # save workflow output path--建议标配
        workflow_output = self._sheet.output
        if re.match(r'tsanger:', workflow_output):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        dump_tool.add_network_detail(self.tool.output_dir, self.option("main_id"))
        network = workflow_output+'/'+self.option("module").strip().replace(",", "_")+'.network.json'
        dump_tool.update_db_record('sg_wgcna_network', self.option('main_id'), output_dir=network, status="end")
        # add result info
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "wgcna_network", 0,  "211247"],
        ])
        super(WgcnaNetworkWorkflow, self).end()

    def run_tool(self):
        # 为了调用人家的package也是没办法
        step3output = self.option('step3output').prop['path']
        for file in glob.glob(step3output+'/*'):
            if 'protein' in file:
                file_ = file.replace('protein', 'gene')
                os.rename(file, file_)
        step2output = self.option('step2output').prop['path']
        for file in glob.glob(step2output + '/*'):
            if 'protein' in file:
                file_ = file.replace('protein', 'gene')
                os.rename(file, file_)
        options = dict(
            module=self.option('module'),
            threshold=self.option('threshold'),
            top=self.option('top'),
            step3output=step3output,
            step2output=step2output,
        )
        self.tool.set_options(options)
        self.tool.run()
