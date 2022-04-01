# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import re
import json
import os
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


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
            dict(name="step3output", type="infile", format='ref_rna_v2.common_dir'),
            dict(name="step2output", type="infile", format='ref_rna_v2.common_dir'),
            # to update status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("ref_rna_v2.wgcna.wgcna_network")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/07 Advanced_Analysis/02 WGCNA')
        self.inter_dirs = []

    def run(self):
        self.tool.on("end", self.set_db)
        self.get_run_log()
        self.run_tool()
        super(WgcnaNetworkWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("whole_transcriptome", table="wgcna_network", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(WgcnaNetworkWorkflow, self).send_log(data)


    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        dump_tool = self.api.api("whole_transcriptome.wgcna")
        # save workflow output path--建议标配
        workflow_output = self._sheet.output
        if re.match(r'tsanger:', workflow_output):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        dump_tool.add_network_detail(self.tool.output_dir, self.option("main_id"))
        # network = workflow_output+'/'+self.option("module").strip().replace(",", "_")+'.network.json'
        network = workflow_output+'/' + 'network.json'
        dump_tool.update_db_record('wgcna_network', self.option('main_id'), output_dir=network, status="end")
        # add result info
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.tool.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.tool.output_dir)
        self.inter_dirs = [
            ["07 Advanced_Analysis", "", "高级分析结果目录",0],
            ["07 Advanced_Analysis/02 WGCNA", "", "WGCNA分析", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "WGCNA可视化分析文件", 0],
            ["./*.json", "", "该模块网络图数据, 文本格式", 0],
            ["./*network.nodes.txt", "", "该模块网络图节点列表", 0],
            ["./*network.edges.txt", "", "该模块网络图边列表", 0],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(WgcnaNetworkWorkflow, self).end()

    def run_tool(self):
        options = dict(
            module=self.option('module'),
            threshold=self.option('threshold'),
            top=self.option('top'),
            step3output=self.option('step3output').prop['path'],
            step2output=self.option('step2output').prop['path'],
        )
        self.tool.set_options(options)
        self.tool.run()
