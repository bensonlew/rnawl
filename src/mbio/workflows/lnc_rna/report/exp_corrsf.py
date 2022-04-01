# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import os
import re
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import json


class ExpCorrsfWorkflow(Workflow):
    """
    表达量相关性
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpCorrsfWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'exp_matrix', 'type': 'infile', 'format': 'lnc_rna.express_matrix'},
            {'name': 'genes_info', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'gt', 'type': 'string', 'default': None},
            {'name': 'anno', 'type': "infile", 'format': "lnc_rna.common"},
            {'name': 'output', 'type': 'string', 'default': None},
            {'name': 'pvalue_cutoff', 'type': 'float', 'default': 0.05},
            {'name': 'qvalue_cutoff', 'type': 'float', 'default': 0.05},
            {'name': 'cor_cutoff', 'type': 'float', 'default': 0.8},
            {'name': 'corr_way', 'type': 'string', 'default': "spearmanr"},
            {'name': 'padjust_way', 'type': 'string', 'default': "fdr_bh"},
            {'name': 'sig_type', 'type': 'int', 'default': 1},
            {'name': 'group_dict', 'type': 'string', 'default': None},
            {'name': 'group_id', 'type': 'string', 'default': None},
            {'name': 'corr_main_id', 'type': 'string', 'default': None},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("lnc_rna.geneset.exp_corrsf")
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/04 GeneSet/10 Exp_Network')
        self.inter_dirs = []

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
        super(ExpCorrsfWorkflow, self).send_log(data)

    def run(self):
        self.tool.on("end", self.set_db)
        self.get_run_log()
        self.run_tool()
        super(ExpCorrsfWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("lnc_rna", table="sg_exp_corrsf", main_id=self.option('corr_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        exp_corrsf = self.api.api("lnc_rna.exp_corrsf")
        # add result info
        # self, upload_dir, work_dir, genes_info = None, main_id = None, corr_way = None,
        # padjust_way = None, params = None, project_sn = None, task_id = None
        genes_info = None
        if self.option('genes_info').is_set:
            genes_info = self.option('genes_info').path
        exp_corrsf.add_ExpCorrsf(self.get_workflow_output_dir(), self.tool.work_dir,
                                 genes_info=genes_info, main_id=self.option('corr_main_id'))
        self.end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if re.match(r'tsanger:', workflow_output):
            workflow_output = workflow_output.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:', '/mnt/ilustre/data/')
        return workflow_output

    def end(self):
        if os.path.exists(os.path.join(self.tool.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.tool.output_dir)
        self.inter_dirs = [
            ["04 GeneSet", "", "基因集分析结果目录", 0],
            ["04 GeneSet/10 Exp_Network", "", "表达相关性分析", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "表达相关性分析文件", 0],
            ["./express_correlation_info.xls", "xls", "表达相关性分析表", 0, "211066"],
            ["./record.json", "json", "表达量相关性分析绘制网络图表", 0, "211067"],
            ['run_parameter.txt', 'txt', '运行参数日志', 0]
        ])
        super(ExpCorrsfWorkflow, self).end()

    def run_tool(self):
        options = dict(
            exp=self.option('exp_matrix'),
            gt=self.option('gt'),
            anno=self.option('anno').prop['path'],
            pvalue_cutoff=self.option('pvalue_cutoff'),
            output=self.option('output'),
            qvalue_cutoff=self.option('qvalue_cutoff'),
            cor_cutoff=self.option('cor_cutoff'),
            corr_way=self.option('corr_way'),
            padjust_way=self.option('padjust_way'),
            sig_type=self.option('sig_type'),
        )
        self.tool.set_options(options)
        self.tool.run()
