# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
import re
import glob
import shutil
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
import os
from mbio.packages.ref_rna_v2.chart import Chart
from mbio.packages.project_demo.interaction_rerun.interaction_delete import InteractionDelete,linkfile,linkdir
from mbio.packages.ref_rna_v2.chart_geneset import ChartGeneset
import glob
import shutil

class ExpCorrsfWorkflow(Workflow):
    """
    表达量相关性
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpCorrsfWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'exp_matrix', 'type': 'infile', 'format': 'denovo_rna_v2.express_matrix'},
            {'name': 'gt', 'type': 'string', 'default': None},
            {'name': 'anno', 'type': "infile", 'format': "ref_rna_v2.common"},
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
        self.tool = self.add_tool("ref_rna_v2.exp_corrsf")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/03 GeneSet/09 Exp_Network')
        self.inter_dirs = []
        if self._sheet.rerun:
            self.logger.info("该交互为重运行项目,先删除mongo库内容")
            interactiondelete = InteractionDelete(bind_object=self, project_type="ref_rna_v2",
                                                  main_id=self.option('corr_main_id'))
            interactiondelete.delete_interactions_records()

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

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        exp_corrsf = self.api.api("ref_rna_v2.exp_corrsf")
        # add result info
        exp_corrsf.add_ExpCorrsf(self.get_workflow_output_dir(), self.tool.work_dir, main_id=self.option('corr_main_id'), )
        self.end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if re.match(r'tsanger:', workflow_output):
            workflow_output = workflow_output.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:', '/mnt/ilustre/data/')
        return workflow_output

    def get_run_log(self):
        get_run_log = GetRunLog("ref_rna_v2", table="sg_exp_corrsf", main_id=self.option('corr_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def chart(self):
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"
        net_json = self.tool.output_dir + "/record.json"
        chart.chart_geneset_corr_net(net_json)
        chart.to_pdf()

        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            print "copy", p, self.tool.output_dir + "/" + os.path.basename(p)
            shutil.copyfile(p, self.tool.output_dir + "/" + os.path.basename(p))


    def end(self):
        # self.chart()
        if os.path.exists(os.path.join(self.tool.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.tool.output_dir)
        self.inter_dirs = [
            ["03 GeneSet", "", "基因集分析结果目录",0],
            ["03 GeneSet/09 Exp_Network", "", "表达相关性分析", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "表达相关性分析文件",0,"211522"],
            ["./express_correlation_info.xls", "xls", "表达相关性分析表", 0, "211066"],
            ["./*.pdf", 'pdf', "相关性网络图",0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
            # ["./record.json", "json", "表达量相关性分析绘制网络图表", 0, "211067"]
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
