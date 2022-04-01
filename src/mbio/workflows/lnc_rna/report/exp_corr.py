# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import os
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
import re
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import shutil


class ExpCorrWorkflow(Workflow):
    """
    表达量相关性
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpCorrWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name='scm', type='string', default='complete'),
            dict(name='scd', type='string', default='correlation'),
            dict(name='corr_method', type='string', default='pearson'),
            dict(name="corr_main_id", type='string'),
            dict(name="log_base", type='int', default=None),
            dict(name="type", type='string'),
            dict(name="rna_type", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("lnc_rna.exp_corr")
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/01 Express/02 Exp_Corr')
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
        super(ExpCorrWorkflow, self).send_log(data)

    def run(self):
        self.tool.on("end", self.set_db)
        self.get_run_log()
        self.run_tool()
        super(ExpCorrWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("lnc_rna", table="sg_exp_corr", main_id=self.option('corr_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("lnc_rna.all_exp")
        # add result info
        all_exp.add_exp_corr2(self.tool.work_dir, main_id=self.option('corr_main_id'), )
        self.end()

    def end(self):
        if os.path.exists(self.tool.output_dir + '/upload'):
            shutil.rmtree(self.tool.output_dir + '/upload')
        os.mkdir(self.tool.output_dir + '/upload')
        shutil.copyfile(self.tool.output_dir + '/sample_correlation.xls', self.tool.output_dir + '/upload/sample_correlation.xls')
        if os.path.exists(os.path.join(self.tool.output_dir + '/upload', os.path.basename(self.run_log))):
            os.remove(os.path.join(self.tool.output_dir + '/upload', os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.tool.output_dir + '/upload', os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.tool.output_dir + '/upload')
        self.inter_dirs = [
            ["01 Express", "", "表达量分析结果目录", 0],
            ["01 Express/02 Exp_Corr", "", "样本间相关性分析", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "样本间相关性分析文件", 0, "211064"],
            ["./sample_correlation.xls", "xls", "样本间相关性系数表", 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0]
        ])
        super(ExpCorrWorkflow, self).end()

    def run_tool(self):
        options = dict(
            exp=self.option('exp_matrix'),
            scm=self.option('scm'),
            scd=self.option('scd'),
            corr_method=self.option('corr_method'),
        )
        if self.option('log_base'):
            options.update(log_base=self.option('log_base'))
        self.tool.set_options(options)
        self.tool.run()
