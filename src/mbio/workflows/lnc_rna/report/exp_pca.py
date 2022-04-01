# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import os
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
import re
from biocluster.core.function import filter_error_info, link, CJsonEncoder


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
            dict(name="rna_type", type="string"),
            dict(name="group_table",type="string"),
            dict(name="analysis_type",type="string",default="pca"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("lnc_rna.exp_pca")
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
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/01 Express/03 Exp_PCA')
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
        super(ExpPcaWorkflow, self).send_log(data)

    def run(self):
        self.get_run_log()
        if not self.ellipse is None:
            self.tool.on("end", self.run_ellipse)
            self.ellipse.on("end", self.set_db)
            self.run_tool()
        else:
            self.tool.on("end", self.set_db)
            self.run_tool()
        super(ExpPcaWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("lnc_rna", table="sg_exp_pca", main_id=self.option('pca_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_ellipse(self):
        options = {}
        if self.option("group_table"):
            options['group_table'] = self.option("group_table")
        pc_map = {'pca':"/PCA.xls",'pcoa': "/Pcoa/pcoa_sites.xls",
                  'dbrda':'/Dbrda/db_rda_sites.xls','nmds':'/Nmds/nmds_sites.xls',
                  'rda_cca': '/Rda'
                  ##'rda_cca'
                  }
        options['analysis'] = self.option('analysis_type')
        options['pc_table'] = self.tool.output_dir + pc_map[self.option('analysis_type')]
        self.ellipse.set_options(options)
        self.ellipse.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("lnc_rna.all_exp")
        # add result info
        all_exp.add_exp_pca2(self.tool.work_dir, main_id=self.option('pca_main_id'), )
        if not self.ellipse is None:
            all_exp.insert_ellipse_table(self.ellipse.work_dir + '/ellipse_out.xls', str(self.option('pca_main_id')))
        else:
            pass
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.tool.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.tool.output_dir)
        self.inter_dirs = [
            ["01 Express", "", "表达量分析结果目录", 0],
            ["01 Express/03 Exp_PCA", "", "样本间PCA分析", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "样本间PCA分析文件", 0,],
            ["./Explained_variance_ratio.xls", 'xls', '主成分解释表',0],
            ['./PCA.xls', 'xls', '样本间PCA详情表',0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0]
        ])
        super(ExpPcaWorkflow, self).end()

    def run_tool(self):
        options = dict(
            exp=self.option('exp_matrix'),
        )
        self.tool.set_options(options)
        self.tool.run()
