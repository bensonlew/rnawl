# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import os
import time
import re
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import json


class WgcnaModuleWorkflow(Workflow):
    """
    666
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(WgcnaModuleWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="main_id", type='string'),
            dict(name="mergeCutHeight", type='string'),
            dict(name="power", type="string"),
            dict(name="networkType", type="string"),
            dict(name="minModuleSize", type="string"),
            dict(name="minKMEtoStay", type="string"),
            dict(name="exp_level", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("lnc_rna.wgcna.wgcna_module")
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/07 Advanced_Analysis/02 WGCNA')
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
        super(WgcnaModuleWorkflow, self).send_log(data)

    def run(self):
        self.tool.on("end", self.set_db)
        self.get_run_log()
        self.run_tool()
        super(WgcnaModuleWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("lnc_rna", table="sg_wgcna_module", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        dump_tool = self.api.api("lnc_rna.wgcna")
        # save workflow output path
        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:',self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            self.workflow_output = self.workflow_output_tmp.replace('sanger:','/mnt/ilustre/data/')
        output_dir = self.workflow_output
        dump_tool.update_db_record('sg_wgcna_module', self.option('main_id'), output_dir=output_dir)
        # add result info
        exp_matrix, gene_id2gene_name = self.option('exp_matrix').split(";")
        dump_tool.add_module_detail(
            self.tool.work_dir,
            gene_id2gene_name,
            exp_matrix,
            self.option('main_id')
        )
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.tool.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.tool.output_dir)
        self.inter_dirs = [
            ["07 Advanced_Analysis", "", "高级分析结果目录", 0],
            ["07 Advanced_Analysis/02 WGCNA", "", "WGCNA分析", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "模块识别文件", 0],
            ['./module_size.stat.xls', 'xls', '模块成员统计表', 0],
            ['./seq_id2gene_name.xls', 'xls', '模块成员基因列表', 0],
            ['./eigengenes.xls', 'xls', '各模块特征基因信息表', 0],
            ['./module_corr.matrix.xls', 'xls', '模块相关性关系表',0],
            ['./membership.xls', 'xls', '模块成员聚类详情表', 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0]
        ])
        result_dir.add_regexp_rules([
            [r'block_.*_dendrogram\..*pdf', 'pdf', '模块分类树图pdf', 0],
            [r'block_.*_dendrogram\..*png', 'png', '模块分类树图png', 0],
            [r'.*\.RData', '', '模块分析R数据, 可用R打开查看', 0],
        ])
        super(WgcnaModuleWorkflow, self).end()

    def run_tool(self):
        options = dict(
            datExpr=self.option('exp_matrix').split(";")[0],
            mergeCutHeight=float(self.option('mergeCutHeight')),
            power=int(self.option('power')),
            networkType=self.option('networkType'),
            minModuleSize=int(self.option('minModuleSize')),
            minKMEtoStay=float(self.option('minKMEtoStay')),
        )
        self.tool.set_options(options)
        self.tool.run()
