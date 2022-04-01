# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import os
import time
import re
from mbio.packages.labelfree.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import json
import glob

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
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("labelfree.wgcna.wgcna_module")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/6_Wgcna')
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
        self.run_tool()
        super(WgcnaModuleWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        dump_tool = self.api.api("labelfree.wgcna")
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

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir+'/'
        chart.output_dir = self.work_dir+'/'
        module_stat = os.path.join(self.tool.work_dir, "module_size.stat.xls")
        module_corr = os.path.join(self.tool.work_dir, "module_corr.matrix.xls")
        module_corr_tree = os.path.join(self.tool.work_dir, "module_corr.tree.txt")
        module_tree = os.path.join(self.tool.work_dir, "unmerged.module_corr.dendrogram.txt")
        if os.path.exists(module_stat):
            chart.chart_wgcna_module_column(module_stat)
        if os.path.exists(module_corr) and os.path.exists(module_corr_tree):
            chart.chart_wgcna_module_corr(module_corr,module_corr_tree)
        if os.path.exists(module_tree):
            chart.chart_wgcna_module_tree(module_tree)
        chart.to_pdf()

        # move pdf to result dir
        for ori_filename,filename in [\
        # ["?????","wgcna_module","dendrogram.pdf"],\
        ["wgcna.module_stat.column.pdf","num.pdf"],#chart_wgcna_module_column\
        ["wgcna.mofule_corr.heat_corr.pdf","relate.pdf"],#chart_wgcna_module_corr\
        ["wgcna.module_tree.heat_corr.pdf","modulecluster.pdf"],#chart_wgcna_module_tree\
        ]:
            if os.path.exists(os.path.join(self.work_dir, ori_filename)):
                os.link(os.path.join(self.work_dir, ori_filename), os.path.join(self.tool.output_dir, filename))

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.tool.output_dir)
        self.inter_dirs = [
            ["6_Wgcna", "", "WGCNA",0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "wgcna_module", 0, "211246"],
            ["*dendrogram.pdf", "", "模块分类树", 0],
            ["num.pdf", "", "模块成员统计图", 0],
            ["relate.pdf", "", "模块相关性图", 0],
            ["modulecluster.pdf", "", "模块聚类图", 0],
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
