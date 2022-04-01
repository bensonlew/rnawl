# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import os
import time
import re
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.project_demo.interaction_rerun.interaction_delete import InteractionDelete,linkfile,linkdir
from mbio.packages.ref_rna_v2.chart_advance import ChartAdvance
import glob

class WgcnaRelateWorkflow(Workflow):
    """
    666
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(WgcnaRelateWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_eigengenes', type='string'),
            dict(name="main_id", type='string'),
            dict(name="trait_path", type='string'),
            dict(name="trait_type", type="string"),
            dict(name="corr_method", type="string"),
            dict(name="exp_level", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
            {"name": "block_Rdata", "type": "infile", "format": "ref_rna_v2.common"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("ref_rna_v2.wgcna.wgcna_relate")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/06 Advanced_Analysis/01 WGCNA')
        self.inter_dirs = []
        if self._sheet.rerun:
            self.logger.info("该交互为重运行项目,先删除mongo库内容")
            interactiondelete = InteractionDelete(bind_object=self, project_type="ref_rna_v2",
                                                  main_id=self.option('main_id'))
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
        super(WgcnaRelateWorkflow, self).send_log(data)

    def run(self):
        self.tool.on("end", self.set_db)
        self.get_run_log()
        self.run_tool()
        super(WgcnaRelateWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("ref_rna_v2", table="sg_wgcna_relate", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()



    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        dump_tool = self.api.api("ref_rna_v2.wgcna")
        # save workflow output path
        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:',self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            self.workflow_output = self.workflow_output_tmp.replace('sanger:','/mnt/ilustre/data/')
        output_dir = self.workflow_output
        dump_tool.update_db_record('sg_wgcna_relate', self.option('main_id'), output_dir=output_dir)
        # add result info
        seq_annot = self.option('exp_eigengenes').split(";")[2]
        dump_tool.add_relate_detail(
            self.tool.work_dir,
            seq_annot,
            self.option('main_id')
        )
        self.end()

    def chart(self):
        chart = ChartAdvance()
        chart.work_dir = self.work_dir + "/"

        relation_corr = self.tool.work_dir + "/module_trait.correlation.xls"
        relation_corr_pvalue = self.tool.work_dir + "/module_trait.correlation_pvalues.xls"
        module_stat = self.tool.work_dir + "/module_size.stat.xls"

        chart.chart_wgcna_relation_corr(relation_corr, relation_corr_pvalue, module_stat)
        gene_trait_corr = self.tool.work_dir + "/gene_trait.correlation.xls"
        seq_annot = self.option('exp_eigengenes').split(";")[2]
        chart.chart_wgcna_relation_ms(gene_trait_corr, seq_annot, relation_corr, relation_corr_pvalue)

        chart.to_pdf()

        pdf_file = glob.glob(self.work_dir + "/wgcna*relation*.pdf")
        for p in pdf_file:
            linkfile(p, self.tool.output_dir + os.path.basename(p))


    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.tool.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.tool.output_dir)
        self.inter_dirs = [
            ["06 Advanced_Analysis", "", "高级分析结果目录",0],
            ["06 Advanced_Analysis/01 WGCNA", "", "WGCNA分析", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "模块分析文件",0],
            ["./block_*_gene_*pdf", "", "基因与表型相关性热图pdf",0,"211623"],
            ["./block_*_gene_*png", "", "基因与表型相关性热图png",0,"211624"],
            ["./module_trait.correlation.xls", "", "模块与表型相关性系数表",0,"211625"],
            ["./module_trait.correlation_pvalues.xls", "", "模块与表型相关显著性统计表",0,"211626"],
            ["./gene_trait.correlation.xls", "", "基因与表型相关性系数表",0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(WgcnaRelateWorkflow, self).end()

    def run_tool(self):
        options = dict(
            datExpr=self.option('exp_eigengenes').split(";")[0],
            MEs=self.option('exp_eigengenes').split(";")[1],
            traits=self.option('trait_path'),
            corType=self.option('corr_method')
        )
        if self.option('block_Rdata').is_set:
            options.update({
                "block_Rdata": self.option('block_Rdata').prop['path']
            })
        self.tool.set_options(options)
        self.tool.run()
