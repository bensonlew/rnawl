# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

from biocluster.workflow import Workflow
from collections import OrderedDict
from bson.objectid import ObjectId
import os
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
import re
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import glob
from mbio.packages.denovo_rna_v2.chart import Chart


class SsrWorkflow(Workflow):
    """
    基因集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SsrWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "unigene_fa", "type": "infile", "format":"denovo_rna_v2.trinity_fasta"},  # 这是fa文件的格式检查
            {"name": "bed", "type": "infile", "format": "gene_structure.bed"},
            {"name": "rept_1", "type": "int", "default": 10},
            {"name": "rept_2", "type": "int", "default": 6},
            {"name": "rept_3", "type": "int", "default": 5},
            {"name": "rept_4", "type": "int", "default": 5},
            {"name": "rept_5", "type": "int", "default": 5},
            {"name": "rept_6", "type": "int", "default": 5},
            {"name": "ssr_distance", "type": "int", "default": 100},
            {"name": "update_info", "type": "string"},
            {"name": "ssr_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.ssr = self.add_tool("denovo_rna_v2.ssr")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/05 Structure_Analysis/03 SSR')
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
        super(SsrWorkflow, self).send_log(data)


    def run(self):
        # super(SsrWorkflow, self).run() 不能把super写在前面，这样框架就不能开启监听，必须先设置绑定关系，最好是写在run函数的最下面
        self.ssr.on("end", self.set_db)
        self.get_run_log()
        self.run_ssr()
        super(SsrWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("denovo_rna_v2", table="sg_ssr", main_id=self.option('ssr_id'), dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_ssr = self.api.api("denovo_rna_v2.ssr")
        self.logger.info("开始进行ssr的导表")
        ssr_statistic_path = self.work_dir + '/Ssr/' + 'ssr_type.txt'
        ssr_detail_path = self.work_dir + '/Ssr/' + 'tmp.txt'
        ssr_class_path = self.work_dir + '/Ssr/' + 'ssr_repeats_class.txt'
        ssr_statistic_dict = OrderedDict()
        ssr_statistic_dict_p = OrderedDict()
        with open(ssr_statistic_path, 'r') as f2:
            header2 = f2.readline()
            for line in f2:
                line = line.strip().split('\t')
                ssr_statistic_dict[line[0]] = int(line[1])

        with open(ssr_statistic_path, 'r') as f2:
            for i in range(4):
                read_over = f2.readline()
            for line in f2:
                line = line.strip().split('\t')
                ssr_statistic_dict_p[line[0] + '_p'] = float("%.6f" % float(line[2]))
        api_ssr.update_db_record('sg_ssr', self.option("ssr_id"), ssr_statistic_dict=ssr_statistic_dict, ssr_statistic_dict_p=ssr_statistic_dict_p, main_id=ObjectId(self.option("ssr_id")))
        api_ssr.add_ssr_detail(self.option("ssr_id"), ssr_detail_path)
        api_ssr.add_ssr_class(self.option("ssr_id"), ssr_class_path)
        self.end()

    def chart(self):
        chart = Chart()
        ssr_statistic_path = self.work_dir + '/Ssr/' + 'ssr_type.txt'
        chart.denovo_chart_ssr_stat(ssr_statistic_path)
        chart.to_pdf()
        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            os.link(p, self.output_dir + "/" + os.path.basename(p))


    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.ssr.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.ssr.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.ssr.output_dir, os.path.basename(self.run_log)))
        os.system('cp {} {}'.format(os.path.join(self.ssr.output_dir,'ssr*'),self.output_dir))
        os.system('cp {} {}'.format(os.path.join(self.ssr.output_dir,'run_parameter.txt'),self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["05 Structure_Analysis", "", "SNP分析结果目录",0],
            ["05 Structure_Analysis/03 SSR", "", "SSR分析", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "SSR分析文件",0,],
            ["ssr_analysis_details", " ", "SSR分析详情表",0],
            ["ssr_repeats", " ", "SSR类型分布统计表",0],
            ["ssr_type", " ", "SSR类型统计表",0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ['SSR.stat.pie.pdf', 'pdf', 'SSR统计图', 0]
        ])
        super(SsrWorkflow, self).end()

    def run_ssr(self):
        # datajson里面的参数是接口传给workflow的，然后这里的opts是workflow传给自己要调用的tool的里面的，二者是独立的，所以如果这里写漏了参数，那么tool的时候就会报错，但是不会追踪到这里，
        # 所以一定要仔细核对参数
        opts = {
            "rept_1": self.option("rept_1"),
            "rept_2": self.option("rept_2"),
            "rept_3": self.option("rept_3"),
            "rept_4": self.option("rept_4"),
            "rept_5": self.option("rept_5"),
            "rept_6": self.option("rept_6"),
            "ssr_distance": self.option("ssr_distance"),
            "unigene_fa": self.option("unigene_fa"),
            "bed": self.option("bed"),

        }
        self.ssr.set_options(opts)
        self.ssr.run()
