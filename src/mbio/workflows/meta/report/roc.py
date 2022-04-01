# -*- coding: utf-8 -*-
# __author__ = 'zhangpeng'
# last modified: Liulinmeng 20181212
# last modified: qingchen.zhang@20190312

""""""

import datetime
from biocluster.workflow import Workflow
import re
import os
import json
import shutil
from mbio.packages.meta.common_function import filter_otu_set
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class RocWorkflow(Workflow):
    """
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RocWorkflow, self).__init__(wsheet_object)
        options = [

            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "level", "type": "int", "default": 9},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 分组表，只接受两个分组
            {"name": "method", "type": "string", "default": "sum"},  # 选定的各物种指标在组内求和、均值、中位数后进行roc计算
            # {"name": "problem_type", "type": "int", "default": 2},
            {"name": "confidence_interval", "type": "float", "default": 0.95},
            {"name": "top_n", "type": "int", "default": -1},  # -1表示前端没传次参数
            {"name": "roc_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "group_id", "type": 'string'},
            {"name": "otu_id", "type": 'string'},
            {"name": "group_detail", "type": "string"},
            {"name": "otuset_id", "type": "string"},
            {"name": "env_id", "type": "string"},  # 环境因子表id
            {"name": "env_labs", "type": "string"},  # 选择的环境因子
            {"name": "norm_method", "type": "string", "default": ""},  # 新增的数据标准化方式
        ]
        self.add_option(options)
        # newtable = os.path.join(self.work_dir, 'otutable1.xls')
        # f2 = open(newtable, 'w+')
        # with open(tablepath, 'r') as f:

        self.set_options(self._sheet.options())
        # self.roc = self.add_tool("meta.beta_diversity.roc")
        self.roc = self.add_tool("statistical.roc")
        # self.samples = re.split(',', self.option("samples"))
        self.output_dir = self.roc.output_dir

    def change_otuname(self, tablepath):
        newtable = os.path.join(self.work_dir, 'otutable1.xls')
        f2 = open(newtable, 'w+')
        with open(tablepath, 'r') as f:
            i = 0
            for line in f:
                if i == 0:
                    i = 1
                    f2.write(line)
                else:
                    line = line.strip().split('\t')
                    line_data = line[0].strip().split(' ')
                    line_he = "".join(line_data)
                    line[0] = line_he
                    # line[0] = line_data[-1]
                    for i in range(0, len(line)):
                        if i == len(line) - 1:
                            f2.write("%s\n" % (line[i]))
                        else:
                            f2.write("%s\t" % (line[i]))
        f2.close()
        return newtable

    def run_roc(self):
        if self.option("otuset_id"):

            try:
                filter_otu_set(self.option("otu_table").prop['path'], self.option("otuset_id"),self.option("otu_id"), self.option("level") ,os.path.join(self.work_dir, "fliter_otuset_id.xls"))
            except:
                self.set_error("筛选asv集失败！")
            newtable = self.change_otuname(os.path.join(self.work_dir, "fliter_otuset_id.xls"))
        else:
            newtable = self.change_otuname(self.option('otu_table').prop['path'])
        options = {
            'abu_table': newtable,
            # 'otutable':self.option('otutable'),
            # 'level': self.option('level'),
            # 'envlabs':self.option('envlabs'),
            'group_table': self.option('group_table'),
            'method': self.option('method'),
            # 'problem_type':self.option('problem_type'),
            'top_n': self.option('top_n'),
            'confidence_interval': str(self.option('confidence_interval')),
            'norm_method': self.option('norm_method'),
        }
        if self.option('env_labs'):
            options['env_labs'] = self.option('env_labs')
        self.roc.set_options(options)
        self.roc.on('end', self.set_db)
        self.output_dir = self.roc.output_dir
        self.roc.run()

    def end(self):
        if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.roc.output_dir))
        '''
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "ROC分析结果目录", 0, "110178"],   #modified by hongdongxuan at 91-92 20170324
            ["./roc_curve.xls", "xls", "ROC曲线结果表", 0, "110180"],
            ["./roc_auc.xls", "xls", "AUC计算结果表", 0, "110182"],
            ["./roc_plot_rocarea.xls", "xls", "置信区间", 0, "110179"], # add 2 lines by hongdongxuan 20170324
            ["./roc_table.xls", "xls", "坐标数据", 0, "110183"]
        ])
        super(RocWorkflow, self).end()
        '''
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.roc.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "ROC分析结果目录", 0, "110178"],
            ["./roc_curve.xls", "xls", "ROC曲线结果表", 0, "110180"],
            ["./roc_auc.xls", "xls", "AUC值计算结果表", 0, "110182"],
            ["./roc_interval.xls", "xls", "ROC曲线置信区间", 0, "110179"],
            ["./best_loc.xls", "xls", "ROC曲线临界值", 0, "110181"],
            ["./roc_curve_smooth.xls", "xls", "ROC曲线平滑处理结果表", 0, "110238"],
            ["./ROC曲线图.pdf", "pdf", "ROC分析结果图", 0, ""]
        ])
        super(RocWorkflow, self).end()

    def set_db(self):
        api_roc = self.api.roc
        datacurve = self.output_dir + '/roc_curve.xls'
        api_roc.add_roc_curve(roc_id=self.option("roc_id"), type="", file=datacurve)

        if os.path.exists(self.output_dir + '/roc_curve_smooth.xls'):
            datacurve_s = self.output_dir + '/roc_curve_smooth.xls'
            api_roc.add_roc_curve(roc_id=self.option("roc_id"), type="smooth", file=datacurve_s)

        if os.path.exists(self.output_dir + '/roc_interval.xls'):
            datainterval = self.output_dir + '/roc_interval.xls'
            api_roc.add_roc_interval(roc_id=self.option("roc_id"), file=datainterval)

        dataauc = self.output_dir + '/roc_auc.xls'
        api_roc.add_roc_auc(roc_id=self.option("roc_id"), file=dataauc, type="")

        if os.path.exists(self.output_dir + '/roc_auc_smooth.xls'):
            dataauc_s = self.output_dir + '/roc_auc_smooth.xls'
            api_roc.add_roc_auc(roc_id=self.option("roc_id"), file=dataauc_s, type="smooth")
        databestloc = self.output_dir + '/best_loc.xls'
        api_roc.add_roc_best_loc(roc_id=self.option("roc_id"), file=databestloc)

        report_files = [datacurve, dataauc, databestloc]
        for f in report_files:
            if not os.path.isfile(f):
                self.logger.error("找不到报告文件:{}".format(f))
                self.set_error("找不到报告文件", code="12703901")

        #self.end()
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("roc_id"), "sg_roc")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("roc_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "roc",
                "interaction": 1,
                "main_table": "sg_roc",
            })
            self.figsave.run()
        else:
            self.end()

    def run(self):
        self.run_roc()
        super(RocWorkflow, self).run()        

