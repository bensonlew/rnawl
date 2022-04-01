# -*- coding: utf-8 -*-
# __author__ = 'zhangpeng'

""""""

import datetime
from biocluster.workflow import Workflow
import re
import os
import json
import shutil
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class NPcaWorkflow(Workflow):
    """
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(NPcaWorkflow, self).__init__(wsheet_object)
        options = [
        
            {"name": "otu_table", "type": "infile","format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "level", "type": "int", "default": 9},
            {"name": "second_group_table", "type": "infile","format": "meta.otu.group_table"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            #{"name": "env_labs", "type": "string", "default": ""},
            #{"name": "PCAlabs", "type": "string", "default": ""},
            {"name": "n_pca_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "second_group_id", "type": 'string'},
            {"name": "otu_id", "type": 'string'},
            {"name": "group_detail", "type": "string"},
            {"name": "second_group_detail", "type": "string"},
            {"name": "group_id","type":"string"},
        ]
        self.add_option(options)

        self.set_options(self._sheet.options())
        self.n_pca = self.add_tool("meta.beta_diversity.n_pca")
        self.output_dir = self.n_pca.output_dir

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
                    #line[0] = line_data[-1]
                    for i in range(0, len(line)):
                        if i == len(line)-1:
                            f2.write("%s\n"%(line[i]))
                        else:
                            f2.write("%s\t"%(line[i]))
        f2.close()
        return newtable



    def run_n_pca(self):
        newtable = self.change_otuname(self.option('otu_table').prop['path'])
        options = {
            'otu_table': newtable,
            'level': self.option('level'),
            'second_group_table':self.option('second_group_table'),
            'group_table':self.option('group_table'),            
        }
        self.n_pca.set_options(options)
        self.n_pca.on('end',self.set_db)
        self.output_dir = self.n_pca.output_dir
        self.n_pca.run()

    def end(self):
        if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "npca分析结果目录", 0, "110119"],
            ["./sites.xls", "xls", "坐标数据", 0, "110120"],
            ["./sd.xls", "xls", "方差大小", 0, "110122"],
            ["./sdmax.xls", "xls", "置信上边界", 0, "110126"],
            ["./sdmin.xls", "xls", "置信下边界", 0, "110121"],
            ["./rotation_mean.xls", "xls", "平均值", 0, "110124"],
            ["./importance.xls", "xls", "百分率返回", 0, "1110123"],
            ["./sitesall.xls", "xls", "所有信息", 0, "110125"],
            ["./多维PCA分析图.pdf", "pdf", "样本多维PCA分析图", 0, ""]
        ])
        super(NPcaWorkflow, self).end()

    def set_db(self):
        api_n_pca = self.api.n_pca
        datasite = self.output_dir + '/rotation_mean.xls'
        datamin = self.output_dir + '/sdmin.xls'
        datamax = self.output_dir + '/sdmax.xls'
        dataimportance = self.output_dir + '/importance.xls'
        datasitesall = self.output_dir + '/sitesall.xls'
        if not os.path.isfile(datasite):
            self.logger.error("找不到报告文件:{}".format(datasite))
            self.set_error("找不到报告文件", code="12702801")
        if not os.path.isfile(datamin):
            self.logger.error("找不到报告文件:{}".format(datamin))
            self.set_error("找不到报告文件", code="12702801")
        if not os.path.isfile(datamax):
            self.logger.error("找不到报告文件:{}".format(datamax))
            self.set_error("找不到报告文件", code="12702801")
        if not os.path.isfile(dataimportance):
            self.logger.error("找不到报告文件:{}".format(dataimportance))
            self.set_error("找不到报告文件", code="12702801")
        if not os.path.isfile(datasitesall):
            self.logger.error("找不到报告文件:{}".format(datasiteall))
            self.set_error("找不到报告文件", code="12702801")
        #api_n_pca.add_n_pca_site(file_path=datasite, table_id=self.option("n_pca_id"))
        #api_n_pca.add_n_pca_min(file_path=datamin,table_id=self.option("n_pca_id"))
        #api_n_pca.add_n_pca_max(file_path=datamax,table_id=self.option("n_pca_id"))
        api_n_pca.add_n_pca_importance(file_path=dataimportance, table_id=self.option("n_pca_id"))
        api_n_pca.add_n_pca_sitesall(file_path=datasitesall, table_id=self.option("n_pca_id"))
        #self.end()
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("n_pca_id"),"sg_npca")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("n_pca_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "beta_npca",
                "interaction": 1,
                "main_table": "sg_npca",
            })
            self.figsave.run()
        else:
            self.end()

    def run(self):
        self.run_n_pca()
        super(NPcaWorkflow, self).run()        

