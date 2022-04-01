# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'

"""样本菌群分型分析模块"""
import os
import re
import json
import shutil
import datetime
import numpy as np
from biocluster.workflow import Workflow
from mbio.packages.metagenomic.id_convert import name2id
from biocluster.core.exceptions import OptionError
from comm_table import CommTableWorkflow
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name
import json


class EnterotypingWorkflow(CommTableWorkflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(EnterotypingWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "anno_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "group_id", "type": "string", "default": ""},
            {"name": "group_detail", "type": "string", "default": ""},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_id", "type": "string", "default": ""},
            {"name": "distance_method", "type": "string", "default": "JSD"},
            {"name": "profile_table", "type": "infile", "format": "meta.otu.otu_table"},
            # {'name': 'save_pdf', 'type': 'int', 'default': 1}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.dist = self.add_tool("meta.beta_diversity.distance_calc")
        #self.abundance = self.add_tool("meta.create_abund_table")
        self.enterotyping = self.add_tool("meta.beta_diversity.mg_enterotyping")
        self.plot_enterotyping = self.add_tool("meta.beta_diversity.mg_plot_enterotyping")
        self.sam = self.add_tool("meta.otu.sort_samples_mg")
        self.add_option(options)
        self.set_options(self._sheet.options())
        # self.enterotyping = self.add_tool("meta.beta_diversity.enterotyping")
        self.a = ''
        self.spe_name = ''
        self.number = ''
        # # self.name = ''
        # group_table_path = os.path.join(self.work_dir, "group_table.xls")
        # self.group_table_path = Meta().group_detail_to_table(self.option("group_detail"), group_table_path)

    # def check_options(self):
    #     if self.option('method') not in ['average', 'single', 'complete', ""]:
    #         raise OptionError('错误的层级聚类方式：%s' % self.option('method'))

    def run_enterotyping(self):
        self.otutable = self.sam.option("out_otu_table")
        self.otutable.get_info()
        if self.otutable.prop['sample_num'] < 10:
            raise OptionError('样品数必须大于等于10', code="12801401")
        if self.otutable.prop['otu_num'] < 2:
            raise OptionError('物种/功能/基因数必须大于等于2', code="12801402")
        if self.option("distance_method") == "JSD":
            self.enterotyping.set_options({
                "otu_table": self.otutable,
                # "group_table": self.group_table_path
            })
        else:
            self.enterotyping.set_options({
                "otu_table": self.otutable,
                "dis_matrix": self.dist.option('dis_matrix')
            })
        self.enterotyping.on("end", self.set_plot_options)
        self.enterotyping.run()

    def run_dist(self):
        otutable = self.sam.option("out_otu_table")
        options = {
            'method': self.option('distance_method'),
            'otutable': otutable
        }
        self.dist.set_options(options)
        self.dist.on("end", self.run_enterotyping)
        self.dist.run()

    def set_plot_options(self):
        all_path = self.enterotyping.output_dir
        print(all_path)
        cluster_path = all_path + "/cluster.txt"
        print(cluster_path)
        up_num = []
        if os.path.exists(cluster_path):
            print("true")
            a = open(cluster_path, "r")
            content = a.readlines()
            for f in content:
                if f.startswith("name") == False:
                    f = f.strip().split("\t")
                    up_num.append(f[1])
            a.close()
        up_num.sort()
        name = os.listdir(all_path)
        self.number = len(name) - 1
        name.remove("ch.txt")  # 1.3
        name.remove("cluster.txt")
        name.sort()
        print(name)
        the_big = len(name)
        print(the_big)
        a = []
        for i in range(1, int(the_big) + 1):
            n = str(i)
            a.append(n)
        a = ','.join(a)
        print(a)
        self.a = a
        spe_name = []
        spe_name_re = []
        for i in name:
            path_c = "/" + i
            print(all_path + path_c)
            if os.path.exists(all_path + path_c):
                b = open(all_path + path_c, "r")
                line = b.readline()
                line = b.readline()
                f = line.strip().split("\t")
                spe_name_re.append(f[0])
                # for f in content:
                #     if f.startswith("taxon_name") == False:
                #         f = f.strip().split("\t")
                #         spe_name_re.append(f[0])
                        # if f[0] not in spe_name_re:
                        #     spe_name_re.append(f[0])
                        #     t = f[0].split(" ")
                        #     spe_name.append(t[-1])
                        #     break
                        # else:
                        #     continue
                b.close()
            else:
                break
        spe_names = ';'.join(spe_name_re)
        print(spe_names)
        self.spe_name = spe_names
        self.run_plot_enterotyping()

    def run_plot_enterotyping(self):
        option = ({
            "otu_table": self.otutable,
            "g": self.a,
            "s": self.a,
            "group": self.option("group_table"),
        })
        if self.option("distance_method") != "JSD":
            option["dis_matrix"] = self.dist.option('dis_matrix')
        self.plot_enterotyping.set_options(option)
        self.plot_enterotyping.on('end', self.set_db)
        self.plot_enterotyping.run()

    def sort_sample(self):
        if self.option("profile_table").is_set:
            otutable = self.option("profile_table")
        else:
            otutable = self.abundance.option('out_table')
        options = {
            'in_otu_table': otutable,
            'group_table': self.option("group_table")
        }
        self.sam.set_options(options)
        self.sam.run()

    def set_db(self):
        self.logger.info("正在写入mongo数据库")
        main_id = self.option('main_id')
        name = os.listdir(self.enterotyping.output_dir)
        name.remove("ch.txt")
        name.remove("cluster.txt")
        api_enterotype = self.api.api('metagenomic.enterotype')
        api_enterotype.add_enterotype(main=False, main_id=main_id, cluster_name=self.a,
                                      spe_name=self.spe_name)
        for i in name:
            api_enterotype.add_enterotype_detail_cluster(self.option("main_id"), self.enterotyping.output_dir + "/" + i,
                                                         name=i.split(".")[0])
        api_enterotype.add_enterotype_detail(self.plot_enterotyping.output_dir + "/summary.txt", "summary",
                                             update_id=main_id)
        for table_type in ["BCA_circle", "BCA_point", "pcoa_circle", "pcoa_point"]:
            if os.path.exists(self.plot_enterotyping.output_dir + "/" + table_type + ".txt"):
                api_enterotype.add_enterotype_detail(
                    self.plot_enterotyping.output_dir + "/" + table_type + ".txt",
                    table_type=table_type,
                    update_id=main_id)
        api_enterotype.add_enterotype_detail(self.enterotyping.output_dir + "/ch.txt", "ch", update_id=main_id,
                                             update_column=False)
        api_enterotype.add_enterotype_detail(self.enterotyping.output_dir + "/cluster.txt", "cluster",
                                             update_id=main_id, update_column=False)
        task_name2id = self.option("task_id")
        self.sample_2_id = name2id(task_name2id, type="task")
        api_group = self.api.api('metagenomic.specimen_group')
        if self.sheet.output:
            if self.sheet.output.endswith("/"):
                group_name = [self.sheet.output.split("/")[-2]]
            else:
                group_name = [self.sheet.output.split("/")[-1]]
            api_group.add_ini_group_table(self.enterotyping.output_dir + "/cluster.txt", self.sample_2_id,
                                          group_name=group_name, task_id=self.option("task_id"), add="Enterotype")
        # if self.option("save_pdf"):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "enterotype")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            submit_loc = json.loads(self.option('params'))["submit_location"]
            # self.logger.info(submit_loc)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "metagenomic",
                "submit_loc": submit_loc,
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        # if self.option("save_pdf"):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        if os.path.exists(self.output_dir + "/summary.txt"):
            os.remove(self.output_dir + "/summary.txt")
        if os.path.exists(self.output_dir + "/enterotyping"):
            shutil.rmtree(self.output_dir + "/enterotyping")
        try:
            shutil.copy2(self.plot_enterotyping.output_dir + "/summary.txt", self.output_dir + "/summary.txt")
            shutil.copytree(self.enterotyping.output_dir, self.output_dir + "/enterotyping")
        except Exception as e:
            self.logger.info("summary.txt copy success{}".format(e))
            raise Exception("summary.txt copy success{}".format(e))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "样本菌群分型分析结果输出目录", 0, "120175"],
            ["./summary.txt", "txt", "summary数据表", 0, "120176"],
            ["./enterotyping", "dir", "分型数据文件夹", 0, "120177"],
            ["./enterotyping/ch.txt", "txt", "CH指数数据表", 0, "120178"],
            ["./enterotyping/cluster.txt", "txt", "cluster数据表", 0, "120179"],
            ["./Enterotypes.pdf", "pdf", "分型图"],
            ["./CH.pdf", "pdf", "CH指数图"],
            ["./summary.pdf", "pdf", "各组分型结果柱形图"]
        ])
        result_dir.add_regexp_rules([
            ["enterotyping/.+\cluster.txt$", "txt", "分型后各组数据表", 0, "120180"]
        ])
        super(EnterotypingWorkflow, self).end()

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        if self.option("profile_table").is_set:
            if self.option("distance_method") == "JSD":
                self.sam.on('end', self.run_enterotyping)
            else:
                self.sam.on('end', self.run_dist)
        else:
            if self.option("distance_method") == "JSD":
                #self.abundance.on('end', self.sort_sample)
                self.run_abundance(self.sort_sample)
                self.sam.on('end', self.run_enterotyping)
            else:
                #self.abundance.on('end', self.sort_sample)
                self.run_abundance(self.sort_sample)
                self.sam.on('end', self.run_dist)
        if self.option("profile_table").is_set:
            self.sort_sample()
        else:
            #self.run_abundance()
            self.abundance.run()
        super(EnterotypingWorkflow, self).run()
