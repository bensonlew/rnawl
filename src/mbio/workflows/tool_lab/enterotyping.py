# -*- coding: utf-8 -*-

"""菌群分型分析"""
import os
import re
import json
import shutil
import datetime
import types
from bson import ObjectId
import numpy as np
from mbio.packages.beta_diversity.filter_newick import get_level_newicktree
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.tool_lab.common_function import rename_name,rename_name_back
import gevent
from mbio.packages.metaasv.common_function import link_dir,link_file


class EnterotypingWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(EnterotypingWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "group_table", "type": "infile", "format": "tool_lab.group_table"},
            {"name": "distance_method", "type": "string", "default": "JSD"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.dist = self.add_tool("meta.beta_diversity.distance_calc")
        self.enterotyping = self.add_tool("tool_lab.mg_enterotyping")
        self.plot_enterotyping = self.add_tool("tool_lab.mg_plot_enterotyping")
        # self.enterotyping = self.add_tool("meta.beta_diversity.enterotyping")
        self.a = ''
        self.spe_name = ''
        self.number = ''

    def run_enterotyping(self):
        self.otutable = self.option("table")
        self.otutable.get_info()
        if self.otutable.prop['sample_num'] < 10:
            raise OptionError('样品数必须大于等于10', code="12701101")
        if self.otutable.prop['otu_num'] < 2:
            raise OptionError('物种/功能/基因数必须大于等于2', code="12701102")
        if self.option("distance_method") == "JSD":
            self.enterotyping.set_options({
                "otu_table": self.data_table,
            })
        else:
            self.enterotyping.set_options({
                "otu_table": self.data_table,
                "dis_matrix": self.dist.option('dis_matrix')
            })
        self.enterotyping.on("end", self.set_plot_options)
        self.enterotyping.run()

    def run_dist(self):
        options = {
            'method': self.option('distance_method'),
            'otutable': self.data_table
        }
        self.dist.set_options(options)
        self.dist.on("end", self.run_enterotyping)
        self.dist.run()

    def filter_otu_sample(self, otu_path, filter_samples, newfile):
        if not isinstance(filter_samples, types.ListType):
            self.logger.error('过滤otu表样本的样本名称应为列表')
            self.set_error("filter_otu_sample错误", code="12701105")
        try:
            with open(otu_path, 'rb') as f, open(newfile, 'wb') as w:
                one_line = f.readline()
                all_samples = one_line.rstrip().split('\t')[1:]
                if not ((set(all_samples) & set(filter_samples)) == set(filter_samples)):
                    self.logger.error('提供的过滤样本存在otu表中不存在的样本all:%s,filter_samples:%s' % (all_samples, filter_samples))
                    self.set_error("含有otu表中不存在的样本", code="12701106")
                if len(all_samples) == len(filter_samples):
                    return otu_path
                samples_index = [all_samples.index(i) + 1 for i in filter_samples]
                w.write('#OTU\t' + '\t'.join(filter_samples) + '\n')
                for line in f:
                    all_values = line.rstrip().split('\t')
                    new_values = [all_values[0]] + [all_values[i] for i in samples_index]
                    w.write('\t'.join(new_values) + '\n')
                return newfile
        except IOError:
            self.set_error('无法打开OTU相关文件或者文件不存在', code="12701107")

    def _get_samplenames(self, groupfile):
        try:
            with open(groupfile, 'rb') as f:
                alllines = f.readlines()
                all_names = [i.split('\t')[0] for i in alllines]
            return all_names[1:]
        except IOError:
            self.set_error('无法打开分组文件或者文件不存在', code="12701108")

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
        spe_name_re = []
        for i in name:
            path_c = "/" + i
            print(all_path + path_c)
            if os.path.exists(all_path + path_c):
                rename_name_back(all_path + path_c,self.name_dict)
                b = open(all_path + path_c, "r")
                line = b.readline()
                line = b.readline()
                f = line.strip().split("\t")
                spe_name_re.append(f[0].split("; ")[-1])
                b.close()
            else:
                break
        spe_names = ';'.join(spe_name_re)
        print(spe_names)
        self.spe_name = spe_names
        self.run_plot_enterotyping()

    def run_plot_enterotyping(self):
        with open(self.work_dir + "/tmp_group.txt","w") as t,open(self.option("group_table").prop["path"]) as f:
            data = f.readlines()
            t.write("#sample\tGroup1\n")
            for i in data[1:]:
                t.write(i)

        option = ({
            "otu_table": self.data_table,
            "g": self.a,
            "s": self.a,
            "group": self.work_dir + "/tmp_group.txt",
        })
        if self.option("distance_method") != "JSD":
            option["dis_matrix"] = self.dist.option('dis_matrix')
        self.plot_enterotyping.set_options(option)
        self.plot_enterotyping.on('end', self.set_db)
        self.plot_enterotyping.run()

    def set_db(self):
        with open(self.enterotyping.output_dir+"/cluster.txt") as f,open(self.work_dir + "/tmp_group.txt") as f2,open(self.plot_enterotyping.output_dir+"/summary.txt","w") as t:
            group_dict = {}
            type_dict = {}
            group_list = []
            data1 = f.readlines()
            data2 = f2.readlines()
            for i in data2[1:]:
                if i.strip().split("\t")[1] in group_dict:
                    group_dict[i.strip().split("\t")[1]].append(i.strip().split("\t")[0])
                else:
                    group_list.append(i.strip().split("\t")[1])
                    group_dict[i.strip().split("\t")[1]] = [i.strip().split("\t")[0]]

            for x in data1[1:]:
                if x.strip().split("\t")[1] in type_dict:
                    type_dict[x.strip().split("\t")[1]].append(x.strip().split("\t")[0])
                else:
                    type_dict[x.strip().split("\t")[1]] = [x.strip().split("\t")[0]]
            t.write("Enterotype"+"\t"+"\t".join(group_list)+"\n")
            for type in type_dict:
                t.write(str(type))
                group_info = {}
                for group in group_list:
                    num = 0
                    for xx in type_dict[type]:
                        if xx in group_dict[group]:
                            num += 1
                        group_info[group] = num
                self.logger.info("group_info{}".format(group_info))
                for xxx in group_list:
                    t.write("\t" + str(group_info[xxx]))
                t.write("\n")

        link_dir(self.enterotyping.output_dir, self.output_dir)
        link_file(self.plot_enterotyping.output_dir + "/summary.txt", self.output_dir + "/summary.txt")
        self.logger.info("正在写入mongo数据库")
        name = os.listdir(self.enterotyping.output_dir)
        name.remove("ch.txt")
        name.remove("cluster.txt")
        if "summary.txt" in name:
            name.remove("summary.txt")
        api_enterotype = self.api.api("tool_lab.enterotyping")
        top1_name = {}
        for file11 in os.listdir(self.output_dir):
            if file11.endswith(".cluster.txt"):
                with open(self.output_dir+"/"+file11)as f:
                    data = f.readlines()
                    top1_name["type"+file11.rstrip(".cluster.txt")] = data[1].strip().split("\t")[0]
        # 更新主表
        if os.path.exists(self.plot_enterotyping.output_dir + "/BCA_circle.txt"):
            self.stats = True
        else:
            self.stats = False
        main_id = api_enterotype.add_enterotype(main=False, main_id=self.option('main_id'), cluster_name=self.a, spe_name=self.spe_name,stats=self.stats)

        for i in name:
            api_enterotype.add_enterotype_detail_cluster(main_id, self.enterotyping.output_dir + "/" + i,
                                                         name=i.split(".")[0], update_type=self.a)

        api_enterotype.add_enterotype_bar(self.output_dir + "/summary.txt", "enterotype", update_id=main_id,
                                          update_column=True,top1_name=str(top1_name))
        api_enterotype.add_enterotype_bar(self.output_dir + "/ch.txt", "ch", update_id=main_id, update_column=True)
        api_enterotype.add_enterotype_detail(self.output_dir + "/cluster.txt", "ch", update_id=main_id,
                                             update_column=True,group_table=self.option("group_table").prop["path"])
        for table_type in ["BCA_circle", "BCA_point", "pcoa_circle", "pcoa_point"]:
            if os.path.exists(self.plot_enterotyping.output_dir + "/" + table_type + ".txt"):
                api_enterotype.add_enterotype_scatter(
                    self.plot_enterotyping.output_dir + "/" + table_type + ".txt",
                    table_type=table_type,
                    update_id=main_id,
                    update_column=True,
                    group_data=self.option("group_table").prop["path"],
                    type_data=self.output_dir + "/cluster.txt",top1_name=str(top1_name))
        gevent.sleep(3)
        """
        api_group = self.api.api("metaasv.group")
        if self.sheet.output:
            if self.sheet.output.endswith("/"):
                group_name = [self.sheet.output.split("/")[-2]]
            else:
                group_name = [self.sheet.output.split("/")[-1]]
            self.logger.info("group: {}".format(group_name))
            api_group.add_ini_group_table(self.output_dir + "/cluster.txt", task_id=self._sheet.id,
                                          group_name=group_name, add="Type")
        """
        self.end()

    def end(self):
        if os.path.exists(self.output_dir + "/summary.txt"):
            os.remove(self.output_dir + "/summary.txt")
        if os.path.exists(self.output_dir + "/enterotyping"):
            shutil.rmtree(self.output_dir + "/enterotyping")
        try:
            shutil.copy2(self.plot_enterotyping.output_dir + "/summary.txt", self.output_dir + "/summary.txt")
            #shutil.copytree(self.enterotyping.output_dir, self.output_dir + "/enterotyping")
        except Exception as e:
            self.logger.error("summary.txt copy failed{}".format(e))
            self.set_error("summary.txt copy failed", code="12701109")
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "样本菌群分型分析结果输出目录", 0, "110129"],
            ["./summary.txt", "txt", "summary数据表", 0, "110134"],
            ["./ch.txt", "txt", "CH指数数据表", 0, "110132"],
            ["./cluster.txt", "txt", "cluster数据表", 0, "110133"]
        ])
        result_dir.add_regexp_rules([
            [".+\cluster.txt$", "txt", "分型后各组数据表", 0, "110131"]
        ])
        super(EnterotypingWorkflow, self).end()

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.data_table_raw = self.option("table").prop["path"]
        self.data_table = self.work_dir + "/input_table.txt"
        self.name_dict = rename_name(self.data_table_raw, self.data_table)
        if self.option("distance_method") == "JSD":
            self.run_enterotyping()
        else:
            self.run_dist()
        super(EnterotypingWorkflow, self).run()
