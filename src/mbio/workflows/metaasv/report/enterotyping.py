# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os
import re
import shutil
import gevent
import types
from mbio.packages.beta_diversity.filter_newick_asv import get_level_newicktree
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir,link_file


class EnterotypingWorkflow(Workflow):
    """
    metaasv 样本菌群分型分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(EnterotypingWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "asv_id", "type": "string"},##主表asv的id
            {"name": "main_id", "type": "string"},##主表id
            {"name": "otu_file", "type": "infile", "format": "meta.otu.otu_table"},##ASV表格
            {"name": "level", "type": "int"},##分类学水平
            {"name": "group_id", "type": "string", "default": ""},
            {"name": "group_detail", "type": "string", "default": ""},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},##分组文件
            {"name": "distance_method", "type": "string", "default": "JSD"}##距离计算方法
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.dist = self.add_tool("meta.beta_diversity.distance_calc")
        self.enterotyping = self.add_tool("meta.beta_diversity.mg_enterotyping")
        self.plot_enterotyping = self.add_tool("meta.beta_diversity.mg_plot_enterotyping")
        self.a = ''
        self.spe_name = ''
        self.number = ''

    def run_enterotyping(self):
        """
        运行样本分析 tool mg_enterotyping
        :return:
        """
        self.otutable = self.option("otu_file")
        self.otutable.get_info()
        if self.otutable.prop['sample_num'] < 10:
            raise OptionError('样品数必须大于等于10')
        if self.otutable.prop['otu_num'] < 2:
            raise OptionError('物种/功能/基因数必须大于等于2')
        if self.option("distance_method") == "JSD":
            self.enterotyping.set_options({
                "otu_table": self.otutable,
            })
        else:
            self.enterotyping.set_options({
                "otu_table": self.otutable,
                "dis_matrix": self.dist.option('dis_matrix')
            })
        self.enterotyping.on("end", self.set_plot_options)
        self.enterotyping.run()

    def run_dist(self):
        """
        计算距离
        :return:
        """
        otutable = self.option("otu_file")
        if 'unifrac' in self.option(
                'distance_method'):  # sanger_bioinfo/src/mbio/workflows/meta/report/distance_calc.py中的解释
            if self.option('level') != 9:
                newicktree = get_level_newicktree(self.option('asv_id'), level=self.option('level'),tempdir=self.work_dir, return_file=False, bind_obj=self)
                all_find = re.findall(r'\'.+?\'', newicktree)
                for n, m in enumerate(all_find):
                    all_find[n] = m.strip('\'')
                all_find = dict((i[1], i[0]) for i in enumerate(all_find))

                def match_newname(matchname):
                    if hasattr(match_newname, 'count'):
                        match_newname.count = match_newname.count + 1
                    else:
                        match_newname.count = 1
                    return 'ASV' + str(match_newname.count)

                newline = re.sub(r'\'.+?\'', match_newname, newicktree)
                temp_tree_file = self.work_dir + '/temp.tree'
                tempfile = open(temp_tree_file, 'w')
                tempfile.write(newline)
                tempfile.close()
                otu_table = self.option('otu_file').path
                temp_otu_file = self.option('otu_file').path + '.temp'
                all_lines = open(otu_table, 'r').readlines()
                if len(all_lines) < 3:
                    self.logger.error('分类水平：%s,asv表数据少于2行：%s' % (self.option('level'), len(all_lines)))
                    self.set_error("asv表数据少于2行", code="12701103")
                new_all = []
                new_all.append(all_lines[0])
                for line in all_lines[1:]:
                    name = line.split('\t')
                    if name[0] not in all_find:
                        self.logger.error('ASV表中存在不是直接通过组合原始表分类名称的ASV名：%s' % name[0])
                        self.set_error("ASV表中有原始表中不存在的ASV名：%s", variables=(name[0]))
                    name[0] = 'ASV' + str(all_find[name[0]] + 1)
                    new_all.append('\t'.join(name))
                otu_file_temp = open(temp_otu_file, 'w')
                otu_file_temp.writelines(new_all)
                otu_file_temp.close()
                options = {
                    'method': self.option('distance_method'),
                    'otutable': temp_otu_file,
                    'newicktree': temp_tree_file
                }
                options['otutable'] = self.filter_otu_sample(options['otutable'],
                                                             self._get_samplenames(self.option('group_table').path),
                                                             options['otutable'] + '.temp')
            else:
                newicktree = get_level_newicktree(self.option('asv_id'), level=self.option('level'),
                                                  tempdir=self.work_dir, return_file=False, bind_obj=self)
                temp_tree_file = self.work_dir + '/temp.tree'
                tempfile = open(temp_tree_file, 'w')
                tempfile.write(newicktree)
                tempfile.close()
                otu_table = self.option('otu_file').path
                temp_otu_file = self.option('otu_file').path + '.temp'
                all_lines = open(otu_table, 'r').readlines()
                new_all = []
                new_all.append(all_lines[0])
                for line in all_lines[1:]:  # ASV表中有复杂的名称ASV名称，包含进化物种类型，进化树种只有ASV名称
                    name = line.split('\t')
                    name[0] = name[0].split(';')[-1].strip()
                    new_all.append('\t'.join(name))
                otu_file_temp = open(temp_otu_file, 'w')
                otu_file_temp.writelines(new_all)
                otu_file_temp.close()
                options = {
                    'method': self.option('distance_method'),
                    'otutable': temp_otu_file,
                    'newicktree': temp_tree_file
                }
                options['otutable'] = self.filter_otu_sample(options['otutable'],
                                                             self._get_samplenames(self.option('group_table').path),
                                                             options['otutable'] + '.temp')
        else:
            options = {
                'method': self.option('distance_method'),
                'otutable': otutable
            }
        self.dist.set_options(options)
        self.dist.on("end", self.run_enterotyping)
        self.dist.run()

    def filter_otu_sample(self, otu_path, filter_samples, newfile):
        if not isinstance(filter_samples, types.ListType):
            self.logger.error('过滤otu表样本的样本名称应为列表')
            self.set_error("filter_otu_sample错误")
        try:
            with open(otu_path, 'rb') as f, open(newfile, 'wb') as w:
                one_line = f.readline()
                all_samples = one_line.rstrip().split('\t')[1:]
                if not ((set(all_samples) & set(filter_samples)) == set(filter_samples)):
                    self.logger.error('提供的过滤样本存在asv表中不存在的样本all:%s,filter_samples:%s' % (all_samples, filter_samples))
                    self.set_error("含有asv表中不存在的样本")
                if len(all_samples) == len(filter_samples):
                    return otu_path
                samples_index = [all_samples.index(i) + 1 for i in filter_samples]
                w.write('ASV\t' + '\t'.join(filter_samples) + '\n')
                for line in f:
                    all_values = line.rstrip().split('\t')
                    new_values = [all_values[0]] + [all_values[i] for i in samples_index]
                    w.write('\t'.join(new_values) + '\n')
                return newfile
        except IOError:
            self.set_error('无法打开OTU相关文件或者文件不存在')

    def _get_samplenames(self, groupfile):
        try:
            with open(groupfile, 'rb') as f:
                alllines = f.readlines()
                all_names = [i.split('\t')[0] for i in alllines]
            return all_names[1:]
        except IOError:
            self.set_error('无法打开分组文件或者文件不存在')

    def set_plot_options(self):
        """
        为画图设置参数
        :return:
        """
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
                b = open(all_path + path_c, "r")
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
        """
        分型画图
        :return:
        """
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

    def set_db(self):
        """
        连接结果文件和导入MongoDB
        :return:
        """
        link_dir(self.enterotyping.output_dir, self.output_dir)
        link_file(self.plot_enterotyping.output_dir + "/summary.txt", self.output_dir + "/summary.txt")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option('main_id')
        name = os.listdir(self.enterotyping.output_dir)
        name.remove("ch.txt")
        name.remove("cluster.txt")
        if "summary.txt" in name:
            name.remove("summary.txt")
        api_enterotype = self.api.api("metaasv.enterotype")
        # 更新主表
        api_enterotype.add_enterotype(main=False, main_id=main_id, cluster_name=self.a,spe_name=self.spe_name)

        for i in name:
            api_enterotype.add_enterotype_detail_cluster(self.option("main_id"), self.enterotyping.output_dir + "/" + i,
                                                         name=i.split(".")[0], update_type=self.a)

        api_enterotype.add_enterotype_bar(self.output_dir + "/summary.txt", "enterotype",update_id=main_id, update_column=True)
        api_enterotype.add_enterotype_bar(self.output_dir + "/ch.txt", "ch",update_id=main_id, update_column=True)
        api_enterotype.add_enterotype_detail(self.output_dir + "/cluster.txt", "ch", update_id=main_id, update_column=True)
        for table_type in ["BCA_circle", "BCA_point", "pcoa_circle", "pcoa_point"]:
            if os.path.exists(self.plot_enterotyping.output_dir + "/" + table_type + ".txt"):
                api_enterotype.add_enterotype_scatter(
                    self.plot_enterotyping.output_dir + "/" + table_type + ".txt",
                    table_type=table_type,
                    update_id=main_id,
                    update_column=True)
        gevent.sleep(3)
        api_group = self.api.api("metaasv.group")
        if self.sheet.output:
            if self.sheet.output.endswith("/"):
                group_name = [self.sheet.output.split("/")[-2]]
            else:
                group_name = [self.sheet.output.split("/")[-1]]
            self.logger.info("group: {}".format(group_name))
            api_group.add_ini_group_table(self.output_dir + "/cluster.txt",task_id=self._sheet.id, group_name=group_name, add="Type")

        self.end()

    def end(self):
        """
        结束和上传结果文件
        :return:
        """
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "样本菌群分型分析结果输出目录", 0, ""],
            ["./summary.txt", "txt", "summary数据表", 0, ""],
            ["./enterotyping", "dir", "分型数据文件夹", 0, ""],
            ["./enterotyping/ch.txt", "txt", "CH指数数据表", 0, ""],
            ["./enterotyping/cluster.txt", "txt", "cluster数据表", 0, ""]
        ])
        result_dir.add_regexp_rules([
            ["enterotyping/.+\cluster.txt$", "txt", "分型后各组数据表", 0, ""]
        ])
        super(EnterotypingWorkflow, self).end()

    def run(self):
        """
        运行
        :return:
        """
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        if self.option("distance_method") == "JSD":
            self.run_enterotyping()
        else:
            self.run_dist()
        super(EnterotypingWorkflow, self).run()
