# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.workflow import Workflow
from mbio.packages.metaasv.filter_newick import get_level_newicktree
from mbio.packages.metaasv.common_function import link_dir, link_file
import re
import os


class HclusterWorkflow(Workflow):
    """
    metaasv 样本层级聚类
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(HclusterWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "dist_method", "type": "string", "default": 'bray_curtis'},
            {"name": "hcluster_method", "type": "string", "default": 'average'},
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "level", "type": 'int', "default": 9},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "asv_id", "type": "string"},
            {"name": "others_value", "type": "float", "default": ""}
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.dist = self.add_tool("meta.beta_diversity.distance_calc")
        self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.hcluster = self.add_tool("meta.beta_diversity.hcluster")

    def run(self):
        """
        运行
        """
        if self.option("others_value") != "":
            self.sort_samples.set_options({
                "in_otu_table": self.option("otu_table"),
                "others": self.option("others_value")
            })
            self.sort_samples.on("end", self.run_dist)
            self.sort_samples.run()
        else:
            self.run_dist()
        super(HclusterWorkflow, self).run()

    def run_dist(self):
        """
        计算距离
        """
        if 'unifrac' in self.option('dist_method'):
            # 查找OTU表对应的进化树
            if self.option('level') != 9:
                newicktree = get_level_newicktree(self.option('asv_id'), level=self.option('level'),
                                                  tempdir=self.work_dir, return_file=False, bind_obj=self)
                all_find = re.findall(r'\'.+?\'', newicktree)  # 找到所有带引号的进化树中复杂的名称
                for n, m in enumerate(all_find):
                    all_find[n] = m.strip('\'')
                all_find = dict((i[1], i[0]) for i in enumerate(all_find))  # 用名称做键，找到的位置数字做值

                def match_newname(matchname):
                    '随着自身被调用，自身的属性count随调用次数增加，返回OTU加次数，用于重命名进化树复杂的名称'
                    if hasattr(match_newname, 'count'):
                        match_newname.count = match_newname.count + 1
                    else:
                        match_newname.count = 1
                    return 'ASV' + str(match_newname.count)  # 后面替换OTU中名称用同样的命名规则
                newline = re.sub(r'\'.+?\'', match_newname, newicktree)  # 替换树种的复杂名称用 OTU 加数字代替 , 选哟注意的是这里的sub查找与findall查到方式是一致的
                temp_tree_file = self.work_dir + '/temp.tree'
                tempfile = open(temp_tree_file, 'w')
                tempfile.write(newline)
                tempfile.close()
                self.logger.info('get_newick:' + temp_tree_file)
                otu_table = self.option('otu_table').path
                temp_otu_file = self.option('otu_table').path + '.temp'
                all_lines = open(otu_table, 'r').readlines()
                if len(all_lines) < 3:
                    self.logger.error('分类水平：%s,otu表数据少于2行：%s' % (self.option('level'), len(all_lines)))
                    self.set_error("asv表数据少于2行")
                new_all = []
                new_all.append(all_lines[0])
                for line in all_lines[1:]:  # 遍历OTU表，将OTU表的复杂OTU名称改为之前find到的复杂名称对应的字典
                    name = line.split('\t')
                    if name[0] not in all_find:
                        self.set_error('ASV表中有原始表不存在的OTU名：%s', variables=(name[0]))
                    name[0] = 'ASV' + str(all_find[name[0]] + 1)
                    new_all.append('\t'.join(name))
                otu_file_temp = open(temp_otu_file, 'w')
                otu_file_temp.writelines(new_all)
                otu_file_temp.close()
                os.rename(otu_table, os.path.join(self.work_dir, "otu_table2.xls"))
                os.rename(temp_otu_file, otu_table)
                input_table = os.path.join(self.work_dir, "otu_table.xls")
                options = {
                    'method': self.option('dist_method'),
                    'otutable': input_table,
                    'newicktree': temp_tree_file
                }
            else:
                newicktree = get_level_newicktree(self.option('asv_id'), level=self.option('level'),
                                                  tempdir=self.work_dir, return_file=False, bind_obj=self)
                temp_tree_file = self.work_dir + '/temp.tree'
                tempfile = open(temp_tree_file, 'w')
                tempfile.write(newicktree)
                tempfile.close()
                otu_table = self.option('otu_table').path
                temp_otu_file = self.option('otu_table').path + '.temp'
                all_lines = open(otu_table, 'r').readlines()
                new_all = []
                new_all.append(all_lines[0])
                for line in all_lines[1:]:  # OTU表中有复杂的名称OTU名称，包含进化物种类型，进化树种只有OTU名称
                    name = line.split('\t')
                    name[0] = name[0].split(';')[-1].strip()
                    new_all.append('\t'.join(name))
                otu_file_temp = open(temp_otu_file, 'w')
                otu_file_temp.writelines(new_all)
                otu_file_temp.close()
                os.rename(otu_table, os.path.join(self.work_dir, "otu_table2.xls"))
                os.rename(temp_otu_file, otu_table)
                input_table = os.path.join(self.work_dir, "otu_table.xls")
                options = {
                    'method': self.option('dist_method'),
                    'otutable': input_table,
                    'newicktree': temp_tree_file
                }
        else:
            options = {
                'method': self.option('dist_method'),
                'otutable': self.option('otu_table')
            }
        self.dist.set_options(options)
        self.dist.on('end', self.run_hcluster)
        self.hcluster.on('end', self.set_db)
        self.dist.run()

    def run_hcluster(self):
        """
        计算层级距离矩阵
        """
        options = {
            'linkage': self.option('hcluster_method'),
            'dis_matrix': self.dist.option('dis_matrix')
        }
        self.hcluster.set_options(options)
        self.hcluster.run()

    def set_db(self):
        """
        保存结果树结果到mongo数据库中
        """
        link_dir(self.hcluster.output_dir, self.output_dir)
        distance_file = os.path.join(self.dist.output_dir, self.option("dist_method")+"_otu_table.xls")
        link_file(distance_file, os.path.join(self.output_dir, "hcluster"+"_"+self.option("dist_method")+".xls"))
        if self.option("others_value") != "" :
            rank_others_file = self.sort_samples.output_dir + "/taxa.percents.table.xls"
            link_file(rank_others_file, os.path.join(self.output_dir, "barplot_table.xls"))
        dir_path = self.output_dir
        api_pca = self.api.api("metaasv.beta_diversity")
        api_pca.add_beta_multi_analysis_result(dir_path, "hcluster",main=False,main_id=str(self.option('main_id')), group_file=self.option("group_file").prop['path'])
        self.logger.info("结束啦")
        self.end()

    def end(self):
        """
        结束和上传结果文件目录
        :return:
        """
        result_dir_hucluster = self.add_upload_dir(self.output_dir)
        result_dir_hucluster.add_relpath_rules([
            [".", "", "样本层级聚类分析结果目录", 0, ""],
            ["./hcluster.tre", "tre", "层级聚类树结果表", 0, ""]
            ])
        result_dir_hucluster.add_regexp_rules([
            [r'%s.*\.xls' % self.option('dist_method'), 'xls', '样本距离矩阵结果表', 0, ""]
        ])
        if self.option("others_value") != "":
            result_dir_hucluster.add_relpath_rules([
                ["barplot_table.xls", "xls", "柱形图结果表",0,""]
            ])
        super(HclusterWorkflow, self).end()
