# -*- coding: utf-8 -*-
# __author__ = 'shenghe'

"""距离矩阵层级聚类"""

import datetime
from biocluster.workflow import Workflow
from mbio.packages.beta_diversity.filter_newick import get_level_newicktree
from bson import ObjectId
import re
import os
import json
import shutil
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class HclusterWorkflow(Workflow):
    """
    报告中调用距离矩阵计算样本层级聚类数使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(HclusterWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "dist_method", "type": "string", "default": 'bray_curtis'},
            {"name": "hcluster_method", "type": "string", "default": 'average'},
            {"name": "level", "type": 'int', "default": 9},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "otu_id", "type": "string"},
            {"name": "others_value", "type": "float", "default": ""}  # by houshuang 20190918
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.dist = self.add_tool("meta.beta_diversity.distance_calc")
        self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")  # by houshuang 20190918, 用于计算柱形图数据
        self.hcluster = self.add_tool("meta.beta_diversity.hcluster")

    def run(self):
        # by houshuang 20190918 >>>
        if self.option("others_value") != "":
            self.sort_samples.set_options({
                "in_otu_table": self.option("otu_table"),
                "others": self.option("others_value")
            })
            self.sort_samples.on("end", self.run_dist)
            self.sort_samples.run()
        else:
            self.run_dist()
        # <<<
        super(HclusterWorkflow, self).run()

    def run_dist(self):
        if 'unifrac' in self.option('dist_method'):
            # 查找OTU表对应的进化树
            if self.option('level') != 9:
                newicktree = get_level_newicktree(self.option('otu_id'), level=self.option('level'),
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
                    return 'OTU' + str(match_newname.count)  # 后面替换OTU中名称用同样的命名规则
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
                    self.set_error("otu表数据少于2行", code="12701801")
                new_all = []
                new_all.append(all_lines[0])
                for line in all_lines[1:]:  # 遍历OTU表，将OTU表的复杂OTU名称改为之前find到的复杂名称对应的字典
                    name = line.split('\t')
                    if name[0] not in all_find:
                        self.set_error('OTU表中有原始表不存在的OTU名：%s', variables=(name[0]), code="12701802")
                    name[0] = 'OTU' + str(all_find[name[0]] + 1)
                    new_all.append('\t'.join(name))
                otu_file_temp = open(temp_otu_file, 'w')
                otu_file_temp.writelines(new_all)
                otu_file_temp.close()
                options = {
                    'method': self.option('dist_method'),
                    'otutable': temp_otu_file,
                    'newicktree': temp_tree_file
                }
            else:
                newicktree = get_level_newicktree(self.option('otu_id'), level=self.option('level'),
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
                options = {
                    'method': self.option('dist_method'),
                    'otutable': temp_otu_file,
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
        params_json = json.loads(self.option('params'))
        api_distance = self.api.distance
        matrix_path = self.dist.output_dir + '/' + os.listdir(self.dist.output_dir)[0]
        final_matrix_path = os.path.join(self.output_dir, os.listdir(self.dist.output_dir)[0])
        shutil.copy2(matrix_path, final_matrix_path)
        if not os.path.isfile(matrix_path):
            self.logger.error("找不到报告文件:{}".format(matrix_path))
            self.set_error("找不到报告文件", code="12701803")
        dist_name = 'Distance_{}_{}'.format(self.option('dist_method'),
                                            datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        matrix_id = api_distance.add_dist_table(matrix_path,
                                                major=True,
                                                name=dist_name,
                                                level=self.option('level'),
                                                otu_id=self.option('otu_id'),
                                                params=params_json)
        # self.add_return_mongo_id('sg_beta_specimen_distance', matrix_id)
        api_newick = self.api.newicktree
        collection = api_newick.db["sg_beta_specimen_distance"]
        newick_fath = self.hcluster.output_dir + "/hcluster.tre"
        final_newick_path = os.path.join(self.output_dir, "hcluster.tre")
        shutil.copy2(newick_fath, final_newick_path)
        if not os.path.isfile(newick_fath):
            self.logger.error("找不到报告文件:{}".format(newick_fath))
            self.set_error("找不到报告文件", code="12701803")
        return_id = api_newick.add_tree_file(newick_fath, major=False, tree_id=self.option('main_id'),
                                             update_dist_id=matrix_id)
        # by houshuang 20190918 >>>
        if self.option("others_value") != "":
            # 链接tool结果文件
            # rank_file = self.sort_samples.output_dir + "/taxa.table.xls"
            rank_others_file = self.sort_samples.output_dir + "/taxa.percents.table.xls"
            # if not os.path.isfile(rank_file):
            #     self.logger.error("找不到报告文件:{}".format(rank_file))
            #     self.set_error("找不到报告文件")
            if not os.path.isfile(rank_others_file):
                self.logger.error("找不到报告文件:{}".format(rank_others_file))
                self.set_error("找不到报告文件", code="12701804")
            # final_rank_file = os.path.join(self.output_dir, "taxa.table.xls")
            # shutil.copy2(rank_file, final_rank_file)
            final_rank_others_file = os.path.join(self.output_dir, "barplot_table.xls")
            shutil.copy2(rank_others_file, final_rank_others_file)
            # 更新物种名至主表
            api_newick.update_newick(final_rank_others_file, return_id)
            # 创建detail表
            api_newick.add_newick_detail(final_rank_others_file, return_id)
        # <<<
        collection.update_one({"_id": ObjectId(matrix_id)}, {'$set': {'newick_tree_id': ObjectId(return_id)}})
        #self.end()
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"),"sg_newick_tree")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "beta_sample_distance",
                "interaction": 1,
                "main_table": "sg_newick_tree",
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        save_params(self.output_dir, self.id)
        result_dir_hucluster = self.add_upload_dir(self.output_dir)
        result_dir_hucluster.add_relpath_rules([
            [".", "", "样本层级聚类分析结果目录", 0, "110049"],
            ["./hcluster.tre", "tre", "层级聚类树结果表", 0, "110050"],  # modified by hongdongxuan 20170321
            ["./样本层级聚类分析图.pdf", "pdf", "样本层级聚类分析图", 0, ""]
            ])
        result_dir_hucluster.add_regexp_rules([
            [r'%s.*\.xls' % self.option('dist_method'), 'xls', '样本距离矩阵结果表', 0, "110048"]
        ])
        # by houshuang 20190918 >>>
        if self.option("others_value") != "":
            result_dir_hucluster.add_relpath_rules([
                ["barplot_table.xls", "xls", "柱形图结果表",0,"110250"]
            ])
        # <<<
        super(HclusterWorkflow, self).end()
