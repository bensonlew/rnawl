# -*- coding: utf-8 -*-
# __author__ = 'shenghe'

"""otu表的样本距离计算"""

import os
import re
from biocluster.workflow import Workflow
from bson import ObjectId
from mbio.packages.beta_diversity.filter_newick import *
import datetime
from mbio.packages.meta.save_params import save_params


class DistanceCalcWorkflow(Workflow):
    """
    报告中调用otu计算样本距离时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(DistanceCalcWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "method", "type": "string", "default": 'bray_curtis'},
            {"name": "update_info", "type": "string"},
            {"name": "otu_id", "type": "string"},
            {"name": "level", "type": "int"},
            {"name": "matrix_id", "type": "string"},
            {"name": "task_type", "type": "string"},
            # {"name": "matrix_out", "type": "outfile", "format": "meta.beta_diversity.distance_matrix"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run(self):
        self.task = self.add_tool("meta.beta_diversity.distance_calc")
        self.logger.info(self.option('otu_file').path)
        if 'unifrac' in self.option('method'):
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
                otu_table = self.option('otu_file').path
                temp_otu_file = self.option('otu_file').path + '.temp'
                all_lines = open(otu_table, 'r').readlines()
                if len(all_lines) < 3:
                    self.logger.error('分类水平：%s,otu表数据少于2行：%s' % (self.option('level'), len(all_lines)))
                    self.set_error("otu表数据少于2行", code="12701001")
                new_all = []
                new_all.append(all_lines[0])
                for line in all_lines[1:]:  # 遍历OTU表，将OTU表的复杂OTU名称改为之前find到的复杂名称对应的字典
                    name = line.split('\t')
                    if name[0] not in all_find:
                        self.set_error('OTU表中存在原表没有的OTU: %s', variables=(name[0]), code="12701002")
                        raise Exception('OTU表中存在不是直接通过组合原始表分类名称的OTU名：%s' % name[0])
                    name[0] = 'OTU' + str(all_find[name[0]] + 1)
                    new_all.append('\t'.join(name))
                otu_file_temp = open(temp_otu_file, 'w')
                otu_file_temp.writelines(new_all)
                otu_file_temp.close()
                options = {
                    'method': self._sheet.option('method'),
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
                otu_table = self.option('otu_file').path
                temp_otu_file = self.option('otu_file').path + '.temp'
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
                    'method': self._sheet.option('method'),
                    'otutable': temp_otu_file,
                    'newicktree': temp_tree_file
                }
        else:
            options = {
                'method': self.option('method'),
                'otutable': self.option('otu_file')
            }
        self.task.set_options(options)
        self.task.on('end', self.set_db)
        self.task.run()
        self.output_dir = self.task.output_dir
        super(DistanceCalcWorkflow, self).run()

    def get_phylo_tree(self):
        tree_path = ''
        return tree_path

    def set_db(self):
        """
        保存结果距离矩阵表到mongo数据库中
        """
        api_distance = self.api.distance
        matrix_path = self.output_dir + '/' + os.listdir(self.output_dir)[0]
        if not os.path.isfile(matrix_path):
            self.logger.error("找不到报告文件:{}".format(matrix_path))
            self.set_error("找不到报告文件", code="12701003")
        params_json = {
            'otu_id': self.option('otu_id'),
            'level_id': self.option('level'),
            'distance_algorithm': self.option('method'),
            'task_type': self.option('task_type'),
            'submit_location': 'beta_sample_distance'
            }
        matrix_id = api_distance.add_dist_table(matrix_path,
                                                major=True,
                                                name='Distance_{}_{}'.format(self.option('method'),
                                                                             datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
                                                level=self.option('level'),
                                                otu_id=self.option('otu_id'),
                                                params=params_json)
        self.add_return_mongo_id('sg_beta_specimen_distance', matrix_id)
        self.logger.info(str(matrix_id))
        self.logger.info('运行self.end')
        self.end()

    def end(self):
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "距离矩阵计算结果输出目录", 0, "110089"],
            ["./%s" % os.listdir(self.output_dir)[0], "xls", "样本距离矩阵文件", 0, "110090"],
            ])
        print self.get_upload_files()
        super(DistanceCalcWorkflow, self).end()
