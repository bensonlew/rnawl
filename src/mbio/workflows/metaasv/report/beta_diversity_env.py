# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'


import os
import re
import numpy as np
import types
from biocluster.workflow import Workflow
from mbio.packages.metaasv.filter_newick import get_level_newicktree
from mbio.packages.metaasv.common_function import link_dir, link_file


class BetaDiversityEnvWorkflow(Workflow):
    """
    metaasv beta多样性分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(BetaDiversityEnvWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "analysis_type", "type": "string", "default": 'rda_cca'},
            {"name": "dist_method", "type": "string", "default": 'bray_curtis'},#距离
            {"name": "update_info", "type": "string"},
            {"name": "asv_id", "type": "string"},#关联asv表ID
            {"name": "main_id", "type": "string"},#主表ID
            {"name": "level", "type": "int"},
            {"name": "env_file", "type": "infile", "format": "meta.otu.group_table"},##环境因子表
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},##分组表
            {"name": "env_labs", "type": "string", "default": ""},##筛选的环境因子
            {"name": "group_detail", "type": "string", "default": ""},
            {"name": "group_id", "type": "string", "default": ""},
            {"name": "env_id", "type": "string", "default": ""},
            {"name": "good_group", "type": "string", "default": "F"},
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.ellipse = self.add_tool("graph.ellipse")
        self.task = self.add_module("meta.beta_diversity.beta_diversity")
        self.rename = {}
        self.skip_ellipse = True

    def run(self):
        # self.IMPORT_REPORT_DATA = True
        # self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info(self.option('otu_file').path)
        options = {
            'analysis': self.option('analysis_type'),
            # 'dis_method': self.option('dist_method'),
            'otutable': self.option('otu_file'),
            # 'scale': self.option('scale'),
            'ellipse': self.option('good_group')
            }
        if self.option('env_file').is_set:
            options['envlabs'] = self.option('env_labs')
            options['envtable'] = self.option('env_file')

        if self.option('analysis_type') in ['dbrda', 'rda_cca']:
            options['group'] = self.option('group_file')

        if self.option('analysis_type') in ['dbrda']:
            options['dis_method'] = self.option('dist_method')
            if 'unifrac' in self.option('dist_method'):  # sanger_bioinfo/src/mbio/workflows/meta/report/distance_calc.py中的解释
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
                    self.logger.info('get_newick:' + temp_tree_file)
                    otu_table = self.option('otu_file').path
                    temp_otu_file = self.option('otu_file').path + '.temp'
                    all_lines = open(otu_table, 'r').readlines()
                    if len(all_lines) < 3:
                        self.logger.error('分类水平：%s,OTU表数据少于2行：%s' % (self.option('level'), len(all_lines)))  #将Otu改成了OTU modified by hongdongxuan 20170310
                        self.set_error("ASV表数据少于2行")
                    self.logger.info(len(all_lines))
                    new_all = []
                    new_all.append(all_lines[0])
                    for line in all_lines[1:]:
                        name = line.split('\t')
                        origin_name = name[0].split("; ")[-1].strip()
                        if name[0] in all_find:
                            name[0] = 'ASV' + str(all_find[name[0]] + 1)
                        new_all.append('\t'.join(name))
                        if name[0] not in self.rename:
                            self.rename[name[0]] = origin_name
                    otu_file_temp = open(temp_otu_file, 'w')
                    otu_file_temp.writelines(new_all)
                    otu_file_temp.close()
                    os.rename(otu_table, os.path.join(self.work_dir, "otu_table2.xls"))
                    os.rename(temp_otu_file, otu_table)
                    input_table = os.path.join(self.work_dir, "otu_file.xls")
                    options['otutable'] = input_table
                    options['phy_newick'] = temp_tree_file
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
                    for line in all_lines[1:]:  # OTU表中有复杂的名称OTU名称，包含进化物种类型，进化树种只有OTU名称
                        name = line.split('\t')
                        origin_name = name[0].split("; ")[-1].strip()
                        name[0] = name[0].split(';')[-1].strip()
                        new_all.append('\t'.join(name))
                        if name[0] not in self.rename:
                            self.rename[name[0]] = origin_name
                    otu_file_temp = open(temp_otu_file, 'w')
                    otu_file_temp.writelines(new_all)
                    otu_file_temp.close()
                    os.rename(otu_table, os.path.join(self.work_dir, "otu_table2.xls"))
                    os.rename(temp_otu_file, otu_table)
                    input_table = os.path.join(self.work_dir, "otu_file.xls")
                    options['otutable'] = input_table
                    options['phy_newick'] = temp_tree_file
        # 计算置信椭圆>>>
        group_detail = eval(self.option('group_detail'))
        if len(group_detail.keys()) > 1:
            for k in group_detail.keys():
                if len(group_detail[k]) > 2:
                    self.skip_ellipse = False
                    break
        if not self.skip_ellipse:
            self.task.on('end', self.run_ellipse)
            self.ellipse.on('end', self.set_db)
        else:
            self.task.on('end', self.set_db)

        self.task.set_options(options)
        self.task.run()
        super(BetaDiversityEnvWorkflow, self).run()

    def run_ellipse(self):
        """
        计算置信椭圆
        :return:
        """
        options = {}
        if self.option("group_file").is_set:
            options['group_table'] = self.option("group_file")
        options['group_id'] = self.option('group_id')
        pc_map = {'pca': "/Pca/pca_sites.xls", 'pcoa': "/Pcoa/pcoa_sites.xls",
                  'dbrda': '/Dbrda/db_rda_sites.xls', 'nmds': '/Nmds/nmds_sites.xls',
                  'rda_cca': '/Rda', 'plsda': "/Plsda/plsda_sites.xls"
                  }
        options['analysis'] = self.option('analysis_type')
        options['meta'] = self.option('analysis_type')
        options['pc_table'] = self.task.output_dir + pc_map[self.option('analysis_type')]
        self.ellipse.set_options(options)
        # self.ellipse.on('end', self.set_db)
        self.ellipse.run()

    def replace_name(self, input):
        """
        对文件替换为原来的名称
        :return:
        """
        out_table = os.path.join(self.work_dir, "out_table.xls")
        with open(input, "r") as f, open(out_table, "w") as w:
            lines = f.readlines()
            w.write(lines[0])
            for line in lines[1:]:
                line = line.strip().split("\t")
                sp_name = line[0].split("; ")[-1].strip()
                if sp_name in self.rename:
                    line[0] = self.rename[sp_name]
                w.write("\t".join(line) + "\n")
        os.remove(input)
        os.rename(out_table, input)

    def set_db(self):
        """
        保存结果距离矩阵表到mongo数据库中
        """
        api_multi = self.api.api("metaasv.beta_diversity_env")
        if self.option("analysis_type") in ["rda_cca"]:
            link_dir(os.path.join(self.task.output_dir,"Rda"), self.output_dir)
        else:
            link_dir(os.path.join(self.task.output_dir, "Dbrda"), self.output_dir)
            self.replace_name(os.path.join(self.output_dir, "db_rda_plot_species_data.xls"))
            self.replace_name(os.path.join(self.output_dir, "db_rda_species.xls"))
        dir_path = self.output_dir
        cond, cons = [], []
        if self.option('env_file').is_set:
            env_labs = self.option("env_labs")
            cond, cons = self.classify_env(self.option('env_file').path, env_labs=env_labs)
            self.logger.info(cond)
            self.logger.info(cons)
        if not os.path.isdir(dir_path):
            self.logger.info("找不到报告文件夹:{}".format(dir_path))
            self.set_error("找不到报告文件夹")

        api_multi.add_beta_multi_analysis_result(dir_path, self.option('analysis_type'),
                                                 main=False,
                                                 remove=cond,
                                                 main_id=str(self.option('main_id'))
                                                 )
        # 导置信椭圆
        api_common = self.api.api("metaasv.common_api")
        ellipse_file = os.path.join(self.ellipse.work_dir, 'ellipse_out.xls')
        if os.path.exists(ellipse_file):
            if os.path.isfile(ellipse_file):
                api_common.insert_ellipse_table(ellipse_file, str(self.option('main_id')), self.option('analysis_type'))
            else:
                self.logger.info("找不到置信椭圆文件:{}".format(ellipse_file))

        self.logger.info('运行self.end')
        self.end()

    def classify_env(self, env_file, env_labs=None):
        """
        获取环境因子中哪些是条件（条件约束）型因子，哪些是数量（线性约束）型的因子
        """
        if isinstance(env_file, types.StringType) or isinstance(env_file, types.UnicodeType):
            if not os.path.exists(env_file):
                self.logger.error('环境因子文件不存在%s' % env_file)
                self.set_error("环境因子文件不存在")
        else:
            self.set_error('提供的环境因子文件名不是字符串')
        if env_labs:
            cond = []  # 记录不全是数字的因子
            cons = []  # 记录全是数字的因子
            env_name = env_labs.split(',')
            with open(env_file, 'r') as f:
                lines = f.readlines()
                head = lines[0].strip("\r\n").split("\t")
                for line in lines[1:]:
                    line = line.strip("\r\n").split("\t")
                    for env in env_name:
                        if env in head:
                            i = head.index(env)
                            try:
                                float(line[i])
                                if not head[i] in cons:
                                    cons.append(head[i])
                            except ValueError:
                                cond.append(head[i])
                                break
                        else:
                            self.set_error("筛选的环境因子不在环境因子表中")
        else:
            frame = np.loadtxt(env_file, dtype=str, comments='')
            len_env = len(frame[0])
            cond = []  # 记录不全是数字的因子
            cons = []  # 记录全是数字的因子
            for n in xrange(len_env - 1):
                env_values = frame[:, n + 1]
                for value in env_values[1:]:
                    try:
                        float(value)
                    except ValueError:
                        cond.append(env_values[0])
                        break
                else:
                    cons.append(env_values[0])
        return cond, cons  # 前者不全是数字分组， 后者是全部都是数字的分组

    def end(self):
        """
        结束和上传结果文件
        :return:
        """
        if self.option('analysis_type') == 'dbrda':
            file_name = "db-RDA分析结果目录"
            code = ""
        elif self.option('analysis_type') == 'rda_cca':
            file_name = "RDA_CCA分析结果目录"
            code = ""
        else:
            file_name = "Beta_diversity分析结果目录"
            code = ""
        repaths = [
            [".", "", file_name, 0, code],
            ["Distance", "", "距离矩阵计算结果输出目录", 0, ""],
            ["Dbrda", "", "db-RDA分析结果目录", 0, ""],
            ["Dbrda/db_rda_sites.xls", "xls", "db_rda样本坐标表", 0, ""],
            ["Dbrda/db_rda_species.xls", "xls", "db_rda物种坐标表", 0, ""],
            ["Dbrda/db_rda_plot_species_data.xls", "xls", "用于绘图的物种坐标表", 0, ""],
            ["Dbrda/db_rda_centroids.xls", "xls", "db_rda哑变量环境因子坐标表", 0, ""],
            ["Dbrda/db_rda_biplot.xls", "xls", "db_rda数量型环境因子坐标表", 0, ""],
            ["Rda", "", "RDA_CCA分析结果目录", 0, ""],
            [r'Rda/dca.xls', 'xls', 'DCA分析结果', 0, ""],
            ]
        regexps = [
            [r'Distance/%s.*\.xls$' % self.option('dist_method'), 'xls', '样本距离矩阵文件', 0, ""],
            [r'.*/.*_importance\.xls$', 'xls', '主成分解释度表', 0, ""],
            [r'Rda/.*_sites\.xls$', 'xls', '样本坐标表', 0, ""],
            [r'Rda/.*_species\.xls$', 'xls', '物种坐标表', 0, ""],
            [r'Rda/.*_biplot\.xls$', 'xls', '数量型环境因子坐标表', 0, ""],
            [r'Rda/.*_centroids\.xls$', 'xls', '哑变量环境因子坐标表', 0, ""],
            [r'.*/.*_envfit\.xls$', 'xls', 'p_value值与r值表', 0, ""],
            [r'.*plot_species_data\.xls$', 'xls', '用于绘图的物种坐标表', 0, ""]
            ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(BetaDiversityEnvWorkflow, self).end()
