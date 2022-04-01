# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os
import re
from biocluster.workflow import Workflow
from mbio.packages.metaasv.filter_newick import get_level_newicktree
from mbio.packages.metaasv.common_function import link_dir, link_file


class BetaDiversityWorkflow(Workflow):
    """
    metaasv 多样性 群落组成分析 不需要环境因子
    包括hcluster、pca、pcoa、nmds
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(BetaDiversityWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "analysis_type", "type": "string", "default": 'pca'},
            {"name": "dist_method", "type": "string", "default": 'bray_curtis'},
            {"name": "update_info", "type": "string"},
            {"name": "asv_id", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "level", "type": "int"},
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_detail", "type": "string", "default": ""},
            {"name": "group_id", "type": "string", "default": ""},
            {"name": "params", "type": "string", "default": ""},
            {"name": "scale", "type": "string", "default": "F"},  # pca是否进行标准化
            {"name": "good_group", "type": "string", "default": "F"},
            {"name": "diff_test_method", "type": "string", "default": ""},  #pca/pcoa/nmds组间差异检验
            {"name": "change_times", "type": "string", "default": ""},  # pca/pcoa/nmds组间差异检验
            {"name": "hcluster_method", "type": "string", "default": 'average'},# 样本层级聚类方法
            {"name": "others_value", "type": "float", "default": ""},  # 是否进行柱形图绘制
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.task = self.add_module("meta.beta_diversity.beta_diversity")
        self.ellipse = self.add_tool("graph.ellipse")
        self.skip_ellipse = True
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False

    def check_option(self):
        """
        参数二次检查
        :return:
        """
        if not self.option("otu_file").is_set:
            raise self.set_error("输入文件不存在".format(self.option("otu_file").prop['path']))
        if not self.option("group_file").is_set:
            raise self.set_error("输入文件不存在".format(self.option("group_file").prop['path']))

    def run_sort_samples(self):
        """
        是否需要对others合并
        """
        self.logger.info("others:{} ".format(self.option("others_value")))
        self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")  # by houshuang 20190918, 用于计算柱形图数据
        self.sort_samples.set_options(dict({
            "in_otu_table": self.option("otu_file").prop["path"],
            "others": float(self.option("others_value"))
        }))
        self.sort_samples.on("end", self.run_beta_diversity)
        self.sort_samples.run()

    def run_beta_diversity(self):
        """
        运行 beta_diversity 的模块
        """
        self.logger.info(self.option('otu_file').path)
        options = {
            'analysis': self.option('analysis_type'),
            'otutable': self.option('otu_file'),
            }
        if self.option('analysis_type') in ['pca', 'pcoa', 'nmds']:
            options['group'] = self.option('group_file')
            options['grouplab'] = self.option('group_file').prop['group_scheme'][0]
            options['scale'] = self.option('scale')
            options['ellipse'] = self.option('good_group')
        if self.option('analysis_type') in ['hcluster']:
            options['linkage'] = self.option('hcluster_method')
        if self.option('analysis_type') in ['pcoa', 'nmds', 'dbrda']:
            options['dis_method'] = self.option('dist_method')
            if 'unifrac' in self.option('dist_method'):  # sanger_bioinfo/src/mbio/workflows/meta/report/distance_calc.py中的解释
                if self.option('level') != 9:
                    newicktree = get_level_newicktree(self.option('asv_id'), level=self.option('level'),
                                                      tempdir=self.work_dir, return_file=False, bind_obj=self)
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
                        self.set_error("OTU表数据少于2行")
                    self.logger.info(len(all_lines))
                    new_all = []
                    new_all.append(all_lines[0])
                    for line in all_lines[1:]:
                        name = line.split('\t')
                        if name[0] in all_find:
                            name[0] = 'ASV' + str(all_find[name[0]] + 1)
                        new_all.append('\t'.join(name))
                    otu_file_temp = open(temp_otu_file, 'w')
                    otu_file_temp.writelines(new_all)
                    otu_file_temp.close()
                    options['otutable'] = temp_otu_file
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
                        name[0] = name[0].split(';')[-1].strip()
                        new_all.append('\t'.join(name))
                    otu_file_temp = open(temp_otu_file, 'w')
                    otu_file_temp.writelines(new_all)
                    otu_file_temp.close()
                    options['otutable'] = temp_otu_file
                    options['phy_newick'] = temp_tree_file
        # 计算分组椭圆
        group_detail = eval(self.option('group_detail'))
        if len(group_detail.keys()) > 1:
            for k in group_detail.keys():
                if len(group_detail[k]) > 2:
                    self.skip_ellipse = False
                    break
        if self.option('analysis_type') in ['hcluster']:
            self.task.on('end', self.set_db)
        if not self.skip_ellipse:
            self.task.on('end', self.run_ellipse)
            self.ellipse.on('end', self.set_db)
        else:
            self.task.on('end', self.set_db)
        if self.option('diff_test_method') != '' and self.option('change_times') != '':
            options['diff_test_method'] = self.option('diff_test_method')
            options['change_times'] = self.option('change_times')
        self.task.set_options(options)
        self.task.run()

    def run(self):
        """
        运行
        """
        if self.option("others_value") != "" and self.option('analysis_type') in ["hcluster"]:
            self.run_sort_samples()
        else:
            self.run_beta_diversity()
        super(BetaDiversityWorkflow, self).run()

    def run_ellipse(self):
        """
        计算置信椭圆
        """
        options = {}
        if self.option("group_file").is_set:
            options['group_table'] = self.option("group_file")
        options['group_id'] = self.option('group_id')
        pc_map = {'pca': "/Pca/pca_sites.xls", 'pcoa': "/Pcoa/pcoa_sites.xls", 'nmds': '/Nmds/nmds_sites.xls'}
        options['analysis'] = self.option('analysis_type')
        options['meta'] = self.option('analysis_type')
        options['pc_table'] = self.task.output_dir + pc_map[self.option('analysis_type')]
        self.ellipse.set_options(options)
        self.ellipse.run()

    def set_db(self):
        """
        保存结果距离矩阵表到mongo数据库中
        """
        if self.option("analysis_type") in ["pca"]:
            link_dir(os.path.join(self.task.output_dir, "Pca") , self.output_dir)
            distance_file = os.path.join(self.task.output_dir, "Distance", "euclidean_otu_file.xls")
            if os.path.exists(distance_file):
                link_file(distance_file, os.path.join(self.output_dir,"PCA_euclidean.xls"))
        elif self.option("analysis_type") in ["pcoa"]:
            link_dir(os.path.join(self.task.output_dir, "Pcoa"), self.output_dir)
            distance_file = os.path.join(self.task.output_dir, "Distance", self.option("dist_method")+"_otu_file.xls")
            link_file(distance_file, os.path.join(self.output_dir, self.option("dist_method")+".xls"))
        elif self.option("analysis_type") in ["nmds"]:
            link_dir(os.path.join(self.task.output_dir, "Nmds"), self.output_dir)
            distance_file = os.path.join(self.task.output_dir, "Distance", self.option("dist_method")+"_otu_file.xls")
            link_file(distance_file, os.path.join(self.output_dir, self.option("dist_method")+".xls"))
        elif self.option('analysis_type') in ['hcluster']:
            link_dir(os.path.join(self.task.output_dir, "Hcluster"), self.output_dir)
            distance_file = os.path.join(self.task.output_dir, "Distance", self.option("dist_method")+"_otu_file.xls")
            link_file(distance_file, os.path.join(self.output_dir, self.option("analysis_type")+"_"+self.option("dist_method")+".xls"))
            if self.option("others_value") != "" :
                rank_others_file = self.sort_samples.output_dir + "/taxa.percents.table.xls"
                link_file(rank_others_file, os.path.join(self.output_dir, "barplot_table.xls"))
        dir_path = self.output_dir
        if self.option("analysis_type") in ["pca", "pcoa", "nmds"]:
            api_pca = self.api.api("metaasv.beta_diversity")
            api_pca.add_beta_multi_analysis_result(dir_path, self.option('analysis_type'),main=False,main_id=str(self.option('main_id')))
        else:
            api_pca = self.api.api("metaasv.beta_diversity")
            api_pca.add_beta_multi_analysis_result(dir_path, self.option('analysis_type'),main=False,main_id=str(self.option('main_id')), group_file=self.option("group_file").prop['path'])

        api_common = self.api.api("metaasv.common_api")
        ### 导入置信椭圆
        self.logger.info("skip_ellipse:{}".format(self.skip_ellipse))
        ellipse_file = os.path.join(self.ellipse.work_dir, 'ellipse_out.xls')
        self.logger.info("ellipse_file:{}".format(ellipse_file))
        if not self.skip_ellipse:
            if os.path.isfile(ellipse_file):
                api_common.insert_ellipse_table(ellipse_file, str(self.option('main_id')), self.option('analysis_type'))

        if self.option("diff_test_method") in ["adonis"]:
            diff_file = os.path.join(self.task.output_dir, "Adonis", "adonis_results.txt")
            link_file(diff_file, os.path.join(self.output_dir, self.option("analysis_type")+"_adonis.xls" ))
        elif self.option("diff_test_method") in ["anosim"]:
            diff_file = os.path.join(self.task.output_dir, "Anosim", "format_results.xls")
            link_file(diff_file, os.path.join(self.output_dir, self.option("analysis_type")+"_anosim.xls" ))
        else:
            diff_file = ''
        # 组间差异检验
        if self.option("analysis_type") in ["pca", "pcoa", "nmds"]:
            correlation_key = self.option("analysis_type") + "_id"
            coll_name = self.option("analysis_type") + "_table"
        if self.option('diff_test_method') != '' and self.option('change_times') != '':
            self.logger.info("检验方法:{}".format(self.option('diff_test_method')))
            if os.path.isfile(diff_file):
                api_common.insert_anosim_detail(diff_file, str(self.option('main_id')), self.option('diff_test_method'), correlation_key=correlation_key,coll_name=coll_name, main_coll=self.option("analysis_type"))

        self.logger.info('数据导入MongoDB完成！')
        self.end()

    def end(self):
        repaths = [
            [".", "", "", 0, ""],
            ["Distance", "", "距离矩阵计算结果输出目录", 0, ""],
            ["Dbrda", "", "db-RDA分析结果目录", 0, ""],
            ["Dbrda/db_rda_sites.xls", "xls", "db_rda样本坐标表", 0, ""],
            ["Dbrda/db_rda_species.xls", "xls", "db_rda物种坐标表", 0, ""],
            ["Dbrda/db_rda_plot_species_data.xls", "xls", "用于绘图的物种坐标表", 0, ""],
            ["Dbrda/db_rda_centroids.xls", "xls", "db_rda哑变量环境因子坐标表", 0, ""],
            ["Dbrda/db_rda_biplot.xls", "xls", "db_rda数量型环境因子坐标表", 0, ""],
            ["Nmds", "", "NMDS分析结果输出目录", 0, ""],
            ["Nmds/nmds_sites.xls", "xls", "样本坐标表", 0, ""],
            ["Nmds/nmds_stress.xls", "xls", "样本特征拟合度值", 0, ""],
            ["Pca", "", "PCA分析结果输出目录", 0, ""],
            ["Pca/pca_importance.xls", "xls", "主成分解释度表", 0, ""],
            ["Pca/pca_rotation_all.xls", "xls", "全部物种主成分贡献度表", 0, ""],
            ["Pca/pca_rotation.xls", "xls", "物种主成分贡献度表", 0, ""],
            ["Pca/pca_sites.xls", "xls", "样本坐标表", 0, ""],
            ["Pcoa", "", "PCoA分析结果目录", 0, ""],
            ["Pcoa/pcoa_eigenvalues.xls", "xls", "矩阵特征值", 0, ""],
            ["Pcoa/pcoa_eigenvaluespre.xls", "xls", "特征解释度百分比", 0, ""],
            ["Pcoa/pcoa_sites.xls", "xls", "样本坐标表", 0, ""],
            ["Anosim", "", "Anosim组间差异检验结果目录",0,""],
            ["Anosim/anosim_results.txt", "txt", "Anosim检验结果表",0,""],
            ["Anosim/format_results.xls", "xls", "Anosim整理结果表",0,""],
            ["Adonis", "", "Adonis组间差异检验结果目录",0,""],
            ["Adonis/adonis_results.txt", "txt", "Adonis检验结果表",0,""]
            # <<<
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
        super(BetaDiversityWorkflow, self).end()
