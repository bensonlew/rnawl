# -*- coding: utf-8 -*-
# __author__ = 'shenghe'

"""beta多元分析"""

import os
import re
import numpy as np
import types
from biocluster.workflow import Workflow
from mbio.packages.beta_diversity.filter_newick import get_level_newicktree
from mbio.packages.meta.common_function import envname_restore
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class BetaMultiAnalysisWorkflow(Workflow):
    """
    报告中调用beta多元分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(BetaMultiAnalysisWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "analysis_type", "type": "string", "default": 'pca'},
            {"name": "dist_method", "type": "string", "default": 'bray_curtis'},
            {"name": "update_info", "type": "string"},
            {"name": "otu_id", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "level", "type": "int"},
            {"name": "multi_analysis_id", "type": "string"},
            {"name": "env_file", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "env_labs", "type": "string", "default": ""},
            {"name": "group_detail", "type": "string", "default": ""},
            {"name": "group_id", "type": "string", "default": ""},
            {"name": "env_id", "type": "string", "default": ""},
            {"name": "params", "type": "string", "default": ""},
            # {"name": "matrix_out", "type": "outfile", "format": "meta.beta_diversity.distance_matrix"}
            {"name": "scale", "type": "string", "default": "F"},  # pca是否进行标准化 ，add by zouxuan
            {"name": "good_group", "type": "string", "default": "F"},
            {"name": "diff_test_method", "type": "string", "default": ""},  # by houshuang 20190924 pca/pcoa/nmds组间差异检验
            {"name": "change_times", "type": "string", "default": ""}  # by houshuang 20190924 pca/pcoa/nmds组间差异检验
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.task = self.add_module("meta.beta_diversity.beta_diversity")  # task改为self.task
        # by houshuang 20190921 置信椭圆 >>>
        self.ellipse = self.add_tool("graph.ellipse")
        self.skip_ellipse = True
        self.rename = {}
        # <<<

    def run(self):
        self.logger.info(self.option('otu_file').path)
        options = {
            'analysis': self.option('analysis_type'),
            # 'dis_method': self.option('dist_method'),
            'otutable': self.option('otu_file'),
            'scale': self.option('scale'),
            'ellipse': self.option('good_group')
            }
        if self.option('env_file').is_set:
            options['envlabs'] = self.option('env_labs')
            options['envtable'] = self.option('env_file')
        else:
            pass
        #if self.option('analysis_type') == 'plsda':
        if self.option('analysis_type') in ['plsda', 'pca', 'pcoa', 'nmds']:
            options['group'] = self.option('group_file')
            options['grouplab'] = self.option('group_file').prop['group_scheme'][0]
        # by houshuang 20190924 计算分组椭圆需要>>>
        if self.option('analysis_type') in ['dbrda', 'rda_cca']:  # 增加dbrda, rda_cca
            options['group'] = self.option('group_file')
        # <<<
        if self.option('analysis_type') in ['pcoa', 'nmds', 'dbrda']:
            options['dis_method'] = self.option('dist_method')
            if 'unifrac' in self.option('dist_method'):  # sanger_bioinfo/src/mbio/workflows/meta/report/distance_calc.py中的解释
                if self.option('level') != 9:
                    newicktree = get_level_newicktree(self.option('otu_id'), level=self.option('level'),
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
                        return 'OTU' + str(match_newname.count)
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
                        self.set_error("OTU表数据少于2行", code="12700601")
                    self.logger.info(len(all_lines))
                    new_all = []
                    new_all.append(all_lines[0])
                    for line in all_lines[1:]:
                        name = line.split('\t')
                        origin_name = name[0].split("; ")[-1].strip()
                        if name[0] in all_find:
                            name[0] = 'OTU' + str(all_find[name[0]] + 1)
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
        # by houshuang 20190921 计算置信椭圆>>>
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
        # pca/pcoa/nmds组间差异检验
        if self.option('diff_test_method') != '' and self.option('change_times') != '':
            options['diff_test_method'] = self.option('diff_test_method')
            options['change_times'] = self.option('change_times')
        # <<<
        self.task.set_options(options)
        self.task.run()
        self.output_dir = self.task.output_dir
        super(BetaMultiAnalysisWorkflow, self).run()

    # by houshuang 20190921, 计算置信椭圆
    def run_ellipse(self):
        options = {}
        if self.option("group_file").is_set:
            options['group_table'] = self.option("group_file")
        options['group_id'] = self.option('group_id')
        pc_map = {'pca': "/Pca/pca_sites.xls", 'pcoa': "/Pcoa/pcoa_sites.xls",
                  'dbrda': '/Dbrda/db_rda_sites.xls', 'nmds': '/Nmds/nmds_sites.xls',
                  'rda_cca': '/Rda', 'plsda': "/Plsda/plsda_sites.xls"
                  ##'rda_cca'
                  }
        options['analysis'] = self.option('analysis_type')
        options['meta'] = self.option('analysis_type')
        options['pc_table'] = self.task.output_dir + pc_map[self.option('analysis_type')]
        self.ellipse.set_options(options)
        self.ellipse.run()

    def replace_name(self, input):
        """
        对将修改的物种名称替换为原来的名称
        qingchen.zhang @20200611
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
        api_multi = self.api.beta_multi_analysis
        if self.option('analysis_type') in ['dbrda']:
            self.replace_name(os.path.join(self.output_dir, "Dbrda" ,"db_rda_plot_species_data.xls"))
            self.replace_name(os.path.join(self.output_dir, "Dbrda" ,"db_rda_species.xls"))
        dir_path = self.output_dir
        cond, cons = [], []
        if self.option('env_file').is_set:
            env_labs = self.option("env_labs")
            cond, cons = self.classify_env(self.option('env_file').path, env_labs=env_labs)
            self.logger.info(cond)
            self.logger.info(cons)
        if not os.path.isdir(dir_path):
            self.logger.info("找不到报告文件夹:{}".format(dir_path))
            self.set_error("找不到报告文件夹", code="12700602")

        api_multi.add_beta_multi_analysis_result(dir_path, self.option('analysis_type'),
                                                 main=False,
                                                 remove=cond,
                                                 main_id=str(self.option('main_id'))
                                                 )
        # by houshuang 20190921 置信椭圆 >>>
        if not self.skip_ellipse:
            ellipse_file = os.path.join(self.ellipse.work_dir, 'ellipse_out.xls')
            if os.path.isfile(ellipse_file):
                api_multi.insert_ellipse_table(ellipse_file, str(self.option('main_id')), self.option('analysis_type'))
            else:
                self.logger.info("找不到分组椭圆文件:{}".format(ellipse_file))
        # 组间差异检验
        if self.option('diff_test_method') != '' and self.option('change_times') != '':
            self.logger.info("检验方法:{}".format(self.option('diff_test_method')))
            if self.option('diff_test_method') == 'anosim':
                out_file = self.task.output_dir.rstrip('/') + '/Anosim/format_results.xls'
            else:
                out_file = self.task.output_dir.rstrip('/') + '/Adonis/adonis_results.txt'
            if os.path.isfile(out_file):
                api_multi.insert_anosim_detail(out_file, str(self.option('main_id')), self.option('diff_test_method'))
            else:
                self.logger.info("找不到组间差异检验结果文件:{}".format(out_file))
        # <<<
        self.logger.info('运行self.end')
        #self.end()
        self.save_pdf()

    def save_pdf(self):
        self.logger.info("get_save_pdf_status：{}".format(get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))))
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status and self.option('analysis_type') in ["pca","pcoa","nmds","anosim","plsda","rda_cca","dbrda"]:
            name = get_name(self.option("main_id"),"sg_beta_multi_analysis")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            submit_loc = "beta_pca"
            if self.option('analysis_type') == "pca":
                submit_loc = "beta_pca"
            elif self.option('analysis_type') == "pcoa":
                submit_loc = "beta_pcoa"
            elif self.option('analysis_type') == "nmds":
                submit_loc = "beta_nmds"
            elif self.option('analysis_type') == "anosim":
                submit_loc = "beta_anosim"
            elif self.option('analysis_type') == "plsda":
                submit_loc = "beta_plsda"
            elif self.option('analysis_type') == "rda_cca":
                submit_loc = "beta_rda_cca"
            elif self.option('analysis_type') == "dbrda":
                submit_loc = "beta_db_rda"
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": submit_loc,
                "interaction": 1,
                "main_table": "sg_beta_multi_analysis",
            })
            self.figsave.run()
        else:
            self.end()

    def classify_env(self, env_file, env_labs=None):
        """
        获取环境因子中哪些是条件（条件约束）型因子，哪些是数量（线性约束）型的因子
        """
        if isinstance(env_file, types.StringType) or isinstance(env_file, types.UnicodeType):
            if not os.path.exists(env_file):
                self.logger.error('环境因子文件不存在%s' % env_file)
                self.set_error("环境因子文件不存在", code="12700603")
        else:
            self.set_error('提供的环境因子文件名不是字符串', code="12700604")
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
                            self.set_error("筛选的环境因子不在环境因子表中", code="12700605")
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

    @envname_restore
    def end(self):
        save_params(self.output_dir, self.id)
        if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            if self.option('analysis_type') == 'plsda':
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir+"/Plsda/"))
            elif self.option('analysis_type') == 'pca':
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir+"/Pca/"))
            elif self.option('analysis_type') == 'pcoa':
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir+"/Pcoa/"))
            elif self.option('analysis_type') == 'nmds':
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir+"/Nmds/"))
            elif self.option('analysis_type') == 'dbrda':
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir+"/Dbrda/"))
            elif self.option('analysis_type') == 'rda_cca':
                os.system("cp -r {}/* {}/".format(pdf_outs+"/RDA/", self.output_dir+"/Rda/"))

        if self.option('analysis_type') == 'plsda':   #add 14 lines by hongdongxuan 20170327
            file_name = "PLS_DA分析结果目录"
            code = "110092"
        elif self.option('analysis_type') == 'pca':
            file_name = "PCA分析结果目录"
            code = "110051"
        elif self.option('analysis_type') == 'pcoa':
            file_name = "PCoA分析结果目录"
            code = "110063"
        elif self.option('analysis_type') == 'nmds':
            file_name = "NMDS分析结果目录"
            code = "110060"
        elif self.option('analysis_type') == 'dbrda':
            file_name = "db-RDA分析结果目录"
            code = "110106"
        elif self.option('analysis_type') == 'rda_cca':
            file_name = "RDA_CCA分析结果目录"
            code = "110097"
        else:
            file_name = "Beta_diversity分析结果目录"
            code = "110091"
        repaths = [
            [".", "", file_name, 0, code],
            ["Distance", "", "距离矩阵计算结果输出目录", 0, "110089"],
            ["Dbrda", "", "db-RDA分析结果目录", 0, "110106"],
            ["Dbrda/db_rda_sites.xls", "xls", "db_rda样本坐标表", 0, "110109"],
            ["Dbrda/db_rda_species.xls", "xls", "db_rda物种坐标表", 0, "110110"],
            ["Dbrda/db_rda_plot_species_data.xls", "xls", "用于绘图的物种坐标表", 0, "110107"],
            ["Dbrda/db_rda_centroids.xls", "xls", "db_rda哑变量环境因子坐标表", 0, "110111"],
            ["Dbrda/db_rda_biplot.xls", "xls", "db_rda数量型环境因子坐标表", 0, "110108"],
            ["Nmds", "", "NMDS分析结果输出目录", 0, "110060"],
            ["Nmds/nmds_sites.xls", "xls", "样本坐标表", 0, "110061"],
            ["Nmds/nmds_stress.xls", "xls", "样本特征拟合度值", 0, "110062"],
            ["Pca", "", "PCA分析结果输出目录", 0, "110051"],
            ["Pca/pca_importance.xls", "xls", "主成分解释度表", 0, "110088"],
            ["Pca/pca_rotation_all.xls", "xls", "全部物种主成分贡献度表", 0, "110081"],
            ["Pca/pca_rotation.xls", "xls", "物种主成分贡献度表", 0, "110082"],
            ["Pca/pca_sites.xls", "xls", "样本坐标表", 0, "110083"],
            ["Pca/pca_envfit_factor_scores.xls", "xls", "哑变量环境因子坐标表", 0, "110085"],
            ["Pca/pca_envfit_factor.xls", "xls", "哑变量环境因子表", 0, "110084"],
            ["Pca/pca_envfit_vector_scores.xls", "xls", "数量型环境因子坐标表", 0, "110087"],
            ["Pca/pca_envfit_vector.xls", "xls", "数量型环境因子表", 0, "110086"],
            ["Pcoa", "", "PCoA分析结果目录", 0, "110063"],
            ["Pcoa/pcoa_eigenvalues.xls", "xls", "矩阵特征值", 0, "110066"],
            ["Pcoa/pcoa_eigenvaluespre.xls", "xls", "特征解释度百分比", 0, "110064"],
            ["Pcoa/pcoa_sites.xls", "xls", "样本坐标表", 0, "110065"],
            ["Plsda", "", "PLS_DA分析结果目录", 0, "110092"],
            ["Plsda/plsda_sites.xls", "xls", "样本坐标表", 0, "110096"],
            ["Plsda/plsda_rotation.xls", "xls", "物种主成分贡献度表", 0, "110095"],
            ["Plsda/plsda_importance.xls", "xls", "主成分组别特征值表", 0, "110093"],
            ["Plsda/plsda_importancepre.xls", "xls", "主成分解释度表", 0, "110094"],
            ["Rda", "", "RDA_CCA分析结果目录", 0, "110097"],
            [r'Rda/dca.xls', 'xls', 'DCA分析结果', 0, "110098"],
            # by houshuang 20190925 组间差异检验结果>>>
            ["Anosim", "", "Anosim组间差异检验结果目录",0,"110241"],
            ["Anosim/anosim_results.txt", "txt", "Anosim检验结果表",0,"110242"],
            ["Anosim/format_results.xls", "xls", "Anosim整理结果表",0,"110243"],
            ["Adonis", "", "Adonis组间差异检验结果目录",0,"110244"],
            ["Adonis/adonis_results.txt", "txt", "Adonis检验结果表",0,"110245"],
            ["Pca/PCA分析散点图.pdf", "pdf", "样本PCA分析散点图", 0, ""],
            ["Pca/PCA分析箱式图.pdf", "pdf", "PCA分析箱式图", 0, ""],
            ["Pcoa/PCoA分析散点图.pdf", "pdf", "样本PCoA分析散点图", 0, ""],
            ["Pcoa/PCoA分析箱式图.pdf", "pdf", "PCoA分析箱式图", 0, ""],
            ["Nmds/NMDS分析散点图.pdf", "pdf", "样本NMDS分析散点图", 0, ""],
            ["Anosim/ANOSIM分析箱式图.pdf", "pdf", "ANOSIM分析箱式图", 0, ""],
            ["Plsda/PLS_DA图.pdf", "pdf", "PLS_DA分析散点图", 0, ""],
            ["Rda/RDA_CCA分析结果图.pdf", "pdf", "RDA/CCA分析结果图", 0, ""],
            ["Dbrda/db_RDA分析结果图.pdf", "pdf", "db_RDA分析结果图", 0, ""],
            # <<<
            ]
        regexps = [
            [r'Distance/%s.*\.xls$' % self.option('dist_method'), 'xls', '样本距离矩阵文件', 0, "110090"],
            [r'.*/.*_importance\.xls$', 'xls', '主成分解释度表', 0, "110104"],
            [r'Rda/.*_sites\.xls$', 'xls', '样本坐标表', 0, "110100"],
            [r'Rda/.*_species\.xls$', 'xls', '物种坐标表', 0, "110103"],
            [r'Rda/.*_biplot\.xls$', 'xls', '数量型环境因子坐标表', 0, "110102"],
            [r'Rda/.*_centroids\.xls$', 'xls', '哑变量环境因子坐标表', 0, "110105"],
            [r'.*/.*_envfit\.xls$', 'xls', 'p_value值与r值表', 0, "110099"],
            [r'.*plot_species_data\.xls$', 'xls', '用于绘图的物种坐标表', 0, "110101"]
            ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(BetaMultiAnalysisWorkflow, self).end()
