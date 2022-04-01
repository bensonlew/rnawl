# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'

"""beta多样性分析"""

import os
import re
import numpy as np
import types
from biocluster.workflow import Workflow
from biocluster.config import Config
from biocluster.core.exceptions import OptionError
from comm_table import CommTableWorkflow
from mbio.packages.meta.common_function import envname_restore
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name
import json

class BetaDiversityWorkflow(CommTableWorkflow):
    """
    交互分析beta多样性
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BetaDiversityWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "analysis_type", "type": "string", "default": 'pca'},
            {"name": "anno_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "profile_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "distance_method", "type": "string", "default": "bray_curtis"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "group_id", "type": "string", "default": ""},
            {"name": "env_file", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_detail", "type": "string", "default": ""},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "env_labs", "type": "string", "default": ""},
            {"name": "group_id", "type": "string", "default": ""},
            {"name": "env_id", "type": "string", "default": ""},
            # {'name': 'save_pdf', 'type': 'int', 'default': 1},  # 是否保存pdf图片导结果目录，add by fangyifei
            {"name": "scale", "type": "string", "default": "F"}  # pca是否进行标准化 ，add by zouxuan

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        #self.abundance = self.add_tool("meta.create_abund_table")
        self.beta = self.add_module("meta.beta_diversity.beta_diversity")
        self.sam = self.add_tool("meta.otu.sort_samples_mg")
        self.ellipse = self.add_tool("graph.ellipse")
        self.skip_ellipse = False

    def run(self):
        # self.IMPORT_REPORT_DATA = True
        # self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info('self._sheet: {}'.format(self._sheet))
        self.logger.info('self._sheet.db_type: {}'.format(self._sheet.db_type))
        if self.option("group_table").is_set:
            if self.option("profile_table").is_set:
                self.sam.on('end', self.run_beta)
            else:
                self.run_abundance(self.sort_sample)
                #self.abundance.on('end', self.sort_sample)
                self.sam.on('end', self.run_beta)
        else:
            if self.option("profile_table").is_set:
                pass
            else:
                self.run_abundance(self.run_beta)
                #self.abundance.on('end', self.run_beta)
        group_detail = eval(self.option('group_detail'))        #zouguanqing
        self.skip_ellipse = True
        if len(group_detail.keys()) >1:
            for k in group_detail.keys():  #guanqing
                if len(group_detail[k]) > 2:
                    #self.big_group_num += 1
                    #if self.big_group_num > 1:
                    self.skip_ellipse = False
                    break
        if not self.skip_ellipse:
            self.beta.on('end', self.run_ellipse)
            self.ellipse.on('end', self.set_db)
        else:
            self.beta.on('end', self.set_db)
        if self.option("profile_table").is_set:
            if self.option("group_table").is_set:
                self.sort_sample()
            else:
                self.run_beta()
        else:
            #self.run_abundance()
            self.abundance.run()
        super(BetaDiversityWorkflow, self).run()

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

    def run_beta(self):
        if self.option("group_table").is_set:
            otutable = self.sam.option("out_otu_table")
        else:
            if self.option("profile_table").is_set:
                otutable = self.option("profile_table")
            else:
                otutable = self.abundance.option('out_table')
        otutable.get_info()
        if otutable.prop['sample_num'] < 2:
            raise OptionError('样品数必须大于等于2', code="12801101")
        if otutable.prop['otu_num'] < 2:
            raise OptionError('物种/功能/基因数必须大于等于2', code="12801102")
        options = {
            'analysis': self.option('analysis_type'),
            # 'dis_method': self.option('dist_method'),
            'otutable': otutable,
            'scale': self.option('scale')
        }
        # options['scale'] = self.option('scale')
        # raise Exception(options['scale'])
        if self.option('env_file').is_set:
            options['envlabs'] = self.option('env_labs')
            options['envtable'] = self.option('env_file')
        else:
            pass
        if self.option('analysis_type') in ['pcoa', 'nmds', 'dbrda']:
            options['dis_method'] = self.option('distance_method')
        # by zhigang.zhao 20200907 分组椭圆需要 >>>
        if self.option('analysis_type') in ['pca', 'pcoa', 'nmds']:
            options['group'] = self.option('group_table')
            if not self.skip_ellipse:
                options['ellipse'] = "T"
        #  <<<
        self.beta.set_options(options)
        #self.beta.on('end', self.set_db)
        self.beta.run()
        self.output_dir = self.beta.output_dir

    def run_ellipse(self):
        options = {}
        if self.option("group_table").is_set:
            options['group_table'] = self.option("group_table")
        options['group_id'] = self.option('group_id')
        pc_map = {'pca':"/Pca/pca_sites.xls",'pcoa': "/Pcoa/pcoa_sites.xls",
                  'dbrda':'/Dbrda/db_rda_sites.xls','nmds':'/Nmds/nmds_sites.xls',
                  'rda_cca': '/Rda'
                  ##'rda_cca'
                  }
        options['analysis'] = self.option('analysis_type')
        options['pc_table'] = self.beta.output_dir + pc_map[self.option('analysis_type')]
        self.ellipse.set_options(options)
        self.ellipse.run()

    def set_db(self):
        """
        保存结果距离矩阵表到mongo数据库中
        """
        api_beta = self.api.api('metagenomic.beta_diversity')
        dir_path = self.output_dir
        # sanger_type, sanger_path = self._sheet.output.split(':')
        # sanger_prefix = Config().get_netdata_config(sanger_type)
        # web_path = os.path.join(sanger_prefix[sanger_type + '_path'], sanger_path)
        # web_path = self._sheet.output.split(":")[-1].lstrip("/")
        web_path = self.output_dir # zzg _sheet.output
        cond, cons = [], []
        if self.option('env_file').is_set:
            cond, cons = self.classify_env(self.option('env_file').path)
            self.logger.info(cond)
            self.logger.info(cons)
        if not os.path.isdir(dir_path):
            self.set_error("找不到报告文件夹:%s", variables=(dir_path), code="12801101")
        file = self.convert(self.option('analysis_type'))
        api_beta.add_beta_diversity(dir_path + '/' + file, self.option('analysis_type'), web_path=os.path.join(web_path, file),
                                    main=None, remove=cond, main_id=str(self.option('main_id')))
        if not self.skip_ellipse:
            api_beta.insert_ellipse_table(self.ellipse.work_dir + '/ellipse_out.xls', str(self.option('main_id')))
        self.logger.info('运行self.end')
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "beta_diversity")
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

    def classify_env(self, env_file):
        """
        获取环境因子中哪些是条件（条件约束）型因子，哪些是数量（线性约束）型的因子
        """
        if isinstance(env_file, types.StringType) or isinstance(env_file, types.UnicodeType):
            if not os.path.exists(env_file):
                raise OptionError('环境因子文件不存在', code="12801103")
        else:
            raise OptionError('提供的环境因子文件名不是一个字符串', code="12801104")
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

    def convert(self, analysis_type):
        if analysis_type == 'plsda':
            file = "PLS_DA"
        if analysis_type == 'pca':
            file = "Pca"
        if analysis_type == 'pcoa':
            file = "Pcoa"
        if analysis_type == 'nmds':
            file = "Nmds"
        if analysis_type == 'dbrda':
            file = 'Dbrda'
        if analysis_type == 'rda_cca':
            file = 'Rda'
        return file

    @envname_restore
    def end(self):
        if self.option('analysis_type') == 'plsda':  # add 14 lines by hongdongxuan 20170327
            file_name = "PLS_DA分析结果目录"
            file_code = "120168"
        elif self.option('analysis_type') == 'pca':
            file_name = "PCA分析结果目录"
            file_code = "120108"
        elif self.option('analysis_type') == 'pcoa':
            file_name = "PCoA分析结果目录"
            file_code = "120109"
        elif self.option('analysis_type') == 'nmds':
            file_name = "NMDS分析结果目录"
            file_code = "120110"
        elif self.option('analysis_type') == 'dbrda':
            file_name = "db-RDA分析结果目录"
            file_code = "120126"
        elif self.option('analysis_type') == 'rda_cca':
            file_name = "RDA_CCA分析结果目录"
            file_code = "120125"
        else:
            file_name = "Beta_diversity分析结果目录"
            file_code = "120168"
        repaths = [
            [".", "", file_name, 0, file_code],
            ["Distance", "", "距离矩阵计算结果输出目录", 0, "120155"],
            ["Dbrda", "", "db-RDA分析结果目录", 0, "120126"],
            ["Dbrda/db_rda_sites.xls", "xls", "样本坐标表", 0, "120162"],
            ["Dbrda/db_rda_envfit.xls", "xls", "p值与r值表", 0, "120159"],
            ["Dbrda/db_rda_species.xls", "xls", "物种_功能坐标表", 0, "120163"],
            # ["Dbrda/db_rda_centroids.xls", "xls", "哑变量环境因子坐标表"],
            ["Dbrda/db_rda_biplot.xls", "xls", "数量型环境因子坐标表", 0, "120158"],
            ["Dbrda/db_rda_plot_species_data.xls", "xls", "绘图物种_功能坐标表", 0, "120161"],
            ["Dbrda/db_rda_importance.xls", "xls", "主成分解释度表", 0, "120160"],
            ["Nmds", "", "NMDS分析结果输出目录", 0, "120110"],
            ["Nmds/nmds_sites.xls", "xls", "样本各维度坐标", 0, "120111"],
            ["Nmds/nmds_stress.xls", "xls", "样本特征拟合度值", 0, "120122"],
            ["Pca", "", "PCA分析结果输出目录", 0, "120108"],
            ["Pca/pca_importance.xls", "xls", "主成分解释度表", 0, "120113"],
            ["Pca/pca_rotation.xls", "xls", "主成分贡献度表", 0, "120115"],
            ["Pca/pca_rotation_all.xls", "xls", "全部主成分贡献度表", 0, "120116"],
            ["Pca/pca_sites.xls", "xls", "样本各成分轴坐标", 0, "120114"],
            ["Pca/pca_envfit_factor_scores.xls", "xls", "哑变量环境因子坐标表", 0, "120117"],
            ["Pca/pca_envfit_factor.xls", "xls", "哑变量环境因子表", 0, "120164"],
            ["Pca/pca_envfit_vector_scores.xls", "xls", "数量型环境因子坐标表", 0, "120119"],
            ["Pca/pca_envfit_vector.xls", "xls", "数量型环境因子表", 0, "120165"],
            ["Pcoa", "", "PCoA分析结果目录", 0, "120109"],
            ["Pcoa/pcoa_eigenvalues.xls", "xls", "矩阵特征值", 0, "120121"],
            ["Pcoa/pcoa_eigenvaluespre.xls", "xls", "特征解释度百分比",0,"120280"],
            ["Pcoa/pcoa_sites.xls", "xls", "样本坐标表", 0, "120123"],
            # ["Plsda", "", "PLS_DA分析结果目录"],
            # ["Plsda/plsda_sites.xls", "xls", "样本坐标表"],
            # ["Plsda/plsda_rotation.xls", "xls", "物种主成分贡献度表"],
            # ["Plsda/plsda_importance.xls", "xls", "主成分组别特征值表"],
            # ["Plsda/plsda_importancepre.xls", "xls", "主成分解释度表"],
            ["Rda", "", "RDA_CCA分析结果目录", 0, "120125"],
            [r'Rda/dca.xls', 'xls', 'DCA分析结果', 0, "120157"],
            ["Pca/PCA.pdf", "pdf", "PCA图"],
            ["Pca/PCA_box.pdf", "pdf", "PCA图箱线图"],
            ["Pcoa/PCoA.pdf", "pdf", "PCoA图"],
            ["Pcoa/PCoA_box.pdf", "pdf", "PCoA图箱线图"],
            ["Nmds/NMDS.pdf", "pdf", "NMDS图"],
            ['Rda/rda_cca.pdf', 'pdf', 'RDA_CCA图'],
            ["Dbrda/db_rda.pdf", "pdf", "db_rda图"],

        ]
        regexps = [
            [r'Distance/%s.*\.xls$' % self.option('distance_method'), 'xls', '样本距离矩阵文件', 0, "120156"],
            [r'Rda/.*_importance\.xls$', 'xls', '主成分解释度表', 0, "120160"],
            [r'Rda/.*_sites\.xls$', 'xls', '样本坐标表', 0, "120162"],
            [r'Rda/.*_species\.xls$', 'xls', '物种_功能坐标表', 0, "120163"],
            [r'Rda/.*_biplot\.xls$', 'xls', '数量型环境因子坐标表', 0, "120158"],
            # [r'Rda/.*_centroids\.xls$', 'xls', '哑变量环境因子坐标表'],
            [r'Rda/.*_envfit\.xls$', 'xls', 'p值与r值表', 0, "120159"],
            [r'Rda/.*_plot_species_data\.xls$', 'xls', '绘图物种_功能坐标表', 0, "120161"]
        ]
        if self.option('analysis_type') == 'pca':
            temp_dir = self.output_dir + "/Pca"
        elif self.option('analysis_type') == 'pcoa':
            temp_dir = self.output_dir + "/Pcoa"
        elif self.option('analysis_type') == 'nmds':
            temp_dir = self.output_dir + "/Nmds"
        elif self.option('analysis_type') == 'dbrda':
            temp_dir = self.output_dir + "/Dbrda"
        else:
            temp_dir = self.output_dir + "/Rda"
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, temp_dir))
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(BetaDiversityWorkflow, self).end()
