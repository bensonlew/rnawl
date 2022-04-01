# -*- coding: utf-8 -*-
# __author__ = 'shenghe'

from biocluster.module import Module
import os
from biocluster.core.exceptions import OptionError
import shutil

class BetaDiversityModule(Module):

    def __init__(self, work_id):
        super(BetaDiversityModule, self).__init__(work_id)
        self.step.add_steps('ChooseAnalysis', 'MultipleAnalysis')
        options = [
            {"name": "analysis", "type": "string",
             "default": "distance,anosim,pca,pcoa,nmds,rda_cca,dbrda,hcluster,plsda"},
            {"name": "dis_method", "type": "string", "default": "bray_curtis"},
            {"name": "dbrda_method", "type": "string", "default": ""},
            # 当设定此值时，dbrda的计算方式将会改变，使用R中自带的距离算法，而不是先计算好距离矩阵，此处的计算方式与一般的距离计算的的值不一致
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "level", "type": "string", "default": "otu"},
            {"name": "phy_newick", "type": "infile", "format": "meta.beta_diversity.newick_tree"},
            {"name": "permutations", "type": "int", "default": 999},
            {"name": "linkage", "type": "string", "default": "average"},
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "envlabs", "type": "string", "default": ""},
            {"name": "pca_envlabs", "type": "string", "default": ""},
            {"name": "dbrda_envlabs", "type": "string", "default": ""},
            {"name": "rda_envlabs", "type": "string", "default": ""},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "grouplab", "type": "string", "default": ""},
            {"name": "anosim_grouplab", "type": "string", "default": ""},
            {"name": "plsda_grouplab", "type": "string", "default": ""},
            {"name": "dis_matrix", "type": "outfile", "format": "meta.beta_diversity.distance_matrix"},
            {"name": "dis_newicktree", "type": "outfile", "format": "meta.beta_diversity.newick_tree"},
            {"name": "scale", "type": "string", "default": "F"},  # pca是否进行标准化 ，add by zouxuan
            {"name": "ellipse", "type": "string", "default": "F"},
            {"name": "diff_test_method", "type": "string", "default": ""},  # by houshuang 20190924 pca/pcoa/nmds组间差异检验
            {"name": "change_times", "type": "string", "default": ""},  # by houshuang 20190924 pca/pcoa/nmds组间差异检验
            {"name": "meta_group_name", "type": "string"},
            {"name": "others_value", "type": "float", "default": ""}
        ]
        self.add_option(options)
        self.matrix = self.add_tool('meta.beta_diversity.distance_calc')
        self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.tools = {}

    def set_otu_table(self):
        """
        根据level返回进行计算的otu表,并设定参数
        :return:
        """
        if self.option('otutable').format == "meta.otu.tax_summary_dir":
            self.otu_table = self.option('otutable').get_table(self.option('level'))
            # self.option('otutable').set_path(otu_table)
            # self.option('otutable').get_info()
            return self.otu_table
        else:
            self.otu_table = self.option('otutable').prop['path']
            return self.otu_table

    def check_options(self):
        for i in ['distance', 'anosim', 'pca', 'pcoa', 'nmds', 'rda_cca', 'dbrda', 'hcluster', 'plsda']:
            if i in self.option('analysis'):
                break
        else:
            raise OptionError('没有选择任何分析或者分析类型选择错误：%s', variables=(self.option('analysis')), code="22700201")
        self.set_otu_table()
        if self.option('permutations') < 0 or self.option('permutations') > 10000:
            raise OptionError('参数permutations：%s 不在范围内(0-10000)', variables=(self.option('permutations')),
                              code="22700202")
        if self.option('linkage') not in ['average', 'single', 'complete']:
            raise OptionError('错误的层级聚类方式：%s', variables=(self.option('linkage')), code="22700203")
        if ('rda_cca' or 'dbrda') in self.option('analysis') and not self.option('envtable').is_set:
            raise OptionError('计算RDA/CCA和dbRDA需要环境因子表', code="22700204")
        if ('anosim' or 'plsda') in self.option('analysis') and not self.option('group').is_set:
            raise OptionError('anosim分析和plsda分析需要相关分组文件', code="22700205")
        if self.option('scale') not in ['T', 'F']:
            raise OptionError('scale只能为T或者F', code="22700206")
        if self.option('ellipse') not in ['T', 'F']:
            raise OptionError('ellipse必须为T或者F', code="22700207")
        return True

    def matrix_run(self):
        """
        运行计算距离矩阵
        :return:
        """
        if self.option('phy_newick').is_set:
            self.matrix.set_options({'method': self.option('dis_method'),
                                     'otutable': self.otu_table,
                                     'newicktree': self.option('phy_newick')})
        else:
            # by houshuang 20190924 组间差异检验，pca距离算法默认为euclidean
            self.matrix.set_options({'method': 'euclidean' if self.option('analysis') == "pca" else self.option('dis_method'),
                                     'otutable': self.otu_table})
        self.matrix.on('end', self.set_output, 'distance')
        self.matrix.run()

    def hcluster_run(self):
        self.tools['hcluster'].set_options({
            'dis_matrix': self.matrix.option('dis_matrix'),
            'linkage': self.option('linkage')
        })
        self.tools['hcluster'].on('end', self.set_output, 'hcluster')
        self.tools['hcluster'].run()

    def anosim_run(self):
        self.tools['anosim'].set_options({
            'dis_matrix': self.matrix.option('dis_matrix'),
            'group': self.option('group'),
            'grouplab': self.option('grouplab') if self.option('grouplab') else self.option('anosim_grouplab'),
            # by houshuang 20190924 增加组间差异检验置换次数和检验方法
            'permutations': self.option('change_times') if self.option('change_times') else self.option('permutations'),
            'diff_test_method': self.option('diff_test_method')
        })
        self.tools['anosim'].on('end', self.set_output, 'anosim')
        self.tools['anosim'].run()

    def anosim_box_run(self):
        """
        add by wangzhaoyue 20170204
        :return:
        """
        self.tools['anosim_box'].set_options({
            'dis_matrix': self.matrix.option('dis_matrix'),
            'grouplab': self.option('grouplab') if self.option('grouplab') else self.option('anosim_grouplab'),
            'group': self.option('group'),
            'permutations': self.option('permutations')
        })
        self.tools['anosim_box'].on('end', self.set_output, 'anosim_box')
        self.tools['anosim_box'].run()

    def box_run(self, rely_obj):
        self.logger.info("box run iiiiiiiiiiiiiiiiiiiiiiing")
        self.tools['box'].set_options({
            'dis_matrix': self.matrix.option('dis_matrix'),
            'grouplab': self.option('grouplab') if self.option('grouplab') else self.option('anosim_grouplab'),
            'group': self.option('group'),
            'permutations': self.option('permutations')
        })
        self.tools['box'].on('end', self.set_output, 'box')
        self.tools['box'].run()

    def pcoa_run(self, rely_obj):
        pcoa_options = {'dis_matrix': self.matrix.option('dis_matrix'), 'ellipse': self.option('ellipse')}
        if self.option('group').is_set:
            pcoa_options['group'] = self.option('group')
        self.tools['pcoa'].set_options(pcoa_options)
        self.tools['pcoa'].on('end', self.set_output, 'pcoa')
        self.tools['pcoa'].run()

    def nmds_run(self, rely_obj):
        nmds_options = {'dis_matrix': self.matrix.option('dis_matrix'), 'ellipse': self.option('ellipse')}
        if self.option('group').is_set:
            nmds_options['group'] = self.option('group')
        self.tools['nmds'].set_options(nmds_options)
        self.tools['nmds'].on('end', self.set_output, 'nmds')
        self.tools['nmds'].run()

    def dbrda_run(self):
        if not self.option('dbrda_method'):
            dbrda_options = {'dis_matrix': self.matrix.option('dis_matrix'), 'envtable': self.option('envtable'),
                             'otutable': self.otu_table, 'method':self.option('dis_method')}  # modify by zhujuan for add db_rda_species.xls 20171011
        else:
            dbrda_options = {'otutable': self.otu_table, 'envtable': self.option('envtable'),
                             'method': self.option('dbrda_method'), 'method':self.option('dis_method')}
        if self.option('envlabs'):
            dbrda_options['envlabs'] = self.option('envlabs')
        else:
            dbrda_options['envlabs'] = self.option('dbrda_envlabs')
        # by houshuang 20190923 分组椭圆 >>>
        if self.option('group').is_set:
            dbrda_options['group_table'] = self.option('group')
        dbrda_options['ellipse'] = self.option('ellipse')
        # <<<
        self.tools['dbrda'].set_options(dbrda_options)
        self.tools['dbrda'].on('end', self.set_output, 'dbrda')
        self.tools['dbrda'].run()

    def rda_run(self):
        rda_options = {'otutable': self.otu_table, 'envtable': self.option('envtable')}
        if self.option('envlabs'):
            rda_options['envlabs'] = self.option('envlabs')
        else:
            rda_options['envlabs'] = self.option('rda_envlabs')
        # by houshuang 20190923 分组椭圆>>>
        if self.option('group').is_set:
            rda_options['group_table'] = self.option('group')
        rda_options['ellipse'] = self.option('ellipse')
        # <<<
        self.tools['rda'].set_options(rda_options)
        self.tools['rda'].on('end', self.set_output, 'rda')
        self.tools['rda'].run()

    def pca_run(self):
        pca_options = {'otutable': self.otu_table,'scale':self.option('scale'), 'ellipse': self.option('ellipse')}
        if self.option('group').is_set:
            pca_options['group_table'] = self.option('group')
        if self.option('envtable').is_set:
            pca_options['envtable'] = self.option('envtable')
            if self.option('envlabs'):
                pca_options['envlabs'] = self.option('envlabs')
            else:
                pca_options['envlabs'] = self.option('pca_envlabs')
        self.tools['pca'].set_options(pca_options)
        self.tools['pca'].on('end', self.set_output, 'pca')
        self.tools['pca'].run()

    def plsda_run(self):
        plsda_options = {'otutable': self.option('otutable'), 'group': self.option('group'),'ellipse': self.option('ellipse')}
        if self.option('grouplab'):
            plsda_options['grouplab'] = self.option('grouplab')
        else:
            plsda_options['grouplab'] = self.option('plsda_grouplab')
        self.tools['plsda'].set_options(plsda_options)
        self.tools['plsda'].on('end', self.set_output, 'plsda')
        self.tools['plsda'].run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'pca':
            self.linkdir(obj.output_dir, 'Pca')
        elif event['data'] == 'rda':
            self.linkdir(obj.output_dir, 'Rda')
        elif event['data'] == 'distance':
            self.linkdir(obj.output_dir, 'Distance')
            self.option('dis_matrix', obj.option('dis_matrix'))
        elif event['data'] == 'hcluster':
            self.linkdir(obj.output_dir, 'Hcluster')
            self.option('dis_newicktree', obj.option('newicktree'))
            if self.option("others_value") != "":
                rank_others_file = self.sort_samples.output_dir + "/taxa.percents.table.xls"
                final_rank_others_file = self.output_dir+"/"+self.option('meta_group_name')+"/Hcluster/barplot_table.xls"
                shutil.copy2(rank_others_file, final_rank_others_file)
        elif event['data'] == 'anosim':
            # pca/pcoa/nmds组间差异检验文件名为Adonis或Anosim
            if self.option('diff_test_method') == 'adonis':
                self.linkdir(obj.output_dir, 'Adonis')
            else:
                self.linkdir(obj.output_dir, 'Anosim')
        elif event['data'] == 'box':
            self.linkdir(obj.output_dir, 'Box')
        elif event['data'] == 'anosim_box':
            self.linkdir(obj.output_dir, 'AnosimBox')  # add by wzy
        elif event['data'] == 'dbrda':
            self.linkdir(obj.output_dir, 'Dbrda')
        elif event['data'] == 'pcoa':
            self.linkdir(obj.output_dir, 'Pcoa')
        elif event['data'] == 'nmds':
            self.linkdir(obj.output_dir, 'Nmds')
        elif event['data'] == 'plsda':
            self.linkdir(obj.output_dir, 'Plsda')
        else:
            pass

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        if self.option('meta_group_name'):
            if not os.path.exists(self.output_dir+"/"+self.option('meta_group_name')):
                os.mkdir(self.output_dir+"/"+self.option('meta_group_name'))
            newdir = os.path.join(self.output_dir, self.option('meta_group_name'), dirname)
        else:
            newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])

    def run(self):
        super(BetaDiversityModule, self).run()
        if self.option("others_value") != "":
            self.sort_samples.set_options({
                "in_otu_table": self.otu_table,
                "others": self.option("others_value")
            })
            self.sort_samples.on("end", self.run_all)
            self.sort_samples.run()
        else:
            self.run_all()

    def run_all(self):
        self.step.ChooseAnalysis.start()
        self.step.update()
        if 'distance' in self.option('analysis'):
            self.tools['distance'] = self.matrix
        if 'anosim' in self.option('analysis'):
            self.tools['anosim'] = self.add_tool('meta.beta_diversity.anosim')
            self.matrix.on('end', self.anosim_run)
            self.tools['anosim_box'] = self.add_tool('meta.beta_diversity.anosim_box')
            self.matrix.on('end', self.anosim_box_run)
            # self.tools['box'] = self.add_tool(
            #     'meta.beta_diversity.distance_box')
            # self.matrix.on('end', self.box_run)  # change by wzy
        if 'pcoa' in self.option('analysis'):
            self.tools['pcoa'] = self.add_tool('meta.beta_diversity.pcoa')
            self.matrix.on('end', self.pcoa_run)
            # by houshuang 20190924 组间差异检验 >>>
            if self.option('diff_test_method') != '' and self.option('change_times') != '':
                self.tools['anosim'] = self.add_tool('meta.beta_diversity.anosim')
                self.matrix.on('end', self.anosim_run)
            # <<<
        if 'nmds' in self.option('analysis'):
            self.tools['nmds'] = self.add_tool('meta.beta_diversity.nmds')
            self.matrix.on('end', self.nmds_run)
            # by houshuang 20190924 组间差异检验 >>>
            if self.option('diff_test_method') != '' and self.option('change_times') != '':
                self.tools['anosim'] = self.add_tool('meta.beta_diversity.anosim')
                self.matrix.on('end', self.anosim_run)
            # <<<
        if 'hcluster' in self.option('analysis'):
            self.tools['hcluster'] = self.add_tool(
                'meta.beta_diversity.hcluster')
            self.matrix.on('end', self.hcluster_run)
        if 'dbrda' in self.option('analysis'):
            self.tools['dbrda'] = self.add_tool('meta.beta_diversity.dbrda')
            if self.option('dbrda_method'):
                pass
            else:
                self.matrix.on('end', self.dbrda_run)
        if 'pca' in self.option('analysis'):
            self.tools['pca'] = self.add_tool('meta.beta_diversity.pca')
            # by houshuang 20190924 组间差异检验 >>>
            if self.option('diff_test_method') != '' and self.option('change_times') != '':
                self.tools['anosim'] = self.add_tool('meta.beta_diversity.anosim')
                self.matrix.on('end', self.anosim_run)
            # <<<
        if 'plsda' in self.option('analysis'):
            self.tools['plsda'] = self.add_tool('meta.beta_diversity.plsda')
        if 'rda_cca' in self.option('analysis'):
            self.tools['rda'] = self.add_tool('meta.beta_diversity.rda_cca')
        self.step.ChooseAnalysis.finish()
        self.step.MultipleAnalysis.start()
        self.step.update()
        print self.tools
        print self.matrix.events["end"]._func
        if self.tools:
            if len(self.tools) == 1:
                self.tools.values()[0].on('end', self.stepend)
            else:
                self.on_rely(self.tools.values(), self.stepend)
            if ('pcoa' in self.option('analysis')) or ('distance' in self.option('analysis')) or ('anosim' in self.option('analysis')) or ('nmds' in self.option('analysis')):
                self.matrix_run()
            if 'dbrda' in self.option('analysis'):
                if self.option('dbrda_method'):
                    self.dbrda_run()
                else:
                    if ('pcoa' in self.option('analysis')) or ('distance' in self.option('analysis')) or ('anosim' in self.option('analysis')) or ('nmds' in self.option('analysis')):
                        pass
                    else:
                        self.matrix_run()
            if 'pca' in self.option('analysis'):
                self.pca_run()
                # by houshuang 20191008 组间差异检验, 兼容大workflow>>>
                if ('pcoa' not in self.option('analysis')) and ('distance' not in self.option('analysis')) and ('anosim' not in self.option('analysis')) and ('nmds' not in self.option('analysis')):
                    if self.option('diff_test_method') != '' and self.option('change_times') != '':
                        self.matrix_run()
                # <<<
            if 'plsda' in self.option('analysis'):
                self.plsda_run()
            if 'rda_cca' in self.option('analysis'):
                self.rda_run()
        else:
            self.matrix_run()


    def stepend(self):
        self.step.MultipleAnalysis.finish()
        self.step.update()
        self.end()

    def end(self):
        repaths = [
            [".", "", "Beta_diversity分析结果文件目录"],
            ["Anosim", "", "anosim&adonis结果输出目录"],
            ["Anosim/anosim_results.txt", "txt", "anosim分析结果"],
            ["Anosim/adonis_results.txt", "txt", "adonis分析结果"],
            ["Anosim/format_results.xls", "xls", "anosim&adonis整理结果表"],
            ["Dbrda", "", "db_rda分析结果目录"],
            ["Dbrda/db_rda_sites.xls", "xls", "db_rda样本坐标表"],
            ["Dbrda/db_rda_species.xls", "xls", "db_rda物种坐标表"],
            ["Dbrda/db_rda_centroids.xls", "xls", "db_rda哑变量环境因子坐标表"],
            ["Dbrda/db_rda_biplot.xls", "xls", "db_rda数量型环境因子坐标表"],
            ["Box", "", "距离统计和统计检验分析结果目录"],
            ["Box/Stats.xls", "xls", "分组统计检验结果"],
            ["Box/Distances.xls", "xls", "组内组间距离值统计结果"],
            ["Distance", "", "距离矩阵计算结果输出目录"],
            ["Hcluster", "", "层次聚类结果目录"],
            ["Hcluster/hcluster.tre", "tre", "层次聚类树"],
            ["Nmds", "", "NMDS分析结果输出目录"],
            ["Nmds/nmds_sites.xls", "xls", "样本坐标表"],
            ["Pca", "", "PCA分析结果输出目录"],
            ["Pca/pca_importance.xls", "xls", "主成分解释度表"],
            ["Pca/pca_rotation.xls", "xls", "物种主成分贡献度表"],
            ["Pca/pca_sites.xls", "xls", "样本坐标表"],
            ["Pca/pca_envfit_factor_scores.xls", "xls", "哑变量环境因子表"],
            ["Pca/pca_envfit_factor.xls", "xls", "哑变量环境因子坐标表"],
            ["Pca/pca_envfit_vector_scores.xls", "xls", "数量型环境因子表"],
            ["Pca/pca_envfit_vector.xls", "xls", "数量型环境因子坐标表"],
            ["Pcoa", "", "pcoa分析结果目录"],
            ["Pcoa/pcoa_eigenvalues.xls", "xls", "矩阵特征值"],
            ["Pcoa/pcoa_sites.xls", "xls", "样本坐标表"],
            ["Plsda", "", "plsda分析结果目录"],
            ["Plsda/plsda_sites.xls", "xls", "样本坐标表"],
            ["Plsda/plsda_rotation.xls", "xls", "物种主成分贡献度表"],
            ["Plsda/plsda_importance.xls", "xls", "主成分解释度表"],
            ["Rda", "", "rda_cca分析结果目录"],
            [r'Rda/dca.xls', 'xls', 'DCA分析结果'],
        ]
        regexps = [
            [r'Distance/%s.*\.xls$' % self.option('dis_method'), 'xls', '样本距离矩阵文件'],
            [r'Rda/.*_importance\.xls$', 'xls', '主成分解释度表'],
            [r'Rda/.*_sites\.xls$', 'xls', '样本坐标表'],
            [r'Rda/.*_species\.xls$', 'xls', '物种坐标表'],
            [r'Rda/.*_biplot\.xls$', 'xls', '数量型环境因子坐标表'],
            [r'Rda/.*_centroids\.xls$', 'xls', '哑变量环境因子坐标表'],
            [r'Rda/.*_envfit\.xls$', 'xls', 'p_value值与r值表']
        ]
        # self.logger.info('shenghe:不能重复添加目录buglog。。。。。。。。。。。。。。。。')
        # self.logger.info(self.upload_dir)
        for i in self.upload_dir:
            self.logger.info(i.path)
        # self.logger.info('shenghe:不能重复添加目录buglog。。。。。。。。。。。。。。。。OVER')
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(BetaDiversityModule, self).end()
