# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.module import Module
import os
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir, link_file


class BetaDiversityModule(Module):
    """
    metaasv多样性分析
    """
    def __init__(self, work_id):
        super(BetaDiversityModule, self).__init__(work_id)
        options = [
            {"name": "analysis", "type": "string","default": "distance,pca,pcoa,nmds,hcluster"},
            {"name": "dis_method", "type": "string", "default": "bray_curtis"},
            {"name": "otutable", "type": "infile", "format": "metaasv.otu_table, metaasv.tax_summary_dir"},
            {"name": "level", "type": "string", "default": "asv"},
            {"name": "phy_newick", "type": "infile", "format": "meta.beta_diversity.newick_tree"},
            {"name": "permutations", "type": "int", "default": 999},
            {"name": "linkage", "type": "string", "default": "average"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "grouplab", "type": "string", "default": ""},
            {"name": "anosim_grouplab", "type": "string", "default": ""},
            {"name": "plsda_grouplab", "type": "string", "default": ""},
            {"name": "dis_matrix", "type": "outfile", "format": "meta.beta_diversity.distance_matrix"},
            {"name": "dis_newicktree", "type": "outfile", "format": "meta.beta_diversity.newick_tree"},
            {"name": "scale", "type": "string", "default": "F"},  # pca是否进行标准化 ，add by zouxuan
            {"name": "ellipse", "type": "string", "default": "F"},
            {"name": "diff_test_method", "type": "string", "default": ""},  # by houshuang 20190924 pca/pcoa/nmds组间差异检验
            {"name": "change_times", "type": "string", "default": ""}  # by houshuang 20190924 pca/pcoa/nmds组间差异检验
        ]
        self.add_option(options)
        self.matrix = self.add_tool('meta.beta_diversity.distance_calc')
        self.analysis = self.option("analysis").split(",")
        self.tools = {}

    def set_otu_table(self):
        """
        根据level返回进行计算的otu表,并设定参数
        :return:
        """
        if self.option('otutable').format == "meta.otu.tax_summary_dir":
            self.otu_table = self.option('otutable').get_table(self.option('level'))
            return self.otu_table
        else:
            self.otu_table = self.option('otutable').prop['path']
            return self.otu_table

    def check_options(self):
        for i in ['pca', 'pcoa', 'nmds', 'hcluster']:
            if i in self.analysis:
                break
        else:
            raise OptionError('没有选择任何分析或者分析类型选择错误：%s')
        self.set_otu_table()
        if self.option('permutations') < 0 or self.option('permutations') > 10000:
            raise OptionError('参数permutations：%s 不在范围内(0-10000)')
        if self.option('linkage') not in ['average', 'single', 'complete']:
            raise OptionError('错误的层级聚类方式：%s')
        if self.option('scale') not in ['T', 'F']:
            raise OptionError('scale只能为T或者F')
        if self.option('ellipse') not in ['T', 'F']:
            raise OptionError('ellipse必须为T或者F')
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
            # 组间差异检验，pca距离算法默认为euclidean
            self.matrix.set_options({'method': 'euclidean' if self.option('analysis') == "pca" else self.option('dis_method'),
                                     'otutable': self.otu_table})
        self.matrix.on('end', self.set_output, 'distance')
        self.matrix.run()

    def hcluster_run(self):
        """
        层级聚类
        :return:
        """
        self.tools['hcluster'].set_options({
            'dis_matrix': self.matrix.option('dis_matrix'),
            'linkage': self.option('linkage')
        })
        self.tools['hcluster'].on('end', self.set_output, 'hcluster')
        self.tools['hcluster'].run()

    def pcoa_run(self):
        """
        PCoA分析
        :return:
        """
        pcoa_options = {'dis_matrix': self.matrix.option('dis_matrix'),
                        'ellipse': self.option('ellipse')}
        if self.option('group').is_set:
            pcoa_options['group'] = self.option('group')
        self.tools['pcoa'].set_options(pcoa_options)
        self.tools['pcoa'].on('end', self.set_output, 'pcoa')
        self.tools['pcoa'].run()

    def nmds_run(self):
        """
        NMDS分析
        :return:
        """
        nmds_options = {'dis_matrix': self.matrix.option('dis_matrix'),
                        'ellipse': self.option('ellipse')}
        if self.option('group').is_set:
            nmds_options['group'] = self.option('group')
        self.tools['nmds'].set_options(nmds_options)
        self.tools['nmds'].on('end', self.set_output, 'nmds')
        self.tools['nmds'].run()


    def pca_run(self):
        """
        PCA分析
        :return:
        """
        pca_options = {'otutable': self.otu_table,
                       'scale':self.option('scale'),
                       'ellipse': self.option('ellipse')}
        if self.option('group').is_set:
            pca_options['group_table'] = self.option('group')
        self.tools['pca'].set_options(pca_options)
        self.tools['pca'].on('end', self.set_output, 'pca')
        self.tools['pca'].run()

    def anosim_run(self):
        """
        Anosim和Adonis组间差异分析
        :return:
        """
        self.tools['anosim'].set_options({
            'dis_matrix': self.matrix.option('dis_matrix'),
            'group': self.option('group'),
            'grouplab': self.option('grouplab') if self.option('grouplab') else self.option('anosim_grouplab'),
            'permutations': self.option('change_times') if self.option('change_times') else self.option('permutations'),
            'diff_test_method': self.option('diff_test_method')
        })
        self.tools['anosim'].on('end', self.set_output, 'anosim')
        self.tools['anosim'].run()

    def set_output(self, event):
        """
        设置结果文件目录
        :param event:
        :return:
        """
        obj = event['bind_object']
        if event['data'] == 'pca':
            link_dir(obj.output_dir, self.output_dir + '/Pca')
        elif event['data'] == 'distance':
            link_dir(obj.output_dir, self.output_dir + '/Distance')
            self.option('dis_matrix', obj.option('dis_matrix'))
        elif event['data'] == 'hcluster':
            link_dir(obj.output_dir, self.output_dir + '/Hcluster')
            self.option('dis_newicktree', obj.option('newicktree'))
        elif event['data'] == 'pcoa':
            link_dir(obj.output_dir, self.output_dir +'/Pcoa')
        elif event['data'] == 'nmds':
            link_dir(obj.output_dir, self.output_dir + '/Nmds')
        elif event['data'] == 'anosim':
            # pca/pcoa/nmds组间差异检验文件名为Adonis或Anosim
            if self.option('diff_test_method') == 'adonis':
                link_dir(obj.output_dir, 'Adonis')
            else:
                link_dir(obj.output_dir, 'Anosim')
        else:
            pass

    def run(self):
        """
        运行
        :return:
        """
        super(BetaDiversityModule, self).run()
        # self.step.ChooseAnalysis.start()
        # self.step.update()
        if 'distance' in self.analysis:
            self.tools['distance'] = self.matrix
        if 'pcoa' in self.analysis:
            self.tools['pcoa'] = self.add_tool('meta.beta_diversity.pcoa')
            self.matrix.on('end', self.pcoa_run)
            # 组间差异检验
            if self.option('diff_test_method') != '' and self.option('change_times') != '':
                self.tools['anosim'] = self.add_tool('meta.beta_diversity.anosim')
                self.matrix.on('end', self.anosim_run)
        if 'nmds' in self.analysis:
            self.tools['nmds'] = self.add_tool('meta.beta_diversity.nmds')
            self.matrix.on('end', self.nmds_run)
            #  组间差异检验
            if self.option('diff_test_method') != '' and self.option('change_times') != '':
                self.tools['anosim'] = self.add_tool('meta.beta_diversity.anosim')
                self.matrix.on('end', self.anosim_run)
        if 'hcluster' in self.analysis:
            self.tools['hcluster'] = self.add_tool('meta.beta_diversity.hcluster')
            self.matrix.on('end', self.hcluster_run)
        if 'pca' in self.analysis:
            self.tools['pca'] = self.add_tool('meta.beta_diversity.pca')
            # 组间差异检验
            if self.option('diff_test_method') != '' and self.option('change_times') != '':
                self.tools['anosim'] = self.add_tool('meta.beta_diversity.anosim')
                self.matrix.on('end', self.anosim_run)

        # self.step.ChooseAnalysis.finish()
        # self.step.MultipleAnalysis.start()
        # self.step.update()
        if self.tools:
            if len(self.tools) == 1:
                self.tools.values()[0].on('end', self.end)
            else:
                self.on_rely(self.tools.values(), self.end)
            if ('pcoa' in self.analysis) or ('distance' in self.analysis) or ('anosim' in self.analysis) or ('nmds' in self.analysis):
                self.matrix_run()
            if 'pca' in self.analysis:
                self.pca_run()
                # 组间差异检验, 兼容大workflow
                if ('pcoa' not in self.analysis) and ('distance' not in self.analysis) and ('anosim' not in self.analysis) and ('nmds' not in self.analysis):
                    if self.option('diff_test_method') != '' and self.option('change_times') != '':
                        self.matrix_run()
        else:
            self.matrix_run()

    def stepend(self):
        self.step.MultipleAnalysis.finish()
        self.step.update()
        self.end()

    def end(self):
        """
        结束
        :return:
        """
        super(BetaDiversityModule, self).end()
