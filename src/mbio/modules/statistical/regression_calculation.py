# -*- coding: utf-8 -*-
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import os
import glob
import types
from mbio.files.meta.otu.otu_table import OtuTableFile


class RegressionCalculationModule(Module):
    """
    module for Regression
    author: gaohao
    """
    MATRIX = ['abund_jaccard', 'binary_chisq', 'binary_chord',
                       'binary_euclidean', 'binary_hamming', 'binary_jaccard',
                       'binary_lennon', 'binary_ochiai',
                       'binary_pearson', 'binary_sorensen_dice',
                       'bray_curtis', 'bray_curtis_faith', 'bray_curtis_magurran',
                       'canberra', 'chisq', 'chord', 'euclidean', 'gower',
                       'hellinger', 'kulczynski', 'manhattan', 'morisita_horn',
                       'pearson', 'soergel', 'spearman_approx', 'specprof',
                       'unifrac', 'unweighted_unifrac', 'unweighted_unifrac_full_tree',
                       'weighted_normalized_unifrac', 'weighted_unifrac']

    def __init__(self, work_id):
        super(RegressionCalculationModule, self).__init__(work_id)
        options = [
            {"name": "diversity_type", "type": "string"},
            {"name": "diversity_analysis_type", "type": "string"},
            {"name": "distance_type", "type": "string"},
            {"name": "tax_abund_file", "type": "infile","format": "meta.otu.otu_table"},
            {"name": "func_abund_file", "type": "infile","format": "meta.otu.otu_table"},
            {"name": "output_tax", "type": "outfile","format": "meta.otu.otu_table"},
            {"name": "output_func", "type": "outfile", "format": "meta.otu.otu_table"},
        ]
        self.tax_beta_diversity = self.add_module('meta.beta_diversity.beta_diversity')
        self.func_beta_diversity = self.add_module('meta.beta_diversity.beta_diversity')
        self.tax_alpha_diversity = self.add_tool('meta.alpha_diversity.estimators')
        self.func_alpha_diversity = self.add_tool('meta.alpha_diversity.estimators')
        self.add_option(options)
        self.beta_list = []


    def tax_diversity_analysis_run(self):
        if self.option('diversity_analysis_type') in ['pca']:
            self.tax_beta_diversity.set_options({
                'otutable': self.option('tax_abund_file'),
                'analysis': self.option('diversity_analysis_type'),
            })
        elif self.option('diversity_analysis_type') in ['pcoa', 'nmds']:
            if self.option('distance_type') in RegressionCalculationModule.MATRIX:
                self.tax_beta_diversity.set_options({
                    'otutable': self.option('tax_abund_file'),
                    'analysis': self.option('diversity_analysis_type'),
                    'dis_method': self.option('distance_type'),
                })
            else:
                raise OptionError('请选择正确的距离方法！', code="24100601")
        else:
            raise OptionError('请选择正确的beta分析类型！', code="24100602")
        self.tax_beta_diversity.on('end',self.set_output, 'beta_tax')
        self.tax_beta_diversity.run()

    def func_diversity_analysis_run(self):
        if self.option('diversity_analysis_type') in ['pca']:
            self.func_beta_diversity.set_options({
                'otutable': self.option('func_abund_file').prop['path'],
                'analysis': self.option('diversity_analysis_type'),
            })
        elif self.option('diversity_analysis_type') in ['pcoa', 'nmds']:
            if self.option('distance_type') in RegressionCalculationModule.MATRIX:
                self.func_beta_diversity.set_options({
                    'otutable': self.option('func_abund_file'),
                    'analysis': self.option('diversity_analysis_type'),
                    'dis_method': self.option('distance_type')
                })
            else:
                raise OptionError('请选择正确的距离方法！', code="24100601")
        else:
            raise OptionError('请选择正确的beta分析类型！', code="24100602")
        self.func_beta_diversity.on('end',self.set_output, 'beta_func')
        #self.beta_list.append(self.beta_diversity)
       # self.on_rely()
        self.func_beta_diversity.run()


    def tax_estimators_run(self):
        if self.option('diversity_analysis_type') in ['shannon','simpson','invsimpson']:
            options = {
               'otu_table':self.option('tax_abund_file'),
               'indices': self.option('diversity_analysis_type'),
               'level': 'otu',
            }
        else:
            raise OptionError('请选择正确的指数！', code="24100603")
        self.tax_alpha_diversity.set_options(options)
        self.tax_alpha_diversity.on("end", self.set_output, 'alpha_tax')
        self.tax_alpha_diversity.run()

    def func_estimators_run(self):
        if self.option('diversity_analysis_type') in ['shannon', 'simpson', 'invsimpson']:
            options = {
                'otu_table':self.option('func_abund_file'),
                'indices': self.option('diversity_analysis_type'),
                'level': 'otu',
            }
        else:
            raise OptionError('请选择正确的指数！', code="24100603")
        self.func_alpha_diversity.set_options(options)
        self.func_alpha_diversity.on("end", self.set_output,'alpha_func')
        self.func_alpha_diversity.run()



    def run(self):
        super(RegressionCalculationModule, self).run()
        if self.option('diversity_type') in ['alpha']:
            if self.option('tax_abund_file').is_set:
                self.tax_estimators_run()
            if self.option('func_abund_file').is_set:
                self.func_estimators_run()
            self.on_rely([self.func_alpha_diversity, self.tax_alpha_diversity], self.end)
        if self.option('diversity_type') in ['beta']:
            if self.option('tax_abund_file').is_set:
               self.tax_diversity_analysis_run()
            if self.option('func_abund_file').is_set:
               self.func_diversity_analysis_run()
            self.on_rely([self.func_beta_diversity,self.tax_beta_diversity], self.end)

    def end(self):
        super(RegressionCalculationModule, self).end()


    def set_output(self, event):
        obj = event['bind_object']
        if self.option('diversity_type') in ['beta']:
            if event['data'] == 'beta_tax':
                self.linkdir(obj.output_dir, 'beta_tax')
                path = self.output_dir + '/beta_tax'
                if self.option('diversity_analysis_type') in ['pca']:
                    #self.logger.info(self.output_dir)
                    self.option('output_tax').set_path(os.path.join(path, "pca_sites.xls"))
                    #self.option('output_tax', os.path.join(self.output_dir,'/beta_tax/', ""))
                if self.option('diversity_analysis_type') in ['pcoa']:
                    self.option('output_tax').set_path(os.path.join(path, "pcoa_sites.xls"))
                if self.option('diversity_analysis_type') in ['nmds']:
                    self.option('output_tax').set_path(os.path.join(path, "nmds_sites.xls"))
            elif event['data'] == 'beta_func':
                self.linkdir(obj.output_dir, 'beta_func')
                path = self.output_dir + '/beta_func'
                if self.option('diversity_analysis_type') in ['pca']:
                    self.option('output_func').set_path(os.path.join(path , "pca_sites.xls"))
                if self.option('diversity_analysis_type') in ['pcoa']:
                    self.option('output_func').set_path(os.path.join(path , "pcoa_sites.xls"))
                if self.option('diversity_analysis_type') in ['nmds']:
                    self.option('output_func').set_path(os.path.join(path , "nmds_sites.xls"))
        if self.option('diversity_type') in ['alpha']:
            if event['data'] == 'alpha_tax':
                self.linkdir(obj.output_dir, 'alpha_tax')
                path=self.output_dir+'/alpha_tax'
                self.option('output_tax').set_path(os.path.join(path, "estimators.xls"))
            elif event['data'] == 'alpha_func':
                self.linkdir(obj.output_dir, 'alpha_func')
                path=self.output_dir+'/alpha_func'
                self.option('output_func').set_path(os.path.join(path, "estimators.xls"))
    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        if self.option('diversity_analysis_type') in ['pca']:
            dir = dirpath +'/Pca'
            allfiles = os.listdir(dir)
            #raise OptionError(dirpath+'/Pca')
            newdir = os.path.join(self.output_dir, dirname)
            if not os.path.exists(newdir):
                os.mkdir(newdir)
            if os.path.exists(newdir + '/pca_sites.xls'):
                os.remove(newdir + '/pca_sites.xls')
            for file in allfiles:
                if file in 'pca_sites.xls':
                    newfile =newdir + '/pca_sites.xls'
                    file=dir+'/pca_sites.xls'
                    os.link(file, newfile)
        if self.option('diversity_analysis_type') in ['pcoa']:
            dir = dirpath +'/Pcoa'
            allfiles = os.listdir(dir)
            newdir = os.path.join(self.output_dir, dirname)
            if not os.path.exists(newdir):
                os.mkdir(newdir)
            if os.path.exists(newdir + '/pcoa_sites.xls'):
                os.remove(newdir + '/pcoa_sites.xls')
            for file in allfiles:
                if file in 'pcoa_sites.xls':
                    file = dir + '/pcoa_sites.xls'
                    newfile =newdir + '/pcoa_sites.xls'
                    os.link(file, newfile)
        if self.option('diversity_analysis_type') in ['nmds']:
            dir = dirpath + '/Nmds'
            allfiles = os.listdir(dir)
            newdir = os.path.join(self.output_dir, dirname)
            if not os.path.exists(newdir):
                os.mkdir(newdir)
            if os.path.exists(newdir+'/nmds_sites.xls'):
                os.remove(newdir+'/nmds_sites.xls')
            for file in allfiles:
                if file in 'nmds_sites.xls':
                    file = dir + '/nmds_sites.xls'
                    newfile = newdir + '/nmds_sites.xls'
                    os.link(file, newfile)

        if self.option('diversity_analysis_type') in ['shannon', 'simpson', 'invsimpson']:
            dir = dirpath
            # raise OptionError(dir)
            allfiles = os.listdir(dir)
            newdir = os.path.join(self.output_dir, dirname)
            if not os.path.exists(newdir):
               os.mkdir(newdir)
            if os.path.exists(newdir + '/estimators.xls'):
               os.remove(newdir + '/estimators.xls')
            for file in allfiles:
               if file in 'estimators.xls':
                  newfile = newdir + '/estimators.xls'
                  file = dir + '/estimators.xls'
                  os.link(file, newfile)
