# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.module import Module
import os
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir, link_file


class BetaDiversityEnvModule(Module):
    """
    metaasv 多样性asv流程
    """
    def __init__(self, work_id):
        super(BetaDiversityEnvModule, self).__init__(work_id)
        self.step.add_steps('ChooseAnalysis', 'MultipleAnalysis')
        options = [
            {"name": "analysis", "type": "string", "default": "rda_cca,dbrda"},
            {"name": "dis_method", "type": "string", "default": "bray_curtis"},
            {"name": "otutable", "type": "infile", "format": "metaasv.otu_table, metaasv.tax_summary_dir"},
            {"name": "level", "type": "string", "default": "asv"},
            {"name": "phy_newick", "type": "infile", "format": "meta.beta_diversity.newick_tree"},
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "envlabs", "type": "string", "default": ""},
            {"name": "dbrda_envlabs", "type": "string", "default": ""},
            {"name": "rda_envlabs", "type": "string", "default": ""},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "dis_matrix", "type": "outfile", "format": "meta.beta_diversity.distance_matrix"},
            {"name": "dis_newicktree", "type": "outfile", "format": "meta.beta_diversity.newick_tree"},
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
        """
        参数二次检查
        :return:
        """
        for i in ['rda_cca', 'dbrda']:
            if i in self.analysis:
                break
        else:
            raise OptionError('没有选择任何分析或者分析类型选择错误：%s')
        self.set_otu_table()
        if ('rda_cca' or 'dbrda') in self.analysis and not self.option('envtable').is_set:
            raise OptionError('计算RDA/CCA和dbRDA需要环境因子表')
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
            self.matrix.set_options({'method': self.option('dis_method'),
                                     'otutable': self.otu_table})
        self.matrix.on('end', self.set_output, 'distance')
        self.matrix.run()

    def dbrda_run(self):
        """
        db-RDA 分析
        :return:
        """
        dbrda_options = {'dis_matrix': self.matrix.option('dis_matrix'),
                         'envtable': self.option('envtable'),
                         'otutable': self.otu_table,
                         'method':self.option('dis_method')}
        if self.option('envlabs'):
            dbrda_options['envlabs'] = self.option('envlabs')
        else:
            dbrda_options['envlabs'] = self.option('dbrda_envlabs')
        self.tools['dbrda'].set_options(dbrda_options)
        self.tools['dbrda'].on('end', self.set_output, 'dbrda')
        self.tools['dbrda'].run()

    def rda_run(self):
        """
        RDA_CCA 分析
        :return:
        """
        rda_options = {'otutable': self.otu_table, 'envtable': self.option('envtable')}
        if self.option('envlabs'):
            rda_options['envlabs'] = self.option('envlabs')
        else:
            rda_options['envlabs'] = self.option('rda_envlabs')
        self.tools['rda'].set_options(rda_options)
        self.tools['rda'].on('end', self.set_output, 'rda')
        self.tools['rda'].run()


    def set_output(self, event):
        """
        设置结果文件目录
        :param event:
        :return:
        """
        obj = event['bind_object']
        if event['data'] == 'rda':
            link_dir(obj.output_dir, self.output_dir + '/Rda')
        elif event['data'] == 'distance':
            link_dir(obj.output_dir, self.output_dir + '/Distance')
            self.option('dis_matrix', obj.option('dis_matrix'))
        elif event['data'] == 'dbrda':
            link_dir(obj.output_dir, self.output_dir + '/Dbrda')
        else:
            pass

    def run(self):
        """
        运行
        :return:
        """
        super(BetaDiversityEnvModule, self).run()
        self.step.ChooseAnalysis.start()
        self.step.update()
        if 'dbrda' in self.analysis:
            self.tools['dbrda'] = self.add_tool('meta.beta_diversity.dbrda')
            self.matrix.on('end', self.dbrda_run)
        if 'rda_cca' in self.analysis:
            self.tools['rda'] = self.add_tool('meta.beta_diversity.rda_cca')
        self.step.ChooseAnalysis.finish()
        self.step.MultipleAnalysis.start()
        self.step.update()
        if self.tools:
            if len(self.tools) == 1:
                self.tools.values()[0].on('end', self.stepend)
            else:
                self.on_rely(self.tools.values(), self.stepend)
            if 'dbrda' in self.analysis:
                self.matrix_run()
            if 'rda_cca' in self.analysis:
                self.rda_run()
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
        super(BetaDiversityEnvModule, self).end()
