# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
from collections import defaultdict
from itertools import chain


class KeggRegulateAgent(Agent):
    """
    Kegg调控统计分析
    version v1.0.1
    author: qiuping
    last_modify: 2016.11.23
    """
    def __init__(self, parent):
        super(KeggRegulateAgent, self).__init__(parent)
        options = [
            {"name": "kegg_table", "type": "infile", "format": "annotation.kegg.kegg_table"},  # 只含有基因的kegg table结果文件
            {"name": "diff_stat", "type": "infile", "format": "rna.diff_stat_table"}
        ]
        self.add_option(options)
        self.step.add_steps("kegg_regulate")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.kegg_regulate.start()
        self.step.update()

    def stepfinish(self):
        self.step.kegg_regulate.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('kegg_table').is_set:
            raise OptionError('必须设置kegg的pathway输入文件')
        if not self.option("diff_stat").is_set:
            raise OptionError("必须设置输入文件diff_stat")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        self._memory = '2G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["pathways", "dir", "kegg调控统计pathway结果图片"],
            ["kegg_regulate_stat.xls", "xls", "kegg调控统计表"],
        ])
        super(KeggRegulateAgent, self).end()


class KeggRegulateTool(Tool):
    def __init__(self, config):
        super(KeggRegulateTool, self).__init__(config)
        self._version = "v1.0.1"
        self.python = '/miniconda2/bin/'
        self.diff_list = self.option('diff_stat').prop['diff_genes']

    def run(self):
        """
        运行
        :return:
        """
        super(KeggRegulateTool, self).run()
        self.run_regulation()

    def run_regulation(self):
        kegg = self.load_package('denovo_rna.express.kegg_regulate')
        ko_genes, path_ko = self.option('kegg_table').get_pathway_koid()
        regulate_gene = self.option('diff_stat').prop['regulate_dict']
        regulate_dict = defaultdict(set)
        path_kos = set(chain(*path_ko.values()))
        for ko in path_kos:
            genes = ko_genes[ko]
            for gene in genes:
                if gene in regulate_gene['up']:
                    regulate_dict['up'].add(ko)
                elif gene in regulate_gene['down']:
                    regulate_dict['down'].add(ko)
                else:
                    # self.logger.info('.....gene:%s' % gene)
                    self.logger.info('......no regulate ko:' + ko)
        up_downs = regulate_dict['up'] & regulate_dict['down']
        regulate_dict['up_down'] = up_downs
        regulate_dict['up'] = regulate_dict['up'] - regulate_dict['down']
        regulate_dict['down'] = regulate_dict['down'] - regulate_dict['up']
        self.logger.info('......dict:%s' % regulate_dict)
        pathways = self.output_dir + '/pathways'
        if not os.path.exists(pathways):
            os.mkdir(pathways)
        try:
            kegg().get_pictrue(path_ko=path_ko, out_dir=pathways, regulate_dict=regulate_dict)
            kegg().get_regulate_table(ko_gene=ko_genes, path_ko=path_ko, regulate_gene=regulate_gene, output=self.output_dir + '/kegg_regulate_stat.xls')
            self.end()
        except Exception:
            import traceback
            self.set_error('kegg调控分析失败：{}'.format(traceback.format_exc()))
