# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess
import os


class DiscomparisonAgent(Agent):
    """
    discomparison:用于检验群落距离矩阵和环境变量距离矩阵之间的相关性
    version: 1.0
    author: wangbixuan
    last_modified: 20160914 by qindanhua
    """
    MATRIX = ['abund_jaccard', 'binary_chisq', 'binary_chord', 'binary_euclidean', 'binary_hamming', 'binary_jaccard',
              'binary_lennon', 'binary_ochiai', 'binary_otu_gain', 'binary_pearson', 'binary_sorensen_dice',
              'bray_curtis', 'bray_curtis_faith', 'bray_curtis_magurran', 'canberra', 'chisq', 'chord', 'euclidean',
              'gower', 'hellinger', 'kulczynski', 'manhattan', 'morisita_horn', 'pearson', 'soergel', 'spearman_approx',
              'specprof', 'unifrac', 'unweighted_unifrac', 'weighted_normalized_unifrac', 'weighted_unifrac']
    MATRIXFACTOR = ['abund_jaccard', 'binary_chisq', 'binary_chord', 'binary_euclidean', 'binary_hamming',
                    'binary_jaccard', 'binary_lennon', 'binary_ochiai', 'binary_otu_gain', 'binary_pearson',
                    'binary_sorensen_dice', 'bray_curtis', 'bray_curtis_faith', 'bray_curtis_magurran', 'canberra',
                    'chisq', 'chord', 'euclidean', 'gower', 'hellinger', 'kulczynski', 'manhattan', 'morisita_horn',
                    'pearson', 'soergel', 'spearman_approx', 'specprof']

    def __init__(self, parent):
        super(DiscomparisonAgent, self).__init__(parent)
        options = [
            {"name": "partialmatrix", "type": "infile", "format": "meta.beta_diversity.distance_matrix"},
            {'name': 'otudistance', 'type': 'infile', 'format': 'meta.beta_diversity.distance_matrix'},
            {'name': 'facdistance', 'type': 'infile', 'format': 'meta.beta_diversity.distance_matrix'}
        ]
        self.add_option(options)
        self.step.add_steps('mantel_test')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.mantel_test.start()
        self.step.update()

    def step_end(self):
        self.step.mantel_test.finish()
        self.step.update()

    def gettable(self):
        """
        get matrix for calculation by level provided
        """
        if self.option("otutable").format == "meta.otu.tax_summary_dir":
            return self.option("otutable").get_table(self.option('level'))
        else:
            return self.option('otutable').prop['path']

    def check_options(self):
        if not self.option('otudistance').is_set:
            raise OptionError('必须提供otu距离表', code="34100301")
        if not self.option('facdistance').is_set:
            raise OptionError('必须提供环境因子距离表', code="34100302")

    def set_resource(self):
        self._cpu = 5
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "Discomparison计算结果输出目录"],
            ["./mantel_results.txt", "txt", "Discomparison结果"]
        ])
        super(DiscomparisonAgent, self).end()


class DiscomparisonTool(Tool):
    def __init__(self, config):
        super(DiscomparisonTool, self).__init__(config)
        self._version = '1.9.1'
        self.cmd_path = 'program/Python/bin/compare_distance_matrices.py'
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + 'gcc/5.1.0/lib64:$LD_LIBRARY_PATH')

    def run(self):
        """
        运行
        """
        super(DiscomparisonTool, self).run()
        self.run_discomparison()
        # self.set_output()
        # self.end()

    def run_discomparison(self):
        """
        run compare_distance_matrices.py
        """
        cmd = self.cmd_path
        if self.option('partialmatrix').is_set:
            cmd += ' -i %s,%s -c %s -o %s --method partial_mantel -n 999' % \
                   (self.option('otudistance').prop['path'], self.option('facdistance').prop['path'],
                    self.option('partialmatrix').prop['path'], self.work_dir)
        else:
            cmd += ' -i %s,%s -o %s --method mantel -n 999' % (self.option('otudistance').prop['path'],
                                                               self.option('facdistance').prop['path'], self.work_dir)
        self.logger.info('运行compare_distance_matrices.py 判断相关性')
        self.logger.info(cmd)
        discomparison_command = self.add_command('distance_comparision', cmd)
        discomparison_command.run()
        self.wait()
        if discomparison_command.return_code == 0:
            self.logger.info('running compare_distance_matrices.py succeed')
            if self.option('partialmatrix').is_set:
                filename = self.work_dir+'/partial_mantel_results.txt'
                linkfile = self.output_dir+'/partial_mantel_results.txt'
            else:
                filename = self.work_dir+'/mantel_results.txt'
                linkfile = self.output_dir+'/mantel_results.txt'
            if os.path.exists(linkfile):
                os.remove(linkfile)
            os.link(filename, linkfile)
            self.end()
        else:
            self.set_error('Error in running compare_distance_matrices.py', code="34100301")
