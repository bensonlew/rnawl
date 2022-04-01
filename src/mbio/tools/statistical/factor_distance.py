# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from mbio.files.meta.otu.otu_table import OtuTableFile
import re


class FactorDistanceAgent(Agent):
    """
    calculate distance of factor matrix according to given matrix type
    author: wangbixuan
    last modified: 20160914 by qindanhua
    """
    MATRIXFACTOR = ['abund_jaccard', 'binary_chisq', 'binary_chord', 'binary_euclidean', 'binary_hamming',
                    'binary_jaccard', 'binary_lennon', 'binary_ochiai', 'binary_otu_gain', 'binary_pearson',
                    'binary_sorensen_dice', 'bray_curtis', 'bray_curtis_faith', 'bray_curtis_magurran', 'canberra',
                    'chisq', 'chord', 'euclidean', 'gower', 'hellinger', 'kulczynski', 'manhattan', 'morisita_horn',
                    'pearson', 'soergel', 'spearman_approx', 'specprof']

    def __init__(self, parent):
        super(FactorDistanceAgent, self).__init__(parent)
        options = [
            {"name": "factor", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "facmatrixtype", "type": "string", "default": "bray_curtis"},
            {"name": "factorselected", "type": "string", "default": ""},
            {"name": "fac_matrix", "type": "outfile", "format": "meta.beta_diversity.distance_matrix"}
        ]
        self.add_option(options)
        self.step.add_steps('factor_distance')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.factor_distance.start()
        self.step.update()

    def step_end(self):
        self.step.factor_distance.finish()
        self.step.update()

    def check_options(self):
        if not self.option('factor').is_set:
            raise OptionError('没有提供环境因子表')
        else:
            self.option('factor').get_info()
            if self.option('factorselected'):
                factors = self.option('factorselected').split(',')
                for f in factors:
                    if f not in self.option('factor').prop['group_scheme']:
                        raise OptionError('该环境因子在输入的环境因子表里不存在：%s' %f)
            else:
                pass
        if self.option('facmatrixtype') not in FactorDistanceAgent.MATRIXFACTOR:
            raise OptionError(' 不支持所选矩阵类型 .')

    def set_resource(self):
        self._cpu = 5
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "环境因子距离矩阵结果输出目录"],
            ])
        result_dir.add_regexp_rules([
            [r'%s.*\.xls' % self.option('facmatrixtype'), 'xls', '环境因子距离矩阵文件']
            ])
        super(FactorDistanceAgent, self).end()


class FactorDistanceTool(Tool):
    def __init__(self, config):
        super(FactorDistanceTool, self).__init__(config)
        self._version = '1.9.1'  # qiime版本
        self.cmd_path = 'program/Python/bin/beta_diversity.py'
        self.biom = None
        # 设置运行环境变量
        # self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + 'gcc/5.1.0/lib64:$LD_LIBRARY_PATH')
        # self.biom = self.biom_fac_table()  # 传入otu表需要转化为biom格式

    def run(self):
        super(FactorDistanceTool, self).run()
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + 'gcc/5.1.0/lib64')
        self.biom = self.biom_fac_table()  # 传入otu表需要转化为biom格式
        self.run_beta_diversity()

    def run_beta_diversity(self):
        cmd = self.cmd_path
        cmd += ' -m %s -i %s -o %s' % (self.option('facmatrixtype'), self.biom, self.work_dir)
        self.logger.info('run beta_diversity.py program')
        fac_matrix_command = self.add_command('fac_matrix', cmd)
        fac_matrix_command.run()
        self.wait()
        if fac_matrix_command.return_code == 0:
            self.logger.info('Succeed on calculating factor matrix')
            filename = self.work_dir + '/' + \
                self.option('facmatrixtype') + '_temp.txt'
            linkfile = self.output_dir+'/factor_out'+'.xls'
            if os.path.exists(linkfile):
                os.remove(linkfile)
            os.link(filename, linkfile)
            self.option('fac_matrix', linkfile)
            self.end()
        else:
            self.set_error('Error in running beta_diversity.py')

    def biom_fac_table(self):
        fac_index = []
        write_line_list = ["#name"]
        with open('transtable.txt', 'w') as w, open(self.option('factor').prop['path'], "r") as f:
            first_line = f.readline().strip().split()
            if self.option('factorselected'):
                selectedf = self.option('factorselected').split(',')
            else:
                selectedf = first_line[1:]
            self.logger.info(selectedf)
            for fac in selectedf:
                if fac in first_line[1:]:
                    fac_index.append(first_line.index(fac))
                    write_line_list.append(fac)
            self.logger.info(fac_index)
            for line in f:
                line = line.strip().split()
                n = 0
                write_line_list[n] += "\t{}".format(line[0])
                for i in fac_index:
                    n += 1
                    # if re.match(r"\D", line[i]):
                    #     continue
                    # else:
                    write_line_list[n] += "\t{}".format(line[i])
            write_line_len = len(write_line_list[0].split("\t"))
            for write_line in write_line_list:
                if len(write_line.split("\t")) < write_line_len:
                    continue
                else:
                    w.write(write_line)
                    w.write("\n")

        trans_newtable = OtuTableFile()
        trans_newtable.set_path('transtable.txt')
        biom_path = os.path.join(self.work_dir, 'temp.biom')
        if os.path.isfile(biom_path):
            os.remove(biom_path)
        # newtable.convert_to_biom(biom_path)
        trans_newtable.check()
        trans_newtable.convert_to_biom(biom_path)
        return biom_path
