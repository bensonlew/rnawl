# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'

from biocluster.module import Module
import os
from biocluster.core.exceptions import OptionError


class NmdsModule(Module):
    """
    小工具nmds模块
    version 1.0
    author zhouxuan
    last_modified:20170615
    """
    def __init__(self, work_id):
        super(NmdsModule, self).__init__(work_id)
        self.step.add_steps('distance_calc', 'nmds')
        options = [
            {"name": "dis_method", "type": "string", "default": "bray_curtis"},
            {"name": "otu_table", "type": "infile", "format": "toolapps.table, meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "T", "type": "string", "default": "column"},
            {"name": "group_table", "type": "infile", "format": "toolapps.group_table"},  # modify by zengjing 20170911
        ]
        self.add_option(options)
        self.matrix = self.add_tool('meta.beta_diversity.distance_calc')
        self.nmds = self.add_tool('meta.beta_diversity.nmds')
        self.METHOD = ['abund_jaccard', 'binary_chisq', 'binary_chord',
                       'binary_euclidean', 'binary_hamming', 'binary_jaccard',
                       'binary_lennon', 'binary_ochiai',
                       'binary_pearson', 'binary_sorensen_dice',
                       'bray_curtis', 'bray_curtis_faith', 'bray_curtis_magurran',
                       'canberra', 'chisq', 'chord', 'euclidean', 'gower',
                       'hellinger', 'kulczynski', 'manhattan', 'morisita_horn',
                       'pearson', 'soergel', 'spearman_approx', 'specprof',
                       'unifrac', 'unweighted_unifrac', 'unweighted_unifrac_full_tree',
                       'weighted_normalized_unifrac', 'weighted_unifrac']
        self.otu_table = ''

    def check_options(self):
        if self.option('T') not in ['row', 'column']:
            raise OptionError("请选择正确的取值方式")
        if self.option('dis_method') not in self.METHOD:
            raise OptionError('请选择正确的距离算法')
        if not self.option('otu_table').is_set:
            raise OptionError('请正确设置输入文件')
        if self.option('T') == 'row':
            self.otu_table = self.t_table(self.option('otu_table').prop['new_table'])
        else:
            self.otu_table = self.option('otu_table').prop['new_table']
        if self.option("group_table").is_set:
            if self.option('T') == 'row':
                for i in self.option('group_table').prop['sample_name']:
                    if i not in self.option('otu_table').prop['row_sample']:
                        raise OptionError('分组文件中的样本{}不存在于表格中，查看是否是数据取值选择错误'.format(i))
            else:
                for i in self.option('group_table').prop['sample_name']:
                    if i not in self.option('otu_table').prop['col_sample']:
                        raise OptionError('分组文件中的样本{}不存在于表格中，查看是否是数据取值选择错误'.format(i))

    def t_table(self, table_file):
        new_table = self.work_dir + '/new_otu_table.xls'
        with open(table_file) as f, open(new_table, 'w') as w:
            table_list = [i.rstrip().split('\t') for i in f.readlines()]
            table_list = map(lambda *a: '\t'.join(a) + '\n', *table_list)
            w.writelines(table_list)
        return new_table

    def matrix_run(self):
        """
        运行计算距离矩阵
        :return:
        """
        self.matrix.set_options({
            'method': self.option('dis_method'),
            'otutable': self.otu_table
        })
        self.matrix.on('end', self.set_output, 'distance')
        self.matrix.run()

    def nmds_run(self):
        self.nmds.set_options({
            'dis_matrix': self.matrix.option('dis_matrix')
        })
        self.nmds.on('end', self.set_output, 'nmds')
        self.nmds.on('end', self.end)
        self.nmds.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'distance':
            self.linkdir(obj.output_dir, 'Distance')
        elif event['data'] == 'nmds':
            self.linkdir(obj.output_dir, 'Nmds')
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
        super(NmdsModule, self).run()
        self.matrix.on('end', self.nmds_run)
        self.matrix_run()

    def end(self):
        repaths = [
            [".", "", "Nmds分析结果文件目录"],
            ["Distance", "", "距离矩阵计算结果输出目录"],
            ["Nmds", "", "NMDS分析结果输出目录"],
            ["Nmds/nmds_sites.xls", "xls", "样本坐标表"],
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
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(NmdsModule, self).end()
