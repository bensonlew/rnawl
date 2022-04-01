# -*- coding: utf-8 -*-
# __author__ = 'shenghe'

from biocluster.module import Module
import os
from biocluster.core.exceptions import OptionError


class BetaDiversityModule(Module):

    def __init__(self, work_id):
        super(BetaDiversityModule, self).__init__(work_id)
        options = [
            {"name": "analysis", "type": "string","default": "pca"},  ### 分析名称
            {"name": "abundance_table", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"}, ### 丰度表
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"}, ###group表
            {"name": "dis_matrix", "type": "outfile", "format": "meta.beta_diversity.distance_matrix"}, ###距离矩阵
            {"name": "dis_method", "type": "string", "default": "bray_curtis"}, ###计算距离矩阵的方法
            {"name": "hcluster_method", "type": "string", "default": "average"}, ###层级聚类的方法
            {"name": "sample_method", "type": "string", "default": "sum"},  # 样本的合并方式, ""为不进行合并
            {'name': 'corr_method', 'type': 'string', 'default': 'pearsonr'}, #相关性heatmap的计算方法
            {"name": "ellipse", "type": "string", "default": "F"},
            {"name": "is_group", "type": "string", "default": "no"},
            {"name": "group_method", "type": "string", "default": ""} ## 这里的group_method有四种["", "sum", "average", "middle"]
        ]
        self.add_option(options)


    def check_options(self):
        """
        进行参数的二次检查
        :return:
        """
        if not self.option("analysis"):
            raise OptionError('没有提供分析类型')
        if not self.option("abundance_table").is_set:
            raise OptionError('没有提供分析所需的otu表')
        # if not self.option("group").is_set:
        #     raise OptionError('没有分组表')
        return True

    def sort_samples_run(self):
        """
        运行计算各分组的结果，之前确定的是按照分组进行计算，后面不再按照分组计算,
        筛选样本的数据也是在这里进行，默认的方法不进行分组计算
        :return:
        """
        opts = ({
            'group_table': self.option('group'),
            'abundance_file': self.option("abundance_table"),
            })
        if self.option("group_method") and self.option("is_group") in ["yes"]:## 为的是防止前端传的是null或者none等字段没有转义
            opts["method"] = self.option("group_method")
        else:
            opts["method"] = ''
        self.sort_samples.set_options(opts)
        self.sort_samples.run()

    def matrix_run(self):
        """
        运行计算距离矩阵
        :return:
        """
        opts = ({
            'method': self.option('dis_method'),
            })
        if self.option("analysis") in ['hcluster', "correlation", "pcoa", "nmds"]:
            opts['otutable'] = self.sort_samples.option("out_otu_table")
        else:
            opts['otutable'] = self.option("abundance_table")
        self.matrix.set_options(opts)
        self.matrix.run()

    def pcoa_run(self):
        """
        运行pcoa分析
        :return:
        """
        pcoa_options = {
            'dis_matrix': self.matrix.option('dis_matrix'),
            'group': self.option("group"),
            'ellipse': self.option('ellipse')
        }
        self.pcoa.set_options(pcoa_options)
        self.pcoa.on('end', self.set_output, 'pcoa'),
        self.pcoa.run()

    def nmds_run(self):
        """
        运行nmds分析
        :return:
        """
        nmds_options = {
            'dis_matrix': self.matrix.option('dis_matrix'),
            'group': self.option("group"),
            'ellipse': self.option('ellipse')
            }
        self.nmds.set_options(nmds_options)
        self.nmds.on('end', self.set_output, 'nmds')
        self.nmds.run()

    def pca_run(self):
        """
        运行pca分析
        :return:
        """
        pca_options = {
            'otutable': self.sort_samples.option("out_otu_table"),
            'group_table': self.option("group"),
            'ellipse': self.option('ellipse')
        }
        self.pca.set_options(pca_options)
        self.pca.on('end', self.set_output, 'pca')
        self.pca.run()

    def hcluster_run(self):
        self.hcluster.set_options({
            'dis_matrix': self.matrix.option('dis_matrix'),
            'linkage': self.option('hcluster_method')
        })
        self.hcluster.on('end', self.set_output, 'hcluster')
        self.hcluster.run()

    def correlation_run(self):
        opts = ({
            'dis_matrix': self.matrix.option('dis_matrix'),
            'hcluster_method': self.option('hcluster_method'),
            'corr_method' : self.option('corr_method')
            })
        if self.option("group").is_set:
            opts['otu_table'] = self.sort_samples.option("out_otu_table").prop['path']
        else:
            opts['otu_table'] = self.option("abundance_table")
        self.correlation.set_options(opts)
        self.correlation.on('end', self.set_output, 'correlation')
        self.correlation.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'pca':
            self.linkdir(obj.output_dir, 'Pca')
        elif event['data'] == 'distance':
            self.linkdir(obj.output_dir, 'correlation')
            self.option('correlation', obj.option('correlation'))
        elif event['data'] == 'hcluster':
            self.linkdir(obj.output_dir, 'Hcluster')
            for file in os.listdir(self.matrix.output_dir):
                file_path = os.path.join(self.matrix.output_dir, file)
                newfile_path = os.path.join(self.output_dir + '/Hcluster', "hcluster_matrix.xls")
                if os.path.exists(newfile_path):
                    os.remove(newfile_path)
                os.link(file_path, newfile_path)
        elif event['data'] == 'pcoa':
            self.linkdir(obj.output_dir, 'Pcoa')
        elif event['data'] == 'nmds':
            self.linkdir(obj.output_dir, 'Nmds')
        elif event['data'] == 'correlation':
            self.linkdir(obj.output_dir, 'Correlation')
        else:
            pass
        self.end()

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
        super(BetaDiversityModule, self).run()
        self.matrix = self.add_tool('meta.beta_diversity.distance_calc')
        if self.option("analysis") in ["pca"]:
            self.sort_samples = self.add_tool("bac_comp_genome.sort_samples")
            self.pca = self.add_tool('meta.beta_diversity.pca')
            self.sort_samples.on("end", self.pca_run)
            self.sort_samples_run()
        elif self.option("analysis") in ["pcoa"]:
            self.pcoa = self.add_tool('meta.beta_diversity.pcoa')
            self.sort_samples = self.add_tool("bac_comp_genome.sort_samples")
            self.matrix.on("end", self.pcoa_run)
            self.sort_samples.on("end", self.matrix_run)
            self.sort_samples_run()
        elif self.option("analysis") in ["nmds"]:
            self.sort_samples = self.add_tool("bac_comp_genome.sort_samples")
            self.nmds = self.add_tool('meta.beta_diversity.nmds')
            self.sort_samples.on("end", self.matrix_run)
            self.matrix.on("end", self.nmds_run)
            self.sort_samples_run()
        elif self.option("analysis") in ["hcluster"]:
            self.hcluster = self.add_tool("meta.beta_diversity.hcluster")
            self.sort_samples = self.add_tool("bac_comp_genome.sort_samples")
            self.sort_samples.on("end", self.matrix_run)
            self.matrix.on("end", self.hcluster_run)
            self.sort_samples_run()
        elif self.option("analysis") in ["correlation"]:
            self.correlation = self.add_tool('bac_comp_genome.correlation')
            self.sort_samples = self.add_tool("bac_comp_genome.sort_samples")
            self.sort_samples.on("end", self.matrix_run)
            self.matrix.on("end", self.correlation_run)
            self.sort_samples_run()


    def end(self):
        super(BetaDiversityModule, self).end()
