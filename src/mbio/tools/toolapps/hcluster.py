# -*- coding: utf-8 -*-
# __author__ = 'zhangpeng'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import re
import subprocess
from biocluster.core.exceptions import OptionError


class HclusterAgent(Agent):
    """
    脚本plot-hcluster_tree_app.pl
    version v2.0
    author: zhangpeng
    last_modified:2017.8.28 zhouxuan
    """

    def __init__(self, parent):
        super(HclusterAgent, self).__init__(parent)
        options = [
            {"name": "otu_table", "type": "infile", "format": "toolapps.table"},  # modify by zhouxuan 20170623
            {"name": "linkage", "type": "string", "default": "average"},
            {"name": "method", "type": "string", "default": "euclidean"},
            {"name": "trans", "type": "string", "default": "column"},
            {"name": "group_table", "type": "infile", "format": "toolapps.group_table"}
        ]
        self.add_option(options)
        self.step.add_steps('hcluster')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.hcluster.start()
        self.step.update()

    def step_end(self):
        self.step.hcluster.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('otu_table').is_set:
            raise OptionError('请设置数据表进行分析')
        if self.option('method') not in ['euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski']:
            raise OptionError('错误的距离方式：%s' % self.option('method'))
        if self.option('trans') not in ['column', 'row']:  # modify by zhouxuan
            raise OptionError('错误的距离方式：%s' % self.option('trans'))
        if self.option('linkage') not in ['average', 'single', 'complete']:
            raise OptionError('错误的层级聚类方式：%s' % self.option('linkage'))
        if self.option('group_table').is_set:  # 当有分组的时候，判断分组中的样本是否存在于数据表中
            if self.option('trans') == 'column':
                for i in self.option('group_table').prop['sample_name']:
                    if i not in self.option('otu_table').prop['col_sample']:
                        raise OptionError('分组文件中的样本不存在于表格中，查看是否是数据取值选择错误')
            else:
                self.logger.info(self.option('group_table').prop['sample_name'])
                self.logger.info(self.option('otu_table').prop['row_sample'])
                for i in self.option('group_table').prop['sample_name']:
                    if i not in self.option('otu_table').prop['row_sample']:
                        raise OptionError('分组文件中的样本不存在于表格中，查看是否是数据取值选择错误')

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        self._memory = '10G'  # 暂时增大内存 by xieshichang 20200628

    def end(self):
        if self.option('group_table').is_set:
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "PCA分析结果输出目录"],
            ])
            result_dir.add_regexp_rules([
                [".+/hcluster.tre", "tre", "层次聚类树"],
                [".+/data_table", 'txt', "数据表"]
            ])
        else:
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "层次聚类结果目录"],
                ["./hcluster.tre", "tre", "层次聚类树"],
                ["./data_table", 'txt', "数据表"]
            ])
        super(HclusterAgent, self).end()


class HclusterTool(Tool):

    def __init__(self, config):
        super(HclusterTool, self).__init__(config)
        self._version = 'v2.1-20140214'  # plot-hcluster_tree.pl版本
        self.cmd_path = os.path.join(
            self.config.SOFTWARE_DIR, 'bioinfo/statistical/scripts/plot-hcluster_tree_app.pl')

    def run(self):
        """
        运行
        """
        super(HclusterTool, self).run()
        if self.option('group_table').is_set:  # 当有分组的时候，判断分组中的样本是否存在于数据表中
            if self.option('trans') == 'column':
                self.data_table = self.option('otu_table').prop['new_table']  # 不进行转置
            else:
                self.data_table = self.work_dir + '/T_table.txt'
                self.t_table(self.option('otu_table').prop['new_table'], self.data_table)
            group_de = self.option('group_table').prop['group_scheme']  # 获取分组文件中的子分组
            self.logger.info(group_de)
            for i in group_de:
                name = []
                name.append(i)
                group_target_path = os.path.join(self.work_dir, i + '_group.xls')
                self.option('group_table').sub_group(group_target_path, name)  # 子分组
                data_target_path = os.path.join(self.work_dir, i + '_data_table.xls')
                self.option('otu_table').get_table_of_main_table(self.data_table, data_target_path, group_target_path)
                self.run_hcluster(data_target_path, i)
        else:
            self.run_hcluster(self.option('otu_table').prop['new_table'])

    def run_hcluster(self, data_table=None, group=None):
        """
        运行plot-hcluster_tree.pl
        """
        if group:
            tmp_name = group + '_distance_matrix.temp'
        else:
            tmp_name = 'distance_matrix.temp'
        real_dis_matrix = os.path.join(self.work_dir, tmp_name)
        self.newname_dict = self.change_sample_name(data_table=data_table, quotes=False, new_path=real_dis_matrix)
        # 修改矩阵的样本名称为不含特殊符号的名称，返回一个旧名称对新名称的字典
        cmd = self.cmd_path
        if self.option('group_table').is_set:
            cmd += ' -i %s -o %s -m %s -l %s -trans col -m_1 %s ' % \
                   (real_dis_matrix, self.work_dir, self.option('linkage'), self.option('method'),
                    self.option('linkage'))
        else:
            if self.option('trans') == 'column':  # 修改行列取值的实际方式，为了和小工具pca保持一致
                cmd += ' -i %s -o %s -m %s -l %s -trans row -m_1 %s ' %\
                       (real_dis_matrix, self.work_dir, self.option('linkage'), self.option('method'), self.option('linkage'))
            else:
                cmd += ' -i %s -o %s -m %s -l %s -trans col -m_1 %s ' % \
                       (real_dis_matrix, self.work_dir, self.option('linkage'), self.option('method'), self.option('linkage'))
        self.logger.info('运行plot-hcluster_tree.pl程序计算Hcluster')
        self.logger.info(cmd)
        try:  # 生成脚本
            subprocess.check_output(cmd, shell=True)
            self.logger.info('生成 hc.cmd.r 文件成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成 hc.cmd.r 文件失败')
            self.set_error('无法生成 hc.cmd.r 文件')
        self.logger.info(self.config.SOFTWARE_DIR + '/program/R-3.3.1/bin/R --restore --no-save < %s/hc.cmd.r' % self.work_dir)
        try:  # 运行脚本
            subprocess.check_output(self.config.SOFTWARE_DIR +
                                    '/program/R-3.3.1/bin/R --restore --no-save < %s/hc.cmd.r' % self.work_dir, shell=True)
            self.logger.info('生成树文件成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成树文件失败')
            raise Exception("数据量太大无法生成树文件")
        filename = self.work_dir + '/hcluster_tree_' + os.path.basename(real_dis_matrix) + '_' + self.option('linkage') + '.tre'  # 结果文件
        if group:
            os.rename(self.work_dir + '/hc.cmd.r', self.work_dir + '/' + group + '_hc.cmd.r')
            os.mkdir(self.output_dir + '/' + group)
            linkfile = self.output_dir + '/' + group +'/hcluster.tre'  # link到结果文件夹的最终文件
        else:
            linkfile = self.output_dir + '/hcluster.tre'  # link到结果文件夹的最终文件
        self.re_recover_name(self.newname_dict, filename, filename + '.temp')  # 换回名称
        # if os.path.exists(linkfile):
        #     os.remove(linkfile)
        os.link(filename + '.temp', linkfile)
        if group:
            os.link(data_table, os.path.join(self.output_dir + '/' + group, "data_table"))
        else:
            os.link(data_table, os.path.join(self.output_dir, "data_table"))
        self.end()

    def change_sample_name(self, data_table, quotes=False, new_path=None):
        """
        修改矩阵的样本名称为不含特殊符号的名称，返回一个旧名称对新名称的字典
        """
        if not new_path:
            new_path = self.work_dir + '/distance_matrix.temp'
        old_matrix = open(data_table, 'rb')
        name_dict = {}
        new_matrix = open(new_path, 'wb')
        frist_line = old_matrix.readline().rstrip().split('\t')[1:]
        if quotes:
            name_dict = {('\"' + name + '\"'): ('name' + str(frist_line.index(name))) for name in frist_line}
        else:
            name_dict = {name: ('name' + str(frist_line.index(name))) for name in frist_line}
        new_matrix.write('\t' + '\t'.join(name_dict.itervalues()) + '\n')
        for line in old_matrix:
            line_split = line.split('\t')
            new_matrix.write('\t'.join(line_split))
        old_matrix.close()
        new_matrix.close()
        return name_dict

    def recover_name(self, namedict, treefile, newfile):
        """
        复原树文件中的名称
        """
        from Bio import Phylo
        from Bio.Phylo.NewickIO import NewickError
        if not isinstance(namedict, dict):
            raise Exception('复原树的枝名称需要旧名称和当前名称的字典')
        namedict = {item[1]: item[0] for item in namedict.iteritems()}
        if not isinstance(treefile, (str, unicode)):
            raise Exception('树文件的路径不是字符串')
        if not isinstance(newfile, (str, unicode)):
            raise Exception('新的树文件的路径不是字符串')
        try:
            tree = Phylo.read(treefile, 'newick')
        except IOError:
            raise Exception('复原树文件时找不到树文件：%s' % treefile)
        except NewickError:
            raise Exception('树文件无法用newick格式解析：%s' % treefile)
        terminals = tree.get_terminals()
        for terminal in terminals:
            if terminal.name not in namedict:
                raise Exception('树的枝名称：%s在旧名称和新名称字典中不存在' % terminal.name)
            terminal.name = namedict[terminal.name]
        Phylo.write(tree, newfile, 'newick')
        return True

    def re_recover_name(self, namedict, treefile, newfile):
        """
        用正则的方式替换复原树文件中的名称
        """
        if not isinstance(namedict, dict):
            raise Exception('复原树的枝名称需要旧名称和当前名称的字典')
        namedict = {item[1]: item[0] for item in namedict.iteritems()}
        if not isinstance(treefile, (str, unicode)):
            raise Exception('树文件的路径不是字符串')
        if not isinstance(newfile, (str, unicode)):
            raise Exception('新的树文件的路径不是字符串')
        try:
            with open(treefile, 'rb') as f, open(newfile, 'wb') as w:
                tree = f.readline()
                for item in namedict.iteritems():
                    tree = re.sub(item[0] + ':', item[1] + ':', tree)
                w.write(tree)
        except IOError, e:
                raise Exception('聚类树文件无法找到或者无法打开：%s' % e)

    def t_table(self, table_file, new_table):  # 表格转置
        """
		转换颠倒表格内容
		"""
        with open(table_file) as f, open(new_table, 'w') as w:
            table_list = [i.rstrip().split('\t') for i in f.readlines()]
            table_list = map(lambda *a: '\t'.join(a) + '\n', *table_list)
            w.writelines(table_list)
