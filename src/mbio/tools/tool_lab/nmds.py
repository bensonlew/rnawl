# -*- coding: utf-8 -*-
# __author__ = 'linna'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import glob
import os
import pandas as pd


class NmdsAgent(Agent):
    METHOD = ['abund_jaccard', 'binary_chisq', 'binary_chord',
              'binary_euclidean', 'binary_hamming', 'binary_jaccard',
              'binary_lennon', 'binary_ochiai',
              'binary_pearson', 'binary_sorensen_dice',
              'bray_curtis', 'bray_curtis_faith', 'bray_curtis_magurran',
              'canberra', 'chisq', 'chord', 'euclidean', 'gower',
              'hellinger', 'kulczynski', 'manhattan', 'morisita_horn',
              'pearson', 'soergel', 'spearman_approx', 'specprof']

    def __init__(self, parent):
        super(NmdsAgent, self).__init__(parent)
        options = [
            {"name": "tooltable", "type": "infile", "format": "tool_lab.table"},
            {"name": "specimen_name", "type": "string", "default": "column"},
            {"name": "group_table", "type": "infile", "format": "toolapps.group_table"},
            {"name": "scale", "type": "string", "default": "False"},
            {"name": "dis_method", "type": "string", "default": "bray_curtis"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数检查
        """
        if not self.option('tooltable').is_set:
            raise OptionError('必须提供数据表', code="32702903")
        self.option('tooltable').get_info()
        if self.option('tooltable').prop['sample_num'] < 3:
            raise OptionError('列数少于3，不可进行分析', code="32702904")
        if self.option('dis_method') not in NmdsAgent.METHOD:
            raise OptionError('错误或者不支持该距离矩阵计算方法', code="32701903")
        # if self.option('group_table').is_set:
            # if self.option('specimen_name') == 'column':
            #     for i in self.option('group_table').prop['sample_name']:
            #         if i not in self.option('tooltable').prop['col_sample']:
            #             raise OptionError('分组文件中的样本不存在于表格中，查看是否数据选择错误', code="32702909")
            # else:
            #     for i in self.option('group_table').prop['sample_name']:
            #         if i not in self.option('tooltable').prop['row_sample']:
            #             raise OptionError('分组文件中的样本不存在于表格中，查看是否数据选择错误', code="32702910")
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 5
        self._memory = '10G'

    def end(self):
        if self.option('group_table').is_set:
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "Nmds分析结果输出目录"],
                ["./nmds_sites.xls", "xls", "样本坐标表"],
                ["./nmds_stress.xls", "xls", "样本特征拟合度值"],
                ["./group.xls", "xls", "分组表"],
            ])
        else:
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "Nmds分析结果输出目录"],
                ["./nmds_sites.xls", "xls", "样本坐标表"],
                ["./nmds_stress.xls", "xls", "样本特征拟合度值"],
            ])
        super(NmdsAgent, self).end()


class NmdsTool(Tool):
    def __init__(self, config):
        super(NmdsTool, self).__init__(config)
        self.cmdpath = 'program/Python/bin/beta_diversity.py'
        self.biom_path = 'program/Python/bin/'
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.perl_path = "program/perl-5.24.0/bin/perl"
        self.cmd_path = self.config.PACKAGE_DIR + '/statistical/ordinationdis.pl'
        self.script_path = "bioinfo/meta/scripts/beta_diver.sh"
        self.R_path = os.path.join(self.config.SOFTWARE_DIR, 'program/R-3.3.1/bin/R')
        self.R_path1 = '/program/R-3.3.1/bin/'
        self.ellipse_path = self.config.PACKAGE_DIR + '/graph/scripts/Ellipse.R'

    def run(self):
        """
        运行
        """
        super(NmdsTool, self).run()
        self.convert_to_biom()
        self.run_beta_diversity()
        self.run_ordination()
        self.linkfile()
        self.end()

    def formattable(self, tablepath):
        """
        转置表格
        """
        this_table = tablepath
        if self.option('specimen_name') != 'column':
            newtable = this_table + '.T'
            self.t_table(this_table, newtable)
            return newtable
        else:
            return this_table

    def t_table(self, table_file, new_table):
        """
        转置表格实现
        """
        with open(table_file) as f, open(new_table, 'w') as w:
            table_list = [i.rstrip().split('\t') for i in f.readlines()]
            table_list = map(lambda *a: '\t'.join(a) + '\n', *table_list)
            w.writelines(table_list)

    def run_beta_diversity(self):
        """
        运行qiime:beta_diversity.py
        """
        cmd = self.cmdpath
        cmd += ' -m %s -i %s -o %s' % (self.option('dis_method'), self.output_dir + '/temp.biom', self.work_dir)
        self.logger.info('运行qiime:beta_diversity.py程序')
        self.logger.info(cmd)
        dist_matrix_command = self.add_command('distance_matrix', cmd)
        dist_matrix_command.run()
        self.wait()
        if dist_matrix_command.return_code == 0:
            self.command_successful()
        else:
            self.set_error("运行qiime:beta_diversity.py出错", code="32701901")

    def command_successful(self):
        self.logger.info('运行qiime:beta_diversity.py完成')
        filename = self.work_dir + '/' + self.option('dis_method') + '_temp.txt'
        linkfile = self.output_dir + '/distance.xls'
        if os.path.exists(linkfile):
            os.remove(linkfile)
        os.link(filename, linkfile)
        self.dis_matrix_check(linkfile)

    def dis_matrix_check(self, linkfile):
        """聚类矩阵结果检查"""
        self.logger.info("开始进行距离矩阵检查")
        dist_dict = dict()
        all_values = []
        with open(linkfile, 'r') as f:
            head = f.readline().rstrip().split('\t')
            head_len = len(head)
            head = head[1:]
            for line in f:
                all_nums = line.rstrip().split('\t')
                if len(all_nums) != head_len:
                    self.set_error('距离矩阵每行数据量格式不正确', code="32701903")
                values = dict(zip(head, all_nums[1:]))
                all_values.extend(all_nums[1:])
                dist_dict[all_nums[0]] = values
            for samp1 in head:
                for samp2 in head:
                    if dist_dict[samp1][samp2] != dist_dict[samp2][samp1]:
                        self.set_error('距离矩阵数据不对称', code="32701904")
            all_values = [float(i) for i in all_values]
            all_plus = sum(all_values)
            if all_plus == 0 and len(all_values) > 1:  # 只有一个样本时距离就为零，不做处理，但是后续不能做任何分析
                self.set_error('所有距离矩阵值全部为零', code="32701905")
            if len(head) != len(set(head)):
                self.set_error('距离矩阵存在重复的样本名', code="32701906")
        self.logger.info("距离矩阵检查正确")

    def convert_to_biom(self):
        """
        转换为biom格式
        """
        cmd = self.biom_path
        biom_filepath = os.path.join(self.work_dir, 'temp.biom')
        cmd += 'biom convert -i %s -o %s --table-type="OTU table" --to-hdf5' % \
               (self.formattable(self.option('tooltable').path), biom_filepath)
        self.logger.info('运行biom.py程序')
        self.logger.info(cmd)
        biom_command = self.add_command('biom', cmd)
        biom_command.run()
        self.wait()
        if biom_command.return_code == 0:
            self.logger.info('运行biom.py完成')
            linkfile = self.output_dir + '/temp.biom'
            if os.path.exists(linkfile):
                os.remove(linkfile)
            os.link(biom_filepath, linkfile)
        else:
            self.set_error("运行biom.py出错", code="32701901")

    def run_ordination(self, cmd1='cmd', cmd2='nmds', cmd3='sample_ellipse', cmd4='group_ellipse'):
        """
        运行ordinationdis_new.pl计算nmds，调用Ellipse.R生成置信椭圆
        """
        if self.option('scale') == "True":
            # 是否标准化参数传给pl脚本值不对导致不能标准化，这里添加处理  modified by wuqin on 20211027
            scale = "T"
        else:
            scale = "F"
        if self.option('group_table').is_set:
            cmd = '{} {} -type nmds -dist {} -outdir {} -scale {} -group {}'.format(
                self.perl_path, self.cmd_path, self.output_dir + '/distance.xls',
                self.work_dir, scale, self.option('group_table').path)
        else:
            dis_table = pd.read_table(self.output_dir + '/distance.xls')
            # 这里没有分组时，调用的脚本不会画椭圆，所以创建一个分组文件，每个样本的分组名都为All
            group_table = pd.DataFrame()
            group_table["name"] = dis_table[list(dis_table.columns)[0]]
            group_table["group"] = "All"
            group_table.to_csv(self.work_dir + "/sample.xls", index=0, sep="\t")
            cmd = '{} {} -type nmds -dist {} -outdir {} -scale {} -group {}'.format(
                self.perl_path, self.cmd_path, self.output_dir + '/distance.xls',
                self.work_dir, scale, self.work_dir + "/sample.xls")
        self.logger.info('运行ordinationdis.pl程序')
        self.logger.info(cmd)
        cmd1 = self.add_command(cmd1, cmd).run()  # 命令1运行
        self.wait(cmd1)
        if cmd1.return_code == 0:
            self.logger.info('生成 cmd.r 文件成功')
        else:
            self.logger.info('生成 cmd.r 文件失败')
            self.set_error('无法生成 cmd.r 文件', code="32702902")
        cmd_ = self.script_path + ' %s %s' % (self.R_path, self.work_dir + "/cmd.r")
        self.logger.info(cmd_)
        cmd2 = self.add_command(cmd2, cmd_).run()  # 命令2 运行
        self.wait(cmd2)
        if cmd2.return_code == 0:
            self.logger.info('nmds计算成功')
        else:
            self.logger.info('nmds计算失败')
            self.set_error('R程序计算nmds失败', code="32702903")
        self.logger.info('运行ordinationdis.pl程序计算nmds完成')

        basename = os.path.basename(self.output_dir + '/distance.xls')
        basename1, basename2 = basename.split('.')
        site_file0 = self.work_dir + '/nmds' + '/' + basename1 + '_' + basename2 + '_' + 'nmds_sites.xls'
        site_table = pd.read_table(site_file0)
        site_table.rename(columns={"Sample_ID": "name", "MDS1": "NMDS1", "MDS2": "NMDS2"}, inplace=True)
        site_table.to_csv(self.output_dir + '/nmds_sites1.xls', sep='\t', index=0)
        target_file = self.output_dir + '/nmds_sites1.xls'
        cmd_s = self.R_path1 + 'Rscript {} -f {} -l {} -o {} '.format(self.ellipse_path, target_file, 0.95,
                                                                      self.output_dir + '/sample_ellipse.xls')
        cmd3 = self.add_command(cmd3, cmd_s).run()  # 命令3 运行
        self.wait(cmd3)
        if cmd3.return_code == 0:
            self.logger.info("运行sample_ellipse完成")
        else:
            self.set_error("运行sample_ellipse出错!", code="32301304")

        if self.option('group_table').is_set:
            group_table = self.option('group_table').path
            if group_table != '':
                lines = open(group_table).readlines()
                open(group_table, 'w').writelines('{}\n'.format(line.strip()) for line in lines)
                cmd_g = self.R_path1 + 'Rscript {} -f {} -l {} -o {} -g {}'.format(self.ellipse_path, target_file,
                                                                                   0.95, self.output_dir +
                                                                                   '/group_ellipse.xls', group_table)
                self.logger.info(cmd_g)
                cmd4 = self.add_command(cmd4, cmd_g).run()  # 命令4 运行
                self.wait(cmd4)
                if cmd4.return_code == 0:
                    self.logger.info("运行group_ellipse完成")
                else:
                    self.set_error("运行group_ellipse出错!", code="32301305")
        else:
            sample_ellipse = self.output_dir + '/sample_ellipse.xls'
            group_ellipse = self.output_dir + '/group_ellipse.xls'
            if os.path.exists(group_ellipse):
                os.remove(group_ellipse)
            os.link(sample_ellipse, group_ellipse)

    def linkfile(self):
        if self.option('group_table').is_set:
            group_file = self.option('group_table').path
            group_link = self.output_dir + '/group.xls'
            if os.path.exists(group_link):
                os.remove(group_link)
            os.link(group_file, group_link)
            basename = os.path.basename(self.output_dir + '/distance.xls')
            basename1, basename2 = basename.split('.')
            stress_file = self.work_dir + '/nmds' + '/' + basename1 + '_' + basename2 + '_' + 'nmds_stress.xls'
            stress_link = self.output_dir + '/nmds_stress.xls'
            if os.path.exists(stress_link):
                os.remove(stress_link)
            os.link(stress_file, stress_link)
            site_file = self.work_dir + '/nmds' + '/' + basename1 + '_' + basename2 + '_' + 'nmds_sites.xls'
            site_link = self.output_dir + '/nmds_sites.xls'
            if os.path.exists(site_link):
                os.remove(site_link)
            os.link(site_file, site_link)
            ellipse_file = self.work_dir + "/nmds" + "/" + "ellipse.xls"
            ellipse_link = self.output_dir + '/ellipse.xls'
            if os.path.exists(ellipse_link):
                os.remove(ellipse_link)
            os.link(ellipse_file, ellipse_link)
        else:
            basename = os.path.basename(self.output_dir + '/distance.xls')
            basename1, basename2 = basename.split('.')
            stress_file = self.work_dir + '/nmds' + '/' + basename1 + '_' + basename2 + '_' + 'nmds_stress.xls'
            stress_link = self.output_dir + '/nmds_stress.xls'
            if os.path.exists(stress_link):
                os.remove(stress_link)
            os.link(stress_file, stress_link)
            site_file = self.work_dir + '/nmds' + '/' + basename1 + '_' + basename2 + '_' + 'nmds_sites.xls'
            site_link = self.output_dir + '/nmds_sites.xls'
            if os.path.exists(site_link):
                os.remove(site_link)
            os.link(site_file, site_link)
            ellipse_file = self.work_dir + "/nmds" + "/" + "ellipse.xls"
            ellipse_link = self.output_dir + '/ellipse.xls'
            if os.path.exists(ellipse_link):
                os.remove(ellipse_link)
            os.link(ellipse_file, ellipse_link)
