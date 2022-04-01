# -*- coding: utf-8 -*-
# __author__ = 'shenghe'  
from biocluster.agent import Agent
from biocluster.tool import Tool
import pandas as pd
import os
import types
import subprocess
import re  # add re by guhaidong 20171025
from mbio.packages.taxon.mask_taxon import mask_taxon  # add by guhaidong 20171025
from biocluster.core.exceptions import OptionError
from mbio.files.meta.otu.otu_table import OtuTableFile


class PcaAgent(Agent):
    """
    脚本ordination.pl
    version v1.0
    author: shenghe
    last_modified:2016.3.24
    """

    def __init__(self, parent):
        super(PcaAgent, self).__init__(parent)
        options = [
            {"name": "otutable", "type": "infile","format": "itraq_and_tmt.ratio_exp"},
            # modify by zhouxuan 20170623 小工具的模块是指定toolapps.table这个文件类型的
            {"name": "level", "type": "string", "default": "otu"},
            {"name": "eigenvalue", "type": "string", "default": "row"},  # column
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "envlabs", "type": "string", "default": ""},
            {"name": "group_table", "type": "infile", "format": "itraq_and_tmt.group_table"},  # modify by zhouxuan 20170823
            {"name": "scale", "type": "string", "default": "F"}  # 是否进行标准化 ，add by zouxuan
        ]
        self.add_option(options)
        self.step.add_steps('PCAanalysis')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.PCAanalysis.start()
        self.step.update()

    def step_end(self):
        self.step.PCAanalysis.finish()
        self.step.update()

    def gettable(self):
        """
        根据level返回进行计算的丰度表
        :return:
        """
        if self.option('otutable').format == "meta.otu.tax_summary_dir":
            return self.option('otutable').get_table(self.option('level'))
        else:
            return self.option('otutable').prop['path']

    def check_options(self):
        """
        重写参数检查
        """
        if self.option('scale') not in ['T', 'F']:
            raise OptionError('scale只能为T或者F', code = "32503001")
        if not self.option('otutable').is_set:
            raise OptionError('必须提供数据表', code = "32503002")
        self.option('otutable').get_info()
        if self.option('otutable').prop['sample_number'] < 3:
            raise OptionError('列数少于3，不可进行分析', code = "32503003")
        if self.option('envtable').is_set:
            self.option('envtable').get_info()
            if self.option('envlabs'):
                labs = self.option('envlabs').split(',')
                for lab in labs:
                    if lab not in self.option('envtable').prop['group_scheme']:
                        raise OptionError('提供的envlabs中有不在环境因子表中存在的因子：%s', variables = (lab), code = "32503004")
            else:
                pass
            if self.option('envtable').prop['sample_number'] < 3:
                raise OptionError('环境因子表的样本数目少于3，不可进行beta多元分析', code = "32503005")
        samplelist = open(self.gettable()).readline().strip().split('\t')[1:]
        if self.option('envtable').is_set:
            self.option('envtable').get_info()
            common_samples = set(samplelist) & set(self.option('envtable').prop['sample'])
            if len(common_samples) < 3:
                raise OptionError("环境因子表和丰度表的共有样本数必需大于等于3个：%s", variables = (str(len(common_samples))), code = "32503006")
        table = open(self.gettable())
        if len(table.readlines()) < 3:
            raise OptionError('数据表特征数量必须大于等于2: %s', variables = (str(len(table.readlines()))), code = "32503007")
        table.close()
        if self.option('group_table').is_set:
            if self.option('eigenvalue') == 'row':
                for i in self.option('group_table').prop['sample_name']:
                    if i not in self.option('otutable').prop['col_sample']:
                        raise OptionError('分组文件中的样本不存在于表格中，查看是否是数据取值选择错误', code = "32503008")
            else:
                for i in self.option('group_table').prop['sample_name']:
                    if i not in self.option('otutable').prop['row_sample']:
                        raise OptionError('分组文件中的样本不存在于表格中，查看是否是数据取值选择错误', code = "32503009")
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        if self.option('group_table').is_set:
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "PCA分析结果输出目录"],
            ])
            result_dir.add_regexp_rules([
                [".+/pca_importance.xls", "xls", "主成分解释度表"],
                [".+/pca_rotation.xls", "xls", "物种主成分贡献度表"],
                [".+/pca_sites.xls", "xls", "样本坐标表"],
                [".+/.+_group.xls$", "xls", "样本坐标表"],
            ])
        else:
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "PCA分析结果输出目录"],
                ["./pca_importance.xls", "xls", "主成分解释度表"],
                ["./pca_rotation.xls", "xls", "物种主成分贡献度表"],
                ["./pca_sites.xls", "xls", "样本坐标表"],
                ["./pca_envfit_factor_scores.xls", "xls", "哑变量环境因子表"],
                ["./pca_envfit_factor.xls", "xls", "哑变量环境因子坐标表"],
                ["./pca_envfit_vector_scores.xls", "xls", "数量型环境因子表"],
                ["./pca_envfit_vector.xls", "xls", "数量型环境因子坐标表"],
            ])
        # print self.get_upload_files()
        super(PcaAgent, self).end()


class PcaTool(Tool):  # PCA需要第一行开头没有'#'的丰度表，filter_otu_sample函数生成的表头没有'#'
    def __init__(self, config):
        super(PcaTool, self).__init__(config)
        self._version = '1.0.1'  # ordination.pl脚本中指定的版本
        # self.cmd_path = os.path.join(
        #     self.config.SOFTWARE_DIR, 'bioinfo/statistical/scripts/ordination.pl')
        # self.cmd_path = 'bioinfo/statistical/scripts/ordination.pl'
        self.perl_path = "program/perl-5.24.0/bin/perl"
        self.cmd_path = self.config.PACKAGE_DIR + '/statistical/ordination.pl'
        self.script_path = "bioinfo/meta/scripts/beta_diver.sh"
        self.R_path = os.path.join(self.config.SOFTWARE_DIR, 'program/R-3.3.1/bin/R')

    def create_otu_and_env_common(self, T1, T2, new_T1, new_T2):
        T1 = pd.read_table(T1, sep='\t', dtype=str)
        T2 = pd.read_table(T2, sep='\t', dtype=str)
        T1_names = list(T1.columns[1:])
        T2_names = list(T2.iloc[0:, 0])
        T1_T2 = set(T1_names) - set(T2_names)
        T2_T1 = set(T2_names) - set(T1_names)
        T1T2 = set(T2_names) & set(T1_names)
        if len(T1T2) < 3:
            return False
        [T1_names.remove(value) for value in T1_T2]
        T1.to_csv(new_T1, sep="\t", columns=[T1.columns[0]] + T1_names, index=False)
        indexs = [T2_names.index(one) for one in T2_T1]
        T2 = T2.drop(indexs)
        T2.to_csv(new_T2, sep="\t", index=False)
        return True

    def get_otu_table(self):
        """
        根据level返回进行计算的丰度表路径
        :return:
        """
        if self.option('otutable').format == "meta.otu.tax_summary_dir":
            otu_path = self.option('otutable').get_table(self.option('level'))
        elif self.option('otutable').format == "toolapps.table":  # add by zhouxuan 20170623
            otu_path = self.option('otutable').prop['new_table']
        else:
            otu_path = self.option('otutable').prop['path']
        otufile = OtuTableFile()
        otufile.set_path(otu_path)
        # 丰度表对象没有样本列表属性
        samplelist = open(otu_path).readline().strip().split('\t')[1:]
        otufile.sub_otu_sample(samplelist, self.work_dir + "/otu_file.no_zero.xls",rmsame=self.option('scale'))
        otu_path = self.work_dir + "/otu_file.no_zero.xls"  # 去除在所有样本中丰度为0的物种 add by zouxuan 20180111
        otu_file = open(otu_path, "r")
        if len(otu_file.readlines()) < 3:
            raise OptionError('数据表特征数量必须大于等于2: %s', variables = (str(len(otu_file.readlines()))), code = "32503010")
        return otu_path
        # if self.option('envtable').is_set:
        #     return self.filter_otu_sample(otu_path, self.option('envtable').prop['sample'],
        #                                   os.path.join(self.work_dir, 'temp_filter.otutable'))
        # else:
        #     return otu_path

    """
    def filter_otu_sample(self, otu_path, filter_samples, newfile):
        if not isinstance(filter_samples, types.ListType):
            raise Exception('过滤丰度表样本的样本名称应为列表')
        try:
            with open(otu_path, 'rb') as f, open(newfile, 'wb') as w:
                one_line = f.readline()
                all_samples = one_line.rstrip().split('\t')[1:]
                if not ((set(all_samples) & set(filter_samples)) == set(filter_samples)):
                    raise Exception('提供的过滤样本存在丰度表中不存在的样本all:%s,filter_samples:%s' % (all_samples, filter_samples))
                if len(all_samples) == len(filter_samples):
                    return otu_path
                samples_index = [all_samples.index(i) + 1 for i in filter_samples]
                w.write('OTU\t' + '\t'.join(filter_samples) + '\n')
                for line in f:
                    all_values = line.rstrip().split('\t')
                    new_values = [all_values[0]] + [all_values[i] for i in samples_index]
                    w.write('\t'.join(new_values) + '\n')
                return newfile
        except IOError:
            raise Exception('无法打开丰度相关文件或者文件不存在')
    """

    def get_new_env(self):
        """
        根据envlabs生成新的envtable
        """
        if self.option('envlabs'):
            new_path = self.work_dir + '/temp_env_table.xls'
            self.option('envtable').sub_group(new_path, self.option('envlabs').split(','))
            return new_path
        else:
            return self.option('envtable').path

    def run(self):
        """
        运行
        """
        super(PcaTool, self).run()
        if self.option('group_table').is_set:
            group_de = self.option('group_table').prop['group_scheme']
            self.logger.info(group_de)
            for i in group_de:
                name = []
                name.append(i)
                target_path = os.path.join(self.work_dir, i + '_group.xls')
                self.logger.info(target_path)
                self.option('group_table').sub_group(target_path, name)
                self.run_ordination(i, target_path, 'cmd_' + i.lower(), 'pca_' + i.lower())
        else:
            self.run_ordination()
        self.end()

    def formattable(self, tablepath):
        # this_table = ''  # 感觉没有用注释掉 zhouxuan
        # with open(tablepath) as table:
        #     if table.read(1) == '#':
        #         newtable = os.path.join(self.work_dir, 'temp_format.table')
        #         with open(newtable, 'w') as w:
        #             w.write(table.read())
        #         this_table = newtable
        this_table = tablepath
        if self.option('eigenvalue') != 'row':
            newtable = this_table + '.T'
            self.t_table(this_table, newtable)
            return newtable
        else:
            return this_table

    def t_table(self, table_file, new_table):
        """
        转换颠倒表格内容
        """
        with open(table_file) as f, open(new_table, 'w') as w:
            table_list = [i.rstrip().split('\t') for i in f.readlines()]
            table_list = map(lambda *a: '\t'.join(a) + '\n', *table_list)
            w.writelines(table_list)

    def dashrepl(self, matchobj):
        """
        add func by guhaidong 20171025
        """
        return self.name_to_name[matchobj.groups()[0]]

    def add_taxon(self, old_result, taxon_result):
        """
        add func by guhaidong 20171025
        description: 将旧注释的名称，根据词典替换成新注释名称
        """
        with open(old_result, "r") as f, open(taxon_result, "w") as w:
            # w.write(old_result)
            for i in f.readlines():
                # line = i.strip()
                new_line = re.sub(r"(name\d+)", self.dashrepl, i)
                w.write(new_line)

    def run_ordination(self, group=None, group_table=None, cmd1='cmd', cmd2='pca'):
        """
        运行ordination.pl
        """
        old_otu_table = self.get_otu_table()  # 根据level返回进行计算的丰度表，特殊的format类型才会进行。小工具不用判断
        self.name_to_name = mask_taxon(old_otu_table, self.work_dir + "/tmp_mask_otu.xls")  # add by guhaidong 20171025
        old_otu_table = self.work_dir + '/tmp_mask_otu.xls'  # add by guhaidong 20171025
        if self.option('envtable').is_set:  # 小工具不用判断
            old_env_table = self.get_new_env()
            self.otu_table = self.work_dir + '/new_otu.xls'
            self.env_table = self.work_dir + '/new_env.xls'
            if not self.create_otu_and_env_common(old_otu_table, old_env_table, self.otu_table, self.env_table):
                self.set_error('环境因子表中的样本与丰度表中的样本共有数量少于2个', code = "32503011")
        else:
            self.otu_table = old_otu_table
        if group:
            target_path = os.path.join(self.work_dir, group + '_datatable.xls')
            self.option('otutable').get_table_of_main_table(old_otu_table, target_path, group_table)
            real_otu_path = target_path
        else:
            real_otu_path = self.formattable(self.otu_table)  # 获取转置文件
        cmd = self.perl_path
        cmd += ' %s -type pca -community %s -outdir %s -scale %s' % (
            self.cmd_path, real_otu_path, self.work_dir, self.option('scale'))
        if self.option('envtable').is_set:
            cmd += ' -pca_env T -environment %s' % self.env_table
        self.logger.info('运行ordination.pl程序计算pca')
        self.logger.info(cmd)
        cmd1 = self.add_command(cmd1, cmd).run()  # 命令1运行
        self.wait(cmd1)
        if cmd1.return_code == 0:
            self.logger.info("生成 cmd.r 文件成功")
        else:
            self.logger.info('生成 cmd.r 文件失败')
            self.set_error('无法生成 cmd.r 文件', code = "32503012")
        cmd_ = self.script_path + ' %s %s' % (self.R_path, self.work_dir + "/cmd.r")
        self.logger.info(cmd_)
        cmd2 = self.add_command(cmd2, cmd_).run()  # 命令2 运行
        self.wait(cmd2)
        if cmd2.return_code == 0:
            self.logger.info("pca计算成功")
        else:
            self.logger.info('pca计算失败')
            self.set_error('R程序计算pca失败', code = "32503013")
        self.logger.info('运行ordination.pl程序计算pca完成')
        if group:
            os.link(self.work_dir + "/cmd.r", self.work_dir + "/" + group + "cmd.r")
            os.remove(self.work_dir + "/cmd.r")
            group_r_path = os.path.join(self.output_dir, group)
            os.mkdir(group_r_path)
            allfiles = self.get_filesname()
            self.linkfile(self.work_dir + '/pca/' + allfiles[0], 'pca_importance.xls', group)
            self.linkfile(self.work_dir + '/pca/' + allfiles[1], 'pca_rotation.xls', group)
            self.linkfile(self.work_dir + '/pca/' + allfiles[2], 'pca_sites.xls', group)
            self.linkfile(self.work_dir + '/pca/' + allfiles[3], 'pca_rotation_all.xls', group)
            os.link(group_table, os.path.join(group_r_path, os.path.basename(group_table)))  # 分组文件link到结果目录下
            os.rename(self.work_dir + '/pca', self.work_dir + '/pca_' + group)
        else:
            allfiles = self.get_filesname()
            self.linkfile(self.work_dir + '/pca/' + allfiles[0], 'pca_importance.xls')
            self.linkfile(self.work_dir + '/pca/' + allfiles[1], 'pca_rotation.xls')
            self.linkfile(self.work_dir + '/pca/' + allfiles[2], 'pca_sites.xls')
            self.linkfile(self.work_dir + '/pca/' + allfiles[3], 'pca_rotation_all.xls')
        if self.option('envtable').is_set:
            if allfiles[4]:
                self.linkfile(self.work_dir + '/pca/' + allfiles[4], 'pca_envfit_factor_scores.xls')
                self.linkfile(self.work_dir + '/pca/' + allfiles[5], 'pca_envfit_factor.xls')
            if allfiles[6]:
                self._magnify_vector(self.work_dir + '/pca/' + allfiles[6], self.work_dir + '/pca/' + allfiles[2],
                                     self.work_dir + '/pca/' + 'pca_envfit_vector_scores_magnify.xls')
                self.linkfile(self.work_dir + '/pca/' + 'pca_envfit_vector_scores_magnify.xls',
                              'pca_envfit_vector_scores.xls')
                self.linkfile(self.work_dir + '/pca/' + allfiles[7], 'pca_envfit_vector.xls')
                # self.end()

    def linkfile(self, oldfile, newname, group=None):
        """
        link文件到output文件夹
        :param oldfile: 资源文件路径
        :param newname: 新的文件名
        :return:
        """
        if group:
            newpath = os.path.join(self.output_dir, group, newname)
        else:
            newpath = os.path.join(self.output_dir, newname)
        if os.path.exists(newpath):
            os.remove(newpath)
        os.link(oldfile, newpath)

    def _magnify_vector(self, vector_file, sites_file, new_vector):
        """
        放大环境因子向量的箭头
        """

        def get_range(path):
            with open(path) as sites:
                sites.readline()
                pc1 = []
                pc2 = []
                for line in sites:
                    split_line = line.rstrip().split('\t')
                    if len(split_line) < 3:
                        self.set_error('未知原因，坐标文件少于3列', code = "32503014")
                    pc1.append(abs(float(split_line[1])))
                    pc2.append(abs(float(split_line[2])))
                    # range_pc = min([max(pc1), max(pc2)])
            return max(pc1), max(pc2)

        range_vector_pc1, range_vector_pc2 = get_range(vector_file)
        range_sites_pc1, range_sites_pc2 = get_range(sites_file)
        magnify_1 = range_sites_pc1 / range_vector_pc1
        magnify_2 = range_sites_pc2 / range_vector_pc2
        magnify = magnify_1 if magnify_1 < magnify_2 else magnify_2
        with open(new_vector, 'w') as new, open(vector_file) as vector:
            new.write(vector.readline())
            for line in vector:
                line_split = line.rstrip().split('\t')
                for i in range(len(line_split) - 1):
                    index = i + 1
                    line_split[index] = str(float(line_split[index]) * magnify)
                new.write('\t'.join(line_split) + '\n')

    def get_filesname(self):
        """
        获取并检查文件夹下的文件是否存在

        :return pca_importance_file, pca_rotation_file,pca_rotation_all_file,
                pca_sites_file, pca_factor_score_file, pca_factor_file,
                pca_vector_score_file, pca_vector_file: 返回各个文件，以及是否存在环境因子，
                存在则返回环境因子结果
        """
        filelist = os.listdir(self.work_dir + '/pca')
        pca_dir = os.path.join(self.work_dir, 'pca')
        pca_importance_file = None
        pca_rotation_file = None
        pca_rotation_all_file = None
        pca_sites_file = None
        pca_factor_score_file = None
        pca_factor_file = None
        pca_vector_score_file = None
        pca_vector_file = None
        for name in filelist:
            if 'pca_importance.xls' in name:
                pca_importance_file = name
            elif 'pca_sites.xls' in name:
                pca_sites_file = name
            elif 'pca_rotation.xls' in name:  # modified by guhaidong 20171025
                # pca_rotation_file = name
                self.add_taxon(os.path.join(pca_dir, name), pca_dir + '/pca_rotation_new.xls')
                pca_rotation_file = 'pca_rotation_new.xls'
            elif 'pca_rotation_all.xls' in name:  # modified by guhaidong 20171025
                # pca_rotation_all_file = name
                self.add_taxon(os.path.join(pca_dir, name), pca_dir + '/pca_rotation_all_new.xls')
                pca_rotation_all_file = 'pca_rotation_all_new.xls'
            elif 'pca_envfit_factor_scores.xls' in name:
                pca_factor_score_file = name
            elif 'pca_envfit_factor.xls' in name:
                pca_factor_file = name
            elif 'pca_envfit_vector_scores.xls' in name:
                pca_vector_score_file = name
            elif 'pca_envfit_vector.xls' in name:
                pca_vector_file = name
        if pca_importance_file and pca_rotation_file and pca_sites_file:
            if self.option('envtable').is_set:
                if pca_factor_score_file:
                    if not pca_factor_file:
                        self.set_error('未知原因，环境因子相关结果丢失或者未生成,factor文件不存在', code = "32503015")
                else:
                    if pca_factor_file:
                        self.set_error('未知原因，环境因子相关结果丢失或者未生成,factor_scores文件不存在', code = "32503016")
                if pca_vector_score_file:
                    if not pca_vector_file:
                        self.set_error('未知原因，环境因子相关结果丢失或者未生成,vector文件不存在', code = "32503017")
                else:
                    if pca_vector_file:
                        self.set_error('未知原因，环境因子相关结果丢失或者未生成,vector_scores文件不存在', code = "32503018")
                    elif not pca_factor_score_file:
                        self.set_error('未知原因，环境因子相关结果全部丢失或者未生成', code = "32503019")
                return [pca_importance_file, pca_rotation_file,
                        pca_sites_file, pca_rotation_all_file, pca_factor_score_file, pca_factor_file,
                        pca_vector_score_file, pca_vector_file]

            else:
                return [pca_importance_file, pca_rotation_file, pca_sites_file, pca_rotation_all_file]
        else:
            self.set_error('未知原因，数据计算结果丢失或者未生成', code = "32503020")
