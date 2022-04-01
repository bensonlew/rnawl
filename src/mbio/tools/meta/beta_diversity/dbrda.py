# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
from biocluster.agent import Agent
from biocluster.tool import Tool

import types
import os
import pandas as pd
import re  # add re by guhaidong 20171031
from mbio.packages.taxon.mask_taxon import mask_taxon,mask_env  # add by guhaidong 20171031 ,zhujuan 20180110
from biocluster.core.exceptions import OptionError
from mbio.files.meta.otu.otu_table import OtuTableFile
from mbio.packages.beta_diversity.dbrda_r import *
from itertools import islice
import subprocess  # by houshuang 20190924


class DbrdaAgent(Agent):
    """
    dbrda_r.py
    version v1.0
    author: shenghe
    last_modified:2016.03.24
    """
    # METHOD = ["manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "morisita", "horn",
    #           "mountford", "raup", "binomial", "chao"]
    METHOD = ["manhattan", "euclidean", "canberra", "bray", "kulczynski", "gower"]
    METHOD_DICT = {"manhattan": "manhattan", "euclidean": "euclidean", "canberra": "canberra",
                   "bray_curtis": "bray", "kulczynski": "kulczynski", "gower": "gower"}

    def __init__(self, parent):
        super(DbrdaAgent, self).__init__(parent)
        options = [
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "method", "type": "string", "default": "bray_curtis",'required':True},
            {"name": "level", "type": "string", "default": "otu"},
            {"name": "dis_matrix", "type": "infile", "format": "meta.beta_diversity.distance_matrix"},
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "envlabs", "type": "string", "default": ""},  # 用逗号分隔的环境因子名称
            {"name": "group_table", "type": "infile", "format": "toolapps.group_table, meta.otu.group_table"},  # by houshuang 20190924
            {"name": "ellipse", "type": "string", "default": "F"}  # by houshuang 20190924
        ]
        self.add_option(options)
        self.step.add_steps('dbRDA')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.dbRDA.start()
        self.step.update()

    def step_end(self):
        self.step.dbRDA.finish()
        self.step.update()

    def gettable(self):
        """
        根据level返回进行计算的丰度表对象，否则直接返回参数otutable对象
        :return:
        """
        if self.option('otutable').format == "meta.otu.tax_summary_dir":
            new_otu = OtuTableFile()
            new_otu.set_path(self.option('otutable').get_table(self.option('level')))
            new_otu.get_info()
            return new_otu
        else:
            return self.option('otutable')

    def check_options(self):
        """
        重写参数检查
        """
        if self.option('envtable').is_set:
            with open(self.option('envtable').prop['path']) as f:  # add by zhujuan 20171115 解决环境因子文件有误的bug
                for line in islice(f, 1, None):
                   line = line.strip('\n').split('\t')
                   for i in range(1, len(line)):
                      if re.match("(-{0,1})[0-9.]+$", line[i]):
                          pass
                      elif re.match("^[a-zA-Z]+[0-9]*$", line[i]):
                          pass
                      else:
                          raise OptionError('环境因子数据中存在非法字符,如含有空格、下划线等,请检查环境因子文件!', code="32701701")
            self.option('envtable').get_info()
            if self.option('envlabs'):
                labs = self.option('envlabs').split(',')
                for lab in labs:
                    if lab not in self.option('envtable').prop['group_scheme']:
                        raise OptionError('提供的envlabs中有不在环境因子表中存在的因子：%s', variables=(lab), code="32701702")
            else:
                pass
            if self.option('envtable').prop['sample_number'] < 3:
                raise OptionError('环境因子表的样本数目少于3，不可进行beta多元分析', code="32701703")
            if self.option('dis_matrix').is_set:
                if not self.option('otutable').is_set:
                    raise OptionError("必须提供丰度表格", code="32701704")
                self.option('dis_matrix').get_info()
                env_collection = set(self.option('envtable').prop['sample'])
                collection = set(self.option('dis_matrix').prop['samp_list']) & env_collection
                if len(collection) < 3:
                    raise OptionError("环境因子表和丰度表的共有样本数必需大于等于3个：%s", variables=(len(collection)), code="32701705")
                pass
                if self.option('envlabs'):  # add by zhujuan 20171123 当环境因子数小于样品数时，才计算
                    if len(collection) <= len(self.option('envlabs').split(',')):
                        raise OptionError("环境因子个数必须小于样品个数，请剔除多余的环境因子!", code="32701706")
                else:
                    envlist = open(self.option('envtable').prop['path']).readline().strip().split('\t')[1:]
                    if len(collection) <= len(envlist):
                        raise OptionError("环境因子个数必须小于样品个数，请剔除多余的环境因子!", code="32701707")
            else:
                if self.option('method') not in DbrdaAgent.METHOD_DICT:
                    raise OptionError('错误或者不支持的距离计算方法', code="32701708")
                self.option('method', DbrdaAgent.METHOD_DICT[self.option('method')])
                if not self.option('otutable').is_set:
                    raise OptionError('没有提供距离矩阵的情况下，必须提供丰度表', code="32701709")
                self.real_otu = self.gettable()
                if self.real_otu.prop['sample_num'] < 3:
                    raise OptionError('丰度表的样本数目少于3，不可进行beta多元分析', code="32701710")
                samplelist = open(self.real_otu.path).readline().strip().split('\t')[1:]
                common_samples = set(samplelist) & set(self.option('envtable').prop['sample'])
                if len(common_samples) < 3:
                    raise OptionError("环境因子表和丰度表的共有样本数必需大于等于3个：%s", variables=(len(common_samples)), code="32701711")
                if self.option('envlabs'):  # add by zhujuan 20171123 当环境因子数小于样品数时，才计算
                    if len(common_samples) <= len(self.option('envlabs').split(',')):
                        raise OptionError("环境因子个数必须小于样品个数，请剔除多余的环境因子!", code="32701712")
                else:
                    envlist = open(self.option('envtable').prop['path']).readline().strip().split('\t')[1:]
                    if len(common_samples) <= len(envlist):
                        raise OptionError("环境因子个数必须小于样品个数，请剔除多余的环境因子!", code="32701703")
                table = open(self.real_otu.path)
                if len(table.readlines()) < 4:
                    raise OptionError('提供的数据表信息少于3行', code="32701714")
                table.close()
                pass
        else:
            raise OptionError('没有提供环境因子表', code="32701715")
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '10G'  # 增加内存 顾海东 20171101 ##fix by qingchen.zhang@20200817

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "db_rda分析结果目录"],
            ["./db_rda_sites.xls", "xls", "db_rda样本坐标表"],
            ["./db_rda_importance.xls", "xls", "db_rda主成分解释度"],  # add by zhujuan 2017.08.21
            ["./db_rda_plot_species_data.xls", "xls", "db_rda物种坐标表"],
            ["./db_rda_centroids.xls", "xls", "db_rda哑变量环境因子坐标表"],
            ["./db_rda_biplot.xls", "xls", "db_rda数量型环境因子坐标表"],
        ])
        # print self.get_upload_files()
        super(DbrdaAgent, self).end()


class DbrdaTool(Tool):
    def __init__(self, config):
        super(DbrdaTool, self).__init__(config)
        self._version = '1.0'
        # 模块脚本路径，并不使用
        self.cmd_path = 'mbio/packages/beta_diversity/dbrda_r.py'
        self.env_table = self.get_new_env()
        if not self.option('dis_matrix').is_set:
            self.otu_table = self.get_otu_table()
            self.name_to_name = mask_taxon(self.otu_table, self.work_dir + "/tmp_mask_otu.xls")  # add by guhaidong 20171031
            self.otu_table = self.work_dir + '/tmp_mask_otu.xls'  # add by guhaidong
            new_otu_table = self.work_dir + '/new_otu_table.xls'
            new_env_table = self.work_dir + '/new_env_table.xls'
            if not self.create_otu_and_env_common(self.otu_table, self.env_table, new_otu_table, new_env_table):
                self.set_error('环境因子表中的样本与丰度表中的样本共有数量少于2个', code="32701701")
            else:
                self.otu_table = new_otu_table
                self.env_table = new_env_table
        else:
            self.otu_table = self.get_otu_table()
            self.name_to_name = mask_taxon(self.otu_table, self.work_dir + "/tmp_mask_otu.xls")  # add by guhaidong 20171031
            self.otu_table = self.work_dir + '/tmp_mask_otu.xls'  # add by guhaidong
            new_otu_table = self.work_dir + '/new_otu_table.xls'
            new_env_table = self.work_dir + '/new_env_table.xls'
            samples = list(
                set(self.option('dis_matrix').prop['samp_list']) & set(self.option('envtable').prop['sample']))
            self.env_table = self.sub_env(samples)
            self.dis_matrix = self.get_matrix(samples)
            if not self.create_otu_and_env_common(self.otu_table, self.env_table, new_otu_table, new_env_table):
                self.set_error('环境因子表中的样本与丰度表中的样本共有数量少于2个', code="32701702")
            else:
                self.otu_table = new_otu_table
        self.env_to_name = mask_env(self.env_table, self.work_dir + '/tmp_mask_env.xls')  # add by zhujuan
        self.env_table = self.work_dir + '/tmp_mask_env.xls'

    def sub_env(self, samples):
        with open(self.env_table) as f, open(self.work_dir + '/sub_env_temp.xls', 'w') as w:
            w.write(f.readline())
            for i in f:
                if i.split('\t')[0] in samples:
                    w.write(i)
        return self.work_dir + '/sub_env_temp.xls'

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

    def get_matrix(self, samples):
        if len(self.option('dis_matrix').prop['samp_list']) == len(self.option('envtable').prop['sample']):
            return self.option('dis_matrix').path
        else:
            # samples = list(set(self.option('dis_matrix').prop['samp_list']) & set(self.option('envtable').prop['sample']))
            self.option('dis_matrix').create_new(samples,
                                                 os.path.join(self.work_dir, 'dis_matrix_filter.temp'))
            return os.path.join(self.work_dir, 'dis_matrix_filter.temp')

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

    def get_otu_table(self):
        """
        根据level返回进行计算的丰度表路径
        :return:
        """
        if self.option('otutable').format == "meta.otu.tax_summary_dir":
            otu_path = self.option('otutable').get_table(self.option('level'))
        else:
            otu_path = self.option('otutable').prop['path']
        # 丰度表对象没有样本列表属性
        return otu_path
        # return self.filter_otu_sample(otu_path, self.option('envtable').prop['sample'],
        #                               os.path.join(self.work_dir + 'temp_filter.otutable'))

    def filter_otu_sample(self, otu_path, filter_samples, newfile):
        if not isinstance(filter_samples, types.ListType):
            self.set_error('OTU表过滤时应提供样本名列表', code="32701703")
        try:
            with open(otu_path, 'rb') as f, open(newfile, 'wb') as w:
                one_line = f.readline()
                all_samples = one_line.rstrip().split('\t')[1:]
                if not ((set(all_samples) & set(filter_samples)) == set(filter_samples)):
                    self.set_error('提供的过滤样本存在丰度表中不存在的样本all:%s,filter_samples:%s', variables=(all_samples, filter_samples), code="32701704")
                if len(all_samples) == len(filter_samples):
                    return otu_path
                samples_index = [all_samples.index(i) + 1 for i in filter_samples]
                w.write('#OTU\t' + '\t'.join(filter_samples) + '\n')
                for line in f:
                    all_values = line.rstrip().split('\t')
                    new_values = [all_values[0]] + [all_values[i] for i in samples_index]
                    w.write('\t'.join(new_values) + '\n')
                return newfile
        except IOError:
            self.set_error('无法打开丰度相关文件或者文件不存在', code="32701705")

    def linkfile(self, oldfile, newname):
        """
        link文件到output文件夹
        :param oldfile: 资源文件路径
        :param newname: 新的文件名
        :return:
        """
        newpath = os.path.join(self.output_dir, newname)
        if os.path.exists(newpath):
            os.remove(newpath)
        os.link(oldfile, newpath)

    def dashrepl(self, matchobj):
        """
        add func by guhaidong 20171031
        """
        return self.name_to_name[matchobj.groups()[0]]

    def dashrepl_env(self, matchobj):
        """
        add by zhujuan 20180110
        将环境因子名称替换成name[0-9]*进行运算，后将结果中的环境因子名改回真实名，避免有些+等特殊字符不能运算
        """
        return self.env_to_name[matchobj.groups()[0]]

    def add_taxon(self, old_result, taxon_result, type):
        """
        add func by guhaidong 20171025, last modify by 20180110
        description: 将旧注释的名称，根据词典替换成新注释名称
        """
        dashrepl = ""
        if type == "tax":
            dashrepl = self.dashrepl
        elif type == "env":
            dashrepl = self.dashrepl_env
        with open(old_result, "r") as f, open(taxon_result, "w") as w:
            # w.write(old_result)
            for i in f.readlines():
                #line = i.strip()
                new_line = re.sub(r"(name\d+)", dashrepl, i)
                w.write(new_line)

    def get_euclidify_params(self, dis_method):
        sqrt_dist = 'F'
        add_dist = 'F'
        euclidi_dis = ['euclidean', 'binary_euclidean','manhattan']
        if dis_method not in euclidi_dis:
            sqrt_dist = 'T'
            #add_dist = 'T'

        return sqrt_dist,add_dist

    def run_dbrda(self):
        """
        运行dbrda.py
        """
        sqrt_dist,add_dist = self.get_euclidify_params(dis_method=self.option('method'))
        self.logger.info('运行dbrda_r.py程序计算Dbrda')
        if self.option('dis_matrix').is_set:
            return_mess = db_rda_dist(dis_matrix=self.dis_matrix, env=self.env_table, species_table=self.otu_table,
                                      output_dir=self.work_dir, sqrt_dist=sqrt_dist, add_dist=add_dist)
        else:
            return_mess = db_rda_new(self.otu_table, self.env_table, self.work_dir,
                                     self.option('method'), sqrt_dist=sqrt_dist, add_dist=add_dist)
        # self.logger.info('运行dbrda_r.py程序计算Dbrda成功')
        if return_mess == 0:
            otu_species_list = self.get_species_name()
            self.linkfile(self.work_dir + '/db_rda_sites.xls', 'db_rda_sites.xls')
            self.linkfile(self.work_dir + '/db_rda_cont.xls', 'db_rda_importance.xls')  ##add by zhujuan 20170821
            self.add_taxon(self.work_dir + '/db_rda_envfit.xls', self.work_dir + '/db_rda_envfit_new.xls', "env")
            self.linkfile(self.work_dir + '/db_rda_envfit_new.xls', 'db_rda_envfit.xls')
            if self.option('dis_matrix').is_set:
                self.add_taxon(self.work_dir + '/db_rda_species.xls', self.work_dir + '/species_new.xls', "tax")  # add by guhaidong 20171031
                self.linkfile(self.work_dir + '/species_new.xls', 'db_rda_species.xls')  # add by guhaidong 20171031
                # self.linkfile(self.work_dir + '/db_rda_species.xls', 'db_rda_species.xls')  # mv by guhaidong
                if len(otu_species_list) == 0:  # add by zhujuan 2017.10.10 for get plot_species.xls
                    # self.linkfile(self.work_dir + '/db_rda_plot_species_data.xls', 'db_rda_plot_species_data.xls')  # mv by guhaidong
                    self.add_taxon(self.work_dir + '/db_rda_plot_species_data.xls', self.work_dir + '/plot_species_data_new.xls', "tax")  # add by guhaidong
                    self.linkfile(self.work_dir + '/plot_species_data_new.xls', 'db_rda_plot_species_data.xls')  # add by guhaidong
                else:
                    new_file_path = self.get_new_species_xls(otu_species_list)
                    self.logger.info(otu_species_list)
                    self.add_taxon(new_file_path, self.work_dir + '/plot_species_data_new.xls' , "tax")  # add by guhaidong
                    self.linkfile(self.work_dir + '/plot_species_data_new.xls', 'db_rda_plot_species_data.xls')  # add by guhaidong
                    # self.linkfile(new_file_path, 'db_rda_plot_species_data.xls')  # mv by guhaidong
            lines = open(self.work_dir + '/env_data.temp').readlines()
            if 'centroids:TRUE' in lines[0]:
                self.add_taxon(self.work_dir + '/db_rda_centroids.xls', self.work_dir + '/db_rda_centroids_new.xls', "env")
                self.linkfile(self.work_dir + '/db_rda_centroids_new.xls', 'db_rda_centroids.xls')
            if 'biplot:TRUE' in lines[1]:
                self.add_taxon(self.work_dir + '/db_rda_biplot.xls', self.work_dir + '/db_rda_biplot_new.xls', "env")
                self.linkfile(self.work_dir + '/db_rda_biplot_new.xls', 'db_rda_biplot.xls')
            self.logger.info('运行dbrda_r.py程序计算Dbrda完成')
            # by houshuang 20190923 增加分组椭圆 >>>
            if self.option("ellipse") == 'T':
                sites_file = os.path.join(self.output_dir, 'db_rda_sites.xls')
                if os.path.isfile(sites_file):
                    cmd = self.config.PACKAGE_DIR + '/statistical/ordination2.pl'
                    cmd += ' -type dbrda -community %s -outdir %s -group %s -sites %s' % (
                        self.otu_table, self.work_dir, self.option('group_table').path, sites_file)
                    self.logger.info(cmd)
                    try:
                        subprocess.check_output(cmd, shell=True)
                        self.logger.info('生成 cmd2.r 文件成功')
                    except subprocess.CalledProcessError:
                        self.logger.info('生成 cmd2.r 文件失败')
                        self.set_error('无法生成 cmd2.r 文件', code="32701707")
                    try:
                        subprocess.check_output(self.config.SOFTWARE_DIR +
                                                '/program/R-3.3.1/bin/R --restore --no-save < %s/cmd2.r' % self.work_dir,
                                                shell=True)
                        self.logger.info('group_ellipse计算成功')
                        self.linkfile(self.work_dir + '/dbrda/ellipse.xls', "ellipse.xls")
                    except subprocess.CalledProcessError:
                        self.logger.info('group_ellipse计算失败')
                        self.set_error('R程序计算group_ellipse失败', code="32701708")
                else:
                    self.logger.info('引用的sites文件不存在')
            # <<<
            self.end()
        else:
            self.set_error('运行dbrda_r.py程序计算Dbrda出错', code="32701706")

    def get_species_name(self):  # modify by zhujuan 1017.10.09
        """
        判断丰度表中的物种数量是否大于30 ，如果大于30，筛选出丰度在前30的物种
        :return: 丰度为前30的物种或者 空的列表
        """
        old_abund_file_path = self.otu_table  # modified by guhaidong 20171102
        df = pd.DataFrame(pd.read_table(old_abund_file_path, sep='\t', index_col=0))
        df['Col_sum'] = df.apply(lambda x: x.sum(), axis=1)
        new_otu_file = df.sort_values(by=['Col_sum'], ascending=0).head(30)
        species_list = list(new_otu_file.index)
        return species_list

    def get_new_species_xls(self, otu_species_list):  # add by zhujuan 1017.10.09
        """
        根据物种列表信息，获取新的species表格文件
        :param otu_species_list: 物种列表信息
        :return: 新的species文件的路径
        """
        if os.path.exists(self.work_dir + "/db_rda_species.xls"):
            old_species_table = self.work_dir + "/db_rda_species.xls"
            new_species_table = self.work_dir + "/db_rda_plot_species_data.xls"
        with open(old_species_table, "rb") as table, open(new_species_table, "a") as w:
            line = table.readlines()
            w.write(line[0])
            for l in line[1:]:
                content = l.strip().split("\t")
                if content[0] in otu_species_list:
                    w.write('\t'.join(content) + "\n")
        return new_species_table

    def run(self):
        """
        运行
        """
        super(DbrdaTool, self).run()
        self.run_dbrda()
