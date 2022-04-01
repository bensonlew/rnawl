#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
import glob
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.alpha_diversity.group_file_split import group_file_spilt
import pandas as pd


class AlphaDiversityModule(Module):
    """
    alpha多样性模块
    version 1.0
    author: qindanhua
    last_modify: 2015.12.29
    """
    ESTIMATORS_E = ['ace', 'bergerparker', 'boneh', 'bootstrap', 'bstick', 'chao', 'coverage', 'default', 'efron',
                    'geometric', 'goodscoverage', 'heip', 'invsimpson', 'jack', 'logseries', 'npshannon', 'nseqs',
                    'qstat', 'shannon', 'shannoneven', 'shen', 'simpson', 'simpsoneven', 'smithwilson', 'sobs', 'solow']
    ESTIMATORS_R = ['ace', 'bootstrap', 'chao', 'coverage', 'default', 'heip', 'invsimpson', 'jack', 'npshannon',
                    'nseqs', 'shannon', 'shannoneven', 'simpson', 'simpsoneven', 'smithwilson', 'sobs']

    def __init__(self, work_id):
        super(AlphaDiversityModule, self).__init__(work_id)
        options = [
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table,meta.otu.tax_summary_dir"},  # 输入文件
            {"name": "estimate_indices", "type": "string", "default": "ace,chao,shannon,simpson,coverage"},
            {"name": "rarefy_indices", "type": "string", "default": "sobs,ace,chao,shannon,simpson,coverage"},  # 指数类型
            {"name": "rarefy_freq", "type": "int", "default": 100},
            {"name": "level", "type": "string", "default": "otu"},  # level水平
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "meta_group_name", "type": "string"},
        ]
        self.add_option(options)
        # self.rank_path = '/mnt/ilustre/users/sanger/app/meta/scripts/'
        self.perl_path = 'Perl/bin/perl'
        self.estimators = self.add_tool('meta.alpha_diversity.estimators')
        self.rarefaction = self.add_module('meta.alpha_diversity.rarefaction')
        self.est_t_test = self.add_tool('statistical.metastat')
        self.step.add_steps('estimators', 'rarefaction','est_t_test')
        self.group_status = 0

    def check_options(self):
        """
        检查参数
        """
        if not self.option("otu_table").is_set:
            raise OptionError("请选择otu表", code="21202701")
        for estimators in self.option('estimate_indices').split(','):
            if estimators not in self.ESTIMATORS_E:
                raise OptionError("请选择正确的指数类型", code="21202702")
        for estimators in self.option('rarefy_indices').split(','):
            if estimators != "":
                if estimators not in self.ESTIMATORS_R:
                    raise OptionError("请选择正确的指数类型", code="21202703")

    def estimators_run(self):
        self.estimators.set_options({
            'otu_table': self.option('otu_table'),
            'indices': self.option('estimate_indices'),
            'level': self.option('level')
            })
        # self.on_rely(estimators, self.rarefaction_run)
        self.step.estimators.start()
        self.estimators.on("end", self.finish_update)
        self.estimators.on("end", self.run_est_t)
        self.estimators.run()
        # self.on_rely(self.estimators, self.finish_update)

    def finish_update(self):
        self.step.estimators.finish()
        self.step.update()

    def rare_finish_update(self):
        self.step.rarefaction.finish()
        self.step.update()

    def est_t_finish_update(self):
        self.step.est_t_test.finish()
        self.step.update()

    def rarefaction_run(self):
        self.logger.info('start rarefaction_run')
        self.rarefaction.set_options({
            'otu_table': self.option('otu_table'),
            'indices': self.option('rarefy_indices'),
            'freq': self.option('rarefy_freq'),
            'level': self.option('level')
            })
        # self.rarefaction.on('end', self.set_output)
        self.step.rarefaction.start()
        self.rarefaction.on("end", self.rare_finish_update)
        self.rarefaction.run()
        self.logger.info('end rarefaction_run')
        # self.on_rely(self.rarefaction, self.finish_update)

    def run_est_t(self):
        if self.group_status:
            self.group_file_dir = self.work_dir + '/two_group_output'
            self.estimators_file = self.estimators.work_dir + '/estimators.xls'
            self.estimators_file2 = self.work_dir + "estimators_ttest.xls"
            group_name = group_file_spilt(self.option('group').prop['path'], self.group_file_dir)
            df = pd.read_table(self.estimators.work_dir+"/estimators.xls",sep='\t',header=0,index_col=0)
            df = df.T
            del_list = []
            for x in df.index.values:
                if "_lci" in x or "_hci" in x:
                    del_list.append(x)
            df = df.drop(del_list)
            df.to_csv(self.estimators_file2,sep="\t")

            name_list = []
            for g in group_name:
                if g[0] > g[1]:
                    gg = g[1] + '|' + g[0]
                    name_list.append(gg)
                else:
                    gg = g[0] + '|' + g[1]
                    name_list.append(gg)
            self.group_name = ",".join(name_list)
            self.logger.info(self.group_name)
            self.logger.info(self.option('est_test_method'))
            options = {
                'est_input': self.estimators_file2,
                'test': 'estimator',
                'est_group': self.group_file_dir,
                'est_test_method': 'mann'
            }
            self.est_t_test.set_options(options)
            self.est_t_test.on('end', self.est_t_finish_update)
            self.est_t_test.run()

    def get_group(self):

        if self.option('group').is_set:
            with open(self.option('group').prop['path'],"r") as f:
                data = f.readlines()
                group_data = {}
                for i in data[1:]:
                    if i.strip().split("\t")[1] in group_data:
                        group_data[i.strip().split("\t")[1]].append(i.strip().split("\t")[0])
                    else:
                        group_data[i.strip().split("\t")[1]] = [i.strip().split("\t")[0]]
            if len(group_data) > 1:
                self.group_status = 1
                for gp in group_data:
                    if len(group_data[gp]) < 2:
                        self.group_status = 0

    def set_output(self):
        self.logger.info('set output')
        for root, dirs, files in os.walk(self.output_dir):
            for names in dirs:
                shutil.rmtree(os.path.join(self.output_dir, names))
            for f in files:
                os.remove(os.path.join(self.output_dir, f))
        if self.option('meta_group_name'):
            if not os.path.exists(self.output_dir + "/" + self.option('meta_group_name')):
                os.mkdir(self.output_dir + "/" + self.option('meta_group_name'))
            esti_path = self.output_dir + "/" + self.option('meta_group_name') + "/Estimators"
            rare_path = self.output_dir + "/" + self.option('meta_group_name') + "/Rarefaction"
            est_t_path = self.output_dir + "/" + self.option('meta_group_name') + "/EstTTest"
        else:
            esti_path = self.output_dir + "/Estimators"
            rare_path = self.output_dir + "/Rarefaction"
            est_t_path = self.output_dir + "/EstTTest"
        if not os.path.exists(rare_path):
            os.mkdir(rare_path)
        if not os.path.exists(esti_path):
            os.mkdir(esti_path)
        # estimators = self.work_dir + '/Estimators/output/estimators.xls'
        estimators = os.path.join(self.estimators.output_dir, 'estimators.xls')  # 防止estimators因为自动重运行而取错路径 by ghd @ 20181012
        os.link(estimators, esti_path + '/estimators.xls')
        for single in glob.glob(self.estimators.work_dir + "/*.summary"):
            # print(single)
            if os.path.exists(single):
                os.link(single, esti_path + '/' + os.path.basename(single))
        # os.system('cp -r %s %s' % (rarefaction, self.output_dir))
        for estimators in self.option('rarefy_indices').split(','):
            # est_path = self.work_dir + '/Rarefaction/output/%s/' % estimators
            if estimators != "":
                est_path = os.path.join(self.rarefaction.output_dir, estimators)  # 防止rarefaction因为重运行而取错路径 by ghd @ 20190124
                os.system('cp -r %s %s' % (est_path, rare_path))
        # self.option('estimators').set_path(self.output_dir+'/estimators')
        # self.option('rarefaction').set_path(self.output_dir+'/rarefaction')
        if self.group_status:
            os.system('cp -r %s %s' % (self.est_t_test.output_dir, est_t_path))
        self.logger.info('done')
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["./Estimators/estimators.xls", "xls", "alpha多样性指数表"],
            ["./Rarefaction", "文件夹", "稀释性曲线输出目录"]
        ])
        # for i in self.option("rarefy_indices").split(","):
        #     # self.logger.info(i)
        #     if i == "sobs":
        #         result_dir.add_relpath_rules([
        #             ["./sobs", "文件夹", "{}指数结果输出目录".format(i)]
        #         ])
        #         result_dir.add_regexp_rules([
        #             [r".*rarefaction\.xls", "xls", "{}指数的simpleID的稀释性曲线表".format(i)]
        #         ])
        #         # self.logger.info("{}指数的simpleID的稀释性曲线表".format(i))
        #     else:
        #         result_dir.add_relpath_rules([
        #             ["./{}".format(i), "文件夹", "{}指数结果输出目录".format(i)]
        #         ])
        #         result_dir.add_regexp_rules([
        #             [r".*{}\.xls".format(i), "xls", "{}指数的simpleID的稀释性曲线表".format(i)]
        #         ])
        # print self.get_upload_files()
        self.end()

    def run(self):
        super(AlphaDiversityModule, self).run()
        self.estimators_run()
        self.rarefaction_run()
        self.get_group()
        if not self.group_status:
            self.on_rely([self.estimators, self.rarefaction], self.set_output)
        else:
            self.on_rely([self.estimators, self.rarefaction, self.est_t_test], self.set_output)
