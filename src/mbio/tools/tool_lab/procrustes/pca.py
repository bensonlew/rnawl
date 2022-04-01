# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.files.meta.otu.otu_table import OtuTableFile
import os
import unittest
import re


class PcaAgent(Agent):
    """
    PCoA
    """
    def __init__(self, parent):
        super(PcaAgent, self).__init__(parent)
        options = [
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir, toolapps.table"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "scale", "type": "string", "default": "F"},
            {"name": "ellipse", "type": "string", "default": "T"},
        ]
        self.add_option(options)
        self.step.add_steps('pca')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.pca.start()
        self.step.update()

    def step_end(self):
        self.step.pca.finish()
        self.step.update()

    def check_options(self):
        if not self.option('otutable').is_set:
            raise OptionError('必须提供数据表', code="32702903")
        with open(self.option('otutable').prop['path'], 'r') as table:
            samplelist = table.readline().strip().split('\t')[1:]
            table.readline()
            if not table.readline():
                raise OptionError('数据表特征数量必须大于等于2: %s', variables=(len(table.readlines())), code="32702908")
        if len(samplelist) < 3:
            raise OptionError('列数少于3，不可进行分析', code="32702904")
        return True

    def set_resource(self):
        """
        运行所需资源
        vcftools一般不耗费内存资源
        vcftools --plink 较耗资源
        """
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        super(PcaAgent, self).end()


class PcaTool(Tool):
    def __init__(self, config):
        super(PcaTool, self).__init__(config)
        self._version = "1.0"
        self.program = {
            'perl': "program/perl-5.24.0/bin/perl",
            'r_path': os.path.join(self.config.SOFTWARE_DIR, 'program/R-3.3.1/bin/R'),
        }
        self.script = {
            'ordination': os.path.join(self.config.PACKAGE_DIR, 'statistical/ordination.pl'),
            'beta_div': "bioinfo/meta/scripts/beta_diver.sh"
        }

    def run_ordination(self, group=None, group_table=None, cmd1='cmd', cmd2='pca'):
        """
        运行ordination.pl
        """
        cmd = self.program['perl'] + ' ' + self.script['ordination']
        if self.option('group_table').is_set and not self.option('group_table').prop['is_empty'] and self.option('ellipse') == 'T':
            group_run = self.option('group_table').prop['path']
            cmd += ' -type pca -community %s -outdir %s -scale %s -group %s' % (self.option('otutable').prop['path'],
                                                                                self.work_dir, self.option('scale'),
                                                                                group_run)
        else:
            cmd += ' -type pca -community %s -outdir %s -scale %s' % (self.option('otutable').prop['path'],
                                                                      self.work_dir, self.option('scale'))
        self.logger.info('运行ordination.pl程序计算pca')
        self.logger.info(cmd)
        cmd1 = self.add_command(cmd1, cmd).run()
        self.wait(cmd1)
        if cmd1.return_code == 0:
            self.logger.info("生成 cmd.r 文件成功")
        else:
            self.logger.info('生成 cmd.r 文件失败')
            self.set_error('无法生成 cmd.r 文件')
        cmd_ = self.script['beta_div'] + ' %s %s' % (self.program['r_path'], self.work_dir + "/cmd.r")
        self.logger.info(cmd_)
        cmd2 = self.add_command(cmd2, cmd_).run()
        self.wait(cmd2)
        if cmd2.return_code == 0:
            self.logger.info("pca计算成功")
        else:
            self.logger.info('pca计算失败')
            self.set_error('R程序计算pca失败', code="32703002")
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
            try:
                self.linkfile(self.work_dir + '/pca/ellipse.xls', 'ellipse.xls')
            except Exception:
                self.set_error('分组椭圆文件丢失', code="32702904")
        else:
            allfiles = self.get_filesname()
            self.linkfile(self.work_dir + '/pca/' + allfiles[0], 'pca_importance.xls')
            self.linkfile(self.work_dir + '/pca/' + allfiles[1], 'pca_rotation.xls')
            self.linkfile(self.work_dir + '/pca/' + allfiles[2], 'pca_sites.xls')
            self.linkfile(self.work_dir + '/pca/' + allfiles[3], 'pca_rotation_all.xls')
        self.logger.info('运行ordination.pl程序计算pca完成')

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
        for name in filelist:
            if 'pca_importance.xls' in name:
                pca_importance_file = name
            elif 'pca_sites.xls' in name:
                pca_sites_file = name
            elif 'pca_rotation.xls' in name:  # modified by guhaidong 20171025
                pca_rotation_file = name
                # self.add_taxon(os.path.join(pca_dir, name), pca_dir + '/pca_rotation_new.xls')
                # pca_rotation_file = 'pca_rotation_new.xls'
            elif 'pca_rotation_all.xls' in name:  # modified by guhaidong 20171025
                pca_rotation_all_file = name
                # self.add_taxon(os.path.join(pca_dir, name), pca_dir + '/pca_rotation_all_new.xls')
                # pca_rotation_all_file = 'pca_rotation_all_new.xls'
        if pca_importance_file and pca_rotation_file and pca_sites_file:
            return [pca_importance_file, pca_rotation_file, pca_sites_file, pca_rotation_all_file]
        else:
            self.set_error('未知原因，数据计算结果丢失或者未生成', code="32702911")

    def run(self):
        super(PcaTool, self).run()
        if self.option('group_table') and self.option('group_table').format == "toolapps.group_table":
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

