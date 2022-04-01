# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import re


class PcoaAgent(Agent):
    """
    PCoA
    """
    def __init__(self, parent):
        super(PcoaAgent, self).__init__(parent)
        options = [
            {"name": "dis_matrix", "type": "infile", "format": "meta.beta_diversity.distance_matrix"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "ellipse", "type": "string", "default": "F"},
        ]
        self.add_option(options)
        self.step.add_steps('pcoa')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.pcoa.start()
        self.step.update()

    def step_end(self):
        self.step.pcoa.finish()
        self.step.update()

    def check_options(self):
        if not self.option('dis_matrix').is_set:
            raise OptionError('必须提供输入距离矩阵表')
        if self.option('ellipse') not in ['T', 'F']:
            raise OptionError('ellipse必须为T或者F')
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
        super(PcoaAgent, self).end()


class PcoaTool(Tool):
    def __init__(self, config):
        super(PcoaTool, self).__init__(config)
        self._version = "1.0"
        self.program = {
            'perl': "program/perl-5.24.0/bin/perl",
            'r_path': os.path.join(self.config.SOFTWARE_DIR, 'program/R-3.3.1/bin/R'),
        }
        self.script = {
            'ordination': os.path.join(self.config.PACKAGE_DIR, 'statistical/ordination.pl'),
            'beta_div': "bioinfo/meta/scripts/beta_diver.sh"
        }

    def run_ordination(self):
        """
        运行ordination.pl
        """
        cmd = self.program['perl'] + ' ' + self.script['ordination']
        if self.option('group').is_set and not self.option('group').prop['is_empty'] and self.option('ellipse')=='T':
            group_run = self.option('group').prop['path']
            cmd += ' -type pcoa -dist %s -outdir %s -group %s' % (
                self.option('dis_matrix').prop['path'], self.work_dir, group_run)
        else:
            cmd += ' -type pcoa -dist %s -outdir %s' % (
                self.option('dis_matrix').prop['path'], self.work_dir)
        self.logger.info('运行ordination.pl程序计算pcoa')
        self.logger.info(cmd)
        cmd1 = self.add_command("cmd.r", cmd).run()
        self.wait(cmd1)
        if cmd1.return_code == 0:
            self.logger.info("生成 cmd.r 文件成功")
        else:
            self.logger.info('生成 cmd.r 文件失败')
            self.set_error('无法生成 cmd.r 文件')
        cmd_ = self.script['beta_div'] + ' %s %s' % (self.program['r_path'], self.work_dir + "/cmd.r")
        self.logger.info(cmd_)
        cmd2 = self.add_command("pcoa", cmd_).run()
        self.wait(cmd2)
        if cmd2.return_code == 0:
            self.logger.info("pcoa计算成功")
        else:
            self.logger.info('pcoa计算失败')
            self.set_error('R程序计算pcoa失败', code="32703002")
        allfile = self.get_filesname()
        sites_file = self.format_header(allfile[1])
        eigenvalues = self.format_eigenvalues(allfile[0], self.work_dir + '/format_eigenvalues.txt')
        eigenvaluespre = self.format_eigenvalues(allfile[2], self.work_dir + '/format_eigenvaluespre.txt')
        self.linkfile(eigenvalues, 'pcoa_eigenvalues.xls')
        self.linkfile(eigenvaluespre, 'pcoa_eigenvaluespre.xls')
        self.linkfile(sites_file, 'pcoa_sites.xls')
        if self.option('group').is_set and not self.option('group').prop['is_empty'] and self.option('ellipse')=='T':
            self.linkfile(self.work_dir + '/pcoa/ellipse.xls', 'ellipse.xls')
        self.logger.info('运行ordination.pl程序计算pcoa完成')

    def format_eigenvalues(self, fp, new_fp):
        with open(fp) as f, open(new_fp, 'w') as w:
            w.write(f.readline())
            for i in f:
                w.write('PC' + i)
        return new_fp

    def format_header(self, old):
        """
        """
        with open(old) as f, open(self.work_dir + '/format_header.temp', 'w') as w:
            headers = f.readline().rstrip().split('\t')[1:]
            for header in headers:
                if header[0] != 'V':
                    self.set_error('Pcoa结果不正确或者不规范', code="32703003")
                num = header[1:]
                if not num.isdigit():
                    self.set_error('Pcoa结果不正确或者不规范', code="32703003")
            news = ['PC' + i[1:] for i in headers]
            news = [''] + news
            new_header = '\t'.join(news) + '\n'
            w.write(new_header)
            for i in f:
                w.write(i)
        return self.work_dir + '/format_header.temp'

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

    def get_filesname(self):
        """
        获取并检查文件夹下的文件是否存在

        :return pcoa_sites_file,pcoa_eigenvalues_file: 返回各个文件
        """
        filelist = os.listdir(self.work_dir + '/pcoa')
        pcoa_eigenvalues_file = None
        pcoa_eigenvaluespre_file = None
        pcoa_sites_file = None
        for name in filelist:
            if name.endswith('pcoa_eigenvalues.xls'):
                pcoa_eigenvalues_file = name
            elif name.endswith('pcoa_sites.xls'):
                pcoa_sites_file = name
            elif name.endswith('pcoa_eigenvaluespre.xls'):
                pcoa_eigenvaluespre_file = name
        if pcoa_eigenvalues_file and pcoa_sites_file:
            return [self.work_dir + '/pcoa/' + pcoa_eigenvalues_file, self.work_dir + '/pcoa/' + pcoa_sites_file,
                    self.work_dir + '/pcoa/' + pcoa_eigenvaluespre_file]
        else:
            self.set_error('未知原因，数据计算结果丢失或者未生成', code="32703004")

    def run(self):
        super(PcoaTool, self).run()
        self.run_ordination()
        self.end()

