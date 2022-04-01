# -*- coding: utf-8 -*-
# __author__ = "guanqing.zou"

from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
import string
import types
import subprocess
from biocluster.core.exceptions import OptionError


class VpaAgent(Agent):
    """
    计算，需要vpa.pl
    author: guanqing.zou
    last_modifued:2018.10.10
    """

    def __init__(self, parent):
        super(VpaAgent, self).__init__(parent)
        options = [
            {"name": "species_table", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "env_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "out_name", "type": "string", "default": 'vpa'}
        ]
        self.add_option(options)
        self.step.add_steps('VpaAnalysis')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.VpaAnalysis.start()
        self.step.update()

    def step_end(self):
        self.step.VpaAnalysis.finish()
        self.step.update()

    def check_options(self):

        if not self.option('species_table').is_set:
             raise OptionError("species_table file must be provided !", code="34101701")
        if not self.option('env_table').is_set:
             raise OptionError("env_table file must be provided !", code="34101702")
        if not self.option('group_table').is_set:
             raise OptionError("group_table file must be provided !", code="34101703")

    def set_resource(self):
        """
        设置内存和CPU
        """
        self._cpu = 10
        self._memory = '15G'
        size_number = os.path.getsize(self.option("species_table").prop['path']) / (1024 * 1024)
        if int(size_number) > 1:
            self._memory = str(50 + int(size_number)*10)+"G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "VPA分析结果目录"],
            ["./env.R2adj.xls", "xls", "R2adj表"]
        ])
        #print self.get_upload_files()
        super(VpaAgent, self).end()


class VpaTool(Tool):
    def __init__(self, config):
        super(VpaTool, self).__init__(config)
        self._version = '1.0.1'
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)

    def run(self):
        """
        运行
        """
        super(VpaTool, self).run()
        self.input_check()  # by xieshichang 20200422
        self.run_vpa()
        self.convert_pdf_to()
        self.change_table_lab()
        self.set_output()
        self.end()

    def input_check(self):
        """
        by xieshichang 20200422
        判断当前所选的样本和环境因子表中，是否在某一环境因子下所有无种的数值相同
        若存在这种情况则无法进行VPA分析
        """
        cmd = '/program/Python/bin/python ' + self.config.PACKAGE_DIR + '/metagenomic/scripts/check_vap_input.py {} {}'.format(
                self.option('species_table').prop['path'], self.option('env_table').prop['path'])
        command = self.add_command('check_vap_input', cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            if os.path.exists('check_vap_input.txt'):
                with open('check_vap_input.txt', 'r') as r:
                    err_info = r.readline()
                self.logger.info(err_info)
                self.set_error('env factors:' + str(err_info) + ' have same values among all selected species, which are not supported in VPA analysis')

    def run_vpa(self):
        """
        运行vpa.pl
        """
        cmd = self.config.SOFTWARE_DIR + '/program/perl/perls/perl-5.24.0/bin/perl ' + self.config.PACKAGE_DIR + '/metagenomic/scripts/vpa.pl '
        cmd += '-spe {0} -env {1} -g {2} -o {3}'.format(self.option('species_table').prop['path'], self.option('env_table').prop['path'], self.option('group_table').prop['path'], self.option('out_name'))

        self.logger.info('开始运行vpa.pl')
        self.logger.info('运行命令：'+ cmd)


        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('生成 VPA.cmd.r 成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成 VPA.cmd.r 失败')
            self.set_error('Unable to generate cmd.r file', code="34101701")
        """
        try:
            subprocess.check_output(self.config.SOFTWARE_DIR +
                                    '/program/R-3.3.1/bin/R --restore --no-save < %s/VPA.cmd.r' % (
                                    self.work_dir), shell=True)
            self.logger.info('R计算成功')
        except subprocess.CalledProcessError:
            self.logger.info('R计算失败')
            self.set_error('R running calculation vpa failed', code="34101702")
        """
        cmd1 = '/program/R-3.3.1/bin/R --restore --no-save -f %s/VPA.cmd.r' % (self.work_dir)
        command1 = self.add_command('vpa_cmd', cmd1).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("R计算成功")
        else:
            self.set_error("R计算成功", code="34101702")
        self.logger.info('运行vpa.pl程序计算完成')

    def convert_pdf_to(self):
        self.image_magick = '/program/ImageMagick/bin/convert'
        old = self.work_dir + '/vpa.pdf'
        png = self.work_dir + '/vpa.png'
        #svg = self.work_dir + '/vpa.svg'
        m = {'png':png}
        for i in ['png']:
            cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' + old + ' ' + m[i]
            command = self.add_command('convert'+i, cmd)
            command.run()
            self.wait()
            if command.return_code == 0:
                self.logger.info('convert %s done' %i)
            else:
                self.logger.info('%s convert failed'%i)

    def change_table_lab(self):

        with open('env.plot.xls') as fr , open('env.plot.xls_new', 'w') as fw :
            line = fr.readline()
            line = line.strip()
            g = line.split(';')[1:]

            if len(g) == 2:
                use_map = {
                    'a' : g[0],
                    'b' : g[0]+'&'+g[1],
                    'c' : g[1],
                    'd' : 'Residuals'
                }

            elif len(g) ==3 :
                use_map = {
                    'a': g[0],
                    'b': g[1],
                    'c': g[2],
                    'd': g[0]+"&"+g[1],
                    'e': g[1]+"&"+g[2],
                    'f': g[0]+"&"+g[2],
                    'g': g[0]+"&"+g[1]+"&"+g[2],
                    'h': 'Residuals'
                }

            elif len(g) == 4:
                A = g[0]
                B = g[1]
                C = g[2]
                D = g[3]
                use_map = {
                    'a': A,
                    'b': B,
                    'c': C,
                    'd': D,
                    'e': A+"&"+B,
                    'f': B+"&"+C,
                    'g': A+"&"+C,
                    'h': A+"&"+D,
                    'i': B+"&"+D,
                    'j': C+"&"+D,
                    'k': A+"&"+B+"&"+D,
                    'l': A+"&"+B+"&"+C,
                    'm': B+"&"+C+"&"+D,
                    'n': A+"&"+C+"&"+D,
                    'o': A+"&"+B+"&"+C+"&"+D,
                    'p': 'Residuals'
                }
            else:
                self.set_error("分组数必须是2-4组", code="34101703")

            line_2 = fr.readline()
            fw.write(line_2)
            for line in fr:
                spline = line.split('\t',1)
                lab = spline[0].split(' ')[0]
                lab_2 = lab[1:2]
                new_lab = use_map[lab_2]
                fw.write(new_lab+'\t'+ spline[1])

        os.rename('env.plot.xls','env.plot.xls_old')
        os.rename('env.plot.xls_new','env.plot.xls')



    def set_output(self):
        files = ['env.R2adj.xls','vpa.png','vpa.pdf']
        for file in files:
            if os.path.exists(self.output_dir + '/' +file):
                os.remove(self.output_dir + '/' +file)
            os.link(self.work_dir + '/' + file,self.output_dir + '/' +file)

