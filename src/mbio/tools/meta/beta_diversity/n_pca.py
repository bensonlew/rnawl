#-*- coding: utf-8 -*-
# __author__ = "zhangpeng"

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import string
import types
import subprocess
from biocluster.core.exceptions import OptionError

class NPcaAgent(Agent):
    """
    需要 .r
    version v1.1
    author: zhangpeng
    last_modifued:2016.12.19
    """

    def __init__(self, parent):
        super(NPcaAgent, self).__init__(parent)
        options = [
            #{"name": "mode", "type": "int", "default":1},
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir", "default":None},
            #{"name": "factor_table", "type": "string", "default":""},
            {"name": "level", "type": "int", "default": 9},
            {"name": "second_group_table", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            #{"name": "method", "type": "string", "default":"sum"},
            #{"name": "name_table", "type": "string", "default":""},
            #{"name": "top_n", "type": "int", "default": 100},
            #{"name": "env_table", "type": "infile", "format": "meta.otu.group_table"},
            #{"name": "env_labs", "type": "string", "default": ""},
            #{"name": "PCAlabs", "type": "string", "default": ""}
            ]
        self.add_option(options)
        self.step.add_steps('NPcaAnalysis')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.NPcaAnalysis.start()
        self.step.update()

    def step_end(self):
        self.step.NPcaAnalysis.finish()
        self.step.update()

    def gettable(self):
        """
        根据输入的otu表和分类水平计算新的otu表
        :return:
        """
        if self.option('otu_table').format == "meta.otu.tax_summary_dir":
            return self.option('otu_table').get_table(self.option('level'))
        else:
            return self.option('otu_table').prop['path']

    def check_options(self):
        if not self.option('otu_table') and self.option('mode') in [1]:
            raise OptionError('必须提供OTU表', code="32702701")
        if self.option('second_group_table').is_set:
            self.option('second_group_table').get_info()
        return True
        if self.option('group_table').is_set:
            self.option('group_table').is_info()
        return True

    
    def set_resource(self):
        """
        设置内存和CPU
        """
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
                [".", "", "npca分析结果目录"],
                ["./sites.xls", "xls", "坐标数据"],
                ["./sd.xls", "xls", "方差大小"],
                ["./sdmax.xls", "xls", "置信上边界"],
                ["./sdmin.xls", "xls", "置信下边界"],
                ["./rotation_mean.xls", "xls", "平均值"],
                #["./linesdata.xls", "xls", "直线信息"],
                #["./sitesall.xls", "xls", "所有信息"],
                ["./importance.xls", "xls", "百分率返回"],
                ["./sitesall.xls", "xls", "所有信息"]
                ])
        print self.get_upload_files()
        super(NPcaAgent, self).end()

class NPcaTool(Tool):
    def __init__(self, config):
        super(NPcaTool, self).__init__(config)
        self._version = '1.0.1'
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)

    def run(self):
        """
        运行
        """
        super(NPcaTool, self).run()
        self.run_npca_perl()

    def get_otu_table(self):
        """
        根据调用的level参数重构otu表
        :return:
        """
        if self.option('otu_table').format == "meta.otu.tax_summary_dir":
            otu_path = self.option('otu_table').get_table(self.option('level'))
        else:
            otu_path = self.option('otu_table').prop['path']
        return otu_path
            
            

                               
    def run_npca_perl(self):
        """
        运行calc_environmentalregression.perl
        """
        old_otu_table = self.get_otu_table()
        self.otu_table = self.work_dir + '/new_otu.xls'      
        cmd = self.config.SOFTWARE_DIR + '/miniconda2/bin/perl ' + self.config.SOFTWARE_DIR + '/bioinfo/meta/scripts/n-pca.pl '
        cmd += '-i %s ' %(self.get_otu_table())
        cmd += '-o %s ' %(self.work_dir + '/npca/')
        if not os.path.exists(self.work_dir + '/npca/'):
            os.mkdir(self.work_dir + '/npca/') 
        cmd += '-m %s ' %(self.option('second_group_table').prop['path'])
        cmd += '-g %s ' %(self.option('group_table').prop['path'])
        self.logger.info('开始运行n-pca.pl计算N维PCA相关数据')
        
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('生成 cmd.r 成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成 cmd.r 失败')
            self.set_error('无法生成 cmd.r 文件', code="32702701")
        try:
            subprocess.check_output(self.config.SOFTWARE_DIR +
                                    '/program/R-3.3.1/bin/R --restore --no-save < %s/cmd.r' % (self.work_dir + '/npca'), shell=True)
            self.logger.info('npca计算成功')
        except subprocess.CalledProcessError:
            self.logger.info('npca计算失败')
            self.set_error('R运行计算npca失败', code="32702702")
            raise "运行R脚本计算npca相关数据失败"
        self.logger.info('运行n-pca.pl程序进行npca计算完成')
        allfiles = self.get_npca_filesname()
        self.linkfile(self.work_dir + '/npca/' + allfiles[0], 'sites.xls')
        self.linkfile(self.work_dir + '/npca/' + allfiles[1], 'sitesall.xls')
        self.linkfile(self.work_dir + '/npca/' + allfiles[2], 'sdmax.xls')
        self.linkfile(self.work_dir + '/npca/' + allfiles[3], 'sdmin.xls')
        self.linkfile(self.work_dir + '/npca/' + allfiles[4], 'rotation_mean.xls')
        self.linkfile(self.work_dir + '/npca/' + allfiles[5], 'importance.xls')
        self.linkfile(self.work_dir + '/npca/' + allfiles[6], 'sd.xls')
        self.end()

    def linkfile(self, oldfile, newname):
        newpath = os.path.join(self.output_dir, newname)
        if os.path.exists(newpath):
            os.remove(newpath)
        os.link(oldfile, newpath)

    def get_npca_filesname(self):
        filelist = os.listdir(self.work_dir + '/npca')
        sites = None
        sd = None
        sdmax = None
        sdmin = None
        rotation_mean = None
        importance = None
        sitesall = None
        for name in filelist:
            if 'sites.xls' in name:
                sites = name
            elif 'sitesall.xls' in name:
                sitesall = name
            elif 'sdmax.xls' in name:
                sdmax = name
            elif 'sdmin.xls' in name:
                sdmin = name
            elif 'rotation_mean.xls' in name:
                rotation_mean = name
            elif 'importance.xls' in name:
                importance = name
            elif 'sd.xls' in name:
                sd = name
        if (sites and sitesall and sdmax and sdmin and rotation_mean and importance and sd):
            return [sites, sitesall, sdmax, sdmin, rotation_mean, importance, sd]
        else:
            self.set_error("未知原因，metagenomeSeq计算结果丢失", code="32702703")
         
