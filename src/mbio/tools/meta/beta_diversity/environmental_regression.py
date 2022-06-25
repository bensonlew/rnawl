# -*- coding: utf-8 -*-
# __author__ = "zhangpeng"

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import string
import types
import subprocess
from biocluster.core.exceptions import OptionError

class EnvironmentalRegressionAgent(Agent):
    """
    需要 .r
    version v1.1
    author: zhangpeng
    last_modifued:2016.12.09
    """

    def __init__(self, parent):
        super(EnvironmentalRegressionAgent, self).__init__(parent)
        options = [
            #{"name": "mode", "type": "int", "default":1},
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir", "default":None},
            #{"name": "factor_table", "type": "string", "default":""},
            {"name": "level", "type": "int", "default": 9},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            #{"name": "method", "type": "string", "default":"sum"},
            #{"name": "name_table", "type": "string", "default":""},
            #{"name": "top_n", "type": "int", "default": 100},
            {"name": "env_table", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "env_labs", "type": "string", "default": ""},
            {"name": "PCAlabs", "type": "string", "default": ""}
            ]
        self.add_option(options)
        self.step.add_steps('EnvironmentalRegressionAnalysis')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.EnvironmentalRegressionAnalysis.start()
        self.step.update()

    def step_end(self):
        self.step.EnvironmentalRegressionAnalysis.finish()
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
            raise OptionError('必须提供OTU表', code="32702101")
        #if not self.option('factor_table') and self.option('mode') == 3:
            #raise OptionError('模式三必须提供factor表')
        if self.option('env_table').is_set:
            self.option('env_table').get_info()
            if self.option('env_labs'):
                labs = self.option('env_labs').split(',')
                if len(labs) > 2:
                    raise OptionError("变量只能有一个", code="32702102")
                #if labs not in self.option('env_table').prop['group_scheme']:
                    #raise OptionError('提供的envlabs中有不在环境因子表中存在的因子：%s' %labs)
            else:
                raise OptionError('必须提供因子', code="32702103")
        else:
            raise OptionError('必须提供环境因子表格', code="32702104")
        if self.option('group_table').is_set:
            self.option('group_table').get_info()
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
                [".", "", "environmentalregression分析结果目录"],
                ["./sites.xls", "xls", "metagenomeseq坐标数据"],
                ["./rotation.xls", "xls", "pca的全部输出"],
                ["./rotationone.xls", "xls", "一个pca的输出"],
                ["./messages.xls", "xls", "主要信息"],
                #["./linesdata.xls", "xls", "直线信息"],
                ["./importance.xls", "xls", "百分率返回"]
                ])
        print self.get_upload_files()
        super(EnvironmentalRegressionAgent, self).end()

class EnvironmentalRegressionTool(Tool):
    def __init__(self, config):
        super(EnvironmentalRegressionTool, self).__init__(config)
        self._version = '1.0.1'
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)

    def run(self):
        """
        运行
        """
        super(EnvironmentalRegressionTool, self).run()
        self.run_environmentalregression_perl()

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
        
    def get_new_env(self):
        """
        根据envlabs生成新的envtable
        """
        if self.option('env_labs'):
            new_path = self.work_dir + '/temp_env_table.xls'
            self.option('env_table').sub_group(new_path, self.option('env_labs').split(','))
            return new_path
        else:
            return self.option('env_table').path
            
            
    def create_otu_and_env_common(self, T1, T2, new_T1, new_T2):
        import pandas as pd
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
                               
    def run_environmentalregression_perl(self):
        """
        运行calc_environmentalregression.perl
        """
        old_env_table = self.get_new_env()
        #print (self.env_table)
        old_otu_table = self.get_otu_table()
        self.otu_table = self.work_dir + '/new_otu.xls'
        self.env_table = self.work_dir + '/temp_env_table.xls'
        #if not self.create_otu_and_env_common(old_otu_table, old_env_table, self.otu_table, self.env_table):
            #self.set_error('环境因子表中的样本与OTU表中的样本共有数量少于2个')       
        cmd = self.config.SOFTWARE_DIR + '/miniconda2/bin/perl ' + self.config.SOFTWARE_DIR + '/bioinfo/meta/scripts/plot-pcoa.pl '
        cmd += '-i %s ' %(self.get_otu_table())
        cmd += '-o %s ' %(self.work_dir + '/environmentalregression/')
        if not os.path.exists(self.work_dir + '/environmentalregression/'):
            os.mkdir(self.work_dir + '/environmentalregression/')
        #cmd += '-i %s ' %(self.get_otu_table())
        #cmd += '-factor %s ' %(self.env_table) 
        cmd += '-fname %s ' %(self.option('env_labs'))
        cmd += '-number %s ' %(self.option('PCAlabs'))
        cmd += '-factor %s ' %(self.env_table)
        #cmd += '-group %s ' %(self.option('group_table').prop['path'])
        """
        with open(self.work_dir + "/environmentalregression/set_env_and_run.cmd", "w") as tmp_file:
            tmp_file.write("#!/usr/bash\n\n")
            tmp_file.write("export PATH=/mnt/ilustre/users/sanger-dev/app/gcc/5.1.0/bin:$PATH")
            tmp_file.write("export LD_LIBRARY_PATH=/mnt/ilustre/users/sanger-dev/app/gcc/5.1.0/lib64:$LD_LIBRARY_PATH")
            tmp_file.write(cmd + "\n")
        cmd = "bash %s/ROC/set_env_and_run.cmd" % self.work_dir
        """
        self.logger.info('开始运行plot-pcoa.pl计算environmentalregression相关数据')
        
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('生成 cmd.r 成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成 cmd.r 失败')
            self.set_error('无法生成 cmd.r 文件', code="32702101")
        try:
            subprocess.check_output(self.config.SOFTWARE_DIR +
                                    '/program/R-3.3.1/bin/R --restore --no-save < %s/cmd.r' % (self.work_dir + '/environmentalregression'), shell=True)
            self.logger.info('environmentalregression计算成功')
        except subprocess.CalledProcessError:
            self.logger.info('environmentalregression计算失败')
            self.set_error('R运行计算environmentalregression失败', code="32702102")
            raise "运行R脚本计算environmentalregression相关数据失败"
        self.logger.info('运行environmentalregression.pl程序进行environmentalregression计算完成')
        allfiles = self.get_environmentalregression_filesname()
        self.linkfile(self.work_dir + '/environmentalregression/' + allfiles[0], 'sites.xls')
        self.linkfile(self.work_dir + '/environmentalregression/' + allfiles[1], 'rotation.xls')
        self.linkfile(self.work_dir + '/environmentalregression/' + allfiles[2], 'rotationone.xls')
        self.linkfile(self.work_dir + '/environmentalregression/' + allfiles[3],  'messages.xls')
        #self.linkfile(self.work_dir + '/environmentalregression/' +  allfiles[4], 'linesdata.xls')
        self.linkfile(self.work_dir + '/environmentalregression/' + allfiles[4], 'importance.xls')
        self.end()

    def linkfile(self, oldfile, newname):
        newpath = os.path.join(self.output_dir, newname)
        if os.path.exists(newpath):
            os.remove(newpath)
        if "sites" in oldfile:
            data = open(oldfile).readlines()[1:]
            with open(oldfile, "w") as tmp_file:
                tmp_file.write("name\tPCA\tfactor\n")
                for s in data:
                    tmp_file.write(s)
        if "rotation" in oldfile:
            data = open(oldfile).readlines()
            with open(oldfile, "w") as tmp_file:
                tmp_file.write("PCA\t")
                for s in data:
                    tmp_file.write(s)
        os.link(oldfile, newpath)

    def get_environmentalregression_filesname(self):
        filelist = os.listdir(self.work_dir + '/environmentalregression')
        sites = None
        rotation = None
        rotationone = None
        messages = None
        #linesdata = None
        importance = None
        for name in filelist:
            if 'sites.xls' in name:
                sites = name
            elif 'rotation.xls' in name:
                rotation = name
            elif 'rotationone.xls' in name:
                rotationone = name
            elif 'messages.xls' in name:
                messages = name
            #elif 'linesdata.xls' in name:
                #linesdata = name
            elif 'importance.xls' in name:
                importance = name
        if (sites and rotation and rotationone and messages and importance):
            return [sites, rotation, rotationone, messages, importance]
        else:
            self.set_error("未知原因，metagenomeSeq计算结果丢失", code="32702103")
         
