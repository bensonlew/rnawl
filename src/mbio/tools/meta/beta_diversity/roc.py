# -*- coding: utf-8 -*-
# __author__ = "zhangpeng"

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import string
import types
import subprocess
from biocluster.core.exceptions import OptionError

class RocAgent(Agent):
    """
    需要calc_roc.pl
    version v1.1
    author: JieYao
    last_modifued:2016.08.22
    """

    def __init__(self, parent):
        super(RocAgent, self).__init__(parent)
        options = [
            {"name": "mode", "type": "int", "default":1},
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir", "default":None},
            #{"name": "factor_table", "type": "string", "default":""},
            {"name": "level", "type": "int", "default": 9},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "method", "type": "string", "default":"sum"},
            #{"name": "name_table", "type": "string", "default":""},
            {"name": "top_n", "type": "int", "default": 100}
            ]
        self.add_option(options)
        self.step.add_steps('RocAnalysis')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.RocAnalysis.start()
        self.step.update()

    def step_end(self):
        self.step.RocAnalysis.finish()
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
            raise OptionError('必须提供OTU表', code="32703601")
        #if not self.option('factor_table') and self.option('mode') == 3:
            #raise OptionError('模式三必须提供factor表')
        if self.option('mode') in [1,2]:
            self.option('otu_table').get_info()
            if self.option('otu_table').prop['sample_num'] < 2:
                raise OptionError('otu表的样本数目小于2，不可进行ROC计算', code="32703602")
            if not self.option('group_table').is_set:
                raise OptionError("必须提供分组表格", code="32703603")
            if self.option('method'):
                if self.option('method') not in ['sum', 'average', 'median']:
                    raise OptionError("丰度计算方法只能选择sum,average,median之一", code="32703604")
                if self.option('mode') == 2:
                    if not os.path.exists(self.option('name_table')):
                        raise OptionError("Mode 2 模式下必须提供物种名列表文件", code="32703605")
                #os.system('cat %s | awk -F "\t" \'{ print $1 }\' > tmp.txt' %(self.gettable()))
                #otu_data = open("tmp.txt", "r").readlines()[1:]
                #os.remove('tmp.txt')
                #otu_data = map(string.rstrip, otu_data)
                #sample_data = open(self.gettable()).readline().strip().split()[1:]
                #group_data = open(self.option('group_table').prop['path']).readlines()[1:]
                #for s in group_data:
                    i#f s.split()[0] not in sample_data:
                        #raise OptionError("分组表中物种%s不在OTU Table中" % s.split()[0])
                    #if s.split()[1] not in ['0','1']:
                        #raise OptionError("分组表中物种分组只能有0和1！")

                if self.option('mode')==2:
                    name_data = open(self.option('name_table'), "r").readlines()[1:]
                    name_data = map(string.strip, name_data)
                    for s in name_data:
                        if s not in otu_data:
                            raise OptionError("物种%s不在OTU Table中", variables=(s), code="32703606")

                #if self.option('mode')==1:
                    #if self.option('top_n')>len(otu_data):
                        #raise OptionError("选择丰度前N高物种时，设定的N多于物种总数：%d>%d" %(self.option('top_n'), len(otu_data)))
        if self.option('mode') == 3:
            os.system('cat %s | awk -F "\t" \'{ print $1 }\' > tmp.txt' %(self.option('factor_table')))
            sample_data = open("tmp.txt", "r").readlines()[1:]
            os.remove('tmp.txt')
            sample_data = map(string.rstrip, sample_data)
            group_data = open(self.option('group_table').prop['path']).readlines()[1:]
            for s in group_data:
                if s.split()[0] not in sample_data:
                    raise OptionError("分组表中物种%s不在Factor Table中", variables=(s.split()[0]), code="32703607")
                if s.split()[1] not in ['0','1']:
                    raise OptionError("分组表中物种分组只能有0和1！", code="32703608")
        return True

    
    def set_resource(self):
        """
        设置内存和CPU
        """
        self._cpu = 10
        self._memory = '15G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
                [".", "", "ROC分析结果目录"],
                ["./roc_curve.xls", "xls", "ROC受试者工作特征曲线数据"],
                ["./roc_auc.xls", "xls", "ROC受试者工作特征曲线-AUC VALUE"]
                ])
        print self.get_upload_files()
        super(RocAgent, self).end()

class RocTool(Tool):
    def __init__(self, config):
        super(RocTool, self).__init__(config)
        self._version = '1.0.1'
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)

    def run(self):
        """
        运行
        """
        super(RocTool, self).run()
        self.run_roc_perl()

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
                               
    def run_roc_perl(self):
        """
        运行calc_roc.perl
        """
        cmd = self.config.SOFTWARE_DIR + '/miniconda2/bin/perl ' + self.config.SOFTWARE_DIR + '/bioinfo/meta/scripts/roc_plus_gq.pl '
        cmd += '-o %s ' %(self.work_dir + '/ROC/')
        if not os.path.exists(self.work_dir + '/ROC/'):
            os.mkdir(self.work_dir + '/ROC/')
        if self.option('mode') in [1,2]:
            cmd += '-i %s ' %(self.get_otu_table())
        else:
            cmd += '-i %s ' %(self.option('factor_table'))
        cmd += '-mode %d ' %(self.option('mode'))
        cmd += '-group %s ' %(self.option('group_table').prop['path'])
        if self.option('mode')==2:
            cmd += '-name %s ' %(self.option('name_table'))
        if self.option('method'):
            cmd += '-method %s ' %(self.option('method'))
        if self.option('mode')==1:
            cmd += '-n %d ' %(self.option('top_n'))
        cmd += ' -ci F -smooth T'
        """
        with open(self.work_dir + "/ROC/set_env_and_run.cmd", "w") as tmp_file:
            tmp_file.write("#!/usr/bash\n\n")
            tmp_file.write("export PATH=/mnt/ilustre/users/sanger-dev/app/gcc/5.1.0/bin:$PATH")
            tmp_file.write("export LD_LIBRARY_PATH=/mnt/ilustre/users/sanger-dev/app/gcc/5.1.0/lib64:$LD_LIBRARY_PATH")
            tmp_file.write(cmd + "\n")
        cmd = "bash %s/ROC/set_env_and_run.cmd" % self.work_dir
        """
        self.logger.info('开始运行calc_roc.pl计算ROC相关数据')
        
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('生成 roc.cmd.r 成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成 roc.cmd.r 失败')
            self.set_error('无法生成 roc.cmd.r 文件', code="32703601")
        try:
            subprocess.check_output(self.config.SOFTWARE_DIR +
                                    '/program/R-3.3.1/bin/R --restore --no-save < %s/roc.cmd.r' % (self.work_dir + '/ROC'), shell=True)
            self.logger.info('ROC计算成功')
        except subprocess.CalledProcessError:
            self.logger.info('ROC计算失败')
            self.set_error('R运行计算ROC失败', code="32703602")
            raise "运行R脚本计算ROC相关数据失败"
        self.logger.info('运行calc_roc.pl程序进行ROC计算完成')
        allfiles = self.get_roc_filesname()
        self.linkfile(self.work_dir + '/ROC/' + allfiles[0], 'roc_curve.xls')
        self.linkfile(self.work_dir + '/ROC/' + allfiles[1], 'roc_auc.xls')
        self.linkfile(self.work_dir + '/ROC/' + allfiles[2], 'roc_plot_rocarea.xls') # guanqing.zou
        self.linkfile(self.work_dir + '/ROC/' + allfiles[3], 'best_loc.xls')
        self.end()

    def linkfile(self, oldfile, newname):
        newpath = os.path.join(self.output_dir, newname)
        if os.path.exists(newpath):
            os.remove(newpath)
        if "roc_curve" in oldfile:
            data = open(oldfile).readlines()[1:]
            with open(oldfile, "w") as tmp_file:
                tmp_file.write("x\ty\tgroup\n")
                for s in data:
                    tmp_file.write(s)
        if "aucvalue" in oldfile:
            data = open(oldfile).readlines()
            with open(oldfile, "w") as tmp_file:
                for i in data:
                    s = i.strip().split()
                    tmp_file.write(s[1]+"\t"+s[2]+"\n")
        os.link(oldfile, newpath)

    def get_roc_filesname(self):
        filelist = os.listdir(self.work_dir + '/ROC')
        roc_curve = None
        roc_auc = None
        plot_rocarea = None
        best_loc = None
        for name in filelist:
            if 'roc_curve.xls' in name:
                roc_curve = name
            elif 'roc_auc.xls' in name:
                roc_auc = name
            elif 'roc_plot_rocarea.xls' in name:    # guanqing.zou
                plot_rocarea = name
            elif 'best_loc.xls' in name:
                best_loc = name
        if (roc_curve and roc_auc and plot_rocarea and best_loc):
            return [roc_curve, roc_auc, plot_rocarea, best_loc]
        else:
            self.set_error("未知原因，ROC计算结果丢失", code="32703603")
    
