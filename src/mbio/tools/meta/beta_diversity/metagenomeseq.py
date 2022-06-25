# -*- coding: utf-8 -*-
# __author__ = "zhangpeng"

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import string
import types
import subprocess
from biocluster.core.exceptions import OptionError

class MetagenomeseqAgent(Agent):
    """
    需要 metagenomeSeq.r
    version v1.1
    author: zhangpeng
    last_modifued:2016.08.22
    """

    def __init__(self, parent):
        super(MetagenomeseqAgent, self).__init__(parent)
        options = [
            #{"name": "mode", "type": "int", "default":1},
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir", "default":None},
            #{"name": "factor_table", "type": "string", "default":""},
            {"name": "level", "type": "int", "default": 9},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            #{"name": "method", "type": "string", "default":"sum"},
            #{"name": "name_table", "type": "string", "default":""},
            #{"name": "top_n", "type": "int", "default": 100}
            ]
        self.add_option(options)
        self.step.add_steps('MetagenomeseqAnalysis')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.MetagenomeseqAnalysis.start()
        self.step.update()

    def step_end(self):
        self.step.MetagenomeseqAnalysis.finish()
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
            raise OptionError('必须提供OTU表', code="32702401")
        #if not self.option('factor_table') and self.option('mode') == 3:
            #raise OptionError('模式三必须提供factor表')
        if self.option('group_table').is_set:
            self.option('group_table').get_info()
            if len(self.option('group_table').prop['sample']) < 2:
                raise OptionError('分组表的样本数目少于2，不可进行特征分析', code="32702402")
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
                [".", "", "Metagenomeseq分析结果目录"],
                ["./diff.xls", "xls", "metagenomeseq数据"],
                ["./list.xls", "xls", "metagenomeseq表头"]
                ])
        print self.get_upload_files()
        super(MetagenomeseqAgent, self).end()

class MetagenomeseqTool(Tool):
    def __init__(self, config):
        super(MetagenomeseqTool, self).__init__(config)
        self._version = '1.0.1'
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)

    def run(self):
        """
        运行
        """
        super(MetagenomeseqTool, self).run()
        self.run_metagenomeseq_perl()

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
                               
    def run_metagenomeseq_perl(self):
        """
        运行calc_metagenomeseq.perl
        """
        cmd = self.config.SOFTWARE_DIR + '/miniconda2/bin/perl ' + self.config.SOFTWARE_DIR + '/bioinfo/meta/scripts/metagenomeSeq.pl '
        cmd += '-o %s ' %(self.work_dir + '/metagenomeseq/')
        if not os.path.exists(self.work_dir + '/metagenomeseq/'):
            os.mkdir(self.work_dir + '/metagenomeseq/')
        cmd += '-i %s ' %(self.get_otu_table())
        cmd += '-group %s ' %(self.option('group_table').prop['path'])
        """
        with open(self.work_dir + "/ROC/set_env_and_run.cmd", "w") as tmp_file:
            tmp_file.write("#!/usr/bash\n\n")
            tmp_file.write("export PATH=/mnt/ilustre/users/sanger-dev/app/gcc/5.1.0/bin:$PATH")
            tmp_file.write("export LD_LIBRARY_PATH=/mnt/ilustre/users/sanger-dev/app/gcc/5.1.0/lib64:$LD_LIBRARY_PATH")
            tmp_file.write(cmd + "\n")
        cmd = "bash %s/ROC/set_env_and_run.cmd" % self.work_dir
        """
        self.logger.info('开始运行metagenomeSeq.pl计算metagenomeSeq相关数据')
        
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('生成 cmd.r 成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成 cmd.r 失败')
            self.set_error('无法生成 cmd.r 文件', code="32702401")
        try:
            subprocess.check_output(self.config.SOFTWARE_DIR +
                                    '/program/R-3.3.1/bin/R --restore --no-save < %s/cmd.r' % (self.work_dir + '/metagenomeseq'), shell=True)
            self.logger.info('metagenomeSeq计算成功')
        except subprocess.CalledProcessError:
            self.logger.info('metagenomeSeq计算失败')
            self.set_error('R运行计算metagenomeSeq失败', code="32702402")
            raise "运行R脚本计算metagenomeSeq相关数据失败"
        self.logger.info('运行metagenomeSeq.pl程序进行metagenomeseq计算完成')
        allfiles = self.get_metagenomeseq_filesname()
        self.linkfile(self.work_dir + '/metagenomeseq/' + allfiles[0], 'diff.xls')
        self.linkfile(self.work_dir + '/metagenomeseq/' + allfiles[1], 'list.xls')
        self.end()

    def linkfile(self, oldfile, newname):
        newpath = os.path.join(self.output_dir, newname)
        if os.path.exists(newpath):
            os.remove(newpath)
        if "diff" in oldfile:
            data = open(oldfile).readlines()[1:]
            with open(oldfile, "w") as tmp_file:
                tmp_file.write("taxa\tazero\tasum\tbzero\tbsum\toddsRatio\tlower\tupper\tfidherP\tadjp\tP.Value\tadj.P.Val\n")
                for s in data:
                    tmp_file.write(s)
        if "list" in oldfile:
            data = open(oldfile).readlines()
            with open(oldfile, "w") as tmp_file:
                tmp_file.write("sample_a_name\tsample_b_name\tcount_a_name\tcount_b_name\n")
                for s in data:
                    tmp_file.write(s)
                    #for i in data:
                    #s = i.strip().split()
                    #tmp_file.write(s[1]+"\t"+s[2]+"\n")
        os.link(oldfile, newpath)

    def get_metagenomeseq_filesname(self):
        filelist = os.listdir(self.work_dir + '/metagenomeseq')
        diff = None
        list = None
        for name in filelist:
            if 'diff.xls' in name:
                diff = name
            elif 'list.xls' in name:
                list = name
        if (diff and list):
            return [diff,list]
        else:
            self.set_error("未知原因，metagenomeSeq计算结果丢失", code="32702403")
    
