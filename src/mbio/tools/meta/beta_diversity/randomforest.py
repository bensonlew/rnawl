# -*- coding: utf-8 -*-
# __author__ = 'JieYao'
# last modify by shaohua.yuan in 20171103
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import types
import subprocess
from biocluster.core.exceptions import OptionError
from mbio.packages.statistical.normalization import Normalization


class RandomforestAgent(Agent):
    """
    需要RandomForest.pl
    version v1.0
    author: JieYao
    last_modified:2016.07.18
    更新：RandomForest_CV_AUC.pl
    2018.04.23 by zhujuan 新增十折交叉验证、AUC验证和预测样本分类功能
    """

    def __init__(self, parent):
        super(RandomforestAgent, self).__init__(parent)
        options = [
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "level", "type": "int", "default": 9},
            {"name": "grouptable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "ntree", "type": "int", "default": 500},
            {"name": "problem_type", "type": "int", "default": 1},
            {"name": "predict_sample", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "method", "type": "string", "default": "CV"},
            {"name": "norm_method", "type": "string", "default": ""},
        ]
        self.add_option(options)
        self.step.add_steps('RandomforestAnalysis')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.RandomforestAnalysis.start()
        self.step.update()

    def step_end(self):
        self.step.RandomforestAnalysis.finish()
        self.step.update()

    def gettable(self):
        """
        根据输入的otu表和分类水平计算新的otu表
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
        if not self.option('otutable').is_set:
            raise OptionError('必须提供otu表', code="32703401")
        self.option('otutable').get_info()
        if self.option('otutable').prop['sample_num'] < 2:
            raise OptionError('otu表的样本数目少于2，不可进行随机森林特征分析', code="32703402")
        if self.option('grouptable').is_set:
            self.option('grouptable').get_info()
            if len(self.option('grouptable').prop['sample']) < 2:
                raise OptionError('分组表的样本数目少于2，不可进行随机森林特征分析', code="32703403")
        samplelist = open(self.gettable()).readline().strip().split('\t')[1:]
        if self.option('grouptable').is_set:
            self.option('grouptable').get_info()
            if len(self.option('grouptable').prop['sample']) > len(samplelist):
                raise OptionError('OTU表中的样本数量:%s少于分组表中的样本数量:%s', variables=(len(samplelist),
                                                                   len(self.option('grouptable').prop['sample'])), code="32703404")
            for sample in self.option('grouptable').prop['sample']:
                if sample not in samplelist:
                    raise OptionError('分组表的样本中存在OTU表中未知的样本%s', variables=(sample), code="32703405")
        table = open(self.gettable())
        if len(table.readlines()) < 4:
            raise OptionError('数据表信息少于3行', code="32703406")
        table.close()
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(RandomforestAgent, self).end()


class RandomforestTool(Tool):
    def __init__(self, config):
        super(RandomforestTool, self).__init__(config)
        self._version = '1.0.1'
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.cmd_path = self.config.PACKAGE_DIR + '/meta/scripts/RandomForest_CV_AUC.pl'
        if self.option('grouptable').is_set:
            self.group_table = self.option('grouptable').prop['path']
        self.otu_table = self.get_otu_table()

    def get_otu_table(self):
        """
        根据调用的level参数重构otu表
        :return:
        """
        if self.option('otutable').format == "meta.otu.tax_summary_dir":
            otu_path = self.option('otutable').get_table(self.option('level'))
        else:
            otu_path = self.option('otutable').prop['path']
        return otu_path

    def run(self):
        """
        运行
        """
        super(RandomforestTool, self).run()
        self.run_RandomForest_perl()

    def formattable(self, tablepath):
        with open(tablepath) as table:
            if table.read(1) == '#':
                newtable = os.path.join(self.work_dir, 'temp_format.table')
                with open(newtable, 'w') as w:
                    w.write(table.read())
                return newtable
        return tablepath

    def deal_index(self, table, i, n=''):
        '''
        i == 0, 改变 table 的行名, 生成文件 'reindex.table', 'reindex.map' 两个文件
        i == 1, 回复改变前的行名
        '''
        cmd = '/miniconda2/bin/perl ' +\
                self.config.PACKAGE_DIR + '/meta/scripts/deal_index.pl -i {} -r {}'.format(table, i)
        command = self.add_command('deal_index' + n, cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info(command.name + ' done!')
        else:
            self.set_error(command.name + ' failed')

    def replace_name(self, infile, outfile):
        """
        fix bug：20200630 qingchen.zhang
        根据前面建好的index名称按照原来的对应关系替换物种名称
        """
        index_map_file = os.path.join(self.work_dir, "reindex.map")
        self.index_map = {}
        self.old_index = {}
        with open(index_map_file, 'r') as f:
            for line in f:
                line = line.strip().split("\t")
                new_name = line[0]
                old_name = line[1]
                if new_name not in self.index_map:
                    self.index_map[new_name] = old_name
                if old_name not in self.old_index:
                    self.old_index[old_name] = new_name
        with open(infile, 'r') as inf, open(outfile, "w") as w:
            lines = inf.readlines()
            w.write(lines[0])
            for lin in lines[1:]:
                lin = lin.strip().split("\t")
                sp_name = lin[0]
                if sp_name in self.old_index:
                    lin[0] = self.old_index[sp_name]
                    w.write("\t".join(lin) + "\n")
        return outfile

    def run_RandomForest_perl(self):
        """
        运行RandomForest.pl
        """
        real_otu_path = self.formattable(self.otu_table)
        if self.option('norm_method'):
            new_otu = os.path.join(self.work_dir, 'normalized.xls')
            Normalization(real_otu_path, norm_meth=self.option('norm_method'), out=new_otu).run()
            real_otu_path = new_otu
        self.deal_index(real_otu_path, 0)  # 输入表行名用行序号代替 by xieshicahng 20200527
        real_otu_path = 'reindex.table'  # 替换完行序号的文件 by xieshichang 20200527
        cmd = self.config.SOFTWARE_DIR + '/miniconda2/bin/perl ' + self.cmd_path
        cmd += ' -i %s -o %s' % (real_otu_path, self.output_dir)
        if self.option('grouptable').is_set:
            cmd += ' -g %s -m %s' % ((self.option('grouptable').prop['group_scheme'])[0], self.group_table)
        cmd += ' -ntree %s' % (str(self.option('ntree')))
        cmd += ' -type %s' % (str(self.option('problem_type')))
        if self.option('predict_sample').is_set:
            replace_file = os.path.join(self.work_dir, "predict_sample.xls")
            predict_file = self.replace_name(self.option('predict_sample').prop['path'], replace_file)
            cmd += ' -pre_i %s' % (predict_file)
        cmd += ' -method %s' % (str(self.option('method')))
        self.logger.info(cmd)
        self.logger.info('运行RandomForest_perl.pl程序进行RandomForest计算')

        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('生成 cmd.r 文件成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成 cmd.r 文件失败')
            self.set_error('无法生成 cmd.r 文件', code="32703401")
        try:
            subprocess.check_output(self.config.SOFTWARE_DIR +
                                    '/program/R-3.3.1/bin/R --restore --max-ppsize=500000 --no-save < %s/randomForest.cmd.r' % (
                                        self.output_dir), shell=True)  # 修改R版本为R-3.3.1上sanger，tsg和tsanger的R版本依旧为R-3.3.3不做修改， add by zhujuan 20180511
            self.logger.info('RandomForest计算成功')
        except subprocess.CalledProcessError:
            self.logger.info('RandomForest计算失败')
            self.set_error('R运行计算RandomForest失败', code="32703402")
        self.logger.info('运行RandomForest_perl.pl程序进行RandomForest计算完成')
        os.rename(self.output_dir + "/randomForest.cmd.r", self.work_dir + "/randomForest.cmd.r")
        n = 1  # 在结果文件中恢复原有的行名称 by xieshichang 20200527
        for f in os.listdir(self.output_dir):  # by xieshichang 20200527
            self.deal_index(os.path.join(self.output_dir, f), 1, str(n))
            n += 1
        self.end()
