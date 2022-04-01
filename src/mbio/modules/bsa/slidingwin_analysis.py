#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = "wangzhaoyue"

from __future__ import division
import os
from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class SlidingwinAnalysisModule(Module):
    """
    version 1.0
    author: wangzhaoyue
    last_modify: 2018.02.23
    """
    def __init__(self, work_id):
        super(SlidingwinAnalysisModule, self).__init__(work_id)
        options = [
            {"name": "pop_vcf", "type": "string"},  # pop.final.vcf文件
            {"name": "wp", "type": "string"},  # 野生型亲本名称
            {"name": "mp", "type": "string"},  # 突变型亲本名称
            {"name": "wb", "type": "string"},  # 野生型混池名称
            {"name": "mb", "type": "string"},  # 突变型混池名称
            {"name": "pdep", "type": "int", "default": 10},  # 亲本测序深度，默认10X
            {"name": "bdep", "type": "int", "default": 10},  # 混池测序深度，默认10X
            {"name": "popt", "type": "string", "default": "F2"},  # 群体类型，默认F2（仅qtlseq.pl使用）
            {"name": "variant_type", "type": "string", "default": "ALL"},  # 只输出变异类型是该类型的结果
            # {"name": "col", "type": "string", "default": "1,2,10"},  # the col of chr pos index
            {"name": "method", "type": "string", "default": "bp"},  # 模式，bp/num
            {"name": "win", "type": "string", "default": "2000000"},  # the window size
            {"name": "step", "type": "string", "default": "10000"},  # the step size
            {"name": "abs", "type": "string", "default": "0"}
        ]
        self.add_option(options)
        self.index_calc = self.add_tool('bsa.index_calc')
        self.slidingwin = self.add_tool('bsa.slidingwin')
        self.step.add_steps('index_calc', 'slidingwin')
        self.tools = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('pop_vcf'):
            raise OptionError('必须输入pop.final.vcf', code="21500201")
        if not self.option('mb'):
            raise OptionError('必须输入突变型混池名称', code="21500202")
        if self.option('variant_type') not in ["ALL", "SNP", "INDEL"]:
            raise OptionError('输出的突变类型只能是"ALL,SNP,INDEL"中的一种', code="21500203")
        if self.option('method') not in ["bp", "num"]:
            raise OptionError('模式只有bp和num两种', code="21500204")
        return True

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def index_calc_run(self):
        self.step.index_calc.start()
        self.step.update()
        opts = {
            "pop_vcf": self.option('pop_vcf'),
            "mb": self.option('mb'),
            "pdep": self.option('pdep'),
            "bdep": self.option('bdep'),
            "popt": self.option('popt'),
            "variant_type": self.option('variant_type'),
        }
        if self.option('wp'):
            opts['wp'] = self.option('wp')
        if self.option('mp'):
            opts['mp'] = self.option('mp')
        if self.option('wb'):
            opts['wb'] = self.option('wb')
        self.index_calc.set_options(opts)
        self.index_calc.on('end', self.finish_update, 'index_calc')
        self.index_calc.run()
        self.tools.append(self.index_calc)

    def slidingwin_run(self):
        if self.option('wb'):
            col = "1,2,3,4"
        else:
            col = "1,2,3"  # 程序以逗号分割，只有长度为3则认为是一个混池
        self.step.slidingwin.start()
        self.step.update()
        opts = {
            "index_file": self.index_calc.output_dir + '/index-calc.result.index',
            "col": col,
            "method": self.option('method'),
            "win": self.option('win'),
            "step": self.option('step'),
            "abs": self.option("abs")
        }
        self.slidingwin.set_options(opts)
        self.slidingwin.on('end', self.finish_update, 'slidingwin')
        self.slidingwin.run()
        self.tools.append(self.slidingwin)

    def run(self):
        """
        运行
        :return:
        """
        self.index_calc.on('end', self.slidingwin_run)
        self.slidingwin.on('end', self.set_output)
        self.index_calc_run()
        super(SlidingwinAnalysisModule, self).run()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        newdir = os.path.join(self.work_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        allfiles = os.listdir(dirpath)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.link(oldfiles[i], newdir)

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for tool in self.tools:
            self.linkdir(tool.output_dir, self.output_dir)
        self.logger.info("设置结果目录成功")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(SlidingwinAnalysisModule, self).end()
