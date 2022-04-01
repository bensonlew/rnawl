#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os,re
import shutil
import glob
from biocluster.core.exceptions import OptionError
from biocluster.module import Module



class AntismashModule(Module):
    """
    次级代谢产物的预测
    last_modify: 2018.04.14
    """

    def __init__(self, work_id):
        super(AntismashModule, self).__init__(work_id)
        options = [
            {"name": "gbk_dir", "type": "infile", "format": "gene_structure.gbk_dir"},  #
            {"name": "sample_name", "type": "string"},
            {"name": "analysis", "type": "string", "default": "complete"}  ###流程分析模式complete，uncomplete
        ]
        self.step.add_steps('antismash')
        self.add_option(options)
        self.modules = []



    def check_options(self):
        """
        检查参数
        """
        if not self.option('gbk_dir').is_set:
            raise OptionError("请设置基因组基因gbk文件夹不存在！", code="21400201")
        if not self.option('sample_name'):
            raise OptionError("请设置样品名称！", code="21400202")


    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def run_antismash_complete(self):
        """
        antismash 完成图运行
        :return:
        """
        if self.option('gbk_dir').is_set:
            files = os.listdir(self.option('gbk_dir').prop['path'])
            n = 0
            for file in os.listdir(self.option('gbk_dir').prop['path']):
                if re.search(r'.gbk', file):
                    self.antismash = self.add_tool('annotation.antismash_v5')
                    self.step.add_steps('antismash{}'.format(n))
                    opts = {
                        "genome_gbk": self.option('gbk_dir').prop['path'] + '/' + file,
                    }
                    self.antismash.set_options(opts)
                    step = getattr(self.step, 'antismash{}'.format(n))
                    step.start()
                    self.step.update()
                    self.antismash.on('end', self.finish_update, 'antismash{}'.format(n))
                    self.modules.append(self.antismash)
                    n += 1
            self.logger.info(self.modules)
            if len(self.modules) >1:
               self.on_rely(self.modules, self.set_output)
            elif len(self.modules) == 1:
                self.modules[0].on('end',self.set_output)
            self.step.update()
            for module in self.modules:
                module.run()

    def run(self):
        super(AntismashModule, self).run()
        self.run_antismash_complete()

    def set_output(self):
        self.logger.info("设置结果目录")
        for module in self.modules:
            files =os.listdir(module.output_dir)
            if len(files) > 1:
                self.logger.info(len(files))
                for file in files:
                    if re.search(r'antismash_anno.xls', file):
                        self.logger.info(module.output_dir + '/' + file)
                        self.logger.info(self.output_dir + '/' + self.option('sample_name') + file)
                        if os.path.exists(self.output_dir + '/' + self.option('sample_name') + '_' + file):
                            os.remove(self.output_dir + '/' + self.option('sample_name') + '_' + file)   #guanqing.zou 20180912
                        os.link(module.output_dir + '/' + file,
                            self.output_dir + '/' + self.option('sample_name') + '_' + file)
                    if re.search(r'gene_antismash.xls', file):
                        self.logger.info(module.output_dir + '/' + file)
                        self.logger.info(self.output_dir + '/' + self.option('sample_name') + file)
                        if os.path.exists(self.output_dir + '/' + self.option('sample_name') + '_' + file):
                            self.add_gene_txt(self.output_dir + '/' + self.option('sample_name') + '_' + file,module.output_dir + '/' + file)
                        else:
                            os.link(module.output_dir + '/' + file,
                                    self.output_dir + '/' + self.option('sample_name') + '_' + file)

        self.end()

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

    def add_gene_txt(self,raw_txt,add_txt):
        with open(raw_txt,"a") as t, open(add_txt) as f:
            add_lines = f.readlines()[1:]
            for i in add_lines:
                t.write(i)

    def end(self):
        super(AntismashModule, self).end()