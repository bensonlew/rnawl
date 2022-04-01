# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict
class AssembleSoapDenovoModule(Module):
    """
    微生物基因组组装MP文库截取
    author: gaohao
    last_modify: 2018.03.25
    """
    def __init__(self, work_id):
        super(AssembleSoapDenovoModule, self).__init__(work_id)
        options = [
            {"name": "config", "type": "infile", "format": "bacgenome.config_file"},#组装的config文件
            {"name": "sample_name", "type": "string"},  # 样品名称
            {"name": "scafSeq", "type": "outfile", "format": "sequence.fasta_dir"},  # 输出文件,sample.scafSeq
            {"name": "kmer", "type": "string"},
            {"name": "D", "type": "int", "default": 1},
            {"name": "d", "type": "string", "default": "3,5,10"},
            {"name": "M", "type": "int", "default": 1},
            {"name": "R", "type": "bool", "default": True},
            {"name": "F", "type": "bool", "default": True},
            {"name": "u", "type": "string", "default": "unmask"},
            {"name": "G", "type": "int", "default": 50}
        ]
        self.sample_path ={}
        self.sample ={}
        self.modules = []
        self.sample_info = {}
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('config').is_set:
            raise OptionError('必须输入config文件夹', code="21400401")
        if not self.option('sample_name'):
            raise OptionError('必须输入sample_name样品名称！', code="21400402")

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def run_bac_assem(self):
        config_file = self.option("config").prop['path']
        n = 0
        kmers = self.option("kmer").split("-")
        kmers_list = list(range(int(kmers[0]),int(kmers[1]), 2))
        d_list = self.option("d").split(",")
        for k in kmers_list:
            for d in d_list:
                self.assemble = self.add_tool('assemble.bac_denova_ass')
                self.step.add_steps('assemble{}'.format(n))
                opts = {
                    "config": config_file,
                    "kmerFreqCutoff": int(d),
                    "sample_name": self.option('sample_name'),
                    "kmer":int(k),
                    "D": self.option("D"),
                    "M": self.option("M"),
                    "R": self.option("R"),
                    "F": self.option("F"),
                    "u": self.option("u"),
                    "G": self.option("G")
                }
                self.assemble.set_options(opts)
                step = getattr(self.step, 'assemble{}'.format(n))
                step.start()
                self.step.update()
                self.assemble.on('end', self.finish_update, 'assemble{}'.format(n))
                self.modules.append(self.assemble)
                n += 1
        self.logger.info(self.modules)
        self.on_rely(self.modules, self.set_output)
        self.step.update()
        for module in self.modules:
            module.run()

    def run(self):
        """
        运行
        :return:
        """
        super(AssembleSoapDenovoModule, self).run()
        self.run_bac_assem()

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
        for module in self.modules:
            self.linkdir(module.output_dir,self.output_dir)
        self.option('scafSeq').set_path(self.output_dir)
        self.end()

    def end(self):
        super(AssembleSoapDenovoModule, self).end()