# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict

class Blast2goModule(Module):
    """
    blast2go注释
    author: gaohao
    last_modify: 2018.09.27
    """
    def __init__(self, work_id):
        super(Blast2goModule, self).__init__(work_id)
        options = [
            {"name": "blastout", "type": "infile", "format": "align.blast.blast_xml"},
            {"name": "blastout2", "type": "infile", "format": "denovo_rna_v2.blast_xml"},
        ]
        self.add_option(options)
        self.modules=[]

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('blastout').is_set:
            raise OptionError("BLAST result file with xml type must be provide!", code="21203001")
        else:
            self.option("blastout2", self.option("blastout").prop['path'])


    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def run_blast2go(self):
        self.blast_nr_out = self.work_dir + '/temp_blast_nr'
        split_file=self.option("blastout2").change_blast_version2(self.blast_nr_out, sub_num=100000)
        a = 0
        for xml_file in split_file:
            self.b2g = self.add_tool('annotation.go.blast2go')
            a += 1
            self.step.add_steps('blast2go{}'.format(a))
            opts = {
                    "blastout": xml_file,
            }
            self.b2g.set_options(opts)
            step = getattr(self.step, 'blast2go{}'.format(a))
            step.start()
            self.step.update()
            self.b2g.on('end', self.finish_update, 'blast2go{}'.format(a))
            self.modules.append(self.b2g)
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
        super(Blast2goModule, self).run()
        self.run_blast2go()

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

    def set_output(self, event):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for module in self.modules:
            self.linkdir(module.output_dir,self.output_dir)
        self.end()

    def end(self):
        super(Blast2goModule, self).end()

