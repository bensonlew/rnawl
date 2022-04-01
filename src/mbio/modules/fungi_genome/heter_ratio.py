# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict

class HeterRatioModule(Module):

    def __init__(self, work_id):
        super(HeterRatioModule, self).__init__(work_id)
        options = [
            {"name": "scaf_fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "heter1", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "heter2", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "heter3", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "heter4", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "heter5", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "heter6", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.modules = []
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('scaf_fa').is_set:
            raise OptionError('必须输入组装文件', code="22101201")
        return True

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def single_heter(self):
        n = 0
        for heter in ['0.01','0.02','0.001','0.002','0.003','0.005']:
            self.single_heter_run = self.add_tool('fungi_genome.heter_ratio')
            self.step.add_steps('heter_ratio{}'.format('_' + heter))
            opts = {
                "scaf_fa": self.option('scaf_fa'),
                "heter": heter,
            }
            self.single_heter_run.set_options(opts)
            step = getattr(self.step, 'heter_ratio{}'.format('_' + heter))
            step.start()
            self.step.update()
            self.single_heter_run.on('end', self.finish_update,'heter_ratio{}'.format('_' + heter))
            self.modules.append(self.single_heter_run)
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
        super(HeterRatioModule, self).run()
        self.single_heter()

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
        self.option('heter1').set_path(self.output_dir + '/0.001_kmer.freq')
        self.option('heter2').set_path(self.output_dir + '/0.002_kmer.freq')
        self.option('heter3').set_path(self.output_dir + '/0.003_kmer.freq')
        self.option('heter4').set_path(self.output_dir + '/0.005_kmer.freq')
        self.option('heter5').set_path(self.output_dir + '/0.01_kmer.freq')
        self.option('heter6').set_path(self.output_dir + '/0.02_kmer.freq')
        self.end()
    def end(self):
        super(HeterRatioModule, self).end()



