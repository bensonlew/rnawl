# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import re
import time,shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir

class CheckmModule(Module):
    """
    多个bins进行分析
    author: gaohao
    last_modify: 2019.03.07
    """
    def __init__(self, work_id):
        super(CheckmModule, self).__init__(work_id)
        options = [
            {"name": "bin_dir", "type": "infile", "format": "sequence.fasta_dir"},#生成bins目录
        ]
        self.modules = []
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('bin_dir').is_set:
            raise OptionError('必须输入bin_dir文件夹')

    def run_checkm(self):
        files = os.listdir(self.option("bin_dir").prop['path'])
        n=0
        for k in files:
            self.checkm = self.add_tool('metagbin.checkm_bin')
            self.step.add_steps('checkm{}'.format(n))
            file = self.option("bin_dir").prop['path'] + "/" + k
            if os.path.exists(self.checkm.work_dir + '/bin'):
                shutil.rmtree(self.checkm.work_dir + '/bin')
            os.mkdir(self.checkm.work_dir + '/bin')
            os.link(file,self.checkm.work_dir + '/bin/' + k)
            opts = {
                "bin_dir": self.checkm.work_dir + '/bin',
                }
            self.checkm.set_options(opts)
            self.modules.append(self.checkm)
            self.checkm.on('end', 'checkm{}'.format(n))
            n +=1
        self.logger.info(self.modules)
        if len(self.modules) > 1:
            self.on_rely(self.modules, self.run_stat)
        else:
            self.modules[0].on('end', self.run_stat)
        for module in self.modules:
            module.run()

    def run_stat(self):
        """
        将所有bins的结果进行合并，
        :return:
        """
        for module in self.modules:
            link_dir(module.output_dir,self.work_dir + "/data")
        with open (self.work_dir + "/all.bin.summary",'w') as f,open(self.work_dir + "/all.marker.xls","w") as g:
            header = ["Bin Id","Genome size","Scaffolds Num","Longest scaffold","N50 (scaffolds)","Mean scaffold length","Completeness","Contamination","Strain heterogeneity","Domain"]
            f.write("\t".join(header) + "\n")
            header2=["Bin Id","markers","marker sets","0","1","2","3","4","5+"]
            g.write("\t".join(header2) + "\n")
            files = os.listdir(self.work_dir + "/data")
            for file in files:
                if re.search(r'.bin.summary.xls',file):
                    with open(self.work_dir + "/data/" + file, 'r') as j:
                        lines = j.readlines()
                        f.write(lines[1])
                elif re.search(r'.marker.xls',file):
                    with open(self.work_dir + "/data/" + file, 'r') as d:
                        lines = d.readlines()
                        g.write(lines[1])
        self.sort_file(self.work_dir + "/all.bin.summary",self.work_dir + "/all.bin.summary.xls")
        self.set_output()

    def sort_file(self,file,output):
        with open (file,"r") as f,open (output,"w") as g:
            dict = {}
            lines =f.readlines()
            g.write(lines[0])
            for lin in lines[1:]:
                line =lin.rstrip('\r\n\t').split("\t")
                dict[lin] =float(line[6])
            list = sorted(dict,key=dict.get,reverse=True)
            for i in list:
                g.write(i)

    def run(self):
        """
        运行
        :return:
        """
        super(CheckmModule, self).run()
        self.run_checkm()
        

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for i in ['all.marker.xls','all.bin.summary.xls']:
            if os.path.exists(self.output_dir + "/" + i):
                os.remove(self.output_dir + "/" + i)
            os.link(self.work_dir + "/" + i,self.output_dir + "/" + i)
        self.end()

    def end(self):
        super(CheckmModule, self).end()