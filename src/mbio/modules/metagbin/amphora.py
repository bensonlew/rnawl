# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir

class AmphoraModule(Module):
    """
    多个bins进行注释
    author: gaohao
    last_modify: 2019.01.14
    """
    def __init__(self, work_id):
        super(AmphoraModule, self).__init__(work_id)
        options = [
            {"name": "bin_dir", "type": "infile", "format": "sequence.fasta_dir"},#生成bins目录
            {"name": "taxon_file", "type": "infile","format": "metagbin.file_table"},  # bin
        ]
        self.modules = []
        self.amphora_anno = self.add_tool('metagbin.amphora_anno')
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('bin_dir').is_set:
            raise OptionError('必须输入bin_dir文件夹')
        if not self.option('taxon_file'):
            raise OptionError('必须输入taxon_file样品名称！')

    def run_amphora(self):
        files = os.listdir(self.option("bin_dir").prop['path'])
        dict = self.get_taxon(self.option("taxon_file").prop['path'])
        n=0
        for k in files:
            m=re.search('(.*).fa',k)
            print m.group(1)
            if m.group(1) in dict:
                self.amphora = self.add_tool('metagbin.amphora')
                opts = {
                    "bin_fa": self.option("bin_dir").prop['path'] + '/' + k,
                    "kingdom": dict[m.group(1)],
                    "bin_name": m.group(1),
                }
                self.amphora.set_options(opts)
                self.modules.append(self.amphora)
                self.amphora.on('end', 'amphora{}'.format(n))
                n +=1
        self.logger.info(self.modules)
        if len(self.modules) > 1:
            self.on_rely(self.modules, self.run_amphora_summ)
        elif len(self.modules) == 1:
            self.modules[0].on("end", self.run_amphora_summ)
        for module in self.modules:
            module.run()

    def run_amphora_summ(self):
        for module in self.modules:
            link_dir(module.output_dir,self.work_dir + "/anno")
        opts = {
            "amphora_dir": self.work_dir + "/anno",
        }
        self.amphora_anno.set_options(opts)
        self.amphora_anno.on("end",self.set_output)
        self.amphora_anno.run()

    def run(self):
        """
        运行
        :return:
        """
        super(AmphoraModule, self).run()
        self.run_amphora()


    def get_taxon(self, file):
        bin_dict={}
        with open (file,'r') as f:
            files=f.readlines()
            for file in files[1:]:
                bins=file.rstrip('\r\n').split('\t')
                bin_dict[bins[0]]=bins[9]
        return bin_dict

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if os.path.exists(self.output_dir + "/summary.anno.xls"):
            os.remove(self.output_dir + "/summary.anno.xls")
        os.link(self.amphora_anno.output_dir + '/summary.anno.xls',self.output_dir + "/summary.anno.xls")
        self.end()

    def end(self):
        super(AmphoraModule, self).end()