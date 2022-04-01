# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import re
import time, shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir


class ScaffoldAbundModule(Module):
    """
    对多个bam文件分批次处理，最后汇总统计
    author: gaohao
    last_modify: 2019.08.07
    """

    def __init__(self, work_id):
        super(ScaffoldAbundModule, self).__init__(work_id)
        options = [
            {"name": "bam_dir", "type": "infile", "format": "metagbin.bam_dir"},  # bam的文件夹
            {"name": "maxbin_depth", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "metabat_depth", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.modules = []
        self.run_summary = self.add_tool("metagbin.coverage_sum")
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('bam_dir').is_set:
            raise OptionError('必须输入bam_dir文件夹')

    def scaffold_coverage(self):
        """
        以20个sam文件为单位；当文件数小于或等于20个时，直接跑；当文件数大于20个，没20个文件为一个单位运行，剩下的不足20个以一个单位运行
        :return:
        """
        files = os.listdir(self.option('bam_dir').prop['path'])
        filess = []
        for file in files:
            if file.endswith('.bam'):
                filess.append(file)
        n = 0
        j= 0
        list = []
        self.logger.info(len(filess))
        if len(filess) > 20:
            for k in filess:
                n += 1
                file = self.option("bam_dir").prop['path'] + "/" + k
                list.append(file)
                if n == 20:
                    j += 1
                    scf_abund = self.add_tool('metagbin.scaffold_abund')
                    des = " ".join(list)
                    opts = {
                        "bam": des,
                        "num": j,
                    }
                    scf_abund.set_options(opts)
                    self.modules.append(scf_abund)
                    list = []
                    n = 0
                elif 20 * j + n == len(filess) and n < 20:
                    scf_abund = self.add_tool('metagbin.scaffold_abund')
                    des = " ".join(list)
                    opts = {
                        "bam": des,
                        "num": j + 1,
                    }
                    scf_abund.set_options(opts)
                    self.modules.append(scf_abund)
            self.logger.info(self.modules)
            if len(self.modules) > 1:
                self.on_rely(self.modules, self.run_summary_coverage)
            for module in self.modules:
                module.run()
        else:
            for k in filess:
                file = self.option("bam_dir").prop['path'] + "/" + k
                list.append(file)
            des = " ".join(list)
            self.scf_abund = self.add_tool('metagbin.scaffold_abund')
            opts = {
                "bam": des,
                "num": 1,
            }
            self.scf_abund.set_options(opts)
            self.scf_abund.on("end", self.set_output)
            self.scf_abund.run()

    def run_summary_coverage(self):
        if os.path.exists(self.work_dir + "/depth"):
            shutil.rmtree(self.work_dir + "/depth")
        for module in self.modules:
            link_dir(module.output_dir, self.work_dir + "/depth")
        opts = {
            "depth_files": self.work_dir + "/depth/",
        }
        self.run_summary.set_options(opts)
        self.run_summary.on("end", self.set_output)
        self.run_summary.run()

    def run(self):
        """
        运行
        :return:
        """
        super(ScaffoldAbundModule, self).run()
        self.scaffold_coverage()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        files = os.listdir(self.option('bam_dir').prop['path'])
        filess = []
        for file in files:
            if file.endswith('.bam'):
                filess.append(file)
        if len(filess) > 20:
            os.system("cut -f 1,3 {} >{}".format(self.run_summary.output_dir + "/depth.txt", self.work_dir + '/maxbin.depth.txt'))
            for i in ['depth.txt','maxbin.depth.txt']:
                if os.path.exists(self.output_dir + "/" + i):
                    os.remove(self.output_dir + "/" + i)
            os.link(self.run_summary.output_dir + "/depth.txt", self.output_dir + "/depth.txt")
            os.link(self.work_dir + '/maxbin.depth.txt', self.output_dir + "/maxbin.depth.txt")
            self.option("metabat_depth", self.output_dir + "/depth.txt")
            self.option("maxbin_depth", self.output_dir + "/maxbin.depth.txt")
            self.end()
        else:
            os.system("cut -f 1,3 {} >{}".format(self.scf_abund.output_dir + "/depth1.txt", self.work_dir + '/maxbin.depth.txt'))
            for i in ['depth.txt','maxbin.depth.txt']:
                if os.path.exists(self.output_dir + "/" + i):
                    os.remove(self.output_dir + "/" + i)
            os.link(self.scf_abund.output_dir + "/depth1.txt", self.output_dir + "/depth.txt")
            os.link(self.work_dir + '/maxbin.depth.txt', self.output_dir + "/maxbin.depth.txt")
            self.option("metabat_depth", self.output_dir + "/depth.txt")
            self.option("maxbin_depth", self.output_dir + "/maxbin.depth.txt")
            self.end()

    def end(self):
        super(ScaffoldAbundModule, self).end()