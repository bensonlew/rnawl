# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# last_modify:2016.09.22

from biocluster.module import Module
import os
import types
from biocluster.core.exceptions import OptionError
from biocluster.config import Config

class CorrNetworkModule(Module):
    def __init__(self, work_id):
        super(CorrNetworkModule, self).__init__(work_id)
        self.step.add_steps('otu_association', 'corr_network_calc')
        options = [
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "lable", "type": "float", "default": 0.03},  # 设定lable（距离水平）的值
            {"name": "method", "type": "string", "default": "pearson"}, # 设定计算相似性的算法
            {"name": "coefficient", "type": "float", "default": 0.08},  # 相关性系数阈值
            {"name": "significance", "type": "float", "default": 0.05}
        ]
        self.add_option(options)
        self.otu_association = self.add_tool("meta.otu.otu_association")
        self.corr_network_calc = self.add_tool("meta.otu.corr_network_calc")
        self._end_info = 0

    def check_options(self):
        if not self.option("otutable").is_set:
            raise OptionError("必须要提供otutable文件")
        self.option('otutable').get_info()
        if self.option('otutable').prop['sample_num'] <= 2:
            raise OptionError('otu表的样本数目少于3，不可进行网络分析')
        return True

    def run(self):
        super(CorrNetworkModule, self).run()
        self.otu_association.on('end', self.corrnetworkcalc_run)
        self.corr_network_calc.on('end', self.end)
        self.otuassociation_run()

    def otuassociation_run(self):
        self.otu_association.set_options({
            "otutable": self.option('otutable').prop["path"],
            "lable": self.option("lable"),
            "method": self.option("method")
        })
        self.otu_association.run()

    def corrnetworkcalc_run(self):
        lable = self.option("lable")
        method = self.option("method")
        # corr_name = "shared."+ str(lable) + "." + str(method) +".otu.corr"
        corr_name = "shared." + str(lable) + "." + str(method) + ".corr"  # modified by hongdongxuan 20170324
        corr_table = os.path.join(self.work_dir, "OtuAssociation/output/", corr_name)
        print corr_table
        self.corr_network_calc.set_options({
            "corr_table": corr_table,
            "coefficient": self.option("coefficient"),
            "significance": self.option("significance")
        })
        self.corr_network_calc.run()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
                    # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_output(self):
        self.linkdir(self.otu_association.output_dir, 'otu_association')
        self.linkdir(self.corr_network_calc.output_dir, 'corr_network_calc')

    def end(self):
        self.set_output()
        super(CorrNetworkModule, self).end()
