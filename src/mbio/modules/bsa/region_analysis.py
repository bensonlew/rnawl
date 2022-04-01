#!/usr/bin/env python
# -*- coding: utf-8 -*-

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import os


class RegionAnalysisModule(Module):
    """
    区域分析的模块，调用的tool有 region_variant， region_gene， region_vcf，三者之间的关系是并发运行，然后都运行完成之后再进行end
    version 1.0
    author: HONGDONG
    last_modify: 20180224
    laste modified by zengjing 增加参数wp、mp、wb、mb、step
    """
    def __init__(self, work_id):
        super(RegionAnalysisModule, self).__init__(work_id)
        options = [
            {"name": "i_c_result", "type": "string"},  # index_calc_result
            {"name": "pop_summary", "type": "string"},  # pop.summary
            {"name": "p_f_vcf", "type": "string"},  # pop.final.vcf
            {"name": "s_w_select", "type": "string"},  # sliding-win.threshold.select
            {"name": "wp", "type": "string", "default": ""},  # 野生型亲本名称
            {"name": "mp", "type": "string", "default": ""},  # 突变型亲本名称
            {"name": "wb", "type": "string", "default": ""},  # 野生型混池名称
            {"name": "mb", "type": "string", "default": ""}  # 突变型混池名称
            # {"name": "step", "type": "int"}  # 滑窗策略较小的数值
        ]
        self.add_option(options)
        self.region_gene = self.add_tool('bsa.region_gene')
        self.region_variant = self.add_tool('bsa.region_variant')
        self.region_vcf = self.add_tool('bsa.region_vcf')

    def check_options(self):
        """
        检查参数
        """
        if not self.option('i_c_result'):
            raise OptionError('必须提供index_calc_result结果表', code="21500101")
        if not self.option('pop_summary'):
            raise OptionError('必须提供pop_summary结果表', code="21500102")
        if not self.option('p_f_vcf'):
            raise OptionError('必须提供index_calc_result结果表', code="21500103")
        if not self.option('s_w_select'):
            raise OptionError('必须提供sliding-win.threshold.select结果表', code="21500104")
        if not self.option("mb"):
            raise OptionError('必须提供突变型混池mb', code="21500105")
        # if not self.option("step"):
        #     raise OptionError('必须提供滑窗策略较小值')

    def region_gene_run(self):
        self.region_gene.set_options({
            'pop_summary': self.option('pop_summary'),
            's_w_select': self.option('s_w_select')
        })
        self.region_gene.on('end', self.set_output, 'region_gene')
        self.region_gene.run()

    def region_variant_run(self):
        self.region_variant.set_options({
            'i_c_result': self.option('i_c_result'),
            's_w_select': self.option('s_w_select')
        })
        self.region_variant.on('end', self.set_output, 'region_variant')
        self.region_variant.run()

    def region_vcf_run(self):
        self.region_vcf.set_options({
            'p_f_vcf': self.option('p_f_vcf'),
            's_w_select': self.option('s_w_select'),
            'wp': self.option('wp'),
            'mp': self.option('mp'),
            'wb': self.option('wb'),
            'mb': self.option('mb')
            # 'step': self.option('step')
        })
        self.region_vcf.on('end', self.set_output, 'region_vcf')
        self.region_vcf.run()

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
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'region_gene':
            self.linkdir(obj.output_dir, 'region_gene')
        elif event['data'] == 'region_variant':
            self.linkdir(obj.output_dir, 'region_variant')
        elif event['data'] == 'region_vcf':
            self.linkdir(obj.output_dir, 'region_vcf')
        else:
            pass

    def run(self):
        super(RegionAnalysisModule, self).run()
        self.region_variant_run()
        self.region_vcf_run()
        self.region_gene_run()
        self.on_rely([self.region_variant, self.region_vcf, self.region_gene], self.end)
