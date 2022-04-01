# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class GenePredictModule(Module):
    """
    宏基因组基因预测模块
    author: guhaidong
    last_modify: 2018.10.23
    """

    def __init__(self, work_id):
        super(GenePredictModule, self).__init__(work_id)
        option = [
            {"name": "input_fasta", "type": "infile", "format": "sequence.fasta_dir"},  # 输入文件夹，去掉小于最短contig长度的序列
            {"name": "min_gene", "type": "string", "default": "100"},  # 输入最短基因长度，如100
            {"name": "out", "type": "outfile", "format": "sequence.fasta"},  # 输出文件，基因预测输出路径
            {"name": "out_fa", "type": "outfile", "format": "sequence.fasta"}, # 输出文件，基因预测只含单拼结果路径
            {"name": "out_fa_mix", "type": "outfile", "format": "sequence.fasta"}, # 输出文件，基因预测混拼结果路径
            {"name": "out_faa", "type": "outfile", "format": "sequence.fasta"}, # 输出文件，基因预测只含单拼蛋白结果路径
            {"name": "out_faa_mix", "type": "outfile", "format": "sequence.fasta"}, # 输出文件，基因预测只含混拼蛋白结果路径
        ]
        self.add_option(option)
        self._is_mix = False  # 根据文件名称判断是否为混拼结果
        self.metagene_tools = []  # run_metagene并行运行tools列表
        self.metagene_tool = self.add_tool('gene_structure.metagene')
        self.metagene_stat = self.add_tool('gene_structure.metagene_stat')
        self.yasuo = self.add_tool('sequence.zip')
        self.len_distribute = self.add_tool('sequence.length_distribute')
        self.step.add_steps("metagene", "metagene_stat", "len_distribute")
        self.sum_tools = []  # 根据此变量移动结果文件
        self.zip_length = []

    def check_options(self):
        """
    检查参数
        :return:
        """
        if not self.option('input_fasta'):
            raise OptionError('必须输入拼接结果路径', code="22200301")
        return True

    '''
    def set_step_start(self, stepx):
        stepx.start()
        self.step.update

    def set_step_finish(self,stepx):
        stepx.finish()
        self.step.update()
    '''

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run_metagene(self):
        opts = {
            'min_gene': self.option('min_gene')
        }
        for f in self.option("input_fasta").fastas_full:
            opts['cut_more_scaftig'] = f
            opts['sample_name'] = os.path.basename(f).split('.contig')[0]
            if opts['sample_name'] in ["newbler", "Megahit_Mix"]:
                self._is_mix = True
            self.metagene_tool = self.add_tool('gene_structure.metagene')
            self.metagene_tool.set_options(opts)
            self.metagene_tools.append(self.metagene_tool)
        self.step.metagene.start()
        self.step.update()
        if len(self.metagene_tools) == 1:
            self.metagene_tools[0].on("end", self.run_metagene_stat)
        else:
            self.on_rely(self.metagene_tools, self.run_metagene_stat)
        for tool in self.metagene_tools:
            tool.run()
        self.step.metagene.finish()
        self.step.update()
        self.sum_tools.append(self.metagene_tools[0])

    def run_metagene_stat(self):
        self.sort_result()
        opts = {
            'contig_dir': self.work_dir + '/Predict'
        }
        self.metagene_stat.set_options(opts)
        self.metagene_stat.on('start', self.set_step, {'start': self.step.metagene_stat})
        self.metagene_stat.on('end', self.set_step, {'end': self.step.metagene_stat})
        self.metagene_stat.on('end', self.run_zip_length)
        self.metagene_stat.run()
        self.sum_tools.append(self.metagene_stat)

    def run_zip_length(self):
        """
        并行运行zip和length_distribute
        :return:
        """
        self.run_zip()
        self.run_len_distribute()
        for tool in self.zip_length:
            tool.run()
        self.on_rely(self.zip_length, self.set_output)
        self.sum_tools.append(self.len_distribute)
        #self.len_distribute.on('end', self.set_step, {'end': self.step.len_distribute})
        #self.step.update()

    def run_zip(self):
        """
        对预测结果进行压缩
        :return:
        """
        if self._is_mix:
            opts = {'file_path': self.metagene_stat.option('fasta').prop['path']}
        else:
            opts = {'file_dir': self.work_dir + '/Predict'}
        self.yasuo.set_options(opts)
        self.zip_length.append(self.yasuo)
        #self.yasuo.on('end', self.run_len_distribute)
        #self.yasuo.run()

    def run_len_distribute(self):
        opts = {
            'len_range': '200,400,500,600,800'
        }
        if self._is_mix:
            opts["fasta_dir"] = self.metagene_stat.output_dir
        else:
            opts["fasta_dir"] = self.work_dir + '/Predict'
        self.len_distribute.set_options(opts)
        self.len_distribute.on('start', self.set_step, {'start': self.step.len_distribute})
        self.zip_length.append(self.len_distribute)
        #self.len_distribute.on('end', self.set_output)
        #self.len_distribute.run()

    def sort_result(self):
        """
        将metagene的预测结果放在一起
        :return:
        """
        if os.path.exists(self.work_dir + '/Predict'):
            shutil.rmtree(self.work_dir + '/Predict')
        os.mkdir(self.work_dir + '/Predict')
        for tool in self.metagene_tools:
            files = os.listdir(tool.output_dir)
            for fasta in files:
                try:
                    os.link(os.path.join(tool.output_dir, fasta), os.path.join(self.work_dir, 'Predict', fasta))
                except OSError:
                    self.logger.info(os.path.join(self.work_dir, 'Predict', fasta) + '已存在')

    def run(self):
        """
        运行module
        :return:
        """
        super(GenePredictModule, self).run()
        self.run_metagene()


    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下所有文件到module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹路径
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.work_dir, dirname)
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
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.link(oldfiles[i], newdir)

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for i, tool in enumerate(self.sum_tools):
            if i == 0:
                pass
                self.linkdir(self.yasuo.output_dir, self.output_dir)
                # if self._is_mix:
                #     pass
                # else:
                #     self.linkdir(self.work_dir + '/Predict', self.output_dir)
            elif i == 1:
                self.linkdir(self.yasuo.output_dir, self.output_dir)
                if os.path.exists(self.output_dir + '/sample.metagene.stat'):
                    os.remove(self.output_dir + '/sample.metagene.stat')
                os.link(tool.output_dir + '/sample.metagene.stat', self.output_dir + '/sample.metagene.stat')
                # if self._is_mix:
                #     self.linkdir(tool.output_dir, self.output_dir)
                # else:
                #     if os.path.exists(self.output_dir + '/sample.metagene.stat'):
                #         os.remove(self.output_dir + '/sample.metagene.stat')
                #     os.link(tool.output_dir + '/sample.metagene.stat', self.output_dir + '/sample.metagene.stat')
            elif i == 2:
                if self._is_mix:
                    for f in os.listdir(tool.output_dir):
                        if f.startswith('Total'):
                            if os.path.exists(self.output_dir + '/len_distribute'):
                                if os.path.exists(self.output_dir + '/len_distribute/' + f):
                                    os.remove(self.output_dir + '/len_distribute/' + f)  # 更正删除操作 顾海东 @20180104
                            else:
                                os.mkdir(self.output_dir + '/len_distribute')
                            os.link(tool.output_dir + '/' + f, self.output_dir + '/len_distribute/' + f)
                else:
                    self.linkdir(tool.output_dir, self.output_dir + '/len_distribute')
        out_path = self.metagene_stat.option('faa').prop['path']
        # self.option('out').set_path(self.metagene_stat.option('fasta').prop['path'])
        self.logger.info(">>>>>>>>>>>>>>START PRINT OUTFILE PATH<<<<<<<<<<<<<<")
        self.option('out').set_path(out_path)
        self.logger.info(self.option('out').prop['path'])
        if self._is_mix:
            self.option('out_fa').set_path(self.metagene_stat.option('fasta_sample').prop['path'])
            self.option('out_faa').set_path(self.metagene_stat.option('faa_sample').prop['path'])
            self.option('out_fa_mix').set_path(self.metagene_stat.option('fasta_mix').prop['path'])
            self.option('out_faa_mix').set_path(self.metagene_stat.option('faa_mix').prop['path'])
            self.logger.info(self.option('out_fa').prop['path'])
            self.logger.info(self.option('out_fa_mix').prop['path'])
        else:
            self.option('out_fa').set_path(self.metagene_stat.output_dir + '/Total.metagene.fa')
            self.logger.info(self.option('out_fa').prop['path'])
        self.logger.info(">>>>>>>>>>>>>PRINT OUTFILE PATH END!!!<<<<<<<<<<<<<")
        self.logger.info("设置基因预测结果目录成功")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [",", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(GenePredictModule, self).end()
