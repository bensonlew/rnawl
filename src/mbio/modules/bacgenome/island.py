#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import shutil
import glob
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.files.sequence.file_sample import FileSampleFile


class IslandModule(Module):
    """
    微生物基因组质控，二代数据：PE文库或MP文库；三代数据：pacbio数据
    last_modify: 2018.03.19
    """

    def __init__(self, work_id):
        super(IslandModule, self).__init__(work_id)
        options = [
            {"name": "fa_dir", "type": "infile", "format": "sequence.fasta_dir"},  #
            {"name": "gbk_dir", "type": "infile", "format": "gene_structure.gbk_dir"},  #
            {"name": "anno", "type": "infile", "format": "sequence.profile_table"},
            {"name": "sample_name", "type": "string"},
            {"name": "analysis", "type": "string", "default": "uncomplete"}  ###流程分析模式complete，uncomplete
        ]
        self.step.add_steps('island_dimob', 'island_islander','isaland_stat')
        self.add_option(options)
        self.diomb = self.add_tool('bacgenome.island_dimob')
        self.islander = self.add_tool('bacgenome.island_islander')
        self.island = self.add_tool('bacgenome.island_stat')
        self.list=[self.diomb,self.islander]


    def check_options(self):
        """
        检查参数
        """
        if not self.option('gbk_dir').is_set:
            raise OptionError("请设置基因组基因gbk文件夹不存在！", code="21401601")
        if not self.option('fa_dir').is_set:
            raise OptionError("请设置基因组学列文件夹！", code="21401602")
        if not self.option('anno').is_set:
            raise OptionError("请设置基因注释总览表！", code="21401603")


    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def run_island_dimob(self):
        """
        island_dimob运行
        :return:
        """
        self.logger.info("正在island_dimob开始")
        self.diomb.set_options({
            'gbk_dir': self.option('gbk_dir')
        })
        self.diomb.on('end', self.set_output, 'island_dimob')
        self.diomb.run()
        self.step.island_dimob.finish()
        self.step.update()
        self.logger.info("运行island_dimob结束")

    def run_island_islander(self):
        """
        island_islander运行
        :return:
        """
        self.logger.info("正在island_islander开始")
        self.islander.set_options({
            'fa_dir':self.option("fa_dir")
        })
        self.islander.on('end', self.set_output, 'island_islander')
        self.islander.run()
        self.step.island_islander.finish()
        self.step.update()
        self.logger.info("island_islander处理结束")

    def run_island_stat(self):
        """
        对合并的island进行统计
        :return:
        """
        self.logger.info("正在对island数据统计开始")
        self.island.set_options({
            'diomb':self.diomb.option('out'),
            'islander': self.islander.option('out'),
            'anno': self.option('anno'),
            'sample_name': self.option("sample_name"),
            'analysis':self.option("analysis"),
        })
        self.island.on("end",self.set_output,'isaland_stat')
        self.island.run()
        self.step.isaland_stat.finish()
        self.step.update()

    def run(self):
        super(IslandModule, self).run()
        self.run_island_dimob()
        self.run_island_islander()
        self.on_rely(self.list,self.run_island_stat)
        self.island.on('end',self.end)

    def set_output(self, event):
        self.logger.info("设置结果目录")
        if event['data'] == 'isaland_stat':
            if os.path.exists(self.output_dir + "/Genomic_Islands"):
                shutil.rmtree(self.output_dir + "/Genomic_Islands")
            shutil.copytree(self.island.output_dir ,self.output_dir + "/Genomic_Islands")

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

    def end(self):
        super(IslandModule, self).end()