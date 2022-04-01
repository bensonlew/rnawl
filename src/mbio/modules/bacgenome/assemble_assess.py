# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict
class AssembleAssessModule(Module):
    """
    微生物基因组组装结果评估
    author: gaohao
    last_modify: 2018.03.26
    """
    def __init__(self, work_id):
        super(AssembleAssessModule, self).__init__(work_id)
        options = [
            {"name": "seq_scaf", "type": "infile", "format": "sequence.fasta"},  # 最佳Gapcloser的scaffold文件
            {'name': 'sample_name', "type": "string"},  # 样本名
            {"name": "scaffold", "type": "outfile", "format": "sequence.fasta"},  # 输出文件
            {"name": "type", "type": "string"}
        ]
        self.step.add_steps('scaf_agp_contig', 'bac_stat', "busco")
        self.scaf_agp_contig = self.add_tool('assemble.scaf_agp_contig')
        self.bac_stat = self.add_tool('bacgenome.bac_assemble_stat')
        self.busco = self.add_tool('bacgenome.busco')
        self.tools = [self.bac_stat, self.busco]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('seq_scaf').is_set:
            raise OptionError('必须输入seq_scaf序列文件', code="21400301")
        if not self.option('sample_name'):
            raise OptionError('必须输入sample_name样品名称！', code="21400302")

    def run_gaploser_scaf(self):
        opts = {
            "seq_scaf": self.option('seq_scaf'),
            "sample_name": self.option('sample_name'),
        }
        self.scaf_agp_contig.set_options(opts)
        self.scaf_agp_contig.on('end', self.set_output, 'scaf_agp_contig')
        self.scaf_agp_contig.run()

    def run_bac_stat(self):
        opts = {
            "scaf_seq": self.scaf_agp_contig.option('scaffold'),
            "cont_seq": self.scaf_agp_contig.option('contig'),
            "sample_name": self.option('sample_name'),
        }
        self.bac_stat.set_options(opts)
        self.bac_stat.on('end', self.set_output, 'bac_stat')
        self.bac_stat.run()

    def run_busco(self):
        opts = {
            "scaf_fa": self.scaf_agp_contig.option('scaffold'),
            "database": "bacteria",
            "sample_name": self.option('sample_name'),
        }
        self.busco.set_options(opts)
        self.busco.on('end', self.set_output, 'busco')
        self.busco.run()

    def run(self):
        """
        运行
        :return:
        """
        super(AssembleAssessModule, self).run()
        if self.option("type") in ['busco']:
            self.on_rely(self.tools, self.end)
            self.scaf_agp_contig.on('end', self.run_bac_stat)
            self.scaf_agp_contig.on('end', self.run_busco)
            self.run_gaploser_scaf()
        else:
            self.bac_stat.on('end', self.end)
            self.scaf_agp_contig.on('end', self.run_bac_stat)
            self.run_gaploser_scaf()


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
                if os.path.exists(newfiles[i]):
                    os.remove(newfiles[i])
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.link(oldfiles[i], newdir)

    def set_output(self, event):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if event["data"] == "scaf_agp_contig":
            self.linkdir(self.scaf_agp_contig.output_dir , self.output_dir + '/assembly')
            self.option('scaffold',self.output_dir + '/assembly/' + self.option('sample_name') + '_scaf.fna')
        if event["data"] == "bac_stat":
            self.linkdir(self.bac_stat.output_dir + '/summary', self.output_dir + '/assembly')
            self.linkdir(self.bac_stat.output_dir + '/len', self.output_dir + '/len')
            self.linkdir(self.bac_stat.work_dir + '/scaffold',self.output_dir + '/scaffold')
            self.linkdir(self.bac_stat.work_dir + '/contig', self.output_dir + '/contig')
        if event["data"] == "busco":
            self.linkdir(self.busco.output_dir, self.output_dir + '/busco')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(AssembleAssessModule, self).end()