# -*- coding: utf-8 -*-
# __author__ = 'liuwentian'
# modified 2018.1204


from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from biocluster.config import Config
import os, re
import json
import time, datetime
from bson.objectid import ObjectId


class GeneSampleSeqModule(Module):
    """
    基因详情页module
    """
    def __init__(self, work_id):
        super(GeneSampleSeqModule, self).__init__(work_id)
        options = [
            {"name": "ref_fa", "type": "string"},  # ref.fa
            {"name": "ref_gff", "type": "string"},  # gff文件
            {"name": "specimen_ids", "type": "string"},  # 样本
            {"name": "pop_final_vcf", "type": "string"},  # pop_final_vcf
            {"name": "gene_id", "type": "string"},  # gene_id
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "project_type", "type": "string"},
            {"name": "ssr_path", "type": "string"},  # 基因组路径
            {"name": "gene_name", "type": "string"},
            {"name": "is_wgs_result", "type": "string"}
        ]
        self.add_option(options)
        self.gene_sample_seq = self.add_tool("dna_evolution.gene_sample_seq")

    def check_options(self):
        if not self.option("ref_fa"):
            raise OptionError("请设置ref_fa文件")
        if not self.option("ref_gff"):
            raise OptionError("请设置ref_gff文件")
        if not self.option("specimen_ids"):
            raise OptionError("请设置specimen_ids")
        if not self.option("pop_final_vcf"):
            raise OptionError("请设置pop_final_vcf文件")
        if not self.option("gene_id"):
            raise OptionError("请设置gene_id")
        if not self.option("main_id"):
            raise OptionError("请设置main_id")
        if not self.option("update_info"):
            raise OptionError("请设置update_info")
        if not self.option("project_type"):
            raise OptionError("请设置project_type")
        if not self.option("ssr_path"):
            raise OptionError("请设置ssr_path")
        if not self.option("gene_name"):
            raise OptionError("请设置gene_name")
        if not self.option("is_wgs_result"):
            raise OptionError("请设置is_wgs_result")

    def run_gene_sample_seq(self):
        """
        """
        self.gene_sample_seq.set_options({
            "ref_fa": self.option("ref_fa"),
            "ref_gff": self.option("ref_gff"),
            "specimen_ids": self.option("specimen_ids"),
            "pop_final_vcf": self.option("pop_final_vcf"),
            "gene_id": self.option("gene_id"),
            "main_id": self.option("main_id"),
            "project_type": self.option("project_type"),
            "ssr_path": self.option("ssr_path"),
            "gene_name": self.option("gene_name"),
            "is_wgs_result": self.option("is_wgs_result"),
            "fa_path": self.parent._sheet.output
        })
        self.gene_sample_seq.on('end', self.set_output, "gene_sample_seq")
        self.gene_sample_seq.on('end', self.end)
        self.gene_sample_seq.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'gene_sample_seq':
            self.linkdir(obj.output_dir, 'gene_sample_seq')
        else:
            pass

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

    def run(self):
        super(GeneSampleSeqModule, self).run()
        self.run_gene_sample_seq()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        super(GeneSampleSeqModule, self).end()
