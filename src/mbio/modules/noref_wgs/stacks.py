# -*- coding: utf-8 -*-
# __author__ = 'wentian.liu'
# modified 2019.1.10


from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from biocluster.config import Config
from bson.objectid import ObjectId
import os, re
import json
import time, datetime
import gevent
import shutil


class StacksModule(Module):
    """
    stacks流程的module
    """
    def __init__(self, work_id):
        super(StacksModule, self).__init__(work_id)
        options = [
            {"name": "fastq_list", "type": 'infile', 'format': "noref_wgs.fastq_list"},  # 03.fq.list
            {"name": "group_list", "type": "string"},  # group.list
            {"name": 'is_split', 'type': 'bool', 'default': True}, # 设定是否需要拆分进行分析，默认是不用拆分的
            {"name": 'split_num', 'type': 'int', "default": 15}, # 将样本分成多少份进行cstack
        ]
        self.add_option(options)
        self.ustacks = self.add_module("noref_wgs.ustacks")
        if self.option('is_split'):
            self.cstackssplit = self.add_module("noref_wgs.cstacks")
        else:
            self.cstacks = self.add_tool("noref_wgs.cstacks")        
        self.sstacks = self.add_module("noref_wgs.sstacks")
        self.genotype = self.add_module("noref_wgs.genotype")
        self.stack_stat = self.add_tool("noref_wgs.stack_stat")
        self.consensus_stat = self.add_tool("noref_wgs.consensus_stat")

    def check_options(self):
        if not self.option("fastq_list").is_set:
            raise OptionError("please input sample's fastq list file!", code="25500403")

    def run_ustacks(self):
        """
        """
        self.ustacks.set_options({
            "fastq_list": self.option("fastq_list").prop['path']
        })
        self.ustacks.on("end", self.set_output, "ustacks")
        if self.option('is_split'):
            self.ustacks.on('end', self.run_cstackssplit)
        else:
            self.ustacks.on('end', self.run_cstacks)
        self.ustacks.run()

    def run_cstacks(self):
        """
        运行cstacks没有进行拆分分析
        """
        if self.option("group_list"):
            self.cstacks.set_options({
                "group_list": self.option("group_list"),
                "ustacks_output": self.ustacks.output_dir
            })
        else:
            self.cstacks.set_options({
                "ustacks_output": self.ustacks.output_dir
            })
        self.cstacks.on("end", self.set_output, "cstacks")
        self.cstacks.on('end', self.run_sstacks)
        self.cstacks.run()

    def run_cstackssplit(self):
        """
        运行cstacks按照某种情形进行分割，然后进行分析，能够解决大内存的问题。
        """
        if self.option("group_list"):
            self.cstackssplit.set_options({
                "group_list": self.option("group_list"),
                "ustacks_path": self.ustacks.output_dir,
                "split_num": self.option("split_num")
            })
        else:
            self.cstackssplit.set_options({
                "ustacks_path": self.ustacks.output_dir,
                "split_num": self.option("split_num")
            })
        self.cstackssplit.on("end", self.set_output, "cstacks")
        self.cstackssplit.on('end', self.run_sstacks)
        self.cstackssplit.run()

    def run_sstacks(self):
        """
        sample_loci_dir:03output
        catalog_dir:04output
        """
        if self.option('is_split'):
            catalog_dir = os.path.join(self.cstackssplit.output_dir, "cstacksmerge")
        else:
            catalog_dir = self.cstacks.output_dir
        self.sstacks.set_options({
            "sample_loci_dir": self.ustacks.output_dir,
            "catalog_dir": catalog_dir
        })
        self.sstacks.on("end", self.set_output, "sstacks")
        self.sstacks.on('end', self.run_genotype)
        self.sstacks.run()

    def run_genotype(self):
        """
        ustacks_output:03output
        cstacks_output:04output
        sstacks_output:05output
        gro_list:04 group.list
        """
        if self.option('is_split'):
            gro_list = os.path.join(self.cstackssplit.work_dir, "group.list")
            cstacks_output = os.path.join(self.cstackssplit.output_dir, "cstacksmerge")
        else:
            cstacks_output = self.cstacks.output_dir
            gro_list = os.path.join(self.cstacks.work_dir, "group.list")
        self.genotype.set_options({
            "ustacks_output": self.ustacks.output_dir,
            "cstacks_output": cstacks_output,
            "sstacks_output": self.sstacks.output_dir,
            "gro_list": gro_list
        })
        self.genotype.on("end", self.set_output, "genotype")
        self.genotype.on('end', self.run_stack_stat)
        self.genotype.run()

    def run_stack_stat(self):
        """
        snp_vcf:populations.snps.vcf
        """
        list_file = os.path.join(self.output_dir, "ustacks/ustacks.list")
        snp_vcf = os.path.join(self.genotype.work_dir, "genotype_input/populations.snps.vcf")
        self.stack_stat.set_options({
            "snp_vcf": snp_vcf,
            "list_file": list_file
        })
        self.stack_stat.on("end", self.set_output, "stack_stat")
        self.stack_stat.on('end', self.run_consensus_stat)
        self.stack_stat.run()

    def run_consensus_stat(self):
        """
        """
        if self.option('is_split'):
            catalog_path = self.output_dir + "/cstacks/cstacksmerge/catalog.tags.tsv.gz"
            sample_list_path = os.path.join(self.cstackssplit.work_dir, "sample.list")
        else:
            sample_list_path = self.cstacks.work_dir + "/sample.list"
            catalog_path = self.output_dir + "/cstacks/catalog.tags.tsv.gz"
        populations_snps_vcf_path = self.output_dir + "/genotype/genotype/populations.snps.vcf"
        ustacks_path = self.output_dir + "/ustacks"
        with open(sample_list_path, "r")as fr:
            lines = fr.readlines()
            sample_num = len(lines)
        self.consensus_stat.set_options({
            "catalog_path": catalog_path,
            "sample_num": sample_num,
            "populations_snps_vcf": populations_snps_vcf_path,
            "ustacks_path": ustacks_path
        })
        self.consensus_stat.on("end", self.set_output, "consensus_stat")
        self.consensus_stat.on('end', self.end)
        self.consensus_stat.run()

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
        if event['data'] == 'stack_stat':
            self.linkdir(obj.output_dir, 'stack_stat')
        elif event['data'] == 'ustacks':
            self.linkdir(obj.output_dir, 'ustacks')
        elif event['data'] == 'cstacks':
            self.linkdir(obj.output_dir, 'cstacks')
        elif event['data'] == 'sstacks':
            self.linkdir(obj.output_dir, 'sstacks')
        elif event['data'] == 'genotype':
            self.linkdir(obj.output_dir, 'genotype')
        elif event['data'] == 'consensus_stat':
            self.linkdir(obj.output_dir, 'consensus_stat')
        else:
            pass

    def run(self):
        super(StacksModule, self).run()
        self.run_ustacks()

    def end(self):
        super(StacksModule, self).end()
