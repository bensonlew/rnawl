# -*- coding: utf-8 -*-
# __author__ = 'wentian.liu'
# modified 2018.0103


from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from biocluster.config import Config
from bson.objectid import ObjectId
import os, re
import json
import time, datetime
import gevent
import shutil


class GenotypeModule(Module):
    """
    群体进化，单倍体图谱接口调用程序
    """
    def __init__(self, work_id):
        super(GenotypeModule, self).__init__(work_id)
        options = [
            {"name": "ustacks_output", "type": "string"},  # 03.ustacks结果目录
            {"name": "cstacks_output", "type": "string"},  # 04.cstacks结果目录
            {"name": "sstacks_output", "type": "string"},  # 05.sstacks结果目录
            {"name": "gro_list", "type": "string"}
        ]
        self.add_option(options)
        self.genotype = self.add_tool("noref_wgs.genotype")
        self.tag_generate = self.add_tool("noref_wgs.tag_generate")
        self.vcf_convert = self.add_tool("noref_wgs.vcf_convert")

    def check_options(self):
        if not self.option("ustacks_output"):
            raise OptionError("请设置ustacks_output", code="25500109")
        if not self.option("cstacks_output"):
            raise OptionError("请设置cstacks_output", code="25500110")
        if not self.option("sstacks_output"):
            raise OptionError("请设置sstacks_output", code="25500111")
        if not self.option("gro_list"):
            raise OptionError("请设置sstacks_output", code="25500112")

    def run_genotype(self):       # 1
        """
        """
        self.genotype.set_options({
            "input_file": os.path.join(self.work_dir, "genotype_input"),
            "gro_list": self.option("gro_list")
        })
        self.genotype.on("end", self.set_output, "genotype")
        self.genotype.on('end', self.run_tag_generate)
        self.genotype.run()

    def run_tag_generate(self):
        """
        :return:
        """
        vcf_path = os.path.join(os.path.join(self.work_dir, "genotype_input"), "populations.snps.vcf")
        tsv_path = os.path.join(self.option("cstacks_output"), "catalog.tags.tsv.gz")
        if not os.path.exists(vcf_path):
            tool_vcf_path = os.path.join(self.genotype.output_dir, "populations.snps.vcf")
            os.link(tool_vcf_path, vcf_path)
        if not os.path.exists(tsv_path):
            tool_tsv_path = os.path.join(self.genotype.output_dir, "catalog.tags.tsv.gz")
            os.link(tool_tsv_path, tsv_path)
        self.tag_generate.set_options({
            "vcf_path": os.path.join(os.path.join(self.work_dir, "genotype_input"), "populations.snps.vcf"),
            "tsv_path": os.path.join(os.path.join(self.work_dir, "genotype_input"), "catalog.tags.tsv.gz")
        })
        self.tag_generate.on("end", self.set_output, "tag_generate")
        self.tag_generate.on('end', self.run_vcf_convert)
        self.tag_generate.run()

    def run_vcf_convert(self):
        """
        :return:
        """
        self.vcf_convert.set_options({
            "vcf_path": os.path.join(os.path.join(self.work_dir, "genotype_input"), "populations.snps.vcf")
        })
        self.vcf_convert.on("end", self.set_output, "vcf_convert")
        self.vcf_convert.on('end', self.end)
        self.vcf_convert.run()

    def cp_output(self):
        """
        复制step03、04、05的结果。
        :return:
        """
        if os.path.exists(os.path.join(self.work_dir, "genotype_input")):
            shutil.rmtree(os.path.join(self.work_dir, "genotype_input"))
        os.mkdir(os.path.join(self.work_dir, "genotype_input"))
        for root, dirs, files in os.walk(self.option("ustacks_output")):
            for name in files:
                if name.endswith(".check"):
                    pass
                else:
                    os.symlink(os.path.join(root, name), os.path.join(os.path.join(self.work_dir, "genotype_input"), name))
        for root, dirs, files in os.walk(self.option("cstacks_output")):
            for name in files:
                os.symlink(os.path.join(root, name), os.path.join(os.path.join(self.work_dir, "genotype_input"), name))
        for root, dirs, files in os.walk(self.option("sstacks_output")):
            for name in files:
                if name.endswith(".check") or name.endswith(".stat"):
                    pass
                else:
                    os.symlink(os.path.join(root, name), os.path.join(os.path.join(self.work_dir, "genotype_input"), name))

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
        if event['data'] == 'tag_generate':
            self.linkdir(obj.output_dir, 'tag_generate')
        elif event['data'] == 'vcf_convert':
            self.linkdir(obj.output_dir, 'vcf_convert')
        elif event['data'] == 'genotype':
            self.linkdir(obj.output_dir, 'genotype')
        else:
            pass

    def run(self):
        super(GenotypeModule, self).run()
        self.cp_output()
        self.run_genotype()

    def end(self):
        super(GenotypeModule, self).end()
