# -*- coding: utf-8 -*-
# __author__ = 'liuwentian'
# modified 2019.02.20


from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from biocluster.config import Config
import os, re
import math
import json
import xlrd
import time, datetime
from bson.objectid import ObjectId


class SvCallModule(Module):
    """
    sv变异统计module
    """
    def __init__(self, work_id):
        super(SvCallModule, self).__init__(work_id)
        options = [
            {"name": "ref_fa", "type": "string"},  # ref.fa
            {"name": "bam_list", "type": "infile", "format": "wgs_v2.bam_dir"},  # 最终的bam文件
            {"name": "snpEff_config", "type": "string"},  # 数据库中取
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "target_path", 'type': "string"}
        ]
        self.add_option(options)
        self.delly_call = self.add_tool("wgs_v2.delly_call")
        self.bcftools_convert = self.add_tool("wgs_v2.bcftools_convert")
        self.snpeff = self.add_tool("wgs.snpeff")
        self.sv_stat_v2 = self.add_tool("wgs_v2.sv_stat_v2")

    def check_options(self):
        if not self.option("ref_fa"):
            raise OptionError("请设置ref_fa")
        if not self.option("bam_list"):
            raise OptionError("请设置bam_list")
        # if not self.option("snpEff_config"):
        #     raise OptionError("请设置snpEff_config")
        # if not self.option("task_id"):
        #     raise OptionError("请设置task_id")

    def run_delly_call(self):
        self.delly_call.set_options({
            "ref_fa": self.option("ref_fa"),
            "bam_list": os.path.join(self.work_dir, "bam.list")
        })
        # self.delly_call.on("end", self.set_output, "delly_call")
        self.delly_call.on('end', self.run_bcftools_convert)
        self.delly_call.run()

    def run_bcftools_convert(self):
        bcf_file = self.delly_call.option("bcf_file").prop["path"]
        self.bcftools_convert.set_options({
            "bcf_file": bcf_file
        })
        self.bcftools_convert.on("end", self.set_output, "bcftools_convert")
        if self.option("snpEff_config"):
            self.bcftools_convert.on('end', self.run_snpeff)
        else:
            self.bcftools_convert.on('end', self.run_sv_stat_v2)
        self.bcftools_convert.run()

    def run_snpeff(self):
        vcf_path = self.bcftools_convert.option("vcf_file").prop["path"]
        self.snpeff.set_options({
            "ref_name": "ref",
            "snpEff_config": self.option('snpEff_config'),
            "filter_recode_vcf": vcf_path,
            "types": "sv"
        })
        self.snpeff.on("end", self.set_output, "snpeff")
        self.snpeff.on('end', self.run_sv_stat_v2)
        self.snpeff.run()

    def run_sv_stat_v2(self):
        if self.option("snpEff_config"):
            sv_vcf = os.path.join(self.snpeff.output_dir, "pop.sv.anno.vcf")
        else:
            sv_vcf = os.path.join(self.bcftools_convert.output_dir, "pop.sort.sv.vcf")
        self.sv_stat_v2.set_options({
            "sv_vcf": sv_vcf
        })
        self.sv_stat_v2.on("end", self.set_output, "sv_stat_v2")
        self.sv_stat_v2.on('end', self.end)
        self.sv_stat_v2.run()

    def make_bam_list(self):
        file_list = os.listdir(self.option('bam_list').prop['path'])
        write_lines = ""
        for i in file_list:
            if i.endswith(".bam"):
                tmp = i.strip().split(".bam")
                sample = tmp[0]
                bam_path = os.path.join(self.option('bam_list').prop['path'], i)
                write_lines += sample + "\t" + bam_path + "\n"
        with open(os.path.join(self.work_dir, "bam.list"), "w")as fw:
            fw.write(write_lines)

    def set_db(self):
        self.logger.info("开始sg_sv_call_stat的导表！")
        api = self.api.api("wgs_v2.sv_call")
        api.add_sg_sv_call_stat(self.option("main_id"), os.path.join(self.sv_stat_v2.output_dir, "stat.txt"))
        self.logger.info("设置sg_sv_call_stat的导表成功！")
        self.logger.info("开始sg_bar的导表！")
        # with open(self.option("bam_list"), "r")as fr:
        #     lines = fr.readlines()
        #     for line in lines:
        #         tmp = line.strip().split("\t")
        #         api.add_sg_sv_len(self.option("main_id"),
        #                           os.path.join(self.sv_stat_v2.output_dir, (tmp[0] + ".sv.length.txt")),
        #                           self.option("task_id"), tmp[0])
        api.add_sg_sv_len(self.option("main_id"), self.option("task_id"), self.sv_stat_v2.output_dir)
        api.update_db_record("sg_task", {"task_id": self.option("task_id")},
                             {"pop_sv_vcf": (self.option('target_path') + "/bcftools_convert/pop.sort.sv.vcf")})
        self.logger.info("设置sg_bar的导表成功！")

    def set_output(self, event):
        obj = event['bind_object']
        # if event['data'] == 'delly_call':
        #     self.linkdir(obj.output_dir, 'delly_call')
        if event['data'] == 'bcftools_convert':
            self.linkdir(obj.output_dir, 'bcftools_convert')
        elif event['data'] == 'snpEff_config':
            self.linkdir(obj.output_dir, 'snpEff_config')
        elif event['data'] == 'sv_stat_v2':
            self.linkdir(obj.output_dir, 'sv_stat_v2')
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
        super(SvCallModule, self).run()
        self.make_bam_list()
        self.run_delly_call()

    def end(self):
        if self.option("main_id"):
            self.set_db()
        self.add_upload_dir(self.output_dir)
        super(SvCallModule, self).end()
