# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# modified 2018.830


from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from biocluster.config import Config
from bson.objectid import ObjectId
import os, re
import json
import time, datetime
import gevent


class HapmapModule(Module):
    """
    群体进化，单倍体图谱接口调用程序
    """
    def __init__(self, work_id):
        super(HapmapModule, self).__init__(work_id)
        options = [
            {"name": "anno_summary_path", "type": "string"},    # 基因组的注释文件
            {"name": "main_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "recode_vcf_path", "type": "infile", "format": "dna_gmap.vcf"},
            {"name": "csv_dir", "type": "string"},     # con从主表里查询到csv_dir，然后传给modu
            {"name": "trait_list", "type": "string"},   # con从主表里查询到，然后用，分割
            {"name": "distance", "type": "int"},    # , "default": 200000
            {"name": "ld_r2", "type": "float"},     # , "default": 0.8
            {"name": "p_value", "type": "float"},   # , "default": 0.05
            {"name": "q_value", "type": "float"},   # , "default": 0.05
            {"name": "update_info", "type": "string"},
            {"name": "is_wgs_result", "type": "string"}
        ]
        self.add_option(options)
        self.vcftools_plink = self.add_tool("dna_evolution.vcftools_plink")
        self.region_hapmap_tools = []
        # self.region_hapmap = self.add_tool("dna_evolution.region_hapmap")

    def check_options(self):
        if not self.option("main_id"):
            raise OptionError("请设置main_id")
        if not self.option("task_id"):
            raise OptionError("请设置task_id")
        # if not self.option("pop_dir"):
        #     raise OptionError("请设置pop_dir")
        if not self.option("recode_vcf_path").is_set:
            raise OptionError("请设置recode_vcf_path")
        if not self.option("anno_summary_path"):
            raise OptionError("请设置anno_summary_path")
        if not self.option("csv_dir"):
            raise OptionError("请设置csv_dir")
        if not self.option("trait_list"):
            raise OptionError("请设置trait_list")
        if not self.option('p_value') and not self.option('q_value'):
            raise OptionError('p:{}/q_value:{}有且只能输入一个'.format(self.option('p_value'), self.option('q_value')))
        if not self.option('distance') and not self.option('ld_r2'):
            raise OptionError('distance:{}/ld_r2:{}有且只能输入一个'.format(self.option('distance'), self.option('ld_r2')))
        if not self.option("is_wgs_result"):
            raise OptionError("请设置is_wgs_result")

    def check_isfile(self, newfile):
        if os.path.exists(newfile):
            return newfile
        else:
            self.set_error("文件不存在：{}".format(newfile))

    # def get_csv_file(self):
    #     self.api = self.api.api("dna_gmap.hapmap")
    #     tuple_ = self.api.get_csv_file("查询主表的trait_list和csv_path")
    #     self.csv_path = tuple_[0]
    #     self.trait_list = tuple_[1]

    def run_vcftools_plink(self):       # 1
        """
        vcftools --vcf --plink --out
        flag为false,不跑
        """
        self.vcftools_plink.set_options({
            "recode_vcf_path": self.option("recode_vcf_path").prop["path"],
            "flag": False
        })
        self.vcftools_plink.on('end', self.run_region_hapmap)
        self.vcftools_plink.run()
        print("\n##### vcftools_plink\n")

    def run_region_hapmap(self):
        """
        循环tool,先核查pop_dir路径下的文件
        MVP.TRAIT3.FarmCPU.csv  MVP.TRAIT3.GLM.csv  MVP.TRAIT3.MLM.csv
        """
        self.check_isfile(os.path.join(self.vcftools_plink.output_dir, "pop.map"))
        self.check_isfile(os.path.join(self.vcftools_plink.output_dir, "pop.ped"))
        pop_path = self.vcftools_plink.output_dir + "/pop"
        self.logger.info('@@@@________pop文件核查通过')
        trait_list = self.option("trait_list").split(',')
        for trait in trait_list:    # 没用遍历的方法，需确定每一个文件是存在的
            csv1 = self.check_isfile(os.path.join(self.option('csv_dir'),
                                                  os.path.join(trait, "MVP." + trait + ".FarmCPU.csv")))
            csv2 = self.check_isfile(os.path.join(self.option('csv_dir'),
                                                  os.path.join(trait, "MVP." + trait + ".GLM.csv")))
            csv3 = self.check_isfile(os.path.join(self.option('csv_dir'),
                                                  os.path.join(trait, "MVP." + trait + ".MLM.csv")))
            for csv_path in [csv1, csv2, csv3]:
                region_hapmap = self.add_tool("dna_evolution.region_hapmap")
                options = {
                    "anno_summary_path": self.option("anno_summary_path"),
                    "recode_vcf_path": self.option("recode_vcf_path").prop["path"],
                    "csv_path": csv_path,
                    "pop_path": pop_path,
                    "task_id": self.option("task_id"),
                    "main_id": self.option("main_id"),
                    "is_wgs_result": self.option("is_wgs_result"),
                    "file_path": self.parent._sheet.output
                }
                for opt in ['distance', 'ld_r2', 'p_value', 'q_value']:
                    if self.option(opt):
                        options[opt] = self.option(opt)
                region_hapmap.set_options(options)
                self.logger.info(trait)
                self.region_hapmap_tools.append(region_hapmap)
        # if self.region_hapmap_tools:
        #     if len(self.region_hapmap_tools) > 1:
        #         self.on_rely(self.region_hapmap_tools, self.hardlink)
        #     elif len(self.region_hapmap_tools) == 1:
        #         self.region_hapmap_tools[0].on('end', self.hardlink, 'trit_mvp')
        # else:
        #     raise Exception("self.region_hapmap_tools列表为空！")
        for j in range(len(self.region_hapmap_tools)):
            self.region_hapmap_tools[j].on('end', self.set_output, 'trit_mvp')
        if len(self.region_hapmap_tools) > 1:
            self.on_rely(self.region_hapmap_tools, self.set_db)
        elif len(self.region_hapmap_tools) == 1:
            self.region_hapmap_tools[0].on('end', self.set_db)
        else:
            self.set_error("region_hapmap_tools列表为空！")
        for tool in self.region_hapmap_tools:
            gevent.sleep(0)
            tool.run()

    def get_header(self):
        """
        得到pop.recode.vcf里是否有indel
        """
        with open(self.option("recode_vcf_path").prop["path"], 'r') as f:
            line = f.readline()
            while line:
                if line.startswith('#'):
                    pass
                else:
                    lines = line.strip().split('\t')
                    list_ = ','.join([lines[3], lines[4]])
                    typ = list_.split(',')
                    for i in typ:
                        if len(i) > 1:
                            return 2
                line = f.readline()
            return 1

    def linkdir(self, dirpath, dirname):
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

    def set_output(self, event):
        obj = event['bind_object']
        if not os.path.exists(self.output_dir + "/region_dir"):
            os.mkdir(self.output_dir + "/region_dir")
        if event['data'] == 'trit_mvp':
            # self.checklink(os.path.join(self.output_dir, "ssr.stat"))
            self.linkdir(obj.output_dir, self.output_dir + "/region_dir")
        # self.end()

    def checklink(self, newfile):
        if os.path.exists(newfile):
            os.remove(newfile)
        else:
            pass
        return True

    def run(self):
        super(HapmapModule, self).run()
        # self.logger.info(self.parent._sheet.output)
        self.run_vcftools_plink()
        # self.run_region_hapmap()
        # self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        super(HapmapModule, self).end()

    def set_db(self):
        """
        单倍体图谱导表
        """
        header = self.get_header()
        hapmap_id = ObjectId(self.option("main_id"))
        task_id = self.option('task_id')
        self.logger.info("开始sg_hapmap的导表！")
        # self.logger.info(self.parent._sheet.output)
        api = self.api.api("dna_evolution.hapmap")
        api.add_sg_hapmap_stat(hapmap_id, header, self.output_dir + "/region_dir")
        api_base = self.api.api("dna_evolution.api_base")
        if header == 1:
            api_base.update_db_record("sg_hapmap", {"_id": hapmap_id}, {"header": 'snp',
                                                                        "region_path": self.output_dir + "/region_dir"})
        else:
            api_base.update_db_record("sg_hapmap", {"_id": hapmap_id},
                                  {"header": 'snp,indel', "region_path": self.output_dir + "/region_dir"})
        self.logger.info("设置sg_hapmap的导表成功！")
        self.end()
        # 需要更新主表sg_hapmap的header和region_path
