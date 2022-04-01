# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# edit by yuguo

"""WGS基础工作流"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from bson.objectid import ObjectId
import os
import datetime
import json
import time
import shutil


class WgsWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        version = 1.0.0
        lasted modifed by HONGDONG 20180418
        写在前面：所有的样本名必须是下划线，不能有中划线（A8_10），有不同批次的时候（A8_10-1）
        """
        self._sheet = wsheet_object
        super(WgsWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'in_fastq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq文件夹
            # {'name': 'speciman_info', 'type': 'string'},  # 样本集信息表，后面改下file格式
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},  # 分组方案表--我们这边导入
            {"name": 'species_id', 'type': 'string'},  # 物种id
            {"name": "genome_version_id", "type": "string"},  # 基因组版本id  2015-11-04||GCF_000004515.4
            {"name": "snp_indel_method", "type": "string", "default": "gatk"},  # call snp indel 的方法：gatk、samtools、freebayes
            {"name": "cnv_method", "type": "string", "default": ""},  # call cnv 的方法, 为空就不进行该分析
            {"name": "sv_method", "type": "string", "default": ""},  # call sv 的方法，为空就不进行该分析
            # {"name": "transgenic_file", "type": "string", "default": ""}  # 转基因分析， 为空就不进行该分析
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.api_base = self.api.api("wgs.api_base")   # 导表的基础模块
        self.wgs_base = self.api.api('wgs.wgs_base')
        self.json_path = self.config.SOFTWARE_DIR + "/database/dna_wgs_geneome/"  # 参考组配置文件
        self.qc_stat = self.add_module("wgs.qc_stat")  # qc_质控
        self.logger.info(self.qc_stat)
        self.mapping = self.add_module("wgs.mapping")  # 样本mapping
        self.mapping_stat = self.add_module("wgs.mapping_stat")  # 样本mapping结果统计
        self.call_snp_indel = self.add_module("wgs.call_snp_indel")  # call_snp_indel
        self.cnv_call = self.add_module("wgs.cnv_call")  # cnv_call
        self.sv_call = self.add_module("wgs.sv_call")  # sv_call
        self.vcf_filter = self.add_module("wgs.vcf_filter")  # vcf_filter
        self.annovar = self.add_module("wgs.annovar")  # annovar
        self.circos = self.add_module("wgs.draw_circos")
        self.step.add_steps("qc_stat", "mapping", "call_snp_indel", "vcf_filter", "call_cnv", "call_sv", "annovar")
        self.target_dir = ""  # 文件上传到磁盘对应的路径，所有的不导表但是需要用于后
        # 面计算都保存磁盘的路径
        # self.sample_info = ""
        self.ref_path = ""
        self.ref_dict = ""
        self.snpeff_path = ""
        self.ref = ""
        self.gff = ""
        self.anno = ""
        self.software = ''
        self.is_cnvnator = False
        self.region = ""
        self.call_type = ''
        self.large_genome = "false"  # 如果是大基因组，并且是wgs的测试方式，要使用samtools mkdup去进行去pcr重复

    def check_options(self):
        """
        检查参数设置--样本信息与分组方案要重新写下
        """
        if not self.option("in_fastq").is_set:
            raise OptionError("必须要输入原始fastq文件夹路径！", code="14501301")
        # if not self.option("speciman_info").is_set:
        #     raise OptionError("必须要输入样本信息表！")
        # elif self.option("group").is_set:
        #     raise OptionError("必须要输入分组方案文件！")
        # if self.option("species_name"):
        #     raise OptionError("必须要输入物种名字！")
        if self.option("snp_indel_method").lower() not in ['gatk', 'samtools', 'freebayes']:
            raise OptionError("snp_indel_method必须为'gatk', 'samtools', 'freebayes'三者之一!", code="14501302")
        return True

    def set_ref_info(self):
        """
        初始化ref，gff，ref_chrlist等信息
        :return:
        """
        self.ref, self.gff, self.anno, self.ref_dict, self.snpeff_path, ref_chrlist, ssr_path, ref_changelog, ref_log = \
            self.wgs_base.set_ref(self.json_path, self.option("genome_version_id"))
        self.ref_path = ssr_path
        os.system("mkdir " + self.output_dir + "/02.reference")
        os.system("mkdir " + self.output_dir + "/08.ssr")
        os.system("cp " + ref_chrlist + " " + self.output_dir + "/02.reference")
        os.system("cp " + self.ref + " " + self.output_dir + "/02.reference")
        os.system("cp " + ref_log + " " + self.output_dir + "/02.reference")
        os.system("cp " + ref_changelog + " " + self.output_dir + "/02.reference")
        os.system("cp " + ssr_path + "/ssr.stat" + " " + self.output_dir + "/08.ssr")
        os.system("cp " + ssr_path + "/ssr.ref.result.xls" + " " + self.output_dir + "/08.ssr")
        os.system("cp " + ssr_path + "/ref.2bit" + " " + self.output_dir + "/02.reference")
        if not os.path.exists(self.output_dir + "/02.reference/ref.gff"):
            os.link(ssr_path + "/ref.gff", self.output_dir + "/02.reference/ref.gff")
        if not os.path.exists(self.output_dir + "/02.reference/ssr.ref.result.xls"):
            os.link(ssr_path + "/ssr.ref.result.xls", self.output_dir + "/02.reference/ssr.ref.result.xls")
        if not os.path.exists(self.output_dir + "/02.reference/anno.summary"):
            os.link(ssr_path + "/anno.summary", self.output_dir + "/02.reference/anno.summary")
        self.is_cnvnator_set(ssr_path + "/total.chrlist")
        self.set_call_type(ref_chrlist)

    def run_qc_stat(self):
        """
        计算所有fastq文件的qc相关文件，并进行统计
        :return:
        """
        self.qc_stat.set_options({
            "fastq_dir": self.option("in_fastq"),
            "task_id": self._sheet.id
        })
        self.qc_stat.on("end", self.set_output, "qc_stat")
        self.qc_stat.on("start", self.set_step, {'start': self.step.qc_stat})
        self.qc_stat.on("end", self.set_step, {'end': self.step.qc_stat})
        self.qc_stat.run()

    def run_mapping(self):
        self.mapping.set_options({
            "fastq_list": self.qc_stat.output_dir + "/clean_data/fastq.list",  # 质控后的fastq列表
            "ref_fa": self.ref,
            "large_genome": self.large_genome
        })
        self.mapping.on("end", self.set_output, "mapping")
        self.mapping.on("start", self.set_step, {'start': self.step.mapping})
        self.mapping.run()

    def run_mapping_stat(self):
        self.mapping_stat.set_options({
            "bam_list": self.mapping.output_dir + "/bam.list",
            "ref_dict": self.ref_dict,
            "step_num": 10000
        })
        self.mapping_stat.on("end", self.set_output, "mapping_stat")
        self.mapping_stat.on("end", self.set_step, {'end': self.step.mapping})
        self.mapping_stat.run()

    def run_call_snp_indel(self):
        self.call_snp_indel.set_options({
            "bam_list": self.mapping.output_dir + "/bam.list",  # bam文件列表，不同call snp方法对应的bam文件不同
            "ref_dict": self.ref_dict,
            "ref_fasta": self.ref,
            "call_type": self.option("snp_indel_method").lower()  # 这里也要转成小写
        })
        self.call_snp_indel.on("end", self.set_output, "snp_indel")
        self.call_snp_indel.on("start", self.set_step, {'start': self.step.call_snp_indel})
        self.call_snp_indel.on("end", self.set_step, {'end': self.step.call_snp_indel})
        self.call_snp_indel.run()

    def run_vcf_filter(self):
        self.vcf_filter.set_options({
            "ref_fasta": self.ref,
            "pop_var_vcf": self.call_snp_indel.output_dir + "/vcf_call/pop.variant.vcf",
            "snpEff_config": self.snpeff_path
        })
        self.vcf_filter.on("end", self.set_output, "vcf_filter")
        self.vcf_filter.on("start", self.set_step, {'start': self.step.vcf_filter})
        self.vcf_filter.on("end", self.set_step, {'end': self.step.vcf_filter})
        self.vcf_filter.run()

    def run_annovar(self):
        self.annovar.set_options({
            "snp_anno_primary_vcf": self.vcf_filter.output_dir + "/eff/snp.anno.primary.vcf",
            "indel_anno_primary_vcf": self.vcf_filter.output_dir + "/eff/indel.anno.primary.vcf",
            "ref_fasta": self.ref,
            "snp_anno_genes": self.vcf_filter.output_dir + "/eff/snp.anno.genes.txt",
            "indel_anno_genes": self.vcf_filter.output_dir + "/eff/indel.anno.genes.txt",
            "anno_summary": self.anno
        })
        self.annovar.on("end", self.set_output, "annovar")
        self.annovar.on("start", self.set_step, {'start': self.step.annovar})
        self.annovar.on("end", self.set_step, {'end': self.step.annovar})
        self.annovar.run()

    def run_cnv_call(self):
        self.cnv_call.set_options({
            "bam_list": self.mapping.output_dir + "/bam.list",  # bam文件路径list.txt文件，第一列样本名，第二列bam路径
            "ref_gff": self.gff
        })
        self.cnv_call.on("end", self.set_output, "cnv_call")
        self.cnv_call.on("start", self.set_step, {'start': self.step.cnv_call})
        self.cnv_call.on("end", self.set_step, {'end': self.step.cnv_call})
        self.cnv_call.run()

    def run_sv_call(self):
        self.sv_call.set_options({
            "bam_list": self.mapping.output_dir + "/bam.list",  # bam文件路径list.txt文件，第一列样本名，第二列bam路径
            "ref_gff": self.gff
        })
        self.sv_call.on("end", self.set_output, "sv_call")
        self.sv_call.on("start", self.set_step, {'start': self.step.sv_call})
        self.sv_call.on("end", self.set_step, {'end': self.step.sv_call})
        self.sv_call.run()

    def run_circos(self):
        self.circos.set_options({
            "genome_version_id": self.option("genome_version_id"),
            "snp": self.vcf_filter.output_dir + "/eff/snp.anno.primary.vcf",
            "indel": self.vcf_filter.output_dir + "/eff/indel.anno.primary.vcf",
            "sv_path": self.sv_call.output_dir + "/anno",
            "cnv_path": self.cnv_call.output_dir + "/anno"
        })
        self.circos.on("end", self.set_output, "circos")
        self.circos.run()

    def run(self):
        # self.start_listener()
        # self.fire("start")
        task_list = [self.mapping_stat, self.annovar]
        # task_list = [self.annovar]
        self.get_target_dir()  # 初始化一下 获取到远程磁盘的路径，用于保存对应文件的路径到主表中
        self.wgs_base.add_sg_task(self._sheet.member_id, self._sheet.member_type, self._sheet.cmd_id)
        self.set_ref_info()
        # hongdong测试修改
        self.logger.info(self.qc_stat)
        self.qc_stat.on('end', self.run_mapping)
        self.mapping.on("end", self.run_mapping_stat)
        self.mapping.on("end", self.run_call_snp_indel)
        self.call_snp_indel.on("end", self.run_vcf_filter)
        self.vcf_filter.on("end", self.run_annovar)
        if self.option("cnv_method") and self.is_cnvnator:
            self.mapping.on("end", self.run_cnv_call)
            task_list.append(self.cnv_call)
            self.step.add_steps('cnv_call')
            # self.run_cnv_call()  # 测试用，需删除
        if self.option("sv_method"):
            self.mapping.on("end", self.run_sv_call)
            task_list.append(self.sv_call)
            self.step.add_steps('sv_call')
            # self.run_sv_call()  # 测试用，需删除
        if self.option("cnv_method") and self.option("sv_method") and self.is_cnvnator:
            task_list = [self.mapping_stat, self.annovar, self.circos]
            self.on_rely([self.vcf_filter, self.cnv_call, self.sv_call], self.run_circos)
        self.on_rely(task_list, self.end)
        self.run_qc_stat()
        super(WgsWorkflow, self).run()
        # self.start_listener()
        # self.fire("start")
        # self.end()

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def move2outputdir(self, olddir, newname):
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下
        :param olddir: 初始路径
        :param newname: 目标路径，可以自定义
        :return:
        """
        start = time.time()
        if not os.path.isdir(olddir):
            self.set_error('需要移动到output目录的文件夹不存在。', code="14501303")
        newdir = os.path.join(self.output_dir, newname)
        if not os.path.exists(newdir):
            os.makedirs(newdir)
        allfiles = os.listdir(olddir)
        oldfiles = [os.path.join(olddir, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        # self.logger.info(newfiles)
        for newfile in newfiles:
            if os.path.isfile(newfile) and os.path.exists(newfile):
                os.remove(newfile)
            elif os.path.isdir(newfile) and os.path.exists(newfile):
                shutil.rmtree(newfile)
        for i in range(len(allfiles)):
            self.move_file(oldfiles[i], newfiles[i])
        end = time.time()
        duration = end - start
        self.logger.info("文件夹{}到{}移动耗时{}s".format(olddir, newdir, duration))

    def move_file(self, old_file, new_file):
        """
        递归移动文件或者文件到指定路径
        :param old_file: 初始路径
        :param new_file: 目的路径
        :return:
        """
        if os.path.isfile(old_file):
            os.link(old_file, new_file)
        else:
            os.mkdir(new_file)
            for file_ in os.listdir(old_file):
                file_path = os.path.join(old_file, file_)
                new_path = os.path.join(new_file, file_)
                self.move_file(file_path, new_path)

    def set_output(self, event):
        obj = event["bind_object"]
        if event['data'] == "qc_stat":
            self.move2outputdir(obj.output_dir, self.output_dir + "/01.fastq_qc")
        if event['data'] == "mapping":
            self.move2outputdir(obj.output_dir + '/sort_bams', self.output_dir + "/03.map_stat/map_bam/sort_bams")
        if event['data'] == "mapping_stat":
            self.move2outputdir(obj.output_dir, self.output_dir + "/03.map_stat")
        if event['data'] == "snp_indel":
            self.move2outputdir(obj.output_dir, self.output_dir + "/04.snp_indel/vcf_call")
        if event['data'] == "vcf_filter":
            self.move2outputdir(obj.output_dir, self.output_dir + "/04.snp_indel")
        if event['data'] == "annovar":
            self.move2outputdir(obj.output_dir, self.output_dir + "/05.annovar")
        if event['data'] == "cnv_call":
            self.move2outputdir(obj.output_dir, self.output_dir + "/06.cnv")
        if event['data'] == "sv_call":
            self.move2outputdir(obj.output_dir, self.output_dir + "/07.sv")
        if event['data'] == "circos":
            self.move2outputdir(obj.output_dir, self.output_dir + "/09.circos")

    def send_files(self):
        self.remove_file()
        repaths = [
            [".", "", "基础分析结果文件夹", 0, "320001"],
            ["01.fastq_qc", "", "数据质控结果目录", 0, "320002"],
            ["01.fastq_qc/rawdata_qc", "", "原始数据统计结果目录", 0, "320003"],
            ["01.fastq_qc/rawdata_qc/qc.xls", "xls", "原始数据统计表", 0, "320004"],
            ["01.fastq_qc/rawdata_qc/atgc", "", "原始数据碱基含量分布", 0, "320005"],
            ["01.fastq_qc/rawdata_qc/qual", "", "原始数据碱基错误率分布", 0, "320006"],
            ["01.fastq_qc/clean_data", "", "质控后fastq文件目录", 0, "320007"],
            ["01.fastq_qc/clean_data/fastq.list", "", "各样品clean data路径", 0, "320008"],
            ["01.fastq_qc/cleandata_qc", "", "质控后数据统计结果目录", 0, "320009"],
            ["01.fastq_qc/cleandata_qc/qc.xls", "xls", "质控后数据统计表", 0, "320010"],
            ["01.fastq_qc/cleandata_qc/atgc", "", "质控后数据碱基含量分布", 0, "320011"],
            ["01.fastq_qc/cleandata_qc/qual", "", "质控后数据碱基错误率分布", 0, "320012"],
            ["02.reference", "", "参考基因组目录", 0, "320013"],
            ["02.reference/ref.fa", "", "参考基因组fa序列", 0, "320014"],
            ["02.reference/ref.chrlist", "", "参考基因组fa文件", 0, "320015"],
            ["02.reference/ref.changelog", "", "参考基因组染色体转换对应表", 0, "320016"],
            ["02.reference/info.log", "", "参考基因组版本信息", 0, "320017"],
            ["02.reference/ref.gff", "", "参考基因组gff文件", 0, "320018"],
            ["02.reference/ref.2bit", "", "参考基因组2bit文件", 0, "320019"],
            ["02.reference/ssr.ref.result.xls", "xls", "参考基因组ssr文件", 0, "320020"],
            ["03.map_stat", "", "基因组比对结果目录", 0, "320021"],
            ["03.map_stat/map_bam", "", "基因组比对bam目录", 0, "320022"],
            ["03.map_stat/map_bam/sort_bams", "", "基因组比对sort_bam目录", 0, "320023"],
            ["03.map_stat/result.stat/Total.mapped.detail.xls", "", "比对结果统计表", 0, "320024"],
            ["03.map_stat/insert", "", "基因组比对插入片段长度结果目录", 0, "320025"],
            ["03.map_stat/depth", "", "基因组比对测序深度结果目录", 0, "320026"],
            ["03.map_stat/coverage", "", "基因组比对覆盖度结果目录", 0, "320027"],
            ["03.map_stat/result.stat", "", "比对结果统计文件夹", 0, "320028"],
            ["04.snp_indel", "", "snp与indel检测结果目录", 0, "320029"],
            ["04.snp_indel/eff", "", "SNP/InDel的eff结果目录", 0, "320030"],
            ["04.snp_indel/eff/snp.anno.primary.vcf", "", "snp功能注释vcf文件", 0, "320031"],
            ["04.snp_indel/eff/indel.anno.primary.vcf", "", "indel功能注释vcf文件", 0, "320032"],
            ["04.snp_indel/variant_stat", "", "SNP/InDel统计", 0, "320033"],
            ["04.snp_indel/anno_stat", "", "SNP/InDel功能信息统计文件夹", 0, "320034"],
            ["04.snp_indel/variant_stat/snp.stat", "", "SNP数据统计表", 0, "320035"],
            ["04.snp_indel/variant_stat/snp.GQ", "", "SNP质量分布统计表", 0, "320036"],
            ["04.snp_indel/variant_stat/snp.depth", "", "SNP深度分布统计表", 0, "320037"],
            ["04.snp_indel/variant_stat/snp.matrix", "", "样本间差异SNP统计表", 0, "320038"],
            ["04.snp_indel/variant_stat/indel.stat", "", "InDel数据统计表", 0, "320039"],
            ["04.snp_indel/variant_stat/indel.len", "", "InDel长度信息统计表", 0, "320040"],
            ["04.snp_indel/variant_stat/indel.matrix", "", "样本间差异InDel统计表", 0, "320041"],
            ["04.snp_indel/variant_stat/indel.GQ", "", "InDel质量分布统计表", 0, "320042"],
            ["04.snp_indel/variant_stat/indel.depth", "", "InDel深度分布统计表", 0, "320043"],
            ["04.snp_indel/anno_stat/snp.stat", "", "SNP功能信息统计表", 0, "320044"],
            ["04.snp_indel/anno_stat/indel.stat", "", "InDel功能信息统计表", 0, "320045"],
            ["05.annovar", "", "基因功能注释结果目录", 0, "320046"],
            ["05.annovar/combine_variants", "", "vcf合并结果目录", 0, "320047"],
            ["05.annovar/combine_variants/pop.final.vcf", "", "vcf合并文件", 0, "320048"],
            ["05.annovar/anno_count", "", "基因功能注释结果目录", 0, "320049"],
            ["05.annovar/anno_count/pop.summary", "", "基因功能注释表", 0, "320050"],
            ["05.annovar/anno_count/pop.stat.csv", "", "基因功能统计表", 0, "320051"],
            ["05.annovar/anno_count/pop.go.stat", "", "GO功能统计表", 0, "320052"],
            ["05.annovar/anno_count/pop.kegg.stat", "", "KEGG功能统计表", 0, "320053"],
            ["05.annovar/anno_count/pop.eggnog.stat", "", "EGGNOG功能统计表", 0, "320054"],
            ["05.annovar/eggnog_anno", "", "eggnog注释结果", 0, "320055"],
            ["05.annovar/eggnog_anno/pop.eggnog.final.stat.detail", "", "eggnog注释表", 0, "320056"],
            ["05.annovar/go_anno", "", "go注释结果", 0, "320057"],
            ["05.annovar/go_anno/pop.go.final.stat.detail", "", "go注释表", 0, "320058"],
            ["05.annovar/kegg_anno", "", "kegg注释结果", 0, "320059"],
            ["05.annovar/kegg_anno/pop.kegg.final.stat.detail", "", "kegg注释表", 0, "320060"],
            ["05.annovar/kegg_anno/pathway_dir", "", "kegg注释表", 0, "320061"],
            ["06.cnv", "", "cnv结果目录", 0, "320062"],
            ["06.cnv/cnv.stat.xls", "xls", "cnv统计表", 0, "320063"],
            ["06.cnv/length", "", "各样本cnv长度统计目录", 0, "320064"],
            ["06.cnv/anno", "", "各样本cnv注释目录", 0, "320065"],
            ["06.cnv/cnv", "", "cnv详情表目录", 0, "320066"],
            ["07.sv", "", "各样品SV详情表目录", 0, "320067"],
            ["07.sv/sv.stat.xls", "xls", "sv统计表", 0, "320068"],
            ["07.sv/length", "", "各样本sv长度统计目录", 0, "320069"],
            ["07.sv/anno", "", "各样本sv注释目录", 0, "320070"],
            ["07.sv/sv", "", "sv详情表目录", 0, "320071"],
            ["08.ssr", "", "ssr分析结果目录", 0, "320072"],
            ["08.ssr/ssr.stat", "", "参考基因组ssr统计表", 0, "320073"],
            ["08.ssr/ssr.ref.result.xls", "xls", "参考基因组ssr详情表", 0, "320074"],
            ["09.circos", "", "circos结果目录", 0, "320075"]
        ]

        regexps = [
            [r"01.fastq_qc/rawdata_qc/atgc/.*\.xls", "xls", "各样本原始数据碱基含量分布", 0, "320076"],
            [r"01.fastq_qc/rawdata_qc/qual/.*\.xls", "xls", "各样本原始数据碱基错误率分布",0,"320077"],
            [r"01.fastq_qc/cleandata_qc/atgc/.*\.xls", "xls", "各样本质控后数据碱基含量分布",0,"320078"],
            [r"01.fastq_qc/cleandata_qc/qual/.*\.xls", "xls", "各样本质控后数据碱基错误率分布",0,"320079"],
            [r"01.fastq_qc/clean_data/.*\.fastq\.gz", "", "样品clean data文件", 0, "320080"],
            [r"03.map_stat/insert/.*\.xls", "xls", "各样本基因组比对插入片段长度", 0, "320081"],
            [r"03.map_stat/depth/.*\.xls", "xls", "各样本基因组比对测序深度", 0, "320082"],
            [r"03.map_stat/coverage/.*\.xls", "xls", "各样本基因组比对覆盖度", 0, "320083"],
            [r"03.map_stat/insert/.*\.xls", "xls", "各样本基因组插入片段长度", 0, "320084"],
            [r"05.annovar/kegg_anno/pathway_dir/.*\.pdf", "", "KEGG pathway通路图", 0, "320085"],
            [r"05.annovar/kegg_anno/pathway_dir/.*\.png", "", "KEGG pathway通路图", 0, "320086"],
            [r"06.cnv/length/.*\.length\.xls", "xls", "各样本cnv长度统计", 0, "320087"],
            [r"06.cnv/anno/.*\.anno\.xls", "xls", "各样本cnv注释文件", 0, "320088"],
            [r"06.cnv/cnv/.*\.cnv\.xls", "xls", "各样本cnv详情表文件", 0, "320089"],
            [r"07.sv/length/.*\.length\.xls", "xls", "各样本sv长度统计文件", 0, "320090"],
            [r"07.sv/anno/.*\.anno\.xls", "xls", "各样本sv注释文件", 0, "320091"],
            [r"06.sv/sv/.*\.sv\.xls", "xls", "各样本sv详情表", 0, "320092"],
            [r"03.map_stat/map_bam/sort_bams/.*\.sort\.bam\.bai$", "", "各样本sort_bam索引文件", 0, "320093"],
            [r"03.map_stat/map_bam/sort_bams/.*\.sort\.bam$", "", "各样本sort_bam文件", 0, "320094"],
            [r"03.map_stat/map_bam/.*\.mkdup\.bai$", "", "各样本mkdup bai文件", 0, "320095"],
            [r"03.map_stat/map_bam/.*\.mkdup\.bam$", "", "各样本mkdup bam文件", 0, "320096"],
            [r"03.map_stat/map_bam/.*\.metric$", "", "各样本metric文件", 0, "320097"],
            [r"09.circos/.*\.png$", "", "各样品circos结果文件", 0, "320098"],
            [r"09.circos/.*\.pdf$", "", "各样品circos结果文件", 0, "320099"],
            [r"09.circos/.*\.svg$", "", "各样品circos结果文件", 0, "320100"]
        ]

        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        # for i in self.get_upload_files():
        #     self.logger.info('upload file:{}'.format(str(i)))

    def run_api(self):
        """
        对所有结果进行导表操作
        :return:
        """
        self.set_software()
        self.import_specimen_info()
        self.import_mapping_results()
        self.import_snp_results()
        self.import_indel_results()
        if self.option("cnv_method") and self.is_cnvnator:
            self.import_cnv_results()
        if self.option("sv_method"):
            self.import_sv_results()
        self.import_annovar_results()
        self.import_ssr_results()
        self.import_genomic_view()
        self.import_files_paths()
        self.set_samples()
        self.send_files()

    def set_software(self):
        if self.option("snp_indel_method").lower() == 'gatk':
            self.software = '{\"method\":\"GATK\"}'
        elif self.option("snp_indel_method").lower() == 'samtools':
            self.software = '{\"method\":\"samtools\"}'
        elif self.option("snp_indel_method").lower() == 'freebayes':
            self.software = '{\"method\":\"freebayes\"}'
        self.logger.info("设置软件名称成功！")

    def import_specimen_info(self):
        self.logger.info("开始进行样本信息导表")
        self.wgs_base.add_sg_specimen()
        self.logger.info("样本信息导表成功！")
        self.wgs_base.update_clean_path(self.output_dir + "/01.fastq_qc/clean_data")
        self.logger.info("更新clean fastq成功！")
        if self.option("group").is_set:
            self.logger.info("开始进行样本分组导表")
            self.wgs_base.add_sg_specimen_group(self.option("group").prop['path'])
            self.logger.info("样本分组导表成功！")
        self.logger.info("开始进行样本质控信息导表")
        self.wgs_base.add_sg_specimen_qc(self.output_dir + "/01.fastq_qc/rawdata_qc/qc.xls",
                                         self.output_dir + "/01.fastq_qc/cleandata_qc/qc.xls")
        self.logger.info("样本质控信息导表成功！")
        self.logger.info("开始进行碱基含量分布导表")
        self.wgs_base.add_qc_atgc_curve(self.output_dir + "/01.fastq_qc/rawdata_qc/atgc", "raw_reads")
        self.wgs_base.add_qc_atgc_curve(self.output_dir + "/01.fastq_qc/cleandata_qc/atgc", "clean_reads")
        self.logger.info("碱基含量分布导表成功！")
        self.logger.info("开始进行碱基错误率分布导表")
        self.wgs_base.add_qc_qual_curve(self.output_dir + "/01.fastq_qc/rawdata_qc/qual", "error_raw_reads")
        self.wgs_base.add_qc_qual_curve(self.output_dir + "/01.fastq_qc/cleandata_qc/qual", "error_clean_reads")
        self.logger.info("碱基错误率分布导表成功！")

    def import_mapping_results(self):
        self.logger.info("开始进行比对结果导表")
        mapping_id = self.wgs_base.add_sg_mapping()
        self.wgs_base.add_sg_mapping_detail(mapping_id,
                                            self.output_dir + "/03.map_stat/result.stat/Total.mapped.detail.xls")
        self.logger.info("比对结果导表成功！")
        self.logger.info("开始进行插入片段分布导表")
        self.wgs_base.insert_sg_mapping_curve(self.output_dir + "/03.map_stat/insert", mapping_id, "insert")
        self.logger.info("插入片段分布导表成功！")
        self.logger.info("开始进行测序深度分布导表")
        self.wgs_base.insert_sg_mapping_curve(self.output_dir + "/03.map_stat/depth", mapping_id, "dep")
        self.logger.info("测序深度分布导表成功！")
        self.logger.info("开始进行基因组覆盖度分布导表")
        self.wgs_base.add_sg_area_detail(mapping_id, self.output_dir + "/03.map_stat/coverage")
        self.logger.info("基因组覆盖度分布导表成功！")

    def import_snp_results(self):
        snp_api = self.api.api('wgs.snp')
        self.logger.info("开始进行snp统计导表")
        call_id = snp_api.add_sg_snp_call(self._sheet.project_sn, self._sheet.id, params=self.software)
        snp_api.add_sg_snp_call_stat(call_id, self.output_dir + "/04.snp_indel/variant_stat/snp.stat")
        self.logger.info("snp统计导表成功！")
        self.logger.info("开始进行snp质量评估导表")
        snp_api.add_snp_qc_curve(self._sheet.id, call_id,
                                 self.output_dir + "/04.snp_indel/variant_stat/snp.GQ", "snp_qc", "snp_qc")
        snp_api.add_snp_qc_curve(self._sheet.id, call_id,
                                 self.output_dir + "/04.snp_indel/variant_stat/snp.depth", "snp_depth", "snp_depth")
        self.logger.info("snp质量评估导表成功！")
        self.logger.info("开始进行snp功能注释导表")
        anno_id = snp_api.add_sg_snp_anno(self._sheet.project_sn, self._sheet.id, params='{\"method\":\"SNPEff\"}')
        # snp功能统计表
        snp_api.add_sg_snp_anno_stat(anno_id, self.output_dir + "/04.snp_indel/anno_stat/snp.stat", "annotion")
        # snp功效统计表
        snp_api.add_sg_snp_anno_stat(anno_id, self.output_dir + "/04.snp_indel/anno_stat/snp.stat", "effect")
        # snp功效与功能累加图与直方图
        snp_api.add_sg_snp_anno_bar(self._sheet.project_sn, self._sheet.id, anno_id,
                                    self.output_dir + "/04.snp_indel/anno_stat/snp.stat")
        self.logger.info("snp功能注释导表成功！")

    def import_indel_results(self):
        indel_api = self.api.api('wgs.indel')
        self.logger.info("开始进行indel统计导表")
        call_id = indel_api.add_sg_indel_call(self._sheet.project_sn, self._sheet.id, params=self.software)
        indel_api.add_sg_indel_call_stat(call_id, self.output_dir + "/04.snp_indel/variant_stat/indel.stat")
        self.logger.info("indel统计表导表成功")
        self.logger.info("开始进行indel长度分布导表")
        indel_api.add_indel_length_bar(self._sheet.id, call_id,
                                       self.output_dir + "/04.snp_indel/variant_stat/indel.len")
        self.logger.info("indel长度分布导表成功！")
        self.logger.info("开始进行indel质量评估导表")
        indel_api.add_indel_qc_curve(self._sheet.id, call_id,
                                     self.output_dir + "/04.snp_indel/variant_stat/indel.GQ", "indel_qc")
        indel_api.add_indel_qc_curve(self._sheet.id, call_id,
                                     self.output_dir + "/04.snp_indel/variant_stat/indel.depth", "indel_depth")
        self.logger.info("indel质量评估导表成功！")
        self.logger.info("开始进行indel功能注释导表")
        anno_id = indel_api.add_sg_indel_anno(self._sheet.project_sn, self._sheet.id, params='{\"method\":\"SNPEff\"}')
        # indel功能统计表
        indel_api.add_sg_indel_anno_stat(anno_id, self.output_dir + "/04.snp_indel/anno_stat/indel.stat", "annotion")
        # indel功效统计表
        indel_api.add_sg_indel_anno_stat(anno_id, self.output_dir + "/04.snp_indel/anno_stat/indel.stat", "effect")
        # indel功效与功能累加图与直方图
        indel_api.add_sg_indel_anno_bar(self._sheet.project_sn, self._sheet.id, anno_id,
                                        self.output_dir + "/04.snp_indel/anno_stat/indel.stat")
        self.logger.info("indel功能注释导表成功！")

    def import_cnv_results(self):
        cnv_api = self.api.api('wgs.cnv')
        self.logger.info("开始进行cnv统计导表")
        call_id = cnv_api.add_sg_cnv_call(self._sheet.project_sn, self._sheet.id, params='{\"method\":\"CNVnator\"}')
        cnv_api.add_sg_cnv_call_stat(call_id, self.output_dir + "/06.cnv/cnv.stat.xls")
        self.logger.info("cnv统计表导表成功！")
        self.logger.info("开始进行cnv长度统计导表")
        cnv_api.add_cnv_length_bar(self._sheet.id, call_id, self.output_dir + "/06.cnv/length")   # 导入length文件夹
        self.logger.info("cnv长度统计导表成功！")

    def import_sv_results(self):
        sv_api = self.api.api('wgs.sv')
        self.logger.info("开始进行sv统计导表")
        call_id = sv_api.add_sg_sv_call(self._sheet.project_sn, self._sheet.id, params='{\"method\":\"BreakDancer\"}')
        sv_api.add_sg_sv_call_stat(call_id, self.output_dir + "/07.sv/sv.stat.xls")
        self.logger.info("sv统计表导表成功！")
        self.logger.info("开始进行sv长度统计导表")
        sv_api.add_sv_length_bar(self._sheet.id, call_id, self.output_dir + "/07.sv/length")   # 导入length文件夹
        self.logger.info("sv长度统计导表成功！")

    def import_annovar_results(self):
        anno_api = self.api.api('wgs.anno')
        self.logger.info("开始进行基因功能注释导表")
        anno_id = anno_api.add_sg_anno(self._sheet.id, self._sheet.project_sn)
        anno_api.add_sg_anno_detail(self.output_dir + "/05.annovar/anno_count/pop.summary", anno_id,
                                    self._sheet.id, self.option("genome_version_id"))
        self.logger.info("基因功能注释导表成功！")
        self.logger.info("开始进行GO/KEGG/EGGNOG统计图表导表")
        anno_api.sg_anno_go_stat(anno_id, self.output_dir + "/05.annovar/go_anno/pop.go.final.stat.detail")
        anno_api.sg_anno_kegg_stat(anno_id, self.output_dir + "/05.annovar/kegg_anno/pop.kegg.final.stat.detail",
                                   self.output_dir + "/05.annovar/kegg_anno/pathway_dir",
                                   self._sheet.output.rstrip('/') + "/05.annovar/kegg_anno/pathway_dir")
        anno_api.sg_anno_eggnog_stat(anno_id, self.output_dir + "/05.annovar/eggnog_anno/pop.eggnog.final.stat.detail")

    def import_ssr_results(self):
        ssr_api = self.api.api("wgs.ssr_analysis")
        self.logger.info("开始进行参考基因组ssr统计导表")
        ssr_id = ssr_api.add_sg_ssr(self._sheet.project_sn, self._sheet.id)
        ssr_api.add_sg_ssr_stat(ssr_id, self.output_dir + "/08.ssr/ssr.stat")
        ssr_api.add_sg_ssr_detail_new(ssr_id, self.output_dir + "/08.ssr/ssr.ref.result.xls",
                                      self._sheet.output.rstrip('/') + "/02.reference/ssr.ref.result.xls")
        self.logger.info("参考基因组ssr统计导表成功！")

    def import_genomic_view(self):
        self.logger.info("开始进行基因组浏览器结果导表")
        self.wgs_base.add_sg_igv(self._sheet.output.rstrip('/') + "/03.map_stat/map_bam/sort_bams/",
                                 self._sheet.output.rstrip('/') +
                                 "/05.annovar/combine_variants/pop.final.vcf",
                                 self._sheet.output.rstrip('/') + "/02.reference/ref.2bit")
        self.logger.info("基因组浏览器结果导表成功！")
        if self.option("cnv_method") and self.option("sv_method") and self.is_cnvnator:
            self.logger.info("开始进行circos结果导表")
            self.wgs_base.add_sg_circos(self._sheet.output.rstrip('/') + "/09.circos/",
                                        self.output_dir + "/09.circos", self.ref_path.rstrip("/") + "/total.chrlist",
                                        params='{\"chrs\": \"\"}')
            self.logger.info("circos结果导表成功！")

    def import_files_paths(self):
        file_paths = {
            "clean_fastq_path": self.target_dir + "/01.fastq_qc/clean_data/",
            "bam_path": self.target_dir + "/03.map_stat/map_bam/sort_bams/",
            "indel_anno_vcf": self.target_dir + "/04.snp_indel/eff/indel.anno.primary.vcf",
            "snp_anno_vcf": self.target_dir + "/04.snp_indel/eff/snp.anno.primary.vcf",
            "pop_final_vcf": self.target_dir + "/05.annovar/combine_variants/pop.final.vcf",
            "cnv_anno_path": self.target_dir + "/06.cnv/anno/",
            "sv_anno_path": self.target_dir + "/07.sv/anno/",
            "genome_version_id": ObjectId(self.option("genome_version_id")),
            "region": self.region
        }
        self.api_base.update_db_record("sg_task", {"task_id": self._sheet.id}, file_paths)

    def set_samples(self):
        samples = self.wgs_base.find_all_sample()
        with open(self.output_dir + "/info.list", "w") as w:
            w.write("PID\t{}\n".format(samples))
            w.write("BID\t{}\n".format(samples))

    def end(self):
        # self.send_files()
        # self.logger.info("上传文件目录完成！")
        self.run_api()
        self.logger.info("全部导表完成！")
        super(WgsWorkflow, self).end()

    def get_target_dir(self):
        """
        获取远程磁盘的路径
        :return:
        """
        self.target_dir = self._sheet.output.rstrip('/')
        self.region = self._sheet.output.split('://')[0]

    def remove_file(self):
        """
        删除不需要上传到磁盘的文件
        :return:
        """
        rm_list = list()
        rm_list.append(self.output_dir + "/04.snp_indel/eff/indel.anno.csv")
        rm_list.append(self.output_dir + "/04.snp_indel/eff/indel.anno.genes.txt")
        rm_list.append(self.output_dir + "/04.snp_indel/eff/snp.anno.csv")
        rm_list.append(self.output_dir + "/04.snp_indel/eff/snp.anno.genes.txt")
        rm_list.append(self.output_dir + "/04.snp_indel/vcf_filter")
        rm_list.append(self.output_dir + "/04.snp_indel/vcf_call")
        rm_list.append(self.output_dir + "/01.fastq_qc/clean_data")
        for files in rm_list:
            if os.path.isfile(files):
                os.remove(files)
                self.logger.info("删除文件{}成功！".format(files))
            elif os.path.isdir(files):
                code = os.system("rm -r {}".format(files))
                if code != 0:
                    self.logger.info("删除文件夹{}失败！".format(files))
            else:
                self.logger.info("文件{}不存在，不用删除！".format(files))

    def is_cnvnator_set(self, total_chrlist):
        """
        参考基因组chr与sca小于1500的话 才会做cnvnator
        :param total_chrlist:
        :return:
        """
        if len(open(total_chrlist, 'rU').readlines()) < 50000:
            self.is_cnvnator = True
        self.logger.info("染色体个数为{}， is_cnvnator为：{}"
                         .format(len(open(total_chrlist, 'rU').readlines()), self.is_cnvnator))

    def set_call_type(self, ref_chrlist):
        """
        检查染色体是不是大于536870912，如果有一条大于536870912，直接用samtool进行call snp不能用gatk跑
        因为染色体很大的时候，就不能建立bai索引，进而导致了call snp失败（默认是使用samtools去call），
        同时没有bai索引的对于wgs实验方法，无法使用picard进行去重复，这个时候使用samtools markdup进行去重复
        :return:
        """
        self.call_type = self.option("snp_indel_method").lower()
        with open(ref_chrlist, 'r') as r:
            data = r.readlines()
            for line in data:
                temp = line.strip().split('\t')
                if int(temp[1]) > 536870912:   # 536870912 301019445
                    self.call_type = 'samtools'
                    self.large_genome = "true"   # 判断是不是大基因组参考组
                    break
        self.logger.info("call type is {}--large_genome is {}".format(self.call_type, self.large_genome))

