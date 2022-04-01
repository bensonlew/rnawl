# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'

"""遗传进化基础工作流"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.dna_evolution.send_email import SendEmail
from bson.objectid import ObjectId
import os
import re
import datetime
import json
import time
import shutil
import gevent


class EvolutionWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        version = 1.0.0
        lasted modifed by HONGDONG 20180918
        写在前面：所有的样本名必须是下划线，不能有中划线（A8_10），有不同批次的时候（A8_10-1）
        """
        self._sheet = wsheet_object
        super(EvolutionWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'in_fastq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq文件夹
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},  # 分组方案表--我们这边导入
            # {"name": 'species_id', 'type': 'string'},  # 物种id
            {"name": "genome_version_id", "type": "string"},  # 基因组版本ObjectID
            {"name": "snp_indel_method", "type": "string", "default": "gatk"},  # gatk、samtools、freebayes
            {"name": "data_type", "type": "string", "default": "raw_path"},  # 数据来源，原始fastq还是wgs路径
            # {"name": "custom_genome", "type": "string", "default": "false"},  # 是否自定义参考基因组
            {"name": "genome_info", "type": "infile", "format": "dna_evolution.genome_path"},  # 参考基因组的文件目录
            {"name": "trait_file", "type": "infile", "format": "dna_gmap.trait"},  # 性状文件
            {"name": "is_gwas", "type": "string"},  # 是个变量名称，其中包含了群体结构，连锁不平衡等接口，为了与之前的兼容，因此保留这个名字
            {'name': 'wgs_path', 'type': 'infile', 'format': 'bsa.bsa_path'},  # 输入wgs分析的路径
            {'name': "email_id", "type": "string"}  # 用于定位项目对应的前端的mysql的id
        ]
        self.add_option(options)
        self.is_structure, self.is_ld_analysis, self.is_sweep, self.is_pop_history, self.is_annovar, self.is_gwas =\
            "", "", "", "", "", ""
        self.set_options(self._sheet.options())
        self.project_type = "dna_evolution"
        self.api_base = self.api.api("dna_evolution.api_base")   # 导表的基础模块
        self.evolution_base = self.api.api('dna_evolution.evolution_base')
        self.json_path = self.config.SOFTWARE_DIR + "/database/dna_geneome/" if self.option("data_type") == "raw_path" \
            else self.config.SOFTWARE_DIR + "/database/dna_wgs_geneome/"  # 参考组配置文件
        self.genome_config = self.add_module("wgs.genome_config_v2")
        self.qc_stat = self.add_module("wgs.qc_stat")  # qc_质控
        self.mapping = self.add_module("wgs.mapping")  # 样本mapping
        self.mapping_stat = self.add_module("wgs.mapping_stat")  # 样本mapping结果统计
        self.call_snp_indel = self.add_module("wgs.call_snp_indel")  # call_snp_indel
        self.vcf_filter = self.add_module("wgs.vcf_filter")  # vcf_filter
        self.annovar = self.add_module("wgs.annovar")  # annovar
        self.vcftools_filter = self.add_tool("dna_evolution.vcftools_filter")
        self.structure = self.add_module('dna_evolution.pop_structure')
        self.pca = self.add_tool('dna_evolution.gcta_pca')
        self.tree = self.add_module('dna_evolution.tree_generic')
        self.ld_analysis = self.add_module('dna_evolution.ld_analysis')
        self.sweep_analysis = self.add_module('dna_evolution.sweep_analysis')
        self.sweep_region = self.add_module('dna_evolution.sweep_region')
        self.pop_history = self.add_module('dna_evolution.pop_history')
        self.gwas_analysis = self.add_module('dna_evolution.gwas_analysis')
        self.hapmap_analysis = self.add_module('dna_evolution.hapmap')
        self.var_diff = self.add_tool("dna_evolution.variant_compare")
        self.go_summary = self.add_tool("dna_evolution.go_summary")
        self.vcftools_plink_tools = []
        self.step.add_steps("vcftools_filter", "structure", "tree", "pca", "ld_analysis", "sweep_analysis",
                            "pop_history", "sweep_region", 'variant_compare')
        self.target_dir = ""  # 文件上传到磁盘对应的路径
        self.base_target_path, self.ref_dict, self.snpeff_path, self.ref, self.gff, self.anno, self.software, self.\
            species_path, self.wgs_task_id= '', '', '', '', '', '', '', '', ''
        self.has_chromosomelist = True  # 默认是有改名文件的
        self.has_genome = True  # 默认是有参考基因组的
        self.chr_set = 0
        self.has_wgs_id = False
        self.genome_version_id = self.option("genome_version_id") if self.option("genome_version_id") else ''
        self.vcf_file, self.group_file, self.pop_final_vcf, self.species, self.genome_version,  = '', '', '', '', ''
        self.api_base_path = ""  # 这个用于后面导表的时候，初始化wgs路径还是原始数据分析后的路径
        self.secret = False
        self.call_type = 'gatk'
        self.wgs_api = self.api.api("wgs.api_base")
        self.group_dict = {}
        self.vcf_id, self.group_id, self.anno_summary, self.anno_id, self.call_id, self.wgs_path, self.region = \
            '', '', '', '', '', '', ''
        self.is_wgs_result = 'no' if self.option("data_type") == "raw_path" else 'yes'
        self.large_genome = "false"  # 如果是大基因组，并且是wgs的测试方式，要使用samtools mkdup去进行去pcr重复

    def check_options(self):
        """
        检查参数设置
        """
        params = self.option("is_gwas").strip().split("|")
        self.is_structure = "true" if "is_structure" in params else "false"
        self.is_ld_analysis = "true" if "is_ld_analysis" in params else "false"
        self.is_sweep = "true" if "is_sweep" in params else "false"
        self.is_pop_history = "true" if "is_pop_history" in params else "false"
        self.is_annovar = "true" if "is_annovar" in params else "false"
        self.is_gwas = "true" if "is_gwas" in params else "false"
        self.logger.info(self.is_gwas)
        if self.option("data_type") == 'raw_path':
            if not self.option("in_fastq").is_set:
                raise OptionError("必须要输入原始fastq文件夹路径！")
            if not self.option("genome_info").is_set and not self.option("genome_version_id"):
                raise OptionError("genome_info或genome_version_id必须存在一个！")
        else:
            if not self.option("wgs_path").is_set:
                raise OptionError('缺少wgs结果目录文件夹')
        if not self.option("trait_file").is_set:
            if self.is_gwas == 'true':
                raise OptionError("缺少性状文件，无法进行gwas的分析！")
        if self.option("snp_indel_method").lower() not in ['gatk', 'samtools', 'freebayes']:
            raise OptionError("snp_indel_method必须为'gatk', 'samtools', 'freebayes'三者之一!")
        if not self.option("group") and self.is_structure == "false":
            raise OptionError("分组方案和群体结构分析必须选择其中一个！")
        return True

    def set_ref_info(self, is_update=False):
        """
        初始化ref，gff，ref_chrlist等信息
        is_update默认是不更新基因组
        :return:
        """
        if is_update:
            wgs_refinfo = self.api.api('wgs.ref_info')
            wgs_refinfo.project_type = "dna_evolution"
            wgs_refinfo.add_sg_genome(self.genome_config.output_dir + "/wgs_genome.json")
            base_path = self.json_path + self.species_path
            self.ref = base_path + '/ref.fa'
            self.gff = base_path + "/ref.gff"
            self.anno = base_path + "/anno.summary"
            self.ref_dict = base_path + "/ref.dict"
            self.snpeff_path = base_path + "/snpEff.config"
            ref_chrlist = base_path + "/ref.chrlist"
            ssr_path = base_path
            ref_changelog = base_path + "/ref.changelog"
            ref_log = base_path + "/info.log"
            self.chr_set = self.api_base.get_file_len(base_path + '/total.chrlist')
            self.anno_summary = base_path + "/anno.summary"
            self.genome_version_id = wgs_refinfo.find_genome_id(self.species, self.genome_version)
            if self.secret:
                wgs_refinfo.update_secret(self.genome_version_id)
        else:
            if self.is_wgs_result == 'yes':
                self.api.del_api("dna_evolution.evolution_base")
            evolution_base = self.api.api('dna_evolution.evolution_base')
            if self.is_wgs_result == 'yes':
                evolution_base.project_type = "dna_wgs"
            self.ref, self.gff, self.anno, self.ref_dict, self.\
                snpeff_path, ref_chrlist, ssr_path, ref_changelog, ref_log, self.chr_set, self.anno_summary = \
                evolution_base.set_ref(self.json_path, self.genome_version_id)
            if self.is_wgs_result == 'yes':
                self.api.del_api("dna_evolution.evolution_base")
        self.logger.info("chr_set:{}".format(self.chr_set))
        os.system("mkdir " + self.output_dir + "/02.reference")
        os.system("cp " + ref_chrlist + " " + self.output_dir + "/02.reference")
        os.system("cp " + self.ref + " " + self.output_dir + "/02.reference")
        os.system("cp " + ref_log + " " + self.output_dir + "/02.reference")
        os.system("cp " + ref_changelog + " " + self.output_dir + "/02.reference")
        self.os_link(ssr_path.rstrip('/') + "/ref.gff", self.output_dir + "/02.reference/ref.gff")
        self.os_link(ssr_path.rstrip('/') + "/anno.summary", self.output_dir + "/02.reference/anno.summary")
        self.set_call_type(ref_chrlist)

    def os_link(self, source, target):
        if not os.path.exists(target):
            os.link(source, target)
        self.logger.info("移动{}到{}成功！".format(source, target))

    def run_genome_config(self):
        options = {
            'reffa': self.option("genome_info").prop['path'].rstrip('/') + "/ref.fa",
            "refgff": self.option("genome_info").prop['path'].rstrip('/') + "/ref.gff",
            "info": self.option("genome_info").prop['path'].rstrip('/') + "/info.log"
        }
        if True:
            options.update({
                "need_rename": 'false',
                "ref_changelog": self.option("genome_info").prop['path'].rstrip('/') + "/ref.changelog"
            })
        if self.has_chromosomelist:
            options.update({"chromosomelist": self.option("genome_info").prop['path'].rstrip('/') + "/chr.list"})
        self.genome_config.set_options(options)
        self.genome_config.on("start", self.set_step, {'start': self.step.genome_config})
        self.genome_config.on("end", self.set_step, {'end': self.step.genome_config})
        self.genome_config.on("end", self.set_output, "genome_config")
        self.genome_config.run()

    def run_qc_stat(self):
        """
        计算所有fastq文件的qc相关文件，并进行统计
        :return:
        """
        self.qc_stat.set_options({
            "fastq_dir": self.option("in_fastq"),
            "task_id": self._sheet.id,
            "project_type": "dna_evolution"
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
            "call_type": self.call_type  # 这里也要转成小写
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
        """
        Annovar/GatkCombineVariants/output/pop.final.vcf
        :return:
        """
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

    def run_go_summary(self):
        if self.option("data_type") == "raw_path":
            pop_summary = self.annovar.output_dir + "/anno_count/pop.summary"
        else:
            pop_summary = self.option("wgs_path").prop['path'] + "/05.annovar/anno_count/pop.summary"
        options = {
            "pop_summary": pop_summary
        }
        self.go_summary.set_options(options)
        self.go_summary.on('end', self.set_output, "go_summary")
        self.go_summary.run()

    def run_variant_compare(self):
        self.make_compare_config()
        options = {
            "filter_recode_vcf": self.vcf_file,
            "variant_compare_config": os.path.join(self.work_dir, "config.txt"),
            "group_table": os.path.join(self.work_dir, "compare_group.txt")
        }
        self.var_diff.set_options(options)
        self.var_diff.on("end", self.set_output, "variant_compare")
        self.var_diff.on("start", self.set_step, {'start': self.step.variant_compare})
        self.var_diff.on("end", self.set_step, {'end': self.step.variant_compare})
        self.var_diff.run()

    def run_vcftools_filter(self):
        """
        vcf进行过滤
        :return:
        """
        options = {
            "vcf_path": self.vcf_file,
            "recode": True,
            "remove_indels": False,
            "remove_filtered_all": False,
            "minDP": 1,
            "maxDP": 300000000,
            "max_missing": 0.3,
            "min_maf": 0.05,
            "max_maf": 1
        }
        self.vcftools_filter.set_options(options)
        self.vcftools_filter.on("end", self.set_output, "vcftools_filter")
        if self.is_structure == "true":
            self.vcftools_filter.on("end", self.run_tree)
            self.vcftools_filter.on("end", self.run_vcftools_plink)
        self.vcftools_filter.on("start", self.set_step, {'start': self.step.vcftools_filter})
        self.vcftools_filter.on("end", self.set_step, {'end': self.step.vcftools_filter})
        self.vcftools_filter.run()

    def run_vcftools_plink(self):
        """
        这个地方写的有点不好，后面有时间的话可以将该步骤加入到对应structure与pca的module中去
        当i == 1的时候，是pop中structure的计算输入，否则就是pca的输入
        :return:
        """
        for i in range(0, 2):
            vcftools_plink = self.add_tool("dna_evolution.vcftools_plink")
            options = {
                "recode_vcf_path": self.vcftools_filter.option("filter_vcf").prop['path'],
                "chr_set": self.chr_set,  # 改物种具体的染色体与sca的个数
                "make_bed": True
            }
            if i == 1:
                options.update({"allow_extra_chr": True})
            vcftools_plink.set_options(options)
            self.vcftools_plink_tools.append(vcftools_plink)
        self.vcftools_plink_tools[0].on('end', self.run_pca)  # pca
        self.vcftools_plink_tools[1].on('end', self.run_structure)  # structure
        for tool in self.vcftools_plink_tools:
            gevent.sleep(1)
            tool.run()

    def run_tree(self):
        """
        进化树计算--默认计算nj树
        :return:
        """
        self.tree.set_options({
            "recode_vcf_path": self.vcftools_filter.option("filter_vcf").prop['path'],
            "tree_type": 'nj',
            "bs_trees": 500
        })
        self.tree.on('end', self.set_output, "tree")
        self.tree.on("start", self.set_step, {'start': self.step.tree})
        self.tree.on("end", self.set_step, {'end': self.step.tree})
        self.tree.run()

    def run_pca(self):
        """
        计算pca分析
        :return:
        """
        self.pca.set_options({
            "bfile_dir": self.vcftools_plink_tools[0].output_dir
        })
        self.pca.on("end", self.set_output, "pca")
        self.pca.on("start", self.set_step, {'start': self.step.pca})
        self.pca.on("end", self.set_step, {'end': self.step.pca})
        self.pca.run()

    def run_structure(self):
        """
        群体结构中--structure分析
        :return:
        """
        self.structure.set_options({
            "pop_fam": self.vcftools_plink_tools[1].output_dir + "/pop.fam",
            "pop_bed": self.vcftools_plink_tools[1].output_dir + "/pop.bed",
            "k_min": 2,
            "k_max": 20
        })
        self.structure.on("end", self.set_output, "structure")
        self.structure.on("start", self.set_step, {'start': self.step.structure})
        self.structure.on("end", self.set_step, {'end': self.step.structure})
        self.structure.run()

    def run_ld_analysis(self):
        """
        进行连锁不平衡分析
        :return:
        """
        self.ld_analysis.set_options({
            "vcf_file": self.vcf_file,
            "min_dp": "1",
            "max_dp": "300000000",
            "max_missing": "0.3",
            "min_maf": "0.05",
            "max_maf": "1",
            "group_dict": self.evolution_base.set_group_dict(self.group_file)
        })
        self.ld_analysis.on("end", self.set_output, "ld_analysis")
        self.ld_analysis.on("start", self.set_step, {'start': self.step.ld_analysis})
        self.ld_analysis.on("end", self.set_step, {'end': self.step.ld_analysis})
        self.ld_analysis.run()

    def run_sweep_analysis(self):
        """
        选择性消除分析
        :return:
        """
        self.sweep_analysis.set_options({
            "vcf_file": self.vcftools_filter.option("filter_vcf").prop['path'],
            "group_file": self.group_file,
            "diff_group": self.evolution_base.make_diff_group(self.group_file),
            "maxDP": 300000000,
            "min_maf": 0.05,
            'max_maf': 1,
            "analysis_method": 'all',
            "window_size": 2000000,
            "window_step": 10000,
            "vcftools_filter": False,
            "max_missing": 0.3,
            'minDP': 1
        })
        self.sweep_analysis.on("end", self.set_output, "sweep_analysis")
        self.sweep_analysis.on("start", self.set_step, {'start': self.step.sweep_analysis})
        self.sweep_analysis.on("end", self.set_step, {'end': self.step.sweep_analysis})
        # self.sweep_analysis.on("end", self.run_sweep_region)
        self.sweep_analysis.run()

    def run_sweep_region(self):
        """
        受选择区域筛选
        :return:
        """
        self.sweep_region.set_options({
            "vcf_file": self.vcftools_filter.option("filter_vcf").prop['path'],
            "pop_summary": self.pop_summary,
            "sweep_dir": self.sweep_analysis.output_dir,
            "diff_group": self.evolution_base.make_diff_group(self.group_file),
            "chr_set": self.chr_set
        })
        self.sweep_region.on("end", self.set_output, "sweep_region")
        self.sweep_region.on("start", self.set_step, {'start': self.step.sweep_region})
        self.sweep_region.on("end", self.set_step, {'end': self.step.sweep_region})
        self.sweep_region.run()

    def run_pop_history(self):
        """
        有效群体大小分析
        :return:
        """
        self.pop_history.set_options({
            "vcf_file": self.vcftools_filter.option("filter_vcf").prop['path'],
            "group_file": self.group_file,
            "vcftools_filter": False
        })
        self.pop_history.on("end", self.set_output, "pop_history")
        self.pop_history.on("start", self.set_step, {'start': self.step.pop_history})
        self.pop_history.on("end", self.set_step, {'end': self.step.pop_history})
        self.pop_history.run()

    def run_gwas_analysis(self):
        """
        全基因组关联分析
        :return:
        """
        self.gwas_analysis.set_options({
            "upload_trait_path": self.option("trait_file").prop['path'],
            "vcf_path": self.vcftools_filter.option("filter_vcf").prop['path']
        })
        self.gwas_analysis.on("end", self.set_output, "gwas_analysis")
        self.gwas_analysis.on("start", self.set_step, {'start': self.step.gwas_analysis})
        self.gwas_analysis.on("end", self.set_step, {'end': self.step.gwas_analysis})
        self.gwas_analysis.run()

    def run_hapmap_analysis(self):
        """
        单倍体型图
        :return:
        """
        self.hapmap_analysis.set_options({
            "anno_summary_path": self.anno_summary,
            "recode_vcf_path": self.gwas_analysis.option("recode_vcf_path"),
            "distance": 200000,
            "ld_r2": 0.8,
            "p_value": 0.05,
            "q_value": 0.05
        })
        self.hapmap_analysis.run()

    def run(self):
        tasks = [self.go_summary, self.var_diff]
        self.get_target_dir()  # 初始化一下 获取到远程磁盘的路径，用于保存对应文件的路径到主表中
        self.set_vcf_group_file()
        if self.option("data_type") == "raw_path" and self.option("genome_info").is_set:
            self.has_genome, self.genome_version_id, self.has_chromosomelist, self.species_path, self.species, self.\
                genome_version, self.secret = \
                self.evolution_base.set_genome_info(self.option("genome_info").prop['path'] + "/info.log",
                                                    self.is_wgs_result)
        self.logger.info("has_genome:{}".format(self.has_genome))
        self.evolution_base.add_sg_task(self._sheet.member_id, self._sheet.member_type, self._sheet.cmd_id)
        if self.is_gwas == 'true':
            self.vcftools_filter.on("end", self.run_gwas_analysis)
            self.step.add_steps("gwas_analysis")
            tasks.append(self.gwas_analysis)
        if self.option("group").is_set:
            self.logger.info("有分组文件，直接用分组文件进行后面计算")
            if self.is_structure == "true":
                tasks.append(self.structure)
                tasks.append(self.tree)
                tasks.append(self.pca)
                self.step.add_steps("tree")
                self.step.add_steps("pca")
                self.step.add_steps("structure")
            if self.is_ld_analysis == 'true':
                self.vcftools_filter.on("end", self.run_ld_analysis)
                self.step.add_steps("ld_analysis")
                tasks.append(self.ld_analysis)
            if self.is_pop_history == 'true':
                self.vcftools_filter.on("end", self.run_pop_history)
                self.step.add_steps("pop_history")
                tasks.append(self.pop_history)
            if self.is_sweep == 'true':
                self.vcftools_filter.on("end", self.run_sweep_analysis)
                self.step.add_steps("sweep_analysis")
                self.step.add_steps("sweep_region")
                tasks.append(self.sweep_analysis)
        else:
            self.logger.info("没有分组文件，将用structure的分组文件进行计算！")
            if self.is_structure == "true":
                tasks.append(self.structure)
                tasks.append(self.tree)
                tasks.append(self.pca)
                self.step.add_steps("tree")
                self.step.add_steps("pca")
                self.step.add_steps("structure")
                if self.is_ld_analysis == 'true':
                    self.structure.on('end', self.run_ld_analysis)
                    self.step.add_steps("ld_analysis")
                    tasks.append(self.ld_analysis)
                if self.is_pop_history == 'true':
                    self.structure.on('end', self.run_pop_history)
                    self.step.add_steps("pop_history")
                    tasks.append(self.pop_history)
                if self.is_sweep == 'true':
                    self.structure.on('end', self.run_sweep_analysis)
                    self.step.add_steps("sweep_analysis")
                    self.step.add_steps("sweep_region")
                    tasks.append(self.sweep_analysis)
            else:
                if self.is_ld_analysis == 'true':
                    self.set_error("没有分组的条件下，如果不做群体结构，不能做连锁不平衡分析！")
                if self.is_pop_history == 'true':
                    self.set_error("没有分组的条件下，如果不做群体结构，不能做种群历史分析！")
                if self.is_sweep == 'true':
                    self.set_error("没有分组的条件下，如果不做群体结构，不能做选择性消除分析！")
        self.on_rely(tasks, self.end)
        if self.option("data_type") != "raw_path":  # 输入的数据是wgs的路径
            self.logger.info("输入的是wgs的路径--后面只将进行遗传进化分析！")
            self.get_wgs_taskid()  # 检查 输入的结果是不是wgs的项目
            self.run_vcftools_filter()
            self.run_variant_compare()
            self.run_go_summary()
        else:
            self.logger.info("输入的是原始fastq的路径--后面先进行变异检测，然后进行遗传进化分析！")
            self.step.add_steps("qc_stat", "mapping", "call_snp_indel", "vcf_filter", "annovar")
            self.set_software()
            self.qc_stat.on('end', self.run_mapping)
            self.mapping.on("end", self.run_mapping_stat)
            self.mapping.on("end", self.run_call_snp_indel)
            # self.mapping_stat.on("end", self.run_send_mail)
            self.call_snp_indel.on("end", self.run_vcf_filter)
            self.vcf_filter.on("end", self.run_annovar)
            self.annovar.on("end", self.run_go_summary)
            self.annovar.on('end', self.run_vcftools_filter)
            self.annovar.on('end', self.run_variant_compare)
            tasks.extend([self.mapping_stat, self.annovar])
            if not self.has_genome:
                self.step.add_steps('genome_config')
                # self.genome_config.on('end', self.run_qc_stat)
                self.run_genome_config()
            else:
                self.set_ref_info()
                self.run_qc_stat()
        # gevent.spawn_later(5, self.end)
        super(EvolutionWorkflow, self).run()

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
            raise Exception('需要移动到output目录的文件夹不存在。')
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
        if event['data'] == 'genome_config':
            has_genome, genome_version_id = self.evolution_base.check_genome_ready(self.species, self.genome_version,
                                                                                   self.is_wgs_result)
            if not has_genome:
                self.move2outputdir(obj.output_dir, self.json_path + self.species_path)
                self.logger.info('指定基因组文件不存在，参考基因组移动到指定位置成功！')
            self.logger.info("开始执行set_ref_info")
            self.set_ref_info(is_update=True)
            self.logger.info("执行set_ref_info成功")
            self.logger.info("开始执行qc_stat模块")
            self.run_qc_stat()
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
        if event['data'] in ["tree", 'pca', 'structure']:
            self.move2outputdir(obj.output_dir, self.output_dir + "/06.pop")
        if event['data'] == "ld_analysis":
            self.move2outputdir(obj.output_dir, self.output_dir + "/07.ld")
        if event['data'] == "sweep_analysis":
            self.move2outputdir(obj.output_dir, self.output_dir + "/08.sweep/sweep_analysis")
        if event['data'] == "sweep_region":
            self.move2outputdir(obj.output_dir, self.output_dir + "/08.sweep/sweep_region")
        if event['data'] == "gwas_analysis":
            self.move2outputdir(obj.output_dir, self.output_dir + "/09.gwas")
        if event['data'] == "pop_history":
            self.move2outputdir(obj.output_dir, self.output_dir + "/10.history")
        if event['data'] == 'variant_compare':
            self.move2outputdir(obj.output_dir, self.output_dir + "/11.variant_compare")
        if event['data'] == 'vcftools_filter':
            self.move2outputdir(obj.output_dir, self.output_dir + "/12.vcf_filter")

    def send_files(self):
        self.remove_file()
        repaths = [
            [".", "", "基础分析结果文件夹"],
            ["01.fastq_qc", "", "数据质控结果目录"],
            ["01.fastq_qc/rawdata_qc", "", "原始数据统计结果目录"],
            ["01.fastq_qc/rawdata_qc/qc.xls", "xls", "原始数据统计表"],
            ["01.fastq_qc/rawdata_qc/atgc", "", "原始数据碱基含量分布"],
            ["01.fastq_qc/rawdata_qc/qual", "", "原始数据碱基错误率分布"],
            ["01.fastq_qc/clean_data/fastq.list", "", "各样品clean data路径"],
            ["01.fastq_qc/cleandata_qc", "", "质控后数据统计结果目录"],
            ["01.fastq_qc/cleandata_qc/qc.xls", "xls", "质控后数据统计表"],
            ["01.fastq_qc/cleandata_qc/atgc", "", "质控后数据碱基含量分布"],
            ["01.fastq_qc/cleandata_qc/qual", "", "质控后数据碱基错误率分布"],
            ["02.reference", "", "参考基因组目录"],
            ["02.reference/ref.fa", "", "参考基因组fa序列"],
            ["02.reference/ref.chrlist", "", "参考基因组fa文件"],
            ["02.reference/ref.changelog", "", "参考基因组染色体转换对应表"],
            ["02.reference/info.log", "", "参考基因组版本信息"],
            ["02.reference/ref.gff", "", "参考基因组gff文件"],
            ["02.reference/ssr.ref.result.xls", "xls", "参考基因组ssr文件"],
            ["03.map_stat", "", "基因组比对结果目录"],
            ["03.map_stat/result.stat/Total.mapped.detail.xls", "", "比对结果统计表"],
            ["03.map_stat/insert", "", "基因组比对插入片段长度结果目录"],
            ["03.map_stat/depth", "", "基因组比对测序深度结果目录"],
            ["03.map_stat/coverage", "", "基因组比对覆盖度结果目录"],
            ["03.map_stat/result.stat", "", "比对结果统计文件夹"],
            ["04.snp_indel", "", "snp与indel检测结果目录"],
            ["04.snp_indel/eff", "", "SNP/InDel的eff结果目录"],
            ["04.snp_indel/eff/snp.anno.primary.vcf", "", "snp功能注释vcf文件"],
            ["04.snp_indel/eff/indel.anno.primary.vcf", "", "indel功能注释vcf文件"],
            ["04.snp_indel/variant_stat", "", "SNP/InDel统计"],
            ["04.snp_indel/anno_stat", "", "SNP/InDel功能信息统计文件夹"],
            ["04.snp_indel/variant_stat/snp.stat", "", "SNP数据统计表"],
            ["04.snp_indel/variant_stat/snp.GQ", "", "SNP质量分布统计表"],
            ["04.snp_indel/variant_stat/snp.depth", "", "SNP深度分布统计表"],
            ["04.snp_indel/variant_stat/snp.matrix", "", "样本间差异SNP统计表"],
            ["04.snp_indel/variant_stat/indel.stat", "", "InDel数据统计表"],
            ["04.snp_indel/variant_stat/indel.len", "", "InDel长度信息统计表"],
            ["04.snp_indel/variant_stat/indel.matrix", "", "样本间差异InDel统计表"],
            ["04.snp_indel/variant_stat/indel.GQ", "", "InDel质量分布统计表"],
            ["04.snp_indel/variant_stat/indel.depth", "", "InDel深度分布统计表"],
            ["04.snp_indel/anno_stat/snp.stat", "", "SNP功能信息统计表"],
            ["04.snp_indel/anno_stat/indel.stat", "", "InDel功能信息统计表"],
            ["05.annovar", "", "基因功能注释结果目录"],
            ["05.annovar/combine_variants", "", "vcf合并结果目录"],
            ["05.annovar/combine_variants/pop.final.vcf", "", "vcf合并文件"],
            ["05.annovar/anno_count", "", "基因功能注释结果目录"],
            ["05.annovar/anno_count/pop.summary", "", "基因功能注释表"],
            ["05.annovar/anno_count/pop.stat.csv", "", "基因功能统计表"],
            ["05.annovar/anno_count/pop.go.stat", "", "GO功能统计表"],
            ["05.annovar/anno_count/pop.kegg.stat", "", "KEGG功能统计表"],
            ["05.annovar/anno_count/pop.eggnog.stat", "", "EGGNOG功能统计表"],
            ["05.annovar/eggnog_anno", "", "eggnog注释结果"],
            ["05.annovar/eggnog_anno/pop.eggnog.final.stat.detail", "", "eggnog注释表"],
            ["05.annovar/go_anno", "", "go注释结果"],
            ["05.annovar/go_anno/pop.go.final.stat.detail", "", "go注释表"],
            ["05.annovar/kegg_anno", "", "kegg注释结果"],
            ["05.annovar/kegg_anno/pop.kegg.final.stat.detail", "", "kegg注释表"],
            ["05.annovar/kegg_anno/pathway_dir", "", "kegg注释表"],
            ['06.pop', '', '群体结构结果目录'],
            ['06.pop/cverror', '', 'cverror分析结果目录'],
            ['06.pop/cverror/best.5.xls', '', ''],
            ['06.pop/cverror/cv.error', '', ''],
            ['06.pop/cverror/group.list', '', ''],
            ['06.pop/vcf2tree', '', 'vcf2tree目录'],
            ['06.pop/tree', '', '进化树分析结果目录'],
            ['06.pop/structure', '', 'structure结果目录'],
            ['07.ld', '', '连锁不平衡结果目录'],
            ['07.ld/vcf_filter_dir', '', 'LD vcf结果目录'],
            ['07.ld/ld_decay_dir', '', 'LD decay结果目录'],
            ['07.ld/ld_draw_dir', '', 'LD draw结果目录'],
            ['08.sweep', '', '选择性消除结果目录'],
            ['08.sweep/sweep_analysis', '', '选择性消除结果目录'],
            ['09.gwas', '', '全基因关联分析结果目录'],
            ['10.history', '', '种群历史结果目录'],
            ['11.variant_compare', '', '变异位点比较分析结果目录'],
            ['12.vcf_filter', '', 'vcftools_filter结果目录']
        ]

        regexps = [
            [r"01.fastq_qc/rawdata_qc/atgc/.*\.xls", "xls", "各样本原始数据碱基含量分布"],
            [r"01.fastq_qc/rawdata_qc/qual/.*\.xls", "xls", "各样本原始数据碱基错误率分布"],
            [r"01.fastq_qc/cleandata_qc/atgc/.*\.xls", "xls", "各样本质控后数据碱基含量分布"],
            [r"01.fastq_qc/cleandata_qc/qual/.*\.xls", "xls", "各样本质控后数据碱基错误率分布"],
            [r"03.map_stat/insert/.*\.xls", "xls", "各样本基因组比对插入片段长度"],
            [r"03.map_stat/depth/.*\.xls", "xls", "各样本基因组比对测序深度"],
            [r"03.map_stat/coverage/.*\.xls", "xls", "各样本基因组比对覆盖度"],
            [r"03.map_stat/insert/.*\.xls", "xls", "各样本基因组插入片段长度"],
            [r"05.annovar/kegg_anno/pathway_dir/.*\.pdf", "", "KEGG pathway通路图"],
            [r"05.annovar/kegg_anno/pathway_dir/.*\.png", "", "KEGG pathway通路图"]
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
        self.import_specimen_info()
        self.import_mapping_results()
        self.import_snp_results()
        self.import_indel_results()
        if self.is_annovar == 'true':
            self.import_annovar_results()
        if self.is_structure == 'true':
            self.import_pop_analysis()
        if self.is_ld_analysis == 'true':
            self.import_ld_analysis()
        if self.is_sweep == 'true':
            self.import_sweep_analysis()
        if self.is_pop_history == 'true':
            self.import_history()
        self.import_variant_compare()
        if self.is_gwas == 'true':
            self.import_gwas()
        self.import_files_paths()
        self.send_files()

    def set_vcf_group_file(self):
        """
        设置初始化一下vcf文件与group分组方案文件
        :return:
        """
        if self.option("group").is_set:
            self.group_file = self.remove_header(self.option("group").prop['path'])
        else:
            self.group_file = self.structure.output_dir + "/cverror/group.list"
        if self.option("data_type") == "raw_path":
            self.vcf_file = self.annovar.output_dir + "/combine_variants/pop.final.vcf"
        else:
            self.vcf_file = self.option("wgs_path").prop['path'] + "/05.annovar/combine_variants/pop.final.vcf"

    def remove_header(self, group_table):
        """
        因为后面计算的时候，不要header头文件，所以这里就将器删除
        :param group_table:
        :return:
        """
        with open(group_table, 'r') as r, open(self.output_dir + "/group.txt", 'w') as w:
            datas = r.readlines()
            for m in datas:
                if re.match('#.*', m):
                    pass
                else:
                    w.write(m)
        return self.output_dir + "/group.txt"

    def get_wgs_taskid(self):
        """
        根据输入的文件夹路径获取wgs的路径，files/m_188/188_5af3dc68d1df4/tsanger_30196/workflow_results
        如果能够获取到wgs task id则直接从wgs项目中拷贝相关的样本信息，否则就要给出样本信息列表，具体参考fastq.txt,
        :return:
        """
        self.logger.info(self.wgs_path)
        task_id = os.path.basename(os.path.dirname(self.wgs_path))
        self.logger.info("WGS的task_id为：{}".format(task_id))
        result = self.wgs_api.col_find_one("sg_task", {"task_id": task_id})
        if result:
            self.has_wgs_id = True
            self.wgs_task_id = task_id
            self.genome_version_id = result['genome_version_id']
            params = self.wgs_api.col_find_one("sg_snp_call", {"task_id": task_id})
            if params:
                self.software = params['params']
            else:
                self.software = '{\"\":\"\"}'
            self.logger.info("设置软件名称成功！")
        else:
            self.set_error("{}在wgs项目的sg_task中没有找到对应信息！".format(task_id))

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
        if not self.has_wgs_id:
            wgs_base = self.api.api('wgs.wgs_base')
            wgs_base.project_type = "dna_evolution"
            wgs_base.add_sg_specimen()
        else:
            sample_info = self.api.api('wgs.wgs_base').get_specimen_other_info(self.wgs_task_id)
            gmap_base = self.api.api("dna_gmap.gmap_base")
            gmap_base.project_type = "dna_evolution"
            gmap_base.add_sg_specimen(sample_info)
        self.api.del_api("wgs.wgs_base")
        wgs_base = self.api.api('wgs.wgs_base')
        wgs_base.project_type = "dna_evolution"
        self.logger.info("样本信息导表成功！")
        # wgs_base.update_clean_path(self.api_base_path + "/01.fastq_qc/clean_data")
        # self.logger.info("更新clean fastq成功！")
        self.logger.info("开始进行样本分组导表")
        self.group_id, self.group_dict = self.api.api("dna_evolution.evolution_base").\
            add_sg_specimen_group(self.group_file, True if self.option("group").is_set else False)
        self.logger.info("样本分组导表成功！")
        self.logger.info("开始进行样本质控信息导表")
        wgs_base.add_sg_specimen_qc(self.api_base_path + "/01.fastq_qc/rawdata_qc/qc.xls",
                                    self.api_base_path + "/01.fastq_qc/cleandata_qc/qc.xls")
        self.logger.info("样本质控信息导表成功！")
        self.logger.info("开始进行碱基含量分布导表")
        wgs_base.add_qc_atgc_curve(self.api_base_path + "/01.fastq_qc/rawdata_qc/atgc", "raw_reads")
        wgs_base.add_qc_atgc_curve(self.api_base_path + "/01.fastq_qc/cleandata_qc/atgc", "clean_reads")
        self.logger.info("碱基含量分布导表成功！")
        self.logger.info("开始进行碱基错误率分布导表")
        wgs_base.add_qc_qual_curve(self.api_base_path + "/01.fastq_qc/rawdata_qc/qual", "error_raw_reads")
        wgs_base.add_qc_qual_curve(self.api_base_path + "/01.fastq_qc/cleandata_qc/qual", "error_clean_reads")
        self.logger.info("碱基错误率分布导表成功！")

    def import_mapping_results(self):
        wgs_base = self.api.api('wgs.wgs_base')
        wgs_base.project_type = "dna_evolution"
        self.logger.info("开始进行比对结果导表")
        mapping_id = wgs_base.add_sg_mapping()
        wgs_base.add_sg_mapping_detail(mapping_id,
                                       self.api_base_path + "/03.map_stat/result.stat/Total.mapped.detail.xls")
        self.logger.info("比对结果导表成功！")
        self.logger.info("开始进行插入片段分布导表")
        wgs_base.insert_sg_mapping_curve(self.api_base_path + "/03.map_stat/insert", mapping_id, "insert")
        self.logger.info("插入片段分布导表成功！")
        self.logger.info("开始进行测序深度分布导表")
        wgs_base.insert_sg_mapping_curve(self.api_base_path + "/03.map_stat/depth", mapping_id, "dep")
        self.logger.info("测序深度分布导表成功！")
        self.logger.info("开始进行基因组覆盖度分布导表")
        coverage_windows = self.api.api('dna_evolution.coverage_window')
        windows_id = coverage_windows.add_sg_coverage_window({"step_num": 100, "submit_location": "coverage_window",
                                                              "task_type": 2, "file_id": str(mapping_id)}, mapping_id)
        coverage_windows.get_sample(self.api_base_path + "/03.map_stat/coverage/", windows_id)
        self.api_base.update_db_record("sg_mapping", {'_id': mapping_id},
                                       {"dep_path": self.base_target_path + "/03.map_stat/map_bam/sort_bams/"})
        self.logger.info("基因组覆盖度分布导表成功！")

    def import_snp_results(self):
        self.api.del_api("dna_evolution.evolution_base")
        snp_api = self.api.api('dna_evolution.evolution_base')
        self.logger.info("开始进行snp统计导表")
        self.call_id = snp_api.add_variant_call(self._sheet.project_sn, self._sheet.id, params=self.software)
        snp_api.add_sg_snp_call_stat(self.call_id, self.api_base_path + "/04.snp_indel/variant_stat/snp.stat")
        self.logger.info("snp统计导表成功！")
        self.logger.info("开始进行snp质量评估导表")
        snp_api.add_snp_qc_curve(self._sheet.id, self.call_id,
                                 self.api_base_path + "/04.snp_indel/variant_stat/snp.GQ", "snp_qc", "snp_qc")
        snp_api.add_snp_qc_curve(self._sheet.id, self.call_id,
                                 self.api_base_path + "/04.snp_indel/variant_stat/snp.depth", "snp_depth", "snp_depth")
        self.logger.info("snp质量评估导表成功！")
        self.logger.info("开始进行snp功能注释导表")
        self.anno_id = snp_api.add_sg_variant_anno(self._sheet.project_sn, self._sheet.id,
                                                   params='{\"method\":\"SNPEff\"}')
        # snp功能统计表
        snp_api.add_sg_snp_anno_stat(self.anno_id, self.api_base_path + "/04.snp_indel/anno_stat/snp.stat",
                                     "annotation")
        # snp功效统计表
        snp_api.add_sg_snp_anno_stat(self.anno_id, self.api_base_path + "/04.snp_indel/anno_stat/snp.stat", "effect")
        # snp功效与功能累加图与直方图
        snp_api.add_sg_snp_anno_bar(self._sheet.id, self.anno_id,
                                    self.api_base_path + "/04.snp_indel/anno_stat/snp.stat")
        self.logger.info("snp功能注释导表成功！")

    def import_indel_results(self):
        indel_api = self.api.api('dna_evolution.evolution_base')
        self.logger.info("开始进行indel统计导表")
        indel_api.add_sg_indel_call_stat(self.call_id, self.api_base_path + "/04.snp_indel/variant_stat/indel.stat")
        self.logger.info("indel统计表导表成功")
        self.logger.info("开始进行indel长度分布导表")
        indel_api.add_indel_length_bar(self._sheet.id, self.call_id,
                                       self.api_base_path + "/04.snp_indel/variant_stat/indel.len")
        self.logger.info("indel长度分布导表成功！")
        self.logger.info("开始进行indel质量评估导表")
        indel_api.add_indel_qc_curve(self._sheet.id, self.call_id,
                                     self.api_base_path + "/04.snp_indel/variant_stat/indel.GQ", "indel_qc")
        indel_api.add_indel_qc_curve(self._sheet.id, self.call_id,
                                     self.api_base_path + "/04.snp_indel/variant_stat/indel.depth", "indel_depth")
        self.logger.info("indel质量评估导表成功！")
        self.logger.info("开始进行indel功能注释导表")
        # indel功能统计表
        indel_api.add_sg_indel_anno_stat(self.anno_id, self.api_base_path + "/04.snp_indel/anno_stat/indel.stat",
                                         "annotation")
        # indel功效统计表
        indel_api.add_sg_indel_anno_stat(self.anno_id, self.api_base_path + "/04.snp_indel/anno_stat/indel.stat",
                                         "effect")
        # indel功效与功能累加图与直方图
        indel_api.add_sg_indel_anno_bar(self._sheet.project_sn, self._sheet.id, self.anno_id,
                                        self.api_base_path + "/04.snp_indel/anno_stat/indel.stat")
        self.logger.info("indel功能注释导表成功！")

    def import_annovar_results(self):
        """
        要修改
        :return:
        """
        anno_api = self.api.api("dna_evolution.region_anno")
        self.logger.info("开始进行基因功能注释导表")
        anno_params_id = anno_api.add_sg_anno_params(self._sheet.id, self._sheet.project_sn, 1,
                                                     {"source": "Whole Genome Annotation"},
                                                     "Whole Genome Annotation", "Whole Genome Annotation",
                                                     ['All Region'],
                                                     self.base_target_path + '/05.annovar/anno_count/pop.summary')
        anno_id = anno_api.add_sg_region_anno(self._sheet.project_sn, self._sheet.id,
                                              {'sg_anno_params_id': str(anno_params_id),
                                               "region_select": ['All Region'], "submit_location": "region_anno"})
        anno_api.add_sg_region_anno_detail(anno_id, self.api_base_path + "/05.annovar/anno_count/pop.summary",
                                           self.genome_version_id)
        self.logger.info("基因功能注释导表成功！")
        self.logger.info("开始进行GO/KEGG/EGGNOG统计图表导表")
        anno_api.sg_region_anno_go_stat(self._sheet.id, anno_id,
                                        self.api_base_path + "/05.annovar/go_anno/pop.go.final.stat.detail",
                                        self.go_summary.output_dir + "/pop.2.enrich",
                                        self.api_base_path + "/05.annovar/anno_count/pop.summary")
        anno_api.sg_region_anno_kegg_stat(self._sheet.id, anno_id,
                                          self.api_base_path + "/05.annovar/kegg_anno/pop.kegg.final.stat.detail",
                                          self.api_base_path + "/05.annovar/kegg_anno/pathway_dir",
                                          self.base_target_path + "/05.annovar/kegg_anno/pathway_dir/")
        anno_api.sg_region_anno_eggnog_stat(self._sheet.id, anno_id,
                                            self.api_base_path + "/05.annovar/eggnog_anno/pop.eggnog.final.stat.detail",
                                            self.api_base_path + "/05.annovar/anno_count/pop.summary")
        self.logger.info("GO/KEGG/EGGNOG统计图表导表成功！")
        self.vcf_id = self.evolution_base.add_sg_variant_compare_filter(
            self.base_target_path + "/05.annovar/combine_variants/pop.final.vcf", self._sheet.id,
            self._sheet.project_sn)
        self.evolution_base.sg_software(self._sheet.project_sn, self._sheet.id, self.call_type)

    def import_pop_analysis(self):
        """
        群体结构导表
        :return:
        """
        self.logger.info("开始进行PCA导表")
        pop_api = self.api.api('dna_evolution.pop_analysis')
        pca_params = {"group_dict": self.group_dict, "group_id": str(self.group_id), "max_maf": '1',
                      "max_missing": '0.3', "maxdp": '',  "min_maf": '0.05', 'mindp': '10', 'vcf_id': str(self.vcf_id),
                      'submit_location': 'pop_pca', 'task_type': 2}
        pca_id = pop_api.add_sg_pop_pca(self._sheet.id, self._sheet.project_sn, pca_params)
        pop_api.add_sg_pop_pca_detail(pca_id, self.output_dir + "/06.pop/pop.pca.eigenvec")
        pop_api.add_sg_scatter(self._sheet.id, pca_id, self.output_dir + "/06.pop/pop.pca.eigenvec")
        pop_api.add_pca_bar(self._sheet.id, pca_id, self.output_dir + '/06.pop/pop.pca.eigenval')
        self.logger.info("PCA导表完成--将进行structure导表")
        structure_params = {"group_dict": self.group_dict, "group_id": str(self.group_id), "max_maf": '1',
                            "max_missing": '0.3', "maxdp": '', "min_maf": '0.05', 'mindp': '10',
                            'vcf_id': str(self.vcf_id), 'kmax': '20', 'kmin': '2',
                            'submit_location': 'pop_structure', 'task_type': 2}
        structure_id = pop_api.add_sg_pop_structure(self._sheet.id, self._sheet.project_sn, structure_params)
        pop_api.add_sg_kvalue_detail(structure_id, self.output_dir + "/06.pop/structure/", self._sheet.id)
        pop_api.add_structure_sg_curve(self._sheet.id, structure_id, self.output_dir + "/06.pop/cverror/cv.error")
        self.logger.info("Structure导表完成--将进行tree导表")
        tree_params = {"group_dict": self.group_dict, "group_id": str(self.group_id), "max_maf": '1',
                       "max_missing": '0.3', "maxdp": '', "min_maf": '0.05', 'mindp': '10', 'vcf_id': str(self.vcf_id),
                       'bootstrap': '1000', 'tree_type': 'nj', "submit_location": 'pop_tree', 'task_type': 2}
        tree_id = pop_api.add_sg_pop_tree(self._sheet.id, self._sheet.project_sn, tree_params)
        pop_api.add_pop_sg_tree(self._sheet.id, tree_id, self.output_dir + "/06.pop/tree/pop.nj.tree")
        self.logger.info("Tree导表完成")

    def import_ld_analysis(self):
        """
        连锁不平衡导表
        :return:
        """
        self.logger.info("开始进行LD导表")
        ld_api = self.api.api('dna_evolution.ld_analysis')
        ld_params = {"group_dict": self.group_dict, "group_id": str(self.group_id), "max_maf": '1',
                     "max_missing": '0.3', "maxdp": '', "min_maf": '0.05', 'mindp': '10', 'vcf_id': str(self.vcf_id),
                     'task_type': 2}
        ld_id = ld_api.sg_ld(self._sheet.id, self._sheet.project_sn, ld_params)
        ld_api.sg_ld_detail(ld_id, self.ld_analysis.work_dir + "/gro_list")
        ld_api.update_ld_analysis(self.target_dir + "/07.ld/ld_draw_dir/ld.png", ld_id)
        # ld_api.add_ld_curve(self.ld_analysis.work_dir + "/gro_list", ld_id, self._sheet.id)
        self.logger.info("LD导表完成")

    def import_sweep_analysis(self):
        """
        选择性消除导表
        :return:
        """
        self.logger.info("开始进行选择性消除导表")
        swan_api = self.api.api('dna_evolution.sweep_analysis')
        diff_group_name = self.evolution_base.make_diff_group(self.group_file)
        params_json = {
            "vcf_id": str(self.vcf_id),
            "group_dict": self.group_dict,
            "diff_group": [diff_group_name],
            "group_id": str(self.group_id),
            "analysis_method": 'all',
            "window_size": 2000,
            "max_missing": '0.3',
            "min_maf": '0.05',
            "max_maf": '1',
            "mindp": '1',
            "maxdp": '',
            "submit_location": 'sweep',
            'task_type': 2
        }
        sweep_id = swan_api.add_sg_sweep([diff_group_name], params_json)
        window_step = 10000
        variant_num_path = os.path.join(self.sweep_analysis.work_dir, "vcf.variant.txt")
        pop1 = diff_group_name.split("_vs_")[0]
        pop2 = diff_group_name.split("_vs_")[1]
        diff_file = pop1 + "-" + pop2 + ".pi_tajimaD_fst.detail"
        pi_tajimad_fst_path = os.path.join("{}/08.sweep/sweep_analysis".format(self.output_dir), diff_file)
        if os.path.exists(pi_tajimad_fst_path):
            swan_api.add_sg_sweep_detail(sweep_id, pi_tajimad_fst_path, variant_num_path, diff_group_name, window_step)
        else:
            self.logger.info("Can not find {}.pi_tajimaD_fst.detail file!".format(diff_group_name))
        sweep_dir = "{}/08.sweep/sweep_analysis/".format(self.target_dir)
        png_path = os.path.join(sweep_dir, pop1 + "-" + pop2 + ".pop1.manhattan.png")
        pdf_path = os.path.join(sweep_dir, pop1 + "-" + pop2 + ".pop1.manhattan.pdf")
        swan_api.add_sg_manhattan_path(sweep_id, diff_group_name, pop1, png_path, pdf_path)
        png_path = os.path.join(sweep_dir, pop1 + "-" + pop2 + ".pop2.manhattan.png")
        pdf_path = os.path.join(sweep_dir, pop1 + "-" + pop2 + ".pop2.manhattan.pdf")
        swan_api.add_sg_manhattan_path(sweep_id, diff_group_name, pop2, png_path, pdf_path)
        vcf_path = os.path.join(self.target_dir + '/12.vcf_filter', "pop.recode.vcf")
        swan_api.update_sweep_path(sweep_id, sweep_dir, vcf_path)
        self.logger.info("选择性消除导表完成！")

    def import_sweep_region(self):
        """
        受选择区域筛选
        :return:
        """
        self.logger.info("开始进行选择性区域筛选导表")
        region_api = self.api.api("dna_evolution.sweep_region")
        sweep_region_id = region_api.add_sg_sweep_region()
        diff_group_name = self.evolution_base.make_diff_group(self.group_file)
        select_region_stat = ''
        region_api.add_sg_sweep_region_detail(sweep_region_id, diff_group_name, select_region_stat)
        region_api.update_diff_group(sweep_region_id, [diff_group_name])
        self.logger.info("选择性区域筛选导表完成！")
        pass

    def import_history(self):
        """
        种群历史--有效群体大小分析
        :return:
        """
        self.logger.info("将种群历史结果导入mongo表")
        history_api = self.api.api("dna_evolution.pop_history")
        params_json = {
            "vcf_id": str(self.vcf_id),
            "group_dict": self.group_dict,
            "group_id": str(self.group_id),
            "analysis_method": 'psmc',
            "max_missing": '0.3',
            "min_maf": '0.05',
            "max_maf": '1',
            "mindp": '1',
            "maxdp": '',
            "submit_location": 'pop_history',
            'task_type': 2
        }
        history_id = history_api.add_sg_history(self._sheet.project_sn, self._sheet.id, params_json)
        psmc_stat = os.path.join(self.output_dir + "/10.history/all_psmc_stat.xls")
        history_api.add_sg_history_detail(history_id, psmc_stat)
        curve_id = history_api.add_sg_curve(history_id)
        history_api.add_sg_curve_detail(curve_id, psmc_stat)
        self.logger.info("种群历史结果导表成功")
        pass

    def import_variant_compare(self):
        """
        突变比较分析结果导表
        :return:
        """
        self.logger.info("开始突变比较分析导表！")
        api_path = self.api.api("dna_evolution.variant_compare")
        chromosome_window = self.api.api("dna_evolution.chromosome_window")
        params_json = {"step_num": 10,
                       "variant_type": "all",
                       "project_sn": self._sheet.project_sn,
                       # "task_id": self._sheet.id,
                       "funtype": "high,moderate,low,modifier",
                       "efftype": "3_prime_UTR_variant,stop_losti,upstream_gene_variant",
                       "task_type": 2,
                       "dep": '0',
                       'location': '',
                       "group_dict": {self.group_dict.keys()[0]: self.group_dict[self.group_dict.keys()[0]]},
                       "group_id": str(self.group_id),
                       "type": "sample",
                       "analysi_type": "3"
                       }
        download_path = self.target_dir + "/11.variant_compare/pop.table"
        main_id = api_path.sg_variant_compare(self._sheet.id, self._sheet.project_sn, params_json)
        window_params = {"file_id": str(main_id), "step_num": 100, "submit_location": "variantcompare_window",
                         "task_type": 2, "variant_type": 'all'}
        chromosome_id = chromosome_window.add_sg_chromosome_window(self._sheet.project_sn, window_params,
                                                                   self._sheet.id, main_id)
        chromosome_window.add_sg_distribution(self.output_dir + "/11.variant_compare", chromosome_id, "all")
        api_path.update_varaint_compare(download_path, self.target_dir + "/11.variant_compare/pop.filtered.vcf",
                                        main_id)
        api_path.add_sg_variant_compare_effect(main_id, self.output_dir + "/11.variant_compare/eff.type")
        api_path.add_variant_compare_effect_bar(main_id, self._sheet.id,
                                                self.output_dir + "/11.variant_compare/eff.type", 'all')
        api_path.sg_variant_compare_impact(main_id,  self.output_dir + "/11.variant_compare/function.type")
        api_path.sg_variant_compare_impact_bar(main_id, self._sheet.id,
                                               self.output_dir + "/11.variant_compare/function.type", 'all')
        api_path.sg_varian_compare_detail(main_id, self.output_dir + "/11.variant_compare/pop.table", 'all')
        api_path.add_variant_compare_effect_bar(main_id, self._sheet.id,
                                                self.output_dir + "/11.variant_compare/eff.type", "snp")
        api_path.add_variant_compare_effect_bar(main_id, self._sheet.id,
                                                self.output_dir + "/11.variant_compare/eff.type", "indel")
        api_path.sg_variant_compare_impact_bar(main_id, self._sheet.id,
                                               self.output_dir + "/11.variant_compare/function.type", "snp")
        api_path.sg_variant_compare_impact_bar(main_id, self._sheet.id,
                                               self.output_dir + "/11.variant_compare/function.type", "indel")
        api_path.sg_varian_compare_detail(main_id, self.output_dir + "/11.variant_compare/pop.table", "snp")
        api_path.sg_varian_compare_detail(main_id, self.output_dir + "/11.variant_compare/pop.table", "indel")
        self.logger.info("突变比较分析导表成功！")

    def import_gwas(self):
        """
        全基因关联分析
        :return:
        """
        self.logger.info("开始进行GWAS结果导表")
        seq_api = self.api.api("dna_evolution.gwas_analysis")
        gwas_id = seq_api.add_sg_gwas_analysis(self._sheet.project_sn, self._sheet.id,
                                               self.option("trait_file").prop['path'])
        bar_dir = self.output_dir + "/trait_annovar"
        trait_stat = os.path.join(bar_dir, "annovar.trait.xls")
        seq_api.add_sg_gwas_analysis_stat(gwas_id=gwas_id, trait_stat=trait_stat, task_id=self._sheet.id,
                                          bar_dir=bar_dir)
        csv_dir = os.path.join(self.output_dir, "09.gwas/gwas_mvp")
        seq_api.add_sg_manhattan(data_path=csv_dir, main_id=gwas_id)
        seq_api.add_sg_scatter(task_id=self._sheet.id, origin_id=gwas_id, gwas_output_dir=csv_dir)
        recode_vcf_path = os.path.join(os.path.join(self.output_dir, "09.gwas/vcftools_filter"), "pop.recode.vcf")
        self.api.api('dna_evolution.api_base').update_db_record(collection="sg_gwas_analysis",
                                                                query_dict={"_id": gwas_id},
                                                                update_dict={"csv_dir": csv_dir,
                                                                             "recode_vcf_path": recode_vcf_path})
        self.logger.info("GWAS结果导表成功")
        pass

    def import_files_paths(self):
        file_paths = {
            # "coverage_path": self.base_target_path + "/03.map_stat/coverage/",
            "pop_final_vcf": self.base_target_path + "/05.annovar/combine_variants/pop.final.vcf",
            "genome_version_id": ObjectId(self.genome_version_id),
            "region": self.region,
            "pop_summary": self.base_target_path + "/05.annovar/anno_count/pop.summary",
            "chr_set": self.chr_set,
            "is_wgs_result": self.is_wgs_result
        }
        self.api_base.update_db_record("sg_task", {"task_id": self._sheet.id}, file_paths)

    def end(self):
        self.logger.info("上传文件目录完成！")
        self.run_api()
        self.logger.info("全部导表完成！")
        super(EvolutionWorkflow, self).end()

    def get_target_dir(self):
        """
        获取远程磁盘的路径
        :return:
        """
        self.target_dir = self._sheet.output.rstrip("/")
        self.region = self._sheet.output.strip().split('://')[0]
        if self.option("data_type") == "raw_path":
            self.pop_final_vcf = self._sheet.output.rstrip('/') + "/05.annovar/combine_variants/pop.final.vcf"
            self.api_base_path = self.output_dir
            self.base_target_path = self._sheet.output.rstrip('/')
        else:
            self.api_base_path = self.option("wgs_path").prop['path']
            if self._sheet.options()['wgs_path'].startswith("/mnt"):
                self.pop_final_vcf = os.path.join(self.option('wgs_path').prop['path'].rstrip('/') +
                                                  "/05.annovar/combine_variants/pop.final.vcf")
                self.base_target_path = self._sheet.options()['wgs_path'].rstrip('/')
            else:
                self.wgs_path = self.get_wgs_path()  # 如果是对象存储的话定位到workflow_results
                if not self.wgs_path.startswith("/mnt"):
                    self.pop_final_vcf = os.path.join(self.option('wgs_path').prop['path'].rstrip('/') +
                                                      "/05.annovar/combine_variants/pop.final.vcf")
                else:  # 这一步是为了兼容开发时工作流进行测试
                    self.pop_final_vcf = self.wgs_path.rstrip('/') + "/05.annovar/combine_variants/pop.final.vcf"
                self.base_target_path = self.wgs_path.rstrip('/')
        self.logger.info("base_path:{}".format(self.api_base_path))
        self.logger.info("pop_final_vcf:{}".format(self.pop_final_vcf))

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
                code = os.system("rm -rf {}".format(files))
                if code != 0:
                    self.logger.info("删除文件夹{}失败！".format(files))
            else:
                self.logger.info("文件{}不存在，不用删除！".format(files))

    def make_compare_config(self):
        """
        用于给变异位点比较生成配置文件
        Variant Type=SNP,INDEL
        Variant Eff=HIGH,MODERATE,LOW,MODIFIER
        Variant Ann=3_prime_UTR_variant,stop_losti,upstream_gene_variant
        Sample Diff=16S15,0,100,0,16S456_1,0,100,0,0
         # sample dep1 dep2 hete sample2 dep3 dep4 hete diff; while 1 is hete or diff
        Group Info=parent,5,100,0,0.1,0,1
        # group1 ad1 ad2 miss1 miss2 maf1 maf2
        :return:
        """
        sample1, sample2 = self.evolution_base.make_compare_sample(self.api_base_path + "/03.map_stat/coverage")
        with open(os.path.join(self.work_dir, "config.txt"), "w") as w, open(
                os.path.join(self.work_dir, "compare_group.txt"), 'w') as w1:
            w.write("Variant Type=SNP,INDEL\n")
            w.write("Variant Eff=HIGH,MODERATE,LOW,MODIFIER\n")
            w.write("Variant Ann=3_prime_UTR_variant,stop_losti,upstream_gene_variant\n")
            w.write("Sample Diff={},0,100,0,{},0,100,0,0".format(sample1, sample2))
            w1.write('{}\tdefault\n'.format(sample1))
            w1.write('{}\tdefault\n'.format(sample2))

    def run_send_mail(self):
        """
        mapping完之后，发送比较结果邮件给产品线
        http://www.sanger.com/task/project_tasks/project_id/64527.html
        :return:
        """
        a = SendEmail("1274095891@qq.com", "smtp.qq.com", "ocfnjnhnsbqubaej", "1274095891@qq.com",
                      "hongdong.xuan@majorbio.com", "项目编号:{}, 任务id:{} 的比对结果信息"
                      .format(self._sheet.project_sn, self._sheet.id), 465)
        a.send_msg("{}".format("http://www.{}.com/task/project_tasks/project_id/{}.html"
                               .format(self._sheet.id.split('_')[0], self.option('email_id'))))
        a.attachment(self.mapping_stat.output_dir + "/result.stat/Total.mapped.detail.xls")
        a.send_email()
        self.logger.info("邮件发送完成")
        pass

    def get_wgs_path(self):
        """
        从mapping_file.txt中获取wgs的路径，用于后面导表使用
        :return:
        """
        file_path = self.work_dir + "/remote_input/wgs_path/mapping_file.txt"
        if not os.path.exists(file_path):
            raise Exception("file:{} is not exists!")
        with open(file_path, 'r') as r:
            data = r.readlines()[0]
            json_data = json.loads(data)
        temp_path = json_data['wgs_path'][0]["file_path"].split("workflow_results")[0].rstrip('/') + '/workflow_results'
        return temp_path

    def set_call_type(self, ref_chrlist):
        """
        检查染色体是不是大于536870912，如果有一条大于536870912，直接用samtool进行call snp不能用gatk跑
        :return:
        """
        self.call_type = self.option("snp_indel_method").lower()
        with open(ref_chrlist, 'r') as r:
            data = r.readlines()
            for line in data:
                temp = line.strip().split('\t')
                if int(temp[1]) > 536870912:
                    self.call_type = 'samtools'
                    self.large_genome = "true"
                    break
        self.logger.info("call type is {}--large_genome is {}".format(self.call_type, self.large_genome))
