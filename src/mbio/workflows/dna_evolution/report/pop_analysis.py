# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180910

import re
import os
from bson.objectid import ObjectId
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class PopAnalysisWorkflow(Workflow):
    """
    遗传进化群体结构--包含了进化树（ml，nj，bayes），structrue，pca, 每次只运行一个分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PopAnalysisWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "vcf_file", "type": "infile", "format": "dna_gmap.vcf", "required": True},
            # pop.recode.vcf或者比较分析的结果
            # {"name": "group_file", "type": "infile", "format": "dna_evolution.group_table"},  # 分组文件可以没有
            {"name": "tree_type", "type": "string", "default": "nj"},  # tree 建树方式
            {"name": "bs_trees", "type": "int", "default": 1000},  # tree params：Bootstrap
            {"name": "recode", "type": "bool", "default": True},
            {"name": "remove_indels", "type": "bool", "default": False},  # --remove-indels
            {"name": "remove_filtered_all", "type": "bool", "default": False},  # --remove-filtered-all
            {"name": "minDP", "type": "int", "default": 1},  # 平均测序深度min
            {"name": "maxDP", "type": "int", "default": 3000000000000},  # 平均测序深度max
            {"name": "max_missing", "type": "float", "default": 0.3},  # 缺失率
            {"name": "min_maf", "type": "float", "default": 0.05},  # 次要等位基因频率min
            {"name": "max_maf", "type": "float", "default": 1},  # 次要等位基因频率max
            {"name": "kmax", "type": "int", "default": 20},  # K值最大值 structure
            {"name": "kmin", "type": "int", "default": 2},  # K值最小值 structure
            {"name": "analysis_type", "type": "string", "default": "structure"},  # 确定分析类型structure or tree or pca
            {"name": "chr_set", "type": "int"},
            {"name": "sample_list", "type": "string"},   # 传入样本列表，用于在vcf过滤的时候保留指定样本
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.vcftools_filter = self.add_tool("dna_evolution.vcftools_filter")
        self.vcftools_plink = None
        self.structure = None
        self.tree = None
        self.pca = None

    def check_options(self):
        if self.option("minDP") >= self.option("maxDP"):
            raise OptionError("测序深度minDP:{}必须小于maxDP:{}".format(self.option("minDP"), self.option("maxDP")))
        if self.option("min_maf") >= self.option("max_maf"):
            raise OptionError("次要等位基因频率min_maf:{}必须小于max_maf:{}"
                              .format(self.option("min_maf"), self.option("max_maf")))
        if self.option("max_maf") > 1:
            raise OptionError("次要等位基因频率最大值:{}必须小于等于1".format(self.option("max_maf")))
        if self.option("max_missing") < 0 or self.option("max_missing") > 1:
            raise OptionError("缺失率max_missing:{}必须在[0,1]".format(self.option("max_missing")))
        if self.option("tree_type").lower() not in ["nj", 'ml', 'bayes']:
            raise OptionError("树的计算方法必须是nj or ml or bayes", code="11111")
        if type(self.option("bs_trees")) != int:
            raise OptionError("bs_trees参数必须是整型", code="11111")
        else:
            if self.option("bs_trees") <= 1:
                raise OptionError("bs_trees参数必须是大于1的整数", code="11111")
        if self.option("kmin"):
            if type(self.option("kmin")) != int:
                raise OptionError("k_value一定要是整型数据", code="123345")
            else:
                if self.option("kmin") > 100 or self.option("kmin") < 1:
                    raise OptionError("k值要大于%s且小于%s", variables=(1, 100), code='12345')
        if self.option("kmax"):
            if type(self.option("kmax")) != int:
                raise OptionError("k_value一定要是整型数据", code="123345")
            else:
                if self.option("kmax") > 100 or self.option("kmax") < 1:
                    raise OptionError("k值要大于%s且小于%s", variables=(1, 100), code='12345')
        if not self.option('main_id'):
            raise OptionError("缺少main_id参数！")
        return True

    def set_group_list(self):
        """
        检查是否sample_list存在，如果存在的生成group.list，然后传给vcftools_filter_run中的keep参数
        :return:
        """
        if not self.option("sample_list"):
            return
        else:
            with open(self.work_dir + "/group.list", 'w') as w:
                for m in self.option("sample_list").split(','):
                    w.write("{}\n".format(m))

    def vcftools_filter_run(self):
        options = {
            "vcf_path": self.option("vcf_file"),
            "recode": self.option("recode"),
            "remove_indels": self.option("remove_indels"),
            "remove_filtered_all": self.option("remove_filtered_all"),
            "minDP": self.option("minDP"),
            "maxDP": self.option("maxDP"),
            "max_missing": self.option("max_missing"),
            "min_maf": self.option("min_maf"),
            "max_maf": self.option("max_maf")
        }
        if self.option("sample_list"):
            options.update({"keep": self.work_dir + "/group.list"})
        self.vcftools_filter.set_options(options)
        self.vcftools_filter.on("end", self.set_output, "vcf_filter")
        self.vcftools_filter.run()

    def vcftools_plink_run(self):
        options = {
            "recode_vcf_path": self.vcftools_filter.option("filter_vcf").prop['path'],
            "chr_set": self.option("chr_set"),  # 改物种具体的染色体与sca的个数
            "make_bed": True
        }
        if self.option("analysis_type") == "structure":
            options.update({"allow_extra_chr": True})
        self.vcftools_plink.set_options(options)
        self.vcftools_plink.on('end', self.set_output, "vcf_plink")
        self.vcftools_plink.run()

    def tree_run(self):
        self.tree.set_options({
            "recode_vcf_path": self.vcftools_filter.option("filter_vcf").prop['path'],
            "tree_type": self.option("tree_type"),
            "bs_trees": self.option("bs_trees")
        })
        self.tree.on('end', self.set_output, "tree")
        # self.tree.on('end', self.end)
        self.tree.run()

    def pca_run(self):
        self.pca.set_options({
            "bfile_dir": self.vcftools_plink.output_dir
        })
        self.pca.on("end", self.set_output, "pca")
        # self.pca.on("end", self.end)
        self.pca.run()

    def structure_run(self):
        self.structure.set_options({
            "pop_fam": self.vcftools_plink.output_dir + "/pop.fam",
            "pop_bed": self.vcftools_plink.output_dir + "/pop.bed",
            "k_min": self.option("kmin"),
            "k_max": self.option("kmax")
        })
        self.structure.on("end", self.set_output, "structure")
        # self.structure.on("end", self.end)
        self.structure.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'vcf_filter':
            self.linkdir(obj.output_dir, self.output_dir + "/vcf_filter")
        if event['data'] == 'vcf_plink':
            self.linkdir(obj.output_dir, self.output_dir + "/vcf_plink")
        if event['data'] == 'tree':
            self.linkdir(obj.output_dir, self.output_dir + "/tree")
            self.end()
        if event['data'] == 'pca':
            self.linkdir(obj.output_dir, self.output_dir + "/pca")
            self.end()
        if event['data'] == 'structure':
            self.linkdir(obj.output_dir, self.output_dir + "/structure")
            self.end()

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
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_db(self):
        api = self.api.api("dna_evolution.pop_analysis")
        self.logger.info("开始导表")
        if self.option("analysis_type") == "pca":
            api.add_sg_pop_pca_detail(self.option("main_id"), self.output_dir + '/pca/pop.pca.eigenvec')
            api.add_sg_scatter(self._sheet.id, self.option('main_id'), self.output_dir + '/pca/pop.pca.eigenvec')
            api.add_pca_bar(self._sheet.id, self.option('main_id'), self.output_dir + '/pca/pop.pca.eigenval')
        elif self.option("analysis_type") == "structure":
            api.add_sg_kvalue_detail(self.option("main_id"), self.output_dir + '/structure/structure/', self._sheet.id)
            api.add_structure_sg_curve(self._sheet.id, self.option("main_id"),
                                       self.output_dir + '/structure/cverror/cv.error')
        else:
            if self.option('tree_type') == 'nj':
                tree_path = self.output_dir + "/tree/tree/pop.nj.tree"
                api.add_pop_sg_tree(self._sheet.id, self.option('main_id'), tree_path)
            elif self.option('tree_type') == 'ml':
                tree_path = self.output_dir + "/tree/tree/pop.phylip.raxml.bestTree"
                api.add_pop_sg_tree(self._sheet.id, self.option('main_id'), tree_path, 'ml')
            else:
                tree_path = self.output_dir + "/tree/tree/bayes.nex.con"
                api.add_pop_sg_tree(self._sheet.id, self.option('main_id'), tree_path, 'bayes')
        self.logger.info("导表结束")

    def run(self):
        self.set_group_list()
        if self.option("analysis_type") == "structure":
            self.vcftools_plink = self.add_tool("dna_evolution.vcftools_plink")
            self.structure = self.add_module('dna_evolution.pop_structure')
            self.vcftools_filter.on("end", self.vcftools_plink_run)
            self.vcftools_plink.on("end", self.structure_run)
        elif self.option("analysis_type") == "pca":
            self.vcftools_plink = self.add_tool("dna_evolution.vcftools_plink")
            self.pca = self.add_tool('dna_evolution.gcta_pca')
            self.vcftools_filter.on("end", self.vcftools_plink_run)
            self.vcftools_plink.on("end", self.pca_run)
        else:
            self.tree = self.add_module('dna_evolution.tree_generic')
            self.vcftools_filter.on("end", self.tree_run)
        self.vcftools_filter_run()
        super(PopAnalysisWorkflow, self).run()

    def end(self):
        """
        这里后面要重新定义下文件名字
        :return:
        """
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        self.set_db()
        super(PopAnalysisWorkflow, self).end()
