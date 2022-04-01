# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
from biocluster.workflow import Workflow
from mbio.api.to_file.denovo import *


class MapAssessmentWorkflow(Workflow):
    """
    报告中计算比对质量评估时使用
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MapAssessmentWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "bed", "type": "string"},  # bed格式文件
            {"name": "bam", "type": "string"},  # bam格式文件,排序过的
            {"name": "fpkm", "type": "string"},  # 基因表达量表
            {"name": "insert_id", "type": "string"},  # 插入的主表的ID
            {"name": "express_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "analysis_type", "type": "string"},  # 分析类型
            {"name": "quality_satur", "type": "int"},  # 测序饱和度分析质量值
            {"name": "quality_dup", "type": "int"},  # 冗余率分析质量值
            {"name": "low_bound", "type": "int"},  # Sampling starts from this percentile
            {"name": "up_bound", "type": "int"},  # Sampling ends at this percentile
            {"name": "step", "type": "int"},  # Sampling frequency
            {"name": "min_len", "type": "int"}  # Minimum mRNA length (bp).
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.map_assess = self.add_module('denovo_rna.mapping.map_assessment')
        self.pca = None
        self.tools = [self.map_assess]

    def run(self):
        options = {
            'bed': self.option('bed'),
            # 'bam': self.option('bam'),
            'analysis': self.option('analysis_type'),
            'quality_satur': self.option('quality_satur'),
            'quality_dup': self.option('quality_dup'),
            'low_bound': self.option('low_bound'),
            'up_bound': self.option('up_bound'),
            'step': self.option('step'),
            # 'fpkm': self.option('fpkm').split(",")[0],
            'min_len': self.option('min_len')
            }
        if self.option("analysis_type") == "correlation":
            options["fpkm"] = self.option('fpkm').split(",")[0]
            self.pca = self.add_tool("meta.beta_diversity.pca")
            self.tools.append(self.pca)
            self.on_rely(self.tools, self.set_db)
            self.pca_run()
        if self.option("analysis_type") in ["saturation", "duplication"]:
            options["bam"] = self.option('bam')
            self.map_assess.on('end', self.set_db)
        print(options)
        self.map_assess.set_options(options)
        # self.map_assess.on('end', self.set_db)
        self.map_assess.run()
        self.output_dir = self.map_assess.output_dir
        super(MapAssessmentWorkflow, self).run()

    def pca_run(self):
        self.pca.set_options({
            "otutable": self.option('fpkm').split(",")[0]
        })
        self.pca.run()

    def set_db(self):
        api_mapping = self.api.denovo_rna_mapping
        print(self.option("insert_id"))
        if self.option("analysis_type") == "duplication":
            dup_path = self.output_dir + "/dup"
            print(dup_path)
            if os.path.isfile(dup_path):
                raise Exception("找不到报告文件夹:{}".format(dup_path))
            api_mapping.add_duplication_detail(dup_path, self.option("insert_id"))
        if self.option("analysis_type") == "saturation":
            satur_path = self.output_dir + "/satur"
            if os.path.isfile(satur_path):
                raise Exception("找不到报告文件夹:{}".format(satur_path))
            api_mapping.add_rpkm_detail(satur_path, self.option("insert_id"))
            # api_mapping.add_rpkm_box(satur_path, self.option("insert_id"))
            api_mapping.add_rpkm_curve(satur_path, self.option("insert_id"))
        if self.option("analysis_type") == "correlation":
            correlation_path = self.output_dir + "/correlation"
            pca_path = self.pca.output_dir + "/pca_importance.xls"
            pca_rotation = self.pca.output_dir + "/pca_rotation.xls"
            pca_sites = self.pca.output_dir + "/pca_sites.xls"
            if os.path.isfile(correlation_path):
                raise Exception("找不到报告文件夹:{}".format(correlation_path))
            api_mapping.add_correlation_detail(correlation_path, self.option("insert_id"), updata_tree=True)
            api_mapping.add_pca(pca_path, self.option("insert_id"))
            api_mapping.add_pca_rotation(pca_rotation, "sg_denovo_correlation_pca_rotation", self.option("insert_id"))
            api_mapping.add_pca_rotation(pca_sites, "sg_denovo_correlation_pca_sites", self.option("insert_id"))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "dir", "结果输出目录"],
            # ["./coverage/", "dir", "基因覆盖度分析输出目录"],
            ["./dup/", "dir", "冗余序列分析输出目录"],
            ["./satur/", "dir", "测序饱和度分析输出目录"],
            ["./bam_stat.xls", "xls", "bam格式比对结果统计表"]
        ])
        result_dir.add_regexp_rules([
            [r".*pos\.DupRate\.xls", "xls", "比对到基因组的序列的冗余统计表"],
            [r".*seq\.DupRate\.xls", "xls", "所有序列的冗余统计表"],
            [r".*eRPKM\.xls", "xls", "RPKM表"],
            [r".*cluster_percent\.xls", "xls", "饱和度作图数据"],
            [r".correlation_matrix*\.xls", "xls", "相关系数矩阵"],
            [r".hcluster_tree*\.xls", "xls", "样本间相关系数树文件"]
        ])
        # print self.get_upload_files()
        super(MapAssessmentWorkflow, self).end()
