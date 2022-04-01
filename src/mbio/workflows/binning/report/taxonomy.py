# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2018/12/21'

from biocluster.workflow import Workflow
import os
import json


class TaxonomyWorkflow(Workflow):
    """
    16s taxonomy
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TaxonomyWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {'name': 'database', 'type': 'string'},  # 数据库选择
            {'name': 'query', 'type': 'infile', 'format': 'sequence.fasta'},  # 查询序列
            {'name': 'revcomp', 'type': 'bool', 'default': False},  # 序列是否翻转
            {'name': 'confidence', 'type': 'float', 'default': 0.7},  # 置信度值
            {'name': 'ref_fasta', 'type': 'infile', 'format': 'sequence.fasta'},  # 参考fasta序列
            {'name': 'ref_taxon', 'type': 'infile', 'format': 'taxon.seq_taxon'},  # 参考taxon文件
            {'name': 'taxon_file', 'type': 'outfile', 'format': 'taxon.seq_taxon'},  # 输出序列的分类信息文件
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.new_taxon = ""
        self.format_taxon = self.add_tool("meta.filecheck.format_taxon")
        self.tax = self.add_module("annotation.meta_tax")

    def run(self):
        # self.run_taxon_format()
        self.set_db()
        super(TaxonomyWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_path = self.api.api("metagenomic.common_api")
        api_path._project_type = "metag_bin" # 需配置
        # 导入16s注释表
        taxonomy_params = {
            "query": "bin1",
            "db": "silva128"
        }
        taxonomy_params = json.dumps(taxonomy_params, sort_keys=True, separators=(',', ':'))
        main_id = api_path.add_main("identif_16s", name="16S_origin", params=taxonomy_params)  # 注意参数中需要task_id 和 project_sn
        api_path.add_main_detail("/mnt/ilustre/users/sanger-dev/sg-users/guhaidong/binning/blasr/blasr.xls", "identif_16s_detail",
                                 main_id, "genome_id,taxon,identify", has_head=True, main_name="s16_id")
        # 导入pocp结果
        pocp_params = {
            "query": "bin1",
            "ref": "user_defined", # "bin2"
            "accession": "XXX"    # ref_path + file_dir_id
        }
        pocp_params = json.dumps(pocp_params, sort_keys=True, separators=(',', ':'))
        main_id = api_path.add_main("identif_genus", name="Genus_origin", params=pocp_params)
        api_path.add_main_detail("/mnt/ilustre/users/sanger-dev/sg-users/guhaidong/binning/pocp/pocp_result.xls", "identif_genus_detail",
                                 main_id, "genome_id,ref,ref_genus,prop", has_head=True, main_name="genus_id")
        # 导入ani结果
        ani_params = {
            "query": "bin1",
            "ref": "bin2"
        }
        ani_params = json.dumps(ani_params, sort_keys=True, separators=(',', ':'))
        main_id = api_path.add_main("identif_species", name="Spe_origin", params=ani_params)
        api_path.add_main_detail("/mnt/ilustre/users/sanger-dev/sg-users/guhaidong/binning/ani/ani_result.xls", "identif_species_detail",
                                 main_id, "genome_id,ref,ref_species,ani", has_head=True, main_name="species_id")
        # 导入共线性分析结果
        col_params = {
            "query": "bin1",
            "ref": "user_defined",
            "ref_path": "",
            "file_dir_id": ""
        }
        col_params = json.dumps(col_params, sort_keys=True, separators=(',', ':'))
        main_id = api_path.add_main("colline", name="Colline_origin", params=col_params, others={"cir_path": "s3://binning/"})  # cir_path为文件路径
        api_path.add_main_detail("/mnt/ilustre/users/sanger-dev/sg-users/guhaidong/binning/colline/400786_block_new.xls", "colline_detail",
                                 main_id, "block,location,strand,start,end,len,genome_id,loc_start,loc_end", has_head=True, main_name="circos_id")
        # 注意结果表的第一列需要去空格，第二列需要去掉前两个字符
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.tax.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果目录"],
        ])
        super(TaxonomyWorkflow, self).end()

    def run_taxon_format(self):
        self.new_taxon = self.format_taxon.output_dir + '/taxon.new.xls'
        if self.option("database") == "custom_mode":
            opts = {"in_taxon_table": self.option('ref_taxon')}
            self.format_taxon.set_options(opts)
            self.format_taxon.on('end', self.run_taxon)
            self.format_taxon.run()
        else:
            self.logger.info("have database")
            self.run_taxon()

    def run_taxon(self):
        opts = {
            "fasta": self.option("query"),
            "revcomp": self.option("revcomp"),
            "confidence": self.option("confidence"),
            "database": self.option("database")
        }
        if self.option("database") == "nt":
            opts.update({
                "query_type": "nucl",
                "blast": "blastn",
                "num_alignment": 100
            })
        elif self.option("database") == "custom_mode":
            opts.update({
                "ref_fasta": self.option("ref_fasta"),
                "ref_taxon": self.new_taxon
            })
        self.tax.set_options(opts)
        self.tax.on("end", self.set_db)
        self.tax.run()
