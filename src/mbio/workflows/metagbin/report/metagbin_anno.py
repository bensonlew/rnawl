# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# __modify__ = '2019.01.31'

from biocluster.workflow import Workflow
import os
import json
import datetime
from mbio.packages.metagbin.common_function import link_dir


class MetagbinAnnoWorkflow(Workflow):
    """
    宏基因组binning模块
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MetagbinAnnoWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "gene_seq", "type": "infile", "format": "sequence.fasta"},  # 序列
            {"name": "gene_gff", "type": "infile", "format": "gene_structure.gff3"},  # 基因预测gff文件
            {"name": "sample", "type": "string", "default": ""},  # 样品名称
            {"name": "database", "type": "string", 'default': ""},  #数据库名称 {\"card\":\"v3.0.9\",\"cazy\":\"V8\",\"kegg\":\"v94.2\",\"nr\":\"v20200604\"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.anno = self.add_module("metagbin.annotation")
        self.bin_name = self.option("sample")
        self.remote_dir = self._sheet.output

    def run(self):
        self.run_anno()
        super(MetagbinAnnoWorkflow, self).run()

    def run_anno(self):
        opts ={
            "gene_seq":self.option("gene_seq"),
            "sample":self.option("sample"),
            "gene_gff":self.option("gene_gff"),
        }
        if self.option("database") != "": ## 新版本数据库增加此参数用于判断数据库版本兼容注释新老版本和新老任务
            database = json.loads(self.option("database"))
            if 'card' in database:
                if database['card'] in ['v3.0.9']: ## 新版本默认跑drug_class，需要兼容老版本
                    database_list = "nr_v20200604,card_v3.0.9,eggnog,kegg_v94.2,cazy_v8"
                else:
                    opts["category"] = "aro_category"
                    database_list = "nr,card,eggnog,kegg,cazy"
            else:
                opts["category"] = "aro_category"
                database_list = "nr,card,eggnog,kegg,cazy"
        else:
            opts["category"] = "aro_category"
            database_list = "nr,card,eggnog,kegg,cazy"
        opts["database_list"] = database_list
        self.anno.set_options(opts)
        self.anno.on('end',self.set_db)
        self.anno.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        self.set_output()
        task_id = '_'.join(self._sheet.id.split("_")[0:2])
        bin_name = self.bin_name
        pi_path = self.api.api("metagbin.common_api")
        anno_nr = self.api.api("metagbin.anno_nr")
        nr_params = {"genome_id": bin_name, "nr": "Diamond", 'task_type': 2, 'submit_location': "anno_nr"}
        main_table_name = 'AnnoNr_' + bin_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        nr_id = anno_nr.add_anno_nr(params=nr_params, name=main_table_name, task_id=task_id)
        anno_nr.add_anno_nr_detail(nr_id, bin_name,self.output_dir + '/' + bin_name + "/NR/" + bin_name + "_anno_nr.xls")
        pi_path.add_sg_status("sg_status", table_name=main_table_name, desc='anno_nr注释', submit_location='anno_nr',params=nr_params, table_id=nr_id, genome_id=bin_name, type_name='anno_nr')
        anno_cog = self.api.api("metagbin.anno_cog")
        cog_params = {"genome_id": bin_name, "cog": "Diamond", 'task_type': 2, 'submit_location': "anno_cog"}
        main_table_name = 'AnnoCog_' + bin_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        cog_id = anno_cog.add_anno_cog(params=cog_params, name=main_table_name, task_id=task_id)
        anno_cog.add_anno_cog_detail(cog_id, bin_name,self.output_dir + '/' + bin_name + "/COG/" + bin_name + "_cog_anno.xls")
        pi_path.add_sg_status("sg_status", table_name=main_table_name, desc='anno_cog注释', submit_location='anno_cog',params=cog_params, table_id=cog_id, genome_id=bin_name, type_name='anno_cog')
        anno_kegg = self.api.api("metagbin.anno_kegg")
        kegg_params = {"genome_id": bin_name, "kegg": "Diamond", 'task_type': 2, 'submit_location': "anno_kegg"}
        main_table_name = 'AnnoKegg_' + bin_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        kegg_id = anno_kegg.add_anno_kegg(self.remote_dir, "/KEGG/", "_kegg_pathway_img",params = kegg_params, task_id=task_id, name=main_table_name)
        anno_kegg.add_anno_kegg_detail(kegg_id, bin_name, self.output_dir + '/' + bin_name + "/KEGG/" + bin_name + "_kegg_anno.xls")
        anno_kegg.add_anno_kegg_level(kegg_id, bin_name, self.output_dir + '/' + bin_name + "/KEGG/" + bin_name + "_kegg_level_stat.xls")
        pi_path.add_sg_status("sg_status", table_name=main_table_name, desc='anno_kegg注释', submit_location='anno_kegg',params=kegg_params, table_id=kegg_id, genome_id=bin_name, type_name='anno_kegg')
        anno_cazy = self.api.api("metagbin.anno_cazy")
        cazy_params = {"genome_id": bin_name, "cazy": "hmmscan", 'task_type': 2, 'submit_location': "annocazy"}
        main_table_name = 'AnnoCazy_' + bin_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        cazy_id = anno_cazy.add_anno_cazy(params=cazy_params, name=main_table_name, task_id=task_id)
        anno_cazy.add_anno_cazy_detail(cazy_id, bin_name, self.output_dir + '/' + bin_name + "/CAZY/" + bin_name + "_anno_cazy.xls")
        pi_path.add_sg_status("sg_status", table_name=main_table_name, desc='anno_cazy注释', submit_location='annocazy',params=cazy_params, table_id=cazy_id, genome_id=bin_name, type_name='annocazy')
        api_card = self.api.api("metagbin.anno_card")
        card_params = {"genome_id": bin_name, "card": "Diamond", 'task_type': 2, 'submit_location': "annocard"}
        main_table_name = 'AnnoCard_' + bin_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        if self.option("database") != "":
            card_id = api_card.add_anno_card(self.output_dir + '/' + bin_name + "/CARD", main=True, name=main_table_name, genome_id=bin_name, params=card_params, task_id=task_id, database_version="new")
        else:
            card_id = api_card.add_anno_card(self.output_dir + '/' + bin_name + "/CARD", main=True, name=main_table_name, genome_id=bin_name, params=card_params, task_id=task_id, database_version="old")
        pi_path.add_sg_status("sg_status", table_name=main_table_name, desc='anno_card注释', submit_location='annocard', type_name='annocard',params=card_params, table_id=card_id, genome_id=bin_name)
        api_sum = self.api.api("metagbin.anno_summary")
        api_stat = self.api.api("metagbin.anno_stat")
        sum_params = {"genome_id": bin_name, 'task_type': 2, 'submit_location': "annosummary"}
        main_table_name = 'AnnoSummary_' + bin_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        sum_id = api_sum.add_anno_summary(params=sum_params,name=main_table_name, task_id=task_id)
        pi_path.add_sg_status("sg_status", table_name=main_table_name, desc='anno_summary注释', submit_location='annosummary',params=sum_params, table_id=sum_id, genome_id=bin_name, type_name='annosummary')
        if self.option("database") != "":
            api_sum.add_anno_summary_detail(sum_id, bin_name, self.output_dir + '/' + bin_name + "/Summary/" +bin_name + "_anno_summary.xls", database_version="new")
        else:
            api_sum.add_anno_summary_detail(sum_id, bin_name, self.output_dir + '/' + bin_name + "/Summary/" +bin_name + "_anno_summary.xls", database_version="old")
        api_stat.add_anno_stat_detail(self.option("main_id"), bin_name, self.output_dir + '/' + bin_name + "/Summary/" + bin_name + "_anno_stat.xls")
        self.end()

    def set_output(self):
        link_dir(self.anno.output_dir,self.output_dir + '/' + self.bin_name)
        self.run_move()

    def run_move(self):
        if os.path.exists(self.output_dir + '/' + self.bin_name + "/annotation/KEGG/kegg_pathway_img.tar.gz"):
            os.remove(self.output_dir + '/' + self.bin_name + "/annotation/KEGG/kegg_pathway_img.tar.gz")

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        repaths = []
        regexps = []
        repaths += [
            [".", "", "基因组注释结果目录"],
            ["/%s" %self.bin_name, "", "%s分析注释结果目录" %self.bin_name, 0],
            ["/%s/NR" %self.bin_name, "", "%s基因组的NR注释结果目录" %self.bin_name, 0],
            ["/%s/COG" %self.bin_name, "", "%s基因组的COG注释结果目录" %self.bin_name, 0],
            ["/%s/KEGG" %self.bin_name, "", "%s基因组的KEGG注释结果目录" %self.bin_name, 0],
            ["/%s/CAZY" %self.bin_name, "", "%s基因组的CAZY注释结果目录" %self.bin_name, 0],
            ["/%s/CARD" %self.bin_name, "", "%s基因组的CARD注释结果目录" %self.bin_name, 0],
            ["/%s/Summary" %self.bin_name, "", "%s基因组的注释结果汇总目录" %self.bin_name , 0],
            ]
        regexps +=[
            [r"/%s/NR/.+_anno_nr.xls" %self.bin_name, "", "%s基因组的NR注释结果表" %self.bin_name, 0],
            [r"/%s/COG/.+_cog_summary.xls" %self.bin_name, "", "%s基因组的COG注释结果分类统计表" %self.bin_name, 0],
            [r"/%s/COG/.+_cog_anno.xls" %self.bin_name, "", "%s基因组的COG注释结果表" %self.bin_name, 0],
            [r"/%s/KEGG/.+_kegg_pathway_img" % self.bin_name, "", "%s基因组的KEGG通路图文件目录" % self.bin_name, 0],
            [r"/%s/KEGG/.+_kegg_anno.xls  " % self.bin_name, "", "%s基因组的KEGG注释结果文件" % self.bin_name, 0],
            [r"/%s/KEGG/.+_kegg_level_stat.xls" % self.bin_name, "", "%s基因组的KEGG注释level统计文件" % self.bin_name, 0],
            [r"/%s/CAZY/.+_anno_cazy.xls" %self.bin_name, "", "%s基因组的CAZY注释分类结果表" %self.bin_name, 0],
            [r"/%s/CAZY/.+_cazy_class_stat.xls" %self.bin_name, "", "%s基因组的CAZY注释class分类结果表" %self.bin_name, 0],
            [r"/%s/CAZY/.+_cazy_family_stat.xls" %self.bin_name, "", "%s基因组的CAZY注释family分类结果表" %self.bin_name, 0],
            [r"/%s/CAZY/.+_cazy_parse_anno.xls" %self.bin_name, "", "%s基因组的CAZY注释结果表" %self.bin_name, 0],
            [r"/%s/CARD/.+_card_align.xls" %self.bin_name, "", "%s基因组的CARD比对结果表" %self.bin_name, 0],
            [r"/%s/CARD/.+_card_anno.xls" %self.bin_name, "", "%s基因组的CARD注释结果表" %self.bin_name, 0],
            [r"/%s/CARD/.+_card_category.xls" %self.bin_name, "", "%s基因组的CARD注释分类结果表"%self.bin_name, 0],
            [r"/%s/Summary/.+_anno_summary.xls" % self.bin_name, "", "%s基因组的注释结果汇总表" % self.bin_name, 0],
            [r"/%s/Summary/.+_anno_stat" % self.bin_name, "", "%s基因组的基因注释统计表" % self.bin_name, 0],
        ]
        result_dir.add_relpath_rules(repaths)
        result_dir.add_regexp_rules(regexps)
        super(MetagbinAnnoWorkflow, self).end()