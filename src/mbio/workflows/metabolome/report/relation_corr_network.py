# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'
# last_modifiy = 2022.1.10

from biocluster.workflow import Workflow
from bson.objectid import ObjectId
import os, glob
import types


class RelationCorrNetworkWorkflow(Workflow):
    """
    相关性热图分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RelationCorrNetworkWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metab_table", "type": "infile", "format": "metabolome.express,sequence.profile_table"},  # 代谢表达量表
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},  # 代谢desc表
            {"name": "metab_set_table", "type": "infile", "format": "sequence.profile_table"},  # 代谢集文件
            {"name": "group_detail", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {"name": "trans_exp_main_id", "type": "string"},  # 转录表达量表id
            {"name": "trans_geneset_main_id", "type": "string"},  # 转录基因集表id
            #参数选项
            {"name": "coefficient", "type": "string", "default": "pearson"},  # 相关性算法
            {"name": "padjust_method", "type": "string", "default": "fdr_bh"},  # 多重检验
            {"name": "log10", "type": "bool", "default": False},
            {"name": "sort", "type": "string", "default": "pvalue"},
            {"name": "top", "type": "int", "default": 200},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录
            {"name": "name", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.option("save_pdf", 0)
        self.metab_tool = self.add_tool("metabolome.select_table")  # 筛选代谢表
        self.trans_tool = self.add_tool("metabolome.relation.trans_select_table")  # 筛选转录表
        self.asso_corr_tool = self.add_tool("metabolome.metabset.asso_corr")
        self.p_adjust_tool = self.add_tool("metabolome.relation.pvalue_adjust")

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run RelationCorrNetwork workflow")
        self.select_tools = [self.metab_tool, self.trans_tool]
        self.on_rely(self.select_tools, self.run_asso_corr)
        self.asso_corr_tool.on('end', self.run_pvalue_adjust)
        self.p_adjust_tool.on('end', self.set_db)
        self.select_meta()
        self.select_trans()
        super(RelationCorrNetworkWorkflow, self).run()

    def select_meta(self):
        self.logger.info("start profile!")
        exp_profile = self.option("metab_table").prop["path"]
        options = {
            "origin_table": exp_profile,
            "select_genes": self.option("metab_set_table"),
            "st": "F",
            "group": self.option("group_detail"),
            "select_columns": "metab_id",
            "merge": "nomerge",
            "group_method": 0,
            "log10": self.option("log10")
        }
        self.metab_tool.set_options(options)
        self.metab_tool.run()

    def select_trans(self):
        self.logger.info("start profile!")
        options = {
            "group": self.option("group_detail"),
            "trans_exp_main_id": self.option("trans_exp_main_id"),
            "trans_geneset_main_id": self.option("trans_geneset_main_id"),
            "task_id" : "_".join(self._sheet.id.split("_")[0:2])
        }
        self.trans_tool.set_options(options)
        self.trans_tool.run()

    def run_asso_corr(self):  # 相关性计算
        self.logger.info("start run asso_corr !")
        metab_table = self.metab_tool.option("select_table")
        trans_table = self.trans_tool.option("select_table")
        metab_des = self.option("metab_desc").prop["path"]
        options = {
            'metab_table': metab_table,
            'asso_table': trans_table,
            'coefficient': self.option("coefficient"),
            'metab_trans': metab_des
        }
        self.asso_corr_tool.set_options(options)
        self.asso_corr_tool.run()

    def run_pvalue_adjust(self):
        self.logger.info("start run pvalue_adjust !")
        options = {
            "pvalue_table": self.asso_corr_tool.output_dir + "/pvalue.xls",
            "corr_table": self.asso_corr_tool.output_dir + "/corr.xls",
            "padjust_method": self.option("padjust_method"),
            "sort": self.option("sort"),
            "top": self.option("top")
        }
        self.p_adjust_tool.set_options(options)
        self.p_adjust_tool.run()

    def set_db(self):
        self.corr_network = self.api.api("metabolome.relation_corr_network")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="14700501")
        corr_result = self.link_file("sorted_express_correlation_info.txt", "sorted_express_correlation_info.txt")
        gene_id2name = None
        if os.path.exists(os.path.join(self.trans_tool.work_dir, "gene_id2name.xls")):
            gene_id2name = os.path.join(self.trans_tool.work_dir, "gene_id2name.xls")
        self.logger.info(gene_id2name)
        self.corr_network.add_relation_corr_network_detail(main_id, result_file=corr_result, gene_id2name=gene_id2name)

        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "relatcorrnetwork",
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()
            
    def link_file(self, oldfile, newfile):
        oldfile = os.path.join(self.p_adjust_tool.output_dir, oldfile)
        newfile = os.path.join(self.output_dir, newfile)
        if os.path.exists(newfile):
            os.remove(newfile)
        os.link(oldfile, newfile)
        return newfile

    def end(self):
        if self.option("save_pdf"):
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        relpath_rules = [
            [".", "", "关联分析相关性结果文件夹", 0, "150069"],
            ["express_correlation_info.txt", "txt", "基因与代谢物相关性系数总表", 0, "150071"],
            ["sorted_express_correlation_info.txt", "txt", "筛选后基因与代谢物相关性系数表", 0, "150071"]
        ]
        regexps = [
        ]
        result_dir.add_relpath_rules(relpath_rules)
        result_dir.add_regexp_rules(regexps)
        super(RelationCorrNetworkWorkflow, self).end()