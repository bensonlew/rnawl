# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re,json
from bson import ObjectId
from mbio.packages.meta.common_function import get_level_id


class FigSaveAgent(Agent):
    """

    """
    def __init__(self, parent):
        super(FigSaveAgent, self).__init__(parent)
        options = [
            {"name": "run_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "table_id", "type": "string"},
            {"name": "table_name", "type": "string"},
            {"name": "submit_loc", "type": "string"},
            {"name": "project", "type": "string"},
            {"name": "proj_dir", "type": "string"},
            {"name": "interaction", "type": "int", "default": 1},
            {"name": "main_table", "type": "string"},
            {"name": "group_id_detail", "type": "string"},
            {"name": "otu_id", "type": "string"},
        ]
        self.add_option(options)

    def check_options(self):
        return True

    def set_resource(self):
        self._cpu = 2
        if not self.option("interaction"):
            self._memory = '20G'
        else:
            self._memory = '10G'

    def end(self):
        super(FigSaveAgent, self).end()

class FigSaveTool(Tool):
    def __init__(self, config):
        super(FigSaveTool, self).__init__(config)
        base_dir = self.config.SOFTWARE_DIR + '/bioinfo/figsave'
        py_ld = self.config.SOFTWARE_DIR + '/program/Python/lib:' + base_dir + "/lib"
        self.python = "/miniconda2/bin/python"
        self.fig_save = base_dir + '/fig_save'
        self.set_environ(PATH=self.fig_save, PYTHONPATH=base_dir + '/lib/python2.7/site-packages',
                         LD_LIBRARY_PATH=py_ld)
        self.common_api = self.api.api("metagenomic.common_api")
        self.common_api._project_type = self.option("project")
        self.run_info = self.common_api.db['sg_task'].find_one({"task_id": self.option("task_id")})
        self.meth2out = {}
        if self.option("group_id_detail"):
            self.group_id_detail = {v: k for k, v in eval(self.option("group_id_detail")).items()}
            self.group_id_detail["All"] ="All"
            self.group_id_detail["all"] = "All"
            self.group_id_detail["ALL"] = "All"
        self.level ={
            "1":"Domain",
            "2": "Kingdom",
            "3": "Phylum",
            "4": "Class",
            "5": "Order",
            "6": "Family",
            "7": "Genus",
            "8": "Species",
            "9": "otu",
        }

    def run(self):
        """
        运行
        :return:
        """
        super(FigSaveTool, self).run()
        self.tmp_dir = os.path.join(self.work_dir, "tmp_out")
        if not os.path.exists(self.tmp_dir):
            os.mkdir(self.tmp_dir)
        cmd = "{} {}/fig2pdf -proj {} -task_id {} -tmp_out {}".format(self.python, self.fig_save,
                                                                         self.option("project"),
                                                                         self.option("task_id"), self.tmp_dir)
        if not self.option("interaction"):
            cmd += " -interaction 0 -png"
            #cmd += " -interaction 0"
        else:
            self.update_status(self.option("table_id"),self.option("main_table"),{"status":"end"})
            cmd += " -sub_loc {} -table_id {} -table_name {}".format(self.option("submit_loc"),
                                                                    self.option("table_id"),
                                                                    self.option("table_name"))
        cmd += ' -pk FigSave.pk '
        #proj_kws = {"proj_type": self.run_info["type"], "mix": 'F'}
        #if "mix" in self.run_info:
        #    proj_kws["mix"] = self.run_info["mix"]
        #cmd += " -proj_kws '{}'".format(json.dumps(proj_kws))
        if self.config.DBVersion:
            cmd += ' -db_version ' + str(self.config.DBVersion)
        if not self.command_done('fig_save'):
            command = self.add_command("fig_save", cmd).run()
            self.wait(command)
        self.set_output()
        self.end()

    def command_done(self, name):
        o_file = name + '.o'
        if os.path.exists(o_file):
            filesize = os.path.getsize(o_file)
            offset = 200 if filesize > 200 else filesize
            with open(o_file, 'rb') as r:
                r.seek(-offset, 2)
                tail_info = r.read()
                if "exitcode:0" in tail_info:
                    return True
        return False

    def update_status(self, obj_id, sub_loc, update_info):
        #sg_status = self.common_api.db["sg_status"].find_one({"submit_location": sub_loc})
        main_col = self.common_api.db[self.option("project")]
        main_col.update_one({"_id": ObjectId(obj_id)}, {"$set": update_info})

    def set_output(self):
        fig_info = os.path.join(self.fig_save, "conf/meta_fig_info.txt")
        with open(fig_info, 'r') as r:
            for line in r:
                line = line.strip().split('\t')
                if len(line) == 8:
                    if not self.option("interaction"):
                        self.meth2out[line[0]] = line[4].replace(' ', '')
                    else:
                        self.logger.info("line：{}".format(line))
                        self.meth2out[line[0]] = line[4].split(' ')[1]
        if not self.option("interaction"):
            with open(os.path.join(self.work_dir, "workflow_tables.txt"), 'r') as r:
                for line in r:
                    line = line.strip().split('\t')
                    table_id = str(line[1])
                    table_name = line[2]
                    self.out(table_id, table_name, line[0])
        else:
            self.out(self.option("table_id"), self.option("table_name"),
                        self.option("submit_loc"))

    def out(self, table_id, table_name, sub_loc):
        image_opt = os.path.join(self.tmp_dir, table_id + "_" + sub_loc + '.sg_image_options.txt')
        if os.path.exists(image_opt):
            with open(image_opt, 'r') as r:
                for line in r:
                    line = line.strip().split('\t')
                    params = eval(line[1])
                    if not self.option("interaction"):
                        if "method" in params:
                            if params["method"] == "get_alpha_rank_abundance":
                                if params["otu_id"] != self.option("otu_id"):
                                    out_path = "OtuTaxon_summary/Otu_Rank_Abundance曲线图.pdf"
                                else:
                                    out_path = "OtuTaxon_summary_depth/Otu_Rank_Abundance曲线图.pdf"
                            elif params["method"] == "showOtuPanCoreCurveByOtuPanCoreId":
                                if params["type"] == 1:
                                    out_path = "Pan_Core/" + self.group_id_detail[
                                        params["color_specimen_ids"].keys()[0]] + "/PanCoreOTU/pan曲线图.pdf"
                                else:
                                    out_path = "Pan_Core/" + self.group_id_detail[
                                        params["color_specimen_ids"].keys()[0]] + "/PanCoreOTU/core曲线图.pdf"
                            elif params["method"] == "showAlphaDiversityBarErrorByAlphaDiversityId":
                                out_path = "Alpha_diversity/" + self.group_id_detail[
                                    params["specimen_specimen_ids"].keys()[0]] + "/Estimators/" + params[
                                               "index_type"] + "指数柱形图.pdf"
                            elif params["method"] == "get_alpha_ttest_bar":
                                out_path = "Alpha_diversity/" + self.group_id_detail[
                                    params["color_specimen_ids"].keys()[0]] + "/EstTTest/" + params[
                                               "index_type"] + "指数检验差异检验柱形图.pdf"
                            elif params["method"] == "showAlphaTtestBoxByAlphaTtestId":
                                out_path = "Alpha_diversity/" + self.group_id_detail[
                                    params["color_specimen_ids"].keys()[0]] + "/EstTTest/" + params[
                                               "index_type"] + "指数检验差异检验箱式图.pdf"
                            elif params["method"] == "showAlphaRarefactionCurveByRarefactionId":
                                out_path = "Alpha_diversity/" + self.group_id_detail[
                                    params["specimen_specimen_ids"].keys()[
                                        0]] + "/Rarefaction/" + "各样本的{}指数稀释性曲线图.pdf".format(params["index_type"])
                            elif params["method"] == "showOtuVennByOtuVennId":
                                out_path = "Composition/" + self.group_id_detail[
                                    params["color_specimen_ids"].keys()[0]] + "/Venn/Venn图.pdf"
                            elif params["method"] == "showOtuGroupAnalyseBarByOtuId":
                                out_path = "Composition/" + self.group_id_detail[
                                    params["specimen_specimen_ids"].keys()[0]] + "/BarPie_" + self.level[
                                               str(params["level_id"])] + "/群落柱形图.pdf"
                            elif params["method"] == "showHcHeatmapByHcHeatmapId":
                                out_path = "Composition/" + self.group_id_detail[
                                    params["specimen_specimen_ids"].keys()[0]] + "/Heatmap/群落Heatmap图.pdf"
                            elif params["method"] == "get_otu_cicros":
                                out_path = "Composition/" + self.group_id_detail[
                                    params["color_specimen_ids"].keys()[0]] + "/Circos/Circos图.pdf"
                            elif params["method"] == "get_beta_specimen_distance_force":
                                out_path = "Beta_diversity/" + self.group_id_detail[
                                    params["color_specimen_ids"].keys()[0]] + "/Hcluster/样本层级聚类分析图.pdf"
                            elif params["method"] == "showBetaMultiAnalysisPcaScatterMarkByMultiAnalysisId":
                                out_path = "Beta_diversity/" + self.group_id_detail[
                                    params["shape_specimen_ids"].keys()[0]] + "/Pca/PCA分析散点图.pdf"
                            elif params["method"] == "get_beta_multi_pca_box_plot":
                                out_path = "Beta_diversity/" + self.group_id_detail[
                                    params["color_specimen_ids"].keys()[0]] + "/Pca/PCA分析箱式图.pdf"
                            elif params["method"] == "showBetaMultiAnalysisPcoaScatterMarkByMultiAnalysisId":
                                out_path = "Beta_diversity/" + self.group_id_detail[
                                    params["shape_specimen_ids"].keys()[0]] + "/Pcoa/PCoA分析散点图.pdf"
                            elif params["method"] == "get_beta_multi_pcoa_box_plot":
                                out_path = "Beta_diversity/" + self.group_id_detail[
                                    params["shape_specimen_ids"].keys()[0]] + "/Pcoa/PCoA分析箱式图.pdf"
                            elif params["method"] == "showBetaMultiAnalysisNmdsScatterMarkByMultiAnalysisId":
                                out_path = "Beta_diversity/" + self.group_id_detail[
                                    params["shape_specimen_ids"].keys()[0]] + "/Nmds/NMDS分析散点图.pdf"
                            elif params["method"] == "showBetaMultiAnosimBoxPlotByMultiAnosimId":
                                out_path = "Beta_diversity/" + self.group_id_detail[
                                    params["group_id"]] + "/AnosimBox/ANOSIM分析箱式图.pdf"
                            elif params["method"] == "get_species_difference_multiple_bar":
                                level = get_level_id("sg_species_difference_check", params["result_table_id"],
                                                     "level_id")
                                out_path = "DiffGroup/" + self.group_id_detail[
                                    params["color_specimen_ids"].keys()[0]] + "/DiffStatMultiple_" + self.level[
                                               str(level)] + "/多组比较差异检验柱形图.pdf"
                            elif params["method"] == "get_species_difference_two_group_bar_error":
                                level = get_level_id("sg_species_difference_check", params["result_table_id"],
                                                     "level_id")
                                out_path = "DiffGroup/" + self.group_id_detail[
                                    params["color_specimen_ids"].keys()[0]] + "/DiffStatTwoGroup_" + self.level[
                                               str(level)] + "/两组比较差异检验柱形图.pdf"
                            elif params["method"] == "lefse_tree":
                                out_path = "DiffGroup/" + self.group_id_detail[
                                    params["color_specimen_ids"].keys()[0]] + "/LEfSe/LEfSe多级物种层级树图.pdf"
                            elif params["method"] == "lefse_bar":
                                out_path = "DiffGroup/" + self.group_id_detail[
                                    params["color_specimen_ids"].keys()[0]] + "/LEfSe/LDA判别结果图.pdf"
                            elif params["method"] == "showBetaMultiAnalysisRdaCcaScatterMarkByMultiAnalysisId":
                                my_param = get_level_id("sg_beta_multi_analysis", params["result_table_id"],
                                                        "params")
                                group_id = eval(my_param)["group_id"]
                                out_path = "Env_analysis/" + self.group_id_detail[group_id] + "/Rda/RDA_CCA分析结果图.pdf"
                            elif params["method"] == "showPersonCorrlationByPearsonCorrelationId":
                                level = get_level_id("sg_species_env_correlation", params["result_table_id"],
                                                     "level_id")
                                my_param = get_level_id("sg_species_env_correlation", params["result_table_id"],
                                                        "params")
                                group_id = eval(my_param)["group_id"]
                                out_path = "Env_analysis/" + self.group_id_detail[group_id] + "/SpearmanCorrelation_" + \
                                           self.level[str(level)] + "/相关性Heatmap图.pdf"
                            elif params["method"] == "get_corr_network":
                                group_id = get_level_id("sg_corr_network", params["result_table_id"], "group_id")
                                out_path = "CorrNetworkSpearmanGenus/" + self.group_id_detail[
                                    str(group_id)] + "/corr_network_calc/单因素相关性网络图.pdf"
                            elif params["method"] == "showCorrNetworkHeatmapByCorrNetworkId":
                                group_id = get_level_id("sg_corr_network", params["result_table_id"], "group_id")
                                out_path = "CorrNetworkSpearmanGenus/" + self.group_id_detail[
                                    str(group_id)] + "/heatmap/物种相关性Heatmap图.pdf"
                            elif params["method"] == "get_specimen_step":
                                out_path = "QC_stat/序列长度分布图.pdf"
                            else:
                                out_path = self.meth2out[params["method"]].format(**params)
                        else:
                            out_path = self.meth2out[params["method"]].format(**params)
                    else:
                        out_path = self.meth2out[params["method"]].format(**params)
                    out_path = out_path.replace('__', '_')
                    out_path = os.path.join(self.output_dir, out_path)
                    if not os.path.exists(os.path.dirname(out_path)):
                        os.makedirs(os.path.dirname(out_path))
                    if os.path.exists(out_path):
                        os.remove(out_path)
                    if os.path.exists(os.path.join(self.tmp_dir, line[0] + '.pdf')):
                        os.link(os.path.join(self.tmp_dir, line[0] + '.pdf'), out_path)
        self.update_status(table_id, sub_loc, {"pdf_saved": 1})
