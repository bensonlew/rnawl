# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os, re
import json
from bson import ObjectId
import tarfile

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
        ]
        self.add_option(options)

    def check_options(self):
        pass
        return True

    def set_resource(self):
        if self.option("interaction") == 0:
            self._cpu = 6
        else:
            self._cpu = 2
        self._memory = '10G'

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
            # [".", "", "结果输出目录"],
        # ])
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
        # self.subloc2col = {}
        # self.get_subloc2col()

    def run(self):
        """
        运行
        :return:
        """
        super(FigSaveTool, self).run()
        self.tmp_dir = os.path.join(self.work_dir, "tmp_out")
        cmd = "{} {}/fig2pdf -proj {} -task_id {} -tmp_out {}".format(self.python, self.fig_save,
                                                                         self.option("project"),
                                                                         self.option("task_id"), self.tmp_dir)
        if not self.option("interaction"):
            cmd += " -interaction 0"
        else:
            self.update_status(self.option("table_id"), self.option("submit_loc"), {"status": "end"})
            cmd += " -sub_loc {} -table_id {} -table_name {}".format(self.option("submit_loc"),
                                                                    self.option("table_id"),
                                                                    self.option("table_name"))
        cmd += " -pk FigSave.pk"
        if self.config.DBVersion:
            cmd += " -db_version {}".format(self.config.DBVersion)
        if not self.command_done('fig_save'):
            command = self.add_command("fig_save", cmd).run()
            self.wait(command)
        self.set_output()
        self.dir_gz()
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
        # print(obj_id)
        # print(sub_loc)
        # print(self.common_api.db)
        # main_col = self.common_api.db[self.subloc2col[sub_loc]]
        main_col = self.common_api.db[self.option("project")]
        main_col.update_one({"_id": ObjectId(obj_id)}, {"$set": update_info})

    """
    def get_subloc2col(self):
        fig_info = os.path.join(self.fig_save, "conf/metag_fig_info.txt")
        with open(fig_info, 'r') as r:
            for line in r:
                line = line.strip().split('\t')
                if len(line) == 6:
                    if not self.option("interaction"):
                        self.meth2out[line[0]] = line[4].replace(' ', '')
                    else:
                        self.meth2out[line[0]] = line[4].split(' ')[1]
                    self.subloc2col[line[1]] = line[2]
    """

    def set_output(self):
        fig_info = os.path.join(self.fig_save, "conf/metag_fig_info.txt")
        with open(fig_info, 'r') as r:
            for line in r:
                line = line.strip().split('\t')
                if len(line) == 7:
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

    def get_specimen_name(self, specimen_id):
        try:
            table = self.common_api.db['data_stat_detail'].find_one({"_id": ObjectId(specimen_id)})
            try:
                sample_name = table["new_name"]
            except:
                sample_name = table["specimen_name"]
        except:
            table = self.common_api.db['data_stat_detail'].find_one({"origin_id": ObjectId(specimen_id)})
            try:
                sample_name = table["new_name"]
            except:
                sample_name = table["specimen_name"]
        return sample_name

    def out(self, table_id, table_name, sub_loc):
        image_opt = os.path.join(self.tmp_dir, table_id + "_" + sub_loc + '.sg_image_options.txt')
        if os.path.exists(image_opt):
            with open(image_opt, 'r') as r:
                for line in r:
                    line = line.strip().split('\t')
                    params = eval(line[1])
                    if "contribute_id" in params.keys():
                        type_id = params['type']
                        if "Heatmap" in params['method']:
                            if type_id == 1:
                                out_path = "Contribute_species_heatmap.pdf"
                            else:
                                out_path = "Contribute_function_heatmap.pdf"
                        else:
                            if type_id == 1:
                                out_path = "Contribute_species_bar.pdf"
                            else:
                                out_path = "Contribute_function_bar.pdf"
                    elif "Random" in params["method"] and "Curve" in params["method"]:
                        random_table = self.common_api.db['randomforest'].find_one({"name": table_name})
                        auth_method = json.loads(random_table["params"])["auth_method"]
                        if auth_method == "AUC":
                            out_path = "AUC_evaluation.pdf"
                        else:
                            out_path = "10fold_evaluation.pdf"
                    elif params["method"] == "showSpecimenGraphicQualityBoxplotBySpecimenName":
                        out_path = "base_quality_distribution.pdf/{}_base_quality.pdf".format(params["sample"])
                    elif params["method"] == "showSpecimenGraphicBoxBySpecimenName":
                        out_path = "base_quality_distribution.pdf/{}_base_distribution.pdf".format(params["sample"])
                    elif params["method"] == "showAnnovfdb3dPieByVfdbId":
                        sample_id = params["specimen_specimen_ids"][params["specimen_specimen_ids"].keys()[0]][
                            params["specimen_specimen_ids"][params["specimen_specimen_ids"].keys()[0]].keys()[0]][0]
                        sample_name = self.get_specimen_name(sample_id)
                        out_path = "vfdb_level2.pdf/{}_vfdb_level2_pie.pdf".format(sample_name)
                    elif params["method"] == "showAssemblestatBarByAssemId":
                        sample_id = params["specimen_name"]
                        if "_" not in sample_id:
                            sample_name = self.get_specimen_name(sample_id)
                            out_path = "contig.length.pdf/{}.contig.length.pdf".format(sample_name)
                        else:
                            out_path = "All.contig.length.pdf"
                    elif params["method"] == "showPredictgeneBarByPredictGeneId":
                        sample_id = params["specimen_name"]
                        if sample_id != "Total":
                            sample_name = self.get_specimen_name(sample_id)
                            out_path = "gene.length.pdf/{}.gene.length.pdf".format(sample_name)
                        else:
                            out_path = "Total.gene.length.pdf"
                    elif sub_loc[0:15] == "composition_bar":
                        graph_type = params["method"][15:18]
                        anno_type = sub_loc[16:]
                        if graph_type == "Bar":
                            out_path = "{}_pdf/bar.pdf".format(anno_type)
                        elif params["method"] == "showCompositionPieByCompositionId":
                            sample_id = params["specimen_id"]
                            if len(sample_id) == 24:
                                sample_name = self.get_specimen_name(sample_id)
                            else:
                                sample_name = sample_id
                            out_path = "nr_pdf/pie.pdf/{}_pie.pdf".format(sample_name)
                        else:
                            sample_id = params["specimen_id"]
                            if len(sample_id) == 24:
                                sample_name = self.get_specimen_name(sample_id)
                            else:
                                sample_name = sample_id
                            out_path = "{}_pdf/one_sample_bar.pdf/{}_bar.pdf".format(anno_type, sample_name)
                    elif params["method"] in ["showAnnogoMultiBarByGoId", "showAnnophiBarByPhiId",
                                              "showAnnomvirBarByMvirId", "showAnnotcdbBarByTcdbId"]:
                        sample_id = params["specimen_specimen_ids"][params["specimen_specimen_ids"].keys()[0]][
                            params["specimen_specimen_ids"][params["specimen_specimen_ids"].keys()[0]].keys()[0]][0]
                        sample_name = self.get_specimen_name(sample_id)
                        method_pdfdir_dic = {"showAnnogoMultiBarByGoId": "go_level_bar.pdf",
                                             "showAnnophiBarByPhiId": "phi_phenotype_bar.pdf",
                                             "showAnnomvirBarByMvirId": "Mvirdb_type_bar.pdf",
                                             "showAnnotcdbBarByTcdbId": "TCDB_class_subclass_bar.pdf"}
                        method_pdfname_dic = {"showAnnophiBarByPhiId": "phi_phenotype_bar",
                                              "showAnnomvirBarByMvirId": "Mvirdb_type_bar",
                                              "showAnnotcdbBarByTcdbId": "TCDB_class_subclass_bar"}
                        if params["method"] == "showAnnogoMultiBarByGoId":
                            level_id = params["search_level"]
                            level_dic = {"60": "2", "61": "3", "62": "4"}
                            out_path = "{}/{}_go_level{}_bar.pdf".format(method_pdfdir_dic[params["method"]],
                                                                         sample_name, level_dic[level_id])
                        elif params["method"] == "showAnnotcdbBarByTcdbId":
                            level = params["level"]
                            out_path = "{}/{}_TCDB_{}_bar.pdf".format(method_pdfdir_dic[params["method"]], sample_name,
                                                                      level)
                        else:
                            out_path = "{}/{}_{}.pdf".format(method_pdfdir_dic[params["method"]], sample_name,
                                                             method_pdfname_dic[params["method"]])
                    else:
                        out_path = self.meth2out[params["method"]].format(**params)
                    out_path = out_path.replace('__', '_')
                    out_path = os.path.join(self.output_dir, out_path)
                    print out_path
                    if not os.path.exists(os.path.dirname(out_path)):
                        os.makedirs(os.path.dirname(out_path))
                    if os.path.exists(out_path):
                        os.remove(out_path)
                    if os.path.exists(os.path.join(self.tmp_dir, line[0] + '.pdf')):
                        os.link(os.path.join(self.tmp_dir, line[0] + '.pdf'), out_path)
        self.update_status(table_id, sub_loc, {"pdf_saved": 1})

    def dir_gz(self):
        file_list = ["base_quality_distribution.pdf", "contig.length.pdf", "gene.length.pdf", "vfdb_level2.pdf", "go_level_bar.pdf", "phi_phenotype_bar.pdf","Mvirdb_type_bar.pdf", "TCDB_class_subclass_bar.pdf"]
        files = os.listdir(self.output_dir)
        os.chdir(self.output_dir)
        for file in files:
            if file in file_list:
                #file_path = self.output_dir + file
                os.system("tar -czvf {} {}".format(os.path.join(file + ".tar.gz"), file))
        dir_names = ["nr_pdf", "cog_pdf", "kegg_pdf", "cazy_pdf", "ardb_pdf", "card_pdf", "vfdb_pdf", "gene_pdf", "personal_pdf"]
        for dir in dir_names:
            if os.path.exists(os.path.join(self.output_dir, dir)):
                os.chdir(os.path.join(self.output_dir, dir))
                if os.path.exists("one_sample_bar.pdf"):
                    os.system("tar -czvf one_sample_bar.pdf.tar.gz one_sample_bar.pdf")
                elif os.path.exists("pie.pdf"):
                    os.system("tar -czvf pie.pdf.tar.gz pie.pdf")

