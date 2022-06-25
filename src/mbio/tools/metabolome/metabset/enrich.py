# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import pandas as pd
import unittest
from mbio.packages.metabolome.common import Relation


class EnrichAgent(Agent):
    """
    代谢集富集分析
    last_modify: 2018.6.11
    """

    def __init__(self, parent):
        super(EnrichAgent, self).__init__(parent)
        options = [
            # {"name": "proteinset_kegg", "type": "string"},  # 基因集与基因的对应文件，行数等于基因集的数目
            # {"name": "kegg_table", "type": "infile", "format": "itraq_and_tmt.kegg_table"},
            {"name": "metabset", "type": "infile", "format": "metabolome.metabset"},  # 代谢集文件
            {"name": "anno_overview", "type": "infile", "format": "metabolome.overview_table"},  # 代谢集总览表
            {"name": "ko_overview", "type": "infile", "format": "sequence.profile_table"},  # 代谢集总览表ko表 add by ghd @20191015
            {"name": "correct", "type": "string", "default": "BH"},  # 多重检验校正方法
            {"name": "bg", "type": "string", "default": "project"},
            # 背景，project:本项目鉴定到的代谢物合集; species:本物种全部代谢物合集; kegg:KEGG数据库全部代谢物合集
            {"name": "species", "type": "string", "default": "all"},
            {"name": "version", "type": "string", "default": ""}, # 用于区分新老版本
        ]
        self.add_option(options)
        self.step.add_steps("enrich")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.enrich.start()
        self.step.update()

    def stepfinish(self):
        self.step.enrich.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('metabset').is_set:
            raise OptionError("必须设置代谢集", code="34701001")
        if not self.option('anno_overview').is_set:
            raise OptionError("必须设置代谢总览表", code="34701002")
        if not self.option("ko_overview").is_set:
            raise OptionError("必须设置代谢ko总览表")
        if self.option("correct") not in ["BH", "BY", "bonferroni", "holm"]:
            raise OptionError("矫正参数不在范围内，错误参数值：%s", variables=(self.option('correct')), code="34701003")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = '8G'

    def end(self):
        super(EnrichAgent, self).end()


class EnrichTool(Tool):
    def __init__(self, config):
        super(EnrichTool, self).__init__(config)
        self.python = '/miniconda2/bin/'
        self.script_path = self.config.PACKAGE_DIR + "/itraq_and_tmt/kegg_enrichment.py"
        if self.option("version") in ['v94.2', '94.2']:
            self.k2e = self.config.SOFTWARE_DIR + "/database/Annotation/all/KEGG/version_202007_meta/K2enzyme.tab"
            self.brite = self.config.SOFTWARE_DIR + "/database/Annotation/all/KEGG/version_202007_meta/br08901.txt"
            self.kegg_organisms = self.config.SOFTWARE_DIR + "/database/metabolome/kegg_v94.2_organisms.xls"
            self.kegg_com_path = self.config.SOFTWARE_DIR + "/database/metabolome/kegg_v94.2_compound_pathway.xls"
        elif self.option("version") in ['kegg', '']:
            self.k2e = self.config.SOFTWARE_DIR + "/bioinfo/rna/scripts/K2enzyme.tab"
            self.brite = self.config.SOFTWARE_DIR + "/bioinfo/rna/scripts/br08901.txt"
            self.kegg_organisms = self.config.SOFTWARE_DIR + "/database/metabolome/kegg_organisms.xls"
            self.kegg_com_path = self.config.SOFTWARE_DIR + "/database/metabolome/kegg_compound_pathway.xls"
        else:
            self.k2e = self.config.SOFTWARE_DIR + "/database/Annotation/all/KEGG/{}/K2enzyme.tab".format(self.option("version"))
            self.brite = self.config.SOFTWARE_DIR + "/database/Annotation/all/KEGG/{}/br08901.txt".format(self.option("version"))
            self.kegg_organisms = self.config.SOFTWARE_DIR + "/database/metabolome/kegg_{}_organisms.xls".format(self.option("version"))
            self.kegg_com_path = self.config.SOFTWARE_DIR + "/database/metabolome/kegg_{}_compound_pathway.xls".format(self.option("version"))
        self.compound_info = self.work_dir + "/2compound.info"
        self.pathway_info = self.work_dir + "/2path.info"
        self.compound_pathway = self.work_dir + "/comp_pathway.info"
        self.deg_path = self.work_dir + "/DE.list.check"

        if self.option('correct') == "BH":
            self.cor_code = 3
        elif self.option('correct') == "bonferroni":
            self.cor_code = 1
        elif self.option('correct') == "holm":
            self.cor_code = 2
        elif self.option('correct') == "BY":
            self.cor_code = 4
        self.bgn = set()
        self.ref_set = set()
        self.metab_set = set()
        self.global_and_overview = set(["map01100", "map01110", "map01120", "map01130", "map01200", "map01210",
                                        "map01212", "map01230", "map01220"])

    def run(self):
        """
        运行
        :return:
        """
        super(EnrichTool, self).run()
        self.read_overview(self.option('anno_overview').path, self.option("ko_overview").path)
        self.read_metaset(self.option("metabset").path)
        self.run_identify()
        self.end()

    def read_overview(self, file_path, ko_path):
        # 读总览表
        f1 = open(self.compound_info, 'w')
        f2 = open(self.pathway_info, 'w')
        f3 = open(self.compound_pathway,"w")
        f1.write("metab\tcompound\n")
        f2.write("metab\tpathway\n")
        f3.write("entry\tpathway\n")
        self.all_pathway = set()
        data = pd.read_table(file_path,sep='\t',header=0)
        ko_data = pd.read_table(ko_path, sep="\t", header=0)
        # with open(file_path, 'r') as f:
        #     lines = f.readlines()
        for i in range(len(data)):
            #for line in lines[1:]:
                #line = line.strip().split("\t")
                #if line[13] == "-":
            map_id = data['pathway_id'][i]
            metab_id = data['metab_id'][i]
            compound_id = data['compound_id'][i]
            if map_id == '-':
                continue
            # elif line[13] in ["01100", "01110", "01120"]:
            #   continue  # 01100 01110 01120 不参与富集分析 @ 20180731
            map_set = set(map_id.split(";")) - self.global_and_overview
            # self.all_pathway.update(map_set)
            self.all_pathway.update(set(map_id.split(";")))  # self.all_pathway作为ratio_in_study分母统计的范围, 需要包含global_and_overview by ghd @20191021
            map_str = ";".join(list(map_set))
            f1.write(metab_id + "\t" + compound_id + "\n")
            f2.write(metab_id + "\t" + map_str + "\n")
            self.ref_set.add(metab_id)
        for i in range(len(ko_data)): # f3 从总览表读取改为从ko总览表读取 by ghd @20191015
            compound_list = ko_data["compound_id"][i].split(";")
            self.bgn.update(compound_list)
            for compound_id in set(compound_list):
                f3.write(compound_id + "\t" + ko_data["pathway_id"][i] + "\n")
        f1.close()
        f2.close()
        f3.close()
        if self.option("bg") == "project":
            self.final_bgn = len(self.bgn)
        elif self.option("bg") == "species":
            if self.option("species") == "all" or self.option("species") == "All":
                self.run_comp_pathway()
            else:
                organism_type = self.option("species").split(";")[0]
                all_file = pd.read_table(self.kegg_organisms, sep="\t", header=0)
                if len(self.option("species").split(";")) > 1:
                    if self.option("version") not in [""]:
                        organism = self.option("species").split(";")[1].strip("_")
                    else:
                        organism = self.option("species").split(";")[1]
                    self.logger.info(organism)
                    all_file = all_file[all_file["second_category"] == organism]
                else:
                    self.logger.info(organism_type)
                    all_file = all_file[all_file["first_category"] == organism_type]
                all_file = all_file.drop("map_list", axis=1).join(
                    all_file["map_list"].str.split(';', expand=True).stack().reset_index( \
                        level=1, drop=True).rename("map_list"))
                map_list_pd = all_file["map_list"].astype('str')
                map_list_pd = pd.DataFrame(map_list_pd).drop_duplicates("map_list").reset_index()
                map_list_pd["map"] = "map"
                map_list_pd["pathway"] = map_list_pd["map"].str.cat(map_list_pd["map_list"], sep="")
                com_path_table = pd.read_table(self.kegg_com_path, sep="\t", header=0)
                com_path_table = com_path_table.drop("pathway", axis=1).join(
                    com_path_table["pathway"].str.split(';', expand=True).stack() \
                        .reset_index(level=1, drop=True).rename("pathway"))
                merge_table = pd.merge(com_path_table, map_list_pd, how='inner', on="pathway")
                merge_table = merge_table[["entry", "pathway"]]
                merge_table.to_csv(self.compound_pathway,sep="\t",index=False)
                merge_table = merge_table.drop_duplicates("entry").reset_index()
                self.final_bgn = len(merge_table)
        elif self.option("bg") == "kegg":
            self.run_comp_pathway()
        else:
            raise Exception("bg name is must be project,species,kegg")

    def run_comp_pathway(self):
        com_path_table = pd.read_table(self.kegg_com_path, sep="\t", header=0)
        com_path_table = com_path_table[com_path_table["pathway"] != "-"]
        self.final_bgn = len(com_path_table)
        ori_comp_pathway = com_path_table.drop("pathway", axis=1).join(\
            com_path_table["pathway"].str.split(';', expand=True).stack().reset_index(level=1,\
            drop=True).rename("pathway"))
        self.logger.info(self.all_pathway)
        self.logger.info(ori_comp_pathway["pathway"])
        comp_pathway = ori_comp_pathway[ori_comp_pathway["pathway"].isin(list(self.all_pathway))]
        comp_pathway.to_csv(self.compound_pathway,sep="\t",index=False)

    def read_metaset(self, file_path):
        f1 = open(self.deg_path, 'w')
        with open(file_path, 'r') as f:
            lines = f.readlines()
            for line in lines[0:]:
                line = line.strip()
                if line in self.ref_set and line != "metab_id":
                    f1.write(line + "\tnone\n")
        f1.close()

    def run_identify(self):
        # self.deg_path, self.pathway_info, self.compound_info, self.bgn, self.k2e, self.brite, self.cor_code
        cmd = self.python + 'python {} ' \
                            '-deg {} -g2p {} -g2k {} -bgn {} -k2e {} -brite {} --FDR -dn 20 -correct {}'.format(
            self.script_path, self.deg_path, self.pathway_info, self.compound_info, self.final_bgn, self.k2e,
            self.brite,
            self.cor_code)
        cmd += ' -metab {}'.format(self.compound_pathway)
        command = self.add_command("cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("kegg富集分析运行完成")
            self.set_output()
        else:
            self.set_error("kegg富集分析运行出错", code="34701001")
            raise Exception("kegg富集分析运行出错")

    def set_output(self):
        all_files = ['DE.list.check.kegg_enrichment.xls']
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            #
            table_path = self.option("anno_overview").prop["path"]
            self.logger.info(table_path)
            self.metab_trans = Relation()
            map_table, id_name_dict = self.metab_trans.get_dataframe_and_dict(table_path, metab_name="metab")
            ####self.metab_trans.add_metabolites_column(id_name_dict, each, link, id_name="Metab_ids")
            os.link(each, link)


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            #"id": "CreatTable" + str(random.randint(1, 10000)),
            "id": "Enrich",
            "type": "tool",
            "name": "metabolome.metabset.enrich",
            "instant": True,
            "options": dict(
                anno_overview="/mnt/ilustre/users/sanger-dev/workspace/20190416/MetabsetEnrich_tsg_33827_84988_847666/anno_overview_input.xls",
                metabset="/mnt/ilustre/users/sanger-dev/workspace/20190416/MetabsetEnrich_tsg_33827_84988_847666/metabset_input.set.xls",
                correct="BH",
                #bg="project",
                bg="species",
                species="All"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
