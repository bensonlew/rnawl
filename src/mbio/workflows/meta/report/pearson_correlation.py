# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
from biocluster.workflow import Workflow
import glob
import os
from bson import ObjectId
import re
from mbio.packages.meta.common_function import env_link
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class PearsonCorrelationWorkflow(Workflow):
    """
    报告中计算alpha多样性指数时使用
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(PearsonCorrelationWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", 'format': "meta.otu.otu_table"},  # 输入的OTU id
            {"name": "otu_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "env_file", "type": "infile", 'format': "meta.otu.group_table"},  # 输入的OTU id
            {"name": "env_id", "type": "string"},
            {"name": "corr_id", "type": "string"},
            {"name": "env_labs", "type": "string"},
            {"name": "level", "type": "int"},
            {"name": "correlation_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "method", "type": "string", "default": "pearsonr"},
            {"name": "env_cluster", "type": "string", "default": "average"},
            {"name": "species_cluster", "type": "string", "default": "average"},
            {"name": "group_detail", "type": "string"},
            {"name": "top_species", "type": "int", "default": 0}  # add new option (flit top N species)
            ]
        self.add_option(options)
        # print(self._sheet.options())
        self.set_options(self._sheet.options())
        self.correlation = self.add_tool('statistical.pearsons_correlation')
        self.params = {}
        self.name_to_name = {}
        self.env_name = {}

    def run_correlation(self):
        env_cluster = self.option("env_cluster")
        species_cluster = self.option("species_cluster")
        if env_cluster == "":
            env_cluster = "average"
        if species_cluster == "":
            species_cluster = "average"
        level = self.convert_level()
        options = {
            'otutable': self.option('otu_file'),
            'envtable': self.option('env_file'),
            "method": self.option('method'),
            "env_cluster": env_cluster,
            "species_cluster": species_cluster,
            "top_species": self.option('top_species'),
            "project":"meta"
            }
        if level in ['otu', 'domain', 'kindom', 'phylum', 'class', 'order','family', 'genus', 'species']:
            options['level'] = level
        self.correlation.set_options(options)
        self.correlation.on("end", self.set_db)
        self.correlation.run()

    def convert_level(self):
        """
        将level由int转为string
        """
        self.logger.info("开始转换level的类型")
        level_dict = {
            1 : "domain",
            2 : "kindom",
            3 : "phylum",
            4 : "class",
            5 : "order",
            6 : "family",
            7 : "genus",
            8 : "species",
            9 : "otu",
        }
        if self.option("level") in level_dict.keys():
            level_string = level_dict[self.option("level")]
            return level_string

    def run(self):
        self.run_correlation()
        super(PearsonCorrelationWorkflow, self).run()

    def get_name(self):
        with open(self.correlation.work_dir + "/name_to_name.xls", "r") as f:
            for line in f:
                line = line.strip().split("\t")
                self.name_to_name[line[0]] = line[1]
        with open(self.correlation.work_dir + "/env_name.xls", "r") as ef:
            self.env_name = eval(ef.readline())
            self.logger.info(self.env_name)

    def dashrepl(self, matchobj):
        return self.name_to_name[matchobj.groups()[0]]

    def dashrepl_env(self, matchobj):
        return self.env_name[matchobj.groups()[0]]

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        new_species_tree = ""
        env_tree = ""
        new_env_tree = ""
        env_list = []
        species_list = []
        api_correlation = self.api.meta_species_env
        corr_path = glob.glob(self.correlation.output_dir+"/*correlation*")
        pvalue_path = glob.glob(self.correlation.output_dir+"/*pvalue*")

        env_tree_path = self.correlation.work_dir + "/env_tree.tre"
        species_tree_path = self.correlation.work_dir + "/species_tree.tre"

        self.get_name()

        if os.path.exists(env_tree_path):
            with open(env_tree_path, "r") as f:
                env_tree = f.readline().strip()
                raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', env_tree)
                env_list = [self.env_name[i[1]] for i in raw_samp]
                new_env_tree = re.sub(r"(colnew\d+)", self.dashrepl_env, env_tree)
                # print(env_list)
        if os.path.exists(species_tree_path):
            with open(species_tree_path, "r") as f:
                species_tree = f.readline().strip()
                raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', species_tree)
                # print(self.name_to_name)
                species_list = [self.name_to_name[i[1]] for i in raw_samp]

                new_species_tree = re.sub(r"(name\d+)", self.dashrepl, species_tree)
                print(new_species_tree)
                # print(species_list)
                # new_species_list = []

        corr_id = ObjectId(self.option("corr_id"))
        api_correlation.add_correlation_detail(corr_path[0], "correlation", corr_id)
        api_correlation.add_correlation_detail(pvalue_path[0], "pvalue", corr_id, species_tree=new_species_tree,
                                               env_tree=new_env_tree, env_list=env_list, species_list=species_list)
        # self.add_return_mongo_id('sg_species_env_correlation', corr_id)
        #self.end()
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("corr_id"), "sg_species_env_correlation")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("corr_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "beta_pearson_correlation",
                "interaction": 1,
                "main_table": "sg_species_env_correlation",
            })
            self.figsave.run()
        else:
            self.end()
    
    def envname_restore(self, output_dir):
        if self.option('env_file').is_set:
            env_id = self.option('env_id')
            self.logger.info('修改结果文件中环境因子名称为真实名称')
            outpath = os.path.join(self.work_dir, 'output_upload')
            env_link(self, env_id, output_dir, outpath)
            self.logger.info('完成名称修改')
            self.output_dir = outpath

    def end(self):
        self.envname_restore(self.correlation.output_dir)
        if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "相关性Heatmap分析结果目录", 0, "110154"],   # add 2 lines by hongdongxuan 20170324
            ["./pearsons_correlation_at_otu_level.xls", "xls", "相关性系数表", 0, "110155"],
            ["./pearsons_pvalue_at_otu_level.xls", "xls", "相关性P值", 0, "110156"],
            ["./相关性Heatmap图.pdf", "pdf", "物种与环境因子相关性Heatmap图", 0, ""]
            # ["./mantel_results.txt", "txt", "mantel检验结果"]
        ])
        # print self.get_upload_files()
        super(PearsonCorrelationWorkflow, self).end()
