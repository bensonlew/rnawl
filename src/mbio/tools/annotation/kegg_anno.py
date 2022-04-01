# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last_modify : 2018.02.26

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
import pandas as pd
import subprocess, os,re


class KeggAnnoAgent(Agent):
    """
    调用kegg注释需要的脚本 kegg_xlm_mongo.py 获取注释信息以及对ko_name进行取名规则处理
    """

    def __init__(self, parent):
        super(KeggAnnoAgent, self).__init__(parent)
        options = [
            {"name": "kegg_xml", "type": "infile", "format": "align.blast.blast_xml"},  # 输入的比对结果xml文件
            {"name": "anno_kegg", "type": "outfile", "format": "sequence.profile_table"},  # 基因具体的注释信息
            {"name": "level_stat", "type": "outfile", "format": "sequence.profile_table"},  # levvl层级信息
            {"name": "pathway_img", "type": "string"},  # pathway通路图的dir文件路径
            {"name": "produce_mark", "type": "bool", "default":False},  # 提取出html的画图信息，供前端在底图上画图
            {"name": "description", "type": "string", "default": "false"},  #
        ]
        self.add_option(options)
        self._memory_increase_step = 10 ##add by qingchen.zhang@20201016

    def check_options(self):
        if not self.option("kegg_xml").is_set:
            raise OptionError("必须设置输入文件", code="31202101")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G' ## 改成10G，原来5G内存，由于升级过程中经常报错，20201118

    def end(self):
        super(KeggAnnoAgent, self).end()


class KeggAnnoTool(Tool):
    def __init__(self, config):
        super(KeggAnnoTool, self).__init__(config)
        self._version = "1.0"
        self.python_path = self.config.SOFTWARE_DIR + "/program/Python/bin/python"
        if self.option("description") in ['false']:## 真菌用,细菌增加了酶等描述信息，真菌未加，为区分开还是分开
            self.python_script = self.config.PACKAGE_DIR + '/fungi_genome/kegg_xlm_mongo.py'
        else:## 细菌用细菌用
            self.python_script = self.config.PACKAGE_DIR + '/bacgenome/kegg_xlm_mongo.py'
        self.python_script2 = self.config.PACKAGE_DIR + '/annotation/mg_annotation/kegg_pathway_img_v94.py'

    def run(self):
        """
        运行
        :return:
        """
        super(KeggAnnoTool, self).run()
        self.run_kegg_anno()
        self.run_level_stat()
        self.run_pathway_img()
        self.set_output()
        self.end()

    def run_kegg_anno(self):
        table = xml2table(self.option('kegg_xml').prop['path'], self.output_dir + '/tmp_kegg_table.xls')
        cmd = '{} {} {} {}'.format(self.python_path, self.python_script, table, self.output_dir)
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('运行kegg_anno完成')
        except subprocess.CalledProcessError:
            self.set_error('运行kegg_anno出错', code="31202101")

    def run_level_stat(self):
        level = self.output_dir + "/kegg_level.xls"
        kegg_level_stat = self.output_dir + "/kegg_level_stat.xls"
        table = pd.DataFrame(pd.read_table(level, sep='\t'))
        uniq_level3 = list(set(table["level3"]))
        with open(kegg_level_stat, 'w') as tab:
            tab.write("Level1\tLevel2\tLevel3\tGene nu\tGene list\tPathway list\tKO list\n")
            for l3 in uniq_level3:
                gene = table[table["level3"] == l3]["#Query"]
                pathway = table[table["level3"] == l3]["pathway_id"]
                ko = table[table["level3"] == l3]["KO"]
                gene_nu = len(gene)
                gene_list = ";".join(gene)
                pathway_list = ";".join(pathway)
                ko_list = ";".join(ko)
                l1 = list(table[table["level3"] == l3]["level1"])[0]
                l2 = list(table[table["level3"] == l3]["level2"])[0]
                level_stat = l1 + "\t" + l2 + "\t" + l3 + "\t" + str(
                    gene_nu) + "\t" + gene_list + "\t" + pathway_list + "\t" + ko_list + "\n"
                tab.write(level_stat)
        uniq_pathay = list(set(table["pathway_id"]))
        pathway_ko = self.work_dir + "/pathwaw_KO.xls"
        with open(pathway_ko, 'w') as tab:
            tab.write("pathway\tKO\n")
            for p in uniq_pathay:
                ko = table[table["pathway_id"] == p]["KO"]
                ko_list = ";".join(ko)
                if p != "-":
                    tab.write(p + "\t" + ko_list + "\n")
        os.remove(level)

    def run_pathway_img(self):
        cmd = "{}/program/Python/bin/python {} -i {} -o {} -p {} -ko {} -KO {}".format("", self.python_script2,
                                                                                       self.option("kegg_xml").prop[
                                                                                           "path"], self.output_dir,
                                                                                       self.work_dir + "/pathwaw_KO.xls",
                                                                                       "pathway", "KO")
        if self.option('produce_mark'):
            cmd += ' -png_file True -html {}/database/Annotation/all/KEGG/version_202007_meta/html/'.format(self.config.SOFTWARE_DIR)
        self.logger.info(cmd)
        command = self.add_command("output_kegg_pathway_img", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            if not self.option('produce_mark'):  # 保留原来的tar.gz 结果
                os.system('tar -zcvf kegg_pathway_img.tar.gz %s' %(self.output_dir + "/pathway_img"))
            self.logger.info("output_kegg_pathway_img succeed")
        else:
            self.set_error("output kegg pathway img failed", code="31202102")

    def set_output(self):
        self.logger.info("set_output")
        if os.path.exists(self.output_dir + "/kegg_pathway_img.tar.gz"):
            os.remove(self.output_dir + "/kegg_pathway_img.tar.gz")
        if not self.option('produce_mark'):  # 保留原来的tar.gz 结果
            os.link(self.work_dir + "/kegg_pathway_img.tar.gz",self.output_dir + "/kegg_pathway_img.tar.gz")
        files = os.listdir(self.output_dir + "/pathway_img")
        for file in files:
            if re.search(r'.png$',file):
                os.remove(self.output_dir + "/pathway_img/" + file)
        try:
            self.option("anno_kegg", self.output_dir + "/gene_kegg_anno.xls")
            self.option("level_stat", self.output_dir + "/kegg_level_stat.xls")
            self.option("pathway_img", self.output_dir + "/pathway_img")
        except Exception as e:
            self.set_error("SET_OUTFILE FAILED %s", variables=(e), code="31202103")
