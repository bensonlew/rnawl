# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import re
from collections import defaultdict
import gridfs
from bson import ObjectId
import time
import json
import unittest
from mbio.packages.metabolome.ipath import Ipath
from biocluster.config import Config
import subprocess


class IpathAgent(Agent):
    """
    代谢集IPATH 分析
    last_modify: 2018.6.5
    """

    def __init__(self, parent):
        super(IpathAgent, self).__init__(parent)
        options = [
            # {"name": "proteinset_kegg", "type": "string"},  # 基因集与基因的对应文件，行数等于基因集的数目
            # {"name": "kegg_table", "type": "infile", "format": "itraq_and_tmt.kegg_table"},
            {"name": "metabset", "type": "infile", "format": "metabolome.mul_metabset"},  # 代谢集文件
            {"name": "anno_overview", "type": "infile", "format": "metabolome.overview_table"},  # 代谢集总览表
            {"name": "main_table_id", "type": "string"}
        ]
        self.add_option(options)
        self.step.add_steps("ipath")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.ipath.start()
        self.step.update()

    def stepfinish(self):
        self.step.ipath.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(IpathAgent, self).end()


class IpathTool(Tool):
    def __init__(self, config):
        super(IpathTool, self).__init__(config)
        # self.python = '/program/Python/bin/'
        # self.map_path = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/map5.r"
        self.ipath_db = self.config.SOFTWARE_DIR + "/database/IPATH"
        self.image_magick = self.config.SOFTWARE_DIR + "/program/ImageMagick/bin/convert"
        self.ipath_input = self.work_dir + '/ipath_input.xls'
        self.gene_ipath_input = self.work_dir + '/gene_ipath_input.xls'
        self.sets = []

    def run(self):
        """
        运行
        :return:
        """
        super(IpathTool, self).run()
        # self.get_kegg_pics()
        self.generate_ipath_file()
        try:
            Ipath1 = Ipath()
            Ipath1.set_db(self.ipath_db)
            Ipath1.set_legend(self.sets)
            Ipath1.get_K_color_width(self.ipath_input)
            Ipath1.map_file(self.image_magick)
        except Exception,e:
            self.logger.info(e)
            self.set_error("ipath运行错误: %s", variables=(e), code="34701201")
        self.logger.info("ipath运行完毕")
        self.set_output()
        self.end()

    def generate_ipath_file(self):
        '''
        根据基因集和注释生成ipath需求的输入文件
        '''
        metab_set_file = self.option("metabset").path
        anno_overview_file = self.option("anno_overview").path
        metab_set1, metab_set2 = self.read_set(metab_set_file)
        bingji_set = metab_set1 | metab_set2
        # self.logger.info("DEBUG:>>>")
        # self.logger.info(len(metab_set1))
        # self.logger.info(len(metab_set2))
        # self.logger.info(len(bingji_set))
        # self.logger.info(len(jiaoji_set))
        jiaoji_set = metab_set1 & metab_set2
        metab_name_dic, compound_id_dic = self.read_overview(anno_overview_file, bingji_set)
        with open(self.ipath_input, 'wb') as w, open(self.gene_ipath_input, 'wb') as w2:
            for one_metab in jiaoji_set:
                if not compound_id_dic.has_key(one_metab) or compound_id_dic[one_metab] == "-":
                    continue
                # self.logger.info("DEBUG1:%s::%s::%s::" % (one_metab, metab_name_dic[one_metab], compound_id_dic[one_metab]))
                w.write("{}\t#0000ff\tW15\n".format(compound_id_dic[one_metab]))
                w2.write("{}\t{}\t#0000ff\tW15\n".format(metab_name_dic[one_metab], compound_id_dic[one_metab]))
            for one_metab in (metab_set1 - metab_set2):
                if not compound_id_dic.has_key(one_metab) or compound_id_dic[one_metab] == "-":
                    continue
                # self.logger.info("DEBUG2:%s::%s::%s::" % (one_metab, metab_name_dic[one_metab], compound_id_dic[one_metab]))
                w.write("{}\t#ff0000\tW15\n".format(compound_id_dic[one_metab]))
                w2.write("{}\t{}\t#ff0000\tW15\n".format(metab_name_dic[one_metab], compound_id_dic[one_metab]))
            for one_metab in (metab_set2 - metab_set1):
                if not compound_id_dic.has_key(one_metab) or compound_id_dic[one_metab] == "-":
                    continue
                # self.logger.info("DEBUG3:%s::%s::%s::" % (one_metab, metab_name_dic[one_metab], compound_id_dic[one_metab]))
                w.write("{}\t#00ff00\tW15\n".format(compound_id_dic[one_metab]))
                w2.write("{}\t{}\t#00ff00\tW15\n".format(metab_name_dic[one_metab], compound_id_dic[one_metab]))

    def read_set(self, file_path):
        # 读代谢集，格式同蛋白组，#代谢集名称    代谢集列表(,分割)
        f = open(file_path, 'r')
        lines = f.readlines()
        f.close()
        metab_set1 = lines[0].strip().split("\t")[1].split(",")
        metab_set2 = []
        self.sets.append(lines[0].strip().split("\t")[0])
        if len(lines) == 2:
            metab_set2 = lines[1].strip().split("\t")[1].split(",")
            self.sets.append(lines[1].strip().split("\t")[0])
        return set(metab_set1), set(metab_set2)

    def read_overview(self, file_path, metab_set):
        # 读总览表
        id2name = {}
        id2compound = {}
        with open(file_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                metab_id = line.split('\t')[0]
                metab_name = line.split('\t')[1]
                compound_id = line.split('\t')[9]
                if metab_id in metab_set:
                    id2name[metab_id] = metab_name
                    id2compound[metab_id] = compound_id
        return id2name, id2compound  # 可能compound里面有"-"

    def set_output(self):
        all_files = ['gene_ipath_input.xls', 'Metabolic_pathways.svg', 'Regulatory_pathways.svg',
        # all_files = ['Metabolic_pathways.svg', 'Regulatory_pathways.svg',
                     'Biosynthesis_of_secondary_metabolities.svg', 'Metabolic_pathways.png',\
                     'Regulatory_pathways.png', 'Biosynthesis_of_secondary_metabolities.png']
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)


"""
class TestFunction(unittest.TestCase):
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir = '/mnt/ilustre/users/sanger-dev/workspace/20180312/ProteinsetKegg_itraq_tmt_7534_948'
        data = {
            "id": "ipath" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "itraq_and_tmt.proteinset.ipath",
            "instant": True,
            "options": dict(
                proteinset_kegg = test_dir + "/" + "multi_proteinset_list",
                kegg_table = test_dir + "/" + "protein_kegg_table.xls"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
"""
