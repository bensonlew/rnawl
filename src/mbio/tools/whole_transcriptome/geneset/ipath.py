# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

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
from mbio.packages.ref_rna_v2.ipath import Ipath
# from mbio.packages.lnc_rna.ipath import Ipath
from biocluster.config import Config
import subprocess

class IpathAgent(Agent):
    """
    基因集IPATH 分析
    last_modify: 2018.3.13
    """
    def __init__(self, parent):
        super(IpathAgent, self).__init__(parent)
        options = [
            {"name": "geneset_kegg", "type": "string"},  # 基因集与基因的对应文件，行数等于基因集的数目
            {"name": "kegg_table", "type": "infile", "format": "ref_rna_v2.kegg_table"},
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
        self._memory = '1G'

    def end(self):
        super(IpathAgent, self).end()


class IpathTool(Tool):
    def __init__(self, config):
        super(IpathTool, self).__init__(config)
        self.python = '/miniconda2/bin/'
        self.map_path = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/map5.r"
        self.ipath_db =  self.config.SOFTWARE_DIR + "/database/IPATH"
        # self.ipath_db =  self.config.SOFTWARE_DIR + "/database/IPATH3"
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
        Ipath1 = Ipath()
        Ipath1.set_db(self.ipath_db)
        Ipath1.set_legend(self.sets)
        Ipath1.get_K_color_width(self.ipath_input)
        Ipath1.map_file()
        self.logger.info("ipath运行完毕")
        self.set_output()
        self.end()

    def generate_ipath_file(self):
        '''
        根据基因集和注释生成ipath需求的输入文件
        '''
        gene_set_file = self.option("geneset_kegg")
        kegg_file = self.option("kegg_table").prop['path']
        # 读入基因集文件
        f = open(gene_set_file, 'r')
        lines = f.readlines()
        f.close()
        gene_set1 = []
        gene_set2 = []
        #self.logger.info("pset lines{}".format("*".join(lines)))
        gene_set1 = lines[0].strip().split("\t")[1].split(",")
        self.sets.append(lines[0].strip().split("\t")[0])
        if len(lines) == 2:
            gene_set2 = lines[1].strip().split("\t")[1].split(",")
            self.sets.append(lines[1].strip().split("\t")[0])
        gene_set1_ko = []
        gene_set2_ko = []

        ko_gene = dict()
        #self.logger.info("pset ko1{}".format(gene_set1))
        #self.logger.info("pset ko1{}".format(gene_set1_ko))
        KO_type = dict()
        # 读入每个基因集的KO
        with open(kegg_file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                acc_id = line.split('\t')[0]
                KO_ids = line.split('\t')[1].split(';')
                for ko in KO_ids:
                    if ko_gene.has_key(ko):
                        if acc_id in gene_set1 or acc_id in gene_set2:
                            ko_gene[ko].append(acc_id)
                        else:
                            pass
                    else:
                        if acc_id in gene_set1 or acc_id in gene_set2:
                            ko_gene.update({ko:[acc_id]})
                        else:
                            pass

                if acc_id in gene_set1:
                    gene_set1_ko.extend(KO_ids)
                else:
                    pass
                if acc_id in gene_set2:
                    gene_set2_ko.extend(KO_ids)
                else:
                    pass
        # self.logger.info("pset ko1{}".format(gene_set1_ko))
        # 根据交集差集些输出文件
        with open(self.ipath_input, 'wb') as w, open(self.gene_ipath_input, 'wb') as w2:
            for i in list(set(gene_set1_ko).intersection(set(gene_set2_ko))):
                w.write("{}\t#0000ff\tW15\n".format(i))
                w2.write("{}\t{}\t#0000ff\tW15\n".format(";".join(set(ko_gene[i])), i))
            for i in list(set(gene_set1_ko).difference(set(gene_set2_ko))):
                w.write("{}\t#ff0000\tW15\n".format(i))
                w2.write("{}\t{}\t#ff0000\tW15\n".format(";".join(set(ko_gene[i])), i))
            for i in list(set(gene_set2_ko).difference(set(gene_set1_ko))):
                w.write("{}\t#00ff00\tW15\n".format(i))
                w2.write("{}\t{}\t#00ff00\tW15\n".format(";".join(set(ko_gene[i])), i))
        return

    def set_output(self):
        all_files = ['gene_ipath_input.xls', 'Metabolic_pathways.svg', 'Regulatory_pathways.svg', 'Biosynthesis_of_secondary_metabolities.svg']
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir = '/mnt/ilustre/users/sanger-dev/workspace/20180312/genesetKegg_ref_rna_v2_7534_948'
        data = {
            "id": "ipath" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v2.geneset.ipath",
            "instant": True,
            "options": dict(
                geneset_kegg = test_dir + "/" + "multi_geneset_list",
                kegg_table = test_dir + "/" + "gene_kegg_table.xls"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
