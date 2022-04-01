# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import re
import unittest


class EnterzToKeggAgent(Agent):
    """
    根据ncbi基因组提取信息  enterz, gene name, gene description
    version v1.0.1
    author: liubinxu
    last_modify: 2019.01.14
    filter_species: Protists|Plants|Fungi|Bacteria|Animals|Archaea
    enterz: ENSG00000284241\tENST00000639108\t3804
    """
    def __init__(self, parent):
        super(EnterzToKeggAgent, self).__init__(parent)
        options = [
            {'name': 'enterz', 'type': 'string', 'default': None},
            {'name': 'output', 'type': 'string', 'default': 'pathway'},
            {'name': 'filter_species', 'type': 'string', 'default': None}
        ]
        self.add_option(options)
        self.step.add_steps("enterz")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.enterz.start()
        self.step.update()

    def stepfinish(self):
        self.step.enterz.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('enterz'):
            raise OptionError('必须输入enterz')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = "30G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(EnterzToKeggAgent, self).end()


class EnterzToKeggTool(Tool):
    def __init__(self, config):
        super(EnterzToKeggTool, self).__init__(config)
        self._version = "v1.0.1"
        self.ncbi_gff = self.config.PACKAGE_DIR + "/ref_genome_db/ncbi_gff.py"
        self.kegg_genes = self.config.SOFTWARE_DIR + "/database/KEGG/ko_genes.list"
        self.path_ko = self.config.SOFTWARE_DIR + "/database/KEGG/pathway"

        self.taxonomy_path = self.config.SOFTWARE_DIR + "/database/KEGG/species/{}.ko.txt".format(self.option("filter_species"))


    def run(self):
        """
        运行
        :return:
        """
        super(EnterzToKeggTool, self).run()
        tran2enterz = self.run_get_enterz()
        ko2path = self.get_ko2path()
        self.logger.info("ko2path {}".format(ko2path.items()[:5]))
        self.get_ko(tran2enterz, ko2path)
        self.set_output()
        self.end()

    def get_ko2path(self):
        '''
        提取ko2path信息
        '''
        filter_path = set()
        if self.option("filter_species"):
            with open(self.taxonomy_path, 'r') as f:
                for line in f:
                    filter_path.add(line.strip().lstrip('ko'))

        ko2path = dict()
        with open(self.path_ko, 'r') as f:
            for line in f:
                cols = line.strip().split("\t")
                if cols[0].lstrip("path:map") in filter_path or len(filter_path) == 0:
                    if cols[1].lstrip("ko:") in ko2path:
                        ko2path[cols[1].lstrip("ko:")].append(cols[0].lstrip("path:map"))
                    else:
                        ko2path[cols[1].lstrip("ko:")] = [cols[0].lstrip("path:map")]
        return ko2path


    def run_get_enterz(self):
        '''
        提取enterz 信息
        '''
        trans_dict = dict()
        with open(self.option("enterz"), 'r') as f:
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) == 3:
                    trans_dict[cols[1]] = [cols[0], cols[2]]
                elif len(cols) == 2:
                    trans_dict[cols[1]] = [cols[0], '']
        return trans_dict

    def get_enterz_ko(self, enterz_set):
        enterz2ko = dict()
        with open(self.kegg_genes, 'r') as f:
            for line in f:
                cols = line.strip().split("\t")
                if cols[1].split(":")[1] in enterz_set:
                    enterz2ko[cols[1].split(":")[1]] = cols[0].split(":")[1]
        self.logger.info("enterz2ko {}".format(enterz2ko.items()[:5]))
        return enterz2ko

    def get_ko(self, tran2enterz, ko2path):
        '''
        提取enterz ko信息
        '''
        enterzs = [tran2enterz[x][1] for x in tran2enterz if tran2enterz[x][1] != '']
        enterz_set = set(enterzs)
        self.logger.info("enterz_set {}".format(list(enterz_set)[:5]))
        enterz2ko = self.get_enterz_ko(enterz_set)
        with open(os.path.join(self.work_dir, self.option("output")), 'w') as pathway_f:
            for tran in tran2enterz:
                gene = tran2enterz[tran][0]
                enterz = tran2enterz[tran][1]
                if enterz != "" and enterz in enterz2ko:
                    ko = enterz2ko[enterz]
                else:
                    ko = ""
                if ko != "" and ko in ko2path:
                    pathways = ";".join(ko2path[ko])
                else:
                    pathways = ""
                if ko != "":
                    pathway_f.write("\t".join([gene, tran, enterz, ko, pathways]) + "\n")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            all_files = [self.option("output")]
            for each in all_files:
                fname = os.path.basename(each)
                link = os.path.join(self.output_dir, fname)
                if os.path.exists(link):
                    os.remove(link)
                os.link(each, link)
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet

        data = {
            "id": "enterz_to_kegg" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_genome_db.enterz_to_kegg",
            "instant": True,
            "options": dict(
                enterz = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/NCBI/Homo_sapiens.GRCh38.biomart_enterz.txt",
                filter_species = "Animals"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
