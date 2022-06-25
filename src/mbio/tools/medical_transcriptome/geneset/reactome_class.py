# -*- coding: utf-8 -*-
# __author__ = 'shijin'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from collections import defaultdict
from collections import OrderedDict
import gridfs
import re
from bson import ObjectId
import time
from biocluster.config import Config
import subprocess
from biocluster.file import getsize, exists
from biocluster.file import download
from biocluster.api.file.lib.transfer import MultiFileTransfer
from mbio.packages.rna.annot_config import AnnotConfig
from mbio.packages.ref_rna_v2.copy_file import CopyFile

from mbio.packages.medical_transcriptome.geneset.reactome_svg import mutiple_svg_convert
import shutil
import unittest
from multiprocessing import Pool
from mbio.packages.ref_rna_v2.copy_file import CopyFile

# from mbio.packages.medical_transcriptome.functions import tryforgood

class ReactomeClassAgent(Agent):
    """
    Reactome分类统计分析，主要用于基因集的重运行步骤
    version v1.0.1
    author: shijin
    last_modify: 2017.8.16
    """
    def __init__(self, parent):
        super(ReactomeClassAgent, self).__init__(parent)
        options = [
            {"name": "geneset_ids", "type": "string"},  # 基因集与基因的对应文件，行数等于基因集的数目
            {"name": "reactome_annot", "type": "string"},  # 基因集与基因的对应文件，行数等于基因集的数目
            {"name": "task_id", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {'name': 'geneset_names', 'type': 'string', 'default': None},
            {"name": "background_links", "type": "string", "default": ""},  # 底图的地址信息，add_info
            {"name": "type", "type": "string", "default": "origin"},  # 取最新的注释表还是原来的注释表
            {"name": "source", "type": "string", "default": ""},
            {"name": "colors", "type": "string", "default": "#FF0000,#0000FF"},
            {'name': 'reactome_version', 'type': 'string', 'default': None}
        ]
        self.add_option(options)
        self.step.add_steps("reactome_regulate")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.reactome_regulate.start()
        self.step.update()

    def stepfinish(self):
        self.step.reactome_regulate.finish()
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
        self._cpu = 11
        self._memory = '30G'

    def end(self):
        super(ReactomeClassAgent, self).end()


class ReactomeClassTool(Tool):
    def __init__(self, config):
        super(ReactomeClassTool, self).__init__(config)
        self._version = "v1.0.1"
        self.python = '/miniconda2/bin/'
        self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.3/bin/Rscript"
        self.map_path = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/map5.r"

        self.path_name = AnnotConfig().get_file_path(
            file ="ReactomePathways.txt",
            db = "reactome",
            version = "72")

        self.path_relation = AnnotConfig().get_file_path(
            file ="ReactomePathwaysRelation.txt",
            db = "reactome",
            version = "72")
        self.svg_path = AnnotConfig().get_file_path(
            file ="svg_convert",
            db = "reactome",
            version = "72")



        self.rgene_complex_path = AnnotConfig().get_file_path(
            file ="Complex_2_hasComponent.tsv",
            db = "reactome",
            version = "72")

        self.rgene_entityset_member = AnnotConfig().get_file_path(
            file ="EntitySet_2_hasMember.tsv",
            db = "reactome",
            version = "72")

        self.rgene_candidate = AnnotConfig().get_file_path(
            file ="CandidateSet_2_hasCandidate.tsv",
            db = "reactome",
            version = "72")

        self.infered_from = AnnotConfig().get_file_path(
            file ="Event_2_inferredFrom.tsv",
            db = "reactome",
            version = "72")

        self.infered_from2 = AnnotConfig().get_file_path(
            file ="PhysicalEntity_2_inferredFrom.tsv",
            db = "reactome",
            version = "72")

        self.path_content = AnnotConfig().get_file_path(
            file ="submap.txt",
            db = "reactome",
            version = "72")

        self.colors = self.option("colors").split(",")


    def run(self):
        """
        运行
        :return:
        """
        super(ReactomeClassTool, self).run()
        self.relation_dict = self.get_path_realtion()

        # 获取节点级别关联信息
        self.complex_child = self.get_entity_member()
        self.complex_child = self.get_rgene_complex(self.complex_child)
        self.complex_child = self.get_rgene_candidate(self.complex_child)
        self.complex_child = self.get_subpath_complex(self.complex_child)
        self.complex_child = self.get_infered_from(self.complex_child)

        # with open(self.work_dir + "/complex_child", 'w') as f:
        #     f.write(self.complex_child)
        self.rpath_dict = dict()

        self.rpath_dict = dict()

        self.rpath_class = dict()
        self.get_rpath_dict()

        self.reac_genes, self.reac_path, self.gene2rgene, self.gene2rpath, self.rpath2gene = self.get_reactome_dicts()
        self.geneset2gene = self.get_geneset_dicts()
        self.get_regulate_table()
        self.get_pics()
        self.set_output()
        self.end()
        '''
        for attr in ["reac_genes", "reac_path", "gene2rgene", "gene2rpath", "rpath2gene", "geneset2gene"]:
            self.logger.info("{} is \n {}\n".format(attr, getattr(self, attr).items()[:5]))

        if not os.path.exists(pathways):
            os.mkdir(pathways)
        self.generate_ko_txt_dir()
        self.generate_new_pics()
        self.end()
        '''

    def get_entity_member(self):
        '''
        获取set的member
        '''
        complex_child = dict()
        with open(self.rgene_entityset_member, 'r') as f:
            f.readline()
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) >= 3:
                    if cols[0] in complex_child:
                        complex_child[cols[0]].add(cols[2])
                    else:
                        complex_child[cols[0]] = set([cols[2]])
        return complex_child

    def set_output(self):
        CopyFile().linkdir('svg', self.output_dir + '/svg')
        if os.path.exists(os.path.join(self.output_dir,"svg.tar.gz")):
            os.remove(os.path.join(self.output_dir,"svg.tar.gz"))
        cmd = "tar -zcvf {} {} ".format(os.path.join(self.output_dir,"svg.tar.gz"), os.path.join(self.output_dir,"svg"))
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('reactome图片打包压缩完成')
        except subprocess.CalledProcessError:
            self.logger.info('reactome图片打包压缩完成')
            self.set_error("reactome图片打包压缩失败失败")


    def get_rgene_candidate(self, complex_child):
        '''
        获取candidate
        '''
        # complex_child = dict()
        with open(self.rgene_candidate, 'r') as f:
            f.readline()
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) >= 3:
                    if cols[0] in complex_child:
                        complex_child[cols[0]].add(cols[2])
                    else:
                        complex_child[cols[0]] = set([cols[2]])
        return complex_child

    def get_rgene_complex(self, complex_child):
        '''
        获取复合节点
        '''

        # complex_child = dict()
        with open(self.rgene_complex_path, 'r') as f:
            f.readline()
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) >= 3:
                    if cols[0] in complex_child:
                        complex_child[cols[0]].add(cols[2])
                    else:
                        complex_child[cols[0]] = set([cols[2]])
        return complex_child

    def get_subpath_complex(self, complex_child):
        '''
        获取复合通路
        '''

        # complex_child = dict()
        '''
        # 只需要该通路含有的基因
        with open(self.path_relation, 'r') as f:
            f.readline()
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) >= 2:
                    par = cols[0].split("-")[-1]
                    chi = cols[1].split("-")[-1]
                    if par in complex_child:
                        complex_child[par].add(chi)
                    else:
                        complex_child[par] = set([chi])
        '''

        with open(self.option("reactome_annot"), "r") as f:
            for line in f:
                line = line.strip("\n").split("\t")
                gene = line[0]
                reactome_link = line[1]
                reactomes = reactome_link.split(";")
                for reactome in reactomes:
                    [rpath, rgene] = reactome.split("&SEL=")
                    par = rpath.split("-")[-1]
                    chi = rgene.split("-")[-1]
                    if par in complex_child:
                        complex_child[par].add(chi)
                    else:
                        complex_child[par] = set([chi])
        return complex_child

    def get_infered_from(self, complex_child):
        '''
        infered from 节点和子通路对应关系, 小鼠出现需要该文件映射id关系
        '''
        with open(self.infered_from, 'r') as f:
            f.readline()
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) >= 3:
                    if cols[0] in complex_child:
                        complex_child[cols[0]].add(cols[2])
                    else:
                        complex_child[cols[0]] = set([cols[2]])
        with open(self.infered_from2, 'r') as f:
            f.readline()
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) >= 4 and cols[3] != "'Complex'" :
                    if cols[0] in complex_child:
                        complex_child[cols[0]].add(cols[2])
                    else:
                        complex_child[cols[0]] = set([cols[2]])
        return complex_child

    def get_path_content(self, complex_child):
        '''
        子节点对应关系
        '''
        with open(self.path_content, 'r') as f:
            f.readline()
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) >= 2:
                    path = cols[0].split("-")[-1]
                    entity = cols[1].split("_")
                    if path in complex_child:
                        complex_child[path].add(entity[1])
                    else:
                        complex_child[path] = set(entity[1])
        return complex_child


    def get_rpath_dict(self):
        with open(self.path_name, 'r') as f:
            for line in f:
                cols = line.strip().split("\t")
                self.rpath_dict[cols[0]] = cols[1]
                self.rpath_class[cols[0]] = self.get_father(cols[0])

    def get_father(self, rpath):
        if rpath in self.relation_dict:
            return self.get_father(self.relation_dict[rpath])
        else:
            return rpath

    def get_path_realtion(self):
        relation_dict = dict()
        with open (self.path_relation, 'r') as f:
            for line in f:
                cols = line.strip().split("\t")
                relation_dict[cols[1]] = cols[0]
        return relation_dict

    def get_reactome_dicts(self):
        '''
        获取reactome 注释信息
        '''
        reac_genes = dict()
        reac_path = dict()
        gene2rgene = dict()
        gene2rpath = dict()
        rpath2gene = dict()

        with open(self.option("reactome_annot"), "r") as f:
            for line in f:
                line = line.strip("\n").split("\t")
                gene = line[0]
                reactome_link = line[1]
                reactomes = reactome_link.split(";")
                for reactome in reactomes:
                    [rpath, rgene] = reactome.split("&SEL=")
                    if rgene in reac_genes:
                        reac_genes[rgene].add(gene)
                    else:
                        reac_genes[rgene]= set([gene])
                    if rpath in reac_path:
                        reac_path[rpath].add(rgene)
                    else:
                        reac_path[rpath]= set([rgene])

                    if gene in gene2rgene:
                        gene2rgene[gene].add(rgene)
                    else:
                        gene2rgene[gene]= set([rgene])
                    if gene in gene2rpath:
                        gene2rpath[gene].add(rpath)
                    else:
                        gene2rpath[gene]= set([rpath])
                    if rpath in rpath2gene:
                        rpath2gene[rpath].add(gene)
                    else:
                        rpath2gene[rpath]= set([gene])

        return reac_genes, reac_path, gene2rgene, gene2rpath, rpath2gene

    def get_geneset_dicts(self):
        gene2set = dict()
        geneset2gene = OrderedDict()

        with open(self.option('geneset_ids'), 'r') as r:
            for line in r:
                cols = line.strip('\n').split('\t')
                geneset2gene[cols[0]] = cols[1].split(",")

        return geneset2gene

    def get_regulate_table(self):
        '''
        获取注释分类结果
        '''
        # 获取有哪些path
        pathids = list()
        for geneset, genes in self.geneset2gene.items():
            for g in genes:
                if g in self.gene2rpath:
                    print g, self.gene2rpath[g]
                    pathids.extend(list(self.gene2rpath[g]))
            # (pathids.extend(list(self.gene2rpath[g])) for g in genes)
        print "pathids"
        print pathids

        pathids_uniq = list(set(pathids))
        class_gene = dict()

        self.path_list = pathids_uniq
        output = self.output_dir + '/reactome_path.xls'

        with open(output, 'w') as f:
            f.write("Pathway ID\tDescription\tcategory")
            for geneset in self.geneset2gene.keys():
                f.write("\t{geneset}_numbers\t{geneset}_genes".format(geneset=geneset))
            f.write("\tlink\n")
            for rpath in pathids_uniq:
                rpath_des = self.rpath_dict[rpath]
                stat_list = list()
                rpath_class = self.rpath_class[rpath]
                rpath_class_des = self.rpath_dict[rpath_class]
                for geneset, genes in self.geneset2gene.items():
                    geneset_genes = set(genes).intersection(set(self.rpath2gene[rpath]))
                    stat_list.append(str(len(geneset_genes)))
                    stat_list.append(";".join(list(geneset_genes)))
                    if rpath_class in class_gene:
                        if geneset in class_gene[rpath_class]:
                            class_gene[rpath_class][geneset].extend(list(geneset_genes))
                        else:
                            class_gene[rpath_class][geneset] = list(geneset_genes)
                    else:
                        class_gene[rpath_class] = {geneset: list(geneset_genes)}

                link = "https://reactome.org/PathwayBrowser/#/{}".format(rpath)
                f.write("\t".join([rpath, rpath_des, rpath_class_des] + stat_list + [link]) + "\n")

        output = self.output_dir + '/reactome_class.xls'
        with open(output, 'w') as f:
            f.write("Pathway ID\tCategory Function description")
            for geneset in self.geneset2gene.keys():
                f.write("\t{geneset}_numbers\t{geneset}_genes".format(geneset=geneset))
            f.write("\n")
            for rclass in class_gene:
                class_des = self.rpath_dict[rclass]
                class_list = list()
                for geneset in self.geneset2gene.keys():
                    class_genes = list(set(class_gene[rclass].get(geneset, [])))
                    class_list.append(str(len(class_genes)))
                    class_list.append(";".join(list(class_genes)))
                f.write("\t".join([rclass, class_des] + class_list) + "\n")

    def get_pics(self):
        if not os.path.exists(self.work_dir + '/svg'):
            os.mkdir(self.work_dir + '/svg')


        paras = list()
        for path in self.path_list:
            svg_path1 = self.svg_path + '/{}.changed.svg'.format(path)
            if os.path.exists(svg_path1):
                # reac_svg = ReactomeSvg(colors = self.colors, complex_child=self.complex_child)
                path_reac_gene = dict()
                for rgene in self.reac_path[path]:
                    path_reac_gene[rgene] = self.reac_genes[rgene]

                paras.append([path, self.colors, self.complex_child, path_reac_gene.copy(), svg_path1, self.geneset2gene, self.work_dir + '/svg'])
                #reac_svg.svg_mark(svg_path1, geneset2gene=self.geneset2gene, rgene2gene=path_reac_gene, path_name=path, output_dir=self.work_dir + '/svg')
            else:
                print "can not find {} svg".format(svg_path1)
        mutiple_svg_convert(paras)

    def get_color(self, ko):
        """
        根据基因的ko号，获取颜色
        :param ko: 基因的ko号
        :return:
        """
        if len(self.category) == 1:
            if ko in self.category[self.category.keys()[0]]:
                return "red"
                # return "blue"
            else:
                return False
        elif len(self.category) == 2:
            # if self.option('source') == 'diff_exp':
            # lst = list(self.category.keys())  # 基因集列表
            lst = self.geneset_list
            # lst.sort()
            if ko in self.category[lst[0]] and ko in self.category[lst[1]]:
                return "pink"
            elif ko in self.category[lst[0]]:
                return "red"
            elif ko in self.category[lst[1]]:
                return "blue"
            else:
                return False


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "test_reactome" + str(random.randint(1, 10000))+"yyyy",
            "type": "tool",
            "name": "medical_transcriptome.geneset.reactome_class",
            "instant": False,
            "options": dict({
                'geneset_ids': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/medical_transcriptome/ref_annot_class2/reactome/ids.list',
                'geneset_names': 'set1,set2',
                "reactome_annot": "/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/medical_transcriptome/ref_annot_class2/reactome/reactome.list"
            })
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
