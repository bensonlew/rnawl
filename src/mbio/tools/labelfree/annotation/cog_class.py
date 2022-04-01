# -*- coding: utf-8 -*-
# __author__ = 'shijin'

from __future__ import division
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from collections import defaultdict
import os

class CogClassAgent(Agent):
    """
    cog功能分类
    """
    def __init__(self, parent):
        super(CogClassAgent, self).__init__(parent)
        options = [
           {"name": "diff_list", "type": "infile", "format": "rna.gene_list"},
            {"name": "cog_table", "type": "infile", "format": "annotation.cog.cog_table"}
        ]
        self.add_option(options)
        self.step.add_steps("cog_class")
        self.on('start', self.start_cog_class)
        self.on("end", self.end_cog_class)

    def start_cog_class(self):
        self.step.cog_class.start()
        self.step.update()

    def end_cog_class(self):
        self.step.cog_class.finish()
        self.step.update()

    def check_options(self):
        if not self.option("diff_list").is_set:
            raise OptionError("参数diff_list不能为空")

    def set_resource(self):
        self._cpu = 4
        self._memory = "4G"

class CogClassTool(Tool):
    def __init__(self, config):
        super(CogClassTool, self).__init__(config)
        self.gene_list = self.option("diff_list").prop["gene_list"]
        self.query_list = []
        self.group_dict = defaultdict(lambda: [])
        self.categories_dict = defaultdict(lambda:{})
        self.func_type = {
            'INFORMATION STORAGE AND PROCESSING': ['A', 'B', 'J', 'K', 'L'],
            'CELLULAR PROCESSES AND SIGNALING': ['D', 'M', 'N', 'O', 'T', 'U', 'V', 'W', 'Y', 'Z'],
            'METABOLISM': ['C', 'E', 'F', 'G', 'H', 'I', 'P', 'Q'],
            'POORLY CHARACTERIZED': ['R', 'S'],

        }
        self.decs_class = {
            'A': 'INFORMATION STORAGE AND PROCESSING',
            'B': 'INFORMATION STORAGE AND PROCESSING',
            'C': 'METABOLISM',
            'D': 'CELLULAR PROCESSES AND SIGNALING',
            'E': 'METABOLISM',
            'F': 'METABOLISM',
            'G': 'METABOLISM',
            'H': 'METABOLISM',
            'I': 'METABOLISM',
            'J': 'INFORMATION STORAGE AND PROCESSING',
            'K': 'INFORMATION STORAGE AND PROCESSING',
            'L': 'INFORMATION STORAGE AND PROCESSING',
            'M': 'CELLULAR PROCESSES AND SIGNALING',
            'N': 'CELLULAR PROCESSES AND SIGNALING',
            'O': 'CELLULAR PROCESSES AND SIGNALING',
            'P': 'METABOLISM',
            'Q': 'METABOLISM',
            'R': 'POORLY CHARACTERIZED',
            'S': 'POORLY CHARACTERIZED',
            'T': 'CELLULAR PROCESSES AND SIGNALING',
            'U': 'CELLULAR PROCESSES AND SIGNALING',
            'V': 'CELLULAR PROCESSES AND SIGNALING',
            'W': 'CELLULAR PROCESSES AND SIGNALING',
            'Y': 'CELLULAR PROCESSES AND SIGNALING',
            'Z': 'CELLULAR PROCESSES AND SIGNALING'
        }
        self.func_decs = {
            'A': 'RNA processing and modification',
            'B': 'Chromatin structure and dynamics',
            'C': 'Energy production and conversion',
            'D': 'Cell cycle control, cell division, chromosome partitioning',
            'E': 'Amino acid transport and metabolism',
            'F': 'Nucleotide transport and metabolism',
            'G': 'Carbohydrate transport and metabolism',
            'H': 'Coenzyme transport and metabolism',
            'I': 'Lipid transport and metabolism',
            'J': 'Translation, ribosomal structure and biogenesis',
            'K': 'Transcription',
            'L': 'Replication, recombination and repair',
            'M': 'Cell wall/membrane/envelope biogenesis',
            'N': 'Cell motility',
            'O': 'Posttranslational modification, protein turnover, chaperones',
            'P': 'Inorganic ion transport and metabolism',
            'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
            'R': 'General function prediction only',
            'S': 'Function unknown',
            'T': 'Signal transduction mechanisms',
            'U': 'Intracellular trafficking, secretion, and vesicular transport',
            'V': 'Defense mechanisms',
            'W': 'Extracellular structures',
            'Y': 'Nuclear structure',
            'Z': 'Cytoskeleton'
        }
        self.summary_file = open(self.output_dir + "/cog_summary.xls", "w")

    def cmd1(self):
        with open(self.option("cog_table").prop["path"], "r") as file:
            for line in file:
                line = line.strip()
                if line.startswith("#"):
                    continue
                tmp = line.split("\t")
                query = tmp[0]
                cog_group = tmp[10]
                cog_categories = tmp[12]
                if query in self.gene_list:
                    # if query not in self.query_list:
                        # self.query_list.append(query)
                    for letter in ["C", "N", "K"]:
                        if not self.categories_dict[cog_categories].has_key(letter):
                            self.categories_dict[cog_categories][letter] = []
                        if cog_group.startswith(letter):
                            if query not in self.categories_dict[cog_categories][letter]:
                                self.categories_dict[cog_categories][letter].append(query)
        name = os.path.splitext(os.path.basename(self.option("diff_list").prop["path"]))[0]
        head = "type\tcategory\t{}_COG\t{}_NOG\t{}_KOG\t{}_COG_LIST\t{}_NOG_LIST\t{}_KOG_LIST\n".format(name, name, name, name, name, name)
        self.summary_file.write(head)
        for g in sorted(self.func_decs.keys()):
            detail = self.func_decs[g]
            category = '[' + g + ']' + ' ' + detail
            try:
                coglist = self.categories_dict[g]["C"]
                cogcount = len(coglist)
            except KeyError:
                cogcount = 0
                coglist = []
            try:
                noglist = self.categories_dict[g]["N"]
                nogcount = len(noglist)
            except KeyError:
                nogcount = 0
                noglist = []
            try:
                koglist = self.categories_dict[g]["K"]
                kogcount = len(koglist)
            except KeyError:
                kogcount = 0
                koglist = []
            self.summary_file.write(  self.decs_class[g]  + '\t' + category + '\t' + str(cogcount) + '\t' + str(nogcount) + '\t'
                                        + str(kogcount) + '\t' + ';'.join(coglist) + '\t' + ';'.join(noglist) + '\t'
                                        + ';'.join(koglist) + '\n')
        self.summary_file.close()


        '''按一级分类排序
        for thekey in ['INFORMATION STORAGE AND PROCESSING', 'CELLULAR PROCESSES AND SIGNALING',
                       'METABOLISM', 'POORLY CHARACTERIZED']:
            for g in self.func_type[thekey]:
                detail = self.func_decs[g]
                category = '[' + g + ']' + ' ' + detail
                try:
                    coglist = self.categories_dict[g]["C"]
                    cogcount = len(coglist)
                except KeyError:
                    cogcount = 0
                    coglist = []
                try:
                    noglist = self.categories_dict[g]["N"]
                    nogcount = len(noglist)
                except KeyError:
                    nogcount = 0
                    noglist = []
                try:
                    koglist = self.categories_dict[g]["K"]
                    kogcount = len(koglist)
                except KeyError:
                    kogcount = 0
                    koglist = []
                self.summary_file.write(thekey + '\t' + category + '\t' + str(cogcount) + '\t' + str(nogcount) + '\t'
                                        + str(kogcount) + '\t' + ';'.join(coglist) + '\t' + ';'.join(noglist) + '\t'
                                        + ';'.join(koglist) + '\n')
        self.summary_file.close()
        '''
    def run(self):
        super(CogClassTool, self).run()
        self.cmd1()
        self.end()
