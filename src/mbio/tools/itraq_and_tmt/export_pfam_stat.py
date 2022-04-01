# coding=utf-8
import os
# import unittest
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from decimal import *
from collections import OrderedDict
import copy
__author__ = 'fengyitong'


class ExportPfamStatAgent(Agent):
    """
    需输入表达量矩阵及样本分组信息
    """
    def __init__(self, parent):
        super(ExportPfamStatAgent, self).__init__(parent)
        options = [
            dict(name="proteinsets", type="string"),
            dict(name="pfam_info", type="infile", format="itraq_and_tmt.common"),
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(ExportPfamStatAgent, self).end()


class ExportPfamStatTool(Tool):
    def __init__(self, config):
        super(ExportPfamStatTool, self).__init__(config)

    def get_proteinset_dict(self):
        name2proteins = OrderedDict()
        for proteinset in self.option('proteinsets').split(','):
            name = os.path.basename(proteinset).split('_protein.list')[0]
            with open(proteinset) as pr:
                proteins = pr.read().split('\n')
                while '' in proteins:
                    proteins.remove('')
            name2proteins[name] = proteins
        return name2proteins

    def pfam_stat(self, name2proteins):
        pfam_df = pd.read_csv(self.option('pfam_info').prop['path'], sep='\t').fillna('')
        pfam_dict = dict()
        for dex in pfam_df.index:
            pfam = pfam_df.loc[dex, 'pfam']
            if pfam not in pfam_dict:
                pfam_dict[pfam] = dict(
                    pfam = pfam,
                    domain = pfam_df.loc[dex, 'domain'],
                    description = pfam_df.loc[dex, 'description'],
                    accession_id = set()
                )
            pfam_dict[pfam]['accession_id'].add(pfam_df.loc[dex, 'accession_id'])
        cols = ['pfam', 'domain', 'description']
        for name in name2proteins:
            cols += ['%s_num'%name, '%s_str'%name]
        pfam_file = os.path.join(self.output_dir, 'pfam_stat.xls')
        pw = open(pfam_file, 'w')
        pw.write('\t'.join(cols)+'\n')
        for pfam, info in pfam_dict.items():
            tmp = copy.copy(info)
            for name, proteins in name2proteins.items():
                tmp[name] = set(proteins) & tmp['accession_id']
                if not tmp[name]:
                    tmp['%s_num'%name] = '0'
                    tmp['%s_str'%name] = ''
                else:
                    tmp['%s_num' % name] = str(len(tmp[name]))
                    tmp['%s_str' % name] = ';'.join(tmp[name])
            # 把都是零的过滤掉
            write = 1
            for name in name2proteins:
                if tmp['%s_str'%name]:
                    break
            else:
                write = 0
            if write:
                write_line = '\t'.join(['{'+x+'}' for x in cols])
                write_line = write_line.format(**tmp)
                pw.write(write_line + '\n')
        pw.close()

    def set_output(self):
        # all ready write results to output
        pass

    def run(self):
        super(ExportPfamStatTool, self).run()
        name2proteins = self.get_proteinset_dict()
        self.pfam_stat(name2proteins)
        self.set_output()
        self.end()



