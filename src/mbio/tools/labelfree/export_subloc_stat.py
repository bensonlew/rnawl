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


class ExportSublocStatAgent(Agent):
    """
    需输入表达量矩阵及样本分组信息
    """
    def __init__(self, parent):
        super(ExportSublocStatAgent, self).__init__(parent)
        options = [
            dict(name="proteinsets", type="string"),
            dict(name="subloc_info", type="infile", format="labelfree.common"),
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 3
        self._memory = '15G'

    def end(self):
        super(ExportSublocStatAgent, self).end()


class ExportSublocStatTool(Tool):
    def __init__(self, config):
        super(ExportSublocStatTool, self).__init__(config)

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

    def subloc_stat(self, name2proteins):
        subloc_df = pd.read_csv(self.option('subloc_info').prop['path'], sep='\t', header=None, usecols=[0,1], names=['accession_id','subloc']).fillna('')
        subloc_dict = dict()
        for dex in subloc_df.index:
            subloc = subloc_df.loc[dex, 'subloc']
            if ':' in subloc:
                subloc = subloc.split(':')[0].strip()
            if subloc not in subloc_dict:
                subloc_dict[subloc] = dict(
                    subloc = subloc,
                    accession_id = set()
                )
            subloc_dict[subloc]['accession_id'].add(subloc_df.loc[dex, 'accession_id'])
        cols = ['subloc']
        for name in name2proteins:
            cols += ['%s_num'%name, '%s_str'%name]
        subloc_file = os.path.join(self.output_dir, 'subloc_stat.xls')
        pw = open(subloc_file, 'w')
        pw.write('\t'.join(cols)+'\n')
        for subloc, info in subloc_dict.items():
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
                if tmp['%s_str' % name]:
                    break
            else:
                write = 0
            if write:
                write_line = '\t'.join(['{' + x + '}' for x in cols])
                write_line = write_line.format(**tmp)
                pw.write(write_line + '\n')
        pw.close()

    def set_output(self):
        # all ready write results to output
        pass

    def run(self):
        super(ExportSublocStatTool, self).run()
        name2proteins = self.get_proteinset_dict()
        self.subloc_stat(name2proteins)
        self.set_output()
        self.end()



