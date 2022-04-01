# -*- coding: utf-8 -*-

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.taxon.name2taxinfo import get_taxinfo
from collections import defaultdict
import re
import shutil
import os
import json
import pandas as pd


class TaxonOutsetAgent(Agent):
    def __init__(self, parent):
        super(TaxonOutsetAgent, self).__init__(parent)
        options = [
            {"name": "result_dir", "type": "infile", "format": "meta_genomic.taxon_dir"},
            {"name": "name2id", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        pass

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(TaxonOutsetAgent, self).end()


class TaxonOutsetTool(Tool):
    """
    对metaphlan3 和 kraken2 的输出结果进行处理后用于导表
    """
    def __init__(self, config):
        super(TaxonOutsetTool, self).__init__(config)
        self.taxons = defaultdict(dict)
        self.name2id = json.loads(self.option("name2id"))
        self.output = os.path.join(self.output_dir, "taxon_")
        self.level_labs = ['D', 'K', 'P', 'C', 'O', 'F', 'G', 'S']
        self.labs_dict = {'D': "Domain", 'K': "Kingdom", 'P': "Phylum", 'C': "Class",
                          'O': "Order", 'F': "Family", 'G': "Genus", 'S': "Species"}

    def set_result(self):
        self.logger.info(self.option("result_dir").prop['level_list'])
        samples = set()
        for level_id, one_list in self.option("result_dir").prop['level_list'].items():
            if not one_list:
                continue
            one_level = None
            print(one_list)
            total = False
            for f in one_list:
                file_name = os.path.basename(f)
                sample_name = re.match(r"(.*)_\w\.xls$", file_name).groups()[0]
                samples.add(sample_name)
                one = pd.read_csv(f, sep='\t')
                if "kraken_assigned_reads" in one.columns:
                    total = True
                one = one.rename(columns={"kraken_assigned_reads": sample_name,
                                          "relative_abundance": sample_name})
                one = one[["taxonomy_id", sample_name, "linkage"]]
                if one_level is None:
                    one_level = one
                else:
                    one_level = pd.merge(one_level, one, how="outer", on=["taxonomy_id", "linkage"])
            if one_level.empty:
                continue
            one_level["level_id"] = level_id
            level_labs = [n for n in self.level_labs[:int(level_id)]]
            lab_names = []
            for lab in level_labs:
                t_lab = lab.lower() + '__'
                lab_name = self.labs_dict[lab.upper()]
                lab_names.append(lab_name)
                one_level[lab_name] = one_level["linkage"].agg(lambda x: t_lab + eval(x)[lab] if lab in eval(x) else "--")
            one_level.drop(["taxonomy_id", "linkage"], axis=1, inplace=True)
            one_level[list(samples)] = one_level[list(samples)].fillna(0)
            one_level = one_level[lab_names + ["level_id"] + list(samples)]
            if total:
                one_level["Total"] = one_level[list(samples)].sum(1)
            else:
                if "Domain" in one_level.columns:
                    one_level.drop(["Domain"], axis=1, inplace=True)
            one_level.to_csv(self.output + level_labs[int(level_id) - 1] + '.xls', sep='\t', index=False)

    def run(self):
        super(TaxonOutsetTool, self).run()
        self.set_result()
        self.end()
