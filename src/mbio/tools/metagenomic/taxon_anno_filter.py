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
from functools import reduce


class TaxonAnnoFilterAgent(Agent):
    def __init__(self, parent):
        super(TaxonAnnoFilterAgent, self).__init__(parent)
        options = [
            {"name": "anno_file", "type": "infile", "format": "meta_genomic.taxon_dir"},
            {"name": "samples", "type": "string"},
            {"name": "level_filter", "type": "string", "default": ""},
            {"name": "sample_filter", "type": "string", "default": ""},
            {"name": "abu_filter", "type": "string", "default": ""},
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
        super(TaxonAnnoFilterAgent, self).end()


class TaxonAnnoFilterTool(Tool):
    """
    用于筛选metaphlan 和 kraken的结果
    """
    def __init__(self, config):
        super(TaxonAnnoFilterTool, self).__init__(config)
        self.taxons = defaultdict(dict)
        self.output = os.path.join(self.output_dir, "taxon_")
        self.all_level = self.option("anno_file").prop['level_list']
        self.level_labs = ['D', 'K', 'P', 'C', 'O', 'F', 'G', 'S']
        self.level_col = ['d__', 'k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
        self.level_col_name = ["Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
        self.default_col = ["level_id", ] + self.level_col
        self.total = False

    def get_lowest_level(self):
        lowest_level = max([int(l) for l in self.all_level])
        lowest_df = []
        for f in self.all_level[str(lowest_level)]:
            one = pd.read_csv(f, sep='\t')
            lowest_df.append(one)
        lowest_df = pd.concat(lowest_df, ignore_index=True)
        if "Total" in lowest_df.columns:
            self.total = True
        lowest_df = lowest_df.rename(columns={self.level_col_name[i]: self.level_col[i] for i in range(8)})
        samples = json.loads(self.option("samples"))
        select_cols = [c for c in lowest_df.columns if c in self.default_col] + samples
        lowest_df = lowest_df[select_cols]
        return lowest_df, lowest_level

    def level_filter(self, df):
        level_filter = json.loads(self.option("level_filter"))
        filter_s = None
        for filt in level_filter:  # {type: 1, level: 1, name: d__Bac}
            level_col = self.level_col[int(filt["level"]) - 1]
            if int(filt["type"]) == 1:
                ft = df[level_col] == filt["name"]
            else:
                ft = df[level_col] != filt["name"]
            if filter_s is None:
                filter_s = ft
            else:
                filter_s = filter_s | ft
        return filter_s

    def sample_filter(self, df):
        sample_filter = json.loads(self.option("sample_filter"))
        filter_s = None
        samples = json.loads(self.option("samples"))
        sp_df = df[samples]
        sp_df = sp_df[samples] / sp_df.sum()
        for filt in sample_filter:  # {type: 1, sample: 1, abu:0.1}
            ft = (sp_df > float(filt["abu"])).astype(int).sum(1) >= int(filt["sample"])
            if int(filt["type"]) == 0:
                ft ^= True
            if filter_s is None:
                filter_s = ft
            else:
                filter_s = filter_s | ft
        return filter_s

    def abu_filter(self, df):
        abu_filter = json.loads(self.option("abu_filter"))
        filter_s = None
        if "Total" in df.columns:
            total = df["Total"]
        else:
            samples = json.loads(self.option("samples"))
            total = df[samples].sum(1)
        total_r = total / total.sum()
        for filt in abu_filter:  # {t_abu: 0.3, type: 1}
            ft = total_r >= float(filt["t_abu"])
            if int(filt["type"]) == 0:
                ft ^= True
            if filter_s is None:
                filter_s = ft
            else:
                filter_s = filter_s | ft
        return filter_s

    def out_level(self, df, level, lowest=False):
        if len(self.all_level[str(level + 1)]) == 0:
            return
        ori_file = self.all_level[str(level + 1)][0]
        outfile = os.path.join(self.output_dir, os.path.basename(ori_file))
        samples = json.loads(self.option("samples"))
        if lowest:
            out_df = df
        else:
            out_df = pd.read_csv(ori_file, sep='\t')
            out_df = out_df.rename(columns={self.level_col_name[i]: self.level_col[i] for i in range(8)})
            select_cols = [c for c in out_df.columns if c in self.default_col] + samples
            out_df = out_df[select_cols]
            level_col = self.level_col[level]
            out_df = out_df[out_df[level_col].isin(df[level_col])]
        if not out_df.empty:
            out_df = out_df.rename(columns={self.level_col[i]: self.level_col_name[i] for i in range(8)})
            if self.total:
                out_df['Total'] = out_df[samples].sum(1)
            out_df.to_csv(outfile, sep='\t', index=False)

    def run(self):
        super(TaxonAnnoFilterTool, self).run()
        lowest_df, lowest_level = self.get_lowest_level()
        all_f = []
        if self.option("level_filter"):
            lvl_f = self.level_filter(lowest_df)
            all_f.append(lvl_f)
        if self.option("sample_filter"):
            sp_f = self.sample_filter(lowest_df)
            all_f.append(sp_f)
        if self.option("abu_filter"):
            ab_f = self.abu_filter(lowest_df)
            all_f.append(ab_f)
        lowest_df = lowest_df[reduce(lambda x, y: x & y, all_f)]
        if lowest_df.empty:
            self.set_error("根据条件筛选之后的注释表为空，请重新填写筛选的条件！")
        for level in range(lowest_level):
            self.out_level(lowest_df, level, level == lowest_level - 1)
        self.end()
