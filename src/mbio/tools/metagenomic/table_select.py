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


class TableSelectAgent(Agent):
    def __init__(self, parent):
        super(TableSelectAgent, self).__init__(parent)
        options = [
            {"name": "table", "type": "infile", "format": "meta_genomic.profile"},
            {"name": "cols", "type": "string"},
            {"name": "out_table", "type": "outfile", "format": "meta_genomic.profile"},
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
        self._cpu = 1
        self._memory = '4G'

    def end(self):
        super(TableSelectAgent, self).end()


class TableSelectTool(Tool):
    """
    对metaphlan3 和 kraken2 的输出结果挑选指定的列
    """
    def __init__(self, config):
        super(TableSelectTool, self).__init__(config)
        self.output = os.path.join(self.output_dir, "new_table.txt")

    def run(self):
        super(TableSelectTool, self).run()
        self.select()
        self.end()

    def select(self):
        cols = json.loads(self.option("cols"))
        key_map = {
            "d__": 'Domain' , "k__": 'Kingdom', 'p__': "Phylum", 'c__': "Class",
            'o__': "Order", 'f__': "Family", 'g__': "Genus", 's__': "Species"
        }
        for i in range(len(cols)):
            cols[i] = key_map[cols[i]] if cols[i] in key_map else cols[i]
        table = pd.read_csv(self.option("table").path, sep='\t', chunksize=500000)
        if os.path.exists(self.output):
            os.remove(self.output)
        for one in table:
            one = one[cols]
            one.to_csv(self.output, sep='\t', quoting=False, index=False, mode="a+")
        self.option("out_table", self.output)
