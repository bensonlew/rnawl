# -*- coding: utf-8 -*-
# __author__ = 'zhengyuan'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from collections import defaultdict
import os

class AddTaxAgent(Agent):

    def __init__(self, parent):
        super(AddTaxAgent, self).__init__(parent)
        options = [
            {"name": "abund", "type": "infile", "format": "sequence.profile_table"},    # 输入文件，gene丰度表
            {"name": "new_abund", "type": "outfile", "format": "sequence.profile_table"}   # 输出文件，将gene转为物种分类
        ]
        self.add_option(options)
        self.step.add_steps("tax")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def check_options(self):
        if not self.option("abund").is_set:
            raise OptionError("必须输入gene表或blast结果表")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def stepstart(self):
        self.step.tax.start()
        self.step.update()

    def stepfinish(self):
        self.step.tax.finish()
        self.step.update()

    def end(self):
        super(AddTaxAgent, self).end()

class AddTaxTool(Tool):

    def __init__(self, config):
        super(AddTaxTool, self).__init__(config)
        self.level_id = 0
        self.ref_db = os.path.join(self.config.SOFTWARE_DIR, "database/Cgc_taxon/gene_nr_anno.xls")

    def get_map(self):
        with open(self.ref_db, "r") as read_file:
            target = [{},{},{},{},{},{},{},{},{}]
            lines = read_file.readlines()
            for line in lines[1:]:
                tmp = line.strip().split("\t")
                if self.level_id == 0:
                    target[0][tmp[0]] = ";".join(tmp[2:10])
                else:
                    target[self.level_id][tmp[0]] = tmp[self.level_id + 1]
        return target

    def anno_prof(self):
        import pandas as pd
        import numpy as np
        gene_to_anno = self.get_map()
        profile = pd.read_table(self.option("abund").path, header=0, low_memory=False)
        profile['annotation'] = profile['gene'].map(gene_to_anno[self.level_id]).replace(np.nan, "no_rank")
        new_file = profile.groupby(profile['annotation']).sum()
        new_file.to_csv(os.path.join(self.output_dir, "new_abund.xls"), sep='\t')

    def set_output(self):
        self.logger.info("设置结果目录")
        self.option("new_abund").set_path(self.output_dir + '/new_abund.xls')
        self.logger.info("设置结果目录成功")

    def run(self):
        super(AddTaxTool, self).run()
        self.anno_prof()
        self.set_output()
        self.end()
