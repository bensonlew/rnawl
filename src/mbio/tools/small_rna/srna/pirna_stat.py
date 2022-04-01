# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'


import glob
import os
import re
import shutil
import subprocess
import unittest
import pandas as pd
from biocluster.agent import Agent
from biocluster.config import Config
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool




class PirnaStatAgent(Agent):
    """
    已知miRNA鉴定
    """

    def __init__(self, parent):
        super(PirnaStatAgent, self).__init__(parent)
        options = [
            {'name': 'pir_list', 'type': 'infile', 'format': 'small_rna.common'},
            {'name': 'known_pirna', 'type': 'outfile', 'format': 'small_rna.common'},

        ]
        self.add_option(options)
        self.step.add_steps("known_pirna")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.known_pirna.start()
        self.step.update()

    def stepfinish(self):
        self.step.known_pirna.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        self.set_resource()
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = '50G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "已知miRNA鉴定结果目录"]
        ])
        result_dir.add_regexp_rules([
            [r"pdfs", "", "已知miRNA二级结构文件目录"],
            ["known_mirna_count.xls", "xls", "已知miRNA定量count表"],
            ["known_mirna_norm.xls", "xls", "已知miRNA定量norm表"],
            ["filtered.fa", "", "过滤FASTA文件"],
        ])
        super(PirnaStatAgent, self).end()


class PirnaStatTool(Tool):
    def __init__(self, config):
        super(PirnaStatTool, self).__init__(config)
        self.file = {
            'known_pirna': os.path.join(self.output_dir, 'known_pirna.txt')
        }

    def run(self):
        """
        运行
        :return:
        """
        super(PirnaStatTool, self).run()
        self.merge()
        self.set_output()
        self.end()

    def merge(self):
        pir_list = list()
        with open(self.option('pir_list').path, 'r') as l:
            for line in l.readlines():
                pir_path = line.strip()
                size = os.path.getsize(pir_path)
                if size == 0:
                    continue
                pir_table = pd.read_table(pir_path, sep='\t', header=0)
                pir_list.append(pir_table)
        pirna_df = pd.concat(pir_list, axis=0)
        df = pirna_df.groupby('sample')
        df_list = list()
        for i in df:
            # i1 = i[1].set_index('pirna_id')
            i1 = i[1]
            i1.rename(columns={'count': i[0]}, inplace=True)
            i1.drop(labels=['sample'], axis=1, inplace=True)
            df_list.append(i1)
        print df_list
        df_merge = reduce(lambda left, right: pd.merge(left, right, how='outer', on=['pirna_id']), df_list)
        df_merge = df_merge.fillna(0)
        df = df_merge.groupby('pirna_id').sum()
        df.to_csv(self.file['known_pirna'], sep='\t', header=True, index=True)

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        self.option('known_pirna').set_path(self.file['known_pirna'])


class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "KnownPirna_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna.srna.known_pirna",
            "instant": False,
            "options": dict(
                species="ssc",
                clean_fa="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/smallrna/test/uniq_1000000.fa",
                config="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/smallrna/test/qc_file.config",
                # #mirdeep2_version='0.1.3',
                # whole=False,
                # database="pmiren"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
