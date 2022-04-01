# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import os,glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
from mbio.packages.whole_transcriptome.utils import runcmd
import unittest
import pandas as pd
import numpy as np
from sklearn import decomposition, preprocessing


class HeatmapAgent(Agent):
    """
    Used for convert RNA-seq expression units convert, count/CPM/RPM/TPM/FPKM/RPKM are implemented.

    RPM/CPM: Reads/Counts of exon model per million mapped reads
    RPM/CPM=Total exon reads/ Mapped reads(Millions)

    RPKM/FPKM: Reads/Fragments Per Kilobase of exon model per Million mapped reads
    RPKM/FPKM=Total exon reads/[Mapped reads(Millions)*Exon length(Kb)]

    TPM is like RPKM and FPKM, except the order of operation is switched.
    """
    def __init__(self, parent):
        super(HeatmapAgent, self).__init__(parent)
        options = [
            {'name': 'otutable', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'log_change', 'type': 'string'},
            {'name': 'scale', 'type': 'string'},
            {'name': 'scale_table', 'type': 'outfile', 'format': 'ref_rna_v2.common'}

        ]
        self.add_option(options)


    def check_options(self):
        # """
        # 重写参数检测函数
        # :return:
        # """
        # if not self.option('exp_matrix').is_set:
        #     raise OptionError('必须输入表达定量文件')
        # if not self.option('convert_type') in ["tpm", "cpm", "fpkm", "TMM", "TMMwzp", "RLF", "uqua", "DESeq2"]:
        #     raise OptionError('不支持该标准化方法')
        # return True
        pass
    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "3G"

    def end(self):

        super(HeatmapAgent, self).end()


class HeatmapTool(Tool):
    def __init__(self, config):
        super(HeatmapTool, self).__init__(config)
        self.file = {
            'scale_table': os.path.join(self.output_dir, 'scale_table.txt')
        }

    def run(self):
        """
        运行
        :return:
        """
        super(HeatmapTool, self).run()
        self.run_heatmap()
        self.set_output()
        self.end()

    def run_heatmap(self):
        table = pd.read_table(self.option('otutable').path, sep='\t', index_col=0)
        table.index.name = 'feature'
        if self.option('log_change') == 'loge':
            table = np.log(table + 1)
        elif self.option('log_change') == 'log2':
            table = np.log2(table + 1)
        elif self.option('log_change') == 'log10':
            table = np.log10(table + 1)
        else:
            table = table
        if self.option('scale') == '0_mean':
            scaler = preprocessing.StandardScaler(with_std=False)
            scaler.fit(table.T)
            table_scale = pd.DataFrame(scaler.transform(table.T)).T
            table_scale.columns = table.columns
            table_scale.insert(0, 'feature', table.index.tolist())
            table_scale = table_scale.set_index('feature')
        elif self.option('scale') == 'zscore':
            scaler = preprocessing.StandardScaler()
            scaler.fit(table.T)
            table_scale = pd.DataFrame(scaler.transform(table.T)).T
            table_scale.columns = table.columns
            table_scale.insert(0, 'feature', table.index.tolist())
            table_scale = table_scale.set_index('feature')
        else:
            table_scale = table
        table_scale.to_csv(self.file['scale_table'], index=True, header=True, sep='\t')



    def set_output(self):
        self.option('scale_table').set_path(self.file['scale_table'])

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "heatmap_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.heatmap",
            "options": dict(
                otutable="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/cluster/otu_table.xls",
                log_change="None",
                scale='zscore'

            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
