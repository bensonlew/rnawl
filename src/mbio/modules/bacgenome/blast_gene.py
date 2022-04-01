# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2019/5/9'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import pandas as pd

class BlastGeneModule(Module):
    """
    16s序列比对和管家基因序列比对
    """

    def __init__(self, work_id):
        super(BlastGeneModule, self).__init__(work_id)
        option = [
            {"name": "fa", "type": "infile", "format": "sequence.fasta", "required": True},
            {"name": "out_table", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(option)
        self.s_blast = self.add_tool('align.blast')

    def check_options(self):
        """
        检查参数
        :return:
        """
        # edit options check
        return True

    def run_16s_blast(self):
        opts = {
            "query": self.option('fa'),
            "query_type": "nucl",
            "database": "16s",
            "outfmt": 6,
            "blast": "blastn",
            "memory": 50,
            "num_alignment": 10
        }
        self.s_blast.set_options(opts)
        self.s_blast.run()

    def process_blast_table(self):
        data = pd.read_table(self.s_blast.option("outtable").prop["path"], header=0, names=["Score", "Evalue", "align_length",
            "Identity (%)", "Similarity (%)", "Query_id", "q.length", "q.start", "q.end", "q.direction", "hit_name",
            "s.length", "s.start",  "s.end", "s.direction", "Organism"])
        data["Coverage (%)"] = abs(data["s.end"] - data["s.start"]) / data["s.length"].astype("float") * 100  # 查询序列coverage改为参考序列coverage
        data.reindex(columns=["Organism", "Identity (%)", "Coverage (%)", "Evalue", "Score"])\
            .to_csv(self.output_dir + "/16s_blast.xls", sep="\t", index=False, header=True)
        self.set_output()

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        self.option("out_table").set_path(self.output_dir + "/16s_blast.xls")
        self.logger.info("设置结果成功")
        self.end()

    def run(self):
        super(BlastGeneModule, self).run()
        self.s_blast.on("end", self.process_blast_table)
        self.run_16s_blast()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [",", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(BlastGeneModule, self).end()
