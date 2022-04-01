# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2019/4/29'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class BlastNtModule(Module):
    """
    比对nt库
    """

    def __init__(self, work_id):
        super(BlastNtModule, self).__init__(work_id)
        option = [
            {"name": "query", "type": "infile", "format": "sequence.fasta", "required": True},
            {"name": "sample_name", "type": "string", "required": True},
            {"name": "nt_table", "type": "outfile", "format": "sequence.profile_table"},  # blast结果
            {"name": "cir_table", "type": "outfile", "format": "sequence.profile_table"},  # 成环结果
            {"name": "fa_out", "type": "outfile", "format": "sequence.fasta"}
        ]
        self.add_option(option)
        self.run_tools = []  # edit tool list
        self.blast = self.add_tool('align.blast')
        # self.blast2 = self.add_tool("align.blast")
        self.circle = self.add_tool("bacgenome.circle_correct")
        self.havecircle_intheory = False
        self.iscircle_infact = False

    def check_options(self):
        """
        检查参数
        :return:
        """
        # edit options check
        return True

    def run(self):
        super(BlastNtModule, self).run()
        # self.on_rely([self.blast, self.blast2], self.set_output)
        super(BlastNtModule, self).run()
        self.blast.on("end", self.run_circle_correct)
        self.circle.on("end", self.set_output)
        self.run_blast()

    # def run_blast(self):
    #     opts = {
    #         'query': self.option('query'),
    #         'query_type': "nucl",
    #         'database': 'nt_20180827',
    #         'outfmt': 5,
    #         'blast': 'blastn',
    #         'memory': 50,
    #         'num_alignment': 1
    #     }
    #     self.blast.set_options(opts)
    #     self.blast.run()

    def run_blast(self):
        opts = {
            'query': self.option('query'),
            'query_type': "nucl",
            'database': 'nt_20180827',
            'outfmt': 6,
            'blast': 'blastn',
            'memory': 50,
            'num_alignment': 1
        }
        self.blast.set_options(opts)
        self.blast.run()

    def run_circle_correct(self):
        opts = {
            'blast_table': self.blast.option("outtable"),
            'query': self.option('query'),
            'sample_name': self.option("sample_name")
        }
        self.circle.set_options(opts)
        self.circle.run()

    # def blast_convert(self):
    #     '''
    #     # 不用这个函数
    #     读取blast结果，判断是否成环，如果结果中有成环的注释信息，做重复序列校正
    #     :return:
    #     '''
    #     if self.havecircle_intheory:
    #         self.run_circle_seq()
    #
    # def run_circle_seq(self):
    #     """
    #     截取两端
    #     :return:
    #     """
    #     opts = {
    #         "blast_table": self.blast.option("outtable"),
    #         "query": self.option("query"),
    #         "sample_name": self.option("sample_name")
    #         # "circle_table"  添加了成环判定的blast表，query根据circle_fasta中的id进行重命名
    #         # "circle_fasta"  环化调整后的fasta，id经过重命名
    #         # 如果一条序列的nt结果非环状，则跳过，只做id 的更新
    #     }


    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        # edit output
        self.option("fa_out").set_path(self.circle.option("fa_out").prop["path"])
        self.option("nt_table").set_path(self.circle.option("table_out").prop["path"])
        self.option("cir_table").set_path(self.circle.option("cir_out").prop["path"])
        self.logger.info("设置结果成功")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [",", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(BlastNtModule, self).end()
