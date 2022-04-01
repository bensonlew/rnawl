# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2019/4/30'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagenomic.common import link_file


class ChrMulCircleModule(Module):
    """
    先做nt比对，然后根据比对结果环化校正，dnaA校正
    """

    def __init__(self, work_id):
        super(ChrMulCircleModule, self).__init__(work_id)
        option = [
            {"name": "scf_input1", "type": "infile", "format": "sequence.fasta_dir", "required": True},  # 主要的拼接结果
            {"name": "scf_input2", "type": "infile", "format": "sequence.fasta_dir", "required": False},  # 次要的拼接结果
            {"name": "samples", "type": "string", "required": True},  # 样品列表
            {"name": "scf_out1", "type": "outfile", "format": "sequence.fasta_dir"},  # 主要的拼接结果
            {"name": "scf_out2", "type": "outfile", "format": "sequence.fasta_dir"},
            {"name": "nt_table", "type": "outfile", "format": "bacgenome.simple_dir"}  # blast结果
        ]
        self.add_option(option)
        self.run_tools = []
        self.sample_path = {}
        self.second_path_list = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        # edit options check
        return True

    def run(self):
        super(ChrMulCircleModule, self).run()
        self.mv_path()
        self.run_circle()

    def mv_path(self):
        for sample in self.option("samples").split(","):
            first_path = os.path.join(self.option("scf_input1").prop["path"], sample)

            path1 = self.have_path(first_path)
            path2 = None
            second_path = ""
            if self.option("scf_input2").is_set:
                second_path = os.path.join(self.option("scf_input2").prop["path"], sample)
                path2 = self.have_path(second_path)
            if path1:
                self.sample_path[sample] = path1
                if path2:
                    self.second_path_list.append(path2)  # 用于放到second路径下的数据
            elif path2:
                self.sample_path[sample] = path2
            if os.path.isfile(os.path.join(first_path, sample + ".scaffold.fna")):
                self.sample_path[sample] = os.path.join(first_path, sample + ".scaffold.fna")
            elif os.path.isfile(os.path.join(first_path, sample + ".scaf.fna")):
                self.sample_path[sample] = os.path.join(first_path, sample + ".scaf.fna")
            elif os.path.isfile(os.path.join(second_path, sample + ".scaffold.fna")):
                self.sample_path[sample] = os.path.join(second_path, sample + ".scaffold.fna")
            elif os.path.isfile(os.path.join(second_path, sample + ".scaf.fna")):
                self.sample_path[sample] = os.path.join(second_path, sample + ".scaf.fna")

    def have_path(self, path):
        if os.path.isfile(path + ".scaffold.fna") and os.path.getsize(path + ".scaffold.fna"):
            return path + ".scaffold.fna"
        elif os.path.isfile(path + ".scaf.fna") and os.path.getsize(path + ".scaf.fna"):
            return path + ".scaf.fna"
        else:
            return False

    def run_circle(self):
        for sample in self.sample_path:
            tool = self.add_module("bacgenome.blast_nt")
            tool.set_options({
                "query": self.sample_path[sample],
                "sample_name": sample
            })
            self.run_tools.append(tool)
        self.on_rely(self.run_tools, self.set_output)
        for tool in self.run_tools:
            tool.run()

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        first_path = os.path.join(self.output_dir, "first")
        second_path = os.path.join(self.output_dir, "second")
        table_path = os.path.join(self.output_dir, "nt")
        for path in [first_path, second_path, table_path]:
            if not os.path.isdir(path):
                os.mkdir(path)
        for tool in self.run_tools:
            link_file(tool.option("fa_out").prop["path"], os.path.join(first_path, tool.option("sample_name") + ".scaffold.fna"))
            link_file(tool.option("nt_table").prop["path"], os.path.join(table_path, tool.option("sample_name") + ".nt.xls"))
            link_file(tool.option("cir_table").prop["path"], os.path.join(table_path, tool.option("sample_name") + ".circle.xls"))
        for path in self.second_path_list:
            file_name = os.path.basename(path)
            file_name.replace(".scaf.fna", ".scaffold.fna")
            link_file(path, os.path.join(second_path, file_name))
        self.option("scf_out1").set_path(first_path)
        self.option("scf_out2").set_path(second_path)
        self.option("nt_table").set_path(table_path)
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
        super(ChrMulCircleModule, self).end()
