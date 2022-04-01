# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

import unittest
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class SeqStructureAgent(Agent):
    def __init__(self, parent):
        super(SeqStructureAgent, self).__init__(parent)
        options = [
            {"name": "seq", "type": "string"},  # consists of gene id/gene name/transcript id, separating by comma
            {"name": "annotation_gtf", "type": "infile", "format": "ref_rna_v2.gtf"},
            {"name": "annotation_gff", "type": "infile", "format": "gene_structure.gff3"},
        ]
        self.add_option(options)
        self.step.add_steps("gene_structure")
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.gene_structure.start()
        self.step.update()

    def step_finish(self):
        self.step.gene_structure.finish()
        self.step.update()

    def check_options(self):
        if not self.option("seq"):
            raise OptionError("必须设置输入至少一个gene id/gene name/transcript id。")
        if not self.option("annotation_gtf").is_set and not self.option("annotation_gff").is_set:
            raise OptionError("必须设置输入gtf/gff格式的注释文件。")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "结果输出目录"]
        # ])
        # result_dir.add_regexp_rules([
        #     [r"disgenet_enrichment.xls$", "xls", "DisGeNET富集分析结果"]
        # ])
        super(SeqStructureAgent, self).end()


class SeqStructureTool(Tool):
    def __init__(self, config):
        super(SeqStructureTool, self).__init__(config)
        self._version = "v1.0"
        software_dir = self.config.SOFTWARE_DIR
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.seq_structure_path = self.config.PACKAGE_DIR + "/tool_lab/seq_structure.r"
        self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.1/bin/Rscript"
        # self.r_path = '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/miniconda/lib/R/bin/Rscript'

    def run(self):
        super(SeqStructureTool, self).run()
        self.run_seq_structure()
        self.set_output()
        self.end()

    def run_seq_structure(self):
        if self.option("annotation_gtf").is_set:
            annotation = self.option("annotation_gtf").prop["path"]
        if self.option("annotation_gff").is_set:
            annotation = self.option("annotation_gff").prop["path"]
        cmd = '{} {}'.format(self.r_path, self.seq_structure_path)
        cmd += ' -s {}'.format(self.option('seq'))
        cmd += ' -a {}'.format(annotation)
        cmd_name = 'plot_gene_structure'
        command = self.add_command(cmd_name, cmd, shell=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} 运行成功".format(cmd_name))
            else:
                self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd))
        else:
            self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd))

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        plots = glob.glob(self.work_dir + '/*.pdf')
        if not plots:
            self.logger.info("None of {} was found in your annotation file".format(self.option("seq")))
        seq_list = self.option('seq').split(",")
        file_list = list()
        for each in plots:
            name = os.path.basename(each)
            link = os.path.join(self.output_dir, name)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)
            seq = name.split(".pdf")[0]
            file_list.append(seq)
        for i in seq_list:
            if file_list and i not in file_list:
                self.logger.info("{} cannot be found in your annotation file".format(i))

