# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.files.sequence.fastq import FastqFile
from mbio.files.sequence.file_sample import FileSampleFile
import os
import re
import unittest

class FileCheckAgent(Agent):
    """
    用于workflow开始之前对所输入的文件进行详细的内容检测
    """

    def __init__(self, parent):
        super(FileCheckAgent, self).__init__(parent)
        options = [
            {"name": "fastq_dir", "type": "infile", 'format': "sequence.fastq_dir"},  # fastq文件夹
            {"name": "fq_type", "type": "string"},  # PE OR SE
            {"name": "group_table", "type": "infile", "format": "sample.group_table"}, # 有生物学重复的时候的分组文件
            {"name": "control_file", "type": "infile", "format": "sample.control_table"}, # 对照组文件，格式同分组文件
            {"name": "assembly_file", "type": "infile", 'format': "denovo_rna_v2.trinity_fasta"},  # trinity组装结果文件
            {"name": "gene_to_trans", "type": "string"},  # 基因和转录本对应关系文件
        ]
        self.add_option(options)
        self.step.add_steps("file_check")
        self.on('start', self.start_file_check)
        self.on('end', self.end_file_check)

    def start_file_check(self):
        self.step.file_check.start()
        self.step.update()

    def end_file_check(self):
        self.step.file_check.finish()
        self.step.update()

    def check_option(self):
        if not self.option('fastq_dir'):
            raise OptionError("必须输入fastq文件参数", code = "33705401")
        if not self.option('fastq_dir').prop['has_list_file']:
            raise OptionError('fastq文件夹中必须含有一个名为list.txt的文件名--样本名的对应文件', code = "33705402")
        if self.option('fq_type') not in ['PE', 'SE']:
            raise OptionError("测试类型只能是PE或者SE", code = "33705403")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = "15G"


class FileCheckTool(Tool):
    """
    检查输入文件的格式是否符合要求
    """
    def __init__(self, config):
        super(FileCheckTool, self).__init__(config)
        self.samples = list()

    def get_unigene(self):
        self.option("assembly_file").set_gene2tran(self.option("gene_to_trans"))
        self.option("assembly_file").get_unigene('unigene.fasta', 'Trinity')

    def check_fastq(self):
        self.logger.info("正在检测fastq_dir文件")
        if not self.option("fastq_dir").prop["has_list_file"]:
            raise OptionError('fastq文件夹中必须含有一个名为list.txt的文件名--样本名的对应文件', code = "33705405")
        self.samples = self.option("fastq_dir").prop["samples"]
        for sample in self.samples:
            if len(sample) >= 20:
                raise OptionError("The length of sample name {} is longer than 20 characters.".format(sample))
            match = re.match(r'[^a-zA-Z]', sample)
            if match is not None:
                raise OptionError('Sample name {} must be start with english letter.'.format(sample))
            for char in sample:
                match = re.search(r'([^a-zA-z0-9_])', char)
                if match is not None:
                    raise OptionError('{} contains special character: {}'.format(sample, match.group()))
        list_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
        col_num = len(open(list_path, "r").readline().split())
        if self.option('fq_type') == "PE" and col_num != 3:
            raise OptionError("PE序列list文件应该包括文件名、样本名和左右端说明三列", code = "33705406")
        elif self.option('fq_type') == "SE" and col_num != 2:
            raise OptionError("SE序列list文件应该包括文件名、样本名两列", code = "33705407")
        self.logger.info("fastq文件检测完毕")

    def get_list_info(self):
        list_path = self.option("fastq_dir").prop["path"] + "/list.txt"
        with open(list_path, "r") as l:
            col_num = len(l.readline().strip().split())
        return col_num

    def check_group(self):
        if self.option('group_table').is_set:
            self.logger.info("正在检测group文件")
            self.option("group_table").get_info()
            gp_sample = self.option("group_table").prop["sample"]
            self.logger.info("group_table中有{}个样本".format(str(len(set(gp_sample)))))
            self.logger.info("fastq_dir中list有{}个样本".format(str(len(set(self.samples)))))
            if set(gp_sample) != set(self.samples):
                self.set_error("group和fastq_dir中list文件不对应")
            for gp in gp_sample:
                if gp not in self.samples:
                    self.set_error("group表出错, 样本%s在fastq文件中未出现", variables = (gp), code = "33705408")
        else:
            self.logger.info("未检测到group文件， 跳过...")
        self.logger.info("group文件检测完毕")

    def check_control(self):
        self.logger.info("正在检测control文件")
        vs_list = self.option("control_file").prop["vs_list"]
        con_samples = []
        if self.option('group_table').is_set:
            group_scheme = self.option('group_table').prop['group_scheme'][0]
            group_name = self.option('group_table').get_group_name(group_scheme)
        for vs in vs_list:
            for i in vs:
                if i not in con_samples:
                    con_samples.append(i)
            for cp in con_samples:
                if self.option('group_table').is_set:
                    if cp not in group_name:
                        self.set_error("control表出错，分组%s在fastq文件中未出现", variables = (cp), code = "33705409")
                else:
                    if cp not in self.samples:
                        self.set_error("control表出错，样本%s在fastq文件中未出现", variables = (cp), code = "33705410")
        self.logger.info("control文件检测完毕")

    def run(self):
        super(FileCheckTool, self).run()
        self.check_fastq()
        self.check_group()
        self.check_control()
        if self.option("assembly_file").is_set:
            self.get_unigene()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "FileCheck_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "whole_transcriptome.smallrna.file_check",
            "instant": False,
            "options": dict(
                fq_type="SE",
                fastq_dir="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/whole_transcriptome/demo_human/smallrna/rawdata",
                group_table="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/whole_transcriptome/demo_human/group.txt",
                control_file="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/whole_transcriptome/demo_human/control.txt"
                #in_gtf="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
