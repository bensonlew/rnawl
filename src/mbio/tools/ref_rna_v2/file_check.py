# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import os
import re
import unittest

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool


class FileCheckAgent(Agent):
    """
    用于workflow开始之前对所输入的文件进行详细的内容检测
    """

    def __init__(self, parent):
        super(FileCheckAgent, self).__init__(parent)
        options = [
            {"name": "sample_num", "type": "string", 'default': "multiple"},
            {"name": "fastq_dir", "type": "infile", 'format': "sequence.fastq_dir"},  # fastq文件夹
            {"name": "fq_type", "type": "string"},  # PE OR SE
            {"name": "group_table", "type": "infile", "format": "sample.group_table"},  # 有生物学重复的时候的分组文件
            {"name": "control_file", "type": "infile", "format": "sample.control_table"},  # 对照组文件，格式同分组文件
            {"name": "in_gtf", "type": "infile", "format": "gene_structure.gtf"},
            {"name": "gtf", "type": "outfile", "format": "gene_structure.gtf"},
            {"name": "bed", "type": "outfile", "format": "gene_structure.bed"},
            {'name': 'is_duplicate', 'type': 'bool', 'default': True}
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
            raise OptionError("必须输入fastq文件参数", code="33705401")
        if not self.option('fastq_dir').prop['has_list_file']:
            raise OptionError('fastq文件夹中必须含有一个名为list.txt的文件名--样本名的对应文件', code="33705402")
        if self.option('fq_type') not in ['PE', 'SE']:
            raise OptionError("测试类型只能是PE或者SE", code="33705403")
        if not self.option("in_gtf").is_set:
            raise OptionError("必须输入GTF文件", code="33705404")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        self._memory = "8G"


class FileCheckTool(Tool):
    """
    检查输入文件的格式是否符合要求
    """

    def __init__(self, config):
        super(FileCheckTool, self).__init__(config)
        self.samples = list()

    def run(self):
        super(FileCheckTool, self).run()
        self.check_fastq()
        if self.option("group_table").is_set:
            self.check_group()
        if self.option("control_file").is_set:
            self.check_control()
        self.gtf_to_bed()
        self.make_true_bed()
        self.filter_bed()
        self.end()

    def check_fastq(self):
        self.logger.info("正在检测fastq_dir文件")
        if not self.option("fastq_dir").prop["has_list_file"]:
            raise OptionError('fastq文件夹中必须含有一个名为list.txt的文件名--样本名的对应文件', code="33705405")
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
            raise OptionError("PE序列list文件应该包括文件名、样本名和左右端说明三列", code="33705406")
        self.logger.info("fastq文件检测完毕")

    def get_list_info(self):
        list_path = self.option("fastq_dir").prop["path"] + "/list.txt"
        with open(list_path, "r") as f:
            col_num = len(f.readline().strip().split())
        return col_num

    def check_group(self):
        if self.option('group_table').is_set:
            self.logger.info("正在检测group文件")
            self.option("group_table").get_info()
            gp_sample = self.option("group_table").prop["sample"]
            self.logger.info("group_table中有{}个样本".format(str(len(set(gp_sample)))))
            self.logger.info("fastq_dir中list有{}个样本".format(str(len(set(self.samples)))))
            if set(gp_sample) !=  set(self.samples):
                self.set_error("group和fastq_dir中list文件不对应")
            for gp in gp_sample:
                if gp not in self.samples:
                    self.set_error("group表出错, 样本%s在fastq文件中未出现", variables=gp, code="33705408")
            if self.option("sample_num") == "multiple" and self.option("is_duplicate"):
                group_dict = self.option("group_table").prop['group_dict']
                group_names = set(group_dict.keys())
                samples = set(self.option("group_table").prop["sample"])
                coflict_names = group_names & samples
                if len(coflict_names) > 0:
                    group_snum = [len(group_dict[g]) for g in group_dict]
                    min_group_num = min(group_snum)
                    if min_group_num == 1:
                        # for g in group_dict:
                        #     if g != group_dict[g][0]:
                        #         self.set_error("当有生物学重复和无生物学重复混组时,无重组组样本名不可与其他组组名相同")
                        pass
                    else:
                        self.set_error("group表出错, 当有生物学重复时,不允许组名和样本名相同")
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
        else:
            group_name = list()
        for vs in vs_list:
            for i in vs:
                if i not in con_samples:
                    con_samples.append(i)
            for cp in con_samples:
                if self.option('group_table').is_set:
                    if cp not in group_name:
                        self.set_error("control表出错，分组%s在fastq文件中未出现", variables=(cp,), code="33705409")
                else:
                    if cp not in self.samples:
                        self.set_error("control表出错，样本%s在fastq文件中未出现", variables=(cp,), code="33705410")
        self.logger.info("control文件检测完毕")

    def gtf_to_bed(self):
        self.logger.info("转换gtf文件为bed文件")
        new_gtf = self.work_dir + "/" + os.path.basename(self.option("in_gtf").prop["path"])
        if os.path.exists(new_gtf):
            os.remove(new_gtf)
        # os.link(self.option("in_gtf").prop["path"], new_gtf)
        # 检查方向, 修改"."为"+"
        with open(self.option("in_gtf").prop["path"], "r") as in_gtf, open(new_gtf, 'w') as out_gtf:
            for line in in_gtf:
                cols = line.split('\t')
                if len(cols) >= 9 and cols[6] == ".":
                    cols[6] = "+"
                    line = "\t".join(cols)
                out_gtf.write(line)

        self.option("gtf").set_path(new_gtf)
        self.option("gtf").to_bed()

    def make_true_bed(self):
        """
        modified by qinjincheng at 20190325
        create a legal bed by ucsc utilities
        :return:
        """
        gtf = self.option('gtf').prop['path']
        genepred = '{}.genepred'.format(self.option('gtf').prop['path'])
        bed = '{}.bed'.format(self.option('gtf').prop['path'])
        gtfToGenePred = 'bioinfo/align/ucsc_tools/gtfToGenePred'
        genePredToBed = 'bioinfo/align/ucsc_tools/genePredToBed'
        cmd1 = '{} -ignoreGroupsWithoutExons {} {}'.format(gtfToGenePred, gtf, genepred)
        self.run_code('gtftogenepred', cmd1)
        cmd2 = '{} {} {}'.format(genePredToBed, genepred, bed)
        self.run_code('genepredtobed', cmd2)

    def run_code(self, cmd_name, cmd, shell=False):
        if shell:
            cmd = self.config.SOFTWARE_DIR + '/' + cmd
        command = self.add_command(cmd_name, cmd, shell=shell)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('succeed in running {}'.format(cmd_name))
        elif command.return_code is None:
            self.logger.warn('fail to run {}, try again'.format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info('succeed in rerunning {}'.format(cmd_name))
            else:
                self.set_error('fail to rerun {}, abord'.format(cmd_name))
        else:
            self.set_error('fail to run {}, abord'.format(cmd_name))

    def filter_bed(self):
        """
        modified by shicaiping at 20181225
        新增对bed文件的过滤，基于以下两点：
        1. 基因单条染色体序列长度超过500000000，无法进行Bamdistribution分析，过滤超过500000000部分的转录本
        2. 过滤exon存在交叉、包含等关系的转录本
        :return:
        """
        with open(self.option("gtf").prop["path"] + ".bed", "r") as f, open(
                self.option("gtf").prop["path"] + ".filter.bed", "w") as w:
            for line in f:
                if int(line.split("\t")[2]) < 500000000:
                    flag = 1
                    start_list = line.strip().split("\t")[11].split(",")[:-1]
                    length_list = line.strip().split("\t")[10].split(",")[:-1]
                    for index, value in enumerate(length_list[:-1]):
                        if (int(length_list[index]) + int(start_list[index])) >= int(start_list[index + 1]):
                            flag = 0
                            break
                    if flag == 1:
                        w.write(line)
        self.option("bed").set_path(self.option("gtf").prop["path"] + ".filter.bed")


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
            "name": "ref_rna_v2.file_check",
            "instant": False,
            "options": dict(
                fq_type="PE",
                fastq_dir="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_rna_v2/check_test",
                # group_table="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_rna_v2/check_test/min_default_group.txt",
                group_table="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_rna_v2/check_test/default_group.txt",
                control_file="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_rna_v2/check_test/control_file.txt",
                # control_file="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_rna_v2/check_test/mix_control_file.txt",
                in_gtf="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus"
                       "/Ensemble_release_89/gtf/Mus_musculus.GRCm38.89.gtf",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
