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
            {"name": "in_gtf", "type": "infile", "format": "gene_structure.gtf"},
            {"name": "gtf", "type": "outfile", "format": "gene_structure.gtf"},
            {"name": "bed", "type": "outfile", "format": "gene_structure.bed"},
            {"name": "is_duplicate", "type": "bool", "default": True},  # 是否有生物学重复
            {"name": "ref_genome", "type": "string"},
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
            raise OptionError("必须输入fastq文件参数", code = "35003201")
        if not self.option('fastq_dir').prop['has_list_file']:
            raise OptionError("fastq文件夹中必须含有一个名为list.txt的文件", code = "35003202")
        if self.option('fq_type') not in ['PE', 'SE']:
            raise OptionError("测试类型只能是PE或者SE", code = "35003203")
        if not self.option("in_gtf").is_set:
            raise OptionError("必须输入GTF文件", code = "35003204")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = "5G"


class FileCheckTool(Tool):
    """
    检查输入文件的格式是否符合要求
    """
    def __init__(self, config):
        super(FileCheckTool, self).__init__(config)
        self.samples = list()

    def check_fastq(self):
        self.logger.info("正在检测fastq_dir文件")
        if not self.option("fastq_dir").prop["has_list_file"]:
            raise OptionError("fastq文件夹中必须含有一个名为list.txt的文件", code = "35003205")
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
            raise OptionError("PE序列list文件应该包括文件名、样本名和左右端说明三列", code = "35003206")
        elif self.option('fq_type') == "SE" and col_num != 2:
            raise OptionError("SE序列list文件应该包括文件名、样本名两列", code = "35003207")
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
                    self.set_error("group表出错, 样本%s在fastq文件中未出现", variables = (gp), code = "35003208")
            if self.option("is_duplicate"):
                group_dict = self.option("group_table").prop['group_dict']
                group_names = set(group_dict.keys())
                samples = set(self.option("group_table").prop["sample"])
                coflict_names = group_names & samples
                if len(coflict_names) > 0:
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
        for vs in vs_list:
            for i in vs:
                if i not in con_samples:
                    con_samples.append(i)
            for cp in con_samples:
                if self.option('group_table').is_set:
                    if cp not in group_name:
                        self.set_error("control表出错，分组%s在fastq文件中未出现", variables = (cp), code = "35003209")
                else:
                    if cp not in self.samples:
                        self.set_error("control表出错，样本%s在fastq文件中未出现", variables = (cp), code = "35003210")
        self.logger.info("control文件检测完毕")

    def gtf_to_bed(self):
        self.logger.info("转换gtf文件为bed文件")
        new_gtf = self.work_dir + "/" + os.path.basename(self.option("in_gtf").prop["path"])
        if os.path.exists(new_gtf):
            os.remove(new_gtf)
        os.link(self.option("in_gtf").prop["path"], new_gtf)
        # 针对原核他们整理gtf文件经常出错，所以加了一个文件检查
        with open(new_gtf, 'r') as ng:
            for line in ng:
                if line.count('"-"') < 2:
                    break
            else:
                self.set_error("""你传入的gtf文件有问题，第九列的参考样式为 'transcript_id "peg.1"; gene_id "peg.1"; gene_name "peg.1"'""")
        self.option("gtf").set_path(new_gtf)
        self.option("gtf").to_bed()

    def filter_bed(self):
        with open(self.option("gtf").prop["path"] + ".bed", "r") as f, open(self.option("gtf").prop["path"] + ".filter.bed", "w") as w:
            for line in f:
                if int(line.split("\t")[2]) < 500000000:
                    w.write(line)
        self.option("bed").set_path(self.option("gtf").prop["path"] + ".filter.bed")

    def get_rna_seq(self):
        with open(self.option("in_gtf").prop["path"], "r") as f, open("rrna.bed", "w") as w:
            for line in f:
                cols = line.split("\t")
                if cols[2].lower() == "rrna":
                    if "16S " in cols[-1] or "23S " in cols[-1]:
                        if cols[-1].startwith("ID="):
                            gene_id = cols[-1].split(";")[0].strip().split("=")[-1]
                        else:
                            gene_id = cols[0] + "_" + cols[3] + "_" + cols[4]
                        w.write("\t".join([
                            cols[0],
                            str(int(cols[3]) - 1),
                            cols[4],
                            gene_id,
                            "0",
                            cols[6]
                        ]))

        self.getfasta_path  = "/bioinfo/rna/bedtools2-master/bin/bedtools"
        out_fa = self.out_dir + "/rrna.fa"
        if os.path.getsize("rrna.bed") > 0:
            cmd = "%s getfasta -fi %s -bed %s -fo %s -name -s" % (self.getfasta_path, self.option("ref_genome"), "rrna.bed", out_fa)
            self.logger.info("开始打印cmd命令!")
            self.logger.info(cmd)
            fa_cmd = self.add_command("rrna_fa",cmd).run()
            self.wait(fa_cmd)
            if fa_cmd.return_code == 0:
                self.logger.info("%s运行完成" % fa_cmd)
            else:
                self.set_error("%s运行出错", variables = (fa_cmd), code = "33705602")
            print "生成rrna的fa序列"

    def run(self):
        super(FileCheckTool, self).run()
        self.check_fastq()
        self.check_group()
        self.check_control()
        self.gtf_to_bed()
        self.filter_bed()
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
            "name": "ref_rna_v2.file_check",
            "instant": False,
            "options": dict(
                fq_type="PE",
                fastq_dir="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/demo_Mouse_small_8samples/rawdata",
                group_table="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/demo_Mouse_small_8samples/default_group.txt",
                control_file="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/demo_Mouse_small_8samples/control_file.txt",
                in_gtf="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/gtf/Mus_musculus.GRCm38.89.gtf",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
