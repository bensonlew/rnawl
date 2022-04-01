# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
from biocluster.config import Config
import ConfigParser
import unittest
from skbio.parse.sequences import parse_fasta

class RepeatStatAgent(Agent):
    """
    repeat比对结果统计
    """
    def __init__(self, parent):
        super(RepeatStatAgent, self).__init__(parent)
        options = [
            {"name": "blast_table", "type": "infile", "format": "align.blast.blast_table"},  # rfam比对结果文件
            {"name": "bowtie_sam", "type": "infile", "format": "small_rna.common"},  # bowtie比对结果sam文件
            {"name": "filter_fa", "type": "outfile", "format": "small_rna.fasta"}, # 鉴定完已知miRNA的过滤文件
            {"name": "config", "type": "infile", "format": "small_rna.common"},  # 样本名转换文件
            {"name": "query", "type": "infile", "format": "small_rna.fasta"},  # 输入文件
            {"name": "gff", "type": "infile", "format": "gene_structure.gff3"},  # 输入文件
        ]
        self.add_option(options)
        self.step.add_steps("repeat_stat")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.repeat_stat.start()
        self.step.update()

    def stepfinish(self):
        self.step.repeat_stat.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        """
        if (not self.option("blast_table").is_set) and (not self.option("bowtie_sam").is_set):
            raise OptionError("必须提供blast比对结果表或者bowtie比对结果表")
        if not self.option("query").is_set:
            raise OptionError("必须提供输入FASTA文件")
        if not self.option("config").is_set:
            raise OptionError("必须提供配置文件")
        self.set_resource()
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(RepeatStatAgent, self).end()


class RepeatStatTool(Tool):
    def __init__(self, config):
        super(RepeatStatTool, self).__init__(config)
        self.python = "program/Python/bin/python"
        self.repeat_stat = self.config.PACKAGE_DIR + "/small_rna/repeat_stat.py"

    def run(self):
        super(RepeatStatTool, self).run()
        if self.option("blast_table").is_set:
            self.run_repeat_stat()
            self.parse_stat()
        elif self.option("bowtie_sam").is_set:
            self.parse_bowtie_result()
            self.bowtie_repeat_stat()
        self.set_output()
        self.end()

    def parse_bowtie_result(self):
        detail_out = self.output_dir + "/repeat_detail.xls"
        blast = self.option("bowtie_sam").prop["path"]
        gff = self.option("gff").prop["path"]

        ssr = {}
        with open(gff, "r") as f:
            head = f.readline()
            for line in f:
                items = line.strip().split("\t")
                info = re.search(r'ID=(\S+);.*;Class=(\S+).*;PercDiv=.*', items[8])
                ID = info.group(1)
                CLASS = info.group(2).split("/")[0]
                ssr[ID] = CLASS

        query = {}
        with open(blast, "r") as f, open (detail_out ,"w") as w:
            w.write("query_name\tHit\tClass\n")
            for line in f:
                items = line.strip().split("\t")
                if items[2] != "*" and items[0] not in query:
                    query[items[0]] = 1
                    columns = [items[0], items[2], ssr[items[2]]]
                    w.write("\t".join(columns) + "\n")
                else:
                    pass

    def bowtie_repeat_stat(self):
        detail_out = self.output_dir + "/repeat_detail.xls"
        summary_out = self.output_dir + "/repeat_summary.xls"
        with open (detail_out, "r") as f, open (summary_out, "w") as w:
            SINE = {}
            LINE = {}
            LTR = {}
            DNA = {}
            Simple_repeat = {}
            Satellite = {}
            Other_repeats = {}
            samples = {}
            f.readline()
            for line in f:
                items = line.strip().split("\t")
                sample = items[0].split("_")[0]
                samples[sample] = sample
                num = int(items[0].split("_x")[1])
                if items[2] == "SINE":
                    if SINE.has_key(sample):
                        SINE[sample] += num
                    else:
                        SINE[sample] = num
                elif items[2] == "LINE":
                    if LINE.has_key(sample):
                        LINE[sample] += num
                    else:
                        LINE[sample] = num
                elif items[2] == "LTR":
                    if LTR.has_key(sample):
                        LTR[sample] += num
                    else:
                        LTR[sample] = num
                elif items[2] == "DNA":
                    if DNA.has_key(sample):
                        DNA[sample] += num
                    else:
                        DNA[sample] = num
                elif items[2] == "srpRNA":
                    if Simple_repeat.has_key(sample):
                        Simple_repeat[sample] += num
                    else:
                        Simple_repeat[sample] = num
                elif items[2] == "Satellite":
                    if Satellite.has_key(sample):
                        Satellite[sample] += num
                    else:
                        Satellite[sample] = num
                else:
                    if Other_repeats.has_key(sample):
                        Other_repeats[sample] += num
                    else:
                        Other_repeats[sample] = num
            w.write("Samples\tSINE\tLINE\tLTR\tDNA\tSimple_repeat\tSatellite\tOther_repeats\n")
            for key in samples:
                sample_name = self.parse_config(file=self.option("config").prop["path"], section="NAME", name=key)
                if not SINE.has_key(key):
                    SINE[key] = 0
                if not LINE.has_key(key):
                    LINE[key] = 0
                if not LTR.has_key(key):
                    LTR[key] = 0
                if not DNA.has_key(key):
                    DNA[key] = 0
                if not Simple_repeat.has_key(key):
                    Simple_repeat[key] = 0
                if not Satellite.has_key(key):
                    Satellite[key] = 0
                if not Other_repeats.has_key(key):
                    Other_repeats[key] = 0
                w.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sample_name, SINE[key], LINE[key], LTR[key], DNA[key], Simple_repeat[key], Satellite[key], Other_repeats[key]))

    def run_repeat_stat(self):
        self.logger.info("开始进行repeat比对结果统计")
        repeat_table = self.option("blast_table").prop["path"]
        detail_out = self.output_dir + "/repeat_detail.xls"
        gff = self.option("gff").prop["path"]
        cmd = "{} {} -i {} -f {} -o {}".format(self.python, self.repeat_stat, repeat_table, gff, detail_out)
        command = self.add_command("repeat_detail", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("repeat统计运行完成")
        else:
            self.set_error("repeat统计运行失败")

    def parse_stat(self):
        detail_out = self.output_dir + "/repeat_detail.xls"
        summary_out = self.output_dir + "/repeat_summary.xls"
        with open (detail_out, "r") as f, open (summary_out, "w") as w:
            SINE = {}
            LINE = {}
            LTR = {}
            DNA = {}
            Simple_repeat = {}
            Satellite = {}
            Other_repeats = {}
            samples = {}
            f.readline()
            for line in f:
                items = line.strip().split("\t")
                sample = items[0].split("_")[0]
                samples[sample] = sample
                num = int(items[0].split("_x")[1])
                if items[4] == "SINE":
                    if SINE.has_key(sample):
                        SINE[sample] += num
                    else:
                        SINE[sample] = num
                elif items[4] == "LINE":
                    if LINE.has_key(sample):
                        LINE[sample] += num
                    else:
                        LINE[sample] = num
                elif items[4] == "LTR":
                    if LTR.has_key(sample):
                        LTR[sample] += num
                    else:
                        LTR[sample] = num
                elif items[4] == "DNA":
                    if DNA.has_key(sample):
                        DNA[sample] += num
                    else:
                        DNA[sample] = num
                elif items[4] == "srpRNA":
                    if Simple_repeat.has_key(sample):
                        Simple_repeat[sample] += num
                    else:
                        Simple_repeat[sample] = num
                elif items[4] == "Satellite":
                    if Satellite.has_key(sample):
                        Satellite[sample] += num
                    else:
                        Satellite[sample] = num
                else:
                    if Other_repeats.has_key(sample):
                        Other_repeats[sample] += num
                    else:
                        Other_repeats[sample] = num
            w.write("Samples\tSINE\tLINE\tLTR\tDNA\tSimple_repeat\tSatellite\tOther_repeats\n")
            for key in samples:
                sample_name = self.parse_config(file=self.option("config").prop["path"], section="NAME", name=key)
                if not SINE.has_key(key):
                    SINE[key] = 0
                if not LINE.has_key(key):
                    LINE[key] = 0
                if not LTR.has_key(key):
                    LTR[key] = 0
                if not DNA.has_key(key):
                    DNA[key] = 0
                if not Simple_repeat.has_key(key):
                    Simple_repeat[key] = 0
                if not Satellite.has_key(key):
                    Satellite[key] = 0
                if not Other_repeats.has_key(key):
                    Other_repeats[key] = 0
                w.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sample_name, SINE[key], LINE[key], LTR[key], DNA[key], Simple_repeat[key], Satellite[key], Other_repeats[key]))

    def parse_config(self, file=None, section=None, name=None):
        config = ConfigParser.ConfigParser()
        config.read(file)
        return config.get(section, name)


    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        self.logger.info("开始生成过滤后的fasta文件")
        input_fa = self.option("query").prop["path"]
        filter_fa = os.path.join(self.output_dir, "filtered.fa")
        rpeat_detail = os.path.join(self.output_dir, "repeat_detail.xls")
        a = dict()
        with open (rpeat_detail, "r") as f1:
            for line in f1:
                if re.compile(r'\S+_\d+_x\d+\s+').match(line):
                    a[line.strip().split()[0]] = 1
        with open(filter_fa, "w") as w:
            for seq_id, seq_sequence in parse_fasta(input_fa):
                if seq_id not in a:
                    w.write('>%s\n%s\n' % (seq_id, seq_sequence))
        self.option('filter_fa', filter_fa)
        self.logger.info("fasta过滤文件处理完毕")

class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "RepeatStat" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna.srna.repeat_stat",
            "instant": False,
            "options": dict(
                bowtie_sam="/mnt/ilustre/users/sanger-dev/workspace/20191115/Single_BowtieRepeat8106/BowtieRepeat/output/bowtie_repeat.sam",
                config="/mnt/ilustre/users/sanger-dev/workspace/20191114/Smallrna_tsg_36156/MirnaQc/output/clean_data/qc_file.config",
                query="/mnt/ilustre/users/sanger-dev/workspace/20191114/Smallrna_tsg_36156/Srna/RfamStat/output/filtered.fa",
                gff="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/Annotation_v2/repeatmasker/repeatmasker_merge.SSR.gff",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()