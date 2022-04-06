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

class RfamStatAgent(Agent):
    """
    Rfam比对结果统计
    """
    def __init__(self, parent):
        super(RfamStatAgent, self).__init__(parent)
        options = [
            {"name": "blast_table", "type": "infile", "format": "align.blast.blast_table"},  # blast比对结果文件
            {"name": "bowtie_sam", "type": "infile", "format": "small_rna.common"},  # bowtie比对结果sam文件
            {"name": "filter_fa", "type": "outfile", "format": "small_rna.fasta"}, # 鉴定完已知miRNA的过滤文件
            {"name": "config", "type": "infile", "format": "small_rna.common"},  # 样本名转换文件
            {"name": "query", "type": "infile", "format": "small_rna.fasta"},  # 输入文件
        ]
        self.add_option(options)
        self.step.add_steps("rfam_stat")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.rfam_stat.start()
        self.step.update()

    def stepfinish(self):
        self.step.rfam_stat.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        """
        if (not self.option("blast_table").is_set) and (not self.option("bowtie_sam").is_set):
            raise OptionError("必须提供blast比对结果表或者bowtie比对结果表")
        if not self.option("query").is_set:
            raise OptionError("必须提供输入FASTA文件")
        self.set_resource()
        if not self.option("config").is_set:
            raise OptionError("必须提供配置文件")
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(RfamStatAgent, self).end()


class RfamStatTool(Tool):
    def __init__(self, config):
        super(RfamStatTool, self).__init__(config)
        self.python = "miniconda2/bin/python"
        self.rfam_stat = self.config.PACKAGE_DIR + "/small_rna/rfam_stat.py"
        self.rfam_seed = self.config.SOFTWARE_DIR + "/database/Rfam_12.3/Rfam.seed"

    def run(self):
        super(RfamStatTool, self).run()
        if self.option("blast_table").is_set:
            self.run_rfam_stat()
            self.parse_stat()
        elif self.option("bowtie_sam").is_set:
            self.parse_bowtie_result()
            self.bowtie_rfam_stat()
        self.set_output()
        self.end()

    def parse_bowtie_result(self):
        rfam_seed = self.config.SOFTWARE_DIR + "/database/Rfam_12.3/Rfam.seed"
        detail_out = self.output_dir + "/rfam_detail.xls"
        blast = self.option("bowtie_sam").prop["path"]

        seed = {}
        rfam_ac = ""
        rfam_id = ""
        rfam_de = ""
        rfam_tp = ""
        with open(rfam_seed, "r") as f:
            lines = f.readlines()
            for line in lines:
                m = re.match(r"#=GF\s+AC\s+(\S+)\s*$", line)
                m1 = re.match(r"#=GF\s+ID\s+(.*)\s*$", line)
                m2 = re.match(r"#=GF\s+DE\s+(.*)\s*$", line)
                m3 = re.match(r"#=GF\s+TP\s+(.*)\s*$", line)
                if m:
                    rfam_ac = m.group(1)
                if m1:
                    rfam_id = m1.group(1)
                if m2:
                    rfam_de = m2.group(1)
                if m3:
                    rfam_tp = m3.group(1)
                seed.update({rfam_ac:{"ID": rfam_id, "DE": rfam_de, "TP": rfam_tp}})

        with open(blast, "r") as f, open (detail_out ,"w") as w:
            w.write("query_name\tHit\tAC\tID\tTP\tDE\n")
            for line in f:
                items = line.strip().split("\t")
                if items[2] != "*":
                    rfam_ac = items[2].split(";")[0]
                    hit = items[2].split(";")[1]
                    query = items[0]
                    if rfam_ac in seed:
                        tp = seed[rfam_ac]["TP"]
                        if re.search('rRNA', tp):
                            columns = [query, hit, rfam_ac, seed[rfam_ac]["ID"], "rRNA", seed[rfam_ac]["DE"]]
                            w.write("\t".join(columns) + "\n")
                        elif re.search('ribozyme', tp):
                            columns = [query, hit, rfam_ac, seed[rfam_ac]["ID"], "rRNA", seed[rfam_ac]["DE"]]
                            w.write("\t".join(columns) + "\n")
                        elif re.search('tRNA', tp):
                            columns = [query, hit, rfam_ac, seed[rfam_ac]["ID"], "tRNA", seed[rfam_ac]["DE"]]
                            w.write("\t".join(columns) + "\n")
                        elif re.search('snRNA', tp):
                            if re.search('snoRNA', tp):
                                columns = [query, hit, rfam_ac, seed[rfam_ac]["ID"], "snoRNA", seed[rfam_ac]["DE"]]
                                w.write("\t".join(columns) + "\n")
                            else:
                                print tp
                                columns = [query, hit, rfam_ac, seed[rfam_ac]["ID"], "snRNA", seed[rfam_ac]["DE"]]
                                w.write("\t".join(columns) + "\n")
                    else:
                        pass
                else:
                    pass

    def bowtie_rfam_stat(self):
        summary_out = self.output_dir + "/rfam_summary.xls"
        detail_out = self.output_dir + "/rfam_detail.xls"
        with open (detail_out, "r") as f, open (summary_out, "w") as w:
            rRNA = {}
            tRNA = {}
            snoRNA = {}
            snRNA = {}
            samples = {}
            f.readline()
            for line in f:
                items = line.strip().split("\t")
                sample = items[0].split("_")[0]
                samples[sample] = sample
                num = int(items[0].split("_x")[1])
                if items[4] == "rRNA":
                    if sample in rRNA:
                        rRNA[sample] += num
                    else:
                        rRNA[sample] = num
                elif items[4] == "tRNA":
                    if sample in tRNA:
                        tRNA[sample] += num
                    else:
                        tRNA[sample] = num
                elif items[4] == "snoRNA":
                    if sample in snoRNA:
                        snoRNA[sample] += num
                    else:
                        snoRNA[sample] = num
                elif items[4] == "snRNA":
                    if sample in snRNA:
                        snRNA[sample] += num
                    else:
                        snRNA[sample] = num
            w.write("Samples\trRNA\ttRNA\tsnoRNA\tsnRNA\n")
            for key in samples:
                sample_name = self.parse_config(file=self.option("config").prop["path"], section="NAME", name=key)
                if not rRNA.has_key(key):
                    rRNA[key] = 0
                if not tRNA.has_key(key):
                    tRNA[key] = 0
                if not snoRNA.has_key(key):
                    snoRNA[key] = 0
                if not snRNA.has_key(key):
                    snRNA[key] = 0
                w.write("{}\t{}\t{}\t{}\t{}\n".format(sample_name, rRNA[key], tRNA[key], snoRNA[key], snRNA[key]))


    def run_rfam_stat(self):
        self.logger.info("开始处理blast运行结果文件")
        rfam_table = self.option("blast_table").prop["path"]
        detail_out = self.output_dir + "/rfam_detail.xls"
        cmd = "{} {} -i {} -db {} -o {}".format(self.python, self.rfam_stat, rfam_table, self.rfam_seed, detail_out)
        command = self.add_command("rfam_detail", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("处理blast运行结果文件完成")
        else:
            self.set_error("开始处理blast运行结果文件失败")

    def parse_stat(self):
        self.logger.info("开始进行rfam比对结果统计")
        detail_out = self.output_dir + "/rfam_detail.xls"
        summary_out = self.output_dir + "/rfam_summary.xls"
        with open (detail_out, "r") as f, open (summary_out, "w") as w:
            rRNA = {}
            tRNA = {}
            snoRNA = {}
            snRNA = {}
            other_sRNA = {}
            samples = {}
            f.readline()
            for line in f:
                items = line.strip().split("\t")
                sample = items[0].split("_")[0]
                samples[sample] = sample
                num = int(items[0].split("_x")[1])
                if items[6] == "rRNA":
                    if sample in rRNA:
                        rRNA[sample] += num
                    else:
                        rRNA[sample] = num
                elif items[6] == "tRNA":
                    if sample in tRNA:
                        tRNA[sample] += num
                    else:
                        tRNA[sample] = num
                elif items[6] == "snoRNA":
                    if sample in snoRNA:
                        snoRNA[sample] += num
                    else:
                        snoRNA[sample] = num
                elif items[6] == "snRNA":
                    if sample in snRNA:
                        snRNA[sample] += num
                    else:
                        snRNA[sample] = num
                # else:
                #     if other_sRNA.has_key(sample_name):
                #         other_sRNA[sample_name] += num
                #     else:
                #         other_sRNA[sample_name] = num
            w.write("Samples\trRNA\ttRNA\tsnoRNA\tsnRNA\n")
            for key in samples:
                sample_name = self.parse_config(file=self.option("config").prop["path"], section="NAME", name=key)
                if not rRNA.has_key(key):
                    rRNA[key] = 0
                if not tRNA.has_key(key):
                    tRNA[key] = 0
                if not snoRNA.has_key(key):
                    snoRNA[key] = 0
                if not snRNA.has_key(key):
                    snRNA[key] = 0
                w.write("{}\t{}\t{}\t{}\t{}\n".format(sample_name, rRNA[key], tRNA[key], snoRNA[key], snRNA[key]))
        self.logger.info("Rfam统计运行完成")

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
        rfam_detail = os.path.join(self.output_dir, "rfam_detail.xls")
        a = dict()
        with open (rfam_detail, "r") as f1:
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
            "id": "RfamStat" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna.srna.rfam_stat",
            "instant": False,
            "options": dict(
                bowtie_sam="/mnt/ilustre/users/sanger-dev/workspace/20191112/Single_BowtieRfam7380/BowtieRfam/output/bowtie_rfam.sam",
                config="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/srna_test/mouse/qc_file.config",
                query="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/srna_test/mouse/input.fa",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()