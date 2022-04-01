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

class IntronStatAgent(Agent):
    """
    repeat比对结果统计
    """
    def __init__(self, parent):
        super(IntronStatAgent, self).__init__(parent)
        options = [
            {"name": "blast_table", "type": "infile", "format": "align.blast.blast_table"},  # rfam比对结果文件
            {"name": "filter_fa", "type": "outfile", "format": "small_rna.fasta"}, # 鉴定完已知miRNA的过滤文件
            {"name": "config", "type": "infile", "format": "small_rna.common"},  # 样本名转换文件
            {"name": "query", "type": "infile", "format": "small_rna.fasta"},  # 输入文件
        ]
        self.add_option(options)
        self.step.add_steps("intron_stat")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.intron_stat.start()
        self.step.update()

    def stepfinish(self):
        self.step.intron_stat.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        """
        if not self.option("blast_table").is_set:
            raise OptionError("必须提供repeat比对结果表")
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
        super(IntronStatAgent, self).end()


class IntronStatTool(Tool):
    def __init__(self, config):
        super(IntronStatTool, self).__init__(config)
        self.python = "program/Python/bin/python"

    def run(self):
        super(IntronStatTool, self).run()
        self.parse_blast()
        self.set_output()
        self.end()

    def parse_blast(self):
        self.logger.info("开始对intron序列比对结果进行统计")
        intron_table = self.option("blast_table").prop["path"]
        detail_out = self.output_dir + "/intron_detail.xls"
        summary_out = self.output_dir + "/intron_summary.xls"
        with open(intron_table, "r") as f, open (detail_out ,"w") as w1, open (summary_out, "w") as w2:
            query = {}
            samples = {}
            intron = {}
            head = f.readline()
            w1.write("query_name\tHit\tEvalue\tIdentity\n")
            for line in f:
                items = line.strip().split("\t")
                sample = items[5].split("_")[0]
                num = int(items[5].split("_x")[1])
                sample_name = self.parse_config(file=self.option("config").prop["path"], section="NAME", name=sample)
                samples[sample_name] = sample
                if not query.has_key(items[5]):
                    query[items[5]] = 1
                    columns = [items[5], items[10], items[1], items[3]]
                    w1.write("\t".join(columns) + "\n")
                    if intron.has_key(sample_name):
                        intron[sample_name] += num
                    else:
                        intron[sample_name] = num
                else:
                    pass
            w2.write("Samples\tIntron\n")
            for key in samples:
                if not intron.has_key(key):
                    intron[key] = 0
                w2.write("{}\t{}\n".format(key, intron[key]))

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
        intron_detail = self.output_dir + "/intron_detail.xls"
        a = dict()
        with open (intron_detail, "r") as f1:
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
            "id": "IntronStat_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna.srna.intron_stat",
            "instant": False,
            "options": dict(
                blast_table="/mnt/ilustre/users/sanger-dev/workspace/20180427/Single_srna15-55-29/Srna/output/rfam_blast/blast.xls",
                config="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/quantifier_test/Uniq.cfg.ini",
                query="/mnt/ilustre/users/sanger-dev/workspace/20180426/Single_KnownMirna5123/KnownMirna/output/filtered.fa",
                gff="/mnt/ilustre/users/sanger-dev/workspace/20180502/Single_srna08-44-49/Srna/Repeatmasker/output/Trinity.SSR.gff",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()