# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from Bio import SeqIO
import unittest


class SplitFastaBaseAgent(Agent):
    """
    SplitFasta:将fasta文件按行数拆分
    version 1.0
    author: qiuping
    last_modify: 2016.11.15
    """

    def __init__(self, parent):
        super(SplitFastaBaseAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "small_rna.fasta"},
            {"name": "bases", "type": "int", "default": 10000000},  # 序列数
        ]
        self.add_option(options)
        self.step.add_steps('splitfasta')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.splitfasta.start()
        self.step.update()

    def step_end(self):
        self.step.splitfasta.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("fasta").is_set:
            raise OptionError("请传入fasta序列文件")
        if not isinstance(self.option('bases'), int):
            raise OptionError("碱基数必须为整数")
        if self.option('bases') <= 0:
            raise OptionError("碱基数小于等于0，请重设！")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '5G'


class SplitFastaBaseTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(SplitFastaBaseTool, self).__init__(config)

    def split_fasta(self):
        """
        """
        bases = 0
        i = 1
        w = open(self.output_dir + '/fasta_1', 'wb')
        for seq_record in SeqIO.parse(self.option('fasta').prop['path'], "fasta"):
            if bases <= self.option('bases'):
                w.write('>{}\n{}\n'.format(seq_record.id, seq_record.seq))
                bases += len(seq_record.seq)
            else:
                i += 1
                w.close()
                bases = 0
                w = open(self.output_dir + '/fasta_%s' % i, 'wb')
                w.write('>{}\n{}\n'.format(seq_record.id, seq_record.seq))
                bases += len(seq_record.seq)
        w.close()

    def run(self):
        super(SplitFastaBaseTool, self).run()
        self.split_fasta()
        self.end()

class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "SplitFasta_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna.srna.split_fasta_base",
            "instant": False,
            "options": dict(
                fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/dna/Mus_musculus.GRCm38.dna_rm.toplevel.clean.fa",
                bases=10000000
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()