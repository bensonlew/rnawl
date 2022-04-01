# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
from Bio import SeqIO
import shutil
# import pandas as pd
__author__ = 'shicaiping'


class SplitRefFastaAgent(Agent):
    """
    split_fasta description
    """
    def __init__(self, parent):
        super(SplitRefFastaAgent, self).__init__(parent)
        options = [
            {'name': 'fasta', 'type': 'infile', 'format': 'ref_rna_v2.common'},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 2
        self._memory = "{}G".format('15')

    def end(self):
        super(SplitRefFastaAgent, self).end()


class SplitRefFastaTool(Tool):
    """
    split_fasta description
    """
    def __init__(self, config):
        super(SplitRefFastaTool, self).__init__(config)

    def split_ref_run(self):
        split_file = "{}/split_file".format(self.output_dir)
        if os.path.isdir(split_file):
            shutil.rmtree(split_file)
        os.mkdir(split_file)

        def write_buffer(filename, data):
            with open(filename, 'wb') as target:
                for d in data:
                    target.write('{}\t{}\t{}\n'.format(d.id, 0, len(d.seq)))
        buffers = list()
        i = 0
        file_size = os.path.getsize(self.option("fasta").prop['path'])/float(1024*1024*1024)
        file_size = round(file_size, 2)
        # 按照大约50个tool投递
        line_limit = 800000*file_size
        for seq_record in SeqIO.parse(self.option("fasta").prop['path'], 'fasta'):
            i += 1
            filename = split_file + '/fasta_{}'.format(i)
            if len(seq_record.seq) >= line_limit:
                write_buffer(filename, [seq_record])
                continue
            buffers.append(seq_record)
            if sum(len(i.seq) for i in buffers) >= line_limit:
                write_buffer(filename, buffers)
                del buffers[:]
        filename = split_file + '/fasta_{}'.format(i+1)
        if buffers:
            write_buffer(filename, buffers)
            del buffers[:]

    def run(self):
        super(SplitRefFastaTool, self).run()
        self.split_ref_run()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "SplitFasta" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v2.split_ref_fasta",
            "instant": False,
            "options": dict(
                fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/GRCm38_Ensembl_96/dna/Mus_musculus.GRCm38.dna.toplevel.fa",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
