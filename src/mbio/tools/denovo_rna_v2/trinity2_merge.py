# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
from Bio import SeqIO
__author__ = 'liubinxu'


class Trinity2MergeAgent(Agent):
    """
    Trinity
    """
    def __init__(self, parent):
        super(Trinity2MergeAgent, self).__init__(parent)
        options = [
            {'type': 'infile', 'name': 'partion_dir', 'format': 'denovo_rna_v2.common_dir'},
            {"name": "trinity_fa", "type": "outfile", "format": "denovo_rna_v2.trinity_fasta"},
            {"name": "gene2trans", "type": "outfile", "format": "denovo_rna_v2.common"}, #组装后gene和转录本的对应列表
            {"name": "min_len", "type": "int", "default": 200}
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", ""]
            ])
        super(Trinity2MergeAgent, self).end()


class Trinity2MergeTool(Tool):
    """
    Trinity
    """
    def __init__(self, config):
        super(Trinity2MergeTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.merge = '/bioinfo/denovo_rna_v2/trinityrnaseq-Trinity-v2.5.0/Trinity2_merge.pl'
        self.perl = self.config.SOFTWARE_DIR + '/program/perl-5.24.0/bin/'
        self.bowtie = self.config.SOFTWARE_DIR + '/bioinfo/denovo_rna_v2/bowtie2-2.3.3.1-linux-x86_64/'
        self.java = self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin'
        self.set_environ(PATH=self.perl)
        self.set_environ(PATH=self.bowtie)
        self.set_environ(PATH=self.java)

    def run_merge(self):
        cmd = '{} {}'.format(self.merge, self.option("partion_dir").prop['path'])
        cmd_name = 'trinity2merge'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} Finished successfully".format(cmd_name))
            else:
                self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "32006301")
        else:
            self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "32006302")

    def rename(self):
        trinity_fa = self.work_dir + '/Trinity.fasta'
        gene2trans = self.work_dir + '/Trinity.fasta.gene_trans_map'

        with open(self.output_dir + '/Trinity.fasta', 'w') as fa_fo, open(self.output_dir + '/Trinity.fasta.gene_trans_map', 'w') as g2t_fo:
            for seq in SeqIO.parse(trinity_fa, "fasta"):
                new_name = seq.id
                if len(seq.seq) > self.option("min_len"):
                    fa_fo.write(">{}\n{}\n".format(new_name, seq.seq))
                    g2t_fo.write("{}\t{}\n".format(new_name.split("_i")[0], new_name))


    def set_output(self):
        '''
        trinity_fa = self.work_dir + '/Trinity.fasta'
        gene2trans = self.work_dir + '/Trinity.fasta.gene_trans_map'
        if os.path.exists(self.output_dir + '/Trinity.fasta'):
            os.remove(self.output_dir + '/Trinity.fasta')
        if os.path.exists(self.output_dir + '/Trinity.fasta.gene_trans_map'):
            os.remove(self.output_dir + '/Trinity.fasta.gene_trans_map')
        os.link(trinity_fa, os.path.join(self.output_dir, 'Trinity.fasta'))
        os.link(gene2trans, os.path.join(self.output_dir, 'Trinity.fasta.gene_trans_map'))
        '''
        self.rename()
        self.logger.info("设置结果目录")

        try:
            self.option('trinity_fa', os.path.join(self.output_dir, 'Trinity.fasta'))
            self.option('gene2trans', os.path.join(self.output_dir, 'Trinity.fasta.gene_trans_map'))
            self.logger.info("设置组装拼接分析结果目录成功")
            self.logger.info("序列文件为 {}".format(self.option('trinity_fa').prop['path']))
            self.logger.info("设置组装拼接分析结果目录成功 {}".format(self.option('gene2trans')))
            self.end()
        except Exception as e:
            self.logger.info("设置组装拼接分析结果目录失败{}".format(e))
            self.set_error("设置组装拼接分析结果目录失败%s", variables = (e), code = "32006303")

    def run(self):
        super(Trinity2MergeTool, self).run()
        self.run_merge()
        self.set_output()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/workspace/20180102/Single_de_tr_data2.2/Trinity2/trinity_out_dir/read_partitions'
        data = {
            "id": "Trinity2Merge" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.trinity2_merge",
            "instant": False,
            "options": dict(
                partion_dir='/mnt/ilustre/users/sanger-dev/workspace/20180102/Single_de_tr_data2.2/Trinity2/trinity_out_dir/read_partitions',
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
