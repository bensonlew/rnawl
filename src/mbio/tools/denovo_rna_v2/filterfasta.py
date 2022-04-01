# coding=utf-8
import os
import glob
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
from Bio import SeqIO
# import pandas as pd
__author__ = 'liubinxu'

class FilterfastaAgent(Agent):
    """
    使用表达量过滤转录组组装结果
    """
    def __init__(self, parent):
        super(FilterfastaAgent, self).__init__(parent)
        options = [
            {'type': 'infile', 'name': 'transcripts', 'format': 'denovo_rna_v2.trinity_fasta'},
            {'type': 'string', 'name': 'matrix'},
            {'type': 'float', 'name': 'min_expr_any', 'default': 0.1},
            {'type': 'int', 'name': 'filter_200', 'default': 0},
            {'type': 'int', 'name': 'filter_500', 'default': 0},
            {'type': 'int', 'name': 'filter_1000', 'default': 0},
            {'type': 'string', 'name': 'matrix_type', 'default': 'matrix'}, # 可以为matrix 或 transrate
            {'default': 0, 'type': 'float', 'name': 'min_percent_dom_iso'},
            {'default': False, 'type': 'bool', 'name': 'highest_iso_only'},
            {'default': True, 'type': 'float', 'name': 'trinity_mode'},
            {'default': '', 'type': 'string', 'name': 'gene_to_trans_map'},
            {'type': 'outfile', 'name': 'exp_filterfasta', 'format': 'denovo_rna_v2.trinity_fasta'},
        ]
        self.add_option(options)

    def check_options(self):
        if not os.path.exists(self.option('matrix')):
            self.set_error('文件不存在：%s', variables=(self.option('matrix')), code="32003701")
        if self.option("matrix_type") == "matrix":
            if self.option("filter_200") != 0 or self.option("filter_500") != 0 or self.option("filter_1000") != 0:
                self.set_error('表达量矩阵为fpkm矩阵不能使用 count过滤', code="32003702")
        if self.option('min_expr_any') < 0 and self.option('min_percent_dom_iso') == 0  and self.option('highest_iso_only') == False:
            self.set_error('条件设置不正确 min_expr_any min_percent_dom_iso highest_iso_only至少设置一个', code="32003703")

    def set_resource(self):
        self._cpu = 1
        self._memory = "{}G".format('6')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", ""]
            ])
        super(FilterfastaAgent, self).end()


class FilterfastaTool(Tool):
    """
    使用表达量过滤转录组组装结果
    """
    def __init__(self, config):
        super(FilterfastaTool, self).__init__(config)
        self.filterfasta = '/bioinfo/denovo_rna_v2/trinityrnaseq-Trinity-v2.5.0/util/filter_low_expr_transcripts.pl'

    def run_filterfasta(self):
        matrix = self.option("matrix")
        if self.option('matrix_type') == 'transrate':
            table = pd.read_table(matrix, header=0)
            table_choose = table[['Name','TPM']]
            table_choose.to_csv('choose_exp.xls',sep="\t",index=False)
            matrix = 'choose_exp.xls'

        cmd = '{} '.format(self.filterfasta)
        cmd += '--{} {} '.format("transcripts", self.option("transcripts").prop['path'])
        cmd += '--{} {} '.format("matrix", matrix)
        cmd += '--{} {} '.format("min_expr_any", self.option("min_expr_any"))
        if self.option("highest_iso_only"):
            cmd += '--{} '.format("highest_iso_only")
        elif self.option("min_percent_dom_iso") > 0:
            cmd += '--{} {} '.format("min_percent_dom_iso", self.option("min_percent_dom_iso"))
        else:
            pass

        if self.option("trinity_mode"):
            cmd += '--{} '.format("trinity_mode")
        elif self.option("gene_to_trans_map"):
            cmd += '--{} {} '.format("gene_to_trans_map", self.option("gene_to_trans_map"))
        else:
            pass

        cmd_name = 'filterfasta'
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
                self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="32003704")
        else:
            self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="32003705")

    def run_filterfasta_bycount(self):
        '''
        根据count标准过滤fasta
        '''
        count_list = []
        matrix = self.option("matrix")
        table = pd.read_table(matrix, header=0)

        count_list.extend(list(table[(table['Length']<=500) & (table['NumReads']>=self.option("filter_200"))]['Name']))
        count_list.extend(list(table[(table['Length']>500) & (table['Length']<1000) & (table['NumReads']>=self.option("filter_500"))]['Name']))
        count_list.extend(list(table[(table['Length']>1000) & (table['NumReads']>=self.option("filter_1000"))]['Name']))

        count_set = set(count_list)

        with open("exp.clean2.fa", 'w') as clean2_w:
            for seq in SeqIO.parse('exp.clean.fa', "fasta"):
                seq_seq = seq.seq
                seq_name = seq.name
                if seq_name in count_set:
                    clean2_w.write('>{}\n{}\n'.format(seq_name, seq_seq))

    def set_output(self):
        '''
        设置输出
        '''
        clean_files = 'exp.clean.fa'
        fname = os.path.basename(clean_files)
        if self.option("filter_200") != 0 or self.option("filter_500") != 0 or self.option("filter_1000") != 0:
            clean_files = 'exp.clean2.fa'
        link = os.path.join(self.output_dir, fname)
        if os.path.exists(link):
            os.remove(link)
        os.link(clean_files, link)
        self.option('exp_filterfasta', link)

    def run(self):
        super(FilterfastaTool, self).run()
        self.run_filterfasta()
        if self.option("filter_200") != 0 or self.option("filter_500") != 0 or self.option("filter_1000") != 0:
            self.run_filterfasta_bycount()
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
        test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_denovo/test_data2'
        data = {
            "id": "FilterFasta" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.filterfasta",
            "instant": True,
            "options": dict(
                transcripts=test_dir + "/" + "Trinity.fasta",
                matrix=test_dir + "/" + "Trinity.fasta_quant.sf",
                min_expr_any="0.1",
                matrix_type="transrate",
                min_percent_dom_iso="0",
                highest_iso_only=False,
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
