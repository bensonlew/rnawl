# coding=utf-8
import os
import re
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
from mbio.packages.denovo_rna_v2.trans_step import step_count
import pandas as pd
__author__ = 'liubinxu'

class TrinityStatAgent(Agent):
    """
    Trinity
    """
    def __init__(self, parent):
        super(TrinityStatAgent, self).__init__(parent)
        options = [
            {'type': 'infile', 'name': 'fasta', 'format': 'denovo_rna_v2.trinity_fasta'},
            {'type': 'string', 'name': 'gene2trans', 'default': ''},
            {'type': 'int', 'name': 'min_len', 'default': 200},
            {'type': 'string', 'name': 'busco_result', 'default': ''},
            {'type': 'string', 'name': 'transrate_result', 'default': ''},
            {'type': 'string', 'name': 'exp_result', 'default': ''}, #transrate 表达统计结果
            {'type': 'outfile', 'name': 'out_dir', "format": "denovo_rna_v2.common"},
            {'type': 'bool', 'name': 'marker', 'default': False},
            #{'type': 'outfile', 'name': 'out_ungene_stat'},
            #{'type': 'outfile', 'name': 'out_unigene', 'format': 'sequence.fasta'},
            #{'type': 'outfile', 'name': 'out_unigene_length'},
            #{'type': 'outfile', 'name': 'out_trans_length'},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fasta").is_set:
            raise OptionError("必须设置参数fasta", code = "32006601")
        if not self.option("gene2trans"):
            raise OptionError("必须设置gene2trans参数", code = "32006602")
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = "{}G".format('5')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ["Trinity_stat.xls", "csv", "Trinity 统计结果"],
            [r"count_stat*.txt", "txt", "组装结果长度分布"],
            ["unigene.fasta", "fasta", "unigene 序列文件"],
            [r"Trinity_t2*", "txt", "转录本标记文件"]
            ])
        super(TrinityStatAgent, self).end()

class TrinityStatTool(Tool):
    """
    Trinity
    """
    def __init__(self, config):
        super(TrinityStatTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.trinity_stat = 'bioinfo/denovo_rna_v2/trinityrnaseq-Trinity-v2.5.0/util/TrinityStats.pl'
        self.e90n50 = 'bioinfo/denovo_rna_v2/trinityrnaseq-Trinity-v2.5.0/util/misc/contig_ExN50_statistic.pl'
        self.stat_data = {}

    def run_trinity_stat(self):
        '''
        统计Trinity结果信息
        '''
        cmd = '{} '.format(self.trinity_stat)
        cmd += '{} '.format(self.option("fasta").prop['path'])
        cmd_name = 'trinitystat'
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
                self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "32006603")
        else:
            self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "32006604")

        # 读入Trinity的统计结果
        with open("Trinity.stat.txt", 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split(':')
                self.stat_data[line[0]] = line[1]

    def run_exp_stat(self):
        '''
        表达量计算
        '''
        matrix = self.option('exp_result')
        table = pd.read_table(matrix, header=0)
        mean_read = table[['NumReads']].mean()['NumReads']
        table_choose = table[['Name','TPM']]
        table_choose.to_csv('choose_exp.xls',sep="\t",index=False)
        self.stat_data['Mean mapped reads'] = mean_read


    def run_e90_n50_stat(self):
        '''
        统计E90N50信息
        '''
        cmd = '{} {} {}'.format(self.e90n50, 'choose_exp.xls', self.option("fasta").prop['path'])
        cmd_name = 'exn50_stat'
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
                self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "32006605")
        else:
            self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "32006606")

        table=pd.read_table('exn50.txt', header=0)
        # e_num = 90
        for e_num in range(90,100):
            try:
                e90n50 = table.loc[table['E']==e_num].iloc[0,1]
                break
            except:
                self.logger.info("E90n50 not exists")

        self.stat_data['E90N50'] = e90n50

    def run_fasta_stat(self):
        '''
        统计fasta文件基本信息
        '''
        form,seq_type,seq_number,bases,longest,shortest = self.option("fasta").get_seq_info()
        self.stat_data['longest'] = longest
        self.stat_data['shortest'] = shortest

    def trinity_fasta_marker(self):
        '''
        获取unigene序列，和相关对应关系列表
        '''
        self.option("fasta").set_gene2tran(self.option("gene2trans"))
        self.option("fasta").get_unigene('unigene.fasta', 'Trinity')
        unigene_list = list()
        with open("unigene.fasta", 'r') as f:
            [unigene_list.append(i) for i in f if i.startswith(">")]
        self.stat_data['Total trinity genes'] = str(len(unigene_list))

    def get_unigene(self):
        '''
        获取unigene序列，和相关对应关系列表
        '''
        #self.option("fasta").set_gene2tran(self.option("gene2trans"))
        self.option("fasta").get_unigene_by_marker('unigene.fasta', self.option("gene2trans"))
        unigene_list = list()
        with open("unigene.fasta", 'r') as f:
            [unigene_list.append(i) for i in f if i.startswith(">")]
        self.stat_data['Total trinity genes'] = str(len(unigene_list))

    def buscle_stat(self):
        '''
        获取buscle的统计结果
        '''
        busco_result = self.option("busco_result")
        with open(busco_result, 'r') as f:
            lines = f.readlines()
            stat_line = lines[7]
            '''
            stat_match = re.search(r'C:(.*)\[S:(.*),D:(.*)\],F:(.*),M:(.*),n:(.*)', stat_line)
            stat_c = stat_match.group(1)
            stat_d = stat_match.group(3)
            self.stat_data['Busco score'] = stat_c + '(' + stat_d + ')'
            '''
            self.stat_data['Busco score'] = stat_line.strip()

    def transrate_stat(self):
        '''
        获取transrate的统计结果
        '''
        transrate_result = self.option("transrate_result")
        with open(transrate_result, 'r') as f:
            lines = f.readlines()
            #head_line = lines[0]
            stat_line = lines[1]
            stat_line = stat_line.strip().split(',')
            self.stat_data['TransRate score'] = stat_line[-4]
            self.stat_data['Mean mapped percent'] = float(stat_line[21] ) * 100


    def length_stat(self):
        '''
        长度分布统计
        '''
        fasta = self.option("fasta").prop['path']
        steps = [500]
        for step in steps:
            files = "trans_length.txt"
            step_count(fasta, files , 10, step, self.work_dir + "/trans_count_stat_" + str(step) + ".txt", self.option("min_len"))
            files = 'unigene_length.txt'
            step_count('unigene.fasta', files , 10, step, self.work_dir + "/unigene_count_stat_" + str(step) + ".txt", self.option("min_len"))

    def save_stat(self):
        '''
        保存统计结果
        '''
        with open('Trinity_stat.xls', 'w') as trinity_stat:
            '''
            trinity_stat.write('Type\tResource\n')
            trinity_stat.write('Total transcripts num\t{}\n'.format(self.stat_data['Total trinity transcripts']))
            trinity_stat.write('Total unigenes num\t{}\n'.format(self.stat_data['Total trinity genes']))
            trinity_stat.write('Total sequence base\t{}\n'.format(self.stat_data['Total assembled bases']))
            trinity_stat.write('Largest\t{}\n'.format(self.stat_data['longest']))
            trinity_stat.write('Smallest\t{}\n'.format(self.stat_data['shortest']))
            trinity_stat.write('Average length\t{}\n'.format(self.stat_data['Average contig']))
            trinity_stat.write('N50\t{}\n'.format(self.stat_data['Contig N50']))
            if 'E90N50' in self.stat_data.keys():
                trinity_stat.write('E90N50\t{}\n'.format(self.stat_data['E90N50']))
            else:
                trinity_stat.write('E90N50\t{}\n'.format(0))
            trinity_stat.write('GC percent\t{}\n'.format(self.stat_data['Percent GC']))
            if 'Mean mapped reads' in self.stat_data.keys():
                trinity_stat.write('Mean mapped reads\t{}\n'.format(self.stat_data['Mean mapped reads']))
            else:
                trinity_stat.write('Mean mapped reads\t{}\n'.format(0))
            trinity_stat.write('TransRate score\t{}\n'.format(self.stat_data['TransRate score']))
            trinity_stat.write('BUSCO score\t{}\n'.format(self.stat_data['Busco score']))
            '''
            trinity_stat.write('Type\tResource\n')
            trinity_stat.write('Total number\t{}\n'.format(self.stat_data['Total trinity transcripts']))
            # trinity_stat.write('Total unigenes num\t{}\n'.format(self.stat_data['Total trinity genes']))
            trinity_stat.write('Total base\t{}\n'.format(self.stat_data['Total assembled bases']))
            trinity_stat.write('Largest length (bp)\t{}\n'.format(self.stat_data['longest']))
            trinity_stat.write('Smallest length (bp)\t{}\n'.format(self.stat_data['shortest']))
            trinity_stat.write('Average length (bp)\t{}\n'.format(self.stat_data['Average contig']))
            trinity_stat.write('N50 length (bp)\t{}\n'.format(self.stat_data['Contig N50']))
            if 'E90N50' in self.stat_data.keys():
                trinity_stat.write('E90N50 length (bp)\t{}\n'.format(self.stat_data['E90N50']))
            else:
                trinity_stat.write('E90N50 length (bp)\t{}\n'.format(0))
            trinity_stat.write('Mean mapped percent (%)\t{}\n'.format(self.stat_data['Mean mapped percent']))
            trinity_stat.write('GC percent (%)\t{}\n'.format(self.stat_data['Percent GC']))
            trinity_stat.write('TransRate score\t{}\n'.format(self.stat_data['TransRate score']))
            trinity_stat.write('BUSCO score\t{}\n'.format(self.stat_data['Busco score']))


    def set_output(self):
        trinity_stat_files = glob.glob('./Trinity_stat.xls')
        length_stat_files = glob.glob('./*count_stat*.txt')
        if self.option('marker'):
            unigene = glob.glob('./unigene.fasta')
            marker = glob.glob('./Trinity_t2*')
            all_files = trinity_stat_files + length_stat_files + unigene + marker
        else:
            all_files = trinity_stat_files + length_stat_files

        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

        #self.option('out_unigene', os.path.join(self.output_dir, 'unigene.fasta'))
        #self.option('outfile', os.path.join(self.output_dir, 'good.Trinity.fasta'))
        #self.option('bad_fa', os.path.join(self.output_dir, 'bad.Trinity.fasta'))
        #self.option('result', os.path.join(self.output_dir, 'assemblies.csv'))


    def run(self):
        super(TrinityStatTool, self).run()
        self.run_trinity_stat()
        self.run_exp_stat()
        self.run_e90_n50_stat()
        self.run_fasta_stat()
        if self.option('marker'):
            self.trinity_fasta_marker()
        else:
            self.get_unigene()
        self.buscle_stat()
        self.transrate_stat()
        self.length_stat()
        self.save_stat()
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
            "id": "TrinityStat" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.trinity_stat",
            "instant": True,
            "options": dict(
                fasta = test_dir + "/" + "Trinity.fasta",
                busco_result = test_dir + "/" + "short_summary_busco_result.txt",
                gene2trans = test_dir + "/" + "assemble/Trinity.fasta_t2g2u2l",
                transrate_result = test_dir + "/" + "assemblies.csv",
                exp_result = test_dir + "/" + "Trinity.fasta_quant.sf",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
