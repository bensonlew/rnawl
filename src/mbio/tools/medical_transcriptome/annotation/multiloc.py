# coding=utf-8
import os
import glob
from collections import Counter
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from Bio import SeqIO
import unittest

# import pandas as pd
__author__ = 'liubinxu'


class MultilocAgent(Agent):
    """
    sub cell location result
    """
    def __init__(self, parent):
        super(MultilocAgent, self).__init__(parent)
        options = [
            {'type': 'outfile', 'name': 'output', 'format': 'itraq_and_tmt.common'},
            {'type': 'string', 'name': 'species', 'default': 'Animals'},
            {'type': 'infile', 'name': 'fa', 'format': 'itraq_and_tmt.common'},
            {'type': 'string', 'name': 'go', 'default': None},
            {'type': 'int', 'name': 'cpu', 'default': 30},
            {'type': 'int', 'name': 'lines' , 'default': 100},
        ]
        self.add_option(options)

    def check_options(self):
        if self.option('species') not in ['Animals','Plants','Fungi']:
            self.option('species','Fungi')
            # raise OptionError('species必须是Fungi Animals Plants中的一种')
        pass

    def set_resource(self):
        self._cpu = self.option('cpu')
        self._memory = "{}G".format('40')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", ""]
            ])
        """
        # more detail
        result_dir.add_regexp_rules([
            [r"*.xls", "xls", "xxx"],
            [r"*.list", "", "xxx"],
            ])
        """
        super(MultilocAgent, self).end()


class MultilocTool(Tool):
    """
    sub cell location result
    """
    def __init__(self, config):
        super(MultilocTool, self).__init__(config)
        self.db = {
            'Fungi': 'fungal',
            'Animals': 'animal',
            'Plants': 'plant',
        }
        self.parafly = '/bioinfo/denovo_rna_v2/trinityrnaseq-Trinity-v2.5.0/trinity-plugins/ParaFly-0.1.0/bin/ParaFly'
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = self.config.SOFTWARE_DIR +  '/miniconda2/bin/python'
        self.multiloc2_prediction = software_dir + '/bioinfo/itraq_and_tmt/MultiLoc2-26-10-2009/src/multiloc2_prediction.py'
        self.blast_path = software_dir + "/bioinfo/itraq_and_tmt/blast-2.2.23/bin/"
        self.svm_path = software_dir + "/bioinfo/itraq_and_tmt/libsvm-3.22/"
        self.set_environ(PATH=self.blast_path)
        self.set_environ(PATH=self.svm_path)
        self.split_fasta = []
        self.parrel_method = []
        self.split_name = dict()

    def run_split_fasta(self):
        '''
        分割fasta序列
        '''
        line = 1
        i = 1
        w = open(self.work_dir + '/fasta_%s' %i, 'wb')
        self.split_fasta.append(self.work_dir + '/fasta_1')
        for seq_record in SeqIO.parse(self.option('fa').prop['path'], "fasta"):
            if line <= self.option('lines'):
                w.write('>{}\n{}\n'.format(seq_record.id, seq_record.seq))
                line += 1
            else:
                i += 1
                w.close()
                line = 1
                w = open(self.work_dir + '/fasta_%s' % i, 'wb')
                w.write('>{}\n{}\n'.format(seq_record.id, seq_record.seq))
                self.split_fasta.append(self.work_dir + '/fasta_%s' %i)
                line += 1

            if self.split_name.has_key(i):
                self.split_name[i].append(seq_record.id)
            else:
                self.split_name[i] = [seq_record.id]
        w.close()

    def run_split_go(self):
        '''
        分割go注释文件
        '''
        go_dict = dict()
        with open(self.option('go'), 'r') as go:
            for line in go.readlines():
                go_dict[line.strip().split("\t")[0]] = line.strip().split("\t")[1]
        for file_key in self.split_name.keys():
            genes = self.split_name[file_key]
            with open(self.work_dir + '/go_%s' % file_key, 'wb') as go_split:
                for gene in genes:
                    if go_dict.has_key(gene):
                        go_split.write("{}\t{}\n".format(gene, go_dict[gene]))

    def run_multiloc2_parrel(self):
        '''
        并行运行
        '''
        cmd = '{} '.format(self.parafly)
        cmd += '-{} {} '.format("c", self.work_dir + "/parrel_cmd")
        cmd += '-{} {} '.format("CPU", self.option("cpu"))
        cmd += '-v -shuffle'
        cmd_name = 'multiloc'
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
                self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))
        else:
            self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))

        cmd_cat = 'cat multiloc_result_* >> multiloc.txt'
        os.system(cmd_cat)

    def run_multiloc2_prediction(self):
        with open('parrel_cmd', 'wb') as f:
            for file_key in self.split_name.keys():
                cmd = '{} {} '.format(self.python_path, self.multiloc2_prediction)
                cmd += '-{}={} '.format("origin", self.db[self.option("species")])
                # cmd += '-{}={} '.format("output", "advanced")
                cmd += '-{}={} '.format("fasta", self.work_dir + '/fasta_%s' %file_key)
                cmd += '-{}={} '.format("go", self.work_dir + '/go_%s' %file_key)
                cmd += '-{}={} '.format("result", self.work_dir + '/multiloc_result_%s' %file_key)
                cmd += '-predictor=HighRes'
                f.write("{}\n".format(cmd))

    def set_output(self):
        first_loc = []
        with open('multiloc.txt', 'r') as multiloc, \
             open('multiloc.xls', 'w') as multiloc_result:
            lines = multiloc.readlines()
            for line in lines[5:]:
                line = line.strip().split('\t')
                if len(line) > 8:
                    multiloc_result.write('\t'.join(line[0:4]) + '\n')
                    first_loc.append(line[1].split(':')[0])
        loc_stat = Counter(first_loc)
        stat_list = sorted(loc_stat.items(), key=lambda ot:ot[1], reverse=True)
        with open('multiloc_stat.xls', 'w') as multiloc_stat:
            multiloc_stat.write('Subcelluar location\tProtein num' + '\n')
            for loc in stat_list:
                multiloc_stat.write(loc[0]+'\t' +str(loc[1]) + '\n')

        fname = 'multiloc.xls'
        link = os.path.join(self.output_dir, fname)
        if os.path.exists(link):
            os.remove(link)
        os.link('multiloc.xls', link)
        self.option('output', link)
        fname = 'multiloc_stat.xls'
        link = os.path.join(self.output_dir, fname)
        if os.path.exists(link):
            os.remove(link)
        os.link('multiloc_stat.xls', link)

    def run(self):
        super(MultilocTool, self).run()
        self.run_split_fasta()
        self.run_split_go()
        self.run_multiloc2_prediction()
        self.run_multiloc2_parrel()
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
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_itraq_and_tmt/'
        data = {
            "id": "Multiloc" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "itraq_and_tmt.annotation.multiloc",
            "instant": True,
            "options": dict(
                fa=test_dir + "/" + "data1.fa",
                go=test_dir + "/" + "data1.GO.list",
                species="Plants",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
