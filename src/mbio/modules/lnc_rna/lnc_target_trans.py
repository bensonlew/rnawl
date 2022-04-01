#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import pandas as pd
import glob
import unittest
import shutil
from biocluster.file import exists
from biocluster.file import download
from biocluster.api.file.lib.transfer import MultiFileTransfer
from Bio import SeqIO


class LncTargetTransModule(Module):
    """
    对所有样本进行定量
    """
    def __init__(self, work_id):
        super(LncTargetTransModule, self).__init__(work_id)
        options = [
            {'name': 'query', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'target', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'diff_summary', 'type': 'string', 'default': None},
            {'name': 'method', 'type': 'string'},
            {'name': 'p', 'type': 'int', "default": 30},
            {'name': 'query_list', 'type': 'string', "default": ''},
            {'name': 'query_list_type', 'type': 'string', "default": "G"},
            {'name': 'target_list', 'type': 'string', "default": ''},
            {'name': 'target_list_type', 'type': 'string', "default": "G"},
            {'type': 'infile', 'name': 'annotation', 'format': 'lnc_rna.common'},
            {'type': 'infile', 'name': 'g2t', 'format': 'lnc_rna.common'},
            {'type': 'infile', 'name': 'gtf', 'format': 'lnc_rna.gtf'},
            {'type': 'infile', 'name': 'lnc_gtf', 'format': 'lnc_rna.gtf'}
        ]
        self.add_option(options)
        self.tools = list()
        self.samples = list()

    def check_options(self):
        if self.option('gtf').is_set or self.option('g2t').is_set:
            pass
        else:
            raise OptionError("gtf g2t 需要至少设置一个")
        pass

    def run(self):
        super(LncTargetTransModule, self).run()
        self.choose_seq()
        target_dict = self.split_seq()
        self.tool_run(target_dict)

    def choose_seq(self):
        '''
        根据列表筛选mRNA 与 靶基因序列
        '''
        t2g = self.get_t2g()
        choose_list = list()

        with open(self.work_dir + ".t2g", 'w') as f:
            for t,g in t2g.items():
                f.write("{}\t{}\n".format(g, t))

        if os.path.isfile(self.option('query_list')):
            with open(self.option('query_list'), 'r') as l_in:
                choose_list = [x.strip() for x in l_in.readlines()]

            if self.option('query_list_type') == 'G':
                choose_list = self.get_translist_from_genelist(choose_list, t2g)
            self.option('query').choose_seq_by_list(choose_list, self.work_dir + "/" +  'query_choose.fa')
            self.query = self.work_dir + "/" + 'query_choose.fa'
        elif self.option("diff_summary"):
            choose_list = self.get_gene_list_summary()
            print "choose_list is {}".format(choose_list)
            if self.option('query_list_type') == 'G':
                choose_list = self.get_translist_from_genelist(choose_list, t2g)
            self.option('query').choose_seq_by_list(choose_list, self.work_dir + "/" + 'query_choose.fa')
            self.query = self.work_dir + "/" + 'query_choose.fa'
        else:
            self.query = self.option('query').prop['path']

        print "**** {}".format(self.query)

        if os.path.isfile(self.option('target_list')):
            with open(self.option('target_list'), 'r') as l_in:
                choose_list = [x.strip() for x in l_in.readlines()]
            if self.option('target_list_type') == 'G':
                choose_list = self.get_translist_from_genelist(choose_list, t2g)
            self.option('target').choose_seq_by_list(choose_list, self.work_dir + "/" + 'target_choose.fa')
            self.target = self.work_dir + "/" + 'target_choose.fa'
        elif self.option("diff_summary"):
            choose_list = self.get_gene_list_summary()
            if self.option('target_list_type') == 'G':
                choose_list = self.get_translist_from_genelist(choose_list, t2g)
            self.option('target').choose_seq_by_list(choose_list, self.work_dir + "/" + 'target_choose.fa')
            self.target = self.work_dir + "/" + 'target_choose.fa'
        else:
            self.target = self.option('target').prop['path']

    def get_translist_from_genelist(self, gene_list, t2g):
        tran_list = []
        for t,g in t2g.items():
            if g in gene_list:
                tran_list.append(t)
        return tran_list

    def split_seq(self):
        split_fasta = self.target
        file_split = dict()
        target_dict = dict()

        num = 0
        print "split_fasta is {}".format(split_fasta)
        for seq in SeqIO.parse(split_fasta, "fasta"):
            file_num = num/50
            num += 1
            if file_num in file_split:
                file_split[file_num].write('>{}\n{}\n'.format(seq.name, seq.seq))
            else:
                file_split[file_num] = open(self.work_dir + "/target.{}".format(file_num), 'w')
                file_split[file_num].write('>{}\n{}\n'.format(seq.name, seq.seq))
        for i in file_split:
            target_dict[i] = self.work_dir + "/target.{}".format(i)
            file_split[i].close()

        return target_dict


    def get_gene_list_summary(self):
        gene_list = []
        if self.option("diff_summary"):
            with open(self.option("diff_summary"), 'r') as diff_f:
                diff_f.readline()
                diff_f.readline()
                for line in diff_f:
                    cols = line.strip().split("\t")
                    gene_list.append(cols[0])
        return gene_list

    def get_t2g(self):
        m_dict = dict()
        if self.option('g2t').is_set:
            g2t = self.option('g2t').prop['path']
            with open(g2t, 'w') as g2t_f:
                for line in g2t_f.readlines():
                    cols = line.strip().split("\t")
                    m_dict[cols[1]] = cols[0]
        else:
            m_dict = self.option('gtf').get_txpt_gene_dic()
            if self.option('lnc_gtf').is_set:
                m_dict.update(self.option('lnc_gtf').get_txpt_gene_dic())
            return m_dict

    def tool_run(self, target_dict):
        self.tools = list()
        # print target_dict
        for num, target_file in target_dict.items():
            tool_opts = {
                "target": target_file,
                "query": self.option("query"),
                "method": self.option("method"),
                "g2t": self.work_dir + ".t2g",
                "lnc_gtf": self.option("lnc_gtf"),
                "annotation": self.option("annotation"),
                "diff_summary": self.option("diff_summary"),
                'query_list_type': self.option("query_list_type"),
                'target_list_type': self.option("target_list_type"),
            }
            tool = self.add_tool('lnc_rna.lnc_target_trans')
            tool.set_options(tool_opts)
            self.tools.append(tool)
        if len(self.tools) == 1:
            self.tools[0].on("end", self.set_output)
            self.tools[0].run()
        elif len(self.tools) == 0:
            self.set_output()
        else:
            self.on_rely(self.tools, self.set_output)
            for tool in self.tools:
                tool.run()

    def set_output(self):
        self.logger.info("合并结果")
        header = True
        with open(self.output_dir + "/all_merge_out.annot.xls", 'w') as f:
            for tool in self.tools:
                target_result = glob.glob(tool.output_dir + "/*_merge_out.annot.xls")
                if os.path.exists(target_result[0]):
                    target_f = open(target_result[0], 'r')
                    if header:
                        [f.write(line) for line in target_f]
                        header = False
                    else:
                        target_f.readline()
                        if self.option("method") == "RNAplex":
                            [f.write(line) for line in target_f if float(line.strip().split()[-1]) < -25]
                        else:
                            [f.write(line) for line in target_f]
                    target_f.close()

        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "表达定量分析结果目录"],
        ])
        super(LncTargetTransModule, self).end()



class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "LncTargetTrans+" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "lnc_rna.quant",
            "instant": False,
            "options": dict(
                transcriptome="/mnt/ilustre/users/sanger-dev/workspace/20190322/Single_assemble_3643_9041/Assemble/output/NewTranscripts/all_transcripts.fa",
                fastq="/mnt/ilustre/users/isanger/workspace/20190213/Single_LncrnaQc_9774/LncrnaQc/output/sickle_dir/fq_list.txt",
                method="Salmon",
                t2g="/mnt/ilustre/users/sanger-dev/workspace/20190322/Single_assemble_3643_9041/Assemble/output/NewTranscripts/trans2gene",
                pool=6,
                thread=6,
                output=None,
                read_len=149,
                read_len_sd=30,
                tpm_threshold=0,
                #map_tool="bowtie2",
            )
        }
        # data['options']['method'] = 'rsem'
        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.run()
        # #
        # data['id'] += '1'
        # data['options']['method'] = 'salmon'
        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.run()
        #
        #data['id'] += 'gdq'
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
