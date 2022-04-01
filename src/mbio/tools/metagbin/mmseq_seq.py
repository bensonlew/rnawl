#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,shutil
import pandas as pd
from Bio import SeqIO
from mbio.packages.metagbin.common_function import link_dir,bin_rename


class MmseqSeqAgent(Agent):
    """
    用于mmseq2对序列去做聚类分析
    version 1.0
    author: gaohao
    last_modify: 2019.07.30
    qingchen.zhang @2019.10.11 modify
    """

    def __init__(self, parent):
        super(MmseqSeqAgent, self).__init__(parent)
        options = [
            {"name": "contig_fa", "type": "infile", "format": "sequence.fasta"},  # 输入序列文件contigs.fa
            {"name": "ref_fa", "type": "outfile", "format": "sequence.fasta"},
            {"name": "cdhit_identity", "type": "float", "default": 0.95},  # 给出cdhit的参数identity
            {"name": "cdhit_coverage", "type": "float", "default": 0.9},  # 给出cdhit的参数coverage
        ]
        self.add_option(options)
        self._memory_increase_step = 50  # 每次重运行增加内存50G by qingchen.zhang

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('contig_fa').is_set:
            raise OptionError("组装的fasta文件不存在！")

    def set_resource(self):
        """
        所需资源
        """
        self.size = os.path.getsize(self.option("contig_fa").prop["path"])
        self._cpu = 12
        num =float(self.size)/float(1000000000)
        if num < 1:
            self._memory = "50G"
        else:
            num2 = int(round(num))*20
            self._memory =str(num2) + 'G'

    def end(self):
        super(MmseqSeqAgent, self).end()


class MmseqSeqTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(MmseqSeqTool, self).__init__(config)
        self.mmseqs ="/bioinfo/metaGenomic/mmseqs/bin/mmseqs"
        self.fa =self.option('contig_fa').prop['path']
        self.ident = self.option('cdhit_identity')
        self.cover = self.option('cdhit_coverage')

    def run_mmseqs(self):
        """
        运行mmseq进行聚类
        """
        cmd = "{} easy-linclust --cluster-mode 2 --threads 10 --min-seq-id {} -c {} --cov-mode 1 {} {} {}".format(self.mmseqs, self.ident, self.cover, self.fa, self.work_dir + "/all",self.work_dir + "/temp")
        self.logger.info(cmd)
        self.logger.info("开始运行run_mmseqs")
        command = self.add_command("run_mmseqs", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_mmseqs完成")
        else:
            self.set_error("运行run_mmseqs运行出错!")

    def run_stat(self):
        """
        对聚类的结果进行整理
        """
        self.logger.info("正在对聚类结果进行统计")
        out_fasta = self.work_dir + '/all_rep_seq.fasta'
        all_cluster = os.path.join(self.work_dir, 'all_cluster.tsv')
        os.system("sed -i '1icluster_id\tquery_id' %s" %all_cluster)
        all_cluster_result = os.path.join(self.work_dir, "homologues_cluster.xls")
        cluster_list = []
        sample_list = []
        for seq_record in SeqIO.parse(out_fasta, 'fasta'):#打开第一次主要是获取到聚类的结果id和样本名称
            seq_id = seq_record.id
            sample = seq_id.split("|")[0]
            if sample not in sample_list:
                sample_list.append(sample)
            if seq_id not in cluster_list:
                cluster_list.append(seq_id)
        sample_list.sort()
        cluster_list.sort()
        data = pd.read_table(all_cluster, sep='\t', header=0)
        cluster = data.groupby('cluster_id').apply(lambda x:','.join(x['query_id']))
        cluster.to_csv('all_sample.xls', sep='\t', header=0)
        cluster.columns = ['cluster_id', 'query_name']
        new_cluster = {}
        with open('all_sample.xls', 'r') as f, open(all_cluster_result, 'w') as w:
            all_sample = "\t".join(sample_list)
            w.write("Cluster_ID\tSample_number\tGene_number\t" + all_sample + "\n")
            lines = f.readlines()
            cluster_num = 0
            for line in lines:
                line = line.strip().split("\t")
                cluster_num += 1
                cluster_id = "CLUSTER" + str(cluster_num)
                new_cluster[line[0]] = "CLUSTER" + str(cluster_num)
                query_id = line[1]
                query_id_list = query_id.split(",")
                all_sample_num =0
                all_gene_num =0
                cluster_sample_list = []
                for sample in sample_list:
                    sample_cluster_list = []
                    gene_num = 0
                    sample_num =0
                    for new_query in query_id_list:
                        if sample in new_query: #将具有相同样本名称的基因id放到一起并以","分割
                            gene_num += 1
                            sample_cluster_list.append(new_query)
                    if len(sample_cluster_list) != 0:
                        sample_name = ','.join(sample_cluster_list)
                        sample_num += 1
                    else:
                        sample_name = '-'
                    cluster_sample_list.append(sample_name)
                    all_gene_num += gene_num
                    all_sample_num += sample_num
                all_sample_id = '\t'.join(cluster_sample_list)#将不同样本名称的基因id放到一起并以"\t"分割
                w.write("{}\t{}\t{}\t{}\n".format(cluster_id, all_sample_num, all_gene_num, all_sample_id))
        self.logger.info("对聚类结果统计完成")
        return new_cluster

    def get_fasta(self, cluster):
        """
        根据OG0000892与cluster对应关系，获得聚类得到的fasta结果,改名换成cluster序列的编号
        :return:
        """
        self.logger.info("正在对代表序列进行改名")
        cluster_fasta = os.path.join(self.output_dir, "/all_cluster.faa")
        with open(cluster_fasta, 'w') as w:
            file_path = self.work_dir + '/all_rep_seq.fasta'
            for seq_record in SeqIO.parse(file_path, 'fasta'):#打开第一次主要是获取到聚类的结果id和样本名称
                seq = seq_record.seq
                seq_id = seq_record.id
                new_seq_id = cluster[seq_id]
                w.write(">{}\n{}\n".format(new_seq_id, seq))

    def set_output(self):
        """
        最后得到两个结果：聚类得到的faa文件和cluster对应的结果表
        """
        self.logger.info("正在生成结果文件目录")
        if os.path.exists(self.output_dir + '/all_cluster.faa'):
            os.remove(self.output_dir + '/all_cluster.faa')
        os.link(self.work_dir + '/all_rep_seq.fasta',self.output_dir + '/all_cluster.faa')
        if os.path.exists(self.output_dir + '/homologues_cluster.xls'):
            os.remove(self.output_dir + '/homologues_cluster.xls')
        os.link(self.work_dir + '/homologues_cluster.xls', self.output_dir + '/homologues_cluster.xls')
        self.option('ref_fa',self.output_dir + '/all_cluster.faa')

    def run(self):
        """
        运行
        """
        self.logger.info("开始运行mmseq的tool")
        super(MmseqSeqTool, self).run()
        self.run_mmseqs()
        cluster = self.run_stat()
        #self.get_fasta(cluster)
        self.set_output()
        self.end()