# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

"""Usearch聚类工具"""
from __future__ import division
import math
from biocluster.tool import Tool
from biocluster.agent import Agent
from Bio import SeqIO
import os
import re
from biocluster.core.exceptions import OptionError


class PanUsearchAgent(Agent):
    """
    Usearch：uparse
    version v7
    author：qingchen.zhang
    细菌比较基因组开发使用, 改写自多样性聚类的tool--Usearch
    """
    def __init__(self, parent=None):
        super(PanUsearchAgent, self).__init__(parent)
        options = [
            {'name': 'fasta', 'type': 'infile', 'format': 'sequence.fasta'},# 输入fasta文件，序列名称格式为'>sampleID_seqID'.
            {'name': 'identity', 'type': 'float', 'default': 0.5},# 相似性值，范围0-1.
            {'name': 'otu_table', 'type': 'outfile','format': 'meta.otu.otu_table'},  # 输出结果otu表
        ]
        self.add_option(options)
        self.step.add_steps('OTUCluster')
        self._memory_increase_step = 50  # 每次重运行增加内存50G by qingchen.zhang
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.OTUCluster.start()
        self.step.update()

    def step_end(self):
        self.step.OTUCluster.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数设置
        """
        if not self.option("fasta").is_set:
            raise OptionError("必须设置输入fasta文件.")
        if self.option("identity") < 0 or self.option("identity") > 1:
            raise OptionError("identity值必须在0-1范围内.")
        return True

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["otu_reps.fasta", "sequence.fasta", "代表序列"],
            ["otu_table.xls", "meta.otu.otu_table", "OTU表"]
        ])
        super(PanUsearchAgent, self).end()

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        total = os.path.getsize(self.option("fasta").prop["path"])
        total = int(math.ceil(total / (1024 * 1024 * 1024)))
        total = int(total * 10)
        self._memory = "{}G".format(total)


class PanUsearchTool(Tool):
    """
    UsearchOTU tool
    """
    def __init__(self, config):

        super(PanUsearchTool, self).__init__(config)
        self._version = "v7.0"
        self.usearch_path = "bioinfo/meta/usearch-v7.0/"
        self.script_path = "bioinfo/meta/scripts/"
        self.qiime_path = "program/Python/bin/"

    def cmd1(self):
        """
        加索引size大小
        :return:
        """
        cmd = self.usearch_path + \
            "uparse -derep_prefix meta.fasta -output meta_derepprefix.fasta -sizeout"
        return cmd

    def cmd2(self):
        """
        排序
        :return:
        """
        cmd = self.usearch_path + \
            "uparse -sortbysize meta_derepprefix.fasta -output meta_derepprefix_sorted.fasta -minsize 2"
        return cmd

    def cmd3(self):
        """
        make db序列
        :return:
        """
        ratio = str(100 - float(self.option('identity')) * 100)
        cmd = self.usearch_path + \
            "uparse -cluster_otus meta_derepprefix_sorted.fasta -otus cluster.fasta -otu_radius_pct " + ratio
        return cmd

    def cmd4(self):
        """
        开始聚类参照db
        :return:
        """
        cmd = self.usearch_path + "uparse -usearch_global meta.fasta -db cluster.fasta -strand plus -id " + \
            str(self.option('identity')) + " -uc map.uc"
        return cmd

    def cmd5(self):
        """
        整理出聚类的id结果
        :return:
        """
        cmd = self.script_path + """uc2otuseqids.pl -i map.uc -o cluster.seqids"""
        return cmd

    def run_stat(self):
        """
        对运行的结果得到的cluster情况进行整理和统计
        :return:
        """
        self.logger.info("开始对聚类结果进行统计")
        cluster_result = os.path.join(self.work_dir, 'cluster.seqids')
        cluster_fasta = os.path.join(self.work_dir, 'cluster.fasta')
        sample_list = []
        for seq_record in SeqIO.parse(cluster_fasta, 'fasta'):#打开第一次主要是获取到聚类的结果id和样本名称
            seq_id = seq_record.id
            sample = seq_id.split("|")[0]
            if sample not in sample_list:
                sample_list.append(sample)
        sample_list.sort()
        all_cluster_result = os.path.join(self.work_dir, "homologues_cluster.xls")
        with open(cluster_result, 'r') as f, open(all_cluster_result, 'w') as w:
            all_sample = "\t".join(sample_list)
            w.write("Cluster_ID\tSample_number\tGene_number\t" + all_sample + "\n")
            lines = f.readlines()
            cluster_num = 0
            for line in lines:
                line = line.strip().split("\t")
                length = len(line)
                cluster_num += 1
                cluster_id = "CLUSTER" + str(cluster_num)
                all_sample_num =0
                all_gene_num =0
                cluster_sample_list = []
                for sample in sample_list:
                    sample_cluster_list = []
                    gene_num = 0
                    sample_num =0
                    for i in range(1, length):
                        if sample in line[i]: #将具有相同样本名称的基因id放到一起并以","分割
                            gene_num += 1
                            sample_cluster_list.append(line[i])
                    if len(sample_cluster_list) != 0:
                        sample_name = ','.join(sample_cluster_list)
                        sample_num += 1
                    else:
                        sample_name = '-'
                    cluster_sample_list.append(sample_name)
                    all_gene_num += gene_num
                    all_sample_num += sample_num
                all_sample_id = '\t'.join(cluster_sample_list)#将不同样本名称的基因id放到一起并以"\t"分割
                w.write("{}\t{}\t{}\t{}\n".format(cluster_id, str(all_sample_num), str(all_gene_num), all_sample_id))
        self.logger.info("对聚类结果统计完成")

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("设置输出结果")
        result_path = os.path.join(self.output_dir, "homologues_cluster.xls")
        if os.path.exists(result_path):
            os.remove(result_path)
        os.link(os.path.join(self.work_dir, "homologues_cluster.xls"), result_path)

    def run(self):
        """
        运行
        :return:
        """
        super(PanUsearchTool, self).run()
        self.logger.info("将输入文件链接到工作目录")
        if os.path.exists(self.work_dir + '/meta.fasta'):
            os.remove(self.work_dir + '/meta.fasta')
        os.link(self.option("fasta").prop['path'], self.work_dir + '/meta.fasta')
        self.logger.info("OK")
        i = 0
        while i < 5:
            i += 1
            self.logger.info("开始运行cmd" + str(i))
            cmd = getattr(self, 'cmd' + str(i))()
            self.logger.info(cmd)
            command = self.add_command('cmd' + str(i), cmd)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("运行cmd" + str(i) + "完成")
            else:
                if command.return_code is None:
                    command.rerun()
                    self.wait()
                    if command.return_code == 0:
                        self.logger.info("重新运行cmd" + str(i) + "完成")
                    elif command.return_code is None:
                        self.logger.info("cmd" + str(i) + "重新运行返回码仍然为None")
                    else:
                        self.logger.info('Run Return Code: {}'.format(command.return_code))
                        self.set_error("cmd %s 运行出错!")
                        break
                else:
                    self.logger.info('Run Return Code: {}'.format(command.return_code))
                    self.set_error("cmd %s 运行出错!")
                    break
        self.run_stat()
        self.set_output()
        self.end()
