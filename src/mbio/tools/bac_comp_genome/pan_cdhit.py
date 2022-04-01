# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import re
from Bio import SeqIO
import subprocess
from biocluster.core.exceptions import OptionError
import shutil


class PanCdhitAgent(Agent):
    """
    细菌比较基因组用于聚类，其基本思路没有变，主要是增加了统计和修改了结果文件的名称
    """

    def __init__(self, parent):
        super(PanCdhitAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入fasta文件
            {"name": "qunum", "type": "int", "default": 0},  # fasta编号
            {"name": "identity", "type": "float", "default": 0.95},  ##给出cdhit的参数identity
            {"name": "coverage", "type": "float", "default": 0.9},  # 给出cdhit的参数coverage
            {"name": "memory_limit", "type": "int", "default": 10000},  # 内存大小，0为无限制
            {"name": "method", "type": "int", "default": 0},  # 1为全局比对，0为局部比对
            {"name": "direction", "type": "int", "default": 1},  # 1为双向比对，0为单向比对
            {"name": "num_threads", "type": "int", "default": 4},  # cpu数
            {"name": "select", "type": "int", "default": 1},  # 1为聚类到最相似的类中，0为聚类到第一个符合阈值的类
            {"name": "compare", "type": "string", "default": ""},  # 比对结果输出路径
            {"name": "pre", "type": "string", "default": "gene.geneset.tmp.fa.div-"},  # 文件前缀
            {"name": "output", "type": "outfile", "format": "sequence.fasta"}, # 输出fasta文件
            {"name": "ana_type", "type": "string", "default": "nucl"}  # 输入分析类型，是对核酸聚类还是随蛋白聚类
        ]
        self.add_option(options)
        self.step.add_steps('cdhitcomparesingle')
        self._memory_increase_step = 50  # 每次重运行增加内存50G by qingchen.zhang
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self.total = 0

    def step_start(self):
        self.step.cdhitcomparesingle.start()
        self.step.update()

    def step_end(self):
        self.step.cdhitcomparesingle.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        """
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query")

    def set_resource(self):
        """
        设置所需资源
        """
        infile = self.option("query").prop['path']
        number = os.path.getsize(infile) / 2000000 ###1G数据量用50G
        if number < 30:
            self.total = 30
        else:
            self.total = number
        self._cpu = self.option("num_threads")
        self._memory = str(self.total) + 'G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(PanCdhitAgent, self).end()


class PanCdhitTool(Tool):
    def __init__(self, config):
        super(PanCdhitTool, self).__init__(config)
        self._version = '1.0'
        self.cdhit_est_path = 'bioinfo/uniGene/cd-hit-v4.6.1-2012-08-27/cd-hit-est'
        self.cdhit_prot_path = 'bioinfo/uniGene/cd-hit-v4.6.1-2012-08-27/cd-hit'
        self.stat = os.path.join(self.work_dir, 'Orthologs_Cluster.txt')
        if self.option("compare") != "":
            self.compare = self.option("compare")
        else:
            self.compare = os.path.join(self.work_dir)

    def run(self):
        """
        运行
        :return:
        """
        super(PanCdhitTool, self).run()
        self.single_compare()
        cluster = self.run_stat()
        self.set_output()

    def word_len(self):
        """
        确定word_length
        :return:
        """
        word_length = 5
        if self.option("identity") >= 0.7:
            word_length = 5
        elif 0.6 <= self.option("identity") < 0.7:
            word_length = 4
        elif 0.5 <= self.option("identity") < 0.6:
            word_length = 3
        elif 0.4 <= self.option("identity") < 0.5:
            word_length = 2
        return word_length

    def single_compare(self):
        """
        用cdhit对序列进行自比挑选出代表序列
        :return:
        """
        infile = self.option("query").prop['path']
        number = os.path.getsize(infile) / 2000000  ###1G数据量用50G
        if number < 30:
            total = 30
        else:
            total = number
        length = self.word_len()
        out_dir = self.compare
        self.total_memory = total * 1000
        if os.path.exists(out_dir):
            pass
        else:
            os.mkdir(out_dir)
        if self.option("ana_type") == "nucl":
            cmd = '%s -i %s -o %s -c %s -aS %s -n %s -G %s -M %s -d %s -r %s -g %s -T %s' % (
                self.cdhit_est_path, self.option("query").prop['path'], out_dir + "/all_cluster.fa", self.option("identity"),
                self.option("coverage"), length, 1, self.total_memory, 0,
                self.option("direction"), self.option("select"), self.option("num_threads"))
        else:
            cmd = '%s -i %s -o %s -c %s -aS %s -n %s -G %s -M %s -d %s -g %s -T %s' % (
                self.cdhit_prot_path, self.option("query").prop['path'], out_dir + "/all_cluster.faa", self.option("identity"),
                self.option("coverage"), length, self.option("method"), self.total_memory, 0, self.option("select"), self.option("num_threads"))
        self.logger.info(cmd)
        command1 = self.add_command('cmd_1', cmd)
        command1.run()
        self.wait(command1)
        if command1.return_code == 0:
            # try:
                # if self.option("ana_type") == "nucl":
                #     self.option('output', out_dir + '/cluster.fa')
                # else:
                #     self.option('output', out_dir + '/cluster.faa')
                self.logger.info("compare single succeed at fist time")
            # except:
            # else:
            #     command1.rerun()
            #     self.wait(command1)
            #     if command1.return_code == 0:
            #         if self.option("ana_type") == "nucl":
            #             self.option('output', out_dir + "/cluster.fa")
            #         else:
            #             self.option('output', out_dir + "/cluster.faa")
            #         self.logger.info("compare single succeed at second time")
            #     else:
            #     self.set_error("cdhit cluster failed")
        else:
            self.set_error("cdhit cluster failed")

    def run_stat(self):
        """
        根据最后生成的序列id的过程文件生成结果统计表
        :return:
        """
        self.logger.info("开始对聚类结果进行统计和整理")
        if self.option("ana_type") == "nucl":
            input_cluster = os.path.join(self.compare, 'all_cluster.fa.clstr')
            input_fasta = os.path.join(self.compare, 'all_cluster.fa')
        else:
            input_cluster = os.path.join(self.compare, 'all_cluster.faa.clstr')
            input_fasta = os.path.join(self.compare, 'all_cluster.faa')
        sample_list = []
        # for seq_record in SeqIO.parse(input_fasta, 'fasta'):#打开第一次主要是获取到样本名称
        #     seq_id = seq_record.id
        #     sample = seq_id.split("|")[0]
        #     if sample not in sample_list:
        #         sample_list.append(sample)
        # sample_list.sort()
        all_sample = os.path.join(self.work_dir, 'homologues_cluster.xls')
        #new_all_sample = self.convert_cluster(, all_sample, sample_list)

        with open(input_cluster, 'r') as f:
            lines = f.readlines()
            cluster = {}
            cluster_sample_list = []
            cluster_index = 0
            for line in lines:
                if line.startswith(">"):
                    cluster_index += 1
                    #cluster_l = line.strip().split(" ")[1]
                    cluster_name = "CLUSTER" + str(cluster_index)
                    cluster_sample_list.append(cluster_name)
                    sample_cluster_list = []
                else:
                    cluster_name = cluster_sample_list[-1]
                    line = line.strip().split('\t')
                    m = re.match(r'(.+)\, \>(.+)\.\.\.(.+)$', line[1])
                    if m:
                        cluster_id = m.group(2)
                        sample_cluster_list.append(cluster_id)
                        sample_name = cluster_id.split("|")[0]
                        if sample_name not in sample_list:
                            sample_list.append(sample_name)
                        cluster[cluster_name] = sample_cluster_list
        sample_list.sort()
        #self.logger.info(cluster)
        with open(all_sample, 'w') as w:
            all_sample = "\t".join(sample_list)
            w.write("Cluster_ID\tSample_number\tGene_number\t" + all_sample + "\n")
            total_cluster = {}
            cluster_num = 0
            for key in cluster_sample_list:
                all_sample_cluster_list = cluster[key]
                all_sample_num =0
                all_gene_num =0
                cluster_num += 1
                cluster_sample_list = []
                cluster_id = "CLUSTER"+ str(cluster_num)
                total_cluster[key] = cluster_id
                for sample in sample_list:
                    sample_c_list = []
                    gene_num = 0
                    sample_num =0
                    for sample_cluster_id in all_sample_cluster_list:
                        if sample in sample_cluster_id:
                            sample_c_list.append(sample_cluster_id)
                    if len(sample_c_list) != 0:
                        sample_name = ','.join(sample_c_list)
                        gene_num = len(sample_c_list)
                        sample_num += 1
                    else:
                        sample_name = '-'
                    all_gene_num += gene_num
                    all_sample_num += sample_num
                    cluster_sample_list.append(sample_name)
                all_sample_id = '\t'.join(cluster_sample_list)#将不同样本名称的基因id放到一起并以"\t"分割
                w.write("{}\t{}\t{}\t{}\n".format(cluster_id, all_sample_num, all_gene_num, all_sample_id))
        self.logger.info("完成统计和整理聚类结果")
        return total_cluster

    def get_fasta(self, cluster):
        """
        根据OG0000892与cluster对应关系，获得聚类得到的fasta结果,改名换成cluster序列的编号
        :return:
        """
        self.logger.info("正在对代表序列进行改名")
        cluster_fasta = os.path.join(self.compare, 'cluster.faa')
        with open(cluster_fasta, 'w') as w:
            file_path = self.work_dir + '/all_rep_seq.fasta'
            for seq_record in SeqIO.parse(file_path, 'fasta'):#打开第一次主要是获取到聚类的结果id和样本名称
                seq = seq_record.seq
                seq_id = seq_record.id
                new_seq_id = cluster[seq_id]
                w.write(">{}\n{}\n".format(new_seq_id, seq))

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        if self.option("ana_type") == "nucl":
            newfile = os.path.join(self.output_dir, 'all_cluster.fa')
        else:
            newfile = os.path.join(self.output_dir, 'all_cluster.faa')
        if os.path.exists(newfile):
            os.remove(newfile)
        oldfiles = os.path.join(self.compare, 'cluster.faa')
        #os.link(oldfiles,newfile)
        cluster_stat = os.path.join(self.output_dir, 'homologues_cluster.xls')
        if os.path.exists(cluster_stat):
            os.remove(cluster_stat)
        os.link(self.work_dir + '/homologues_cluster.xls', cluster_stat)
        self.end()