# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict
class GenomeAssess2Module(Module):
    """
    微生物基因组质控数据库比对没有加入
    author: guhaidong
    last_modify: 2019.04.16
    """
    def __init__(self, work_id):
        super(GenomeAssess2Module, self).__init__(work_id)
        options = [
            {"name": "seq", "type": "infile", "format": "sequence.fasta", "required": True},  # 参考序列
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # 原始序列文件夹
            {"name": "kmer_type", "type": "int", "default": 17},  # kmer的大小
            {"name": "bases", "type": "int"},  #计算时使用base数量
            {"name": "sample_name", "type": "string", "required": True},
        ]
        self.gc_depth = self.add_module('bacgenome.gc_depth')
        self.genome_size = self.add_tool('bacgenome.genome_size')
        self.cat_reads = self.add_tool('bacgenome.cat_reads')
        self.pca = self.add_tool("bacgenome.kmer_pca")
        self.sixteens = self.add_module("bacgenome.blast_gene")
        self.hgene = self.add_tool("bacgenome.hgene_blastx")
        self.path = self.work_dir + '/list.txt'
        self.add_option(options)
        self.start_times = 0
        self.end_times = 0
        self.list = [self.gc_depth,self.genome_size, self.pca, self.sixteens, self.hgene]
        self.list2 = [self.pca, self.sixteens, self.hgene]
        self.fq1 = ""
        self.fq2 = ""

    def check_options(self):
        """
        检查参数
        :return:
        """
        return True

    def run_gc(self):
        opts = {
            "seq": self.option('seq'),
            "fastq_list": self.path,
        }
        self.gc_depth.set_options(opts)
        self.gc_depth.on('end', self.set_output, 'gc_depth')
        self.gc_depth.run()
        self.start_times += 1

    def run_pca(self, coverage=None):
        opts = {
            "scaf": self.option("seq")
        }
        if coverage:
            opts["cov"] = coverage
        self.pca.set_options(opts)
        self.pca.on("end", self.set_output, "pca")
        self.pca.run()

    def run_sixteens(self):
        self.sixteens.set_options({
            "fa": self.option("seq")
        })
        self.sixteens.on("end", self.set_output, "sixteens")
        self.sixteens.run()

    def run_hgene(self):
        self.hgene.set_options({
            "seq_fa": self.option("seq")
        })
        self.hgene.on("end", self.set_output, "hgene")
        self.hgene.run()

    def run_cat_reads(self):
        # 不用这个了
        opts = {
            "map_dir": self.option('fastq_dir'),
        }
        self.cat_reads.set_options(opts)
        self.cat_reads.run()

    def run_size(self):
        self.logger.info("=========genome_assess2 debug===========")
        self.logger.info(self.fq1)
        self.logger.info(self.fq2)
        self.genome_size.set_options({
            "fasta1": self.fq1,  # self.cat_reads.option('fasta1'),
            "fasta2": self.fq2,  # self.cat_reads.option('fasta2'),
            "bases": self.option('bases'),
            "sample_name":self.option('sample_name'),
        })
        self.genome_size.on('end', self.set_output, 'genome_size')
        self.genome_size.run()
        self.start_times += 1

    def get_info(self):
        list_path = os.path.join(self.option('fastq_dir').prop['path'],'list.txt')
        with open (list_path,'r') as f,open(self.path,'w') as file:
            sample_path = defaultdict(list)
            lines=f.readlines()
            for line in lines:
                tmp =line.rstrip('\r\n').split('\t')
                if tmp[1] in sample_path.keys():
                    if tmp[2] == 'l':
                        sample_path[tmp[1]].insert(0, self.option('fastq_dir').prop['path'] + '/' + tmp[0])

                    else:
                        sample_path[tmp[1]].append(self.option('fastq_dir').prop['path'] + '/' + tmp[0])

                else:
                    sample_path[tmp[1]].append(self.option('fastq_dir').prop['path'] + '/' + tmp[0])
                if tmp[2] == 'l':
                    self.fq1 = self.option("fastq_dir").prop["path"] + "/" + tmp[0]
                elif tmp[2] == 'r':
                    self.fq2 = self.option("fastq_dir").prop["path"] + "/" + tmp[0]

            for sample in sample_path:
                 file.write(self.option('sample_name') + '\t' + sample + '\t' + sample_path[sample][0] + ';' + sample_path[sample][1] + "\n")
            file.close()

    def run(self):
        """
        运行
        :return:
        """
        # self.cat_reads.on('end', self.run_size)
        super(GenomeAssess2Module, self).run()
        coverage = os.path.join(os.path.dirname(self.option("seq").prop["path"]), self.option("sample_name") + ".abund")
        coverage = coverage if os.path.isfile(coverage) else None
        if self.option("fastq_dir").is_set:
            self.get_info()
            # self.on_rely(self.list,self.end)
            self.run_size()
            self.run_gc()
        # else:
        #     self.on_rely(self.list2, self.end)
        self.run_pca(coverage=coverage)
        self.run_sixteens()
        self.run_hgene()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        newdir = os.path.join(self.work_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        allfiles = os.listdir(dirpath)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.link(oldfiles[i], newdir)
    def set_output(self, event):
        """
        将结果文件复制到output文件夹下面
        :return
        """
        self.logger.info("设置结果目录")
        self.end_times += 1
        if event["data"] == "gc_depth":
            for j in 1000, 3000, 5000, 8000, 10000:
                g1 = self.gc_depth.output_dir + "/" + "depth_gc_" + str(j) + "/"
                if not os.path.exists(g1):
                    break
                self.linkdir(g1, self.output_dir +  "/" + "depth_gc_" + str(j) + "/")
        if event["data"] == "genome_size":
            self.linkdir(self.genome_size.output_dir  + "/kmer_frequency/", self.output_dir  + "/kmer_frequency/")
            self.linkdir(self.genome_size.output_dir + "/genome_size/", self.output_dir + "/genome_size/")
        if event["data"] == "pca":
            self.linkdir(self.pca.output_dir, self.output_dir + "/kmer_pca")
        if event["data"] == "sixteens":
            self.linkdir(self.sixteens.output_dir, self.output_dir + "/organism")
        elif event["data"] == "hgene":
            self.linkdir(self.hgene.output_dir, self.output_dir + "/organism")
        if self.option("fastq_dir").is_set:
            check_end_list = self.list
        else:
            check_end_list = self.list2
        for tool in check_end_list:
            if not tool.is_end:
                return
        self.end()

    def end(self):
        if self.is_end:
            return
        super(GenomeAssess2Module, self).end()