# -*- coding: utf-8  -*-
# __author__ = 'gaohao'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict

class GenomeAssessModule(Module):
    """
    微生物基因组质控
    author: gaohao
    last_modify: 2018.03.21
    """
    def __init__(self, work_id):
        super(GenomeAssessModule, self).__init__(work_id)
        options = [
            {"name": "seq", "type": "infile", "format": "sequence.fasta"},  # 参考序列
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # 原始序列文件夹
            {"name": "kmer_type", "type": "int", "default": 17},  # kmer的大小
            {"name": "bases", "type": "int"},##计算时使用base数量
            {"name": "sample_name", "type": "string"},
        ]
        self.step.add_steps('gc_depth', 'cat_reads','genome_size')
        #self.gc_depth = self.add_tool('bacgenome.gc_depth')
        #self.gc_depth = self.add_tool('bacgenome.gc_depth')
        self.gc_depth = self.add_module('bacgenome.gc_depth')
        self.genome_size = self.add_tool('bacgenome.genome_size')
        self.cat_reads = self.add_tool('bacgenome.cat_reads')
        self.path = self.work_dir + '/list.txt'
        self.add_option(options)
        self.start_times = 0
        self.end_times = 0
        self.list = [self.gc_depth,self.genome_size]

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('bases'):
            raise OptionError('请提供总的bases数量！', code="21401501")
        if not self.option('seq').is_set:
            raise OptionError('请必须添加序列文件或序列文件夹！', code="21401502")
        if not self.option('fastq_dir').is_set:
            raise OptionError('必须输入原始序列文件夹！', code="21401503")
        return True

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def run_gc(self):
        self.get_info()
        opts = {
            "seq": self.option('seq'),
            "fastq_list": self.path,
        }
        self.gc_depth.set_options(opts)
        self.gc_depth.on('end', self.set_output, 'gc_depth')
        self.gc_depth.run()
        self.step.gc_depth.finish()
        self.step.update()
        self.start_times += 1

    def run_cat_reads(self):
        opts = {
            "map_dir": self.option('fastq_dir'),
        }
        self.cat_reads.set_options(opts)
        self.cat_reads.run()
        self.step.cat_reads.finish()
        self.step.genome_size.start()
        self.step.update()

    def run_size(self):
        self.genome_size.set_options({
            "fasta1": self.cat_reads.option('fasta1'),
            "fasta2": self.cat_reads.option('fasta2'),
            "bases": self.option('bases'),
            "sample_name":self.option('sample_name'),
        })
        self.genome_size.on('end', self.set_output, 'genome_size')
        self.genome_size.run()
        self.step.genome_size.finish()
        self.step.update()
        self.start_times += 1

    def get_info(self):
        list_path = os.path.join(self.option('fastq_dir').prop['path'],'list.txt')
        with open (list_path,'r') as f,open(self.path,'w') as file:
            sample_path = defaultdict(list)
            lines=f.readlines()
            for line in lines:
                tmp =line.rstrip('\r\n').split('\t')
                if re.search(r'PE',tmp[1]):
                    if tmp[1] in sample_path.keys():
                        if tmp[2] == 'l':
                            sample_path[tmp[1]].insert(0, self.option('fastq_dir').prop['path'] + '/' + tmp[0])
                        else:
                            sample_path[tmp[1]].append(self.option('fastq_dir').prop['path'] + '/' + tmp[0])
                    else:
                        sample_path[tmp[1]].append(self.option('fastq_dir').prop['path'] + '/' + tmp[0])

            for sample in sample_path:
                 file.write(self.option('sample_name') + '\t' + sample + '\t' + sample_path[sample][0] + ';' + sample_path[sample][1] + "\n")
            file.close()

    def run(self):
        """
        运行
        :return:
        """
        super(GenomeAssessModule, self).run()
        self.on_rely(self.list,self.end)
        self.run_gc()
        self.cat_reads.on('end', self.run_size)
        self.run_cat_reads()

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
                if os.path.exists(g1):
                    self.linkdir(g1, self.output_dir +  "/" + "depth_gc_" + str(j) + "/")
        if event["data"] == "genome_size":
            self.linkdir(self.genome_size.output_dir  + "/kmer_frequency/", self.output_dir  + "/kmer_frequency/")
            self.linkdir(self.genome_size.output_dir + "/genome_size/", self.output_dir + "/genome_size/")
            self.logger.info("设置结果目录成功")

    def end(self):
        super(GenomeAssessModule, self).end()