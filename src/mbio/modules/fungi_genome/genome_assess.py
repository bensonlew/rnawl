# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict
import shutil

class GenomeAssessModule(Module):
    """
    真菌基因组评估
    author: gaohao
    last_modify: 2018.06.22
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
        self.step.add_steps('gc_depth', 'cat_reads','genome_size','heter_stat','heter_ratio')
        #self.gc_depth = self.add_tool('fungi_genome.gc_depth')  #20200303
        self.gc_depth = self.add_module('bacgenome.gc_depth')
        self.genome_size = self.add_tool('bacgenome.genome_size')
        self.cat_reads = self.add_tool('bacgenome.cat_reads')
        self.heter_ratio = self.add_module('fungi_genome.heter_ratio')
        self.heter_stat = self.add_tool('fungi_genome.heter_stat')
        self.path = self.work_dir + '/list.txt'
        self.add_option(options)
        self.start_times = 0
        self.end_times = 0
        self.list = [self.genome_size,self.heter_ratio]
        self.list2 = [self.gc_depth,self.heter_stat]


    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('bases'):
            raise OptionError('请提供总的bases数量！', code="22101101")
        if not self.option('seq').is_set:
            raise OptionError('请必须添加序列文件或序列文件夹！', code="22101102")
        if not self.option('fastq_dir').is_set:
            raise OptionError('必须输入原始序列文件夹！', code="22101103")
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
            "windl": "10000, 20000, 30000, 50000, 100000"
        }
        self.gc_depth.set_options(opts)
        self.gc_depth.on('end', self.set_output, 'gc_depth')
        self.gc_depth.run()
        self.step.gc_depth.finish()
        self.step.update()

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

    def run_heter(self):
        self.heter_ratio.set_options({
            "scaf_fa": self.option('seq'),
        })
        self.heter_ratio.on('end', self.set_output, 'heter_ratio')
        self.heter_ratio.run()
        self.step.heter_ratio.finish()
        self.step.update()

    def run_heter_stat(self):
        self.heter_stat.set_options({
            "frequency": self.genome_size.work_dir + "/" + self.option('sample_name') + ".frequency.xls",
            "heter1": self.heter_ratio.option('heter1'),
            "heter2":self.heter_ratio.option('heter2'),
            "heter3": self.heter_ratio.option('heter3'),
            "heter4": self.heter_ratio.option('heter4'),
            "heter5": self.heter_ratio.option('heter5'),
            "heter6": self.heter_ratio.option('heter6'),
        })
        self.heter_stat.on('end', self.set_output, 'heter_stat')
        self.heter_stat.run()
        self.step.heter_stat.finish()
        self.step.update()

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
        self.on_rely(self.list2,self.end)
        self.on_rely(self.list,self.run_heter_stat)
        self.run_heter()
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
        if event["data"] == "gc_depth":
            for j in 10000, 20000, 30000, 50000, 100000:
                g1 = self.gc_depth.output_dir + "/" + "depth_gc_" + str(j) + "/"
                if not os.path.exists(g1):   #zouguanqing 20181018
                    continue
                if os.path.exists(self.output_dir +  "/" + "depth_gc_" + str(j) + "/"):
                    shutil.rmtree(self.output_dir +  "/" + "depth_gc_" + str(j) + "/")
                self.linkdir(g1, self.output_dir +  "/" + "depth_gc_" + str(j) + "/")
            if os.path.exists(self.output_dir + '/' + self.option('sample_name') + "_gc_depth.xls"):
                os.remove(self.output_dir + '/' + self.option('sample_name') + "_gc_depth.xls")    
            #os.link(self.gc_depth.work_dir + '/cov.depth.txt',self.output_dir + '/' + self.option('sample_name') + "_gc_depth.xls")
            os.link(self.gc_depth.work_dir + '/GcDepthStep/depth.final.out',self.output_dir + '/' + self.option('sample_name') + "_gc_depth.xls") #20200303

        if not os.path.exists(self.output_dir  + "/kmer_frequency/"):
            os.mkdir(self.output_dir  + "/kmer_frequency/")
        if event["data"] == "genome_size":
            if os.path.exists(self.output_dir + "/genome_size/"):
                shutil.rmtree(self.output_dir + "/genome_size/")
            self.linkdir(self.genome_size.output_dir + "/genome_size/", self.output_dir + "/genome_size/")
            if os.path.exists(self.output_dir + "/kmer_frequency/" + self.option('sample_name') + ".frequency.xls"):
                os.remove(self.output_dir + "/kmer_frequency/" + self.option('sample_name') + ".frequency.xls")
            os.link(self.genome_size.output_dir  + "/kmer_frequency/" + self.option('sample_name') + ".frequency.xls", self.output_dir  + "/kmer_frequency/" + self.option('sample_name') + ".frequency.xls")
        if event["data"] == "heter_stat":
            if os.path.exists(self.output_dir  + "/kmer_frequency/" + self.option('sample_name') + "_Kmer_frequency.xls"):
                os.remove(self.output_dir  + "/kmer_frequency/" + self.option('sample_name') + "_Kmer_frequency.xls")
            os.link(self.heter_stat.output_dir  + "/all_heter_kmer_freq.xls", self.output_dir  + "/kmer_frequency/" + self.option('sample_name') + "_Kmer_frequency.xls")
        self.logger.info("设置结果目录成功")

    def end(self):
        super(GenomeAssessModule, self).end()