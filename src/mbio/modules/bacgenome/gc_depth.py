# -*- coding: utf-8 -*-
# __author__ = 'guanqing.zou'
# last_modify:201808

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class GcDepthModule(Module):
    def __init__(self, work_id):
        super(GcDepthModule, self).__init__(work_id)
        options = [
            {"name": "seq_dir", "type": "infile","format": "sequence.fasta_dir"},
            {"name": "fastq_list", "type": "infile","format": "meta.profile"},  # 输入列表文件
            #{"name": "fastq_list", "type": "string"},
            {"name": "seq", "type": "infile", "format": "sequence.fasta"},  # 参考序列
            #{"name": "seq", "type": "string","default":''},
            {"name": "fq_file_id", "type": "int", "default": 0},
            {"name": "windl", "type": "string", "default": "1000,3000,5000,8000,10000"},  # 滑动窗口大小

        ]
        self.add_option(options)
        self.align = self.add_tool("align.bowtie2")
        self.gcdepth = self.add_tool("bacgenome.gc_depth_step")
        #self.fastq_list = self.option("fastq_list").prop['path']
        #self.fastq_list = self.option("fastq_list")
        self.ref = self.work_dir + "/" + 'ref.fna'

        

    def check_options(self):
        if not self.option("fastq_list").is_set:
            raise OptionError("必须添加fastq_list的列表文件！", code="21401901")
        

    def run_align(self):
        self.fastq_list = self.option("fastq_list").prop['path']
        self.ref = self.work_dir + "/" + 'ref.fna'
        if self.option('seq').is_set:
            os.system('cp {} {}'.format(self.option('seq').prop['path'],self.ref))
        if self.option('seq_dir').is_set:
            self.dir = self.option('seq_dir').prop['path']
            os.system('cat {}/* > {}'.format(self.dir,self.ref))
        #with open(self.fastq_list,'r') as f:
        f = open(self.option("fastq_list").prop['path'])
        lines = f.readlines()
        line = lines[self.option('fq_file_id')]
        line = line.strip().split('\t')
        (self.read1,self.read2)=line[2].split(';')
        self.sam_prex = os.path.basename(self.read1).split('.')[0]
        self.align.set_options({
            "ref_fasta":self.ref,
            "fastq1": self.read1,
            "fastq2":self.read2
        })
        self.align.on('end',self.run_gc_depth)
        self.align.run()

    def run_gc_depth(self):
        self.logger.info(self.align.option('sam_file').prop['path']+'/{}.pair.sam'.format(self.sam_prex))
        self.gcdepth.set_options({
            'sam':self.align.option('sam_file').prop['path']+'/{}.pair.sam'.format(self.sam_prex),
            'ref':self.ref,
            'windl': self.option('windl')
        })
        self.gcdepth.on('end',self.set_output)
        self.logger.info(self.align.option('sam_file').prop['path'])
        self.gcdepth.run()


    def set_output(self):
        #for j in 1000, 3000, 5000, 8000, 10000:
        for j in map(int, self.option('windl').split(',')):
            g1 = self.gcdepth.output_dir + "/" + "depth_gc_" + str(j) + "/"
            if not os.path.exists(g1):
                break
            if os.path.exists(self.output_dir +  "/" + "depth_gc_" + str(j) + "/"):
                shutil.rmtree(self.output_dir +  "/" + "depth_gc_" + str(j) + "/")
            self.linkdir(g1, self.output_dir +  "/" + "depth_gc_" + str(j) + "/")
        self.end()

    def run(self):
        super(GcDepthModule, self).run()
        self.run_align()

    def end(self):
        super(GcDepthModule, self).end()



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