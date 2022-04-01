# -*- coding: utf-8 -*-
# __author__ = 'guanqing.zou hao.gao'
# version 1.0
# last_modify: 2018.8.14
import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
class GcDepthBwaAgent(Agent):
    """
    细菌基因组depth_gc评估
    """
    def __init__(self, parent):
        super(GcDepthBwaAgent, self).__init__(parent)
        options = [
            {"name": "seq_dir", "type": "infile","format": "sequence.fasta_dir"},  # 参考序列文件夹
            {"name": "fastq_list", "type": "infile","format": "meta.profile"},  # 输入列表文件
            {"name": "seq", "type": "infile", "format": "sequence.fasta"},  # 参考序列
            {"name": "windl", "type": "int", "default": "1000,3000,5000,8000,10000"},  # 滑动窗口大小
            {"name": "fq_file_id", "type": "int", "default": 0}
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fastq_list").is_set:
            raise OptionError("必须添加fastq_list的列表文件！")
        if not self.option('seq_dir').is_set and not self.option('seq').is_set:
            raise OptionError('请必须添加序列文件或序列文件夹！')

    def set_resource(self):
        self._cpu = 5
        self._memory = '20G'

    def end(self):
        super(GcDepthBwaAgent, self).end()

class GcDepthBwaTool(Tool):
    def __init__(self, config):
        super(GcDepthBwaTool, self).__init__(config)
        self.fastq_list = self.option("fastq_list").prop['path']
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/"
        #self.bwt_builder= "bioinfo/uniGene/soap2.21release/2bwt-builder"
        #self.soap = "bioinfo/uniGene/soap2.21release/soap"
        #self.soap_coverage = "bioinfo/uniGene/soap2.21release/2.7.7/soap.coverage"
        self.R_path = 'program/R-3.3.1/bin/Rscript'

        self.ref = self.work_dir + "/" + 'ref.fna'
        if self.option('seq').is_set:
            os.system('cp {} {}'.format(self.option('seq').prop['path'],self.ref))

        self.sh ="../../../../../.." + self.config.PACKAGE_DIR + '/sequence/scripts/depth_gc_cat.sh'
        self.depth = self.work_dir + "/" + 'cov.depth.txt'
        #self.file = self.work_dir + "/" + 'soap.list'
        #self.bwa = 'bioinfo/align/bwa-0.7.9a/bwa'
        self.bowtie = 'bioinfo/align/bowtie2'
        self.samtools = "bioinfo/align/samtools-1.4/bin/samtools"
        self.bac_depth = self.config.PACKAGE_DIR + '/bacgenome/bac_depth.sh'

    def run_ref_cat(self):
        self.dir = self.option("seq_dir").prop['path']
        cmd =self.sh + ' ' + self.dir + '/*' + ' ' + self.ref
        self.logger.info(cmd)
        command = self.add_command("run_ref_cat", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("生成合并sh运行完成")
        else:
            self.set_error("生成合并sh运行出错!")

    def run_bwa_depth(self):
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/align/bwa-0.7.9a/bwa')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/align/samtools-1.4/bin')
        self.depth = self.work_dir + '/depth.txt'
        self.logger.info('bwa 比对开始')
        with open(self.fastq_list,'r') as f:
            lines = f.readlines()
            line = lines[self.option('fq_file_id')]
            line = line.strip().split('\t')
            (self.read1,self.read2)=line[2].split(';')
        cmd1 = '/program/sh {} {} {} {} {}'.format(self.bac_depth, self.ref, self.read1, self.read2, self.depth)
        command1 = self.add_command('bwa-depth', cmd1)
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("%运行成功"%cmd1)
        else:
            self.set_error("%s运行失败"%cmd1)

    '''
    def run_ref_index(self):
        cmd = '{} index {}'.format(self.bwa, self.ref)
        #cmd = '{} {}'.format(self.bwt_builder,ref)
        self.logger.info(cmd)
        command = self.add_command("bwa-index", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("bwa index运行完成")
        else:
            self.set_error("bwa index运行出错!")

    def run_aln(self):
        self.logger.info('bwa 比对开始')
        with open(self.fastq_list,'r') as f:
            lines = f.readlines()
            line = lines[self.option('fq_file_id')]
            line = line.strip().split('\t')
            PEname=line[1]
            (self.read1,self.read2)=line[2].split(';')
            self.sai1 = self.work_dir + '/out.1.sai'
            self.sai2 = self.work_dir + '/out.2.sai'
            cmd = '{} aln {} {} -f {}'.format(self.bwa, self.ref, self.read1, self.sai1)
            cmd1 = '{} aln {} {} -f {}'.format(self.bwa, self.ref, self.read2, self.sai2)
            command = self.add_command('aln1', cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("%s运行完成"%cmd)
            else:
                self.set_error("%s 运行出错!" %cmd)

            command1 = self.add_command('aln2', cmd1).run()
            self.wait(command1)
            if command.return_code == 0:
                self.logger.info("%s运行完成"%cmd1)
            else:
                self.set_error("%s 运行出错!" %cmd1)

    def run_sampe(self):
        self.out_sam = self.work_dir + '/map.sam'
        cmd_sam = '{} sampe -a 500 -f {} {} {} {} {} {}'.format(self.bwa, self.out_sam, self.ref, self.sai1, self.sai2, self.read1, self.read2)
        self.logger.info('开始运行bwa sampe')
        self.logger.info('{}'.format(cmd_sam))
        command = self.add_command('bwa-sampe', cmd_sam)
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("%s运行完成"%cmd_sam)
        else:
            self.set_error("%s运行失败"%cmd_sam)

        self.logger.info('开始生成bam文件')
        self.out_bam = self.work_dir + '/map.bam'
        cmd_bam = "{} view -F 12 -S {} -b -o {}".format(self.samtools, self.out_sam, self.out_bam)
        command1 = self.add_command('to-bam', cmd_bam)
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("%s运行完成"%cmd_bam)
        else:
            self.set_error("%s运行失败"%cmd_bam)

    def run_depth(self):
        self.logger.info('对bam文件进行排序')
        self.sort_bam = self.work_dir + '/sort'
        cmd = '{} sort {} -o {}'.format(self.samtools, self.out_bam, self.sort_bam)
        command = self.add_command('sort-bam',cmd)
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("%s运行完成"%cmd)
        else:
            self.set_error("%s运行失败"%cmd)

        self.logger.info('开始计算覆盖度')
        self.tmp_sh = self.work_dir + '/tmp.sh'
        self.depth = self.work_dir + '/depth.txt'
        os.system('echo {} depth {} \> {} > {}'.format(self.samtools, self.sort_bam, self.depth, self.tmp_sh))
        cmd1 = '/program/sh {}'.format(self.tmp_sh)
        command1 = self.add_command('depth', cmd1)
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("%运行成功"%cmd1)
        else:
            self.set_error("%s运行失败"%cmd1)
        '''

    def change_depth_format(self):
        self.logger.info('开始转化depth文件格式')
        self.depth_final = self.work_dir + '/depth.final.out'
        fw = open(self.depth_final, 'w')
        with open(self.depth) as f:
            line1 = f.next()
            line1 = line1.strip()
            spline1= line1.split('\t')
            k = spline1[0]
            fw.write('>'+k+'\n'+spline1[2])
            for i in f:
                i = i.strip()
                spi = i.split('\t')
                if spi[0] == k:
                    fw.write(' '+spi[2])
                else:
                    k = spi[0]
                    fw.write('\n>'+k+"\n"+spi[2])
        fw.close()
        self.logger.info('转化depth文件格式完成，生成%s'%self.depth_final)


    def run_depth_gc(self):
        wind_arrry = (self.option("windl")).split(',')
        self.logger.info("run_depth_gc运行开始")
        #ref = ''
        #if self.option('seq').is_set:
        #    ref = self.option('seq').prop['path']
        #elif self.option('seq_dir').is_set:
        #    ref = self.ref
        for wind in wind_arrry:
            out_dir = self.work_dir + "/" + 'depth_gc_' + str(wind) + "/"
            cmd = '{} {}gc_depth_v2.0.pl {} {} -windl {} -step 500 -outdir {}'.format(self.perl_path,self.perl_script,self.ref,self.depth_final,wind,out_dir)
            gc_depth =  'depth_gc_' + str(wind)
            command = self.add_command(gc_depth, cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info(gc_depth + "运行完成")
            else:
                self.set_error(gc_depth + "运行出错!")
            cmd_1 = '%s %s' % (self.R_path, out_dir + "depth_gc.cmd.r")
            self.logger.info(cmd_1)
            gc_depth_r = 'depth_gc_' + str(wind) +'_r'
            command1 = self.add_command(gc_depth_r, cmd_1)
            command1.run()
            self.wait(command1)
            if command1.return_code == 0:
                self.logger.info("R程序计算depth_gc成功")
            else:
                self.set_error("R程序计算depth_gc失败")
                raise Exception("R程序计算depth_gc失败")
        self.logger.info("run_depth_gc运行结束")

    def set_output(self):
        for j in 1000, 3000, 5000, 8000, 10000:
            g1 = self.output_dir + "/" + "depth_gc_" + str(j) + "/"
            l1 = self.work_dir + "/" + "depth_gc_" + str(j) + "/"
            if os.path.exists(g1):
                pass
            else:
                os.mkdir(g1)
            for type in ['png', 'svg']:
                if os.path.exists(g1 + 'gc_depth.' + type):
                    os.remove(g1 + 'gc_depth.' + type)
                os.link(l1 + 'gc_depth.' + type, g1 + 'gc_depth.' + type)


    def run(self):
        super(GcDepthBwaTool, self).run()
        if self.option('seq_dir').is_set:
            self.run_ref_cat()
        self.run_bwa_depth()
        #self.run_ref_index()
        #self.run_aln()
        #self.run_sampe()
        #self.run_depth()
        self.change_depth_format()
        self.run_depth_gc()
        self.set_output()
        self.end()


