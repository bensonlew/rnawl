# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
# version 1.0
# last_modify: 2018.8.14
import os
import re,shutil
from Bio import SeqIO
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class GcDepthBinAgent(Agent):
    """
    细菌基因组depth_gc评估
    """
    def __init__(self, parent):
        super(GcDepthBinAgent, self).__init__(parent)
        options = [
            {"name": "sam", "type": "infile", "format": "align.bwa.sam_dir"},   #输入文件 Sam文件夹
            {"name": "windl", "type": "int", "default": "1000"},  # 滑动窗口大小
            {"name": "ref","type":"infile","format":"sequence.fasta"},#组装结果的基因组非binning的结果
            {"name": "fastq1", "type": "infile", "format": "sequence.fastq"},  # 输入文件 mapping后的fastq1序列
            {"name": "fastq2", "type": "infile", "format": "sequence.fastq"},  # 输入文件 mapping后的fastq2序列
            {"name": "fastqs", "type": "infile", "format": "sequence.fastq"},  # 输入文件 mapping后的fastqs序列
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sam"):
            raise OptionError("必须添加sam文件！")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(GcDepthBinAgent, self).end()

class GcDepthBinTool(Tool):
    def __init__(self, config):
        super(GcDepthBinTool, self).__init__(config)

        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.perl_script = self.config.PACKAGE_DIR + "/metagbin/"
        self.R_path = 'program/R-3.3.1_gcc5.1/bin/Rscript'

        self.depth = self.work_dir + "/" + 'cov.depth.txt'
        self.samtools = "bioinfo/align/samtools-1.4/bin/samtools"
        self.bac_depth = self.config.PACKAGE_DIR + '/bacgenome/bac_depth.sh'

    def run_pick(self):
        self.logger.info('开始生成bam文件')
        n = 0
        self.logger.info('%s'%(self.option('sam').prop['path']))
        for file in os.listdir(self.option('sam').prop['path']):
            n = n+1
            sam_path = os.path.join(self.option('sam').prop['path'], file)
            self.out_bam = self.work_dir + '/' + str(n) + '_map.bam'
            self.logger.info(self.out_bam)
            cmd_bam = "{} view -F 4 -S {} -b -o {}".format(self.samtools, sam_path, self.out_bam)
            self.logger.info(cmd_bam)
            to_bam = str(n) + 'to-bam'
            command1 = self.add_command(to_bam, cmd_bam).run()
            self.wait(command1)
            if command1.return_code == 0:
                self.logger.info("%s运行完成"%cmd_bam)
            else:
                self.set_error("%s运行失败", variables=(cmd_bam))
            self.logger.info('对bam文件进行排序')
            self.sort_bam = self.work_dir + '/' + str(n) +'_sort.bam'
            cmd = '{} sort {} -o {}'.format(self.samtools, self.out_bam, self.sort_bam)
            to_sort = str(n) + 'sort-bam'
            command = self.add_command(to_sort,cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("%s运行完成"%cmd)
            else:
                self.set_error("%s运行失败", variables=(cmd))

    def run_merge(self):
        self.logger.info('开始将所有bam文件进行合并')
        if os.path.exists(self.work_dir + '/merge.bam'):
            os.remove(self.work_dir + '/merge.bam')
        self.merge_bam = self.work_dir + '/merge.bam'
        self.logger.info(self.merge_bam)
        cmd_bam = "{} merge {} {} {}".format(self.samtools, self.merge_bam,self.work_dir + '/1_sort.bam', self.work_dir + '/2_sort.bam' )
        self.logger.info(cmd_bam)
        to_merge = 'to_merge'
        command1 = self.add_command(to_merge, cmd_bam).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("%s运行完成"%to_merge)
        elif command1.return_code == 1:
            self.add_state('merge.bam文件已经存在', 'memory is low!')
        else:
            self.set_error("运行失败")

    def run_depth(self):
        self.logger.info('开始计算覆盖度')
        self.tmp_sh = self.work_dir + '/tmp.sh'
        self.depth = self.work_dir + '/depth.txt'
        if len(os.listdir(self.option('sam').prop['path'])) >= 2:
            self.merge_bam = self.work_dir + '/merge.bam'
        else:
            self.merge_bam = self.work_dir + '/1_sort.bam'
        self.tmp_sh = "{} {} {} {}".format(self.config.PACKAGE_DIR + '/bacgenome/samtools_depth.sh', self.config.SOFTWARE_DIR ,self.merge_bam, self.depth)
        cmd1 = '/program/sh {}'.format(self.tmp_sh)
        self.logger.info('1111111111')
        command1 = self.add_command('depth', cmd1).run()
        self.logger.info('bbbbbbbbb')
        self.wait(command1)
        self.logger.info('ccccccccc')
        if command1.return_code == 0:
            self.logger.info("%s运行成功"%cmd1)
        else:
            self.set_error("%s运行失败", variables=(cmd1))

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
        self.ref=self.option('ref').prop['path']
        self.logger.info("run_depth_gc运行开始")
        wind = self.option('windl')
        out_dir = self.work_dir + "/" + 'depth_gc_' + str(wind) + "/"
        cmd = '{} {}gc_depth.pl {} {} -windl {} -step 500 -outdir {}'.format(self.perl_path,self.perl_script,self.ref,self.depth_final,wind,out_dir)
        gc_depth =  'depth_gc_' + str(wind)
        command = self.add_command(gc_depth, cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info(gc_depth + "运行完成")
        else:
            self.set_error("%s运行出错!", variables=(gc_depth))
        cmd_1 = '%s %s' % (self.R_path, out_dir + "depth_gc.cmd.r")
        self.logger.info(cmd_1)
        gc_depth_r = 'depth_gc_' + str(wind) +'_r'
        command1 = self.add_command(gc_depth_r, cmd_1).run()
        #command1.run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("R程序计算depth_gc成功")
        else:
            self.set_error("R程序计算depth_gc失败")
            self.set_error("R程序计算depth_gc失败")
        self.logger.info("run_depth_gc运行结束")

    def run_cal_coverage(self):
        """
        计算组装结果的reads的测序深度
        :return:
        """
        value_1 = 0
        value_2 = 0
        value_s = 0
        value_r = 0
        self.logger.info('正在开始计算coverage')
        coverage_path = self.output_dir + '/Genome_coverage.xls'
        with open(coverage_path, 'w') as w:
            #w.write("#Coverage\n")
            if self.option('fastqs').is_set:
                for seq_record in SeqIO.parse(self.option('fastq1').prop['path'], 'fastq'):
                    length = int(len(seq_record.seq))
                    value_1 += length
                self.logger.info("value_1: %s" %value_1)
                for seq_record in SeqIO.parse(self.option('fastq2').prop['path'], 'fastq'):
                    length = int(len(seq_record.seq))
                    value_2 += length
                self.logger.info("value_2: %s" %value_2)
                for seq_record in SeqIO.parse(self.option('fastqs').prop['path'], 'fastq'):
                    length = int(len(seq_record.seq))
                    value_s += length
                self.logger.info("value_s: %s" %value_s)
                for seq_record in SeqIO.parse(self.option('ref').prop['path'], 'fasta'):
                    length = int(len(seq_record.seq))
                    value_r += length
                self.logger.info("value_r: %s" %value_r)
                coverage = int((value_1 + value_2 + value_s)/value_r)
                self.logger.info("assembly_coverage: %s" %coverage)
                w.write('#Assembly_coverage\n{}\n'.format(coverage))
            else:
                for seq_record in SeqIO.parse(self.option('fastq1').prop['path'], 'fastq'):
                    length = int(len(seq_record.seq))
                    value_1 += length
                self.logger.info("value_1: %s" %value_1)
                for seq_record in SeqIO.parse(self.option('fastq2').prop['path'], 'fastq'):
                    length = int(len(seq_record.seq))
                    value_2 += length
                self.logger.info("value_2: %s" %value_2)
                for seq_record in SeqIO.parse(self.option('ref').prop['path'], 'fasta'):
                    length = int(len(seq_record.seq))
                    value_r += length
                self.logger.info("value_r: %s" %value_r)
                coverage = int((value_1 + value_2)/value_r)
                self.logger.info("assembly_coverage: %s" %coverage)
                w.write('#Assembly_coverage\n{}\n'.format(coverage))

    def set_output(self):
        l1 = self.work_dir + "/" + "depth_gc_" + str(1000) + "/"
        for type in ['png', 'svg']:
            pre_name=os.path.basename(self.option('ref').prop['path']).split('_')
            sample_name = "G_" + pre_name[1]
            if os.path.exists(self.output_dir + '/gc_depth.' + type):
                os.remove(self.output_dir + '/gc_depth.' + type)
            os.link(l1 + 'gc_depth.' + type, self.output_dir + '/gc_depth.' + type)

    def run(self):
        super(GcDepthBinTool, self).run()
        self.run_pick()
        if len(os.listdir(self.option('sam').prop['path'])) >= 2:
            self.run_merge()
        self.run_depth()
        self.change_depth_format()
        self.run_depth_gc()
        self.run_cal_coverage()
        self.set_output()
        self.end()


