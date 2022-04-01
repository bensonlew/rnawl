# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2018.1.2
import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
class GcDepthAgent(Agent):
    """
    细菌基因组depth_gc评估
    """
    def __init__(self, parent):
        super(GcDepthAgent, self).__init__(parent)
        options = [
            {"name": "seq_dir", "type": "infile","format": "sequence.fasta_dir"},  # 参考序列文件夹
            {"name": "fastq_list", "type": "infile","format": "meta.profile"},  # 输入列表文件
            {"name": "seq", "type": "infile", "format": "sequence.fasta"},  # 参考序列
            {"name": "windl", "type": "int", "default": "1000,3000,5000,8000,10000"}  # 滑动窗口大小
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fastq_list").is_set:
            raise OptionError("必须添加fastq_list的列表文件！", code="31401801")
        if not self.option('seq_dir').is_set and not self.option('seq').is_set:
            raise OptionError('请必须添加序列文件或序列文件夹！', code="31401802")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(GcDepthAgent, self).end()

class GcDepthTool(Tool):
    def __init__(self, config):
        super(GcDepthTool, self).__init__(config)
        self.fastq_list = self.option("fastq_list").prop['path']
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/"
        self.bwt_builder= "bioinfo/uniGene/soap2.21release/2bwt-builder"
        self.soap = "bioinfo/uniGene/soap2.21release/soap"
        self.soap_coverage = "bioinfo/uniGene/soap2.21release/2.7.7/soap.coverage"
        self.R_path = 'program/R-3.3.1_gcc5.1/bin/Rscript'
        self.ref =self.work_dir + "/" + 'ref.fna'
        self.sh ="../../../../../.." + self.config.PACKAGE_DIR + '/sequence/scripts/depth_gc_cat.sh'
        self.depth = self.work_dir + "/" + 'cov.depth.txt'
        self.file = self.work_dir + "/" + 'soap.list'

    def run_ref_cat(self):
        self.dir = self.option("seq_dir").prop['path']
        cmd =self.sh + ' ' + self.dir + '/*' + ' ' + self.ref
        self.logger.info(cmd)
        command = self.add_command("run_ref_cat", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("生成合并sh运行完成")
        else:
            self.set_error("生成合并sh运行出错!", code="31401801")


    def run_ref_index(self):
        ref = ''
        if self.option('seq').is_set:
            ref =self.option('seq').prop['path']
        elif self.option('seq_dir').is_set:
            ref =self.ref
        cmd = '{} {}'.format(self.bwt_builder,ref)
        self.logger.info(cmd)
        command = self.add_command("bwt-builder", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("2bwt-builder运行完成")
        else:
            self.set_error("2bwt-builder运行出错!", code="31401802")

    def run_soap(self):
        self.logger.info("总run_soap运行开始")
        index =''
        if self.option('seq').is_set:
            index =self.option('seq').prop['path'] + '.index'
        elif self.option('seq_dir').is_set:
            index =self.ref + '.index'
        file =open(self.file,'w')
        with open(self.fastq_list,'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split('\t')
                PEname=line[1]
                (read1,read2)=line[2].split(';')
                single = self.work_dir + "/" + PEname + '.aln.single'
                file.write(single + '\n')
                soap = self.work_dir + "/" + PEname + '.aln.soap'
                file.write(soap + '\n')
                cmd = '{} -a {} -b {} -D {} -g 3 -x 580 -s 35 -l 32 -v 2 -m 380 -p 20 -2 {} -o {} '.format(self.soap,read1,read2,index,single,soap)
                path= PEname.lower() + 'run_soap'
                command = self.add_command(path, cmd).run()
                self.wait(command)
                if command.return_code == 0:
                    self.logger.info("%s运行完成" %path)
                else:
                    self.set_error("%s 运行出错!" , variables=(path), code="31401803")
        self.logger.info("总run_soap运行结束")

    def run_coverage(self):
        detail = self.work_dir + "/" + 'cov.detail.txt'
        ref = ''
        if self.option('seq').is_set:
            ref =self.option('seq').prop['path']
        elif self.option('seq_dir').is_set:
            ref =self.ref
        with open(self.file,'r') as f:
            lines = f.readlines()
            input = ' '.join(lines)
        cmd = '%s -cvg -p 20 -i %s -refsingle %s -o %s -depthsingle %s' %(self.soap_coverage,input,ref,detail,self.depth)
        command = self.add_command("run_coverage", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_coverage运行完成")
        else:
            self.set_error("run_coverage运行出错!", code="31401804")

    def run_depth_gc(self):
        wind_arrry = (self.option("windl")).split(',')
        self.logger.info("run_depth_gc运行开始")
        ref = ''
        if self.option('seq').is_set:
            ref = self.option('seq').prop['path']
        elif self.option('seq_dir').is_set:
            ref = self.ref
        for wind in wind_arrry:
            out_dir = self.work_dir + "/" + 'depth_gc_' + str(wind) + "/"
            cmd = '{} {}gc_depth_v2.0.pl {} {} -windl {} -step 500 -outdir {}'.format(self.perl_path,self.perl_script,ref,self.depth,wind,out_dir)
            gc_depth =  'depth_gc_' + str(wind)
            command = self.add_command(gc_depth, cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info(gc_depth + "运行完成")
            else:
                self.set_error('%s运行出错!', variables=(gc_depth), code="31401805")
            cmd_1 = '%s %s' % (self.R_path, out_dir + "depth_gc.cmd.r")
            self.logger.info(cmd_1)
            gc_depth_r = 'depth_gc_' + str(wind) +'_r'
            command1 = self.add_command(gc_depth_r, cmd_1)
            command1.run()
            self.wait(command1)
            if command1.return_code == 0:
                self.logger.info("R程序计算depth_gc成功")
            else:
                self.set_error("R程序计算depth_gc失败", code="31401806")
                self.set_error("R程序计算depth_gc失败", code="31401807")
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
        super(GcDepthTool, self).run()
        if self.option('seq_dir').is_set:
            self.run_ref_cat()
            self.run_ref_index()
            self.run_soap()
            self.run_coverage()
            self.run_depth_gc()
            self.set_output()
            self.end()
        elif self.option('seq').is_set:
            self.run_ref_index()
            self.run_soap()
            self.run_coverage()
            self.run_depth_gc()
            self.set_output()
            self.end()


