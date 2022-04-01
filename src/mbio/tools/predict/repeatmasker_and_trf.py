# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# version 1.0
# last_modify: 2017.12.27

import os
import re
import glob
import subprocess
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from operator import itemgetter


class RepeatmaskerAndTrfAgent(Agent):
    """
    RepeatModeler 构建denovo数据库
    RepeatMasker 查找散在重复序列(interspersed repeats)
    TRF 寻找串联重复序列 (tandem repeat)
    """

    def __init__(self, parent):
        super(RepeatmaskerAndTrfAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的scaffold文件
            {"name": "analysis_type", "type": "string", "default": "1"},
            # 采用的trf预测方法默认trf:1，RepeatMasker：2，RepeatModeler：3
            {"name": "match", "type": "int", "default": 2},  # 选择1时， 匹配权重，设置为2
            {"name": "mismatch", "type": "int", "default": 7},  # 选择1时， 错配，设置为7
            {"name": "delta", "type": "int", "default": 7},  # 选择1时， 插入，设置为7
            {"name": "pm", "type": "int", "default": 80},  # 选择1时， match probability，设置为80
            {"name": "pi", "type": "int", "default": 500},  # 选择1时， indel probability，设置为500
            {"name": "gff", "type": "outfile", "format": "gene_structure.gff3"},  # 预测生成的GFF格式
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome", code="33300701")
        if self.option("analysis_type") not in ["1", "2", "3"]:
            raise OptionError("不存在该预测方法，请重新设置！", code="33300702")

    def set_resource(self):
        if self.option("analysis_type") == "1":
            self._cpu = 2
            self._memory = '2G'
        else:
            self._cpu = 5
            self._memory = '20G'

    def end(self):
        super(RepeatmaskerAndTrfAgent, self).end()


class RepeatmaskerAndTrfTool(Tool):
    def __init__(self, config):
        super(RepeatmaskerAndTrfTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        self.sample_name = (self.genome_fasta).split('/')[-1].rsplit('.')[0]
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/sequence/scripts/"
        self.repeatmodeler_path = "/bioinfo/Genomic/Sofware/RepeatModeler-open-1.0.10/"
        self.repeatmasker_path = "/bioinfo/Genomic/Sofware/RepeatMasker/"
        self.trf_path = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/"
        perl5lib_path="{0}/program/perl-5.24.0/lib/site_perl/5.24.0:{0}/program/perl/perls/perl-5.24.0/lib:{0}/program/perl/perls/perl-5.24.0/lib/5.24.0/x86_64-linux-thread-multi:{0}/program/perl/perls/perl-5.24.0/lib/site_perl/5.24.0/x86_64-linux-thread-multi:{0}/program/perl/perls/perl-5.24.0/lib".format(self.config.SOFTWARE_DIR)
        self.set_environ(PERL5LIB = perl5lib_path)  #add for nb zouguanqing 20181022
        self.set_environ(PATH = "{0}/program/perl-5.24.0/bin".format(self.config.SOFTWARE_DIR))

    def run_repeatmasker(self):
        if self.option("analysis_type") == "2":
            cmd = "{}RepeatMasker -parallel 5 -engine ncbi -nolow -no_is -norna -dir {}  {}".format(
                self.repeatmasker_path, self.work_dir, self.genome_fasta)
        elif self.option("analysis_type") == "3":
            repeatmodeler_consensi = glob.glob(self.work_dir + "/RM*/consensi.fa.classified")[0]
            cmd = "{}RepeatMasker -parallel 5 -lib {} -engine ncbi -nolow -no_is -norna -dir {}  {}".format(
                self.repeatmasker_path, repeatmodeler_consensi, self.work_dir,
                self.genome_fasta)
        command = self.add_command("repeatpredict", cmd).run()
        self.wait(command)
        self.logger.info(command.return_code)
        if command.return_code == 0:
            self.logger.info("Repeatmasker运行完成")
            self.run_repeat_to_gff()
        else:
            self.set_error("Repeatmasker运行出错!", code="33300701")

    def run_repeatmodeler(self):
        cmd1 = "{}BuildDatabase -engine ncbi -name {} {}".format(self.repeatmodeler_path, self.sample_name,
                                                                 self.genome_fasta)
        cmd2 = "{}RepeatModeler -engine ncbi -database {} ".format(self.repeatmodeler_path, self.sample_name)
        command = self.add_command("build database", cmd1).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("build database运行完成")
            command1 = self.add_command("build database2", cmd2).run()
            self.wait(command1)
            if command1.return_code == 0:
                self.logger.info("build database运行完成")
                self.run_repeatmasker()
            else:
                self.logger.info("RepeatModeler运行完成")
        else:
            self.set_error("BuildDatabase运行出错!", code="33300702")

    def run_trf(self):
        cmd = "{}trf  {} {} {} {} 80 10 {} {} -d -h -ngs > {}".format(self.trf_path, self.genome_fasta,
                                                                      self.option("match"), self.option("mismatch"),
                                                                      self.option("delta"), self.option("pm"),
                                                                      self.option("pi"), self.sample_name + ".TRF.dat")
        # command = self.add_command("repeat_predict", cmd).run()
        # self.wait(command)
        # self.logger.info(command.return_code)
        # if command.return_code == 0:
        #    self.logger.info("TRF运行完成")
        #    self.run_repeat_to_gff()
        # else:
        #    self.set_error("TRF运行:没有预测到repeat!!")
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info("TRF运行完成")
            self.run_repeat_to_gff()
        except subprocess.CalledProcessError:
            self.set_error("TRF运行失败!!", code="33300703")

    def run_repeat_to_gff(self):
        if self.option("analysis_type") == "1":
            file = self.sample_name + ".TRF.dat"
        elif self.option("analysis_type") in ["2", "3"]:
            file = (self.genome_fasta).split('/')[-1] + ".out"
            os.link(self.work_dir + "/" + file, self.work_dir + "/" + self.sample_name + ".out.old")
            os.link(self.work_dir + "/" + (self.genome_fasta).split('/')[-1] + ".tbl",
                    self.work_dir + "/" + self.sample_name + ".tbl")
            file_old = self.work_dir + "/" + self.sample_name + ".out.old"
            file = self.work_dir + "/" + self.sample_name + ".out"
            self.sort_out(file_old, file)   #zouguanqing 20180706 原.out文件没有排序，改名为out.old。排序生成.out文件
        else:
            self.logger.info('分析类型设置值不是1,2或3')
            self.set_error('分析类型设置值不是1,2或3', code="33300704")

        cmd = '{} {}repeat_to_gff.pl {}'.format(self.perl_path, self.perl_script, file)
        self.logger.info(cmd)
        command = self.add_command("run_repeat_to_gff", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("生成gff文件运行完成")
            self.set_output()
        else:
            self.set_error("生成gff文件运行出错!", code="33300705")

    def set_output(self):
        gff_path = ""
        if self.option("analysis_type") in ["1"]:
            gff = self.work_dir + "/" + self.sample_name + ".TRF.dat.gff"
            gff_path = self.output_dir + "/" + self.sample_name + ".TRF.gff"
            if os.path.exists(gff_path):   #guanqing.zou 20180820
                os.remove(gff_path)
            os.link(gff, gff_path)
            if os.path.exists(self.output_dir + "/" + self.sample_name + ".TRF.dat"):    #guanqing.zou 20180820
                os.remove(self.output_dir + "/" + self.sample_name + ".TRF.dat")
            os.link(self.work_dir + "/" + self.sample_name + ".TRF.dat.new", self.output_dir + "/" + self.sample_name + ".TRF.dat")
        if self.option("analysis_type") in ["2", "3"]:
            gff = self.work_dir + "/" + self.sample_name + ".out.gff"
            gff_path = self.output_dir + "/" + self.sample_name[:-5] + "_Rep.gff"
            if os.path.exists(gff_path):
                os.remove(gff_path)
            os.link(gff, gff_path)
            tbl = self.work_dir + "/" + self.sample_name + ".tbl"
            tbl_new_path = self.output_dir + "/" + self.sample_name[:-5] + "_Rep.tbl"
            if os.path.exists(tbl_new_path):
                os.remove(tbl_new_path)
            os.link(tbl, tbl_new_path)
            out = self.work_dir + "/" + self.sample_name + ".out"
            out_new_path = self.output_dir + '/' + self.sample_name[:-5] + "_Rep.out"
            if os.path.exists(out_new_path):
                os.remove(out_new_path)
            os.link(out, out_new_path)
        self.option('gff', gff_path)

    def run(self):
        super(RepeatmaskerAndTrfTool, self).run()
        if self.option("analysis_type") == "1":
            self.run_trf()
        elif self.option("analysis_type") == "2":
            self.run_repeatmasker()
        elif self.option("analysis_type") == "3":
            self.run_repeatmodeler()
        self.end()


    #zouguanqing 20180706 根据scaffold编号和start 排序
    def sort_out(self, infile ,newfile):
        compile = re.compile(r'{}'.format('\s*'))
        compile1 = re.compile(r'[S|s]caffold')
        fr = open(infile)
        lines = fr.readlines()
        rlist = []
        for line in lines[3:]:
            line = line.strip()
            spl = compile.split(line)
            if compile1.match(spl[4]):
                spl.append(int(spl[4][8:])) #scaffold的编号
            else:
                self.set_error("out 文件格式不通过", code="33300706")
            spl.append(int(spl[6]))  #重复序列start
            rlist.append(spl)

        sl = sorted(rlist, key=itemgetter(-2,-1))  #根据scaffold编号和start 排序
        fw = open(newfile, 'w')
        fw.write(''.join(lines[:2]))
        id = 1
        for i in sl:
            fw.write(' '.join(i[:-4]) + ' ' + str(id) + '\n')
            id += 1


