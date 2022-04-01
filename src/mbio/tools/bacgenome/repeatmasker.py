# -*- coding: utf-8 -*-
# last_modify: 2020.10.23

import os
import re
import glob
import subprocess
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from operator import itemgetter


class RepeatmaskerAgent(Agent):
    """
    RepeatMasker 查找散在重复序列(interspersed repeats)
    """

    def __init__(self, parent):
        super(RepeatmaskerAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},
            {"name": "gff", "type": "outfile", "format": "gene_structure.gff3"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome", code="33300701")

    def set_resource(self):
        self._cpu = 5
        self._memory = '20G'

    def end(self):
        super(RepeatmaskerAgent, self).end()

class RepeatmaskerTool(Tool):
    def __init__(self, config):
        super(RepeatmaskerTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        self.sample_name = (self.genome_fasta).split('/')[-1].rsplit('.')[0]
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/sequence/scripts/"
        self.repeatmasker_path = "/bioinfo/Genomic/Sofware/RepeatMasker/"
        perl5lib_path="{0}/program/perl-5.24.0/lib/site_perl/5.24.0:{0}/program/perl/perls/perl-5.24.0/lib:{0}/program/perl/perls/perl-5.24.0/lib/5.24.0/x86_64-linux-thread-multi:{0}/program/perl/perls/perl-5.24.0/lib/site_perl/5.24.0/x86_64-linux-thread-multi:{0}/program/perl/perls/perl-5.24.0/lib".format(self.config.SOFTWARE_DIR)
        self.set_environ(PERL5LIB = perl5lib_path)
        self.set_environ(PATH = "{0}/program/perl-5.24.0/bin".format(self.config.SOFTWARE_DIR))

    def run_repeatmasker(self):
        cmd = "{}RepeatMasker -parallel 5 -engine ncbi -nolow -no_is -norna -dir {}  {}".format(
                self.repeatmasker_path, self.work_dir, self.genome_fasta)
        command = self.add_command("repeatpredict", cmd).run()
        self.wait(command)
        self.logger.info(command.return_code)
        if command.return_code == 0:
            self.logger.info("Repeatmasker运行完成")
            self.run_repeat_to_gff()
        else:
            self.set_error("Repeatmasker运行出错!", code="33300701")

    def run_repeat_to_gff(self):
        if os.path.exists(self.work_dir + "/" + (self.genome_fasta).split('/')[-1] + ".tbl"):
            file = (self.genome_fasta).split('/')[-1] + ".out"
            if os.path.exists(self.work_dir + "/" + self.sample_name + ".out"):
                os.remove(self.work_dir + "/" + self.sample_name + ".out")
            os.link(self.work_dir + "/" + file, self.work_dir + "/" + self.sample_name + ".out")
            if os.path.exists(self.work_dir + "/" + self.sample_name + ".tbl"):
                os.remove(self.work_dir + "/" + self.sample_name + ".tbl")
            os.link(self.work_dir + "/" + (self.genome_fasta).split('/')[-1] + ".tbl",
                self.work_dir + "/" + self.sample_name + ".tbl")
            file = self.work_dir + "/" + self.sample_name + ".out"
            cmd = '{} {}repeat_to_gff.pl {}'.format(self.perl_path, self.perl_script, file)
            self.logger.info(cmd)
            command = self.add_command("run_repeat_to_gff", cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("生成gff文件运行完成")
                self.set_output()
            else:
                self.set_error("生成gff文件运行出错!", code="33300705")
        else:
            pass

    def set_output(self):
        gff = self.work_dir + "/" + self.sample_name + ".out.gff"
        gff_path = self.output_dir + "/" + "Rep.gff"
        if os.path.exists(gff):
            if os.path.exists(gff_path):
                os.remove(gff_path)
            os.link(gff, gff_path)
            tbl = self.work_dir + "/" + self.sample_name + ".tbl"
            tbl_new_path = self.output_dir + "/" + "Rep.tbl"
            if os.path.exists(tbl_new_path):
                os.remove(tbl_new_path)
            os.link(tbl, tbl_new_path)
            out = self.work_dir + "/" + self.sample_name + ".out"
            out_new_path = self.output_dir + '/' + "Rep.out"
            if os.path.exists(out_new_path):
                os.remove(out_new_path)
            os.link(out, out_new_path)
            self.option('gff', gff_path)
        else:
            pass

    def run(self):
        super(RepeatmaskerTool, self).run()
        self.run_repeatmasker()
        self.end()
    '''
    def sort_out(self, infile ,newfile):
        compile = re.compile(r'{}'.format('\s*'))
        compile1 = re.compile(r'[S|s]caffold')
        compile2 = re.compile(r'[C|c]hromosome')
        compile3 = re.compile(r'[P|p]la')
        fr = open(infile)
        lines = fr.readlines()
        rlist = []
        for line in lines[3:]:
            line = line.strip()
            spl = compile.split(line)
            if compile1.match(spl[4]):
                spl.append(int(spl[4][8:])) #scaffold的编号
            elif compile2.match(spl[4]):
                spl.append(1)  #chr的编号
            elif compile3.match(spl[4]):
                spl.append(1)  #Plasmid编号
            else:
                self.set_error("out 文件格式不通过", code="33300706")
            spl.append(int(spl[5]))  #重复序列start
            rlist.append(spl)

        sl = sorted(rlist, key=itemgetter(-2,-1))  #根据编号和start 排序
        fw = open(newfile, 'w')
        fw.write(''.join(lines[:2]))
        id = 1
        for i in sl:
            fw.write(' '.join(i[:-4]) + ' ' + str(id) + '\n')
            id += 1
    '''

