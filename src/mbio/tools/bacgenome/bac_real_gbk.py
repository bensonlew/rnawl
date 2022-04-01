#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re
import time
import subprocess
import shutil


class BacRealGbkAgent(Agent):
    """
    用于gbk文件生成
    last_modify: 2018.03.29
    之所以修改此脚本的原因是：应产品线第三次升级的需求，gff文件、ptt文件的结果与原来的一致
    且如果是完成图，生成合并的gff文件，为打通小工具使用
    """

    def __init__(self, parent):
        super(BacRealGbkAgent, self).__init__(parent)
        options = [
            {"name": "gen_gff", "type": "infile", "format": "gene_structure.gff3"},  #
            {"name": "rrna_gff", "type": "infile", "format": "gene_structure.gff3"},  #
            {"name": "trna_gff", "type": "infile", "format": "gene_structure.gff3"},  #
            {"name": "pro_fa", "type": "infile", "format": "sequence.fasta"},  #
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},  #
            {"name": "anno", "type": "infile", "format": "sequence.profile_table"},  #
            {"name": "sample_name", "type": "string"},
            {"name": "analysis", "type": "string", "default": "uncomplete"},  ###流程分析模式complete，uncomplete
            {"name": "gbk_dir", "type": "outfile", "format": "gene_structure.gbk_dir"},  # 输出文件
            {"name": "all_gbk", "type": "outfile", "format": "gene_structure.gbk"},  # 输出文件
        ]
        self.add_option(options)
        self.list =[]

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('gen_gff').is_set:
            raise OptionError("请设置基因组基因预测基因gff文件！", code="31400401")
        if not self.option('rrna_gff').is_set:
            raise OptionError("请设置基因组rRNA预测基因gff文件！", code="31400402")
        if not self.option('trna_gff').is_set:
            raise OptionError("请设置基因组tRNA预测基因gff文件！", code="31400403")
        if not self.option('pro_fa').is_set:
            raise OptionError("请设置基因组基因蛋白序列文件！", code="31400404")
        if not self.option('genome_fa').is_set:
            raise OptionError("请设置基因组组装序列文件！", code="31400405")
        if not self.option('anno').is_set:
            raise OptionError("请设置基因组基因注释汇总表！", code="31400406")
        if not self.option('analysis'):
            raise OptionError("请提供分析流程类型！", code="31400407")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["./fastq_stat.xls", "xls", "fastq信息统计表"]
        ])
        super(BacRealGbkAgent, self).end()


class BacRealGbkTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(BacRealGbkTool, self).__init__(config)
        self.genegff =self.option('gen_gff').prop['path']
        self.rrnagff = self.option('rrna_gff').prop['path']
        self.trnagff = self.option('trna_gff').prop['path']
        self.genomefa = self.option('genome_fa').prop['path']
        self.pro = self.option('pro_fa').prop['path']
        self.anno = self.option('anno').prop['path']
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/"
        self.list = []
        self.python = self.config.SOFTWARE_DIR + '/program/Python/bin/python'
        self.python_scrit = self.config.PACKAGE_DIR + "/bacgenome/gbk_rm_ty_tbl.py"
        self.tbl2asn = "bioinfo/Genomic/Sofware/tbl2asn/linux64.tbl2asn"
        self.template_sbt = self.config.PACKAGE_DIR + "/bacgenome/template.sbt"
        self.all_gbk = self.work_dir + '/' +self.option('sample_name') + '.gbk'


    def run_tbl(self):
        cmd = "{} {}tbl_generation.pl uncomplete {} {} {} {} {} {} ".format(self.perl_path, self.perl_script,self.genegff,self.trnagff,self.rrnagff,self.pro ,self.genomefa,self.anno)
        self.logger.info(cmd)
        self.logger.info("开始运行run_tbl")
        command = self.add_command("run_tbl", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_tbl完成")
        else:
            self.set_error("运行run_tbl运行出错!", code="31400401")

    def run_gbk(self):
        self.get_list()
        for file in self.list:
            fsa =self.work_dir + '/' + file + '/' + file + '.fsa'
            tbl = self.work_dir + '/' + file + '/' + file + '.tbl'
            ##guanqing.zou 20180904
            with open(tbl) as f:
                lines = f.readlines()
            if len(lines)==1:
                cmd = "{} -t {} -i {}  -V vb".format(self.tbl2asn, self.template_sbt, fsa)
            else:
                cmd = "{} -t {} -i {} -k {} -V vb".format(self.tbl2asn, self.template_sbt, fsa, tbl)

            self.logger.info(cmd)
            path = file.lower() + '_gbk_run'
            self.logger.info("开始运行%s" %path)
            # command = self.add_command(path, cmd)
            # command.run()
            # self.wait(command)
            # if command.return_code in [0,1]:  #1  是This copy of tbl2asn is more than a year old. 的返回码，有结果。
            try:
                subprocess.check_output(self.config.SOFTWARE_DIR + '/'+cmd, shell=True)
                self.logger.info("运行%s完成" % path)
            except Exception as e:
            #else:
                #self.set_error("运行%s运行出错!" , variables=( path), code="31400402")
                self.logger.info(e)
            os.system("sed  -i '/^LOCUS/s/circular//' %s"%(self.work_dir + '/' + file + '/' + file + '.gbf'))
            os.system("sed  -i '/^LOCUS/s/linear//' %s"%(self.work_dir + '/' + file + '/' + file + '.gbf'))
            os.system('sed \'s/ACCESSION/ACCESSION   %s/g\' %s -i' % (file, self.work_dir + '/' + file + '/' + file + '.gbf'))
            os.system('sed \':a;N;$!ba;s/-\\n[ ]*/M/g\' %s -i' % ( self.work_dir + '/' + file + '/' + file + '.gbf'))
            os.system('%s %s %s %s %s' % (self.python, self.python_scrit, self.work_dir + '/' + file + '/' + file + '.tbl', self.work_dir + '/' + file + '/' + file + '.gbf',
                                              self.work_dir + '/' + file + '/' + file + '.gbk'))
            os.system('cat %s >> %s' % (self.work_dir + '/' + file + '/' + file + '.gbk', self.all_gbk))
            self.logger.info("完成运行%s" %path)

    def run_tbl_complete(self):
        cmd = "{} {}tbl_generation.pl complete {} {} {} {} {} {} ".format(self.perl_path, self.perl_script,self.genegff,self.trnagff,self.rrnagff,self.pro ,self.genomefa,self.anno)
        self.logger.info(cmd)
        self.logger.info("开始运行run_tbl")
        command = self.add_command("run_tbl", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_tbl完成")
        else:
            self.set_error("运行run_tbl运行出错!", code="31400403")

    def run_gbk_complete(self):
        self.get_list()
        for file in self.list:
            fsa =self.work_dir + '/' + file + '/' + file + '.fsa'
            tbl = self.work_dir + '/' + file + '/' + file + '.tbl'
            cmd = "{} -t {} -i {} -k {} -V vb".format(self.tbl2asn, self.template_sbt, fsa, tbl)
            self.logger.info(cmd)
            path = file.lower() + '_gbk_run'
            self.logger.info("开始运行%s" %path)
            # command = self.add_command(path, cmd)
            # command.run()
            # self.wait(command)
            #if command.return_code in [0, 1] :
            try:
                subprocess.check_output(self.config.SOFTWARE_DIR + '/'+ cmd, shell=True)
                self.logger.info("运行%s完成" % path)
            # else:
            except Exception as e:
                #self.set_error("运行%s运行出错!" , variables=( path), code="31400404")
                self.logger.info(e)
            os.system("sed  -i '/^LOCUS/s/circular//' %s"%(self.work_dir + '/' + file + '/' + file + '.gbf'))
            os.system("sed  -i '/^LOCUS/s/linear//' %s"%(self.work_dir + '/' + file + '/' + file + '.gbf'))
            os.system('sed \'s/ACCESSION/ACCESSION   %s/g\' %s -i' %(file,self.work_dir + '/' + file + '/' + file + '.gbf'))
            os.system('sed \':a;N;$!ba;s/-\\n[ ]*/M/g\' %s -i' % (self.work_dir + '/' + file + '/' + file + '.gbf'))
            os.system('%s %s %s %s %s' %(self.python, self.python_scrit, self.work_dir + '/' + file + '/' + file + '.tbl', self.work_dir + '/' + file + '/' + file + '.gbf', self.work_dir + '/all.gbk'))
            os.system('cp %s %s' % (self.work_dir + '/all.gbk', self.work_dir + '/' + file + '/' + file + '.gbf'))
            self.logger.info("完成运行%s" %path)

    def run_bio_gbk_un(self):
        if os.path.exists(self.work_dir + 'all.antismash.gbk'):
            os.remove(self.work_dir + 'all.antismash.gbk')
        cmd = "{} {}GBK_generation_un.pl {} {} {} {} {} {} ".format(self.perl_path, self.perl_script,
                                                                            self.genegff, self.trnagff, self.rrnagff,
                                                                            self.pro, self.genomefa, self.anno)
        self.logger.info(cmd)
        self.logger.info("开始运行run_bio_gbk")
        command = self.add_command("run_bio_gbk", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_bio_gbk完成")
        else:
            self.set_error("运行run_bio_gbk运行出错!", code="31400405")

    def run_bio_gbk(self):
        if os.path.exists(self.work_dir + 'all.antismash.gbk'):
            os.remove(self.work_dir + 'all.antismash.gbk')
        cmd = "{} {}GBK_generation.pl {} {} {} {} {} {} ".format(self.perl_path, self.perl_script,
                                                                            self.genegff, self.trnagff, self.rrnagff,
                                                                            self.pro, self.genomefa, self.anno)
        self.logger.info(cmd)
        self.logger.info("开始运行run_bio_gbk")
        command = self.add_command("run_bio_gbk", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_bio_gbk完成")
        else:
            self.set_error("运行run_bio_gbk运行出错!", code="31400406")

    def get_list(self):
        with open(self.genomefa,'r') as f:
            lines=f.readlines()
            for line in lines:
                if re.search(r'^>',line):
                    line = line.replace('>', '')
                    tmp = line.rstrip('\n').split(' ')[0]
                    print(tmp)
                    self.list.append(tmp)


    def set_output(self):
        self.logger.info("set output")
        if self.option('analysis') in ['uncomplete']:
            if os.path.exists(self.output_dir + '/seq_gbk'):
                shutil.rmtree(self.output_dir + '/seq_gbk')
            os.mkdir(self.output_dir + '/seq_gbk')
            if os.path.exists(self.output_dir + '/analysis_gbk'):
                shutil.rmtree(self.output_dir + '/analysis_gbk')
            os.mkdir(self.output_dir + '/analysis_gbk')
            if os.path.exists(self.output_dir + '/gbk'):
                shutil.rmtree(self.output_dir + '/gbk')
            os.mkdir(self.output_dir + '/gbk')
            if os.path.exists(self.output_dir + '/gbk/Scaffold'):
                shutil.rmtree(self.output_dir + '/gbk/Scaffold')
            os.mkdir(self.output_dir + '/gbk/Scaffold')
            for file in self.list:
                if os.path.exists(self.output_dir + '/seq_gbk/' + file + '.gbk'):
                    os.remove(self.output_dir + '/seq_gbk/' + file + '.gbk')
                os.link(self.work_dir + '/gbk/' + file + '.gbk',
                            self.output_dir + '/seq_gbk/' + file + '.gbk')
            self.option('gbk_dir').set_path(self.output_dir + '/seq_gbk')
            if os.path.exists(self.output_dir + '/analysis_gbk/' + 'all.antismash.gbk'):
                os.remove(self.output_dir + '/analysis_gbk/' + 'all.antismash.gbk')
            os.link(self.work_dir + '/all.antismash.gbk', self.output_dir + '/analysis_gbk/' + 'all.antismash.gbk')
            self.option('all_gbk').set_path(self.output_dir + '/analysis_gbk/' + 'all.antismash.gbk')
            if os.path.exists(self.output_dir + '/gbk/Scaffold/' + self.option('sample_name') + '.gbk'):
                os.remove(self.output_dir + '/gbk/Scaffold/' + self.option('sample_name') + '.gbk')
            os.link(self.all_gbk, self.output_dir + '/gbk/Scaffold/' + self.option('sample_name') + '.gbk')

        elif self.option('analysis') in ['complete']:
            if os.path.exists(self.output_dir + '/analysis_gbk'):
                shutil.rmtree(self.output_dir + '/analysis_gbk')
            os.mkdir(self.output_dir + '/analysis_gbk')
            if os.path.exists(self.output_dir + '/seq_gbk'):
                shutil.rmtree(self.output_dir + '/seq_gbk')
            os.mkdir(self.output_dir + '/seq_gbk')
            for file in self.list:
                if not os.path.exists(self.output_dir + '/seq_gbk' + '/' + file):
                    os.mkdir(self.output_dir + '/seq_gbk' + '/' + file)
                if os.path.exists(self.output_dir + '/seq_gbk/' + file+'/'+ self.option('sample_name') + '_' + file + '.gbk'):
                    os.remove(self.output_dir + '/seq_gbk/' + file + '/' + self.option('sample_name') + '_' + file + '.gbk')
                os.link(self.work_dir + '/' + file + '/' + file + '.gbf',
                            self.output_dir + '/seq_gbk/' + file + '/' + self.option('sample_name') + '_' + file + '.gbk')
                if os.path.exists(self.output_dir + '/analysis_gbk/' +  file + '.gbk'):
                    os.remove(self.output_dir + '/analysis_gbk/' +  file + '.gbk')
                os.link(self.work_dir + '/gbk/' + file + '.gbk',self.output_dir + '/analysis_gbk/' + file + '.gbk')
            self.option('gbk_dir').set_path(self.output_dir + '/analysis_gbk')

        ##链接gff和ptt的结果
        if self.option("analysis") == 'uncomplete':
            if os.path.exists(self.output_dir + '/'+  self.option('sample_name') + '.gff'):
                os.remove(self.output_dir + '/'+  self.option('sample_name') + '.gff')
            os.link(self.work_dir + '/gff/out.gff',self.output_dir + '/'+ self.option('sample_name') + '.gff')

            if os.path.exists(self.output_dir + '/'+  self.option('sample_name') + '.ptt'):
                os.remove(self.output_dir + '/'+  self.option('sample_name') + '.ptt')
            os.link(self.work_dir + '/ptt/out.ptt',self.output_dir + '/'+  self.option('sample_name') + '.ptt')
            if os.path.exists(self.output_dir + '/'+  self.option('sample_name') + '.all.gff'):
                os.remove(self.output_dir + '/'+  self.option('sample_name') + '.all.gff')
            if os.path.exists(self.work_dir + '/gff/out.gff'):
                os.link(self.work_dir + '/gff/out.gff', self.output_dir + '/'+  self.option('sample_name') + '.all.gff')

        else:
            map = {
                'Chr':'Chromosome','Chr2':'Chromosome2',
                "Chr3":'Chromosome3','Chr4':'Chromosome4',
                "p":'Plasmid',"pA":'PlasmidA',"pB":'PlasmidB',
                "pC":'PlasmidC','pD':'PlasmidD','pE':'PlasmidE',
                "pF":'PlasmidF','pG':'PlasmidG','pH':'PlasmidH'
                }
            gff_files =os.listdir(self.work_dir + '/gff')
            for f in gff_files:
                pre = f.rstrip('.gff')
                if pre in map.keys():
                    new_f = map[pre]+'.gff'
                else:
                    new_f = f
                target_f =  self.output_dir + '/'+  self.option('sample_name') + '_' + new_f
                if os.path.exists(target_f):
                    os.remove(target_f)
                os.link(self.work_dir+'/gff/'+f, target_f)

            ptt_files =os.listdir(self.work_dir + '/ptt')
            for f in ptt_files:
                pre = f.split('.ptt')[0]
                if pre in map.keys():
                    new_f = map[pre]+'.ptt'
                else:
                    new_f = f
                target_f =  self.output_dir + '/'+  self.option('sample_name') + '_' + new_f
                if os.path.exists(target_f):
                    os.remove(target_f)
                os.link(self.work_dir+'/ptt/'+f, target_f)
            if os.path.exists(self.output_dir + '/'+  self.option('sample_name') + '.all.gff'):
                pass
            else:
                self.merge_gff()

    def run_gff_ptt(self):  # zouguanqing 20190409
        if self.option("analysis") == 'complete':
            gff_src = "generation_gff_complete.pl"
            ptt_src = "generation_ptt_complete.pl"
        else:
            gff_src = "generation_gff.pl"
            ptt_src = "generation_ptt.pl"

        cmd = "{} {}{} {} {} {} {} {}".format(self.perl_path, self.perl_script, gff_src,
                                                            self.genegff, self.trnagff, self.rrnagff,self.anno, 'gff')
        self.logger.info(cmd)
        self.logger.info("开始生成gff")
        command = self.add_command("run_bio_gff", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gff生成完成")
        else:
            self.set_error("gff生成出错!")

        cmd1 = "{} {}{} {} {} {}".format(self.perl_path, self.perl_script, ptt_src, self.genegff, self.anno, "gff")

        self.logger.info(cmd1)
        self.logger.info("开始生成ptt")
        command1 = self.add_command("run_bio_ptt", cmd1)
        command1.run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("ptt生成完成")
        else:
            self.set_error("ptt生成出错!")

    def merge_gff(self):
        """
        如果是完成图，需要将完成图的gff文件合并起来
        """
        gff_files = os.listdir(self.work_dir + '/gff')
        target_f =  self.output_dir + '/'+  self.option('sample_name') + '.all.gff'
        if os.path.exists(target_f):
            os.remove(target_f)
        allf = open(target_f, 'w')
        allf.writelines("##gff-version 3\n")
        for f in gff_files:
            gff_path = os.path.join(self.work_dir, "gff" ,f)
            for line in open(gff_path):
                if not line.startswith("#"):
                    allf.writelines(line)
        allf.close()

    def run(self):
        """
        运行
        """
        super(BacRealGbkTool, self).run()
        if self.option('analysis') in ['uncomplete']:
            self.run_tbl()
            self.run_gbk()
            self.run_bio_gbk_un()
            self.run_gff_ptt()
            time.sleep(5)
            self.set_output()
            self.end()
        elif self.option('analysis') in ['complete']:
            self.run_tbl_complete()
            self.run_gbk_complete()
            self.run_bio_gbk()
            self.run_gff_ptt()
            self.merge_gff()
            time.sleep(5)
            self.set_output()
            self.end()

