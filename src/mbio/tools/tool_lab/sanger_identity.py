#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'zzg'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os, shutil
from Bio import SeqIO
import glob
import subprocess


class SangerIdentityAgent(Agent):
    """
    用于一代测序菌种鉴定
    version 1.0
    author: zzg
    last_modify: 2021.4.14
    """

    def __init__(self, parent):
        super(SangerIdentityAgent, self).__init__(parent)
        options = [
            {"name": "file_dir", "type": "infile", 'format': "denovo_rna_v2.common_dir"},
            {"name": "method", "type": "string", "default": "false"} # NT数据库比对,默认否
        ]
        self.add_option(options)
        self.list = []

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option("file_dir").is_set:
            raise OptionError("必须输入文件！")
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 16
        self._memory = '150G'

    def end(self):
        super(SangerIdentityAgent, self).end()


class SangerIdentityTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(SangerIdentityTool, self).__init__(config)
        self.Rscript = "program/miniconda3/envs/R_v4/bin/Rscript"  #  dev: bioinfo/tool_lab/miniconda2/envs/R4version/bin/Rscript
        self.software_dir = self.config.SOFTWARE_DIR
        self.perl_script = self.config.PACKAGE_DIR + "/tool_lab/bac_identity/gi_tax.pl"
        self.set_environ(
            LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/library/glibc-2.14/lib/')
        self.perl_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/perl"
        self.perl5path = self.config.SOFTWARE_DIR + "/miniconda2/bin/perl"
        self.set_environ(PERL5LIB=self.perl5path)
        self.script = self.config.PACKAGE_DIR + '/tool_lab/sanger_identity.r'
        self.assemble_sh = self.config.PACKAGE_DIR + '/tool_lab/assemble.sh'
        self.assemble_script = self.config.SOFTWARE_DIR + '/bioinfo/tool_lab/pairASS/'
        self.blastn = "bioinfo/align/ncbi-blast-2.10.1+/bin/blastn"
        self.perl = 'program/perl-5.24.0/bin/perl'
        self.nt = os.path.join(self.software_dir, "database/Tool_lab/bac_identity/db_blast/nt")
        self.nt_tax = os.path.join(
            self.software_dir, "database/Tool_lab/bac_identity/NT_tax.xls")
        self.readme = os.path.join(self.software_dir, "database/Tool_lab/bac_identity/README.txt")

    def get_info(self):
        if os.path.exists(self.work_dir + "/tmp_dir"):
            shutil.rmtree(self.work_dir + "/tmp_dir")
        os.mkdir(self.work_dir + "/tmp_dir")
        all_file = os.listdir(self.option("file_dir").prop["path"])
        self.list_file = self.work_dir + "/tmp_dir/list.txt"
        with open(self.option("file_dir").prop["path"] + "/list.txt","r") as f,open(self.list_file,"w") as t:
            t.write("\t"+"样品名称"+"\t"+"测序引物"+"\t"+"测序状态"+"\t"+"测序结果"+"\n")
            data = f.readlines()
            self.sample = {}
            for i in data[1:]:
                if i.split("\t")[0] not in self.sample:
                    for ii in all_file:
                        if ii.endswith("ab1"):
                            if i.split("\t")[0] in ii and i.split("\t")[1] in ii.split("-")[-1]:
                                self.sample[i.split("\t")[0]] = [{i.split("\t")[1]: ii}]
                                os.link(self.option("file_dir").prop["path"] + "/" + ii,
                                        self.work_dir + "/tmp_dir/" + ii)
                else:
                    for ii in all_file:
                        if ii.endswith("ab1"):
                            if i.split("\t")[0] in ii and i.split("\t")[1] in ii.split("-")[-1]:
                                self.sample[i.split("\t")[0]].append({i.split("\t")[1]: ii})
                                os.link(self.option("file_dir").prop["path"] + "/" + ii,
                                        self.work_dir + "/tmp_dir/" + ii)
            num =1
            for x in self.sample:
                t.write(str(num) + "\t" + x + "\t" + self.sample[x][0].keys()[0] + "\t" + "测序完成" + "\t" + "正常" + "\n")
                num +=1
                t.write(str(num) + "\t" + x + "\t" + self.sample[x][1].keys()[0] + "\t" + "测序完成" + "\t" + "正常" + "\n")

    def run_identity(self):
        os.chdir(self.work_dir + "/tmp_dir")
        cmd = '{} {}'.format(self.Rscript, self.script)
        cmd += " {}".format(self.list_file)
        self.logger.info(cmd)
        self.logger.info("开始运行run_identity")
        command = self.add_command("run_identity", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_identity完成")
        else:
            self.set_error("运行run_identity运行出错!")

    def get_all_contig(self):
        #os.system('cd {}'.format(self.work_dir + "/tmp_dir"))
        #os.system(self.assemble_sh + " " + self.software_dir + "/bioinfo/tool_lab/pairASS/")
        #self.tmp_sh = "{} {}".format(self.assemble_sh, self.software_dir + "/bioinfo/tool_lab/pairASS/")
        #cmd1 = '/program/sh {}'.format(self.tmp_sh)
        #self.logger.info('start')
        #command1 = self.add_command('assemble_sh', cmd1).run()
        #self.logger.info('running')
        #self.wait(command1)
        #self.logger.info('end')
        #if command1.return_code == 0:
        #    self.logger.info("%s运行成功" % cmd1)
        #else:
        #    self.set_error("%s运行失败", variables=(cmd1))
        self.assemble(self.work_dir + "/tmp")
        allcontigs = open(os.path.join(self.work_dir, "all.contigs"), "w")
        sample_list = []
        sample_num = len(self.sample.keys())
        for i in self.sample:
            sample_list.append(i)
        file_list = []
        self.logger.info(sample_list)
        self.logger.info(sample_num)
        for root, dir_name, file in os.walk(self.work_dir + "/tmp" + "/contigs"):
            # if len(file) < sample_num:
            #    self.logger.info("assemble出现错误，组装的样本小于列表的样本")
            for name in file:
                if name == "readme.txt":
                    continue
                # self.logger.info(name)
                file_list.append(name.split('.')[0])
                # if name[:-4] not in sample_list:
                #     nofound_file.append(file[:-4])
                #     # self.set_error("样本{}没有生成".format(file))
            self.logger.info(file_list)
            nofound_file = list(set(sample_list) - set(file_list))
            # if len(nofound_file) > 0:
            #    self.set_error("样本{}没有生成".format(",".join(nofound_file)))
            for name in file:
                if name == "readme.txt":
                    os.remove(os.path.join(root, name))
                    continue
                with open(os.path.join(root, name), "r") as contig, open(
                        os.path.join(root, name.split(".")[0] + ".seq"), "w") as seq:
                    n = 0
                    text = ""
                    while 1:
                        line = contig.readline()

                        if not line:
                            break
                        if line[0] == ">":
                            text += ">{}\n".format(name.split(".")[0])
                            n += 1
                            continue
                        # seq.write(">{}\n".format(name.split(".")[0]))
                        # allcontigs.write(">{}\n".format(name.split(".")[0]))
                        text += line
                    if n > 1:
                        mj = name.split(".")[0]
                        primer = self.info[mj].primer[0] if float(self.info[mj].identity[0]) > float(
                            self.info[mj].identity[1]) else self.info[mj].primer[1]
                        fafile = glob.glob(os.path.join(self.fasta_dir + '/' + mj + '-' + primer + '*.fa'))[0]
                        fa_record = SeqIO.read(fafile, "fasta")
                        SeqIO.write(fa_record, self.fasta_dir + "/temp.txt", "fasta")
                        with open(self.fasta_dir + "/temp.txt", "r") as temp_seq:
                            text = temp_seq.read()
                    seq.write(text)
                allcontigs.write(text)
                os.remove(os.path.join(root, name))

    def assemble(self,path):
        """
        tmp 文件夹组装序列
        :param path:
        :return:
        """
        if os.path.exists(self.work_dir + "/tmp"):
            shutil.rmtree(self.work_dir + "/tmp")
        os.mkdir(self.work_dir + "/tmp")
        os.mkdir(self.work_dir + "/tmp/contigs")
        file_dict = {}
        sample_list = []
        for file in os.listdir(self.work_dir + "/tmp_dir"):
            if file.endswith(".fa"):
                with open(file) as f:
                    data = f.readlines()
                    for i in data:
                        if i.startswith(">"):
                            sample = i.strip().lstrip(">")
                if sample in file_dict:
                    file_dict[sample].append(file)
                else:
                    file_dict[sample] = [file]
        for sample in file_dict:
            if len(file_dict[sample]) == 2:
                os.mkdir(self.work_dir + "/tmp/" + sample)
                os.system("cat {} {} > {}".format(self.work_dir + "/tmp_dir/"+file_dict[sample][0],self.work_dir + "/tmp_dir/"+file_dict[sample][1],self.work_dir + "/tmp/" + sample + "/" + sample + ".fasta.screen"))
                os.system("{}/phrap {} -new_ace".format(self.assemble_script,self.work_dir + "/tmp/" + sample + "/" + sample + ".fasta.screen"))
                self.logger.info("{}/tagRepeats.perl {} {}".format(self.assemble_script,self.assemble_script,sample + ".fasta.screen.ace"))
                os.chdir(self.work_dir + "/tmp/" + sample + "/")
                #cmd = "export LC_ALL=C && {} {}/tagRepeats.perl {} {}".format(self.perl_path,self.assemble_script,self.assemble_script,sample + ".fasta.screen.ace")
                #subprocess.call(cmd, shell=True)
                os.system("export LC_ALL=C && {} {}/tagRepeats.perl {} {}".format(self.perl_path,self.assemble_script,self.assemble_script,sample + ".fasta.screen.ace"))

                if os.path.exists(self.work_dir + "/tmp/" + sample + "/" + sample + ".contigs"):
                    if os.path.getsize(self.work_dir + "/tmp/" + sample + "/" + sample + ".contigs") != 0:
                        os.system("{}/cross_match {} {} -tags -minmatch 10".format(self.assemble_script,self.work_dir + "/tmp/" + sample + "/" + sample + ".contigs", self.assemble_script + "/lib/screenLibs/repeats.fasta"))
                        shutil.copyfile(self.work_dir + "/tmp/" + sample + "/" + sample + ".contigs", self.work_dir + "/tmp/contigs/" + sample + ".contigs")

    def run_nt_blastn(self):
        """
        菌鉴流程的模块
        blastn -db nt -query contigs/all.contigs -out NT_blast.m8.txt -evalue 1e-5 -outfmt 6 -max_target_seqs 10
        perl /gi_tax.pl NT_blast.m8.txt /NT_tax.xls  NT.fasta.blast.taxonomy.xls
        """
        """
        with open(self.work_dir + "/all.contigs","w") as t:
            for ii in os.listdir(self.work_dir + "/tmp_dir"):
                for x in self.sample:
                    if ii.endswith("fa"):
                        if (x in ii) and (self.sample[x][0].keys()[0] in ii.split("-")[-1]):
                            with open(ii) as f:
                                data = f.readlines()
                                t.write(">" + str(ii) + "\n" + "\n".join(data[1:]))
                        elif (x in ii) and (self.sample[x][1].keys()[0] in ii.split("-")[-1]):
                            with open(ii) as f:
                                data = f.readlines()
                                t.write(">" + str(ii) + "\n" + "\n".join(data[1:]))
        """
        cmd = "{} -db {} -query {}".format(self.blastn, self.nt, os.path.join(self.work_dir, "all.contigs"))
        out = os.path.join(self.work_dir, "NT_blast.m8.txt")
        out_tax = os.path.join(self.work_dir, "NT.fasta.blast.taxonomy.xls")
        cmd += " -out {} -evalue 1e-5 -outfmt 6 -max_target_seqs 10 -num_threads 16".format(out)
        self.logger.info(cmd)
        command = self.add_command("blastn", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行blastn完成")
        else:
            self.logger.info("运行错误，请重新检查输入参数")
        cmd1 = "{} {}".format(self.perl, self.perl_script)
        cmd1 += " {} {} {}".format(out, self.nt_tax, out_tax)
        command1 = self.add_command("gitax", cmd1).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("运行gi_tax完成")
        else:
            self.set_error("运行错误，请重新检查输入参数")

    def get_result(self):
        if self.option("method") == "true":
            if self.contig_size != 0:
                with open(self.work_dir + "/NT_blast.m8.txt") as f, open(
                        self.work_dir + "/NT.fasta.blast.taxonomy.xls") as g, open(
                    self.output_dir + "/sangerSeq_NTblast.xls", "w") as t:
                    t.write(
                        "Sample Name\tSubject\tTaxonomy\tIdentity(%)\tAlign Len (bp)\tMismatches\tGap_Open\tQ_Start\tQ_End\tS_Start\tS_End\tEvalue\tBit Score\n")
                    data1 = f.readlines()
                    data2 = g.readlines()
                    sample_list1 = []
                    sample_list2 = {}
                    tax_dict = {}

                    for x in data2:
                        tax_dict[x.split("\t")[2]] = x.split("\t")[1]
                    for i in data1:
                        t.write(i.split("\t")[0] + "\t" + i.split("\t")[1] + "\t" + tax_dict[
                            i.split("\t")[1]] + "\t" + "\t".join(i.strip().split("\t")[2:]) + "\n")
                        if i.split("\t")[0] in sample_list1:
                            pass
                        else:
                            sample_list1.append(i.split("\t")[0])
                            sample_list2[i.split("\t")[0]] = tax_dict[i.split("\t")[1]]
            else:
                sample_list2 = {}

        with open(self.output_dir +"/sangerSeq_stat.xls","w") as v:
            if self.option("method") == "true":
                v.write("Sample Name\tF Quality\tF Raw Len (bp)\tF Identity(%)\tR Quality\tR Raw Len (bp)\tR Trim Len (bp)\tR Identity(%)\tNT Blast\n")
            else:
                v.write(
                    "Sample Name\tF Quality\tF Raw Len (bp)\tF Identity(%)\tR Quality\tR Raw Len (bp)\tR Trim Len (bp)\tR Identity(%)\n")
            for x in self.sample:
                for ii in os.listdir(self.work_dir + "/tmp_dir"):
                    if ii.endswith("report"):
                        if "F" in self.sample[x][0].keys()[0]:
                            if (x in ii) and (self.sample[x][0].keys()[0] in ii.split("-")[-1]):
                                with open(self.work_dir + "/tmp_dir/" + ii) as f:
                                    data = f.readlines()
                                    for xx in data:
                                        if xx.startswith("rawseqlength"):
                                            raw_len = xx.strip().split(" ")[1]
                                        elif xx.startswith("trimseqlength"):
                                            trim_len = xx.strip().split(" ")[-1]
                                        elif xx.startswith("identity"):
                                            identity = xx.strip().split(" ")[-1].strip("()")
                                    if identity >= "99.5":
                                        quality = "正常"
                                    else:
                                        quality = "污染"
                                    v.write(x+"\t"+quality+"\t"+raw_len+"\t"+trim_len+"\t"+identity+"\t")
                for ii in os.listdir(self.work_dir + "/tmp_dir"):
                    if ii.endswith("report"):
                        if "F" in self.sample[x][0].keys()[0]:
                            if (x in ii) and (self.sample[x][1].keys()[0] in ii.split("-")[-1]):
                                with open(self.work_dir + "/tmp_dir/" + ii) as f:
                                    data = f.readlines()
                                    for xx in data:
                                        if xx.startswith("rawseqlength"):
                                            raw_len = xx.strip().split(" ")[1]
                                        elif xx.startswith("trimseqlength"):
                                            trim_len = xx.strip().split(" ")[-1]
                                        elif xx.startswith("identity"):
                                            identity = xx.strip().split(" ")[-1].strip("()")
                                    if identity >= "99.5":
                                        quality = "正常"
                                    else:
                                        quality = "污染"
                                    if self.option("method") == "true":
                                        if x in sample_list2:
                                            v.write(quality + "\t" + raw_len + "\t" + trim_len + "\t" + identity + "\t" +sample_list2[x].split("s__")[-1] + "\n")
                                        else:
                                            v.write(quality + "\t" + raw_len + "\t" + trim_len + "\t" + identity + "\t" + "-" + "\n")
                                    else:
                                        v.write(quality+"\t"+raw_len+"\t"+trim_len+"\t"+identity+"\n")
                for ii in os.listdir(self.work_dir + "/tmp_dir"):
                    if ii.endswith("report"):
                        if "R" in self.sample[x][0].keys()[0]:
                            if (x in ii) and (self.sample[x][1].keys()[0] in ii.split("-")[-1]):
                                with open(self.work_dir + "/tmp_dir/" + ii) as f:
                                    data = f.readlines()
                                    for xx in data:
                                        if xx.startswith("rawseqlength"):
                                            raw_len = xx.strip().split(" ")[1]
                                        elif xx.startswith("trimseqlength"):
                                            trim_len = xx.strip().split(" ")[-1]
                                        elif xx.startswith("identity"):
                                            identity = xx.strip().split(" ")[-1].strip("()")
                                    if identity >= "99.5":
                                        quality = "正常"
                                    else:
                                        quality = "污染"
                                    v.write(x+"\t"+quality+"\t"+raw_len+"\t"+trim_len+"\t"+identity+"\t")
                for ii in os.listdir(self.work_dir + "/tmp_dir"):
                    if ii.endswith("report"):
                        if "R" in self.sample[x][0].keys()[0]:
                            if (x in ii) and (self.sample[x][0].keys()[0] in ii.split("-")[-1]):
                                with open(self.work_dir + "/tmp_dir/" + ii) as f:
                                    data = f.readlines()
                                    for xx in data:
                                        if xx.startswith("rawseqlength"):
                                            raw_len = xx.strip().split(" ")[1]
                                        elif xx.startswith("trimseqlength"):
                                            trim_len = xx.strip().split(" ")[-1]
                                        elif xx.startswith("identity"):
                                            identity = xx.strip().split(" ")[-1].strip("()")
                                    if identity >= "99.5":
                                        quality = "正常"
                                    else:
                                        quality = "污染"
                                    if self.option("method") == "true":
                                        if x in sample_list2:
                                            v.write(quality + "\t" + raw_len + "\t" + trim_len + "\t" + identity + "\t" + sample_list2[x].split("s__")[-1] +"\n")
                                        else:
                                            v.write(quality + "\t" + raw_len + "\t" + trim_len + "\t" + identity + "\t" + "-" + "\n")
                                    else:
                                        v.write(quality+"\t"+raw_len+"\t"+trim_len+"\t"+identity+"\n")

    def run(self):
        """
        运行
        """
        super(SangerIdentityTool, self).run()
        self.get_info()
        self.run_identity()
        if self.option("method") == "true":
            self.get_all_contig()
            self.contig_size = os.path.getsize(os.path.join(self.work_dir, "all.contigs"))
        else:
            self.contig_size = 0
        if self.contig_size != 0:
            self.run_nt_blastn()
        self.get_result()
        self.end()