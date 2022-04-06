# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
#from biocluster.config import config
from collections import namedtuple, defaultdict
from Bio import SeqIO
import time
import re
import time
import unittest
import os
import glob
import sys
import shutil


class BacIdentityAgent(Agent):
    def __init__(self, parent):
        super(BacIdentityAgent, self).__init__(parent)
        options = [
            {"name": "fasta_dir", "type": "infile", "format": "denovo_rna_v2.common_dir"},
            {"name": "threshold", "type": "float", "default": 99.5},
            {"name": "list_xls", "type": "infile", "format": "denovo_rna_v2.common"},
            {"name": "sample_info","type":"infile","format":"denovo_rna_v2.common"}
        ]
        self.add_option(options)
        self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_option(self):
        """
        参数检查
        """
        if not self.option("fasta_dir"):
            raise OptionError("没有找到fasta_dir")
        if not self.option("list_xls"):
            raise OptionError("没有找到list.xls")
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 8
        self._memory = "125G"

    def end(self):
        super(BacIdentityAgent, self).end()


class BacIdentityTool(Tool):
    def __init__(self, config):
        super(BacIdentityTool, self).__init__(config)
        self.mj_info = namedtuple('mj_info', [
                                  'mj', 'primer', 'rawLength', 'trimLength', 'identity', 'quality', 'specieses','asses'])
        self.info = {}#各个样本的信息，为可命名元祖的list
        self.Rscript = "bioinfo/tool_lab/miniconda2/envs/R4version/bin/Rscript"
        self.software_dir = self.config.SOFTWARE_DIR
        # self.env = 'export LD_LIBRARY_PATH={}:$LD_LIBRARY_PATH'.format(self.config.SOFTWARE_DIR + '/library/glibc-2.14/lib/')
        # self.perl_script = self.config.PACKAGE_DIR + "/tool_lab/bac_identity/gi_tax.pl"
        self.set_environ(
            LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/library/glibc-2.14/lib/')
        self.set_environ(
            BLASTDB=self.config.SOFTWARE_DIR + '/database/Tool_lab/bac_identity/NT_20211102/')
        self.script = self.config.PACKAGE_DIR + '/tool_lab/bac_identity/sangerSeq2_new.R'
        self.assemble_sh = 'bioinfo/tool_lab/assemble.sh.new'
        self.blastn = "bioinfo/align/ncbi-blast-2.10.1+/bin/blastn"
        # self.perl = 'program/perl-5.24.0/bin/perl'
        # 20211102 更新数据库，改用taxonkit搜索物种名
        self.python_script = self.config.PACKAGE_DIR + "/tool_lab/bac_identity/getNtTaxon.py"
        self.python = "miniconda2/bin/python"
        self.nt = os.path.join(self.software_dir, "database/Tool_lab/bac_identity/NT_20211102/nt")
        # self.nt_tax = os.path.join(
        #     self.software_dir, "database/Tool_lab/bac_identity/NT_tax.xls")
        self.readme = os.path.join(self.software_dir,"database/Tool_lab/bac_identity/README.txt")
        self.pictrue_dir = "/mnt/ilustre/users/sanger-dev/junjian/picture"
        
        # self.mj_primer_quatity= {}

    def run(self):
        """
        运行
        """
        super(BacIdentityTool, self).run()
        self.run_sangerSeq2_R()
        self.run_collate_Rresult()
        self.run_assemble()
        self.run_nt_blastn()
        self.check_M8()
        self.run_getreuslt()
        self.set_output()
        self.end()

    def run_sangerSeq2_R(self):
        """
        Rscript sangerSeq2_new.R list.xls
        """
        self.fasta_dirname = self.option("fasta_dir").prop["path"].split('/')[-1]
        self.fasta_dir = os.path.join(self.work_dir, self.fasta_dirname)
        try:
            shutil.copytree(self.option("fasta_dir").prop["path"],self.fasta_dir)
        except Exception as e:
            self.set_error(e)
        os.chdir(self.fasta_dir)
        cmd = '{} {}'.format(self.Rscript, self.script)
        cmd += " {}".format(self.option("list_xls").prop["path"])
        self.logger.info(cmd)
        self.logger.info("开始鉴定污染")
        command1 = self.add_command("sangerseq_r", cmd,ignore_error=True)
        command1.run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("鉴定完成")
        else: #判定报错信息，获取报错信息
            if os.path.exists(self.work_dir+'/sangerseq_r.o'):
                with open(self.work_dir+'/sangerseq_r.o','r') as error_info:
                    while 1:
                        line = error_info.readline()
                        if not line:
                            break
                        if 'Error in rawToChar(x, ...) : embedded nul in string' in line:
                            self.set_error("ab1文件传输方式请设置成二进制")
            self.set_error("鉴定污染出现未知错误")
            # self.set_error("鉴定出现错误")

    def run_collate_Rresult(self):
        outReport1 = open(self.work_dir+"/sanger1.report", 'w')
        outReport1.write(
            'sample\tprimer\trawSeqLength\ttrimSeqLength\tidentity\tquality\n')
        san = open(self.work_dir+'/sanger.report', 'w')
        san.write('\t'.join(["#Majorbio-No", "Primer", "rawSeqLength",
                             "trimSeqLength", "QC-Identity", "Quality-Rank", "Date"])+'\n')
        inResult = os.path.join(self.fasta_dir, "result.txt")
        sn_primer = {}
        with open(self.option("list_xls").prop['path'],'r') as lx:
            while 1:
                line = lx.readline()
                if not line:
                    break
                fd = line.rstrip().split('\t')
                if sn_primer.has_key(fd[1]):
                    sn_primer[fd[1]].append(fd[2])
                else:
                    sn_primer[fd[1]] = [fd[2]]
        with open(inResult, "r") as rt:
            rt.readline()
            while 1:
                line = rt.readline()
                if not line:
                    break
                term = line.split()
                info = term[1].split("-")
                mj = info[0]
                # self.mj_primer_quatity[mj] = {}
                if not self.info.has_key(mj):
                    self.info[mj] = self.mj_info(mj, [],[],[],[],[],[],[])
                tem = "-".join(info[1:])
                primer = tem.split('.')[0].split("_")[0]
                self.info[mj].primer.append(primer)
                identity = term[6].lstrip("(").rstrip("%)")
                self.info[mj].identity.append(identity)
                # ratio = info[5].split("/")
                trimLength = term[10]
                self.info[mj].trimLength.append(trimLength)
                rawLength = term[12].rstrip()
                self.info[mj].rawLength.append(rawLength)
                threshold = self.option("threshold")
                quality = ""
                if float(identity) < 99:
                    quality = "污染"
                elif float(identity) < threshold:
                    quality = "疑似污染"
                else:
                    quality = "正常"
                self.info[mj].quality.append(quality)
                # self.mj_primer_quatity[mj][primer] = quality

                today = time.strftime("%Y%m%d", time.localtime())
                outReport1.write(
                    '\t'.join([mj, primer, rawLength, trimLength, identity, quality])+'\n')
        for sn in sorted(self.info.keys()): 
            if sn_primer[sn][0] == self.info[sn].primer[0] and sn_primer[sn][1] == self.info[sn].primer[1]:
                san.write(
                    '\t'.join([self.info[sn].mj, ",".join(self.info[sn].primer), ",".join(self.info[sn].rawLength), ",".join(self.info[sn].trimLength), ",".join(self.info[sn].identity), ",".join(self.info[sn].quality), today]))
                san.write("\n")
            elif sn_primer[sn][1] == self.info[sn].primer[0] and sn_primer[sn][1] == self.info[sn].primer[0]:
                self.info[sn].primer.reverse()
                self.info[sn].rawLength.reverse()
                self.info[sn].trimLength.reverse()
                self.info[sn].identity.reverse()
                self.info[sn].quality.reverse()
                san.write(
                    '\t'.join([self.info[sn].mj, ",".join(self.info[sn].primer), ",".join(self.info[sn].rawLength), ",".join(self.info[sn].trimLength), ",".join(self.info[sn].identity), ",".join(self.info[sn].quality), today]))
                san.write("\n")
            else:
                self.set_error("list信息和鉴定结果样品信息不一致")
        outReport1.close()
        san.close()

    def run_assemble(self):
        # fasta_dirname = os.path.basename(self.option("fasta_dir").prop["path"])
        # fasta_dir = os.path.join(self.work_dir, fasta_dirname)
        os.system('cd {}'.format(self.fasta_dir))
        command1 = self.add_command("assemble", self.assemble_sh).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("assemble完成")
        else:
            self.set_error("assemble出现错误")
        allcontigs = open(os.path.join(self.work_dir, "all.contigs"), "w")
        sample_num = 0
        sample_list = []
        with open(self.option("sample_info").prop['path'],'r') as lx:
            sample_num = len(lx.readlines())
        with open(self.option("sample_info").prop['path'],'r') as lx:
            while 1:
                line = lx.readline()
                if not line:
                    break
                if line[0] == "#":
                    continue
                self.logger.info(line.rstrip().split('\t'))
                sample_list.append(line.rstrip().split('\t')[7])
        # nofound_file = []
        file_list = []
        self.logger.info(sample_list)
        self.logger.info(sample_num)
        for root, dir_name, file in os.walk(self.fasta_dir + "/contigs"):
            if len(file) < sample_num:
                self.logger.info("assemble出现错误，组装的样本小于列表的样本")
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
            if len(nofound_file) > 0:
                self.logger.info("样本{}没有生成".format(",".join(nofound_file)))
            for name in file:
                if name == "readme.txt":
                    os.remove(os.path.join(root,name))
                    continue
                with open(os.path.join(root, name), "r") as contig, open(os.path.join(root, name.split(".")[0]+".seq"), "w") as seq:
                    n = 0
                    text =""
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
                        primer = self.info[mj].primer[0] if float(self.info[mj].identity[0]) > float(self.info[mj].identity[1]) else self.info[mj].primer[1]
                        fafile = glob.glob(os.path.join(self.fasta_dir+'/'+mj+'-'+primer+'*.fa'))[0]
                        fa_record = SeqIO.read(fafile,"fasta")
                        SeqIO.write(fa_record, self.fasta_dir+"/temp.txt","fasta")
                        with open(self.fasta_dir+"/temp.txt","r") as temp_seq:
                            text = temp_seq.read()
                    seq.write(text)
                allcontigs.write(text)
                os.remove(os.path.join(root, name))
            for nf_sn in nofound_file:
                    # for pm in self.info[nf_sn].primer:
                    #     if os.path.exists('{}-{}.fa')
                if len(self.info[nf_sn].primer) == 2:
                    pm1 = self.info[nf_sn].primer[0]
                    pm2 = self.info[nf_sn].primer[1]
                    file1 = self.fasta_dir + '/{}-{}.fa'.format(nf_sn,pm1)
                    file2 = self.fasta_dir + '/{}-{}.fa'.format(nf_sn,pm2)
                    if os.path.exists(file1) and os.path.exists(file2):
                        contig_file = file1 if os.path.getsize(file1) >= os.path.getsize(file2) else file2
                        fa_record = SeqIO.read(contig_file,"fasta")
                        SeqIO.write(fa_record, self.fasta_dir+"/temp.txt","fasta")
                        with open(self.fasta_dir+"/temp.txt","r") as temp_seq:
                            text = temp_seq.read()
                        allcontigs.write(text)
        allcontigs.close()

    def run_nt_blastn(self):
        """
        blastn -db nt -query contigs/all.contigs -out NT_blast.m8.txt -evalue 1e-5 -outfmt 6 -max_target_seqs 10
        perl /gi_tax.pl NT_blast.m8.txt /NT_tax.xls  NT.fasta.blast.taxonomy.xls    
        """
        cmd = "{} -db {} -query {}".format(self.blastn, self.nt, os.path.join(self.work_dir, "all.contigs"))
        out = os.path.join(self.work_dir, "NT_blast.m8.txt")
        out_tax = os.path.join(self.work_dir, "NT.fasta.blast.taxonomy.temp.xls")
        cmd += " -out {} -evalue 1e-5 -max_target_seqs 10 -num_threads 2 -outfmt \"6 -delim \\\"\\t\\\" qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids \"".format(out)
        self.logger.info(cmd)
        command = self.add_command("blastn", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行blastn完成")
        else:
            self.logger.info("运行错误，请重新检查输入参数")
        # cmd1 = "{} {}".format(self.perl, self.perl_script)
        # cmd1 += " {} {} {}".format(out, self.nt_tax, out_tax)
        #更换了获取物种名的方法
        cmd1 = "{} {}".format(self.python, self.python_script)
        cmd1 += " -i {} -o {}".format(out,out_tax)
        command1 = self.add_command("gettax", cmd1).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("运行getNtTax完成")
        else:
            self.set_error("运行错误，请重新检查输入参数")
        
        
    def check_M8(self):
        NT_fbt = os.path.join(self.work_dir, "NT.fasta.blast.taxonomy.temp.xls")
        new_NT = os.path.join(self.work_dir, "NT.fasta.blast.taxonomy.xls")
        with open(NT_fbt,"r") as old,open(new_NT,"w") as new:
            n = 0
            sn_last = ""
            while 1:
                line = old.readline()
                if not line:
                    break
                sn = line.rstrip().split("\t")[0]
                if sn != sn_last:
                    sn_last = sn
                    n = 0
                if n < 10:
                    new.write(line)
                n += 1

    def run_getreuslt(self):
        m8file = open(os.path.join(self.work_dir, "blastM8.xls"), "w")
        taxfile = open(os.path.join(self.work_dir, "NT.tax.xls"), "w")
        m8file.write(
            '#Query\tTaxonomy\tSubject\tIdentity%\tAlign_length\tMismatches\tGap_open\tQ_start\tQ_end\tS_start\tS_end\tEvalue\tBit_score\n')
        taxfile.write('#Sample\tAssession ID\tTaxonomy\tBase calling\n')
        sn_species = {}
        sn_ass = {}
        sn_score={}
        with open(os.path.join(self.work_dir, "NT.fasta.blast.taxonomy.xls"), "r") as m8:
            while 1:
                line = m8.readline()
                if not line:
                    break
                m8file.write(line)
                temp = line.rstrip().split("\t")
                if len(temp) == 13:
                    species = temp[1]
                    ass = temp[2]
                    score = float(temp[-1])
                    if sn_score.has_key(temp[0]):
                        if sn_score[temp[0]].has_key(species):
                            if sn_score[temp[0]][species] < score:
                                sn_score[temp[0]][species] = score
                        else:
                            sn_score[temp[0]][species] = score
                    else:
                        sn_score[temp[0]] = {}
                        sn_score[temp[0]][species] = score
                    self.info[temp[0]].specieses.append(species)
                    self.info[temp[0]].asses.append(ass)
                    # self.info[temp[0]].specieses.append(species)
                    # self.info[temp[0]].asses.append(ass)
        for sn in sorted(self.info.keys()):
            if len(self.info[sn].specieses) == 0:
                sn_species[sn] = "NO-TaxID"
                sn_ass[sn] = "NO-TaxID"
                # self.info[sn].species = "NO-TaxID"
                # self.info[sn].ass = "NO-TaxID"
                continue
            n = 0
            genus = ""
            score = 0
            sn_times = {}
            tax_ass={}
            while 1:
                if n >= len(self.info[sn].specieses):
                    if len(sn_times.keys()) >= 2:
                        fn=""
                        fass=''
                        fs=0
                        for i in sn_times.keys():
                            if sn_times[i] > fs:
                                fs = sn_times[i]
                                fn = i.split(";")[-1].split("s__")[-1]
                                fass = tax_ass[i]
                        sn_species[sn] = fn
                        sn_ass[sn] = fass
                        break
                    elif len(sn_times.keys()) == 1:
                        sn_species[sn] = sn_times.keys()[0].split(";")[-1].split("s__")[-1]
                        sn_ass[sn]=tax_ass[sn_times.keys()[0]]
                        break
                elif sn_score[sn][self.info[sn].specieses[n]] < score:
                    if len(sn_times.keys()) >= 2:
                        fn=""
                        fass=''
                        fs=0
                        for i in sn_times.keys():
                            self.logger.info("{}:{}=={}".format(sn,i,sn_times[i]))
                            if sn_times[i] > fs:
                                fs = sn_times[i]
                                fn = i.split(";")[-1].split("s__")[-1]
                                fass = tax_ass[i]
                        sn_species[sn] = fn
                        sn_ass[sn] = fass
                        break
                    elif len(sn_times.keys()) == 1:
                        sn_species[sn] = sn_times.keys()[0].split(";")[-1].split("s__")[-1]
                        sn_ass[sn]=tax_ass[sn_times.keys()[0]]
                        break
                
                if re.search('s__uncultured|s__.*_bacterium|s__bacterium|NO-TaxID',self.info[sn].specieses[n]) and re.search("norank",";".join(self.info[sn].specieses[n].split(";")[-6:])):
                    if genus == "":
                        if re.search("norank",self.info[sn].specieses[0]):
                            if not re.search("p__norank",";".join(self.info[sn].specieses[0].split(";")[-6:])):
                                tax_info = self.info[sn].specieses[0].split(";")[-6:]
                                tem_info = ";".join(self.info[sn].specieses[0].split(";")[-6:])
                                tem_level = tax_info.index(tem_info.split("norank")[0].split(';')[-2])
                                genus = ";".join(tax_info[tem_level:])
                                n += 1
                                continue
                            else:
                                genus = ";".join(self.info[sn].specieses[0].split(";")[-8:])
                                n += 1
                                continue
                        else:
                            genus = ";".join(self.info[sn].specieses[0].split(";")[-2:])
                            n += 1
                            continue
                    if re.search("norank",genus) and not re.search("norank",";".join(self.info[sn].specieses[n].split(";")[-6:])):
                        genus = ";".join(self.info[sn].specieses[n].split(";")[-2:])
                    elif re.search("norank",genus) and not re.search("p__norank",self.info[sn].specieses[n]):
                        tax_info = self.info[sn].specieses[n].split(";")[-6:]
                        tem_info = ";".join(self.info[sn].specieses[n].split(";")[-6:])
                        tem_level = tax_info.index(tem_info.split("norank")[0].split(';')[-2])
                        if tem_level > (7-len(genus.split(';'))):
                            genus = ";".join(tax_info[tem_level:])
                    n += 1
                    if n == len(self.info[sn].specieses):
                        # self.info[sn].species = genus
                        sn_species[sn] = genus
                        # self.info[sn].ass = self.info[sn].asses[n-1]
                        sn_ass[sn] = self.info[sn].asses[n-1]
                        break
                    continue
                else:
                    # genus = self.info[sn].specieses[n].species[-1].split("g__")[1].split(";")[0]
                    # species = self.info[sn].specieses[n].split(";")[-1].split("s__")[-1]
                    # self.info[sn].species = species
                    # sn_species[sn] = species
                    # self.info[sn].ass = self.info[sn].asses[n]
                    # sn_ass[sn] = self.info[sn].asses[n]
                    # tax_ass[self.info[sn].specieses[n]] = self.info[sn].asses[n]
                    '''
                    增加score的判定因素
                    '''
                    if sn_score[sn][self.info[sn].specieses[n]] >= score:
                        score = sn_score[sn][self.info[sn].specieses[n]]
                        if sn_times.has_key(self.info[sn].specieses[n]):
                            sn_times[self.info[sn].specieses[n]] += 1
                        else:
                            tax_ass[self.info[sn].specieses[n]] = self.info[sn].asses[n]
                            sn_times[self.info[sn].specieses[n]] = 1
                    n+=1
                    continue
                # if re.search('s__uncultured|s__.*_sp\.|s__.*_bacterium|s__bacterium|NO-TaxID', i):
                #             continue
        m8file.close()
        for sn in sorted(self.info.keys()):
            # ass = self.info[sn].ass
            ass =sn_ass[sn]
            # tax = self.info[sn].species
            tax = sn_species[sn]
            base_calling = []
            for pr in self.info[sn].primer:
                index= self.info[sn].primer.index(pr)
                base_calling.append("{}:{}".format(pr,self.info[sn].quality[index]))
            if "污染" in self.info[sn].quality:
                taxfile.write("{}\t-\t-\t{}\n".format(sn,",".join(base_calling)))
            else:
                taxfile.write("{}\t{}\t{}\t{}\n".format(sn,ass,tax,",".join(base_calling)))
        taxfile.close()

    def set_output(self):
        '''
        将结果文件赋值到output文件夹下面
        :return:
        '''
        # fasta_dir = os.path.basename(self.option("fasta_dir"))
        sample_info = self.option("sample_info").prop["path"]
        # fasta_dirname = os.path.basename(self.option("fasta_dir").prop["path"])
        # fasta_dir = os.path.join(self.work_dir, fasta_dirname)
        dir_name = ""
        with open(sample_info,"r") as sinfo:
            sinfo.readline()
            line = sinfo.readline()
            temp= line.split("\t")
            temp= line.split("\t")
            # client = temp[5].replace(' ','_')
            # dir_name = "{}-{}-{}-{}-{}sample".format(temp[0],temp[1],client,temp[4],temp[6])
            # result_name = "{}-{}-{}-{}-{}sample".format(temp[0],temp[1],client,temp[4],temp[6])
            # pictrue_path = os.path.join(self.pictrue_dir,"{}.jpg".format(temp[4]))
            dir_name = "result"
            result_name = "result"
            if len(glob.glob(os.path.join(self.pictrue_dir+'/'+temp[4]+'='+'*.jpg'))) > 0:
                pictrue_path = glob.glob(os.path.join(self.pictrue_dir+'/'+temp[4]+'='+'*.jpg'))[0]
                new_pictrue_name = "{}.jpg".format(temp[4])
            else:
                time.sleep(60)
                self.logger.info("1 min")
                if len(glob.glob(os.path.join(self.pictrue_dir+'/'+temp[4]+'='+'*.jpg'))) > 0:
                    pictrue_path = glob.glob(os.path.join(self.pictrue_dir+'/'+temp[4]+'='+'*.jpg'))[0]
                    new_pictrue_name = "{}.jpg".format(temp[4])
                else:
                    time.sleep(180)
                    self.logger.info("3 min")
                    if len(glob.glob(os.path.join(self.pictrue_dir+'/'+temp[4]+'='+'*.jpg'))) > 0:
                        pictrue_path = glob.glob(os.path.join(self.pictrue_dir+'/'+temp[4]+'='+'*.jpg'))[0]
                        new_pictrue_name = "{}.jpg".format(temp[4])
                    else:
                        time.sleep(300)
                        self.logger.info("5 min")
                        if len(glob.glob(os.path.join(self.pictrue_dir+'/'+temp[4]+'='+'*.jpg'))) > 0:
                            pictrue_path = glob.glob(os.path.join(self.pictrue_dir+'/'+temp[4]+'='+'*.jpg'))[0]
                            new_pictrue_name = "{}.jpg".format(temp[4])
                        else:
                            time.sleep(600)
                            self.logger.info("10 min")
                            if len(glob.glob(os.path.join(self.pictrue_dir+'/'+temp[4]+'='+'*.jpg'))) > 0:
                                pictrue_path = glob.glob(os.path.join(self.pictrue_dir+'/'+temp[4]+'='+'*.jpg'))[0]
                                new_pictrue_name = "{}.jpg".format(temp[4])
                            else:
                                self.set_error("找不到订单号{}的胶图,请检查胶图是否上传，如已上传，请10分钟后重新运行".format(temp[4]))
        if len(self.output_dir) > 0:
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        ot_dir = os.path.join(self.output_dir, result_name)
        os.mkdir(ot_dir)
        os.mkdir(ot_dir + "/" + dir_name)
        rstdir = ot_dir + "/" + dir_name
        # os.mkdir(self.output_dir + "/" + dir_name + "/" + "3730")
        dir3730 = ot_dir + "/" + dir_name + "/" + "3730"
        os.mkdir(ot_dir + "/" + dir_name + "/" + "NCBI")
        dirNCBI = ot_dir + "/" + dir_name + "/" + "NCBI"
        # os.mkdir(self.output_dir + "/" + dir_name + "/" + "contigs")
        dircontigs = ot_dir + "/" + dir_name + "/" + "contigs"
        self.logger.info("设置结果目录")
        try:
            os.link(os.path.join(self.work_dir, "sanger.report"),
                    os.path.join(ot_dir, "sanger.report"))
            shutil.copyfile(pictrue_path,
                    os.path.join(ot_dir, new_pictrue_name))
            # os.link(os.path.join(self.work_dir, "sanger1.report"),
            #         os.path.join(self.output_dir, "sanger1.report"))
            shutil.copytree(os.path.join(self.fasta_dir, "contigs"),
                            dircontigs)
            os.link(os.path.join(self.work_dir, "blastM8.xls"),
                    os.path.join(dirNCBI, "blastM8.xls"))
            os.link(os.path.join(self.work_dir,"NT.tax.xls"),
                    os.path.join(dirNCBI, "NT.tax.xls"))
            shutil.copytree(self.option("fasta_dir").prop["path"],
                            dir3730)
            os.link(self.readme, os.path.join(rstdir,"README.txt"))
            os.link(sample_info, os.path.join(rstdir,"sample.detail.xls"))
            
        except Exception as e:
            self.set_error("设置结果目录失败{}".format(e))

class TestFunction(unittest.TestCase):
    """
    测试脚本用
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'bac_identity_' + str(random.randint(1, 100000)),
            "type": "tool",
            "name": "tool_lab.bac_identity",
            "options": {
                 "fasta_dir":"/mnt/ilustre/users/sanger-dev/tsanger/workspace/20210730/BacIdentity_175774_20210730_091444/BacDatapre/fasta_dir",
                 "sample_info":"/mnt/ilustre/users/sanger-dev/tsanger/workspace/20210730/BacIdentity_175774_20210730_091444/BacDatapre/sample.detail.xls",
                 "list_xls":"/mnt/ilustre/users/sanger-dev/tsanger/workspace/20210730/BacIdentity_175774_20210730_091444/BacDatapre/list.xls"
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
