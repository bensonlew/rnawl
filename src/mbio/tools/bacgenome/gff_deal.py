# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
from biocluster.core.exceptions import OptionError
import subprocess


class GffDealAgent(Agent):
    """
    last_modify: 2018.04.23
    """

    def __init__(self, parent):
        super(GffDealAgent, self).__init__(parent)
        options = [
            {"name": "fna", "type": "infile", "format": "sequence.fasta"},
            {"name": "all_gff", "type": "infile", "format": "gene_structure.gff3"},
            {"name": "analysis", "type": "string", "default": "complete"},  ###流程分析模式complete，uncomplete
            {'name': 'sample_name', "type": "string"},  # 样本名
            {"name": "rrna_gff", "type": "outfile", "format": "gene_structure.gff3" },
            {"name": "trna_gff", "type": "outfile", "format": "gene_structure.gff3"},
            {"name": "gene_gff", "type": "outfile", "format": "gene_structure.gff3"},
            {"name": "faa", "type":"outfile", "format": "sequence.fasta"},
            {"name": "ffn", "type":"outfile", "format": "sequence.fasta"},
            {"name": "txt_info","type":"string","default":""}  #{"chr_num","plasmid_num"}
            #{"name": "gene_statistics", "type": "outfile", "format": "sequence.profile_table"},
            #{"name": "fna_summary", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("fna").is_set:
            raise OptionError("请传入基因组文件！")
        if not self.option("all_gff").is_set:
            raise OptionError("请传入gff文件！")
        if not self.option("sample_name"):
            raise OptionError("请传入样本名称！")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 1
        self._memory = '2G'

    def end(self):
        super(GffDealAgent, self).end()

class GffDealTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(GffDealTool, self).__init__(config)
        self.gffread_path = "bioinfo/rna/cufflinks-2.2.1/"
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/sequence/scripts/"
        self.transeq = "/bioinfo/seq/EMBOSS-6.6.0/emboss/transeq"
        self.LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/bioinfo/seq/EMBOSS-6.6.0/lib"
        self.set_environ(LD_LIBRARY_PATH=self.LD_LIBRARY_PATH)
        self.seqkit = "/bioinfo/seq/seqkit"

    def get_fna_len(self):
        pat = re.compile('\s+')
        fna = self.option("fna").prop['path']
        self.each_scf = {}
        self.each_scf_len = {}
        self.scf_order = []
        with open(fna) as fr:
            for line in fr:
                if line[0] == '>':
                    k = pat.split(line)[0][1:]
                    if k not in self.each_scf_len.keys():
                        self.scf_order.append(k)
                        self.each_scf[k] = ''
                    else:
                        raise(k + ' duplication name')
                else:
                    self.each_scf[k] += line.strip()

        fw = open('sequence_len.xls','w')
        fw.write('Sequence id\tlength\tseq_id\n')
        for k in self.each_scf:
            seq_len = len(self.each_scf[k])
            self.each_scf_len[k] = seq_len
            fw.write(k+'\t'+str(seq_len)+'\t'+k+'\n')
        fw.close()


    def get_gene_trna_rrna_gff(self):
        self.gene_info = {}

        gff = self.option("all_gff").prop["path"]
        with open(gff) as fr,open('gene.gff','w') as gene, open('trna.gff','w') as trna, open('rrna.gff','w') as rrna, open('gene_desc.txt','w') as g_desc:
            fw_map = {'CDS':gene, 'tRNA':trna, 'rRNA':rrna}
            gene.write('Gene ID\tSequence id\tStart\tEnd\tStrand\tGene Length(bp)\tProtein Length\tA.start\tA.end\tInitiator Codon\tTerminator Codon\n')
            trna.write('Gene ID\tSequence id\tStart\tEnd\ttRNA Type\tAnti Codon\tIntron Begin\tBounds End\tScore\tA.start\tA.end\tInitiator Codon\tTerminator Codon\n')
            rrna.write('Gene ID\tSequence id\tStart\tEnd\tE-values\tStrand\tPhase\tAttributes\tA.start\tA.end\tInitiator Codon\tTerminator Codon\n')
            g_desc.write('Gene ID\tgene_ori_name\tgene_ori_desc\n')
            for line in fr:
                if line.startswith('#'):
                    continue
                line = line.strip()
                sp_line = line.split('\t')
                type = sp_line[2]
                if type not in ['CDS','tRNA','rRNA']:
                    continue
                seq_index = self.scf_order.index(sp_line[0])
                add_len = 0
                for pre in self.scf_order[:seq_index]:
                    add_len += self.each_scf_len[pre]

                start = sp_line[3]
                end =  sp_line[4]
                if int(start) > self.each_scf_len[sp_line[0]] or int(end) > self.each_scf_len[sp_line[0]]:  #去除有问题的基因
                    continue
                a_start = int(start) + add_len
                a_end = int(end) + add_len

                if type == 'CDS': #,'tRNA','rRNA']:
                    f_name = re.findall(';Parent=([^;]*)',sp_line[8])
                    if len(f_name) == 0:
                        self.logger.info('cannot match Parent : %s'%line)
                        continue
                    else:
                        name = f_name[0]
                    seq_id = sp_line[0] + "_ORF"
                    strand = sp_line[6]
                    nuc_len = abs(int(end)-int(start)+1)
                    pro_len = nuc_len/3
                    if self.option('analysis') == 'complete':
                        fw_map[type].write('\t'.join([name,seq_id,start,end,strand,str(nuc_len),str(pro_len),start, end,'-','-\n']))
                    else:
                        fw_map[type].write('\t'.join([name,seq_id,start,end,strand,str(nuc_len),str(pro_len),str(a_start), str(a_end),'-','-\n']))
                    self.gene_info[name] = {'start':start, 'end':end, 'strand':strand,'seq_id':seq_id}
                    gene = re.findall('gene=([^;]*)',sp_line[8])
                    desc = re.findall('product=([^;]*)',sp_line[8])
                    if gene:
                        gene = gene[0]
                    else:
                        gene = ''
                    if desc:
                        desc = desc[0]
                    else:
                        desc = ''
                    g_desc.write('\t'.join([name,gene,desc])+'\n')
                elif type == 'tRNA':
                    f_name =  re.findall(';Parent=([^;]*)',sp_line[8])
                    if len(f_name) == 0:
                        self.logger.info('cannot match Parent : %s'%line)
                        continue
                    else:
                        name = f_name[0]
                    seq_id = sp_line[0] + '_tRNA'
                    find_result = re.findall('Note=tRNA-([^\(]*)\s*\(([\w]*)\)',sp_line[8])
                    if len(find_result) != 0:
                        t_type = find_result[0][0]
                        anti_code = find_result[0][1]
                    else:
                        find_result = re.findall('product=([^;]*)',sp_line[8])
                        if len(find_result) !=0:
                            t_type = find_result[0]
                            anti_code = '***'
                        else:
                            self.logger.info('cannot match producter and anti code infomation %s'%line)
                            continue
                    if self.option('analysis') == 'complete':
                        trna.write('\t'.join([name,seq_id,start,end,t_type,anti_code,'-','-','-',start,end,'-','-\n']))
                    else:
                        trna.write('\t'.join([name,seq_id,start,end,t_type,anti_code,'-','-','-',str(a_start),str(a_end),'-','-\n']))
                elif type == 'rRNA':
                    f_name =  re.findall(';Parent=([^;]*)',sp_line[8])
                    if len(f_name) == 0:
                        self.logger.info('cannot match Parent : %s'%line)
                        continue
                    else:
                        name = f_name[0]
                    seq_id = sp_line[0] + '_rRNA'

                    attribute = re.findall('product=([^;]*)', sp_line[8])[0]
                    strand = sp_line[6]
                    if '5S' in attribute:
                        attribute2 = 'Name=5S_rRNA;product=5S ribosomal RNA'
                    elif '16S' in attribute:
                        attribute2 = 'Name=16S_rRNA;product=16S ribosomal RNA'
                    elif '23S' in attribute:
                        attribute2 = 'Name=23S_rRNA;product=23S ribosomal RNA'
                    else:
                        attribute2 = attribute
                    if self.option('analysis') == 'complete':
                        rrna.write('\t'.join([name,seq_id,start,end,'-',strand,'-',attribute2,start, end,'-' , '-\n']))
                    else:
                        rrna.write('\t'.join([name,seq_id,start,end,'-',strand,'-',attribute2,str(a_start),str(a_end),'-' , '-\n']))



    def extract_gene_ffn(self):
        gff = self.option("all_gff").prop["path"]
        fasta = self.option('fna').prop["path"]

        self.logger.info("开始处理fasta，使每行碱基数一致")
        new_all_fasta = self.work_dir +'/new_all.fasta'
        cmd = "{}  seq -w 80 -o {}  {}".format(self.seqkit,new_all_fasta, fasta )
        command = self.add_command("samelen", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("seqkit 成功")
        else:
            raise Exception("seqkit 失败")

        cmd = "{}gffread {} -g {} -x cds.ffn_ori".format(self.gffread_path, gff, new_all_fasta)  #-y cds.faa_ori
        self.logger.info("开始运行cufflinks的gffread，合成、提取CDS")
        command = self.add_command("gffread", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("CDS提取完成")
            with open('cds.ffn','w') as fw, open('cds.ffn_ori') as fr:
                skip = 0
                for line in fr:
                    if line[0] == '>':
                        skip = 0
                        name = line.split()[0][1:]
                        if name in self.gene_info.keys():
                            info = self.gene_info[name]
                        else:
                            skip = 1
                            continue
                        loc = info['seq_id']+'_'+name
                        if info['strand'] == '-':
                            line = '>' + name +' ' + ' '.join([info['end'] , info['start'] ,loc,info['end'],info['start']])
                        else:
                            line = '>' + name +' ' + ' '.join([info['start'] , info['end'] ,loc ,info['start'],info['end']])
                        fw.write(line+'\n')
                    else:
                        if skip != 1:
                            fw.write(line)

            nul_seq = 'cds.ffn'
            prot_seq = 'cds.faa_ori'
            cmd = '{} -trim -table 11 -sequence {} -outseq {}'.format(self.transeq, nul_seq, prot_seq)
            command = self.add_command("transeq", cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("翻译蛋白序列运行完成")

            with open('cds.faa','w') as fw, open('cds.faa_ori') as fr:
                for line in fr:
                    if line[0] == '>':
                        name = line.split()[0][1:-2] ##末尾有_1
                        info = self.gene_info[name]
                        loc = info['seq_id']+'_'+name
                        if info['strand'] == '-':
                            line = '>' + name +' ' + ' '.join([info['end'] , info['start'] ,loc ,info['end'],info['start']])
                        else:
                            line = '>' + name +' ' + ' '.join([info['start'] , info['end'] ,loc ,info['start'],info['end']])
                        fw.write(line+'\n')
                    else:
                        line = line.replace('.','*')
                        fw.write(line)

        else:
            self.set_error("运CDS提取出错")

    # def gene_stat(self):
    #     all_gene_fnn = self.work_dir + '/cds.ffn'
    #     all_fasta = self.option('fna').prop["path"]
    #     out_path = self.output_dir + "/" + "CDS_predict"
    #     cmd = '{} {}dnabac_sample_stat.pl {} {} {} {}'.format(self.perl_path, self.perl_script, all_fasta, all_gene_fnn,
    #                                                           self.option("sample_name"), out_path)
    #
    #     command = self.add_command("sample_gene_info_stat", cmd).run()
    #     self.wait(command)
    #     if command.return_code == 0:
    #         self.logger.info("样品编码基因预测统计完成")
    #     else:
    #         self.set_error("样品编码基因预测统计出错!")

    def run_summary(self):
        with open(self.work_dir+'/trna.gff') as trna,  open(self.work_dir+'/rrna.gff') as rrna, open(self.work_dir+'/gene.gff') as gene:
            trna_num = len(trna.readlines())-1
            rrna_num = len(rrna.readlines())-1
            gene_num = len(gene.readlines())-1
        scaf_num = 0
        seq = ''
        plasmid_num = 0
        with open(self.option('fna').prop["path"]) as fna:
            for line in fna:
                if line[0] ==  '>':
                    if 'plasmid' in line:
                        plasmid_num +=1
                    scaf_num +=1
                else:
                    line = line.strip()
                    seq += line
        chr_num = scaf_num - plasmid_num

        if self.option('txt_info'):
            txt_info = eval(self.option('txt_info'))
            plasmid_num = txt_info['plasmid_num']
            chr_num = txt_info['chr_num']

        base_num = len(seq)
        gc_num = seq.count('g') + seq.count('G') + seq.count('C')+seq.count('c')
        gc_percent = round(gc_num*1.0/base_num,4)*100
        if self.option('analysis') in ['uncomplete']:
            with open('project.summary', 'w') as fw:
                fw.write('Sample\tGenome Size\tscaffold no\tGC Content(%)\tCDS No.\trRNA No.\ttRNA No.\n')
                fw.write(self.option('sample_name')+'\t'+'\t'.join([str(i) for i in [base_num, scaf_num, gc_percent,gene_num,rrna_num,trna_num]]) + '\n')
        else:
            with open('project.summary', 'w') as fw:
                fw.write('Sample\tGenome Size\tChrom No.\tPlas No.\tGC Content(%)\tCDS No.\trRNA No.\ttRNA No.\n')
                fw.write(self.option('sample_name')+'\t'+'\t'.join([str(i) for i in [base_num,chr_num,plasmid_num, gc_percent,gene_num,rrna_num,trna_num]]) + '\n')


    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        self.logger.info("设置结果目录")
        file_names = ['gene.gff','trna.gff','rrna.gff','cds.faa','cds.ffn', 'project.summary','sequence_len.xls']
        for f in file_names:
            if os.path.exists(self.output_dir + '/' + self.option('sample_name') + '.' + f ):
                os.remove(self.output_dir + '/' + self.option('sample_name') + '.' + f)
            os.link(self.work_dir + '/' + f,self.output_dir + '/' + self.option('sample_name') + '.' + f)

        self.option('ffn',self.output_dir + '/' + self.option('sample_name') + '.cds.ffn')
        self.option('faa',self.output_dir + '/' + self.option('sample_name') + '.cds.faa')
        self.option('trna_gff',self.output_dir + '/' + self.option('sample_name') + '.trna.gff')
        self.option('rrna_gff',self.output_dir + '/' + self.option('sample_name') + '.rrna.gff')
        self.option('gene_gff',self.output_dir + '/' + self.option('sample_name') + '.gene.gff')
        self.logger.info("设置结果目录成功")


    def run(self):
        super(GffDealTool, self).run()
        self.get_fna_len()
        self.get_gene_trna_rrna_gff()
        self.extract_gene_ffn()
        self.run_summary()
        self.set_output()
        self.end()
