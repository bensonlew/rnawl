# -*- coding: utf-8 -*-
# __author__ = 'litangjian, liubinxu'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import shutil
import pandas as pd
from biocluster.core.exceptions import OptionError
from mbio.packages.denovo_rna.gene_structure.pfam_domtblout import pfam_out
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Seq
from mbio.packages.denovo_rna_v2.trans_step import step_count
from mbio.packages.rna.annot_config import AnnotConfig


class PredictAgent(Agent):
    """
    transdecoder:Predict,cds预测软件
    hmmscan：Pfam数据库比对工具
    version 1.0
    author: qindanhua
    last_modify: 2016.07.11
    """

    def __init__(self, parent):
        super(PredictAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "ref_rna_v2.fasta"},  # 输入文件
            {"name": "search_pfam", "type": "bool", "default": True},  # 是否比对Pfam数据库，cds预测涉及到这个吗？
            {"name": "Markov_length", "type": "int", "default": 3000},  # 马尔科夫训练长度
            {"name": "bed", "type": "outfile", "format": "ref_rna_v2.bed"},  # 输出结果
            {"name": "blast_cds", "type": "outfile", "format": "ref_rna_v2.fasta"},  # 输出结果
            {"name": "blast_pep", "type": "outfile", "format": "ref_rna_v2.fasta"},  # 输出结果
            {"name": "cds", "type": "outfile", "format": "ref_rna_v2.fasta"},  # 输出结果
            {"name": "pep", "type": "outfile", "format": "ref_rna_v2.fasta"},  # 输出结果
            {"name": "species_type", "type": "string"},  # 最开始workflow的参数是选择植物还是动物，这个决定我们是用planttfdb还是animaltfdb比对最后的结果
            {"name": "isoform_unigene", "type": "string"},  #实际为一个文件，只是不检查
            {"name": "blast_orf", "type": "string"},
            {"name": "pfam.domtblout", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "pfam.tblout", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "pfam_domain", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "pfam_version", "type": "string", "default": "32"},
        ]
        self.add_option(options)
        self._memory_increase_step = 50
        self.step.add_steps('predict')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.predict.start()
        self.step.update()

    def step_end(self):
        self.step.predict.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        """
        if not self.option("fasta").is_set:
            raise OptionError("请传入fasta序列文件", code = "33702401")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            # ["./estimators.xls", "xls", "alpha多样性指数表"]
        ])
        result_dir.add_regexp_rules([
            [r"transdecoder.pep$", "fasta", "蛋白质序列文件"],
            [r"transdecoder.cds$", "fasta", "cds序列文件"],
            [r"transdecoder.bed$", "bed", "orf位置信息bed格式文件"]
        ])
        if self.option("search_pfam"):
            result_dir.add_regexp_rules([
                ["./pfam_domain", "", "Pfam比对蛋白域结果信息"]
            ])

        if self.option("species_type").lower() == "animal":
            result_dir.add_regexp_rules([
                ["./animal_transcription_analysis_details", "", "动物转录因子分析详情"]
            ])

        if self.option("species_type").lower() == "plant":
            result_dir.add_regexp_rules([
                ["./plant_transcription_analysis_details", "", "植物转录因子分析详情"]
            ])

        # print self.get_upload_files()
        super(PredictAgent, self).end()


class PredictTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(PredictTool, self).__init__(config)
        # self.transdecoder_path = "bioinfo/gene-structure/TransDecoder-3.0.0/"
        self.transdecoder_path = "bioinfo/ref_rna_v2/miniconda2/bin/"
        self.hmmscan_path = "bioinfo/align/hmmer-3.1b2-linux-intel-x86_64/binaries/"
        self.pfam_db = AnnotConfig().get_file_path(
            file ="Pfam-A",
            db = "pfam",
            version = self.option("pfam_version"))
        # self.pfam_db = self.config.SOFTWARE_DIR + "/database/Pfam/Pfam-A.hmm"  # hmm参考库
        self.fasta_name = self.option("fasta").prop["path"].split("/")[-1]  # 得到fasta文件名
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        # TransDecoder-3.0.0的运行依赖perl的环境，所以程序必须设置一个环境变量，直接在我们的服务器程序里测试单条是不行的，我们没有perl_env
        self.perl_lib = self.config.SOFTWARE_DIR + '/program/perl/perls/perl-5.24.0/bin'
        self.set_environ(PATH=self.perl_lib)
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.tfdb_path = self.config.SOFTWARE_DIR + "/database/TFDB/"

    def td_predict(self, hmm_out=None):
        pre_dir = self.work_dir + "/" + os.path.basename(self.option("fasta").prop["path"]) + ".transdecoder_dir"
        if not os.path.exists(pre_dir):
            os.mkdir(pre_dir)
            for file in os.listdir(self.work_dir + "/../Predict/" + os.path.basename(self.option("fasta").prop["path"]) + ".transdecoder_dir"):
                os.link(self.work_dir + "/../Predict/" + os.path.basename(self.option("fasta").prop["path"]) + ".transdecoder_dir/" + file, pre_dir + "/" + file)
        if hmm_out is None:
            hmm_out = ""
        else:
            hmm_out = "--retain_pfam_hits {}".format(hmm_out)
        cmd = "{}TransDecoder.Predict -t {} -T {} {}".format(
            self.transdecoder_path, self.option("fasta").prop["path"], self.option("Markov_length"), hmm_out)
        print(cmd)
        self.logger.info("开始预测编码区域")
        self.logger.info(cmd)
        command = self.add_command("transdecoder_predict", cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("预测编码区域完成！")
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("预测过程过程出错", code = "33702402")


    def merge_and_change_id(self):
        # 获取预测结果

        predict_orfdict = dict()
        trans_predict_num = dict()
        old2new = dict()


        with open(os.path.basename(self.option("fasta").prop["path"]) + ".transdecoder.bed", 'r') as f_bed:
            f_bed.readline()
            for line in f_bed:
                cols = line.strip().split("\t")
                if cols[0] in trans_predict_num:
                    trans_predict_num[cols[0]] += 1
                else:
                    trans_predict_num[cols[0]] = 1

                new_name = cols[0] + "_orfp" +  str(trans_predict_num[cols[0]])
                old_id = cols[3].split(";")[0].split("=")[1]
                cds_type = cols[3].split(";")[2].split(":")[1]
                cds_type = cds_type.split("_len")[0]
                # print "old_id is {}".format(old_id)
                old2new[old_id] = new_name
                predict_orfdict[new_name] = {
                    "tran_id": cols[0],
                    "new_id": new_name,
                    "start": str(int(cols[6]) + 1),
                    "end": cols[7],
                    "strand": cols[5],
                    "des": cols[3],
                    "cds_type": cds_type
                }

        with open(self.option("blast_orf"), 'r') as f_blast:
            for line in f_blast:
                cols = line.strip().split("\t")
                predict_orfdict[cols[0]] = {
                    "tran_id": cols[0].split("_orf")[0],
                    "new_id": cols[0],
                    "start": cols[2],
                    "end": str(int(cols[3])),
                    "strand": cols[6],
                    "des": cols[1],
                    "cds_type": cols[5]
                }

        tran2gene = self.get_unigene_pep(predict_orfdict)


        # 合并提取的orf和预测的ORF 替换预测id

        record = SeqIO.parse(open(self.option("blast_cds").prop['path']), "fasta")
        SeqIO.write(record, "all_predicted.cds.fa", "fasta")
        with open("all_predicted.cds.fa", 'a') as f_cds:
            record = SeqIO.parse(os.path.basename(self.option("fasta").prop["path"]) + ".transdecoder.cds", "fasta")
            for seq in record:
                seq.id = old2new[seq.id]
                SeqIO.write(seq, f_cds, "fasta")

        record = SeqIO.parse("all_predicted.cds.fa", "fasta")
        with open("all_predicted.cds.unigene.fa", 'w') as f_cds:
            for seq in record:
                tran_id = predict_orfdict[seq.id]['tran_id']
                if tran2gene[tran_id][1] == "yes":
                    f_cds.write(">{}\n{}\n".format(seq.id, str(seq.seq)))

        record = SeqIO.parse(open(self.option("blast_pep").prop['path']), "fasta")
        SeqIO.write(record, "all_predicted.pep.fa", "fasta")
        with open("all_predicted.pep.fa", 'a') as f_cds:
            record = SeqIO.parse(os.path.basename(self.option("fasta").prop["path"]) + ".transdecoder.pep", "fasta")
            for seq in record:
                seq.id = old2new[seq.id]
                SeqIO.write(seq, f_cds, "fasta")

        # 生成ORF文件
        with open("all_predicted.xls", 'w') as f_detail,  open("all_predicted.bed", 'w') as f_bed:
            f_detail.write("protein_id\ttranscript_id\tgene_id\tis_gene\tposition\tcds_len\ttype\n")
            for pep_id, pep in sorted(predict_orfdict.items(), key=lambda x:x[0]):
                f_detail.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    pep_id,
                    pep['tran_id'],
                    tran2gene[pep['tran_id']][0],
                    tran2gene[pep['tran_id']][1],
                    pep['start'] + "-" + pep['end'] + "(" + pep['strand'] + ")",
                    str(int(pep['end']) - int(pep['start']) + 1),
                    pep['cds_type']

                ))
                f_bed.write("\t".join([
                    pep['tran_id'],
                    "0",
                    tran2gene[pep['tran_id']][2],
                    "ID={}".format(tran2gene[pep['tran_id']][0]),
                    "0",
                    pep['strand'],
                    str(int(pep['start']) - 1),
                    pep['end'],
                    "0",
                    "1",
                    tran2gene[pep['tran_id']][2],
                    "0"
                ]) + "\n")

        with open(self.option("pfam.tblout").prop['path']) as tbout, \
             open(self.work_dir + "/best_1_domain.tblout", "w") as best_1_domain:
            best_1_domain_line = set()
            best_1_domain.write("\t".join(["target_name", "accession", "query_name", "E_value", "score",
                                           "description_of_target"]) + "\n")

            for line in tbout:
                if line.startswith("#"):
                    continue
                line = line.strip().split()
                pep_id = line[2]
                if pep_id in old2new:
                    pep_id = old2new[pep_id]
                if pep_id not in predict_orfdict:
                    continue
                tran_id = predict_orfdict[pep_id]['tran_id']

                if tran_id not in best_1_domain_line:
                    best_1_domain.write("\t".join([line[0], line[1].split(".")[0], tran_id, line[4], line[5], (" ".join(line[18:]))]) + "\n")
                    best_1_domain_line.add(tran_id)

        with open(self.option("pfam_domain").prop['path'], 'r') as tbout, \
             open(self.work_dir + "/pfam_domain.out", "w") as domain:
            header = tbout.readline()
            domain.write(header)
            for line in tbout:
                line = line.split("\t")
                pep_id = line[1]
                if pep_id in old2new:
                    pep_id = old2new[pep_id]
                if pep_id not in predict_orfdict:
                    continue
                tran_id = predict_orfdict[pep_id]['tran_id']
                line[0] = tran_id
                line[1] = pep_id
                line[2] = line[2].split(".")[0]
                domain.write("\t".join(line))


    def get_unigene_pep(self, predict_orfdict):
        '''
        生成基因转录本对应关系文件
        '''
        bed = '{}.transdecoder.bed'.format(self.fasta_name)
        trans2pep = dict()
        trans2pep_len = dict()
        with open(bed, 'rb') as bed_file:
            for line in bed_file.readlines():
                cols = line.strip().split("\t")
                if len(cols) > 4 and cols[3].startswith("ID="):
                    tran_id = cols[0]
                    pep_id = cols[3].split("=")[1].split(";")[0]
                    pep_len = cols[3].split(":")[-1].split("_")[0]


        for pep_id,pep in predict_orfdict.items():
            tran_id = pep['tran_id']
            pep_len = abs(int(pep['end']) - int(pep['start']))
            if trans2pep.has_key(tran_id):
                trans2pep[tran_id] += ";" + pep_id
            else:
                trans2pep[tran_id] = pep_id

            if trans2pep_len.has_key(tran_id):
                if trans2pep_len[tran_id] > int(pep_len):
                    pass
                else:
                    trans2pep_len[tran_id] = int(pep_len)
            else:
                 trans2pep_len[tran_id] = int(pep_len)

        gene_list = list()
        gene_trans = dict()
        gene_unigene = dict()
        trans_len = dict()
        with open(self.option("isoform_unigene"), 'rb') as isoform_unigene:
            for line in isoform_unigene.readlines():
                cols = line.split("\t")
                '''
                if cols[2] == "yes":
                    gene_list.append(cols[1])
                    gene_unigene[cols[1]] = cols[0]
                '''
                trans_len[cols[0]] = int(cols[3])
                if cols[1] in gene_trans:
                    gene_trans[cols[1]].append(cols[0])
                else:
                    gene_trans[cols[1]] = [cols[0]]

        '''
        # 优先考虑蛋白长度
        for gene, trans in gene_trans.items():
            pep_lens = [[tran, trans2pep_len[tran]] for tran in trans if tran in trans2pep_len]
            if len(pep_lens) > 0:
                unigene = max(pep_lens, key=lambda x:x[1])[0]
                gene_unigene[gene] = unigene
            else:
                tran_lens = [[tran, trans_len[tran]] for tran in trans if tran in trans_len]
                unigene = max(tran_lens, key=lambda x:x[1])[0]
                gene_unigene[gene] = unigene
        '''

        

        '''
        for gene, unigene in gene_unigene.items():
            trans = gene_trans[gene]
            lens = [[tran, trans2pep_len[tran]] for tran in trans if tran in trans2pep_len]
            if len(lens) != 0:
                print "lens is {}".format(lens)
                unigene = max(lens, key=lambda x:x[1])[0]
                gene_unigene[gene] = unigene
        '''

        print gene_unigene
        tran2gene = dict()
        with open(self.option("isoform_unigene"), 'rb') as isoform_unigene, open("isoform_unigene.out", 'wb') as isoform_unigene2:
            for line in isoform_unigene.readlines():
                cols = line.strip().split("\t")
                if len(cols) <= 4:
                    cols.append("")

                '''
                if cols[1] in gene_unigene  and  cols[0] == gene_unigene[cols[1]]:
                    tran2gene[cols[0]] = [cols[1], "yes", cols[3]]
                    cols[2] = "yes"
                else:
                    tran2gene[cols[0]] = [cols[1], "no", cols[3]]
                    cols[2] = "no"
                '''
                tran2gene[cols[0]] = [cols[1], cols[2], cols[3]]
                
                if trans2pep.has_key(cols[0]):
                    cols[4] = trans2pep[cols[0]]
                else:
                    cols[4] = ""
                isoform_unigene2.write("{}\n".format("\t".join(cols)))
        return tran2gene





    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        # self.get_unigene_pep()
        self.logger.info("将结果文件link到output文件夹下面")
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        # pep = '{}.transdecoder.pep'.format(self.fasta_name)
        # bed = '{}.transdecoder.bed'.format(self.fasta_name)
        # cds = '{}.transdecoder.cds'.format(self.fasta_name)

        pep = 'all_predicted.pep.fa'
        bed = 'all_predicted.bed'
        cds = 'all_predicted.cds.fa'
        detail = 'all_predicted.xls'
        t2g = 'isoform_unigene.out'

        os.link(self.work_dir+"/isoform_unigene.out", self.output_dir+"/all_tran2gen.txt")

        os.link(self.work_dir+"/"+pep, self.output_dir+"/"+pep)
        self.option('pep').set_path(self.output_dir+"/"+pep)

        os.link(self.work_dir+"/"+bed, self.output_dir+"/"+bed)
        self.option('bed').set_path(self.output_dir+"/"+bed)

        os.link(self.work_dir+"/"+cds, self.output_dir+"/"+cds)
        self.option('cds').set_path(self.output_dir+"/"+cds)


        os.link(self.work_dir+"/"+detail, self.output_dir+"/"+detail)
        os.link(self.work_dir+"/"+t2g, self.output_dir+"/"+t2g)

        os.link(self.work_dir+"/"+'total_len_cds_list_only_unigene.txt', self.output_dir+"/"+'cds_len_unigene.txt')
        os.link(self.work_dir+"/"+'total_len_cds_list_only_transcript.txt', self.output_dir+"/"+'cds_len_transcript.txt')

        self.logger.info("link所有的非判断文件")

        if self.option("search_pfam"):
            if os.path.exists(self.work_dir + "/pfam_domain.out"):
                os.link(self.work_dir + "/pfam_domain.out", self.output_dir + "/pfam_domain")
            else:
                os.link(os.path.dirname(self.option("pfam.domtblout").prop['path']) + "/pfam_domain", self.output_dir + "/pfam_domain")
            self.logger.info("这次分析是用pfam库作为比对分析")

        if self.option("species_type").lower() == "animal":
            os.link(self.work_dir + "/" + 'merge_only_unigene_animal', self.output_dir+"/"+"unigene_tf_detail.xls")
            os.link(self.work_dir + "/" + 'merge_only_unigene_animal.stat', self.output_dir+"/"+"unigene_tf_stat.xls")
            os.link(self.work_dir + "/" + 'merge_only_transcript_animal', self.output_dir+"/"+"transcript_tf_detail.xls")
            os.link(self.work_dir + "/" + 'merge_only_transcript_animal.stat', self.output_dir+"/"+"transcript_tf_stat.xls")
            self.logger.info("这次分析是用分析的动物转录因子")

        if self.option("species_type").lower() == "plant":
            os.link(self.work_dir + "/" + 'merge_only_unigene_plant', self.output_dir+"/"+"unigene_tf_detail.xls")
            os.link(self.work_dir + "/" + 'merge_only_unigene_plant.stat', self.output_dir+"/"+"unigene_tf_stat.xls")
            os.link(self.work_dir + "/" + 'merge_only_transcript_plant', self.output_dir+"/"+"transcript_tf_detail.xls")
            os.link(self.work_dir + "/" + 'merge_only_transcript_plant.stat', self.output_dir+"/"+"transcript_tf_stat.xls")
            self.logger.info("这次分析是用分析的植物转录因子")

        self.logger.info("成功将结果文件link到output文件夹下面")

    def orf_stat(self):
        if self.option("search_pfam") is True:
            # self.td_predict(self.option("pfam.domtblout").prop['path'])
            self.logger.info("123456")
            # 统计cds的长度区间
            if not os.path.exists(self.work_dir + "/all_predicted.cds.fa"):
                cds_seq = self.work_dir + "/../Predict/" + self.fasta_name + ".transdecoder.pep"
            else:
                cds_seq = self.work_dir + "/all_predicted.cds.fa"
                cds_seq_unigene = self.work_dir + "/all_predicted.cds.unigene.fa"

            if os.path.exists(self.work_dir + "/total_len_cds_list_only_transcript.txt"):
                os.remove(self.work_dir + "/total_len_cds_list_only_transcript.txt")
            if os.path.exists(self.work_dir + "/total_len_cds_list_only_unigene.txt"):
                os.remove(self.work_dir + "/total_len_cds_list_only_unigene.txt")
            step_count(cds_seq, self.work_dir + "/trans_cds.txt" , 10, 200, self.work_dir + "/total_len_cds_list_only_transcript.txt")
            step_count(cds_seq_unigene, self.work_dir + "/trans_cds.txt" , 10, 200, self.work_dir + "/total_len_cds_list_only_unigene.txt")


            '''

            with open(cds_seq, "r") as cds_seq, open(self.work_dir + "/total_len_cds_list_only_unigene.txt","w") as total_cds_unigene, open(self.work_dir + "/total_len_cds_list_only_transcript.txt", "w") as total_cds_transcript:



                line_new = ""
                for line in cds_seq:
                    line_new += line.strip("\n")
                #print line_new
                line_split = line_new.strip().split(">")[1:]
                cds_whole_list = list()
                for line1 in line_split:
                    line1 = line1.split()
                    query_name, cds_len, cds_start, cds_end, cds_seq_concrete = line1[0].split("::")[1], int(line1[4].split(":")[
                        1]), \
                    line1[6].split(":")[1].split("(")[0].split("-")[0], line1[6].split(":")[1].split("(")[0].split("-")[1], \
                    line1[6].split(")")[1]
                    # 下面是统计整个的相关信息
                    # cds_single_list = [query_name, cds_len, cds_start, cds_end, cds_seq_concrete]
                    cds_single_list = [query_name, cds_len]
                    cds_whole_list.append(cds_single_list)
                total_cds_len_list = pd.DataFrame(cds_whole_list)  # 这时候是以行为标准写入的
                # 这个是把序列的信息也导入这个list
                # data.rename(columns={0:"query_name", 1:"cds_len", 2:"cds_start", 3:"cds_end", 4:"cds_seq_concrete"},inplace=True)
                # print(data)
                # 我们现在的结果是只需要统计长度的信息
                total_cds_len_list.rename(columns={0:"protein", 1:"cds_len"},inplace=True)
                total_cds_len_list.to_csv(self.work_dir + "/total_cds_len_list", sep = '\t', index=False)

                # 下面和那个带着yes或者no的文件进行比对，得到长度的分段统计
                isoform_unigene = pd.read_table("isoform_unigene.out", header=None, sep='\t')
                isoform_unigene.columns = ["transcript", "gene", "yes_no", "length", "protein"]
                # 筛选出只有yes的文件
                isoform_only_unigene = isoform_unigene[isoform_unigene.yes_no.isin(['yes'])][["transcript", "gene", "yes_no", "protein"]]
                cds_only_unigene = pd.merge(total_cds_len_list, isoform_only_unigene, on="protein")
                #print cds_only_unigene
                # 按照区间范围统计cds的长度范围,文本的读入都是字符串，所以前面需要类型转换下，或者df[数字列].astype(int, inplace=True)
                #cds_0_2_200 = cds_only_unigene[(cds_only_unigene["cds_len"]>=0) & (cds_only_unigene["cds_len"]<=200)]
                # 这个起始和终止坐标是有规律的，所以可以这么写，2个参数循环是嵌套循环的[[x,y] for x in range(2) for y in range(2)][[0, 0], [0, 1], [1, 0], [1, 1]]
                len_cds_list_only_unigene = [cds_only_unigene[cds_only_unigene["cds_len"].between(x+1, x+200,inclusive=True)].shape[0] for x in
                 range(0,1800,200)]
                # 长度为1800到无穷大的需要单独考虑计算
                cds_1800_2_inf = cds_only_unigene[cds_only_unigene["cds_len"]>=1801]
                len_1800_2_inf = cds_1800_2_inf.shape[0]
                len_cds_list_only_unigene.append(len_1800_2_inf)
                total_len_cds_list_only_unigene = sum(len_cds_list_only_unigene)
                total_cds_unigene.write("0~200" + "\t" + str(len_cds_list_only_unigene[0]) + "\n" +
                                        "201_400" + "\t" + str(len_cds_list_only_unigene[1]) + "\n" +
                                        "401_600" + "\t" + str(len_cds_list_only_unigene[2]) + "\n" +
                                        "601_800" + "\t" + str(len_cds_list_only_unigene[3]) + "\n" +
                                        "801_1000" + "\t" + str(len_cds_list_only_unigene[4]) + "\n" +
                                        "1001_1200" + "\t" + str(len_cds_list_only_unigene[5]) + "\n" +
                                        "1201_1400" + "\t" + str(len_cds_list_only_unigene[6]) + "\n" +
                                        "1401_1600" + "\t" + str(len_cds_list_only_unigene[7]) + "\n" +
                                        "1601_1800" + "\t" + str(len_cds_list_only_unigene[8]) + "\n" +
                                        ">1800" + "\t" + str(len_cds_list_only_unigene[9]) + "\n" +
                                        "total" + "\t" + str(total_len_cds_list_only_unigene))


                # 这里计算所有的转录本的cds长度，所有的转录本均计算，不对yes或者no筛选
                # isoform_unigene = pd.read_table("isoform_unigene.xls", header=None, sep = '\t')
                # isoform_unigene.columns = ["query_name", "unigene", "yes_no"]
                cds_only_transcript = pd.merge(total_cds_len_list, isoform_unigene, on="protein")
                # 按照区间范围统计cds的长度范围,文本的读入都是字符串，所以前面需要类型转换下，或者df[数字列].astype(int, inplace=True)
                #cds_0_2_200 = cds_only_transcript[(cds_only_transcript["cds_len"]>=0) & (cds_only_transcript["cds_len"]<=200)]
                # 这个起始和终止坐标是有规律的，所以可以这么写，2个参数循环是嵌套循环的[[x,y] for x in range(2) for y in range(2)][[0, 0], [0, 1], [1, 0], [1, 1]]
                len_cds_list_only_transcript = [cds_only_transcript[cds_only_transcript["cds_len"].between(x+1, x+200,inclusive=True)].shape[0] for x in
                 range(0,1800,200)]
                # 长度为1800到无穷大的需要单独考虑计算
                cds_1800_2_inf = cds_only_transcript[cds_only_transcript["cds_len"]>=1801]
                len_1800_2_inf = cds_1800_2_inf.shape[0]
                len_cds_list_only_transcript.append(len_1800_2_inf)
                total_len_cds_list_only_transcript = sum(len_cds_list_only_transcript)
                total_cds_transcript.write("0~200" + "\t" + str(len_cds_list_only_transcript[0]) + "\n" +
                "201_400" + "\t" + str(len_cds_list_only_transcript[1]) + "\n"
                "401_600" + "\t" + str(len_cds_list_only_transcript[2]) + "\n"
                "601_800" + "\t" + str(len_cds_list_only_transcript[3]) + "\n"
                "801_1000" + "\t" + str(len_cds_list_only_transcript[4]) + "\n"
                "1001_1200" + "\t" + str(len_cds_list_only_transcript[5]) + "\n"
                "1201_1400" + "\t" + str(len_cds_list_only_transcript[6]) + "\n"
                "1401_1600" + "\t" + str(len_cds_list_only_transcript[7]) + "\n"
                "1601_1800" + "\t" + str(len_cds_list_only_transcript[8]) + "\n"
                ">1800" + "\t" + str(len_cds_list_only_transcript[9]) + "\n"
                "total" + "\t" + str(total_len_cds_list_only_transcript))


            with open(self.option("pfam.tblout").prop['path']) as tbout, \
                 open(self.work_dir + "/best_1_domain.tblout", "w") as best_1_domain:
                # 我们只保留best 1 domain最高的那个数据对应的PF号，原始软件跑出来的数据就是从小到大的E_value排好序的，所以我们通过集合的不重复性进行元素的去重，也不是一整行重复，只是第一个accession的
                # 那一列重复，我们取::分割的第一个数据就可以判断
                best_1_domain_line = set()
                best_1_domain.write("\t".join(["target_name", "accession", "query_name", "E_value", "score",
                                               "description_of_target"]) + "\n")
                for line in tbout:
                    if line.startswith("#"):
                        continue
                    line = line.strip().split()
                    if "::" in line[2]:
                        if line[2].split("::")[0] not in best_1_domain_line:
                            best_1_domain.write("\t".join([line[0], line[1].split(".")[0], line[2].split("::")[1], line[4], line[5], (" ".join(line[18:]))]) + "\n")
                            best_1_domain_line.add(line[2].split("::")[0])
                    if line[2].split("_orf")[0] not in best_1_domain_line:
                            best_1_domain.write("\t".join([line[0], line[1].split(".")[0], line[2].split("_orf")[0], line[4], line[5], (" ".join(line[18:]))]) + "\n")
                            best_1_domain_line.add(line[2].split("_orf")[0])
            # 在测试的时候单条命令都可以通过，但是如果best_1_domain.tblout这个读在上面with的缩进里面，就会出现在流程里面报错，no columns to parse,所以即使pandas读有自己的读写
            # 但是我们还是要在with打开的写结束后再读取，这样更安全
        '''

    def tf_predict(self):
        if self.option("species_type").lower() == "animal":
            # 通过pandas的merge函数通过PF号来连接动物的2个文件
            pd_best = pd.read_table(self.work_dir + "/best_1_domain.tblout", header=0, sep = '\t')
            pd_animal = pd.read_table(self.tfdb_path + "animaltfdb.txt", header=0, sep = '\t')
            merge_animal = pd.merge(pd_best, pd_animal, on = "accession")
            merge_animal.to_csv(self.work_dir + "/merge_animal", sep = '\t', index=False)

            pd_merge_animal = pd.read_table(self.work_dir + "/merge_animal", header= 0, sep = '\t')
            pd_gene2tran = pd.read_table("isoform_unigene.out", header=None, sep = '\t')

            # 给传过来的文件增加一个行名字，作为匹配行，得到动物的所有转录本的转录因子的信息，要求这个文件传入三列，如果是四列文件，需要删除不需要的那一列
            pd_gene2tran.columns = ["query_name", "unigene", "yes_no", "len", "pep"]

            # 对动物的转录因子和unigene和transcript合并
            merge_gene2tran_animal = pd.merge(pd_merge_animal, pd_gene2tran, on="query_name")
            merge_gene2tran_animal.to_csv(self.work_dir + "/merge_gene2tran_animal", sep = '\t', index=False)

            # 这里筛选出第三列为yes的，也就是只有unigene的转录因子的信息
            pd_only_unigene = pd_gene2tran[pd_gene2tran.yes_no.isin(['yes'])]

            #　得到动物unigene转录因子信息
            merge_only_unigene_animal = pd.merge(pd_merge_animal, pd_only_unigene, on="query_name")
            merge_only_unigene_animal.to_csv(self.work_dir + "/merge_only_unigene_animal", sep = '\t', index=False)
            tf_trans = dict()
            for tran, fam in  zip(merge_only_unigene_animal["unigene"], merge_only_unigene_animal["Family"]):
                if fam in tf_trans:
                    tf_trans[fam].append(tran)
                else:
                    tf_trans[fam] = [tran]
            with open("merge_only_unigene_animal.stat", 'w') as f:
                f.write("tf_family\tnum\tseq_list\n")
                for fam in sorted(tf_trans.keys(), key=lambda x:len(tf_trans[x]), reverse=True):
                    f.write("{}\t{}\t{}\n".format(fam, len(tf_trans[fam]), ",".join(tf_trans[fam])))
        
                

            # 得到动物所有转录本的转录因子信息
            merge_only_transcript_animal = pd.merge(pd_merge_animal, pd_gene2tran, on="query_name")
            merge_only_transcript_animal.to_csv(self.work_dir + "/merge_only_transcript_animal", sep = '\t', index=False)

            tf_trans = dict()
            for tran, fam in  zip(merge_only_transcript_animal["query_name"], merge_only_transcript_animal["Family"]):
                if fam in tf_trans:
                    tf_trans[fam].append(tran)
                else:
                    tf_trans[fam] = [tran]
            with open("merge_only_transcript_animal.stat", 'w') as f:
                f.write("tf_family\tnum\tseq_list\n")
                for fam in sorted(tf_trans.keys(), key=lambda x:len(tf_trans[x]), reverse=True):
                    f.write("{}\t{}\t{}\n".format(fam, len(tf_trans[fam]), ",".join(tf_trans[fam])))


        elif self.option("species_type").lower() == "plant":
            # 通过pandas的merge函数通过PF号来连接植物的2个文件
            pd_best = pd.read_table(self.work_dir + "/best_1_domain.tblout", header=0, sep = '\t')
            pd_plant = pd.read_table(self.tfdb_path + "planttfdb.txt", header=0, sep = '\t')
            merge_plant = pd.merge(pd_best, pd_plant, on="accession")
            merge_plant.to_csv(self.work_dir + "/merge_plant", sep = '\t', index=False)


            pd_merge_plant = pd.read_table(self.work_dir + "/merge_plant", header= 0, sep = '\t')
            pd_gene2tran = pd.read_table("isoform_unigene.out", header=None, sep = '\t')

            # 给传过来的文件增加一个行名字，作为匹配行，得到动物的所有转录本的转录因子的信息，要求这个文件传入三列，如果是四列文件，需要删除不需要的那一列
            pd_gene2tran.columns = ["query_name", "unigene", "yes_no", "len", "pep"]


            # 对植物的转录因子和unigene和transcript合并
            merge_gene2tran_plant = pd.merge(pd_merge_plant, pd_gene2tran, on="query_name")
            merge_gene2tran_plant.to_csv(self.work_dir + "/merge_gene2tran_plant", sep = '\t', index=False)
        

            # 这里筛选出第三列为yes的，也就是只有unigene的转录因子的信息
            pd_only_unigene = pd_gene2tran[pd_gene2tran.yes_no.isin(['yes'])]

            #　得到植物unigene转录因子信息
            merge_only_unigene_plant = pd.merge(pd_merge_plant, pd_only_unigene, on="query_name")
            merge_only_unigene_plant.to_csv(self.work_dir + "/merge_only_unigene_plant", sep = '\t', index=False)

            tf_trans = dict()
            for tran, fam in  zip(merge_only_unigene_plant["unigene"], merge_only_unigene_plant["Family"]):
                if fam in tf_trans:
                    tf_trans[fam].append(tran)
                else:
                    tf_trans[fam] = [tran]
            with open("merge_only_unigene_plant.stat", 'w') as f:
                f.write("tf_family\tnum\tseq_list\n")
                for fam in sorted(tf_trans.keys(), key=lambda x:len(tf_trans[x]), reverse=True):
                    f.write("{}\t{}\t{}\n".format(fam, len(tf_trans[fam]), ",".join(tf_trans[fam])))
        

            # 得到植物所有转录本的转录因子信息
            merge_only_transcript_plant = pd.merge(pd_merge_plant, pd_gene2tran, on="query_name")
            merge_only_transcript_plant.to_csv(self.work_dir + "/merge_only_transcript_plant", sep = '\t', index=False)
            tf_trans = dict()
            for tran, fam in  zip(merge_only_transcript_plant["query_name"], merge_only_transcript_plant["Family"]):
                if fam in tf_trans:
                    tf_trans[fam].append(tran)
                else:
                    tf_trans[fam] = [tran]
            with open("merge_only_transcript_plant.stat", 'w') as f:
                f.write("tf_family\tnum\tseq_list\n")
                for fam in sorted(tf_trans.keys(), key=lambda x:len(tf_trans[x]), reverse=True):
                    f.write("{}\t{}\t{}\n".format(fam, len(tf_trans[fam]), ",".join(tf_trans[fam])))



        else:
            self.logger.info("我们目前只对动植物进行转录因子预测")

    def run(self):
        super(PredictTool, self).run()
        """将Blast和Pfam搜索结果整合到编码区域选择
        TransDecoder借助上面生成的输出结果来确定将这些被blast命中的和结构域命中的多肽保留在报告的编码区集合中。像这样运行TransDecoder.Predict：
        TransDecoder.Predict -t target_transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6
        最终的编码区预测结果将包含与编码区域一致的序列字符以及blast得到的直系同源结果或pfam结构域的内容。这里选择是否比对pfam库,
        这一步是为了做注释使用，hmmscan其实做转录因子其实是不需要做domtblout的结果输出的，我们流程默认的Predict是要整合pfam数据库的结果"""
        if self.option("search_pfam") is True:
            self.td_predict(self.option("pfam.domtblout").prop['path'])
            self.merge_and_change_id()
            # self.get_unigene_pep()
            self.orf_stat()
        else:
            # 未测试过不预测orf的情况
            self.td_predict()
            self.merge_and_change_id()
            # self.get_unigene_pep()
            self.orf_stat()

        if self.option("search_pfam") is True:
            self.tf_predict()

        self.set_output()
        self.end()
