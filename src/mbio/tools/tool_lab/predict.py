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
from mbio.packages.lnc_rna.copy_file import CopyFile


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
            {"name": "genetic_code", "type": "string", "default": "genetic_code"},
            {"name": "single_best_only", "type": "string", "default": "yes"},
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
        self._cpu = 20
        self._memory = '40G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            # ["./estimators.xls", "xls", "alpha多样性指数表"]
        ])
        result_dir.add_regexp_rules([
            [r"transdecoder.pep$", "fasta", "蛋白质序列文件"],
            [r"transdecoder.cds$", "fasta", "cds序列文件"],
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
        self.transdecoder_path = "miniconda2/bin/"
        self.hmmscan_path = "bioinfo/align/hmmer-3.1b2-linux-intel-x86_64/binaries/"
        self.pfam_db = self.config.SOFTWARE_DIR + "/database/Pfam/Pfam-A.hmm"  # hmm参考库
        self.fasta_name = self.option("fasta").prop["path"].split("/")[-1]  # 得到fasta文件名
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        # TransDecoder-3.0.0的运行依赖perl的环境，所以程序必须设置一个环境变量，直接在我们的服务器程序里测试单条是不行的，我们没有perl_env
        self.perl_lib = self.config.SOFTWARE_DIR + '/miniconda2/bin'
        self.set_environ(PATH=self.perl_lib)
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.tfdb_path = self.config.SOFTWARE_DIR + "/database/TFDB/"
        self.predicted = True

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
        elif command.return_code == 2:
            self.logger.info("序列数量少， 无法建立训练模型")
            self.predicted = False
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("预测过程过程出错", code = "33702402")


    def merge_and_change_id(self):
        # 获取预测结果a

        predict_orfdict = dict()
        trans_predict_num = dict()
        old2new = dict()


        if self.predicted:
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

        # tran2gene = self.get_unigene_pep(predict_orfdict)


        # 合并提取的orf和预测的ORF 替换预测id

        record = SeqIO.parse(open(self.option("blast_cds").prop['path']), "fasta")
        SeqIO.write(record, "all_predicted.cds.fa", "fasta")
        if self.predicted:
            with open("all_predicted.cds.fa", 'a') as f_cds:
                record = SeqIO.parse(os.path.basename(self.option("fasta").prop["path"]) + ".transdecoder.cds", "fasta")
                for seq in record:
                    seq.id = old2new[seq.id]
                    SeqIO.write(seq, f_cds, "fasta")


        record = SeqIO.parse(open(self.option("blast_pep").prop['path']), "fasta")
        SeqIO.write(record, "all_predicted.pep.fa", "fasta")
        if self.predicted:
            with open("all_predicted.pep.fa", 'a') as f_cds:
                record = SeqIO.parse(os.path.basename(self.option("fasta").prop["path"]) + ".transdecoder.pep", "fasta")
                for seq in record:
                    seq.id = old2new[seq.id]
                    SeqIO.write(seq, f_cds, "fasta")

        # 生成ORF文件
        with open("all_predicted.xls", 'w') as f_detail,  open("all_predicted.bed", 'w') as f_bed:
            f_detail.write("protein_id\tseq_id\tposition\tcds_len\ttype\n")
            for pep_id, pep in sorted(predict_orfdict.items(), key=lambda x:x[0]):
                f_detail.write("{}\t{}\t{}\t{}\t{}\n".format(
                    pep_id,
                    pep['tran_id'],
                    pep['start'] + "-" + pep['end'] + "(" + pep['strand'] + ")",
                    str(int(pep['end']) - int(pep['start']) + 1),
                    pep['cds_type']

                ))
                f_bed.write("\t".join([
                    pep['tran_id'],
                    "0",
                    pep['tran_id'],
                    "ID={}".format(pep['tran_id']),
                    "0",
                    pep['strand'],
                    str(int(pep['start']) - 1),
                    pep['end'],
                    "0",
                    "1",
                    pep['tran_id'],
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



    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        self.logger.info("将结果文件link到output文件夹下面")

        # pep = '{}.transdecoder.pep'.format(self.fasta_name)
        # bed = '{}.transdecoder.bed'.format(self.fasta_name)
        # cds = '{}.transdecoder.cds'.format(self.fasta_name)

        pep = 'all_predicted.pep.fa'
        cds = 'all_predicted.cds.fa'
        detail = 'all_predicted.xls'

        for each in [pep, cds, detail]:
            CopyFile().linkdir(os.path.join(self.work_dir, each), os.path.join(self.output_dir, each))




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
        else:
            # 未测试过不预测orf的情况
            self.td_predict()
            self.merge_and_change_id()


        self.set_output()
        self.end()
