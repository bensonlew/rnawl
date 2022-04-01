# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2019/4/29'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from Bio import SeqIO
import pandas as pd


class CircleCorrectAgent(Agent):
    """
    DemoAgent:
    version 1.0
    """

    def __init__(self, parent):
        super(CircleCorrectAgent, self).__init__(parent)
        options = [
            {"name": "blast_table", "type": "infile", "format": "align.blast.blast_table", "required": True},
            {"name": "query", "type": "infile", "format": "sequence.fasta", "required": True},
            {"name": "sample_name", "type": "string"},
            {"name": "fa_out", "type": "outfile", "format": "sequence.fasta"},
            {"name": "table_out", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "cir_out", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        # modified check_option
        pass

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = "10G"


class CircleCorrectTool(Tool):
    def __init__(self, config):
        super(CircleCorrectTool, self).__init__(config)
        self.blast = "bioinfo/align/ncbi-blast-2.3.0+/bin"
        self.dnaA = self.config.SOFTWARE_DIR + "/database/bacgenome/dnaA/dnaA.faa"
        self.blast_data = ""
        self.circle_status = {}
        self.fastas = []
        self.circle_table_list = []
        self.chrome_info = {}

    def run_read_fa(self):
        predict_dir = os.path.join(self.work_dir, "predict_dir")
        if not os.path.isdir(predict_dir):
            os.mkdir(predict_dir)
        self.parse_blast()
        index_count = 0
        tmp_fa_dict = {}
        tmp_data_record = []
        for seq_record in SeqIO.parse(self.option("query").prop["path"], "fasta"):
            index_count += 1
            seq_record_id = seq_record.id
            # seq.description 可能会用到
            # new_record_id = "Scaffold%s" % index_count
            # self.logger.info(self.blast_data)
            # self.blast_data.loc[self.blast_data["Query_id"] == seq_record_id, "Query_id"] = new_record_id
            # self.logger.info(self.blast_data.loc[self.blast_data["Query_id"] == new_record_id, "Query_id"])
            if len(seq_record) >= 10000:
                cut_len = len(seq_record) / 3 if len(seq_record) / 3 < 3000 else 3000
                left_seq = seq_record[:cut_len]  # 左端reads
                cov_start = self.compare_seq(left_seq, seq_record) # 右端序列从这个位置开始比对上的
                if cov_start:
                    seq_record = seq_record[:cov_start]  # 从cov_len+1位置开始取
                    query = os.path.join(predict_dir, seq_record_id + ".fa")
                    SeqIO.write(seq_record, query, "fasta")
                    dnaA_start, direction = self.run_dnaA(query, query + ".out")  # 增加比对方向，反向比对结果进行反向互补处理
                    if direction == 1:
                        # dnaA起始位置前需要保留100个碱基
                        seq_record = seq_record[dnaA_start-101:] + seq_record[:dnaA_start-101]  # 如果起始点少于101，可以向序列末端借碱基，代码是一样的
                    elif direction == -1:
                        # 反向互补需要对Seq进行处理，不能是Seq_record
                        if len(seq_record) - dnaA_start >= 100:
                            seq_record.seq = seq_record.seq[:dnaA_start+100].reverse_complement() + seq_record.seq[dnaA_start+100:].reverse_complement()
                        else:
                            seq_need_len = dnaA_start + 100 - len(seq_record)  # 末端不够100bp作为起始序列，到序列前端借
                            seq_record.seq = seq_record.seq[:seq_need_len].reverse_complement() + seq_record.seq[seq_need_len:].reverse_complement()
                    tmp_fa_dict[seq_record_id] = seq_record
                    tmp_data_record.append([seq_record_id, len(seq_record)])
                    self.blast_data.loc[self.blast_data["Query_id"]==seq_record_id, "chrome_info"] = "Circular"
                    self.circle_status[seq_record_id] = "Circular"
                else:
                    tmp_fa_dict[seq_record_id] = seq_record
                    tmp_data_record.append([seq_record_id, len(seq_record)])
                    self.blast_data.loc[self.blast_data["Query_id"]==seq_record_id, "chrome_info"] = "Linear"
                    self.circle_status[seq_record_id] = "Linear"
            else:
                tmp_fa_dict[seq_record_id] = seq_record
                tmp_data_record.append([seq_record_id, len(seq_record)])
                self.blast_data.loc[self.blast_data["Query_id"]==seq_record_id, "chrome_info"] = "Linear"
                self.circle_status[seq_record_id] = "Linear"

        number = 1
        for seq_id in pd.DataFrame(tmp_data_record).sort_values(by=[1], ascending=False)[0]:
            new_seq_id = "Scaffold%s" % number
            number += 1
            seq_record = tmp_fa_dict[seq_id]
            seq_record.id = new_seq_id
            self.fastas.append(seq_record)
            self.circle_table_list.append(["-", new_seq_id, new_seq_id, self.circle_status[seq_id]])
            if seq_id in self.blast_data["Query_id"].values:
                self.blast_data.loc[self.blast_data["Query_id"] == seq_id, "Query_id"] = new_seq_id

    def compare_seq(self, left_seq, ref_seq):
        compare_dir = os.path.join(self.work_dir, "compare_dir")
        if not os.path.isdir(compare_dir):
            os.mkdir(compare_dir)
        cmd_name = left_seq.id
        left_seq.id = left_seq.id + "_l"
        ref_seq.id = ref_seq.id
        left_path = os.path.join(compare_dir, left_seq.id + ".fa")
        ref_path = os.path.join(compare_dir, ref_seq.id + ".fa")
        SeqIO.write(left_seq, left_path, "fasta")
        SeqIO.write(ref_seq,ref_path, "fasta")
        cmd = os.path.join(self.blast, "makeblastdb")
        cmd += " -dbtype nucl -in %s -parse_seqids -out %s " % (ref_path, ref_path)
        makedb = self.add_command("mk%s" % cmd_name, cmd).run()
        self.wait(makedb)
        location = None
        if makedb.return_code == 0:
            out_path = os.path.join(self.work_dir, "compare_dir", cmd_name + ".out")
            location = self.run_blast(ref_path, left_path, out_path, len(ref_seq), len(left_seq))
        else:
            self.set_error("序列左端构建的blast数据库失败")
        return location

    def run_blast(self, ref, query, out, r_len, q_len):
        out = str(out)
        # if os.path.isfile(out):  # 为加快测试速度临时增加
        #     return self.get_loc(out, r_len, q_len)
        cmd = os.path.join(self.blast, "blastn") + " -query %s -db %s -out %s -evalue 1e-5 -outfmt 6 -num_threads 2" % (query, ref, out)
        name = os.path.basename(out)
        command = self.add_command("blast%s" % name, cmd).run()
        self.wait(command)
        if command.return_code == 0:
            return self.get_loc(out, r_len, q_len)  # 没有则返回None
        else:
            self.set_error("序列右端比对左端失败")

    def get_loc(self, blast_out, r_len, q_len):
        """
        get seq location, 用于截取序列， 获取比对上的start位置
        :param blast_out:
        :return:
        """
        if os.path.getsize(blast_out) == 0:
            return None
        data = pd.read_table(blast_out, header=0, names=["que", "ref", "iden", "cov", "mis", "gap", "qs", "qe", "rs", "re", "eva", "sco"])
        if data.empty:
            return None
        cov_max = data["cov"].max()
        data = data.loc[(data["rs"] < data["re"]) & (data["rs"] > q_len) & (data["re"] > r_len*0.65) & (data["cov"] > cov_max*0.65),]
        if data.empty:
            return None
        return data.loc[data["re"]==data["re"].max(), "rs"].iloc[0]
        # return data.loc[0, "rs"]

    def get_dnaA_loc(self, blast_out):
        if os.path.getsize(blast_out) == 0:
            return None, 0
        data = pd.read_table(blast_out, header=None)
        direction = 1 if data.loc[0, 6] < data.loc[0,7] else -1
        return data.loc[0, 6], direction

    def run_dnaA(self, query, out):
        out = str(out)
        # if os.path.isfile(out):
        #     return self.get_dnaA_loc(out)  # 为加快测试计算速度，临时增加
        cmd = os.path.join(self.blast, "blastx") + " -query %s -db %s -out %s -evalue 1e-5 -outfmt 6 -max_hsps 1 -num_threads 2 -max_target_seqs 1" % (query, self.dnaA, out)
        name = os.path.basename(out)
        command = self.add_command("dnaa%s" % name, cmd).run()
        self.wait(command)
        if command.return_code == 0:
            return self.get_dnaA_loc(out)
        else:
            self.set_error("dnaA比对出错")


    def parse_blast(self):
        """
        解析参数传过来的blast结果
        :return:
        """
        data = pd.read_table(self.option("blast_table").prop["path"], header=0, names=["Score", "Evalue", "align_length",
            "Identity (%)", "Similarity (%)", "Query_id", "q.length", "q.start", "q.end", "q.direction", "hit_name",
            "s.length", "s.start",  "s.end", "s.direction", "Description"])
        data["chrome_info"] = data.apply(self.is_chrom, axis=1)
        data["Coverage (%)"] = (data["q.end"] - data["q.start"]) / data["q.length"] * 100
        data["Subject ID"] = data.apply(lambda x: x["Description"].split(" ")[0], axis=1)
        data["Sample Name"] = self.option("sample_name")
        self.blast_data = data.reindex(columns=["Sample Name", "Query_id", "align_length", "q.length", "q.start", "q.end",
                                     "q.direction", "Subject ID", "s.length", "s.start", "s.end", "s.direction",
                                     "Description", "Identity (%)", "Coverage (%)", "Evalue", "Score", "chrome_info"])
        self.chrome_info = data.reindex(columns=["Query_id", "chrome_info"]).set_index("Query_id").to_dict()["chrome_info"]

    def is_chrom(self, df):
        if "plasmid" in df["Description"]:
            return "plasmid"
        elif "complete genome" in df["Description"]:
            return "chr"
        else:
            return "linear"

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        fasta_result = os.path.join(self.output_dir, self.option("sample_name") + ".scaffold.fna")
        SeqIO.write(self.fastas, fasta_result, "fasta")
        self.option("fa_out").set_path(fasta_result)
        table_result = os.path.join(os.path.join(self.output_dir, self.option("sample_name") + ".xls"))
        table_result2 = os.path.join(os.path.join(self.output_dir, self.option("sample_name") + "cir.xls"))
        self.blast_data.to_csv(table_result, sep="\t")
        pd.DataFrame(self.circle_table_list).to_csv(table_result2, sep="\t", header=False, index=False)
        self.option("table_out").set_path(table_result)
        self.option("cir_out").set_path(table_result2)

    def run(self):
        super(CircleCorrectTool, self).run()
        self.run_read_fa()
        self.set_output()
        self.end()