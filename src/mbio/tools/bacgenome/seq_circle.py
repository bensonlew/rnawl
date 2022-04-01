# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# __last_modify__ = '2020/07/08'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from Bio import SeqIO
import pandas as pd


class SeqCircleAgent(Agent):
    """
    DemoAgent:主要判断基因组的序列是环形还是线形；主要是将一条序列的左端和右端比较，左端的起始位置必须小于或等于10，右端的终点位置在[len(seq)-10,len(seq)]
    version 1.0
    """

    def __init__(self, parent):
        super(SeqCircleAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta", "required": True},
            {"name": "sample_name", "type": "string"},
            {"name": "fa_out", "type": "outfile", "format": "sequence.fasta"},
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


class SeqCircleTool(Tool):
    def __init__(self, config):
        super(SeqCircleTool, self).__init__(config)
        self.blast = "bioinfo/align/ncbi-blast-2.3.0+/bin"
        self.blast_data = ""
        self.circle_status = {}
        self.fastas = []
        self.circle_table_list = []
        self.chrome_info = {}

    def run_read_fa(self):
        predict_dir = os.path.join(self.work_dir, "predict_dir")
        if not os.path.isdir(predict_dir):
            os.mkdir(predict_dir)
        index_count = 0
        tmp_fa_dict = {}
        tmp_data_record = []
        for seq_record in SeqIO.parse(self.option("query").prop["path"], "fasta"):
            index_count += 1
            seq_record_id = seq_record.id
            if len(seq_record) >= 10000:
                cut_len = len(seq_record) / 3 if len(seq_record) / 3 < 3000 else 3000
                left_seq = seq_record[:cut_len]  # 左端reads
                right_seq = seq_record[len(seq_record)-cut_len:]
                cov_start = self.compare_seq(left_seq, right_seq) # 右端序列从这个位置开始比对上的
                if cov_start:
                    seq_record = seq_record[:(len(seq_record)-cut_len+cov_start)]  # 去除尾部的重复序列
                    query = os.path.join(predict_dir, seq_record_id + ".fa")
                    SeqIO.write(seq_record, query, "fasta")
                    tmp_fa_dict[seq_record_id] = seq_record
                    tmp_data_record.append([seq_record_id, len(seq_record)])
                    self.circle_status[seq_record_id] = "Circular"
                else:
                    tmp_fa_dict[seq_record_id] = seq_record
                    tmp_data_record.append([seq_record_id, len(seq_record)])
                    self.circle_status[seq_record_id] = "Linear"
            else:
                tmp_fa_dict[seq_record_id] = seq_record
                tmp_data_record.append([seq_record_id, len(seq_record)])
                self.circle_status[seq_record_id] = "Linear"

        number = 1
        for seq_id in pd.DataFrame(tmp_data_record).sort_values(by=[1], ascending=False)[0]:
            new_seq_id = "Scaffold%s" % number
            number += 1
            seq_record = tmp_fa_dict[seq_id]
            seq_record.id = new_seq_id
            self.fastas.append(seq_record)
            self.circle_table_list.append([new_seq_id, new_seq_id, self.circle_status[seq_id]])

    def compare_seq(self, left_seq, right_seq):
        compare_dir = os.path.join(self.work_dir, "compare_dir")
        if not os.path.isdir(compare_dir):
            os.mkdir(compare_dir)
        cmd_name = left_seq.id.lower()
        left_seq.id = "left_seq_l"
        right_seq.id = "right_seq_r"
        left_path = os.path.join(compare_dir, left_seq.id + ".fa")
        ref_path = os.path.join(compare_dir, right_seq.id + ".fa")
        SeqIO.write(left_seq, left_path, "fasta")
        SeqIO.write(right_seq,ref_path, "fasta")
        cmd = os.path.join(self.blast, "makeblastdb")
        cmd += " -dbtype nucl -in %s -parse_seqids -out %s " % (ref_path, ref_path)
        makedb = self.add_command("mk%s" % cmd_name, cmd).run()
        self.wait(makedb)
        location = None
        if makedb.return_code == 0:
            out_path = os.path.join(self.work_dir, "compare_dir", cmd_name + ".out")
            location = self.run_blast(ref_path, left_path, out_path, len(right_seq), len(left_seq))
        else:
            self.set_error("序列左端构建的blast数据库失败")
        return location

    def run_blast(self, ref, query, out, r_len, q_len):
        out = str(out)
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
        data = data.loc[(data["qs"] < data["qe"]) & (data["qs"] <= 10) & (data["re"] > data["rs"]) & (r_len-10 <= data["re"]) & (data["re"] <=r_len)]
        if data.empty:
            return None
        cov_max = data["cov"].max()
        if cov_max >= 1000:
            return data.loc[data["cov"]==cov_max, "rs"].iloc[0]

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        fasta_result = os.path.join(self.output_dir, self.option("sample_name") + ".scaffold.fna")
        SeqIO.write(self.fastas, fasta_result, "fasta")
        self.option("fa_out").set_path(fasta_result)
        table_result2 = os.path.join(os.path.join(self.output_dir, self.option("sample_name") + ".cir.log"))
        pd.DataFrame(self.circle_table_list).to_csv(table_result2, sep="\t", header=False, index=False)
        self.option("cir_out").set_path(table_result2)


    def run(self):
        super(SeqCircleTool, self).run()
        self.run_read_fa()
        self.set_output()
        self.end()