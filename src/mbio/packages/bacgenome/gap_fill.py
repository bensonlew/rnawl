# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import os
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd
from biocluster.file import download
from bson.objectid import ObjectId

my_output_fa_record = []
my_output_log_record = pd.DataFrame()

def parse_blast(table, aln_len, identity):
    # 读取存储blast表
    # index=["query", "ref", "identity", "aln_len", "mis", "gap", "qstart", "qend", "rstart", "rend", "evalue", "score"]
    data = pd.read_table(table, header=None, names=["query", "ref", "identity", "aln_len", "mis", "gap", "qstart", "qend", "rstart", "rend", "evalue", "score"])
    blast_data = data.loc[(data["aln_len"]>=aln_len) & (data["identity"]>=identity),]
    return blast_data

def parse_log(log_file):
    have_circle_info = False
    genome_info = {}
    parent_info = {}
    with open(log_file, "r") as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip().split("\t")
            if genome_info.has_key(line[0]):
                genome_info[line[0]].append(line[1])
            else:
                genome_info[line[0]] = [line[1]]
    data = pd.read_table(log_file, header=None)
    if len(data.columns) == 3:
        data.columns = ["genome", "seq", "parent"]
    else:
        have_circle_info = True
        data.columns = ["genome", "seq", "parent", "is_circle"]
    parent_info = data.set_index(["genome", "seq"]).to_dict(orient="index")
    return parent_info, genome_info, have_circle_info

def seq_is_circled(blast_data, seq_id):
    """
    判断一个定义染色体/质粒只包含一条序列时，是否成环
    一个定义染色体/质粒包含多条序列时，序列间补洞后，对于独立的序列，判断是否成环
    """
    left_data = blast_data.loc[blast_data["ref"] == seq_id + "_left", "query"]
    right_data = blast_data.loc[blast_data["ref"] == seq_id + "_right", "query"]
    left_query_set = set(left_data.tolist())
    right_query_set = set(right_data.tolist())
    gap_seq_list = list(left_query_set & right_query_set)
    if len(gap_seq_list) == 0:
        return False
    else:
        return True

def get_ref2(df, ref1_direction, ref2_blast):
    """
    八种情况符合补gap：
    1.ref1的右端序列正向比对到query, query比对区域的右端正向比对到ref2的左端序列
        此时query和ref2都是正向的
    2.ref1的右端序列正向比对到query, query比对区域的右端反向比对到ref2的右端序列
        此时query是正向， ref2是反向
    3.ref1的右端序列反向比对到query, query比对区域的左侧正向比对到ref2的右端序列
        此时query和ref2都是反向
    4.ref1的有段序列反向比对到query, query比对区域的左侧反向比对到ref2的左端序列
        此时query是反向, ref2是正向
    5.ref1的左端序列正向比对到query, query比对区域的左侧正向比对到ref2的右端序列
        此时query和ref2都是正向
    6.ref1的左端序列正向比对到query, query比对区域的左侧反向比对到ref2的左端序列
        此时query是正向，ref2是反向
    7.ref1的左端序列反向比对到query, query比对区域的右侧正向比对到ref2的左端序列
        此时query和ref2都是反向
    8.ref1的左端序列反向比对到query, query比对区域的右侧反向比对到ref2的右端序列
        此时query是反向，ref2是正向
    :param df:
    :param ref1_direction:
    :param ref2_blast:
    :return:
    """
    query_direction = True
    if ref1_direction=="right" and df["rstart"] < df["rend"]:
        """
        情况1,2
        找query 比对区域end右侧比对上ref2的结果
        ref2的左端正向比对到query接受（情况1）
        ref2的右端反向比对到query接受（情况2）
        """
        query_direction = "True" # query的相对ref1比对方向
        filtered = ref2_blast.loc[(ref2_blast["qstart"] > df["qend"]), ]
        if filtered.empty:
            return pd.Series()
        filtered = filtered.loc[( filtered["ref"].str.contains("left") & (filtered["rstart"] < filtered["rend"]) )|
                                ( filtered["ref"].str.contains("right") & (filtered["rstart"] > filtered["rend"]) ),]
        if filtered.empty:
            return pd.Series()
        ref2_result = filtered.loc[filtered["qstart"]==filtered["qstart"].min(),].iloc[0]  # 优先选择gap区最短的
        ref2_direction = "True" if ref2_result["rstart"] < ref2_result["rend"] else "False"
        return pd.Series({"q_direction": query_direction, "pos1": df["qend"], "pos2": ref2_result["qstart"], "ref2_direction": ref2_direction})
    elif ref1_direction=="right" and df["rstart"] > df["rend"]:
        """
        ref1的右端序列反向比对上query（情况3,4）
        此时首先找query 比对区域start左侧比对上ref2的结果
        ref2的右端正向比对到query接受（情况3）
        或者ref2的左端反向比对到query接受（情况4）
        """
        query_direction = "False"
        filtered = ref2_blast.loc[(ref2_blast["qend"] < df["qstart"]),]
        if filtered.empty:
            return pd.Series()
        filtered = filtered.loc[( filtered["ref"].str.contains("right")&(filtered["rstart"] < filtered["rend"])  ) |
                                ( filtered["ref"].str.contains("left")&(filtered["rstart"] > filtered["rend"])  ),]
        if filtered.empty:
            return pd.Series()
        ref2_result = filtered.loc[filtered["qend"]==filtered["qend"].max(),].iloc[0]
        ref2_direction = "True" if ref2_result["rstart"] > ref2_result["rend"] else "False"
        return pd.Series({"q_direction": query_direction, "pos1": df["qstart"], "pos2": ref2_result["qend"], "ref2_direction": ref2_direction})
    elif ref1_direction=="left" and df["rstart"] < df["rend"]:
        """
        ref1的左端序列正向比对上query（情况5,6）
        此时首先找query比对区域start左侧比对上ref2的结果
        ref2的右端正向比对到query接受（情况5）
        ref2的左端反向比对到qeury接受（情况6）
        """
        query_direction = "True"
        filtered = ref2_blast.loc[(ref2_blast["qend"] < df["qstart"]),]
        if filtered.empty:
            return pd.Series()
        filtered = filtered.loc[( filtered["ref"].str.contains("right")&(filtered["rstart"]<filtered["rend"]) ) |
                                ( filtered["ref"].str.contains("left")&(filtered["rstart"]>filtered["rend"])), ]
        if filtered.empty:
            return pd.Series()
        ref2_result = filtered.loc[filtered["qend"]==filtered["qend"].max(),].iloc[0]
        ref2_direction = "True" if ref2_result["rstart"] < ref2_result["rend"] else "False"
        return pd.Series({"q_direction": query_direction, "pos1": ref2_result["qend"], "pos2": df["qstart"], "ref2_direction": ref2_direction})
    elif ref1_direction=="left" and df["rstart"] > df["rend"]:
        """
        ref1的左端序列反向比对上query (情况7,8)
        此时首先找query比对区域end右侧比对上ref2的结果
        ref2的左端正向比对到query接受（情况7）
        ref2的右端反向比对到query接受（情况8）
        """
        query_direction = "False"
        filtered = ref2_blast.loc[(ref2_blast["qstart"] > df["qend"]),]
        if filtered.empty:
            return pd.Series()
        filtered = filtered.loc[( filtered["ref"].str.contains("left")&(filtered["rstart"]<filtered["rend"]) ) |
                                ( filtered["ref"].str.contains("right")&(filtered["rstart"]>filtered["rend"]) ),]
        if filtered.empty:
            return pd.Series()
        ref2_result = filtered.loc[filtered["qstart"]==filtered["qstart"].min(),].iloc[0]
        ref2_direction = "True" if ref2_result["rstart"] > ref2_result["rend"] else "False"
        return pd.Series({"q_direction": query_direction, "pos1": ref2_result["qstart"], "pos2": df["qend"], "ref2_direction": ref2_direction})

    print "get_ref2判断出错\n"
    print ref1_direction
    print df["rstart"]
    print df["rend"]
    return pd.Series()

def get_blast_info(ref1_direction, ref1_blast, ref2_blast, query_seq):
    ref2_direction = True  # True为正向 False为反向
    filter_blast = ref1_blast.apply(get_ref2, args=(ref1_direction, ref2_blast,), axis=1)
    filter_blast = filter_blast.dropna(axis=0)
    if filter_blast.empty:
        return False,ref2_direction
    elif len(filter_blast) == 1:
        filter_blast = filter_blast.iloc[0]
    else:
        filter_blast = filter_blast.loc[(filter_blast["pos1"] - filter_blast["pos2"]).abs()==(filter_blast["pos1"] - filter_blast["pos2"]).abs().min(),].iloc[0]
    print "test get_blast_info\n"
    # filter_blast = filter_blast.loc[(filter_blast["pos1"] - filter_blast["pos2"])==(filter_blast["pos1"] - filter_blast["pos2"]).abs().min(),]
    # print filter_blast
    # filter_blast = filter_blast.iloc[0]
    print filter_blast
    if filter_blast["q_direction"] == "True":
        new_query_seq = query_seq[int(filter_blast["pos1"])-1:int(filter_blast["pos2"]):1]
    else:
        new_query_seq = query_seq[int(filter_blast["pos1"])-1:int(filter_blast["pos2"])-2:-1]
        new_query_seq.seq = new_query_seq.seq.complement()  # 增加互补
    return new_query_seq, filter_blast["ref2_direction"]

def find_query(ref, end_ref, seq_list, seq_out_list, blast_data, ref_seq, query_seq):
    """
    避免query两端的来源相同，
    除非对应的是原始的右端序列
    :return:
    """
    have_end_ref=False
    find_query_data = blast_data.loc[blast_data["ref"] == ref, "query"].drop_duplicates()
    if "_right" in ref:
        ref_id = ref[:-6]
        ref1_direction = "right"
    elif "_left" in ref:
        ref_id = ref[:-5]
        ref1_direction = "left"
    end_ref_id = end_ref[:-5]
    sequence = ref_seq[ref_id].seq
    per_query_sequence = sequence
    per_query_seq_out_list = []
    if find_query_data.empty:
        ######return [ref_id], "uncircled", sequence # 直接退出
        return [], "uncircled", sequence
    else:
        for query in find_query_data.tolist():
            find_ref = blast_data.loc[(blast_data["query"]==query) & (blast_data["ref"]!=ref), "ref"].drop_duplicates()
            query_ref1_data = blast_data.loc[(blast_data["query"]==query) & (blast_data["ref"]==ref)]
            if find_ref.empty:
                continue # 没找到query的另一端，直接退出
            else:
                for ref_i in find_ref:
                    query_ref2_data = blast_data.loc[(blast_data["query"]==query) & (blast_data["ref"]==ref_i)]
                    query_cut_seq,ref2_direction = get_blast_info(ref1_direction, query_ref1_data, query_ref2_data, query_seq[query])
                    if not query_cut_seq:
                        continue  # 说明此比对到query另一端的位置不合适，舍弃掉
                    query_cut_seq = query_cut_seq.seq
                    if ref_i == end_ref:
                        if ref_id == end_ref_id:
                            have_end_ref = True
                            sequence_for_one_circled = sequence + query_cut_seq
                        else:
                            if ref1_direction == "right":
                                sequence = sequence + query_cut_seq
                            else:
                                sequence = query_cut_seq + sequence
                            ######return [ref_id], "multi_circled", sequence # 作为被递归调用的程序，发现直接成环，返回所需数据
                            return [], "multi_circled", sequence
                    else:
                        if "_left" in ref_i:
                            ref_i_id = ref_i[:-5]
                            new_ref = ref_i_id + "_right"
                        elif "_right" in ref_i:
                            ref_i_id = ref_i[:-6]
                            new_ref = ref_i_id + "_left"
                        if ref_i_id in seq_out_list:
                            continue # 放弃已有的序列
                        #######
                        if ref_i_id not in seq_list:
                            continue # 序列不在所比较的基因组中，放弃
                        if ref_i_id == ref_id:
                            continue # 不允许非起点序列的左端与右端同时比对一条query
                        else:
                            #####tmp_seq_out_list, status, tmp_sequence = find_query(new_ref, end_ref, seq_out_list + [ref_i_id], blast_data, ref_seq, query_seq)
                            #######添加seq_list
                            tmp_seq_out_list, status, tmp_sequence = find_query(new_ref, end_ref, seq_list, seq_out_list + [ref_id], blast_data, ref_seq, query_seq)
                            tmp_sequence = tmp_sequence if ref2_direction else tmp_sequence.reverse_complement()  # 改反向互补 tmp_sequence[::-1]
                            if status == "multi_circled":
                                if ref1_direction == "right":
                                    sequence += query_cut_seq + tmp_sequence
                                else:
                                    sequence = tmp_sequence + query_cut_seq + sequence
                                return [ref_i_id] + tmp_seq_out_list, status, sequence  # 调用的子程序反馈了成环
                            else:
                                if len(tmp_seq_out_list) >= len(per_query_seq_out_list):
                                    if ref1_direction == "right":
                                        per_query_sequence += query_cut_seq + tmp_sequence
                                    else:
                                        per_query_sequence = tmp_sequence + query_cut_seq + per_query_sequence
                                    per_query_seq_out_list = [ref_i_id] + tmp_seq_out_list # 贪婪算法，永远选择包含参考序列最多的
        if have_end_ref:
            ######return [end_ref_id], "single_circled", sequence_for_one_circled  # 单个染色体成环和多个染色体但不成环，优先选择前者
            return [], "single_circled", sequence_for_one_circled
        else:
            return per_query_seq_out_list, "uncircled", per_query_sequence

def find_gap_multi(seq, seq_list, seq_out_list, blast_data, refseqs, query_seq):
    """
    情况一：一条序列的一端，通过query找到另一条序列的一端 status=multi_circled
    情况二：一条序列的一端，没能通过query找到其他序列 status=uncircled
    情况三：一条序列的一端，通过query只找到了自身的另一端 status=single_circled
    """
    global my_output_fa_record
    start_part = seq + "_right"
    end_part = seq + "_left"
    #######
    tmp_seq_out_list, status, sequence = find_query(start_part, end_part, seq_list, seq_out_list, blast_data, refseqs, query_seq)
    ######
    tmp_seq_out_list.insert(0, seq)
    # seq_out_list 表示之前用过的ref
    # tmp_seq_out_list 表示本次用的那些ref
    new_seq = SeqRecord(sequence, id=seq, description=seq)
    print new_seq
    my_output_fa_record.append(new_seq)
    return {"status": status, "seq_list": tmp_seq_out_list}

def find_gap(genome, seq_list, blast_data, refseqs, query_seq, parent_info):
    global my_output_log_record
    seq_out_list = []
    max_len = 0
    for seq in seq_list:
        seq_record = refseqs[seq]
        if len(seq_record.seq) > max_len:
            max_len = len(seq_record.seq)
            seq_list.insert(0, seq) # 根据长度，将对应id放到首位，后面有重复，不需要处理，根据seq_out_list跳过
    for seq in seq_list:
        if seq in seq_out_list:
            continue
        print seq + "    start"
        print seq_out_list
        ##### result = find_gap_multi(seq, seq_out_list + [seq], blast_data, refseqs, query_seq)
        #######
        result = find_gap_multi(seq, seq_list, seq_out_list, blast_data, refseqs, query_seq)
        if result["status"] == "multi_circled":
            seq_out_list += result["seq_list"]
            suggest_circle = "True"
        elif result["status"] == "uncircled":
            seq_out_list += result["seq_list"]
            suggest_circle = "False"
        elif result["status"] == "single_circled":
            seq_out_list += result["seq_list"]
            suggest_circle = "True"
        # 处理log
        print seq + "   end"
        print seq_out_list
        parent_list = []
        for ids in result["seq_list"]:
            parent_str = parent_info[(genome, ids)]["parent"]
            parent_list += parent_str.split(",")
        my_output_log_record = my_output_log_record.append({"loc": genome, "seq": seq, "from": ",".join(parent_list), "suggest_circle": suggest_circle}, ignore_index=True)


def parse(fasta, query, blast_table, log_file, output_prefix, aln_len=1000, identity=95):
    global my_output_log_record, my_output_fa_record
    seq_records = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    query_records = SeqIO.to_dict(SeqIO.parse(query, "fasta"))
    blast_data = parse_blast(blast_table, aln_len, identity)
    parent_info, genome_info, have_circle_info = parse_log(log_file)
    for genome in genome_info:
        seq_list = genome_info[genome]
        if len(seq_list) == 1 and have_circle_info:
            if parent_info[(genome, seq_list[0])]["is_circle"]:
                my_output_log_record = my_output_log_record.append({"loc": genome, "seq": seq_list[0], "from": parent_info[(genome, seq_list[0])]["parent"], "suggest_circle": "True"}, ignore_index=True)
                my_output_fa_record.append(seq_records[seq_list[0]])
            else:
                find_gap(genome, seq_list, blast_data, seq_records, query_records, parent_info)  # 为了将它拼成环
        elif len(seq_list) == 1 and not have_circle_info:
            find_gap(genome, seq_list, blast_data, seq_records, query_records, parent_info)  # 没有是否成环的信息，尝试拼成环
        else:
            find_gap(genome, seq_list, blast_data, seq_records, query_records, parent_info)


    # 写log记录
    my_output_log_record.reindex(columns=["loc", "seq", "from", "suggest_circle"]).to_csv(output_prefix + ".log", sep="\t", header=False, index=False)
    # 生成fasta
    SeqIO.write(my_output_fa_record, output_prefix + ".fa", "fasta")

def _main():
    parser = argparse.ArgumentParser(description='filter fa by json file')
    parser.add_argument('-fa', '--fa', help="fasta file")
    parser.add_argument('-query', '--query', help="query fasta file")
    parser.add_argument('-blast', '--blast', help="blast table in format 6")
    parser.add_argument('-log', '--log', help="log table")
    parser.add_argument('-aln_len', '--aln_len', help="aln_len", type=int)
    parser.add_argument('-identity', '--identity', help="identity", type=float)
    parser.add_argument('-o_prefix', '--o_prefix', help="output prefix, output a new fasta file and a log file")
    args = parser.parse_args()
    parse(args.fa, args.query, args.blast, args.log, args.o_prefix, args.aln_len, args.identity)


if __name__ == "__main__":
    _main()
    '''
        python get_fa_by_json.py -fa ~/app/bioinfo/bin/download/test.test.fa -json "{\"chr1\": [1,2,3,5,9,55], \"plasmidB\": [4,6,8]}" -seq_prefix Scaffold -o_prefix test
    '''