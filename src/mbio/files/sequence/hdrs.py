# -*- coding: utf-8 -*-

# __author__ = 'linfang.jin'
# time: 2017/1/13 17:26

import re
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from mbio.files.sequence.fasta import FastaFile
import os
from Bio import SeqIO


class HdrsFile(File):
    '''
    定义与fa文件配套的Hdrs文件
    author:linfang.jin
    date: 2017.01.15
    '''

    def __init__(self):
        super(HdrsFile, self).__init__()

    def check(self):
        """
         检测文件是否满足要求,发生错误时应该触发FileError异常
        :return:
        """
        if super(HdrsFile, self).check():
            if not self.check_line():
                raise FileError("文件格式错误")
            # if not self.check_consistence():
            #     raise FileError("hdrs文件（{}）和他的fasta文件（{}）内容上不一致".format(self.path, self.fasta.path))
        return True

    def set_fasta(self, fasta_obj):
        if not isinstance(fasta_obj, FastaFile):
            raise FileError("传入的文件不是FastaFile对象")
        self._properties['co_fasta_file'] = fasta_obj
        return self

    def fasta(self):
        if self._properties.has_key("co_fasta_file"):
            return self._properties['co_fasta_file']
        else:
            raise Exception("这个fasta文件还没有被设置过与其对应的hdrs文件")

    def check_consistence(self):
        records = self.get_line_dic()
        fa_seqs = self.get_fa_stat()
        if records == fa_seqs:
            return True
        else:
            return False


    def get_fa_stat(self):
        fa_file_path = self.fasta().path
        records = set()
        org_name = os.path.basename(fa_file_path).strip().split(".")[0]
        for seq in SeqIO.parse(fa_file_path,"fasta"):
            seq_id = seq.id
            seq_len = len(seq.seq)
            nonNlen = seq_len - seq.seq.count("N") - seq.seq.count("n")
            d = {"id": seq_id, "seq_len": seq_len, "seq_non_n_len":nonNlen, "org": org_name}
            records.add(d)
        return records

    def get_line_dic(self):
        if super(HdrsFile, self).check():
            try:
                record_set = set()
                with open(self.path) as fr:
                    line = fr.readline().strip()
                    m = re.match(r">(\S+)\s+/len=(\d+)\s+/nonNlen=(\d+)\s+/org=(\S+)", line)
                    if m:
                        d = {"id": m.group(1), "seq_len": m.group(2), "seq_non_n_len": m.group(3), "org": m.group(4)}
                        record_set.add(d)
            except Exception as e:
                raise Exception("打开hdrs文件运行出错: {}".format(e))
        return record_set

    def check_line(self):
        """
        检查每行是否满足hdrs的语法要求
        :return:
        """
        try:
            with open(self.path) as fr:
                line = fr.readline().strip()
                m = re.match(r">(\S+)\s+/len=(\d+)\s+/nonNlen=(\d+)\s+/org=(\S+)", line)
                if not m:
                    raise FileError("hdrs文件（{}）的这一行：{}不满足要求".format(self.path, line))
        except Exception as e:
            raise Exception("打开hdrs文件运行出错: {}".format(e))
        return True

    def write_hdrs(self):
        '''
        生成fasta文件的.hdrs文件: 每一行的格式为">seq_id len=seq_length /nonNlen=non_N_seq_length  /org=org_name"
        供asprofile使用
        :param hdrs_file_path:
        :return: hdrs_file_path
        author: linfang.jin
        date:2017.01.15
        '''

        try:
            fw = open(self.path, 'w')
            org_name = os.path.basename(self.fasta.path).strip().split(".")[0]
            for seq in SeqIO.parse(self.fasta.path, "fasta"):
                seq_id = seq.id
                seq_len = len(seq.seq)
                nonNlen = seq_len - seq.seq.count("N") - seq.seq.count("n")
                newline = ">{} len={} /nonNlen={}  /org={}\n".format(str(seq_id), str(seq_len), str(nonNlen), org_name)
                fw.write(newline)
            fw.close()
        except Exception as e:
            raise FileError("{}:生成{}的hdrs文件({})过程中出错".format(e, self.fasta.path, self.path))
