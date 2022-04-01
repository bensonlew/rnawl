# -*- coding: utf-8 -*-
# __author__ = 'xuting'

"""移动可变碱基, 滤去找不到barcode的序列和嵌合barcode的序列"""
from __future__ import division
import os
import errno
import multiprocessing
import time
from collections import defaultdict
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from mbio.packages.datasplit.miseq_split import code2index, str_check


class RawValidAgent(Agent):
    def __init__(self, parent=None):
        super(RawValidAgent, self).__init__(parent)
        self._run_mode = "ssh1"
        options = [
            {'name': 'sample_info', 'type': "infile", 'format': 'datasplit.miseq_split'},  # 样本拆分信息表
            {'name': 'unzip_path', 'type': "string"}  # bcl2fastq软件拆分出来的fastq解压后的输出目录
        ]
        self.add_option(options)

    def check_option(self):
        if not self.option('sample_info').is_set:
            raise OptionError("参数sample_info不能为空")
        if not self.option('unzip_path'):
            raise OptionError("参数unzip_path不能为空")
        return True

    def set_resource(self):
        self._cpu = 32
        self._memory = ''


class RawValidTool(Tool):
    def __init__(self, config):
        super(RawValidTool, self).__init__(config)
        self._version = 1.0
        self.option('sample_info').get_info()
        self.total = defaultdict(int)
        self.noFindex = defaultdict(int)
        self.noRindex = defaultdict(int)
        self.chimeric = defaultdict(int)
        self.Findex = dict()
        self.Rindex = dict()
        self.Fvarbase = dict()
        self.Rvarbase = dict()
        self.vilid = defaultdict(int)
        self.dump_database_to_memery()
        self.l = multiprocessing.Lock()

    def dump_database_to_memery(self):
        """
        一直访问数据库读取数据会导致运行变慢
        所以先把里面的内容读取到内存中以后再开始运行程序
        """
        for c_id in self.option('sample_info').prop["child_ids"]:
            index_code = self.option('sample_info').child_sample(c_id, "index")
            (self.Findex[c_id], self.Rindex[c_id], self.Fvarbase[c_id], self.Rvarbase[c_id]) = code2index(index_code)

    def make_ess_dir(self):
        var2endDir = os.path.join(self.work_dir, "var2end")
        dir_list = [var2endDir]
        for name in dir_list:
            try:
                os.makedirs(name)
            except OSError as exc:
                if exc.errno == errno.EEXIST and os.path.isdir(name):
                    pass
                else:
                    raise OSError("创建目录失败")

    def WriteStat(self, p):
        p_id = p["sample_id"]
        stat_path = os.path.join(self.work_dir, "output", "stat.xls")
        self.logger.debug(self.noFindex)
        self.logger.debug(self.noRindex)
        self.logger.debug(self.vilid)
        self.logger.debug(self.chimeric)
        self.logger.debug(self.total)
        with open(stat_path, 'ab') as a:
            valid_rate = "{:.2f}%".format((self.vilid[p_id] / self.total[p_id]) * 100)
            str_ = (p['sample_id'] + "\t" + p["library_name"] + "\t" + str(self.total[p_id])
                    + "\t" + str(self.vilid[p_id]) + "\t" + valid_rate + "\t" + str(self.chimeric[p_id])
                    + "\t" + str(self.noRindex[p_id]) + "\t" + str(self.noFindex[p_id]) + "\n")
            a.write(str_)

    def MoveVarBase(self):
        stat_path = os.path.join(self.work_dir, "output", "stat.xls")
        with open(stat_path, 'wb') as w:
            w.write("#library_id\tlibrary_name\ttotal_reads\tvalid_reads\tvalid_rate\tchimeric_reads\tno_Rbarcode\tno_Fbarcode\n")
        process_list = list()
        for p in self.option('sample_info').prop["parent_sample"]:
            if p["has_child"]:
                process = multiprocessing.Process(target=self._MoveVarBase, args=(p,))
                process_list.append(process)
        for my_p in process_list:
            my_p.daemon = True
            my_p.start()
            time.sleep(2)
        for my_p in process_list:
            my_p.join()
        self.logger.info("Library with child done!")

        with open(stat_path, 'ab') as a:
            for p in self.option('sample_info').prop["parent_sample"]:
                if not p["has_child"]:
                    file_r1 = os.path.join(self.option('unzip_path'), p['sample_id'] + "_r1.fastq")
                    num_lines = sum(1 for line in open(file_r1))
                    total_reads = int(num_lines / 4)
                    str_ = (p['sample_id'] + "\t" + p["library_name"] + "\t" + str(total_reads)
                            + "\t" + str(total_reads) + "\t" + "NA\tNA\tNA\tNA\n"
                            )
                    a.write(str_)

    def _MoveVarBase(self, p):
        """
        对每一个文库进行处理，方便多进程
        """
        var2endDir = os.path.join(self.work_dir, "var2end")
        p_id = p['sample_id']
        file_r1 = os.path.join(self.option('unzip_path'), p['sample_id'] + "_r1.fastq")
        file_r2 = os.path.join(self.option('unzip_path'), p['sample_id'] + "_r2.fastq")
        file3 = os.path.join(var2endDir, p['sample_id'] + ".valid_r1.fastq")
        file4 = os.path.join(var2endDir, p['sample_id'] + ".valid_r2.fastq")
        file5 = os.path.join(var2endDir, p['sample_id'] + ".discard_r1.fastq")
        file6 = os.path.join(var2endDir, p['sample_id'] + ".discard_r2.fastq")
        with open(file_r1, 'rb') as r1, open(file_r2, 'rb') as r2, open(file3, 'wb') as w1, open(file4, 'wb') as w2, open(file5, 'wb') as w3, open(file6, 'wb') as w4:
            c_id_list = list()
            c_id_list = self.option('sample_info').find_child_ids(p_id)
            for head1 in r1:
                self.total[p_id] += 1
                if self.total[p_id] % 10000 == 0:
                    self.l.acquire()
                    self.logger.info("Processing library {}: {}".format(p_id, self.total[p_id]))
                    self.l.release()
                head1 = head1.rstrip("\r\n")
                seq1 = r1.next().rstrip("\r\n")
                desc1 = r1.next().rstrip("\r\n")
                qual1 = r1.next().rstrip("\r\n")
                head2 = r2.next().rstrip("\r\n")
                seq2 = r2.next().rstrip("\r\n")
                desc2 = r2.next().rstrip("\r\n")
                qual2 = r2.next().rstrip("\r\n")
                flag = 0
                for c_id in c_id_list:
                    my_Flen = self.Fvarbase[c_id] + len(self.Findex[c_id])
                    my_Rlen = self.Rvarbase[c_id] + len(self.Rindex[c_id])
                    seq1RealIndex = seq1[self.Fvarbase[c_id]:my_Flen]
                    seq2RealIndex = seq2[self.Rvarbase[c_id]:my_Rlen]
                    fIndexMiss = str_check(seq1RealIndex, self.Findex[c_id])

                    if fIndexMiss == 0:
                        flag = 1
                        rIndexMiss = str_check(seq2RealIndex, self.Rindex[c_id])
                        if rIndexMiss == 0:  # 正确匹配
                            self.vilid[p_id] += 1
                            newSeq1 = seq1[self.Fvarbase[c_id]:] + seq1[0:self.Fvarbase[c_id]]
                            newSeq2 = seq2[self.Rvarbase[c_id]:] + seq2[0:self.Rvarbase[c_id]]
                            newQual1 = qual1[self.Fvarbase[c_id]:] + "#" * self.Fvarbase[c_id]
                            newQual2 = qual2[self.Rvarbase[c_id]:] + "#" * self.Rvarbase[c_id]
                            w1.write(head1 + "\n" + newSeq1 + "\n" + desc1 + "\n" + newQual1 + "\n")
                            w2.write(head2 + "\n" + newSeq2 + "\n" + desc2 + "\n" + newQual2 + "\n")
                            flag = 2
                        else:
                            for tmp_c_id in c_id_list:
                                if tmp_c_id != c_id:
                                    my_Rlen = self.Rvarbase[tmp_c_id] + len(self.Rindex[tmp_c_id])
                                    seq2RealIndex = seq2[self.Rvarbase[tmp_c_id]:my_Rlen]
                                    rIndexMiss = str_check(seq2RealIndex, self.Rindex[tmp_c_id])
                                    if rIndexMiss == 0:  # 嵌合
                                        self.chimeric[p_id] += 1
                                        w3.write(head1 + "\tFbarcode=" + self.Findex[c_id] + "\tRbarcode=" +
                                                 self.Rindex[tmp_c_id] + "\n" + seq1 + "\n" + desc1 + "\n" +
                                                 qual1 + "\n")
                                        w4.write(head2 + "\tFbarcode=" + self.Findex[c_id] + "\tRbarcode=" +
                                                 self.Rindex[tmp_c_id] + "\n" + seq2 + "\n" + desc2 + "\n" +
                                                 qual2 + "\n")
                                        flag = 2
                                        break
                    else:  # 尝试用r2去匹配左barcode
                        seq1RealIndex = seq1[self.Rvarbase[c_id]:my_Rlen]
                        seq2RealIndex = seq2[self.Fvarbase[c_id]:my_Flen]
                        fIndexMiss = str_check(seq2RealIndex, self.Findex[c_id])
                        if fIndexMiss == 0:
                            flag = 1
                            rIndexMiss = str_check(seq1RealIndex, self.Rindex[c_id])
                            if rIndexMiss == 0:
                                self.vilid[p_id] += 1
                                newSeq1 = seq1[self.Rvarbase[c_id]:] + seq1[0:self.Rvarbase[c_id]]
                                newSeq2 = seq2[self.Fvarbase[c_id]:] + seq2[0:self.Fvarbase[c_id]]
                                newQual1 = qual1[self.Rvarbase[c_id]:] + "#" * self.Rvarbase[c_id]
                                newQual2 = qual2[self.Fvarbase[c_id]:] + "#" * self.Fvarbase[c_id]
                                w1.write(head1 + "\n" + newSeq1 + "\n" + desc1 + "\n" + newQual1 + "\n")
                                w2.write(head2 + "\n" + newSeq2 + "\n" + desc2 + "\n" + newQual2 + "\n")
                                flag = 2
                                break
                            else:
                                for tmp_c_id in c_id_list:
                                    if tmp_c_id != c_id:
                                        my_Rlen = self.Rvarbase[tmp_c_id] + len(self.Rindex[tmp_c_id])
                                        seq1RealIndex = seq1[self.Rvarbase[tmp_c_id]:my_Rlen]
                                        rIndexMiss = str_check(seq1RealIndex, self.Rindex[tmp_c_id])
                                        if rIndexMiss == 0:
                                            self.chimeric[p_id] += 1
                                            w3.write(head1 + "\tFbarcode=" + self.Findex[c_id] + "\tRbarcode=" +
                                                     self.Rindex[tmp_c_id] + "\n" + seq1 + "\n" + desc1 + "\n" +
                                                     qual1 + "\n")
                                            w4.write(head2 + "\tFbarcode=" + self.Findex[c_id] + "\tRbarcode=" +
                                                     self.Rindex[tmp_c_id] + "\n" + seq2 + "\n" + desc2 + "\n" +
                                                     qual2 + "\n")
                                            flag = 2
                                            break
                if flag == 1:
                    self.noRindex[p_id] += 1
                    w3.write(head1 + "\tnoRbarcode" + "\n" + seq1 + "\n" + desc1 + "\n" + qual1 + "\n")
                    w4.write(head2 + "\tnoRbarcode" + "\n" + seq2 + "\n" + desc2 + "\n" + qual2 + "\n")
                elif flag == 0:
                    self.noFindex[p_id] += 1
                    w3.write(head1 + "\tnoFbarcode" + "\n" + seq1 + "\n" + desc1 + "\n" + qual1 + "\n")
                    w4.write(head2 + "\tnoFbarcode" + "\n" + seq2 + "\n" + desc2 + "\n" + qual2 + "\n")
        self.WriteStat(p)

    def run(self):
        super(RawValidTool, self).run()
        self.make_ess_dir()
        self.MoveVarBase()
        self.end()
