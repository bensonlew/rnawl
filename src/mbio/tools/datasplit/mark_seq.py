# -*- coding: utf-8 -*-
# __author__ = 'xuting'

"""配对左右barcode和引物, 正确配对之后切去, 然后将序列名标记上样本的id"""
from __future__ import division
import os
import errno
import re
import multiprocessing
import time
from collections import defaultdict
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from mbio.packages.datasplit.miseq_split import reverse_complement, code2index, code2primer, str_check


class MarkSeqAgent(Agent):
    def __init__(self, parent=None):
        super(MarkSeqAgent, self).__init__(parent)
        self._run_mode = "ssh1"
        options = [
            {'name': 'sample_info', 'type': "infile", 'format': 'datasplit.miseq_split'},  # 样本拆分信息表
            {'name': 'merge_path', 'type': "string"}  # 经过merge之后序列的目录
        ]
        self.add_option(options)

    def check_option(self):
        if not self.option('sample_info').is_set:
            raise OptionError("参数sample_info不能为空")
        if not self.option('merge_path'):
            raise OptionError("参数merge_path不能为空")
        return True

    def set_resource(self):
        self._cpu = 32
        self._memory = ''


class MarkSeqTool(Tool):
    def __init__(self, config):
        super(MarkSeqTool, self).__init__(config)
        self._version = 1.0
        self.option('sample_info').get_info()
        self.Findex = dict()
        self.Rindex = dict()
        self.Fvarbase = dict()
        self.Rvarbase = dict()
        self.fPrimer = dict()
        self.rPrimer = dict()
        self.leftChompLen = dict()
        self.rightChompLen = dict()
        self.primerMissMatch = defaultdict(int)
        self.valid = defaultdict(int)
        self.RbarcodeMissMatch = defaultdict(int)
        self.noBarcode = defaultdict(int)
        self.total = defaultdict(int)
        self.dump_database_to_memery()
        self.l = multiprocessing.Lock()

    def dump_database_to_memery(self):
        for c_id in self.option('sample_info').prop["child_ids"]:
            index_code = self.option('sample_info').child_sample(c_id, "index")
            (self.Findex[c_id], self.Rindex[c_id], self.Fvarbase[c_id], self.Rvarbase[c_id]) = code2index(index_code)
            primercode = self.option('sample_info').child_sample(c_id, "primer")
            (self.fPrimer[c_id], self.rPrimer[c_id]) = code2primer(primercode)
            self.leftChompLen[c_id] = len(self.Findex[c_id]) + len(self.fPrimer[c_id])
            self.rightChompLen[c_id] = len(self.Rindex[c_id])

    def make_ess_dir(self):
        marked_dir = os.path.join(self.work_dir, "markedSeq")
        dir_list = [marked_dir]
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
        self.logger.debug(self.total)
        self.logger.debug(self.valid)
        self.logger.debug(self.primerMissMatch)
        self.logger.debug(self.RbarcodeMissMatch)
        self.logger.debug(self.noBarcode)
        with open(stat_path, 'ab') as a:
            str_ = (p['sample_id'] + "\t" + p["library_name"] + "\t" + str(self.total[p_id]) + "\t"
                    + str(self.valid[p_id]) + "\t" + str(self.primerMissMatch[p_id]) + "\t"
                    + str(self.RbarcodeMissMatch[p_id]) + "\t" + str(self.noBarcode[p_id]) + "\n")
            a.write(str_)
        self.logger.info(p_id + "end!")

    def MarkSeq(self):
        stat_path = os.path.join(self.work_dir, "output", "stat.xls")
        with open(stat_path, 'wb') as w:
            w.write("#library_id\tlibrary_name\ttotal\tvilid\tprimerMissMatch\tRbarcodeMissMatch\tnoBarcode\n")
        process_list = list()
        for p in self.option('sample_info').prop["parent_sample"]:
            if p["has_child"]:
                process = multiprocessing.Process(target=self._MarkSeq, args=(p,))
                process_list.append(process)
        for my_p in process_list:
            my_p.daemon = True
            my_p.start()
            time.sleep(2)
        for my_p in process_list:
            my_p.join()

    def _MarkSeq(self, p):
        """
        对每一个文库进行处理，方便多进程
        """
        p_id = p['sample_id']
        file_ = os.path.join(self.option('merge_path'), p_id + '.merged.fastq')
        fileMarked = os.path.join(self.work_dir, "markedSeq", p_id + '.marked.fastq')
        fileDiscarded = os.path.join(self.work_dir, "markedSeq", p_id + '.discarded.fastq')
        self.l.acquire()
        self.logger.info("processing Library: {}".format(p_id))
        self.l.release()
        with open(file_, 'rb') as r, open(fileMarked, 'wb') as w1, open(fileDiscarded, 'wb') as w2:
            c_id_list = list()
            c_id_list = self.option('sample_info').find_child_ids(p_id)
            for head in r:
                self.total[p_id] += 1
                if (self.total[p_id] % 10000) == 0:
                    self.l.acquire()
                    self.logger.info("Processing Library {}: {}".format(p_id, self.total[p_id]))
                    self.l.release()
                head = head.rstrip('\r\n')
                head = re.sub('^@', '', head)
                ori_seq = r.next().rstrip('\r\n')
                rev_ori_seq = ori_seq[::-1]
                rev_ori_seq = reverse_complement(rev_ori_seq)
                desc = r.next().rstrip('\r\n')
                qual = r.next().rstrip('\r\n')
                rev_qual = qual[::-1]
                flag = 0
                for c_id in c_id_list:
                    seq1RealIndex = ori_seq[0:len(self.Findex[c_id])]
                    seq2RealIndex = rev_ori_seq[0:len(self.Rindex[c_id])]
                    fIndexMiss = str_check(seq1RealIndex, self.Findex[c_id])
                    if fIndexMiss == 0:
                        rIndexMiss = str_check(seq2RealIndex, self.Rindex[c_id])
                        if rIndexMiss == 0:
                            tmp = len(self.Findex[c_id]) + len(self.fPrimer[c_id])
                            seqPrimer = ori_seq[len(self.Findex[c_id]):tmp]
                            PrimerMiss = str_check(seqPrimer, self.fPrimer[c_id])
                            if PrimerMiss <= int(self.option('sample_info').child_sample(c_id, "primer_miss")):
                                sp_name = self.option('sample_info').child_sample(c_id, "sample_name")
                                new_head = "@{}_{}\t{}\torig_bc={}\tnew_bc={}\tbc_diffs=0".format(c_id, sp_name, head, self.Findex[c_id], self.Findex[c_id])
                                new_seq = ori_seq[self.leftChompLen[c_id]:-self.rightChompLen[c_id]]
                                new_qual = qual[self.leftChompLen[c_id]:-self.rightChompLen[c_id]]
                                str_ = "{}\n{}\n{}\n{}\n".format(new_head, new_seq, desc, new_qual)
                                w1.write(str_)
                                flag = 1
                                self.valid[p_id] += 1
                                break
                            else:
                                self.primerMissMatch[p_id] += 1
                                new_head = "@{}\tPrimer".format(head)
                                str_ = "{}\n{}\n{}\n{}\n".format(new_head, ori_seq, desc, qual)
                                w2.write(str_)
                                flag = 1
                                break
                        else:
                            self.RbarcodeMissMatch[p_id] += 1
                            new_head = "@{}\tRBarcode".format(head)
                            str_ = "{}\n{}\n{}\n{}\n".format(new_head, ori_seq, desc, qual)
                            w2.write(str_)
                            flag = 1
                            break
                    else:
                        seq1RealIndex = rev_ori_seq[0:len(self.Findex[c_id])]
                        seq2RealIndex = ori_seq[0:len(self.Rindex[c_id])]
                        fIndexMiss = str_check(seq1RealIndex, self.Findex[c_id])
                        if fIndexMiss == 0:
                            rIndexMiss = str_check(seq2RealIndex, self.Rindex[c_id])
                            if rIndexMiss == 0:
                                tmp = len(self.Findex[c_id]) + len(self.fPrimer[c_id])
                                seqPrimer = rev_ori_seq[len(self.Findex[c_id]):tmp]
                                if PrimerMiss <= int(self.option('sample_info').child_sample(c_id, "primer_miss")):
                                    sp_name = self.option('sample_info').child_sample(c_id, "sample_name")
                                    new_head = "@{}_{}\t{}\torig_bc={}\tnew_bc={}\tbc_diffs=0".format(c_id, sp_name, head, self.Findex[c_id], self.Findex[c_id])
                                    new_seq = rev_ori_seq[self.leftChompLen[c_id]:-self.rightChompLen[c_id]]
                                    new_qual = rev_qual[self.leftChompLen[c_id]:-self.rightChompLen[c_id]]
                                    str_ = "{}\n{}\n{}\n{}\n".format(new_head, new_seq, desc, new_qual)
                                    w1.write(str_)
                                    flag = 1
                                    self.valid[p_id] += 1
                                    break
                                else:
                                    self.primerMissMatch[p_id] += 1
                                    new_head = "@{}\tPrimer".format(head)
                                    str_ = "{}\n{}\n{}\n{}\n".format(new_head, ori_seq, desc, qual)
                                    w2.write(str_)
                                    flag = 1
                                    break
                            else:
                                self.RbarcodeMissMatch[p_id] += 1
                                new_head = "@{}\tPrimer".format(head)
                                str_ = "{}\n{}\n{}\n{}\n".format(new_head, ori_seq, desc, qual)
                                w2.write(str_)
                                flag = 1
                                break
                if flag == 0:
                    self.noBarcode[p_id] += 1
                    new_head = "@{}\tnoBarcode".format(head)
                    str_ = "{}\n{}\n{}\n{}\n".format(new_head, ori_seq, desc, qual)
                    w2.write(str_)
        self.WriteStat(p)

    def run(self):
        super(MarkSeqTool, self).run()
        self.make_ess_dir()
        self.MarkSeq()
        self.logger.info("end")
        self.end()
