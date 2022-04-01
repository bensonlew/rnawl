# -*- coding: utf-8 -*-
# __author__ = 'xuting'

"""对经过trimmomatic的序列进行截短， 方便后面的merge"""
import os
import errno
import multiprocessing
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class ChompSeqAgent(Agent):
    def __init__(self, parent=None):
        super(ChompSeqAgent, self).__init__(parent)
        self._run_mode = "ssh1"
        options = [
            {'name': 'sample_info', 'type': "infile", 'format': 'datasplit.miseq_split'},  # 样本拆分信息表
            {'name': 'trim_path', 'type': "string"}  # 经过质控之后序列的目录
        ]
        self.add_option(options)

    def check_option(self):
        if not self.option('sample_info').is_set:
            raise OptionError("参数sample_info不能为空")
        if not self.option('trim_path'):
            raise OptionError("参数trim_path不能为空")
        return True

    def set_resource(self):
        self._cpu = 32
        self._memory = ''


class ChompSeqTool(Tool):
    def __init__(self, config):
        super(ChompSeqTool, self).__init__(config)
        self._version = 1.0
        self.option('sample_info').get_info()
        self.l = multiprocessing.Lock()

    def make_ess_dir(self):
        chompDir = os.path.join(self.work_dir, "chomped")
        dir_list = [chompDir]
        for name in dir_list:
            try:
                os.makedirs(name)
            except OSError as exc:
                if exc.errno == errno.EEXIST and os.path.isdir(name):
                    pass
                else:
                    raise OSError("创建目录失败")

    def ChompSeq(self):
        process_list = list()
        for p in self.option('sample_info').prop["parent_sample"]:
            if p["has_child"]:
                process = multiprocessing.Process(target=self._ChompSeq, args=(p,))
                process_list.append(process)
        for my_p in process_list:
            my_p.daemon = True
            my_p.start()
        for my_p in process_list:
            my_p.join()

    def _ChompSeq(self, p):
        p_id = p['sample_id']
        chompDir = os.path.join(self.work_dir, "chomped")
        file_r1 = os.path.join(self.option("trim_path"), p_id + ".trimmo_r1.fastq")
        file_r2 = os.path.join(self.option("trim_path"), p_id + ".trimmo_r2.fastq")
        file_chomp1 = os.path.join(chompDir, p_id + ".chomped_r1.fastq")
        file_chomp2 = os.path.join(chompDir, p_id + ".chomped_r2.fastq")
        chomp_len = self.option("sample_info").parent_sample(p_id, "chomp_len")
        try:
            chomp_len = int(chomp_len)
        except:
            chomp_len = "default"
        if chomp_len == "default":
            insert_len = int(p["insert_len"])
            if insert_len <= 220:
                chomp_len = 160
            elif insert_len <= 320:
                chomp_len = 200
            elif insert_len <= 380:
                chomp_len = 250
            else:
                chomp_len = 99999
        with open(file_r1, 'rb') as r1, open(file_r2, 'rb') as r2, open(file_chomp1, 'wb') as w1, open(file_chomp2, 'wb') as w2:
            for head1 in r1:
                head1 = head1.rstrip("\r\n")
                seq1 = r1.next().rstrip("\r\n")
                desc1 = r1.next().rstrip("\r\n")
                qual1 = r1.next().rstrip("\r\n")
                head2 = r2.next().rstrip("\r\n")
                seq2 = r2.next().rstrip("\r\n")
                desc2 = r2.next().rstrip("\r\n")
                qual2 = r2.next().rstrip("\r\n")

                newseq1 = seq1[0:chomp_len]
                newQual1 = qual1[0:chomp_len]
                newseq2 = seq2[0:chomp_len]
                newQual2 = qual2[0:chomp_len]
                str1 = "{}\n{}\n{}\n{}\n".format(head1, newseq1, desc1, newQual1)
                str2 = "{}\n{}\n{}\n{}\n".format(head2, newseq2, desc2, newQual2)
                self.l.acquire()
                w1.write(str1)
                w2.write(str2)
                self.l.release()

    def run(self):
        super(ChompSeqTool, self).run()
        self.make_ess_dir()
        self.ChompSeq()
        self.end()
