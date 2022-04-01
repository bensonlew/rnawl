# -*- coding: utf-8 -*-
# __author__ = 'xuting'

"""将标记好的文库拆分到每一个样本，并从qual_control(trimmomatic)的结果中提取相应的raw序列"""
import errno
import multiprocessing
import os
import re
import time
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class SeqExtractAgent(Agent):
    def __init__(self, parent=None):
        super(SeqExtractAgent, self).__init__(parent)
        self._run_mode = "ssh1"
        options = [
            {'name': 'sample_info', 'type': "infile", 'format': 'datasplit.miseq_split'},  # 样本拆分信息表
            {'name': 'LenControled_path', 'type': "string"},  # 经过长度控制之后序列的文库
            {'name': 'rawValid_path', 'type': "string"}  # 经过长度控制之后序列的文库
        ]
        self.add_option(options)

    def check_option(self):
        if not self.option('sample_info').is_set:
            raise OptionError("参数sample_info不能为空")
        if not self.option('LenControled_path'):
            raise OptionError("参数LenControled_path不能为空")
        return True

    def set_resource(self):
        self._cpu = 32
        self._memory = ''


class SeqExtractTool(Tool):
    def __init__(self, config):
        super(SeqExtractTool, self).__init__(config)
        self._version = 1.0
        self.option('sample_info').get_info()
        self.id2Name = dict()
        self.GetId2Name()
        self.l = multiprocessing.Lock()

    def GetId2Name(self):
        for c_id in self.option('sample_info').prop["child_ids"]:
            self.id2Name[c_id] = self.option('sample_info').child_sample(c_id, "sample_name")

    def make_ess_dir(self):
        dir_list = list()
        for p in self.option('sample_info').prop["parent_sample"]:
            if p["has_child"]:
                p_id = p["sample_id"]
                libraryDir = os.path.join(self.work_dir, 'child_sample', "library_" + p_id, "final")
                dir_list.append(libraryDir)
                libraryDir = os.path.join(self.work_dir, 'child_sample', "library_" + p_id, "rawValid")
                dir_list.append(libraryDir)
        for name in dir_list:
            try:
                os.makedirs(name)
            except OSError as exc:
                if exc.errno == errno.EEXIST and os.path.isdir(name):
                    pass
                else:
                    raise OSError("创建目录失败")

    def SeqExtract(self):
        process_list = list()
        for p in self.option('sample_info').prop["parent_sample"]:
            if p["has_child"]:
                process = multiprocessing.Process(target=self._SeqExtract, args=(p,))
                process_list.append(process)
        for my_p in process_list:
            my_p.daemon = True
            my_p.start()
            time.sleep(2)
        for my_p in process_list:
            my_p.join()

    def _SeqExtract(self, p):
        p_id = p['sample_id']
        libraryDir = os.path.join(self.work_dir, 'child_sample', "library_" + p_id)
        finalDir = os.path.join(libraryDir, "final")
        rawValidDir = os.path.join(libraryDir, "rawValid")
        file_ = os.path.join(self.option("LenControled_path"), p_id + ".lenControl.fastq")
        rawfile1 = os.path.join(self.option("rawValid_path"), p_id + ".valid_r1.fastq")
        rawfile2 = os.path.join(self.option("rawValid_path"), p_id + ".valid_r2.fastq")
        seqName2sampleId = dict()
        finalHandle = dict()
        rawValidHandle1 = dict()
        rawValidHandle2 = dict()
        count = 0
        for c_id in self.option('sample_info').find_child_ids(p_id):
            finalName = os.path.join(finalDir, c_id + ".final.fastq")
            rawValidName1 = os.path.join(rawValidDir, c_id + ".rawValid_r1.fastq")
            rawValidName2 = os.path.join(rawValidDir, c_id + ".rawValid_r2.fastq")
            finalHandle["finalHandle" + c_id] = open(finalName, 'wb')
            rawValidHandle1["rawValidHandle1_" + c_id] = open(rawValidName1, 'wb')
            rawValidHandle2["rawValidHandle2_" + c_id] = open(rawValidName2, 'wb')
        with open(file_, 'rb') as r:
            for head in r:
                count += 1
                if count % 50000 == 0:
                    self.l.acquire()
                    self.logger.info("Library:{}, processing sequence {}".format(p_id, count))
                    self.l.release()
                head = head.rstrip("\r\n")
                head = re.sub("@", "", head)
                tmp_head = re.split("\t", head)
                tmp_head.pop(0)
                tmp_head = "\t".join(tmp_head)
                seq = r.next()
                desc = r.next()
                qual = r.next()
                seqName = re.split("\t", head)[1]
                sampleId = re.split("\t", head)[0]
                sampleId = re.split("_", sampleId)[0]
                seqName2sampleId[seqName] = sampleId
                newHead = "@{}_{}\t{}\n".format(self.id2Name[sampleId], count, tmp_head)
                finalHandle["finalHandle" + sampleId].write(newHead + seq + desc + qual)
        for myHandle in finalHandle.values():
            myHandle.close()
        with open(rawfile1, 'rb') as r1, open(rawfile2, 'rb') as r2:
            for head1 in r1:
                head1 = head1.rstrip("\r\n")
                head1 = re.sub("@", "", head1)
                seq1 = r1.next()
                desc1 = r1.next()
                qual1 = r1.next()
                head2 = r2.next()
                seq2 = r2.next()
                desc2 = r2.next()
                qual2 = r2.next()
                seqName = re.split("\t", head1)[0]
                if seqName in seqName2sampleId:
                    my_sampleId = seqName2sampleId[seqName]
                    rawValidHandle1["rawValidHandle1_" + my_sampleId].write(head1 + "\n" + seq1 + desc1 + qual1)
                    rawValidHandle2["rawValidHandle2_" + my_sampleId].write(head2 + "\n" + seq2 + desc2 + qual2)
        for myHandle in rawValidHandle1.values():
            rawValidHandle1.close()
        for myHandle in rawValidHandle2.values():
            rawValidHandle2.close()

    def run(self):
        super(SeqExtractTool, self).run()
        self.make_ess_dir()
        self.SeqExtract()
        self.end()
