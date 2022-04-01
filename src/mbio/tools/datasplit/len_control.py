# -*- coding: utf-8 -*-
# __author__ = 'xuting'

"""进行最后的长度过滤， 滤去短于某特定长度的序列"""
import os
import errno
import re
import multiprocessing
from collections import defaultdict
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class LenControlAgent(Agent):
    def __init__(self, parent=None):
        super(LenControlAgent, self).__init__(parent)
        self._run_mode = "ssh1"
        options = [
            {'name': 'sample_info', 'type': "infile", 'format': 'datasplit.miseq_split'},  # 样本拆分信息表
            {'name': 'marked_path', 'type': "string"}  # 标记号样本的序列
        ]
        self.add_option(options)

    def check_option(self):
        if not self.option('sample_info').is_set:
            raise OptionError("参数sample_info不能为空")
        if not self.option('marked_path'):
            raise OptionError("参数marked_path不能为空")
        return True

    def set_resource(self):
        self._cpu = 32
        self._memory = ''


class LenControlTool(Tool):
    def __init__(self, config):
        super(LenControlTool, self).__init__(config)
        self._version = 1.0
        self.option('sample_info').get_info()
        self.minLen = dict()
        self.total = defaultdict(int)
        self.valid = defaultdict(int)
        self.GetMinLen()
        self.l = multiprocessing.Lock()

    def GetMinLen(self):
        for c_id in self.option('sample_info').prop['child_ids']:
            self.minLen[c_id] = self.option('sample_info').child_sample(c_id, "filter_min")

    def make_ess_dir(self):
        LenControledPath = os.path.join(self.work_dir, "LenControled")
        dir_list = [LenControledPath]
        for name in dir_list:
            try:
                os.makedirs(name)
            except OSError as exc:
                if exc.errno == errno.EEXIST and os.path.isdir(name):
                    pass
                else:
                    raise OSError("创建目录失败")

    def LenControl(self):
        stat_path = stat_path = os.path.join(self.work_dir, "output", "stat.xls")
        with open(stat_path, 'wb') as w:
            w.write("#library_id\tlibrary_name\ttotal\tvilid\n")
        process_list = list()
        for p in self.option('sample_info').prop["parent_sample"]:
            if p["has_child"]:
                process = multiprocessing.Process(target=self._LenControl, args=(p,))
                process_list.append(process)
        for my_p in process_list:
            my_p.daemon = True
            my_p.start()
        for my_p in process_list:
            my_p.join()

    def _LenControl(self, p):
        p_id = p['sample_id']
        file_ = os.path.join(self.option("marked_path"), p_id + ".marked.fastq")
        fileLenControled = os.path.join(self.work_dir, 'LenControled', p_id + '.lenControl.fastq')
        filediscarded = os.path.join(self.work_dir, 'LenControled', p_id + '.discarded.fastq')
        with open(file_, 'rb') as r, open(fileLenControled, 'wb') as w1, open(filediscarded, 'wb') as w2:
            for head in r:
                self.total[p_id] += 1
                head = head.rstrip('\r\n')
                head = re.sub('^@', '', head)
                sample_id = re.split('_', head)[0]
                seq = r.next().rstrip('\r\n')
                desc = r.next()
                qual = r.next()
                if len(seq) >= int(self.minLen[sample_id]):
                    self.valid[p_id] += 1
                    w1.write("@" + head + "\n" + seq + "\n" + desc + qual)
                else:
                    w2.write("@" + head + "\n" + seq + "\n" + desc + qual)
        self.WriteStat(p)

    def WriteStat(self, p):
        stat_path = stat_path = os.path.join(self.work_dir, "output", "stat.xls")
        p_id = p["sample_id"]
        with open(stat_path, 'ab') as a:
            str_ = "{}\t{}\t{}\t{}\n".format(p_id, p["library_name"], self.total[p_id], self.valid[p_id])
            self.l.acquire()
            a.write(str_)
            self.l.release()

    def run(self):
        super(LenControlTool, self).run()
        self.make_ess_dir()
        self.LenControl()
        self.end()
