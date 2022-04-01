# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 20181107

import os
from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class RnaQcModule(Module):
    """
    常规RNA质控
    方法：SeqPrep、sickle
    """
    def __init__(self, work_id):
        super(RnaQcModule, self).__init__(work_id)
        options = [
            {"name": "sample_path", "type": "infile", "format": "datasplit.list_file", "required": True},  #第一列fastq路径，第二列样本名，第三列序列类型 l or r
            {"name": "quality", "type": "int", "default": 20},
            {"name": "length", "type": "int", "default": 30},
            {"name": "adapter_a", "type": "string", "default": "AGATCGGAAGAGCACACGTC"},
            {"name": "adapter_b", "type": "string", "default": "AGATCGGAAGAGCGTCGTGT"},
        ]
        self.add_option(options)
        self.sample_info = {}
        self.end_times = 0

    def check_options(self):
        pass

    def get_info(self):
        self.sample_list = []
        with open(self.option("sample_path").prop["path"], "r") as f:
            for line in f:
                item = line.strip().split('\t')
                if item[1] in self.sample_list:
                    if item[2] == 'l':
                        self.sample_info[item[1]].insert(0, item[0])
                    else:
                        self.sample_info[item[1]].append(item[0])
                else:
                    self.sample_list.append(item[1])
                    self.sample_info[item[1]] = []
                    self.sample_info[item[1]].append(item[0])
            for key in self.sample_info.keys():
                if len(self.sample_info[key]) > 2:
                    self.set_error('常规RNA质控的样本：{}有重名，请改样本名或分开质控！'.format(key))
                elif len(self.sample_info[key]) < 2:
                    self.set_error('常规RNA质控的样本：{}对应的R1,R2序列不全,请核实！'.format(key))

    def run_seq_prep(self):
        """
        tool:seq_prep
        """
        for sample in self.sample_info.keys():
            options = {
                "fastq_r": self.sample_info[sample][1],
                "fastq_l": self.sample_info[sample][0],
                "quality": self.option("quality"),
                "length": self.option("length"),
                "adapter_a": self.option("adapter_a"),
                "adapter_b": self.option("adapter_b")
            }
            self.seq_prep = self.add_tool("datasplit.seq_prep")
            self.seq_prep.set_options(options)
            self.seq_prep.on("end", self.run_sickle, sample)
            self.seq_prep.run()

    def run_sickle(self, event):
        """
        tool:sickle_mg
        """
        obj = event["bind_object"]
        options = {
            "fq_type": 'PE',
            "fastq_l": obj.option("seqprep_l"),
            "fastq_r": obj.option("seqprep_r"),
            "length": self.option("length")
        }
        self.sickle = self.add_tool('datasplit.sickle_mg')
        self.sickle.set_options(options)
        self.sickle.on("end", self.set_output, event["data"])
        self.sickle.run()

    def set_output(self, event):
        obj = event["bind_object"]
        for f in os.listdir(obj.output_dir):
            old = os.path.join(obj.output_dir, f)
            if f.startswith("sickle_l.fastq.gz"):
                new = os.path.join(self.output_dir, event["data"] + ".clean.1.fastq.gz")
                if os.path.exists(new):
                    os.remove(new)
                os.link(old, new)
            if f.startswith("sickle_r.fastq.gz"):
                new = os.path.join(self.output_dir, event["data"] + ".clean.2.fastq.gz")
                if os.path.exists(new):
                    os.remove(new)
                os.link(old, new)
        self.end_times += 1
        if self.end_times == len(self.sample_list):
            self.end()

    def run(self):
        super(RnaQcModule, self).run()
        self.get_info()
        self.run_seq_prep()

    def end(self):
        super(RnaQcModule, self).end()
