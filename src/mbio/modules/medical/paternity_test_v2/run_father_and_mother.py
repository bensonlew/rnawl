#!/usr/bin/env python
# -*- coding: utf-8 -*-
# kefei.huang
# last modified by hongdongxuan @ 20171124
from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class RunFatherAndMotherModule(Module):
    def __init__(self, work_id):
        super(RunFatherAndMotherModule, self).__init__(work_id)
        options = [
            {"name": "sample_id", "type": "string"},  # 输入F/M/S的样本ID
            {"name": "fastq_path", "type": "string"},  # fastq所在路径
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},  # 参考序列
            {"name": "targets_bedfile", "type": "string"},  # 位点信息
            {"name": "batch_id", "type": "string"},
            {"name": "board_batch", "type": "string"},
            {"name": "analysis_type", "type": "string"}
        ]
        self.add_option(options)
        self.modules_ = []

    def check_options(self):
        if not self.option("sample_id"):
            raise OptionError("必须输入样本编号")
        if not self.option("ref_fasta").is_set:
            raise OptionError("必须输入参考基因组的fastq文件")
        if not self.option('fastq_path'):
            raise OptionError('必须提供fastq文件所在的路径')
        if not self.option('targets_bedfile'):
            raise OptionError('必须提供target_bedfile文件')
        if not self.option("batch_id"):
            raise OptionError("必须输入batch_id")
        if not self.option("board_batch"):
            raise OptionError("必须输入board_batch")
        if not self.option('analysis_type'):
            raise OptionError("必须指定测序类型")
        if self.option('analysis_type') != "dcpt" and self.option('analysis_type') != "pt":
            raise OptionError("测序类型必须是dcpt/pt")
        return True

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def dcpt_run(self):
        n = 0
        all_samples = self.option("sample_id").split(",")
        for m in all_samples:
            single_run = self.add_module("medical.paternity_test_v2.fastq_call_snp")
            self.step.add_steps('call_snp{}'.format(n))
            single_run.set_options({
                "sample_id": m,
                "fastq_path": self.option("fastq_path"),
                "ref_fasta": self.option("ref_fasta"),
                "targets_bedfile": self.option("targets_bedfile"),
                "batch_id": self.option("batch_id"),
                "board_batch": self.option("board_batch"),
                "analysis_type": self.option("analysis_type")
            })
            step = getattr(self.step, 'call_snp{}'.format(n))
            step.start()
            single_run.on('end', self.finish_update, 'call_snp{}'.format(n))
            self.modules_.append(single_run)
            n += 1
        self.on_rely(self.modules_, self.end)
        for module in self.modules_:
            module.run()

    def run(self):
        super(RunFatherAndMotherModule, self).run()
        self.dcpt_run()

    def end(self):
        super(RunFatherAndMotherModule, self).end()
