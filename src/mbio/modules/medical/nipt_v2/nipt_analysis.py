# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'

from biocluster.module import Module
import os
from biocluster.core.exceptions import OptionError


class NiptAnalysisModule(Module):
    def __init__(self, work_id):
        super(NiptAnalysisModule, self).__init__(work_id)
        self.step.add_steps('fastq2bed', 'bed_analysis', 'fastqc', 'identification')
        options = [
            {"name": "fastq_path", "type": "infile", "format": "sequence.fastq_dir"},
            {"name": "sample_id", "type": "string"},
            {"name": "bw", "type": "int", "default": 10},
            {"name": "bs", "type": "int", "default": 1},
            {"name": "ref_group", "type": "int", "default": 1},
            {"name": "single", "type": "string", "default": "false"}
        ]
        self.add_option(options)
        self.fastq2bed = self.add_tool("medical.nipt_v2.fastq_process")
        self.bin_step = [{'bin': 10, 'step': 1}, {'bin': 5, 'step': 5}, {'bin': 1, 'step': 1},
                         {'bin': 500, 'step': 500}]
        self.tools = []
        self.fastqc = None

    def check_options(self):
        """
         重写参数检测函数
         :return:
         """
        if not self.option('fastq_path').is_set:
            raise OptionError('必须提供fastq文件所在的路径！')
        if not self.option('sample_id'):
            raise OptionError('必须提供要进行分析的样本id！')
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def fastq2bed_run(self):
        self.fastq2bed.set_options({
            "sample_id": self.option('sample_id'),
            "fastq_path": self.option("fastq_path"),
            "single": self.option("single")
        })
        self.fastq2bed.on('end', self.set_output, 'fastq2bed')
        self.fastq2bed.on('start', self.set_step, {'start': self.step.fastq2bed})
        self.fastq2bed.on('end', self.set_step, {'end': self.step.fastq2bed})
        self.fastq2bed.on('end', self.bed_run)
        self.fastq2bed.run()

    def bed_run(self):
        bed_dir = self.fastq2bed.output_dir
        # bed_dir = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/tool_test/nipt"
        n = 0
        for m in self.bin_step:
            bed_analysis = self.add_tool("medical.nipt_v2.bed_analysis")
            self.step.add_steps('bed_analysis{}'.format(n))
            bed_analysis.set_options({
                "bed_file": bed_dir + '/' + self.option('sample_id') + '.bed.2',
                # "bed_file": bed_dir + '/' + 'WS17100541-2.bed',
                "bw": m['bin'],
                'bs': m['step'],
                'ref_group': self.option('ref_group'),
                "single_chr": "false" if m['bin'] != 500 else 'true'
            })
            step = getattr(self.step, 'bed_analysis{}'.format(n))
            step.start()
            bed_analysis.on('end', self.finish_update, 'bed_analysis{}'.format(n))
            self.tools.append(bed_analysis)
            n += 1
        for j in range(len(self.tools)):
            self.tools[j].on('end', self.set_output, 'bed_analysis')
        if len(self.tools) > 1:
            self.on_rely(self.tools, self.fastqc_run)
        elif len(self.tools) == 1:
            self.tools[0].on('end', self.fastqc_run)
        for tool in self.tools:
            tool.run()

    def fastqc_run(self):
        self.fastqc = self.add_tool("medical.nipt_v2.fastqc")
        bed_dir = self.fastq2bed.output_dir
        self.fastqc.set_options({
            "sample_id": self.option('sample_id'),
            "fastq_path": self.option('fastq_path'),
            "bam_file": bed_dir + '/' + self.option('sample_id') + '.map.valid.bam',
        })
        self.fastqc.on('end', self.set_output, 'fastqc')
        self.fastqc.on('start', self.set_step, {'start': self.step.fastqc})
        self.fastqc.on('end', self.set_step, {'end': self.step.fastqc})
        self.fastqc.on('end', self.end)
        self.fastqc.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'fastq2bed' or event['data'] == 'bed_analysis' or event['data'] == 'fastqc':
            allfiles = os.listdir(obj.output_dir)
            oldfiles = [os.path.join(obj.output_dir, i) for i in allfiles]
            newfiles = [os.path.join(self.output_dir, i) for i in allfiles]
            for i in range(len(allfiles)):
                os.link(oldfiles[i], newfiles[i])

    def run(self):
        super(NiptAnalysisModule, self).run()
        self.fastq2bed_run()

    def end(self):
        repaths = [
            [".", "", "无创产前筛查结果输出目录"],
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        super(NiptAnalysisModule, self).end()
