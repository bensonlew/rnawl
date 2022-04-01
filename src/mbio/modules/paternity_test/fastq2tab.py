# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# last_modify:2016.11.28

from biocluster.module import Module
import os
import types
from biocluster.core.exceptions import OptionError


class Fastq2tabModule(Module):
    def __init__(self, work_id):
        super(Fastq2tabModule, self).__init__(work_id)
        self.step.add_steps('fastq2bam', 'bam2tab')
        options = [
            {"name": "sample_id", "type": "string"},  # 输入F/M/S的样本ID
            {"name": "fastq_path", "type": "infile","format":"sequence.fastq_dir"},  # fastq所在路径
            {"name": "cpu_number", "type": "int", "default": 4}, #cpu个数
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},  # 参考序列
            {"name": "targets_bedfile", "type": "infile","format":"denovo_rna.gene_structure.bed"}  # 位点信息
        ]
        self.add_option(options)
        self.fastq2bam = self.add_tool("paternity_test.family2bam")
        self.bam2tab = self.add_tool("paternity_test.bam2tab")
        self._end_info = 0

    def check_options(self):
        """
         重写参数检测函数
         :return:
         """
        if not self.option("sample_id"):
            raise OptionError("必须输入样本编号")
        if not self.option("ref_fasta").is_set:
            raise OptionError("必须输入参考基因组的fastq文件")
        if not self.option('fastq_path'):
            raise OptionError('必须提供fastq文件所在的路径')
        if not self.option('targets_bedfile'):
            raise OptionError('必须提供target_bedfile文件')
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def fastq2bam_run(self):
        self.fastq2bam.set_options({
            "fastq": self.option("sample_id"),
            "ref_fasta": self.option("ref_fasta").prop['path'],
            "targets_bedfile": self.option("targets_bedfile").prop['path'],
            "seq_path": self.option("fastq_path").prop['path'],
            "cpu_number": self.option("cpu_number")
        })
        self.fastq2bam.on('end', self.set_output, 'fastq2bam')
        self.fastq2bam.on('start', self.set_step, {'start': self.step.fastq2bam})
        self.fastq2bam.on('end', self.set_step, {'end': self.step.fastq2bam})
        self.fastq2bam.run()

    def bam2tab_run(self):
        # bam_dir = os.path.join(self.work_dir, "Family2bam/output")
        bam_dir = self.fastq2bam.output_dir
        self.bam2tab.set_options({
            "sample_id": self.option("sample_id"),
            "bam_dir": bam_dir,
            "ref_fasta": self.option("ref_fasta").prop['path'],
            "targets_bedfile": self.option("targets_bedfile").prop['path']
        })
        self.bam2tab.on('end', self.set_output, 'bam2tab')
        self.bam2tab.on('start', self.set_step, {'start': self.step.bam2tab})
        self.bam2tab.on('end', self.set_step, {'end': self.step.bam2tab})
        self.bam2tab.run()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
                    # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'fastq2bam':
            self.linkdir(obj.output_dir, 'fastq2bam')
        elif event['data'] == 'bam2tab':
            self.linkdir(obj.output_dir, 'bam2tab')
        else:
            pass

    def run(self):
        super(Fastq2tabModule, self).run()
        self.fastq2bam.on('end', self.bam2tab_run)
        self.fastq2bam_run()
        self.bam2tab.on('end', self.end)

    def end(self):
        repaths = [
            [".", "", "无创亲子鉴定结果输出目录"],
        ]
        regexps = [
            [r"fastq2bam/.*\.filter\.bam$", "bam", "mem.sort.hit.filter.bam文件"],
            [r"bam2tab/.*\.mem\.sort\.hit\.vcf\.tab", "tab", "mem.sort.hit.vcf.tab文件"],
            [r"bam2tab/.*\.qc", "qc", "质控数据"]
        ]

        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        print regexps
        super(Fastq2tabModule, self).end()