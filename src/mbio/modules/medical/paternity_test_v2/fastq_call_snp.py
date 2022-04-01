#!/usr/bin/env python
# -*- coding: utf-8 -*-
# kefei.huang
# 20171120

import os
from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class FastqCallSnpModule(Module):
    """
        这个程序的作用是，把需要的文件读进来，然后启动call_snp的模块。
        目前只允许F和M多重PCR的样品进来
        以后可能需要运行S杂捕的样品进来
        预期在不久的将来，可能要放进F/M/S的多重样品
        author: kefei.huang
        last modified by hongdong 20171124
        last modified by kefeihuang 20171218
    """
    def __init__(self, work_id):
        super(FastqCallSnpModule, self).__init__(work_id)
        self.step.add_steps('fastq2bam', 'bam2tab', 'bam2tab_dc', 'pol_degree')
        options = [
            {"name": "sample_id", "type": "string"},  # 输入F/M/S的样本ID
            {"name": "fastq_path", "type": "infile", "format": "sequence.fastq_dir"},  # fastq所在路径
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},  # 参考序列
            {"name": "targets_bedfile", "type": "string"},  # 位点信息
            {"name": "batch_id", "type": "string"},
            {"name": "board_batch", "type": "string"},  # 版号170928_TPNB500180_0112_AHMWJKAFXX
            {"name": "analysis_type", "type": "string"},
            {"name": "is_update", "type": "string", "default": "true"}  # 用于设定该流程中样本call 完snp 要不要更新状态
        ]
        self.add_option(options)
        self.fastq2bam = self.add_tool("medical.paternity_test_v2.fastq2bam")
        self.bam2tab = self.add_tool("medical.paternity_test_v2.bam2tab")
        self.bam2tab_dc = self.add_tool("medical.paternity_test_v2.bam2tab_dc")
        self.pol_degree = self.add_tool("medical.paternity_test_v2.pol_degree")

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
        self.logger.info(self.option('fastq_path'))
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def fastq2bam_run(self):
        self.logger.info(self.option("fastq_path"))
        is_dcpt = True if self.option('analysis_type') == "dcpt" else False
        self.fastq2bam.set_options({
            "fastq": self.option("sample_id"),
            "is_dcpt": is_dcpt,
            "ref_fasta": self.option("ref_fasta"),
            "targets_bedfile": self.option("targets_bedfile"),
            "seq_path": self.option("fastq_path"),
            "cpu_number": 5,
            "memory": 10
        })
        self.fastq2bam.on('end', self.set_output, 'fastq2bam')
        self.fastq2bam.on('start', self.set_step, {'start': self.step.fastq2bam})
        self.fastq2bam.on('end', self.set_step, {'end': self.step.fastq2bam})
        self.fastq2bam.run()

    def bam2tab_run(self):
        bam_path = os.path.join(self.fastq2bam.output_dir, self.option("sample_id") + ".mem.sort.hit.filter.bam")
        self.logger.info(bam_path)
        self.bam2tab.set_options({
            "sample_id": self.option("sample_id"),
            "bam_path": bam_path,
            "ref_fasta": self.option("ref_fasta"),
            "targets_bedfile": self.option("targets_bedfile"),
            "batch_id": self.option('batch_id'),
            "board_batch": self.option('board_batch'),
            "cpu_number": 4,
            "memory": 20,
            "is_update": self.option('is_update')
        })
        self.bam2tab.on('end', self.set_output, 'bam2tab')
        self.bam2tab.on('start', self.set_step, {'start': self.step.bam2tab})
        self.bam2tab.on('end', self.set_step, {'end': self.step.bam2tab})
        # self.bam2tab.on('end', self.end)
        self.bam2tab.run()

    def bam2tab_dc_run(self):
        bam_path = os.path.join(self.fastq2bam.output_dir, self.option("sample_id") + ".mem.sort.hit.filter.bam")
        self.logger.info(bam_path)
        self.bam2tab_dc.set_options({
            "sample_id": self.option("sample_id"),
            "bam_path": bam_path,
            "ref_fasta": self.option("ref_fasta"),
            "targets_bedfile": self.option("targets_bedfile"),
            "batch_id": self.option('batch_id'),
            "board_batch": self.option('board_batch'),
            "cpu_number": 4,
            "memory": 10,
            "is_update": self.option('is_update')
        })
        self.bam2tab_dc.on('end', self.set_output, 'bam2tab_dc')  # 如果没有这个文件夹就会创建文件夹
        self.bam2tab_dc.on('start', self.set_step, {'start': self.step.bam2tab_dc})
        self.bam2tab_dc.on('end', self.set_step, {'end': self.step.bam2tab_dc})
        # self.bam2tab_dc.on('end', self.end)
        self.bam2tab_dc.run()

    def pol_degree_run(self):
        if self.option('analysis_type') == "dcpt":
            ref_dp_xls = os.path.join(self.bam2tab_dc.output_dir, self.option("sample_id") + ".ref.dp.xls")
            dp = self.api.api('medical.paternity_test_v2').query_dp(self.option("sample_id"))
            print 'dp:', dp
            self.pol_degree.set_options({
                "ref_dp_xls": ref_dp_xls,
                "sample_id": self.option("sample_id"),
                "dp": dp                                          # 多重需要提供dp的值覆盖tool中的默认值
            })
        elif self.option('analysis_type') == "pt":
            ref_dp_xls = os.path.join(self.bam2tab.output_dir, self.option("sample_id") + ".ref.dp.xls")
            self.pol_degree.set_options({
                "ref_dp_xls": ref_dp_xls,
                "sample_id": self.option("sample_id"),
            })
        else:
            return

        self.pol_degree.on('end', self.set_output, 'pol_degree')
        self.pol_degree.on('start', self.set_step, {'start': self.step.pol_degree})
        self.pol_degree.on('end', self.set_step, {'end': self.step.pol_degree})
        self.pol_degree.on('end', self.end)
        self.pol_degree.run()

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
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'fastq2bam':
            self.linkdir(obj.output_dir, 'fastq2bam')  # obj.output_dir是tools的内置参数
        elif event['data'] == 'bam2tab' or event['data'] == 'bam2tab_dc':
            self.linkdir(obj.output_dir, "bam2tab")
        elif event['data'] == 'pol_degree':
            self.linkdir(obj.output_dir, 'pol_degree')

    def run(self):
        super(FastqCallSnpModule, self).run()
        if self.option('analysis_type') == "pt":
            self.fastq2bam.on('end', self.bam2tab_run)
            self.bam2tab.on('end', self.pol_degree_run)
        elif self.option('analysis_type') == "dcpt":
            self.fastq2bam.on('end', self.bam2tab_dc_run)
            self.bam2tab_dc.on('end', self.pol_degree_run)
        else:
            pass
        self.fastq2bam_run()

    def end(self):
        repaths = [
            [".", "", "无创亲子鉴定结果输出目录"],
        ]
        regexps = [
            [r"fastq2bam/.*\.filter\.bam$", "bam", "mem.sort.hit.filter.bam文件"],
            [r"bam2tab/.*\.mem\.sort\.hit\.vcf\.tab", "tab", "mem.sort.hit.vcf.tab文件"],
            [r"bam2tab/.*\.qc", "qc", "质控数据"]
        ]
        sdir = self.add_upload_dir(self.output_dir)  # 把必要文件放到一个统一的目录下。
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        # print regexps
        super(FastqCallSnpModule, self).end()
