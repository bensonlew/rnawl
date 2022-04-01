# -*- coding: utf-8 -*-
# __author__ = 'moli.zhou'
# last_modify:2016.12.2

from biocluster.module import Module
import os
import types
from biocluster.core.exceptions import OptionError


class RetabModule(Module):
    def __init__(self, work_id):
        super(RetabModule, self).__init__(work_id)
        self.step.add_steps('fastq2bam', 're_bam2tab')
        options = [
            {"name": "sample_id1", "type": "string"},  # 输入F/M/S的样本ID
            {"name": "sample_id2", "type": "string"},  # 输入F/M/S的样本ID
            {"name": "fastq_path", "type": "string"},  # fastq所在路径
            {"name": "cpu_number", "type": "int", "default": 4}, #cpu个数
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},  # 参考序列
            {"name": "targets_bedfile", "type": "string"}  # 位点信息
        ]
        self.add_option(options)
        self.re_bam2tab = self.add_tool("paternity_test.re_bam2tab")
        # self.fastq2bam = self.add_tool("paternity_test.family2bam")
        # self.step.add_steps('fastq2bam')
        self._end_info = 0
        self.tools = []
        self.sum_tools = []

    def check_options(self):
        """
         重写参数检测函数
         :return:
         """
        if not self.option("sample_id1"):
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

    # def fastq2bam_run(self):
    #     self.fastq2bam.set_options({
    #         "fastq": self.option('sample_id1'),
    #         "ref_fasta": self.option("ref_fasta"),
    #         "targets_bedfile": self.option("targets_bedfile"),
    #         "seq_path": self.option("fastq_path"),
    #         "cpu_number": self.option("cpu_number")
    #     })
    #     self.fastq2bam.on('end', self.set_output, 'fastq2bam{}'.format(n))
    #     self.fastq2bam.on('start', self.set_step, {'start': self.step.fastq2bam})
    #     self.fastq2bam.on('end', self.set_step, {'end': self.step.fastq2bam})
    #     self.fastq2bam.run()

    def finish_update(self,event):
        step = getattr(self.step,event['data'])
        step.finish()
        self.step.update()

    def fastq2bam_run(self):
        n = 1
        samples = [self.option('sample_id1'), self.option('sample_id2')]
        for f in samples:
            fastq2bam = self.add_tool("paternity_test.family2bam")
            self.step.add_steps('Fastq2bam{}'.format(n))
            fastq2bam.set_options({
                "fastq": f,
                "ref_fasta": self.option("ref_fasta"),
                "targets_bedfile": self.option("targets_bedfile"),
                "seq_path": self.option("fastq_path"),
                "cpu_number": self.option("cpu_number")
            })
            step = getattr(self.step, 'Fastq2bam{}'.format(n))
            step.start()
            fastq2bam.on('end',self.finish_update,'Fastq2bam{}'.format(n))
            self.tools.append(fastq2bam)
            n += 1
        self.tools[0].on('end', self.set_output)
        self.tools[1].on('end', self.set_output)

        for tool in self.tools:
            tool.run()

    def re_bam2tab_run(self):
        self.re_bam2tab.set_options({
            "bam_file1": self.work_dir + '/Family2bam/output/' + self.option("sample_id1") + '.mem.sort.hit.filter.bam',
            "bam_file2": self.work_dir + '/Family2bam1/output/' + self.option("sample_id2") + '.mem.sort.hit.filter.bam',
            "ref_fasta": self.option("ref_fasta"),
            "targets_bedfile": self.option("targets_bedfile")
        })
        self.re_bam2tab.on('end', self.set_output, 're_bam2tab')
        self.re_bam2tab.on('start', self.set_step, {'start': self.step.re_bam2tab})
        self.re_bam2tab.on('end', self.set_step, {'end': self.step.re_bam2tab})
        self.re_bam2tab.run()

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
        elif event['data'] == 're_bam2tab':
            self.linkdir(obj.output_dir, 're_bam2tab')
        else:
            pass

    def run(self):
        super(RetabModule, self).run()
        self.fastq2bam_run()
        self.on_rely(self.tools, self.re_bam2tab_run)
        self.re_bam2tab.on('end', self.end)

    def end(self):
        repaths = [
            [".", "", "无创亲子鉴定结果输出目录"],
        ]
        regexps = [
            [r"ReBam2tab/.*tab$", "tab", "mem.sort.hit.vcf.tab文件"],
        ]

        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        print regexps
        super(RetabModule, self).end()
