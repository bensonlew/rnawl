# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
# last_modify:2017.06.19

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class QcAndStatModule(Module):
    def __init__(self, work_id):
        super(QcAndStatModule, self).__init__(work_id)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # 输入的fastq文件夹其中包含list文件
            {"name": "insert_size", "type": "infile", "format": "sequence.profile_table"},  # 关于各个样本的insert_size的文件
            {"name": "stat_dir", "type": "infile", "format": "sequence.baif_dir"},  # 输入的碱基质量统计结果文件夹
            {"name": "sickle_dir", "type": "outfile", "format": "sequence.fastq_dir"},  # 设置结果文件后面要用
            {"name": "after_qc_dir", "type": "outfile", "format": "sequence.fastq_dir"}, # 质控结果
            {"name": "before_qc_stat", "type": "outfile", "format": "sequence.profile_table"},  # 原始序列统计信息文件
            {"name": "after_qc_stat", "type": "outfile", "format": "sequence.profile_table"},  # 质控后的高质量序列信息
            {'name': 'fastp', 'type': 'string', 'default': ''}, #宏基因组质控 add by qingchen.zhang 20190507
            {'name': 'qc_quality', 'type': 'int', 'default': 20},
            {'name': 'qc_length', 'type': 'int', 'default': 50},
        ]
        self.add_option(options)
        self.qc = self.add_module("sequence.hiseq_qc_mg")  # 使用统一的质控模块，用参数控制使用该模块
        self.qc_stat = self.add_tool("sequence.raw_qc_stat")  # 质控前后序列信息统计
        self.remove_short_reads = self.add_tool("sequence.remove_short_reads")  # 去除一些比较短的reads
        self.step.add_steps('qc_', 'qc_stat', 'remove_reads')
        self.new_sickle = ''

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def check_options(self):
        if not self.option("fastq_dir").is_set:
            raise OptionError("必须设置参数fastq_dir", code="22701001")
        #if not self.option("stat_dir").is_set:
            #raise OptionError("必须设置参数stat_dir", code="22701002")
        if not self.option("insert_size").is_set:
            raise OptionError("必须设置参数insert_size", code="22701003")
        return True

    def run_qc(self):
        self.qc.set_options({
            "fastq_dir": self.option("fastq_dir"),
            "fq_type": "PSE",
            "quality_q": self.option('qc_quality'),
            "length_q": self.option('qc_length'),
           # "quality_s": 20,
            #"length_s": 50,
        })
        self.qc.on('start', self.set_step, {'start': self.step.qc_})
        self.qc.on('end', self.set_step, {'end': self.step.qc_})
        self.qc.on('end', self.set_step, {'start': self.step.remove_reads})
        self.qc.on('end', self.set_output, "qc")
        self.qc.on('end', self.run_remove)
        self.qc.run()

    def run_fastp_qc(self):
        """
        新软件fastp质控
        :return:
        """
        self.fastp_qc = self.add_module("meta.qc.fastp_qc")
        opts = ({
            "fastq_dir": self.option("fastq_dir"),
            "fq_type": "PE",
            "qualified_quality_phred": str(self.option('qc_quality')),
            "length_required": str(self.option('qc_length'))
        })
        self.fastp_qc.set_options(opts)
        self.fastp_qc.on('end', self.run_fastp_stat)
        self.fastp_qc.run()

    def run_fastp_stat(self):
        """
        fastp质控统计
        :return:
        """
        self.fastp_stat = self.add_tool('meta.qc.fastp_stat')
        opts = ({
            "json_dir": self.fastp_qc.output_dir + '/qc_stat',
        })
        self.fastp_stat.set_options(opts)
        self.fastp_stat.on('end', self.set_output)
        self.fastp_stat.run()

    def run_remove(self):
        self.remove_short_reads.set_options({
            'fastq_dir': self.qc.option('sickle_dir')
        })
        self.remove_short_reads.on('start', self.set_step, {'start': self.step.remove_reads})
        self.remove_short_reads.on('end', self.set_step, {'end': self.step.remove_reads})
        self.remove_short_reads.on('end', self.set_step, {'start': self.step.qc_stat})
        self.remove_short_reads.on('end', self.set_output, "remove_reads")
        self.remove_short_reads.on('end', self.run_qc_stat)
        self.remove_short_reads.run()

    def run_qc_stat(self):
        # self.new_sickle = self.get_new_sickle_dir(self.qc.option('sickle_dir').prop['path'])
        # 处理qc获得的文件的名称，以便后续使用
        self.qc_stat.set_options({
            "base_info_dir": self.option('stat_dir'),
            "sickle_dir":  self.remove_short_reads.option('reasult_dir'),
        })
        self.qc_stat.on('start', self.set_step, {'start': self.step.qc_stat})
        self.qc_stat.on('end', self.set_step, {'end': self.step.qc_stat})
        self.qc_stat.on('end', self.set_output, "qc_stat")
        self.qc_stat.run()

    def get_new_sickle_dir(self, dir_path):  # 这个新文件夹中没有list.txt文件
        new_path = os.path.join(self.work_dir, "new_sickle_dir")
        if os.path.exists(new_path):
            pass
        else:
            os.mkdir(new_path)
        with open(dir_path + "/list.txt", 'r') as r:
            for line in r:
                line = line.strip("\n").split("\t")
                if len(line) == 3:
                    if line[2] == "l":
                        os.link(os.path.join(dir_path, line[0]), os.path.join(new_path, line[1] + ".sickle.1.fq"))
                    else:
                        os.link(os.path.join(dir_path, line[0]), os.path.join(new_path, line[1] + ".sickle.2.fq"))
                else:
                    os.link(os.path.join(dir_path, line[0]), os.path.join(new_path, line[1] + ".sickle.s.fq"))
        return new_path

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'qc':
            self.linkdir(obj.option("sickle_dir").prop["path"], "sickle_dir")
        if event['data'] == 'remove_reads':
            self.linkdir(obj.option('reasult_dir').prop['path'], 'after_qc_dir')
        elif event['data'] == 'qc_stat':
            self.linkdir(obj.output_dir, 'qc_stat')
            self.end()
        if self.option('fastp') != "":
            self.linkdir(self.fastp_qc.option('sickle_dir').prop['path'], 'after_qc_dir')
            self.linkdir(self.fastp_stat.output_dir, 'qc_stat')
            self.end()

    def run(self):
        super(QcAndStatModule, self).run()
        if self.option('fastp') != "":
            self.run_fastp_qc()
        else:
            self.run_qc()

    def end(self):
        try:
            #self.option("sickle_dir", self.output_dir + "/sickle_dir")
            self.option("after_qc_stat", self.output_dir + "/qc_stat/reads.cleanData.stat")
        except Exception as e:
            raise Exception("设置输出文件出错——{}".format(e))
        file_path = self.get_new_reads_stat()
        try:
            self.option("before_qc_stat", file_path)
            self.option('after_qc_dir', self.output_dir + '/after_qc_dir')
        except Exception as e:
            raise Exception("设置输出文件出错——{}".format(e))
        super(QcAndStatModule, self).end()

    def get_new_reads_stat(self):
        old_file_path = self.output_dir + "/qc_stat/reads.rawData.stat"
        new_file_path = self.output_dir + "/qc_stat/new_reads.rawData.stat"  # 含有样品信息的raw.data.stat add by zhujuan 20170925
        if os.path.isfile(new_file_path):
            os.remove(new_file_path)  # 避免重复写入文件 add by guhaidong 20180102
        size = {}
        with open(self.option('insert_size').prop['path'], 'r') as r1:
            for line in r1:
                line = line.strip("\n").split("\t",1)
                size[line[0]] = line[1]
        with open(old_file_path, 'r') as r2 ,open(new_file_path, 'a') as w:
            w.write("#Sample\tSource type\tInsert size(bp)\tRead length(bp)\tRaw reads\tRaw bases(bp)\n")
            for line in r2:
                line = line.strip("\n").split("\t")
                if line[0] in size.keys():
                    w.write(line[0]+"\t"+size[line[0]]+"\t"+line[1]+"\t"+line[2]+"\n")
                else:
                    continue
        return new_file_path

    def linkdir(self, dirpath, dirname):  # 暂时无用
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
                file_name = os.listdir(oldfiles[i])
                os.mkdir(newfiles[i])
                for file_name_ in file_name:
                    os.link(os.path.join(oldfiles[i], file_name_), os.path.join(newfiles[i], file_name_))
