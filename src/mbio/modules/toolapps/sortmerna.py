# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class SortmernaModule(Module):
    """
    对每个文件质控后的reads进行分析
    """
    def __init__(self, work_id):
        super(SortmernaModule, self).__init__(work_id)
        options = [
            {'name': 'in_fastq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq文件夹
            {'name': 'qc', 'type': 'bool', 'default': False},  # 是否需要质控
            {'name': 'qc_quality', 'type': 'int', 'default': 20},  # 质控质量值标准
            {'name': 'qc_length', 'type': 'int', 'default': 30},  # 质控最短序列长度
            {'name': 'rm_host', 'type': 'bool', 'default': False},  # 是否需要去除宿主
            {'name': 'ref_database', 'type': 'string', 'default': ''},  # 宿主参考序列库中对应的物种名，eg：E.coli ,B.taurus
            {'name': 'second_ref', 'type': 'string', 'default': ''}, #参考数据库二级下拉框数据
            {'name': 'ref_undefined', "type": 'infile', 'format': 'sequence.fasta_dir'},
            {'name': 'ref_undefined_name', 'type': 'string', 'default': 'undefined'},  # 自定义参考宿主名称，适应页面参数
            {'name': 'clean_fq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq文件夹
            {"name": "database", "type": 'string', "default": "rfam_5.8s,rfam_5s,arc_16s,arc_23s,bac_16s,bac_23s,euk_18s,euk_28s"}
        ]
        self.add_option(options)
        self.step.add_steps('data', 'reads_qc', 'rm_host', 'clean_data', 'sormerna') #module没法调用
        self.sequence = self.add_module('metagenome.reads_unzip')
        self.qc = self.add_module('meta.qc.fastp_qc')
        self.rm_host = self.add_module('meta.qc.bwa_remove_host')
        self.clean_data = self.add_module("sequence.metagbin_clean_fq")
        #self.sortmerna = self.add_module('metagenome.sortmerna')
        self.fastq_split = self.add_tool('tool_lab.split_fastq_dir')
        self.sortmerna_list = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('in_fastq').is_set and not self.option('clean_fq').is_set:
            raise OptionError('请输入原始序列文件夹或者质控优质序列文件夹！')
        if self.option('rm_host'):
            if self.option('ref_database') == '' and not self.option('ref_undefined').is_set:
                raise OptionError('已选择去宿主，需输入参考数据库或参考序列')
            if self.option('ref_database') not in ['', 'Custom'] and self.option('ref_undefined').is_set:
                raise OptionError('去宿主不可同时提供参考数据库及参考序列')

    def run_sequence(self):
        opts = {
            'fastq_dir': self.option('in_fastq'),
        }
        self.sequence.set_options(opts)
        self.sequence.on('end', self.set_output, 'sequence')
        self.sequence.run()

    def run_qc(self):
        opts = {
            'fastq_dir': self.sequence.output_dir + '/data',
        }
        self.qc.set_options(opts)
        self.qc.on('end', self.set_output, 'reads_qc')
        self.qc.run()

    def run_rm_host(self):
        ref_database = self.option('ref_database') + "," + self.option('second_ref')
        opts = {
            'fq_type': 'PE',
            'ref_database': ref_database,
            'ref_undefined': self.option('ref_undefined'),
        }
        if self.option('ref_database') == 'Custom':
            opts['ref_database'] = ""
        if self.option('qc'):
            opts['fastq_dir'] = self.qc.output_dir + '/after_qc_dir'
        else:
            opts['fastq_dir'] = self.option('in_fastq')
        self.rm_host.set_options(opts)
        self.rm_host.on('end', self.set_output, 'rm_host')
        self.rm_host.run()

    def run_clean_sequence(self):
        opts = {
            'fastq_dir': self.option('clean_fq'),
        }
        self.clean_data.set_options(opts)
        self.clean_data.on('end', self.set_output, 'clean_data')
        self.clean_data.run()

    def run_split(self):
        if self.option("qc"):
            if self.option('rm_host'):
                fastq_dir = self.rm_host.output_dir
            else:
                fastq_dir = self.qc.output_dir + '/after_qc_dir'
        else:
            fastq_dir = self.clean_data.output_dir + '/data'

        opts = {
            "in_fastq": fastq_dir,
        }
        self.fastq_split.set_options(opts)
        self.fastq_split.run()

    def cat_all(self):
        result_dir = {}
        with open(self.fastq_split.work_dir+"/sample_split_file.txt") as f:
            data = f.readlines()
            for x in data:
                if x.strip().split("\t")[0] in result_dir:
                    result_dir[x.strip().split("\t")[0]][x.strip().split("\t")[1]] = ""
                else:
                    result_dir[x.strip().split("\t")[0]]={x.strip().split("\t")[1]:""}
        for sample in result_dir:
            for dir in result_dir[sample]:
                for tool in self.sortmerna_list:
                    if os.path.exists(tool.output_dir + "/" + dir + ".unalign.1.fq.gz"):
                        result_dir[sample][dir] = tool.output_dir
        for sample in result_dir:
            tmp_list = []
            for dir in result_dir[sample]:
                tmp_list.append(result_dir[sample][dir]+"/"+dir)
            for name in [".align.1.fq.gz",".align.2.fq.gz",".unalign.1.fq.gz",".unalign.2.fq.gz"]:
                tmp_list2 = []
                for xx in tmp_list:
                    tmp_list2.append(xx+name)
                os.system("cat {} > {}".format(" ".join(tmp_list2),self.output_dir + "/" + sample + name))

        with open(self.output_dir + "/result.txt","w") as t:
            t.write("Sample\tclean_reads\trRNA_reads\trRNA_ratio\tnon_rRNA_reads\tnon_rRNA_ratio\n")
            for sample in result_dir:
                align = os.popen("zcat {} | wc -l".format(self.output_dir + "/" + sample+".align.1.fq.gz"))
                unalign = os.popen("zcat {} | wc -l".format(self.output_dir + "/" + sample+".unalign.1.fq.gz"))
                align_num = int(align.read()) /2
                unalign_num =  int(unalign.read()) /2
                clean_base = align_num + unalign_num
                rRNA_ratio = round((float(align_num)/float(clean_base))*100,2)
                non_rRNA_ratio = 100 - rRNA_ratio
                t.write(sample + "\t" + str(clean_base) + "\t" + str(align_num) + "\t" + str(rRNA_ratio) + "\t" + str(
                    unalign_num) + "\t" + str(non_rRNA_ratio) + "\n")

    def run_sortmerna(self):
        for dir in os.listdir(self.fastq_split.output_dir):
            sortmerna = self.add_module('metagenome.sortmerna')
            sortmerna.set_options({
                'fa_dir': self.fastq_split.output_dir + "/" + dir,
                'database': self.option("database"),
            })
            self.sortmerna_list.append(sortmerna)
        if len(self.sortmerna_list) > 1:
            self.on_rely(self.sortmerna_list, self.end)
        else:
            self.sortmerna_list[0].on('end', self.end)
        for tool in self.sortmerna_list:
            tool.run()

    def run(self):
        """
        运行
        :return:
        """
        super(SortmernaModule, self).run()
        if self.option('qc'):
            self.sequence.on('end', self.run_qc)
            if self.option('rm_host'):
                self.qc.on('end', self.run_rm_host)
                self.rm_host.on("end", self.run_split)
                self.fastq_split.on("end", self.run_sortmerna)
            else:
                self.qc.on('end', self.run_split)
                self.fastq_split.on("end", self.run_sortmerna)
            self.run_sequence()
        else:
            self.clean_data.on("end",self.run_split)
            self.fastq_split.on("end", self.run_sortmerna)
            self.run_clean_sequence()

    def set_output(self, event):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        pass
        """
        obj = event['bind_object']
        if event['data'] == 'sortmerna':
            if len(os.listdir(self.output_dir)) >=1:
                for i in os.listdir(self.output_dir):
                    os.remove(self.output_dir + "/" + i)
            for i in os.listdir(obj.output_dir):
                os.link(obj.output_dir + "/" + i, self.output_dir + "/" + i)
            self.end()
        """

    def end(self):
        self.cat_all()
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "sortmerna结果输出目录"],
        ])
        super(SortmernaModule, self).end()