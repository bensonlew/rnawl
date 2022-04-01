# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
from biocluster.core.exceptions import OptionError
import subprocess
from Bio import SeqIO


class PredictStatAgent(Agent):
    """
    根据gff文件夹统计出tRNA个数、

    """

    def __init__(self, parent):
        super(PredictStatAgent, self).__init__(parent)
        options = [
            {"name": "merge_dir", "type": "infile", "format": "paternity_test.data_dir"},  # 输入合并文件夹
            {"name": "sample", "type": "string"} #样本名称
        ]
        self.add_option(options)
        self.step.add_steps('predict_stat')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.predict_stat.start()
        self.step.update()

    def step_end(self):
        self.step.predict_stat.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("merge_dir").is_set:
            raise OptionError("请传入merge_dir文件夹路径")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(PredictStatAgent, self).end()

class PredictStatTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(PredictStatTool, self).__init__(config)
        self.sh_path = "../../../../../.." + self.config.PACKAGE_DIR + '/sequence/scripts/'
        self.merge_file = self.option('merge_dir').prop['path']
        self.s16_path = ''
        for file in os.listdir(self.merge_file):
            if re.search(r'16S.fna', file):
                self.s16_path = os.path.join(self.merge_file, file)

    def stat_table(self):
        """
        统计文件
        :return:
        """
        self.logger.info('正在进行结果统计')
        stat_dir = self.option('merge_dir').prop['path']
        stat_file = self.work_dir + '/' + self.option('sample') + '_sample_stat.xls'
        s5 = 0
        s16 = 0
        s23 = 0
        trna_num = 0
        rrna_num = 0
        cds_num = 0
        for file in os.listdir(stat_dir):
            if re.search(r'CDS\.gff', file):
                gff_cds = os.path.join(self.option('merge_dir').prop['path'], file)
                with open(gff_cds, 'r') as f:
                    lines = f.readlines()
                    cds_num = len(lines) - 1
            if re.search(r'tRNA\.gff', file):
                gff_trna = os.path.join(self.option('merge_dir').prop['path'], file)
                with open(gff_trna, 'r') as f:
                    lines = f.readlines()
                    trna_num = len(lines) - 1
            if re.search(r'rRNA\.gff', file):
                gff_rrna = os.path.join(self.option('merge_dir').prop['path'], file)
                with open(gff_rrna, 'r') as f:
                    lines = f.readlines()
                    rrna_num = len(lines) - 1
                    for line in lines:
                        line = line.strip().split('\t')
                        if re.search(r'5S', line[7]):
                            s5 += 1
                        if re.search(r'16S', line[7]):
                            s16 += 1
                        if re.search(r'23S', line[7]):
                            s23 += 1

        with open(stat_file, 'w') as w:
            w.write('Sample\tCDS_num\trRNA_num\t16sRNA\t23sRNA\t5sRNA\ttRNA\n')
            w.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(self.option('sample'),
                        cds_num, rrna_num, s16, s23,s5,trna_num))
        self.logger.info('统计文件数结束啦')

    def run_seq(self):
        """
        对合并后的16S序列重新进行排序
        :return:
        """
        seq_id = []
        seq_length = []
        total_seq = {}
        input_dir = self.option('merge_dir').prop['path']
        for file in os.listdir(input_dir):
            if re.search(r'16S.fna', file):
                s16_path = os.path.join(input_dir, file)
                for seq_record in SeqIO.parse(s16_path, 'fasta'):
                    id = seq_record.id
                    length = int(len(seq_record.seq))
                    seq_length.append(length)
                    total_seq[id] = seq_record.seq
        new_seq_length = sorted(seq_length, reverse=True)
        new_path = self.output_dir + '/' + self.option('sample') + '_16S.fna'
        if os.path.exists(new_path):
            os.remove(new_path)
        with open(new_path, 'w') as w:
            for seq in new_seq_length:
                for _id in total_seq.keys():
                    length = int(len(total_seq[_id]))
                    if length == seq:
                        if _id not in seq_id:
                            seq_id.append(_id)
                            w.write('>{}\n{}\n'.format(_id,total_seq[_id]))


    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        file_path = self.output_dir + '/' + self.option('sample') + '_sample_stat.xls'
        if os.path.exists(file_path):
            os.remove(file_path)
        os.link(self.work_dir + '/' + self.option('sample') + '_sample_stat.xls', file_path)

    def run(self):
        super(PredictStatTool, self).run()
        self.stat_table()
        if self.s16_path !='':
            self.run_seq()
        self.set_output()
        self.end()