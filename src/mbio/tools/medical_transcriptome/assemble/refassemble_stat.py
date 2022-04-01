# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue,shicaiping'

import glob
import os
import re
import shutil

from Bio import SeqIO

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.ref_rna_v2.trans_step import step_count, class_code_count, gene_trans_exon, count_trans_or_exons


class RefassembleStatAgent(Agent):
    """
    有参组装数据统计
    """

    def __init__(self, parent):
        super(RefassembleStatAgent, self).__init__(parent)
        options = [
            {"name": "all_files_dir", "type": "string"},  # 输入的文件夹
            {"name": "assemble_method", "type": "string", "default": "cufflinks"},  # 选择拼接软件
        ]
        self.add_option(options)
        self.step.add_steps("RefassembleStat")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.RefassembleStat.start()
        self.step.update()

    def stepfinish(self):
        self.step.RefassembleStat.finish()
        self.step.update()

    def check_options(self):
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(RefassembleStatAgent, self).end()


class RefassembleStatTool(Tool):
    def __init__(self, config):
        super(RefassembleStatTool, self).__init__(config)

    def run(self):
        """
        运行
        :return:
        """
        super(RefassembleStatTool, self).run()
        self.get_numberlist()
        self.trans_stat()
        self.count_genes_trans_exons()
        self.set_output()
        self.end()

    def get_numberlist(self):
        """
        获得文件夹中每个gtf文件中转录本和基因的数量
        :return:
        """
        file_list = []
        numberlist_path = os.path.join(self.work_dir, "number_list.txt")
        with open(numberlist_path, "w+") as w:
            w.write("#file_names\ttrans\tgenes\n")
            a = os.listdir(self.option('all_files_dir'))
            for f in a:
                file_list.append(f)
                if f.endswith("_out.gtf") or f.endswith("merged.gtf"):
                    files = os.path.join(self.option('all_files_dir'), f)
                    r = open(files)
                    list1 = set("")
                    list2 = set("")
                    for line in r:
                        m = re.match("#.*", line)
                        if not m:
                            arr = line.strip().split("\t")
                            nine_array = arr[-1].strip().split(";")
                            gene_id = nine_array[0].strip().split("\"")
                            trans_id = nine_array[1].strip().split("\"")
                            if len(trans_id) == 3 and trans_id[1] not in list1:
                                list1.add(trans_id[1])
                            if len(gene_id) == 3 and gene_id[1] not in list2:
                                list2.add(gene_id[1])
                    num_count = f + "\t" + str(len(list1)) + "\t" + str(len(list2)) + "\n"
                    w.write(num_count)
                    r.close()

    def trans_stat(self):
        """
        统计组装后的序列的长度分布以及class_code分布情况
        :return:
        """
        # 避免tool反复运行时trans_count*统计文件反复追加
        code_num_filepath = os.path.join(self.work_dir, "code_num.txt")
        files = glob.glob(r"trans_count_stat_*.txt")
        for fp in files:
            os.remove(os.path.join(self.work_dir, fp))
        all_file = os.listdir(self.option('all_files_dir'))
        for f in all_file:
            # if f.endswith("merged.fa"):
            if f.endswith("all_transcripts.fa"):
                files = os.path.join(self.option('all_files_dir'), f)
                steps = [200, 300, 600, 1000]
                for step in steps:
                    step_count(files, self.option('all_files_dir') + "/" + f + ".txt", 10, step,
                               self.work_dir + "/trans_count_stat_" + str(step) + ".txt")
            if f.endswith("add_code_merged.gtf"):
                files = os.path.join(self.option('all_files_dir'), f)
                class_code_count(files, code_num_filepath)
                if len(open(code_num_filepath).readline().split('\t')) == 3:
                    self.logger.info("完成class_code统计")
                else:
                    self.set_error("class_code统计失败！", code="33703801")
        transcript_id_set = {r.id for r in
                             SeqIO.parse(os.path.join(self.option('all_files_dir'), 'all_transcripts.fa'), 'fasta')}
        transcript_id_to_code = dict()
        for line in open(code_num_filepath):
            eles = line.strip().split('\t')
            if len(eles) == 3:
                for transcript_id in eles[1].split(','):
                    transcript_id_to_code[transcript_id] = eles[0]
        code_to_transcript_id_set = {'=': set()}
        for transcript_id in transcript_id_set:
            if transcript_id in transcript_id_to_code:
                code = transcript_id_to_code[transcript_id]
                if code in code_to_transcript_id_set:
                    code_to_transcript_id_set[code].add(transcript_id)
                else:
                    code_to_transcript_id_set[code] = {transcript_id}
            else:
                code_to_transcript_id_set['='].add(transcript_id)
        with open(code_num_filepath, 'w') as fw:
            for code, id_set in code_to_transcript_id_set.items():
                fw.write('{}\t{}\t{}\n'.format(code, ','.join(id_set), len(id_set)))

    def count_genes_trans_exons(self):
        all_file = os.listdir(self.option('all_files_dir'))
        for f in all_file:
            if f.endswith("old_genes.gtf") or f.endswith("old_trans.gtf") or f.endswith(
                    "new_genes.gtf") or f.endswith("new_transcripts.gtf"):
                files = os.path.join(self.option('all_files_dir'), f)
                gene_trans = os.path.join(self.work_dir, f + ".trans")
                trans_exon = os.path.join(self.work_dir, f + ".exon")
                gene_trans_exon(files, self.option("assemble_method").lower(), gene_trans, trans_exon)
        gtf_files = os.listdir(self.work_dir)
        for f in gtf_files:
            if f.endswith('old_genes.gtf.trans') or f.endswith('new_genes.gtf.trans') or f.endswith(
                    'new_transcripts.gtf.exon') or f.endswith('old_trans.gtf.exon'):
                files = os.path.join(self.work_dir, f)
                steps = [1, 5, 10, 20]
                for step in steps:
                    middle_txt = os.path.join(self.work_dir, f + "_" + str(step) + ".test.txt")
                    final_txt = os.path.join(self.work_dir,
                                             f + "_" + str(step) + ".txt")
                    count_trans_or_exons(files, step, middle_txt, final_txt)

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        all_file = os.listdir(self.work_dir)
        for fp in all_file:
            if fp.endswith(".txt"):
                if fp != 'log.txt':
                    os.link(self.work_dir + '/' + fp, self.output_dir + "/" + fp)
