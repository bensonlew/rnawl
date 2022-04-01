# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.ref_rna.trans_step import step_count, class_code_count, gene_trans_exon, count_trans_or_exons
import re
import glob

class RefassembleStatAgent(Agent):
    """
    有参组装数据统计
    version v1.0.1
    author: wangzhaoyue
    last_modify: 2016.09.09
    """
    def __init__(self, parent):
        super(RefassembleStatAgent, self).__init__(parent)
        options = [
            {"name": "all_files_dir", "type": "string"},  # 输入的文件夹
            {"name": "assemble_method", "type": "string", "default": "cufflinks"},  # 选择拼接软件
            # {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},  # 参考基因文件
            # {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},  # 参考基因的注释文件
            # {"name": "cpu", "type": "int", "default": 10},  #RefassembleStat软件所分配的cpu数量
            # {"name": "F", "type": "int", "default": 0.1},  # min-isoform-fraction
            # {"name": "fr_stranded", "type": "string", "default": "fr-unstranded"},  # 是否链特异性
            # {"name": "strand_direct", "type": "string", "default": "none"},  # 链特异性时选择正负链
            # {"name": "sample_gtf", "type": "outfile", "format": "gene_structure.gtf"},  # 输出的转录本文件
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
        """
        重写参数检测函数
        :return:
        """
        # if not self.option('sample_bam'):
        #     raise OptionError('必须输入样本文件为bam格式')
        # if not self.option('ref_fa'):
        #     raise OptionError('必须输入参考序列ref.fa')
        # if not self.option('ref_gtf'):
        #     raise OptionError('必须输入参考序列ref.gtf')
        # if self.option("fr_stranded") != "fr-unstranded" and not self.option("strand_direct").is_set:
        #     raise OptionError("当链特异性时必须选择正负链")
        # return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        # result_dir.add_regexp_rules([
        #     ["_out.gtf", "gtf", "样本拼接之后的gtf文件"]
        # ])
        super(RefassembleStatAgent, self).end()


class RefassembleStatTool(Tool):
    def __init__(self, config):
        super(RefassembleStatTool, self).__init__(config)
        self._version = "v1.0.1"
        # self.RefassembleStat_path = '/bioinfo/rna/RefassembleStat-2.2.1/'

    def run(self):
        """
        运行
        :return:
        """
        super(RefassembleStatTool, self).run()
        self.get_numberlist()
        self.trans_stat()
        self.count_genes_trans_exons()
        # self.run_RefassembleStat()
        # self.run_gffread()
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
        files = glob.glob(r"trans_count_stat_*.txt")
        for file in files:
            os.remove(os.path.join(self.work_dir, file)) ## modified by shicaiping at 20180726
        all_file = os.listdir(self.option('all_files_dir'))
        for f in all_file:
            if f.endswith("merged.fa"):
                files = os.path.join(self.option('all_files_dir'), f)
                steps = [200, 300, 600, 1000]
                for step in steps:
                    step_count(files, self.option('all_files_dir') + "/" + f + ".txt", 10, step,
                               self.work_dir + "/trans_count_stat_" + str(step) + ".txt")
            # if f.endswith("change_id_merged.gtf"):
            if f.endswith("add_code_merged.gtf"):
                files = os.path.join(self.option('all_files_dir'), f)
                code_count = os.path.join(self.work_dir, "code_num.txt")
                class_code_count(files, code_count)
                if len(open(code_count).readline().split('\t')) == 3:
                    self.logger.info("完成class_code统计")
                else:
                    raise Exception("class_code统计失败！")

    def count_genes_trans_exons(self):
        all_file = os.listdir(self.option('all_files_dir'))
        for f in all_file:
            if f.endswith("old_genes.gtf") or f.endswith("old_trans.gtf") or f.endswith(
                    "new_genes.gtf") or f.endswith("new_transcripts.gtf"):
                files = os.path.join(self.option('all_files_dir'), f)
                gene_trans = os.path.join(self.work_dir, f + ".trans")
                trans_exon = os.path.join(self.work_dir, f + ".exon")
                gene_trans_exon(files, self.option("assemble_method"), gene_trans, trans_exon)
                # count_trans_or_exons(gene_trans, final_trans_file)
                # count_trans_or_exons(trans_exon, final_exon_file)
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
        try:
            all_file = os.listdir(self.work_dir)
            for files in all_file:
                if files.endswith(".txt"):
                    if files != 'log.txt':
                        os.link(self.work_dir + '/' + files, self.output_dir + "/" + files)
            self.logger.info("设置组装拼接分析结果目录成功")

        except Exception as e:
            self.logger.info("设置组装拼接分析结果目录失败{}".format(e))
            self.set_error("设置组装拼接分析结果目录失败{}".format(e))
