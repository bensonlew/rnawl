# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2019/5/7'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import pandas as pd
from biocluster.core.exceptions import OptionError
from Bio import SeqIO
from mbio.packages.metagenomic.common import link_file


class GapFaAgent(Agent):
    """
    DemoAgent:
    version 1.0
    """

    def __init__(self, parent):
        super(GapFaAgent, self).__init__(parent)
        options = [
            {"name": "circle_fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "chr_fa", "type": "infile", "format": "sequence.fasta", "required": True},
            {"name": "fill_fa", "type": "infile", "format": "sequence.fasta", "required": True},
            {"name": "pre_log", "type": "infile", "format": "sequence.profile_table"},  # 可能没有文件
            # chromsome/plasmid ID  seq_id  seq_ids_contained_in_this_seq_id(equal seq_id at first time running this tool)  is_circle or not
            # seq_id可作为索引
            {"name": "overlap", "type": "int", "default": 1000},
            {"name": "identity", "type": "float", "default": 0.95},
            {"name": "out_fa", "type": "outfile", "format": "sequence.fasta"},
            {"name": "log", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = "10G"


class GapFaTool(Tool):
    def __init__(self, config):
        super(GapFaTool, self).__init__(config)
        self.all_circled = False
        self.blast_path = "bioinfo/align/ncbi-blast-2.3.0+/bin"
        self.dnaA = self.config.SOFTWARE_DIR + "/database/bacgenome/dnaA/dnaA.faa"
        self.blast_table = os.path.join(self.work_dir, "blast.xls")
        self.script = self.config.PACKAGE_DIR + '/bacgenome/gap_fill2.py'
        self.python_path = '/program/Python/bin/python'

    def make_db(self):
        # 骨架序列建库
        # 将参考序列建库
        # 去除低于overlap长度的序列
        cmd = "%s/makeblastdb -dbtype nucl -in %s -parse_seqids -out chr_ref" % (self.blast_path, self.option("chr_fa").prop["path"])
        makedb = self.add_command("mkdb", cmd).run()
        self.wait(makedb)
        if makedb.return_code == 0:
            self.logger.info("构建参考库成功")
        else:
            self.set_error("构建参考库失败")

    def blast(self):
        # blast比对
        cmd = "%s/blastn -query %s -db chr_ref -out %s -evalue 1e-5 -outfmt 6 -max_hsps 10 -num_threads 2" % (self.blast_path, self.option("fill_fa").prop["path"], self.blast_table)
        command = self.add_command("blast", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("blast比对成功")
        else:
            self.set_error("blast比对失败")

    def gap_fill(self):
        cmd = "%s %s -fa %s -query %s -blast %s -log %s -aln_len %s -identity %s -o_prefix %s" % (self.python_path,
        self.script, self.option("chr_fa").prop["path"], self.option("fill_fa").prop["path"], self.blast_table, self.option("pre_log").prop["path"],
         self.option("overlap"), self.option("identity"), self.work_dir + "/tmp_result")
        command = self.add_command("gap_fill", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gapfill成功")
        else:
            self.set_error("gapfill失败")

    def get_dnaA_loc(self, blast_out):
        if os.path.getsize(blast_out) == 0:
            return None, 0
        data = pd.read_table(blast_out, header=None)
        direction = 1 if data.loc[0, 6] < data.loc[0,7] else -1
        return data.loc[0, 6], direction

    def run_dnaA(self, query, out):
        out = str(out)
        cmd = os.path.join(self.blast_path, "blastx") + " -query %s -db %s -out %s -evalue 1e-5 -outfmt 6 -max_hsps 1 -num_threads 2 -max_target_seqs 1" % (query, self.dnaA, out)
        name = os.path.basename(out).lower()
        command = self.add_command("dnaa%s" % name, cmd).run()
        self.wait(command)
        if command.return_code == 0:
            return self.get_dnaA_loc(out)
        else:
            self.set_error("dnaA比对出错")

    def start_loc(self):
        predict_dir = os.path.join(self.work_dir, "predict_dir")
        if not os.path.isdir(predict_dir):
            os.mkdir(predict_dir)
        new_fasta = []
        fasta_path = os.path.join(self.work_dir, "tmp_result.fa")
        self.out = self.work_dir + "/last_result.fa"
        if os.path.exists(self.out):
            os.remove(self.out)
        if self.option("circle_fa").is_set:
            os.system("cat {} {} >{}".format(self.option("circle_fa").prop['path'], fasta_path, self.out))
        else:
             os.link(self.work_dir + "/tmp_result.fa", self.out)
        new_fasta_path = os.path.join(self.output_dir, "result.fa")
        log = pd.read_table(self.work_dir + "/tmp_result.log", header=None, index_col=1)
        log.columns = ["genome", "parent", "is_circle"]
        for seq_record in SeqIO.parse(self.work_dir + "/last_result.fa", "fasta"):
            seq_id = seq_record.id
            if seq_id in list(log['genome']):
                new_fasta.append(seq_record)
            elif log.loc[seq_id, "is_circle"]:
                self.logger.info("check log is true or false: %s" % log.loc[seq_id, "is_circle"])
                self.logger.info(seq_id)
                query = os.path.join(predict_dir, seq_id + ".fasta")
                SeqIO.write(seq_record, query, "fasta")
                dnaA_start, direction = self.run_dnaA(query, query + ".out")  # 增加比对方向，反向比对结果进行反向互补处理
                if direction == 1:
                    # dnaA起始位置前需要保留100个碱基
                    seq_record = seq_record[dnaA_start-101:] + seq_record[:dnaA_start-101]  # 如果起始点少于101，可以向序列末端借碱基，代码是一样的
                elif direction == -1:
                    # 反向互补需要对Seq进行处理，不能是Seq_record
                    if len(seq_record) - dnaA_start >= 100:
                        seq_record.seq = seq_record.seq[:dnaA_start+100].reverse_complement() + seq_record.seq[dnaA_start+100:].reverse_complement()
                    else:
                        seq_need_len = dnaA_start + 100 - len(seq_record)  # 末端不够100bp作为起始序列，到序列前端借
                        seq_record.seq = seq_record.seq[:seq_need_len].reverse_complement() + seq_record.seq[seq_need_len:].reverse_complement()
                new_fasta.append(seq_record)
            else:
                new_fasta.append(seq_record)
        SeqIO.write(new_fasta, new_fasta_path, "fasta")
        link_file(self.work_dir + "/tmp_result.log", self.output_dir + "/result.log")
        link_file(self.work_dir + "/blast.xls", self.output_dir + "/blast.xls")

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        self.option("out_fa").set_path(self.output_dir + "/result.fa")
        self.option("log").set_path(self.output_dir + "/result.log")

    def run(self):
        super(GapFaTool, self).run()
        self.make_db()
        self.blast()
        self.gap_fill()
        self.start_loc()
        self.set_output()
        self.end()