# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2019/5/9'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import pandas as pd
from Bio import SeqIO


class HgeneBlastxAgent(Agent):
    """
    DemoAgent:
    version 1.0
    """

    def __init__(self, parent):
        super(HgeneBlastxAgent, self).__init__(parent)
        options = [
            {"name": "seq_fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "out_table", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)


    def check_options(self):
        """
        检查参数是否正确
        """
        pass

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 4
        self._memory = "20G"


class HgeneBlastxTool(Tool):
    def __init__(self, config):
        super(HgeneBlastxTool, self).__init__(config)
        self.diamond = "/bioinfo/align/diamond-0.8.35/diamond"
        self.core_gene = self.config.SOFTWARE_DIR + "/database/bacgenome/House-keeping_gene/Core_gene"
        self.all_core_gene = self.config.SOFTWARE_DIR + "/database/bacgenome/House-keeping_gene/all_Coregene"


    def run_get_gene(self):
        """
        description
        :return:
        """
        cmd = "{} blastx -d {} -q {} -o {}".format(self.diamond, self.core_gene, self.option("seq_fa").prop["path"], self.work_dir + "/matches.m8")
        command = self.add_command("run_get_gene", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_get_gene完成")
        else:
            self.set_error("运行run_get_gene出错！")

    def run_get_seq(self):
        data = pd.read_table(self.work_dir + "/matches.m8", header=None, names=["query", "ref", "ident", "aln", "mis", "gap", "qstart", "qend", "rstart", "rend", "evalue", "score"])
        hgene = data["ref"].tolist()
        hgene_set = set(hgene)
        seq_list = []
        seq_records = SeqIO.parse(self.option("seq_fa").prop["path"], "fasta")
        seq_dict = SeqIO.to_dict(seq_records)
        for gene in hgene_set:
            gene_blast = data.loc[data["ref"]==gene,]
            best_one = gene_blast.sort_values(["aln", "score"], ascending=False).iloc[0,]
            query = best_one["query"]
            qstart = best_one["qstart"]
            qend = best_one["qend"]
            query_seq = seq_dict[query]
            if qstart < qend:
                query_seq = query_seq[qstart-1:qend]
            else:
                query_seq = query_seq[qstart:qend-1:-1]
            query_seq.id = gene
            seq_list.append(query_seq)
        SeqIO.write(seq_list, os.path.join(self.work_dir, "find_gene.fa"), "fasta")

    def run_get_all(self):
        cmd = "{} blastx -d {} -q {} -o {} --max-target-seqs 10 -f 6 qseqid sseqid pident qlen mismatch gapopen qstart qend sstart send evalue bitscore".format(self.diamond, self.all_core_gene, os.path.join(self.work_dir, "find_gene.fa"), self.work_dir + "/matches2.m8")
        command = self.add_command("get_all", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_get_all完成")
        else:
            self.set_error("运行run_get_all出错!")

    def process_result(self):
        data = pd.read_table(self.work_dir + "/matches2.m8", header=None, names=["query", "ref", "Identity (%)", "qlen", "mis", "gap", "qstart", "qend", "rstart", "rend", "Evalue", "Score"])
        data["House-keeping gene"] = data.apply(lambda x: x["query"].split(" ")[0], axis=1)
        data["Organism"] = data["ref"]
        data['Coverage (%)'] = (data["qend"]-data["qstart"]) / data['qlen'] * 100
        data.reindex(columns=["House-keeping gene", "Organism", "Identity (%)", "Coverage (%)", "Evalue", "Score"]).\
            to_csv(self.output_dir + "/hgene_blast.xls", sep="\t", header=True, index=False)

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        self.option("out_table").set_path(self.output_dir + "/hgene_blast.xls")

    def run(self):
        super(HgeneBlastxTool, self).run()
        self.run_get_gene()
        self.run_get_seq()
        self.run_get_all()
        try:
            self.process_result()
            self.set_output()
        except Exception,e:
            self.logger.error(e)
            self.logger.error("没有比对到结果，正常结束")
        self.end()