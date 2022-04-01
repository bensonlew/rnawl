# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.files.sequence.fastq import FastqFile
from mbio.files.sequence.file_sample import FileSampleFile
import os
import re
from BCBio import GFF
from gtfparse import read_gtf
from Bio import SeqIO
import pandas as pd

import unittest


class FileCheckAgent(Agent):
    """
    用于workflow开始之前对所输入的文件进行详细的内容检测
    """

    def __init__(self, parent):
        super(FileCheckAgent, self).__init__(parent)
        options = [
            {"name": "in_gtf", "type": "infile", "format": "ref_genome_db.gtf"},
            {"name": "gff", "type": "infile", "format": "ref_genome_db.gtf"},
            {"name": "genome", "type": "infile", "format": "ref_genome_db.common"},
            {"name": "bed", "type": "outfile", "format": "ref_genome_db.bed"},
            {"name": "out_gtf", "type": "outfile", "format": "ref_genome_db.gtf"},
            {"name": "cds", "type": "infile", "format": "ref_genome_db.common"},
            {"name": "pep", "type": "infile", "format": "ref_genome_db.common"},
            {"name": "g2t2p", "type": "infile", "format": "ref_genome_db.common"},
            {"name": "trans2desc", "type": "infile", "format": "ref_genome_db.common"},
            {"name": "trans2name", "type": "infile", "format": "ref_genome_db.common"},
            {"name": "go", "type": "infile", "format": "ref_genome_db.common"},
            {"name": "kegg", "type": "infile", "format": "ref_genome_db.common"},
            {"name": "biomart", "type": "infile", "format": "ref_genome_db.common"},
            {"name": "biomart_type", "type": "string", 'default': ""},
            {"name": "ncbi_gene", "type": "infile", 'format': "ref_genome_db.common"},
            {"name": "level_file", "type": "infile", 'format': "ref_genome_db.common"},
            {"name": "g2t2p_output", "type": "outfile", "format": "ref_genome_db.common"},
        ]
        self.add_option(options)
        self.step.add_steps("file_check")
        self.on('start', self.start_file_check)
        self.on('end', self.end_file_check)

    def start_file_check(self):
        self.step.file_check.start()
        self.step.update()

    def end_file_check(self):
        self.step.file_check.finish()
        self.step.update()

    def check_option(self):
        if not self.option("in_gtf").is_set and not self.option("gff").is_set:
            raise OptionError("GTF/GFF文件必须至少输入其中一个")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = "15G"


class FileCheckTool(Tool):
    """
    检查输入文件的格式是否符合要求
    """

    def __init__(self, config):
        super(FileCheckTool, self).__init__(config)
        self.gffread_path = "bioinfo/rna/cufflinks-2.2.1/gffread"

    def gff_to_gtf(self):
        self.logger.info("转换gff文件为gtf文件")
        gff = self.option("gff").prop["path"]
        tmp_gtf = os.path.join(self.work_dir, os.path.basename(gff) + ".tmp.gtf")
        gtf = os.path.join(self.work_dir, os.path.basename(gff) + ".gtf")
        genome = self.option("genome").prop["path"]
        to_gtf_cmd = '%s %s -g %s -T -O -o %s' % (self.gffread_path, gff, genome, tmp_gtf)
        command = self.add_command("gff_to_gtf", to_gtf_cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行gff_to_gtf完成")
        else:
            self.set_error("运行gff_to_gtf运行出错!")
            return False
        ## 过滤gtf文件，去掉不含gene_id、transcript_id、链信息的注释行
        with open(tmp_gtf, "r") as f1, open(gtf, "w") as w1:
            trans_pattern = re.compile(r'transcript_id \"([^\"]+)')
            trans_ids = dict()
            for line in f1:
                items = line.strip().split("\t")
                trans_id = trans_pattern.match(items[8]).group(1)
                if items[2] == 'exon' and trans_id not in trans_ids:
                    trans_ids[trans_id] = 1
                if len(items) != 9:
                    continue
                if not re.search(r'transcript_id', items[8]):
                    continue
                if not re.search(r'gene_id', items[8]):
                    continue
                if not re.search(r'[+|-]', items[6]):
                    continue
                if items[2] == 'CDS' and trans_id not in trans_ids:
                    w1.write(line.replace('CDS', 'exon'))
                if not line.strip().endswith(";"):
                    w1.write(line.strip() + ";\n")
                else:
                    w1.write(line)
        self.option("out_gtf", gtf)

    def check_input_file(self):
        ## check chromosome list
        gtf = self.option("out_gtf").prop["path"]
        genome = self.option("genome").prop["path"]
        genome_chr = list()
        for seq_record in SeqIO.parse(genome, "fasta"):
            if seq_record.id not in genome_chr:
                genome_chr.append(seq_record.id)
        df = read_gtf(gtf)
        if not (set(genome_chr) & set(df['seqname'])):
            self.set_error('GTF文件和参考基因组文件不匹配')
        if set(df['seqname']).difference(set(genome_chr)):
            self.set_error('{} 存在于gtf文件，但不存在参考基因组文件'.format(set(df['seqname']).difference(set(genome_chr))))
        if 'gene_id' not in df:
            self.set_error("{}文件中不存在gene_id注释信息，请检查文件是否规范！".format(os.path.basename(gtf)))
        if self.option("level_file").is_set:
            level = self.option("level_file").prop["path"]
            level_pd = pd.read_table(level, header=None)
            level_pd_chr = set(level_pd[0].astype("string"))
            if level_pd_chr.difference(set(genome_chr)):
                self.set_error(
                    "{}文件中存在错误的染色体名称：{}".format(os.path.basename(level), level_pd_chr.difference(set(genome_chr))))
            level_pd_lev = set(level_pd[1])
            for level_l in [string.lower() for string in level_pd_lev]:
                if level_l not in ["chromosome", "scaffold", "contig"]:
                    self.set_error(
                        "{}文件中存在不规范的组装水平：{}，只允许包含chromosome、scaffold、contig。".format(os.path.basename(level), level_l))
        gene_ids = set(df['gene_id'])
        if "" in gene_ids:
            gene_ids.remove("")
        if 'transcript_id' in df:
            transcript_ids = set(df['transcript_id'])
            if '' in transcript_ids:
                transcript_ids.remove("")
        else:
            transcript_ids = set()
        if 'protein_id' in df:
            protein_ids = set(df['protein_id'])
            if '' in protein_ids:
                protein_ids.remove("")
        else:
            protein_ids = set()
        ## check gene ids. transcript ids and protein ids
        if self.option("g2t2p").is_set:
            g2t2p = self.option("g2t2p").prop["path"]
            g2t2p_gene_ids  = set()
            g2t2p_transcript_ids = set()
            g2t2p_protein_ids = set()
            # 改写兼容第一行没有蛋白ID的表
            with open(g2t2p) as f:
                for line in f:
                    lines = line.strip().split("\t")
                    g2t2p_gene_ids.add(lines[0])
                    g2t2p_transcript_ids.add(lines[1])
                    if len(lines) > 2:
                        g2t2p_protein_ids.add(lines[2])

            # # g2t2p_pd = pd.read_table(g2t2p, header=None)
            # if len(lines) < 2:
            #     self.set_error('{}文件至少应包含两列，第一列为gene_id, 第二列为transcript_id！'.format(os.path.basename(g2t2p)))
            # g2t2p_gene_ids = set(g2t2p_pd[0])
            # if "" in g2t2p_gene_ids:
            #     g2t2p_gene_ids.remove("")
            # g2t2p_transcript_ids = set(g2t2p_pd[1])
            # if "" in g2t2p_transcript_ids:
            #     g2t2p_transcript_ids.remove("")
            # if len(g2t2p_pd.columns) == 3:
            #     g2t2p_protein_ids = set(g2t2p_pd[2])
            #     if "" in g2t2p_protein_ids:
            #         g2t2p_protein_ids.remove("")
            # else:
            #     g2t2p_protein_ids = set()
            if not gene_ids.intersection(g2t2p_gene_ids):
                self.set_error('{}文件与{}文件中gene id完全不一致，请检查！'.format(os.path.basename(gtf), os.path.basename(g2t2p)))
            if g2t2p_gene_ids.difference(gene_ids):
                self.logger.info(
                    "有{}个gene id仅存在于{}文件中.".format(len(g2t2p_gene_ids.difference(gene_ids)), os.path.basename(g2t2p)))
                # self.set_error('{}文件与{}文件中gene id不完全一致，如{}，请检查！'.format(os.path.basename(gtf), os.path.basename(g2t2p), g2t2p_gene_ids.difference(gene_ids)))
            if gene_ids.difference(g2t2p_gene_ids):
                self.logger.info(
                    "有{}个gene id仅存在于{}文件中.".format(len(gene_ids.difference(g2t2p_gene_ids)), os.path.basename(gtf)))
                self.set_error('{}文件与{}文件中gene id不完全一致，如{}，请检查！'.format(os.path.basename(gtf), os.path.basename(g2t2p),
                                                                        list(gene_ids.difference(g2t2p_gene_ids))[0]))
            if transcript_ids:
                if not transcript_ids.intersection(g2t2p_transcript_ids):
                    self.set_error(
                        '{}文件与{}文件中transcript id完全不一致，请检查！'.format(os.path.basename(gtf), os.path.basename(g2t2p)))
                if g2t2p_transcript_ids.difference(transcript_ids):
                    self.logger.info(
                        "有{}个transcript id仅存在于{}文件中.".format(len(g2t2p_gene_ids.difference(transcript_ids)),
                                                             os.path.basename(g2t2p)))
                if transcript_ids.difference(g2t2p_transcript_ids):
                    self.logger.info(
                        '{}文件与{}文件中transcript id不完全一致，如{}，请检查！'.format(os.path.basename(gtf), os.path.basename(g2t2p),
                                                                       list(transcript_ids.difference(g2t2p_transcript_ids))[0]))
                    self.set_error(
                        '{}文件与{}文件中transcript id不完全一致，如{}，请检查！'.format(os.path.basename(gtf), os.path.basename(g2t2p),
                                                                       list(transcript_ids.difference(g2t2p_transcript_ids))[0]))
            if protein_ids and g2t2p_protein_ids:
                if not transcript_ids.intersection(g2t2p_transcript_ids):
                    self.set_error(
                        '{}文件与{}文件中protein id完全不一致，请检查！'.format(os.path.basename(gtf), os.path.basename(g2t2p)))
                if g2t2p_protein_ids.difference(protein_ids):
                    self.logger.info("有{}个protein id仅存在于{}文件中.".format(len(g2t2p_protein_ids.difference(protein_ids)),
                                                                       os.path.basename(g2t2p)))
                if protein_ids.difference(g2t2p_protein_ids):
                    self.logger.info("有{}个protein id仅存在于{}文件中.".format(len(protein_ids.difference(g2t2p_protein_ids)),
                                                                       os.path.basename(gtf)))
                    self.set_error(
                        '{}文件与{}文件中protein id不完全一致，如{}，请检查！'.format(os.path.basename(gtf), os.path.basename(g2t2p),
                                                                    list(protein_ids.difference(g2t2p_protein_ids))[0]))
            with open(g2t2p, 'r') as f, open(os.path.join(self.output_dir, "g2t2p.txt"), 'w') as fw:
                for line in f:
                    lines = line.strip().split("\t")
                    if lines[0] in gene_ids:
                        fw.write(line)
            self.option("g2t2p_output").set_path(os.path.join(self.output_dir, "g2t2p.txt"))
        else:
            g2t2p_gene_ids = set()
            g2t2p_transcript_ids = set()
            g2t2p_protein_ids = set()
        total_gene_ids = gene_ids.union(g2t2p_gene_ids)
        total_transcript_ids = transcript_ids.union(g2t2p_transcript_ids)
        total_protein_ids = protein_ids.union(g2t2p_protein_ids)
        if self.option("cds").is_set:
            cds = self.option("cds").prop['path']
            for seq_record in SeqIO.parse(cds, "fasta"):
                cds_id = seq_record.id
                if total_transcript_ids:
                    if cds_id not in total_transcript_ids:
                        if '.' in cds_id:
                            cds_id_else = cds_id[:cds_id.rfind('.')]
                            if cds_id_else not in total_transcript_ids:
                                self.set_error("{}文件中存在错误的cds_id, 如{}，请核实后重新上传！".format(os.path.basename(cds), cds_id))
        if self.option("pep").is_set:
            pep = self.option("pep").prop['path']
            for seq_record in SeqIO.parse(pep, "fasta"):
                pep_id = seq_record.id
                if total_protein_ids:
                    if pep_id not in total_protein_ids:
                        if '.' in pep_id:
                            pep_id_else = pep_id[:pep_id.rfind('.')]
                            if pep_id_else not in total_protein_ids:
                                self.set_error(
                                    "{}文件中存在错误的protein_id, 如{}，请核实后重新上传！".format(os.path.basename(pep), pep_id))
        if self.option("trans2desc").is_set:
            trans2desc = self.option("trans2desc").prop['path']
            trans2desc_pd = pd.read_table(trans2desc, header=None)
            if len(trans2desc_pd.columns) < 2:
                self.set_error("{}文件至少需要包含两列，第一列为transcript_id，第二列为描述信息。".format(os.path.basename(trans2desc)))
            trans2desc_ids = set(trans2desc_pd[0])
            if "" in trans2desc_ids:
                trans2desc_ids.remove("")
            for id in trans2desc_ids:
                if total_transcript_ids:
                    if id not in total_transcript_ids:
                        self.set_error(
                            "{}文件中存在错误的transcript_id, 如{}，请核实后重新上传！".format(os.path.basename(trans2desc), id))
        if self.option("trans2name").is_set:
            trans2name = self.option("trans2name").prop['path']
            trans2name_pd = pd.read_table(trans2name, header=None)
            if len(trans2name_pd.columns) < 2:
                self.set_error("{}文件至少需要包含两列，第一列为transcript_id，第二列为gene name。".format(os.path.basename(trans2name)))
            trans2name_ids = set(trans2name_pd[0])
            if "" in trans2name_ids:
                trans2name_ids.remove("")
            for id in trans2name_ids:
                if total_transcript_ids:
                    if id not in total_transcript_ids:
                        self.set_error(
                            "{}文件中存在错误的transcript_id, 如{}，请核实后重新上传！".format(os.path.basename(trans2name), id))
        if self.option("biomart").is_set:
            biomart = self.option("biomart").prop['path']
            biomart_pd = pd.read_table(biomart, header=None)
            biomart_gene_ids = set(biomart_pd[0])
            if "" in biomart_gene_ids:
                biomart_gene_ids.remove("")
            biomart_trans_ids = set(biomart_pd[1])
            if "" in biomart_trans_ids:
                biomart_trans_ids.remove("")
            if self.option("biomart_type") == "type1":
                if len(biomart_pd.columns) != 18:
                    self.set_error("当biomart type为type1时，需要有18列相关信息")
                biomart_pep_ids = set(biomart_pd[6])
            elif self.option("biomart_type") == "type2":
                if len(biomart_pd.columns) != 16:
                    self.set_error("当biomart type为type1时，需要有16列相关信息")
                biomart_pep_ids = set(biomart_pd[4])
            elif self.option("biomart_type") == "type3":
                if len(biomart_pd.columns) != 14:
                    self.set_error("当biomart type为type1时，需要有14列相关信息")
                biomart_pep_ids = set(biomart_pd[2])
            else:
                biomart_pep_ids = set()
            if "" in biomart_pep_ids:
                biomart_pep_ids.remove("")
            if biomart_gene_ids:
                if not biomart_gene_ids.intersection(total_gene_ids):
                    self.set_error("{}文件中gene id与参考基因组注释文件不匹配，请核实后重新上传！".format(os.path.basename(biomart)))
            else:
                self.set_error("{}文件中缺少gene id列.".format(os.path.basename(biomart)))
            if biomart_trans_ids:
                if total_transcript_ids:
                    if not biomart_trans_ids.intersection(total_transcript_ids):
                        self.set_error("{}文件中transcript id与参考基因组注释文件不匹配，请核实后重新上传！".format(os.path.basename(biomart)))
            else:
                self.set_error("{}文件中缺少transcript id列.".format(os.path.basename(biomart)))
            if biomart_pep_ids:
                if total_protein_ids:
                    if not biomart_pep_ids.intersection(total_protein_ids):
                        self.set_error("{}文件中transcript id与参考基因组注释文件不匹配，请核实后重新上传！".format(os.path.basename(biomart)))
        if self.option("ncbi_gene").is_set:
            ncbi_gene = self.option("ncbi_gene").prop["path"]
            ncbi_gene_pd = pd.read_table(ncbi_gene, header=None)
            ncbi_gene_ids = set(ncbi_gene_pd[0])
            if "" in ncbi_gene_ids:
                ncbi_gene_ids.remove("")
            ncbi_transcript_ids = set(ncbi_gene_pd[1])
            if "" in ncbi_transcript_ids:
                ncbi_transcript_ids.remove("")
            if ncbi_gene_ids:
                if not ncbi_gene_ids.intersection(total_gene_ids):
                    self.set_error("{}文件中gene id与参考基因组注释文件不匹配，请核实后重新上传！".format(os.path.basename(ncbi_gene)))
            if ncbi_transcript_ids:
                if total_transcript_ids:
                    if not ncbi_transcript_ids.intersection(total_transcript_ids):
                        self.set_error("{}文件中transcript id与参考基因组注释文件不匹配，请核实后重新上传！".format(os.path.basename(ncbi_gene)))

    def gtf_to_bed(self):
        self.logger.info("转换gtf文件为bed文件")
        if self.option("in_gtf").is_set:
            new_gtf = self.work_dir + "/" + os.path.basename(self.option("in_gtf").prop["path"])
            in_gtf = self.option("in_gtf").prop["path"]
            # if os.path.exists(new_gtf):
            #     os.remove(new_gtf)
            # os.link(self.option("in_gtf").prop["path"], new_gtf)
            ## 过滤gtf文件，去掉不含gene_id、transcript_id、链信息的注释行
            with open(in_gtf, "r") as f1, open(new_gtf, "w") as w1:
                for line in f1:
                    items = line.strip().split("\t")
                    if len(items) != 9:
                        continue
                    if not re.search(r'transcript_id', items[8]):
                        continue
                    if not re.search(r'gene_id', items[8]):
                        continue
                    if not re.search(r'[+|-]', items[6]):
                        continue
                    if not line.strip().endswith(";"):
                        w1.write(line.strip() + ";\n")
                    else:
                        w1.write(line)
            self.option("out_gtf").set_path(new_gtf)
        self.option("out_gtf").to_bed()

    def make_true_bed(self):
        '''
        modified by qinjincheng at 20190402
        create a legal bed by ucsc utilities
        '''
        gtf = self.option('out_gtf').prop['path']
        genepred = '{}.genepred'.format(self.option('out_gtf').prop['path'])
        bed = '{}.bed'.format(self.option('out_gtf').prop['path'])
        gtfToGenePred = 'bioinfo/align/ucsc_tools/gtfToGenePred'
        genePredToBed = 'bioinfo/align/ucsc_tools/genePredToBed'
        cmd1 = '{} -ignoreGroupsWithoutExons {} {}'.format(gtfToGenePred, gtf, genepred)
        self.run_code('gtftogenepred', cmd1)
        cmd2 = '{} {} {}'.format(genePredToBed, genepred, bed)
        self.run_code('genepredtobed', cmd2)

    def run_code(self, cmd_name, cmd, shell=False):
        if shell:
            cmd = self.config.SOFTWARE_DIR + '/' + cmd
        command = self.add_command(cmd_name, cmd, shell=shell)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('succeed in running {}'.format(cmd_name))
        elif command.return_code is None:
            self.logger.warn('fail to run {}, try again'.format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info('succeed in rerunning {}'.format(cmd_name))
            else:
                self.set_error('fail to rerun {}, abord'.format(cmd_name))
        else:
            self.set_error('fail to run {}, abord'.format(cmd_name))

    def filter_bed(self):
        ## modified by shicaiping at 20181225
        ## 新增对bed文件的过滤，基于以下两点：
        # ① 基因单条染色体序列长度超过500000000，无法进行Bamdistribution分析，过滤超过500000000部分的转录本
        # ② 过滤exon存在交叉、包含等关系的转录本
        with open(self.option("out_gtf").prop["path"] + ".bed", "r") as f, open(
                self.option("out_gtf").prop["path"] + ".filter.bed", "w") as w:
            for line in f:
                if int(line.split("\t")[2]) < 500000000:
                    flag = 1
                    start_list = line.strip().split("\t")[11].split(",")[:-1]
                    length_list = line.strip().split("\t")[10].split(",")[:-1]
                    for index, value in enumerate(length_list[:-1]):
                        if ((int(length_list[index]) + int(start_list[index])) >= int(start_list[index + 1])):
                            flag = 0
                            break
                    if flag == 1:
                        w.write(line)
        self.option("bed").set_path(self.option("out_gtf").prop["path"] + ".filter.bed")

    def set_output(self):
        bed = self.option("out_gtf").prop["path"] + ".filter.bed"
        gtf = self.option("out_gtf").prop["path"]
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(bed))):
            os.remove(os.path.join(self.output_dir, os.path.basename(bed)))
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(gtf))):
            os.remove(os.path.join(self.output_dir, os.path.basename(gtf)))
        os.link(bed, os.path.join(self.output_dir, os.path.basename(bed)))
        os.link(gtf, os.path.join(self.output_dir, os.path.basename(gtf)))

    def run(self):
        super(FileCheckTool, self).run()
        if not self.option("in_gtf").is_set:
            self.gff_to_gtf()
        self.gtf_to_bed()
        self.make_true_bed()
        self.filter_bed()
        self.check_input_file()
        self.set_output()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "FileCheck_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_genome_db_v2.file_check",
            "instant": False,
            "options": dict(
                in_gtf="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Mus_musculus/ensembl/Mus_musculus.GRCm38.96.gtf",
                genome="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Mus_musculus/ensembl/Mus_musculus.GRCm38.dna.toplevel.fa",
                g2t2p="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Mus_musculus/ensembl/g2t2p.txt",
                cds="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Mus_musculus/ensembl/Mus_musculus.GRCm38.cds.all.fa",
                pep="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Mus_musculus/ensembl/Mus_musculus.GRCm38.pep.all.fa",
                go="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Mus_musculus/ensembl/go.txt",
                biomart="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Homo_sapiens/ensembl/biomart.txt",
                biomart_type="type1",
                ncbi_gene="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Mus_musculus/ensembl/ncbi_entrez.txt"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
