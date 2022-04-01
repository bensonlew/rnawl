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
from Bio.SeqIO.FastaIO import FastaWriter
import pandas as pd
from pandas.api.types import CategoricalDtype
import unittest


class CheckGenomeAgent(Agent):
    """
    用于workflow开始之前对所输入的参考基因组及其注释文件进行详细的内容检测
    """

    def __init__(self, parent):
        super(CheckGenomeAgent, self).__init__(parent)
        options = [
            {"name": "in_type", "type": "string", "default": ""}, # gtf or gff
            {"name": "in_file", "type": "string", "default": ""},
            {"name": "genome", "type": "string", "default": ""},
            {"name": "gffread", "type": "string", "default": "false"},      # switch off filtering by gffread
            {"name": "out_file", "type": "outfile", "format": "prok_rna.common"},
        ]
        self.add_option(options)
        self.step.add_steps("check_genome")
        self.on('start', self.start_file_check)
        self.on('end', self.end_file_check)

    def start_file_check(self):
        self.step.check_genome.start()
        self.step.update()

    def end_file_check(self):
        self.step.check_genome.finish()
        self.step.update()

    def check_option(self):
        if self.option("in_file") is None or not os.path.exists(self.option("in_file")):
            raise OptionError("gff/gtf文件不存在")
        pass

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = "15G"


class CheckGenomeTool(Tool):
    """
    检查输入文件的格式是否符合要求
    """

    def __init__(self, config):
        super(CheckGenomeTool, self).__init__(config)
        self.gffread_path = "bioinfo/rna/cufflinks-2.2.1/gffread"

    def check_genome(self):
        self.logger.info("开始检查参考基因组文件")
        genome = self.option("genome")
        genome_chr = list()
        genome_range = os.path.join(self.work_dir, "genome_range.txt")
        with open(genome_range, "w") as w:
            for seq_record in SeqIO.parse(genome, "fasta"):
                if seq_record.id not in genome_chr:
                    genome_chr.append(seq_record.id)
                    w.write(seq_record.id + ":0.." + str(len(seq_record.seq)) + "\n")
        # 过滤参考基因组fa文件空行，并且按照固定长度输出
        filter_fa = os.path.join(self.work_dir, 'ref_filter.fa')
        output_fasta = open(filter_fa, "w")
        writer = FastaWriter(output_fasta, wrap=70)
        writer.write_file(SeqIO.parse(genome, "fasta"))
        self.option("genome",filter_fa)
        with open(self.option("in_file"), "r") as f:
            for line in f:
                if line.startswith("#") or not len(line.strip()):
                    continue
                if line.strip().split()[0] not in genome_chr:
                    self.set_error('{}文件和{}不匹配, {}在{}中不存在'.format(os.path.basename(self.option("in_file")),
                                                                  os.path.basename(genome), line.strip().split()[0],
                                                                  os.path.basename(genome)))

    def filter_gtf(self):
        self.logger.info("开始过滤gtf文件")
        gtf = self.option("in_file")
        filter_gtf = os.path.join(self.work_dir, os.path.basename(gtf) + ".filter.gtf")
        genome = self.option("genome")
        if self.option('gffread') == 'true':
            if os.path.exists(genome):
                to_gtf_cmd = '%s %s -g %s -T -r %s -R -F -o %s' % (
                self.gffread_path, gtf, genome, os.path.join(self.work_dir, "genome_range.txt"), filter_gtf)
            else:
                to_gtf_cmd = '%s %s -T -r %s -R -F -o %s' % (
                self.gffread_path, gtf, os.path.join(self.work_dir, "genome_range.txt"), filter_gtf)
            command = self.add_command("filter_gtf", to_gtf_cmd, ignore_error=True)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("运行filter_gtf完成")
            else:
                error_file = os.path.join(self.work_dir, 'filter_gtf.o')
                if os.path.exists(error_file):
                    self.logger.info("解析filter_gtf.o")
                    error_flag = False
                    with open(error_file, "r") as f:
                        for line in f:
                            if 'error' in line or 'Error' in line:
                                error_flag = True
                                self.set_error(line.strip())
                    if not error_flag:
                        self.set_error("gffread过滤GTF文件报错，请核实格式是否符合规范!")
                else:
                    self.set_error("gffread过滤GTF文件报错，请核实格式是否符合规范!")
            # 当gtf文件不符合规范过滤后为空时，放弃过滤
            if os.path.getsize(filter_gtf) == 0 or float(os.path.getsize(filter_gtf))/float(os.path.getsize(gtf)) <= 0.8:
                os.system('cp %s %s' % (gtf, filter_gtf))
        if self.option('gffread') == 'false':
            os.system('cp %s %s' % (gtf, filter_gtf))
        # 替换ID中空格为‘-’
        target_ids = dict()
        with open(filter_gtf, "r") as f:
            for line in f:
                if line.startswith("#") or len(line.strip().split("\t")) < 9:
                    continue
                else:
                    features = line.strip().split("\t")[8].split(";")
                    for feature in features:
                        if '_id \"' in feature and ' ' in feature.split("_id ")[1]:
                            target_ids[feature.strip().split('_id ')[1]] = feature.strip().split('_id ')[1].replace(" ", "-")
        for id in target_ids:
            self.logger.info("替换{}为{}".format(id, target_ids[id]))
            os.system("sed -i 's/%s/%s/g' %s " % (id, target_ids[id], filter_gtf))

    def filter_gff(self):
        self.logger.info("开始过滤gff文件")
        gff = self.option("in_file")
        filter_gff = os.path.join(self.work_dir, os.path.basename(gff) + ".filter.gff")
        genome = self.option("genome")
        if self.option('gffread') == 'true':
            if os.path.exists(genome):
                to_gff_cmd = '%s %s -g %s -r %s -R -F -o %s' % (
                self.gffread_path, gff, genome, os.path.join(self.work_dir, "genome_range.txt"), filter_gff)
            else:
                to_gff_cmd = '%s %s -r %s -R -F -o %s' % (
                self.gffread_path, gff, os.path.join(self.work_dir, "genome_range.txt"), filter_gff)
            command = self.add_command("filter_gff", to_gff_cmd, ignore_error=True)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("运行filter_gff完成")
            else:
                error_file = os.path.join(self.work_dir, 'filter_gff.o')
                if os.path.exists(error_file):
                    self.logger.info("解析filter_gff.o")
                    error_flag = False
                    with open(error_file, "r") as f:
                        for line in f:
                            if 'error' in line or 'Error' in line:
                                self.logger.info("找到错误行")
                                if 'duplicate/invalid \'tRNA\' feature' in line:
                                    self.logger.info("兼容内部产生错误GFF格式文件")
                                    error_flag = True
                                    pass
                                else:
                                    error_flag = True
                                    self.set_error(line.strip())
                    if not error_flag:
                        self.set_error("gffread过滤GFF文件报错，请核实格式是否符合规范!")
                else:
                    self.set_error("gffread过滤GFF文件报错，请核实格式是否符合规范!")
            # 当gff文件不符合规范过滤后为空时，放弃过滤
            if os.path.getsize(filter_gff) == 0 or float(os.path.getsize(filter_gff))/float(os.path.getsize(gff)) <= 0.8:
                os.system('cp %s %s' % (gff, filter_gff))
        if self.option('gffread') == 'false':
            os.system('cp %s %s' % (gff, filter_gff))
        # 替换ID中空格为‘-’
        target_ids = dict()
        chrom_list = list()
        processed_gff = os.path.join(self.work_dir, os.path.basename(gff) + ".processed.gff")
        with open(filter_gff, "r") as f, open(processed_gff, 'w') as p:
            line2skip = 0
            line_index = 0
            fasta_seq = list()
            for line in f:
                line_index += 1
                if line.startswith('##FASTA'):  # dealing with gff contains fasta sequence
                    fasta_seq.append(line.strip())
                elif line.startswith("#") or len(line.strip().split("\t")) < 9:
                    if line_index == line2skip + 1:
                        line2skip += 1
                    p.write(line)
                elif not line.startswith("#") and len(line.strip().split("\t")) >=9:
                    if line.strip().split("\t")[0] not in chrom_list:
                        chrom_list.append(line.strip().split("\t")[0])
                    features = line.strip().split("\t")[8].split(";")
                    for feature in features:
                        if 'ID=' in feature and ' ' in feature:
                            target_ids[feature.strip().split('=')[1]] = feature.strip().split('=')[1].replace(" ", "-")
                elif line.strip():
                    fasta_seq.append(line.strip())

        for id in target_ids:
            self.logger.info("替换{}为{}".format(id, target_ids[id]))
            os.system("sed -i 's/%s/%s/g' %s " % (id, target_ids[id], filter_gff))

        self.gff_processing(filter_gff, processed_gff, line2skip, chrom_list, fasta_seq)

    def gff_processing(self, filter_gff, processed_gff, line2skip=0, chrom_list=None, fasta_seq=None):
        # Sorting gff by gene length, before removing any ID-duplicated records
        if line2skip == 0:
            line2skip = None
        gff_pd = pd.read_table(filter_gff, header=None, skiprows=line2skip, sep='\t', dtype={0: 'str'})
        gff_pd.rename(columns={0: 'chrom', 1: 'source', 2: 'type', 3: 'start', 4: 'end', 5: 'score', 6: 'strand',
                               7: 'phase', 8: 'attributes'}, inplace=True)
        gff_pd['length'] = gff_pd.apply(lambda x: abs(x['end']-x['start']), axis=1)
        if chrom_list:
            cat_order = CategoricalDtype(chrom_list, ordered=True)
            gff_pd['chrom'] = gff_pd['chrom'].astype(cat_order)
        gff_pd.sort_values(by=['chrom', 'start', 'length'], ascending=[True, True, False], inplace=True)
        gff_pd.drop('length', axis=1, inplace=True)
        gff_pd.dropna(axis=0, how='any', inplace=True)
        locus2id = dict()       # replace any duplicated locus_tag for multiple gene IDs
        id2locus = dict()
        cds_list = list()
        gene_list = list()
        removed_gff = os.path.join(self.work_dir, os.path.basename(self.option('in_file')) + ".removed.gff")
        # ID-locus tag modification
        with open(processed_gff, 'a') as p, open(removed_gff, 'w') as r:
            for row in gff_pd.iterrows():
                retain = True
                items = [str(int(i)) if type(i) == float else str(i) for i in row[1].values.tolist()]
                attr = items[8]
                if items[6] not in ['+', '-']:  # Improper strand
                    items[6] = '+'
                if not attr.endswith(';'):
                    attr += ';'
                locus = re.match(r'.*?locus_tag=(.*?);.*', attr)
                if items[2].lower() == 'gene':
                    id_re = re.match(r'.*?ID=(.*?);.*', attr)
                    if id_re and id_re.group(1) not in gene_list:   # remove records with duplicated ID
                        gene_list.append(id_re.group(1))
                    else:
                        retain = False
                elif items[2].lower() == 'cds':
                    id_re = re.match(r'.*?Parent=(.*?);.*', attr)
                    if not id_re:
                        id_re = re.match(r'.*?ID=(.*?);.*', attr)
                    product_re = re.match(r'.*?product=(.*?);.*', attr)
                    if id_re and id_re.group(1) not in cds_list:
                        cds_list.append(id_re.group(1))
                    else:
                        retain = False
                    if product_re:
                        product = product_re.group(1).replace('=', '-')  # avoid multiple '=' in product
                        items[8] = items[8].replace('product=' + product_re.group(1), 'product=' + product)
                else:
                    id_re = None
                if locus and id_re:
                    locus_tag = locus.group(1)
                    gene_id = id_re.group(1)
                    if gene_id in id2locus.keys() and locus_tag != id2locus[gene_id][0]:  # ID with modified locus_tag
                        items[8] = items[8].replace('locus_tag=' + locus_tag, 'locus_tag=' + id2locus[gene_id])
                    elif gene_id not in id2locus.keys() and locus_tag in locus2id.keys():  # Duplicate locus_tag
                        if len(gene_id.split(locus_tag)) > 1:
                            locus_processed = locus_tag + gene_id.split(locus_tag)[1]  # remove potential prefix
                        else:
                            locus_processed = gene_id
                        id2locus[gene_id] = locus_processed
                        locus2id[locus_processed] = gene_id
                        items[8] = items[8].replace('locus_tag=' + locus_tag, 'locus_tag=' + locus_processed)
                    elif gene_id not in id2locus.keys() and locus_tag not in locus2id.keys():  # new locus_tag
                        id2locus[gene_id] = locus_tag
                        locus2id[locus_tag] = gene_id
                if retain:
                    p.write('\t'.join([str(i) for i in items]) + '\n')
                else:
                    r.write('\t'.join([str(i) for i in items]) + '\n')
            if fasta_seq:
                r.write('\n'.join(fasta_seq) + '\n')

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


    def set_output(self):
        if self.option("in_type").lower() == "gtf":
            filter_gtf = os.path.join(self.work_dir, os.path.basename(self.option("in_file")) + ".filter.gtf")
            if os.path.exists(os.path.join(self.output_dir, os.path.basename(filter_gtf))):
                os.remove(os.path.join(self.output_dir, os.path.basename(filter_gtf)))
            os.link(filter_gtf, os.path.join(self.output_dir, os.path.basename(filter_gtf)))
            self.option("out_file", os.path.join(self.output_dir, os.path.basename(filter_gtf)))
        else:
            # filter_gff = os.path.join(self.work_dir, os.path.basename(self.option("in_file")) + ".filter.gff")
            filter_gff = os.path.join(self.work_dir, os.path.basename(self.option("in_file")) + ".processed.gff")
            if os.path.exists(os.path.join(self.output_dir, os.path.basename(filter_gff))):
                os.remove(os.path.join(self.output_dir, os.path.basename(filter_gff)))
            os.link(filter_gff, os.path.join(self.output_dir, os.path.basename(filter_gff)))
            self.option("out_file", os.path.join(self.output_dir, os.path.basename(filter_gff)))


    def run(self):
        super(CheckGenomeTool, self).run()
        if os.path.exists(self.option("genome")):
            self.check_genome()
        if self.option("in_type").lower() == "gtf":
            self.filter_gtf()
        else:
            self.filter_gff()
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
            "id": "CheckGenome_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "prok_rna.check_genome",
            "instant": False,
            "options": dict(
                in_type="gtf",
                in_file="/mnt/ilustre/users/sanger-dev/wpm2/workspace/20210610/Prokrna_gbfl_1pln6ksbdhrmsc3qrd1u4m/Download/output/GCF_000009345.1_ASM934v1/GCF_000009345.1_ASM934v1_genomic.gtf",
                genome="/mnt/ilustre/users/sanger-dev/wpm2/workspace/20210610/Prokrna_gbfl_1pln6ksbdhrmsc3qrd1u4m/Download/output/GCF_000009345.1_ASM934v1/GCF_000009345.1_ASM934v1_genomic.fna",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
