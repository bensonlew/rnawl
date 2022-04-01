# coding=utf-8
# -*- coding: utf-8 -*-
import os
import subprocess
import shlex
import pandas as pd
import argparse
from multiprocessing import Pool
__author__ = 'gdq'


def run_cmd(cmd):
    # subprocess.check_call should not be used
    # for cannot terminate it when error happens during multiprocessing.
    subprocess.call(cmd)


class QuantExpToolbox(object):
    """
    A toolbox contains several expression quantification tools for RNAseq data.
    salmon: only  quasi-mapping-based mode of Salmon used here.
            http://salmon.readthedocs.io/en/latest/salmon.html
    fastq：a ";" separated string with 3 columns, but 3th column is optional.
    ------------------------------------
    sample1;lfq,lfq2,lfq3;rfq,rfq2,rfq3
    ------------------------------------
    """
    def __init__(self, transcriptome, fastq, method='salmon', libtype=None,
                 pool=6, threads=10, output=None,
                 salmon=None, rsem=None, kallisto=None,
                 read_len=149, read_len_sd=30,
                 map_tool='bowtie2', map_tool_path=None, samtools=None,
                 ):
        self.method = method
        self.transcriptome = transcriptome
        self.output = os.getcwd() if output is None else output
        self.salmon = salmon
        self.rsem = rsem
        self.kallisto = kallisto
        self.pool = pool
        self.threads = threads
        self.fastq = self.parse_fastq(fastq)
        self.libtype = libtype
        self.read_len = read_len
        self.read_len_sd = read_len_sd
        self.map_tool = map_tool
        self.map_tool_path = map_tool_path
        if samtools is None:
            self.samtools = "samtools"
        else:
            self.samtools = samtools

    @staticmethod
    def parse_fastq(fastq):
        fastq_info = dict()
        sample_fastqs_list = fastq.strip().split()
        for each_sample in sample_fastqs_list:
            fastq_detail = each_sample.split(';')
            sample = fastq_detail[0]
            fastq_info[sample] = list()
            for each in fastq_detail[1:]:
                fastqs = each.split(',')
                fastq_info[sample].append(fastqs)
        return fastq_info

    def build_index(self, kmer=31):
        if self.method.lower() == 'salmon':
            cmd = '{}/salmon index '.format(self.salmon)
            cmd += '-t {} '.format(self.transcriptome)
            cmd += '-i {}/transcripts_index '.format(self.output)
            cmd += '-k {}'.format(kmer)
            print(cmd)
            subprocess.check_call(shlex.split(cmd))
        elif self.method.lower() == 'kallisto':
            cmd = '{}/kallisto index '.format(self.kallisto)
            cmd += '-i {}/transcripts_index '.format(self.output)
            cmd += '-k {} '.format(kmer)
            cmd += self.transcriptome
            print(cmd)
            subprocess.check_call(shlex.split(cmd))
        elif self.method.lower() == 'rsem':
            '''
            # http://www.bioinfo-scrounger.com/archives/482 RSEM的中文使用测试，对输出文件也有详细的解释
            --gtf <file>
            If this option is on, RSEM assumes that 'reference_fasta_file(s)' contains the sequence of a genome, 
            and will extract transcript reference sequences using the gene annotations specified in <file>, 
            which should be in GTF format.
            If this option is off, RSEM will assume 'reference_fasta_file(s)' contains the reference transcripts.
             In this case, RSEM assumes that name of each sequence in the Multi-FASTA files is its transcript_id.
            (Default: off)
            # rsem跑比对定量的命令格式如下，并且可以bowtie, bowtie2一起跑
            rsem-prepare-reference --gtf mm9.gtf \
                        --star \
                        --star-path /sw/STAR \
                        -p 8 \
                        /data/mm9/chr1.fa,/data/mm9/chr2.fa,...,/data/mm9/chrM.fa \
                        /ref/mouse_0
            '''
            if self.map_tool not in ['bowtie', 'bowtie2', 'star']:
                raise Exception('{} is not supported for rsem'.format(self.map_tool))
            cmd = '{}/rsem-prepare-reference '.format(self.rsem)
            cmd += "--{} ".format(self.map_tool)
            if self.map_tool_path is not None:
                cmd += "--{}-path {} ".format(self.map_tool, self.map_tool_path)
            cmd += '-p {} '.format(self.threads)
            cmd += '{} '.format(self.transcriptome)
            index_path = '{}/transcripts_index'.format(self.output)
            if not os.path.exists(index_path):
                os.mkdir(index_path)
            cmd += '{}/ref_index'.format(index_path)
            print(cmd)
            subprocess.check_call(shlex.split(cmd))
        else:
            raise Exception('unexpected method: ' + self.method)

    def get_salmon_cmd(self):
        """
        http://salmon.readthedocs.io/en/latest/salmon.html (alignments should not be sorted by target or position.
        If your reads or alignments do not appear in a random order with respect to the target transcripts,
        please randomize / shuffle them before performing quantification with Salmon.)
        fastq: {sample name: [[fq1,fq2],[fq1,fq2]], }
        :return: cmd list
        """
        cmd_list = list()
        for sample in self.fastq:
            fq_list = self.fastq[sample]
            cmd = '{}/salmon quant '.format(self.salmon)
            cmd += '-i {}/transcripts_index '.format(self.output)
            cmd += '-l A '
            if len(fq_list) == 2:
                cmd += '-1 {} -2 {} '.format(' '.join(fq_list[0]), ' '.join(fq_list[1]))
            else:
                cmd += '-r {} '.format(' '.join(fq_list[0]))
            cmd += '-o {}/{}_quant '.format(self.output, sample)
            cmd += '--gcBias '
            cmd += '-p {} '.format(self.threads)
            print("此处是salmon")
            print(cmd)
            cmd_list.append(shlex.split(cmd))
            print(cmd_list)
        else:
            return cmd_list

    def get_kallisto_cmd(self):
        cmd_list = list()
        for sample in self.fastq:
            fq_list = self.fastq[sample]
            cmd = '{}/kallisto quant '.format(self.kallisto)
            cmd += '-i {}/transcripts_index '.format(self.output)
            cmd += '-o {}/{}_quant '.format(self.output, sample)
            # cmd += '--plaintext '
            if len(fq_list) == 1:
                cmd += '--single '
                cmd += '-l {} '.format(self.read_len)
                cmd += '-s {} '.format(self.read_len_sd)
            if self.libtype is not None:
                if self.libtype == 'fr':
                    cmd += '--fr-stranded '
                elif self.libtype == 'rf':
                    cmd += '--rf-stranded '
                else:
                    print('the library type {} is invalid and will be ignored'.format(self.libtype))
                    raise Exception('library type can only be either "fr" or "fr"')
            cmd += '-t {} '.format(self.threads)
            if len(fq_list) == 1:
                cmd += ' '.join(fq_list[0])
            else:
                fq_list = list(zip(fq_list[0], fq_list[1]))
                fq_list = [x for y in fq_list for x in y]
                cmd += ' '.join(fq_list)
            print("此处是kallisto")
            print(cmd)
            cmd_list.append(shlex.split(cmd))
            print(cmd_list)
        else:
            return cmd_list

    def get_rsem_cmd(self):
        cmd_list = list()
        for sample in self.fastq:
            fq_list = self.fastq[sample]
            cmd = '{}/rsem-calculate-expression '.format(self.rsem)
            cmd += '-p {} '.format(self.threads)
            if self.libtype is not None:
                if self.libtype == 'fr':
                    cmd += '--strandedness forward '
                elif self.libtype == 'rf':
                    cmd += '--strandedness reverse '
            cmd += '--sort-bam-memory-per-thread 2G '
            cmd += '--estimate-rspd '
            cmd += '--ci-memory 2048 '
            cmd += "--{} ".format(self.map_tool)
            if self.map_tool_path is not None:
                cmd += "--{}-path {} ".format(self.map_tool, self.map_tool_path)
            if len(fq_list) == 1:
                cmd += "--fragment-length-mean {} ".format(self.read_len)
                cmd += "--fragment-length-sd {} ".format(self.read_len_sd)
                cmd += "{} ".format(','.join(fq_list[0]))
            else:
                cmd += '--paired-end '
                cmd += '{} {} '.format(','.join(fq_list[0]), ','.join(fq_list[1]))
            cmd += "{} ".format('{}/transcripts_index/ref_index'.format(self.output))
            out_path = '{}/{}_quant'.format(self.output, sample)
            if not os.path.exists(out_path):
                os.mkdir(out_path)
            cmd += "{}/{}".format(out_path, sample)
            print("此处是rsem")
            print(cmd)
            cmd_list.append(shlex.split(cmd))
            print(cmd_list)
        else:
            return cmd_list

    def run_quant(self):
        if self.method.lower() == 'salmon':
            cmd_list = self.get_salmon_cmd()
        elif self.method.lower() == 'kallisto':
            cmd_list = self.get_kallisto_cmd()
        elif self.method.lower() == 'rsem':
            cmd_list = self.get_rsem_cmd()
        else:
            raise Exception(self.method + ' is not supported')
        if len(self.fastq) <= self.pool:
            pool = Pool(len(self.fastq))
        else:
            pool = Pool(self.pool)
        pool.map(run_cmd, cmd_list)
        pool.close()
        pool.join()
        if len(cmd_list) >= 2:
            self.generate_exp_table()

    def run_sort_bam(self):
        if self.method.lower() == 'rsem':
            samples = sorted(self.fastq.keys())
            bam_list = [os.path.join(self.output, x + '_quant', x + '.transcript.bam') for x in samples]
            cmd_list = list()
            sorted_bam_list = list()
            for each in bam_list:
                out_file = each[:-4] + '.sorted.bam'
                sorted_bam_list.append(out_file)
                cmd = '{} sort --threads {} -o {} {}'.format(self.samtools, self.threads, out_file, each)
                cmd_list.append(shlex.split(cmd))
            bam_list_file = os.path.join(self.output, 'bam.list')
            print(cmd_list)
            with open(bam_list_file, 'w') as f:
                for each in sorted_bam_list:
                    f.write(each + '\n')
            if len(self.fastq) <= self.pool:
                pool = Pool(len(self.fastq))
            else:
                pool = Pool(self.pool)
            pool.map(run_cmd, cmd_list)
            pool.close()
            pool.join()

    def generate_exp_table(self):
        def merge_file(results, target_cols, new_col_names, out):
            for each_col, new_name in zip(target_cols, new_col_names):
                column_list = list()
                for sample, quant in results:
                    tmp_col = pd.read_table(quant, index_col=0, header=0)[each_col]
                    tmp_col.name = sample
                    tmp_col.index.name = 'seq_id'
                    column_list.append(tmp_col)
                result_table = pd.concat(column_list, axis=1)
                result_table.to_csv(out+'.'+new_name+'.matrix', sep='\t')

        samples = sorted(self.fastq.keys())
        join = os.path.join
        if self.method.lower() == 'salmon':
            isoform_results = [join(self.output, x+'_quant', 'quant.sf') for x in samples]
            # gene_results = [join(self.output, x+'_quant', 'quant.genes.sf') for x in samples]
            merge_file(zip(samples, isoform_results), ['TPM', 'NumReads'], ['tpm', 'count'],
                       join(self.output, 'transcript'))
            # merge_file(zip(samples, gene_results), ['TPM', 'NumReads'], ['tpm', 'count'],
            #            join(self.output, 'gene'))  # 此处注释了3行
        elif self.method.lower() == 'rsem':
            iso_results = [join(self.output, x+'_quant', x+'.isoforms.results') for x in samples]
            # gene_results = [join(self.output, x+'_quant', x+'.genes.results') for x in samples]
            merge_file(zip(samples, iso_results), ['TPM', 'expected_count'], ['tpm', 'count'],
                       join(self.output, 'transcript'))
            # merge_file(zip(samples, gene_results), ['TPM', 'expected_count'], ['tpm', 'count'],
            #            join(self.output, 'gene'))  # 此处注释了3行
            # get mapping rate
            cnt_files = [join(self.output, x+'_quant', x+'.stat', x+'.cnt') for x in samples]
            with open(join(self.output, 'alignment_rate.txt'), 'w') as f:
                f.write('sample\ttotal_reads\taligned_reads\taligned_rate\n')
                for s, each in zip(samples, cnt_files):
                    with open(each) as f2:
                        tmp_list = f2.readline().strip('\n').split()
                    map_num = int(tmp_list[1]) + int(tmp_list[2])
                    map_rate = float(map_num)/int(tmp_list[3])
                    f.write('{}\t{}\t{}\t{}\n'.format(s, tmp_list[3], map_num, map_rate))
        elif self.method.lower() == 'kallisto':
            isoform_results = [join(self.output, x+'_quant', 'abundance.tsv') for x in samples]
            merge_file(zip(samples, isoform_results), ['tpm', 'est_counts'], ['tpm', 'count'],
                       join(self.output, 'transcript'))
        else:
            pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', type=str, metavar="transcript", required=True,
                        help="transcripts fasta file")
    parser.add_argument('-fq', type=str, metavar="fastq", required=True,
                        help="str contains the fastq path info, such as:  "
                             "'sample1;lfq,lfq2,lfq3;rfq,rfq2,rfq3 sample2;...'")
    parser.add_argument('-m', type=str, metavar="method", default='salmon',
                        help="salmon[default], kallisto and rsem are supported.")
    parser.add_argument('-o', type=str, metavar="out_dir", default=None,
                        help="Output directory. Default: current working directory.")
    parser.add_argument('-strand', type=str, metavar="library_type", default=None,
                        help="'fr' or 'rf'. Default: None. Strand-specific protocol type.")
    parser.add_argument('-pool', type=int, metavar="pool_size", default=6,
                        help="Process number for batch analysis of many samples. Default: 6")
    parser.add_argument('-thread', type=str, metavar="thread_number", default=10,
                        help="Threads number for analysis of every single sample.")
    parser.add_argument('-salmon', type=str,
                        default="/mnt/ilustre/users/sanger-dev/app/bioinfo/rna/Salmon-0.8"
                                ".2_linux_x86_64/bin/",
                        help="where is salmon installed")
    parser.add_argument('-rsem', type=str,
                        default="/mnt/ilustre/users/sanger-dev/app/bioinfo/rna/RSEM-1.2.31/bin",
                        help="where is rsem installed")
    parser.add_argument('-kallisto', type=str,
                        default="/mnt/ilustre/users/sanger-dev/app/bioinfo/rna/kallisto_linux-v0"
                                ".43.1/",
                        help="where is kallisto installed")
    parser.add_argument('-rl', type=int, metavar="read_length", default=149,
                        help="Only used for single-end mode of kallisto/rsem. "
                             "Estimated average fragment length. Default: 149")
    parser.add_argument('-sd', type=int, metavar="read_length_sd", default=30,
                        help="Only used for single-end mode of kallisto/rsem. "
                             "Use '-sd' along with '-rl. "
                             "Estimated standard deviation of fragment length. Default: 30")
    parser.add_argument('-mapper', type=str, metavar="align_tool", default='bowtie2',
                        help='RSEM argument. Only bowtie, bowtie2 and star are supported. '
                             'Default: bowtie2')
    parser.add_argument('-mapper_path', type=str, metavar="mapper_path",
                        # default='/mnt/ilustre/users/sanger-dev/app/bioinfo/align/bowtie2-2.2.9/',
                        default='/mnt/ilustre/users/sanger-dev/app/bioinfo/ref_rna_v2/miniconda2/bin/',
                        help='RSEM argument. Absolute path of alignment tool/mapper')

    args = parser.parse_args()
    toolbox = QuantExpToolbox(transcriptome=args.t,
                              fastq=args.fq,
                              method=args.m,
                              output=args.o,
                              pool=args.pool,
                              threads=args.thread,
                              libtype=args.strand,
                              salmon=args.salmon,
                              kallisto=args.kallisto,
                              rsem=args.rsem,
                              read_len=args.rl,
                              read_len_sd=args.sd,
                              map_tool=args.mapper,
                              map_tool_path=args.mapper_path,
                              samtools=None,
                              )
    toolbox.build_index()
    toolbox.run_quant()
    toolbox.run_sort_bam()

