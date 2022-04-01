#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import pandas as pd
import glob
import unittest
import shutil
from biocluster.api.file.lib.transfer import MultiFileTransfer

class QuantModule(Module):
    """
    对所有样本进行定量
    """
    def __init__(self, work_id):
        super(QuantModule, self).__init__(work_id)
        options = [
            dict(name="transcriptome", type="infile", format="denovo_rna_v2.trinity_fasta"),
            dict(name="fastq", type="infile", format="denovo_rna_v2.fastq_list"),
            dict(name="bamlist", type="outfile", format="denovo_rna_v2.bamlist"),
            dict(name="method", type="string", default="rsem"),
            dict(name="libtype", type="string", default=None),
            dict(name="t2g", type="string"),
            dict(name="pool", type="int", default=6),
            dict(name="thread", type="int", default=8),
            dict(name="output", type="string", default=None),
            dict(name="read_len", type="int", default=149),
            dict(name="read_len_sd", type="int", default=30),
            dict(name="map_tool", type="string", default="bowtie2"),
            dict(name="gene_count", type="outfile", format="denovo_rna_v2.express_matrix"),  # for output
            dict(name="transcript_count", type="outfile", format="denovo_rna_v2.express_matrix"),  # for output
            dict(name="transcript_tpm", type="outfile", format="denovo_rna_v2.express_matrix"),  # for output
            dict(name="transcript_fpkm", type="outfile", format="denovo_rna_v2.express_matrix"),  # for output
            dict(name="gene_fpkm", type="outfile", format="denovo_rna_v2.express_matrix"),  # for output
            dict(name="gene_tpm", type="outfile", format="denovo_rna_v2.express_matrix"),  # for output
            dict(name="align_rate", type="string"),  # for output
        ]
        self.add_option(options)
        self.tools = list()
        self.samples = list()

    def check_options(self):
        if self.option("method").lower() not in ["rsem", "salmon", "kallisto"]:
            raise OptionError("Method is incorrect", code = "22001101")
        if not self.option("fastq").is_set:
            raise OptionError("fastq file not exist", code = "22001102")
        if not self.option("transcriptome").is_set:
            raise OptionError("transcriptome not exist", code = "22001103")
        if self.option("map_tool").lower() not in ["bowtie", "bowtie2", "star"]:
            raise OptionError("map_tool/aligner is incorrect", code = "22001104")
        if self.option("libtype") is not None:
            if self.option("libtype").lower() not in ["fr", "rf"]:
                raise OptionError("libtype argument is not in [None, 'fr', 'rf']", code = "22001105")

    def run(self):
        super(QuantModule, self).run()
        self.tool_run()

    def tool_run(self):
        fastq_dict = self.option('fastq').to_dict()
        transfer = MultiFileTransfer()
        for sample, fq_list in sorted(fastq_dict.items()):
            if not os.path.exists(''.join(fq_list[0])):
                local_fq = self.work_dir + "/fastq/"
                self.logger.info(''.join(fq_list[0]))
                self.logger.info(local_fq)
                transfer.add_download(''.join(fq_list[0]), local_fq)
                transfer.perform()
                tool_fq_arg = sample + ";" + ''.join(local_fq + os.path.basename(''.join(fq_list[0])))
            else:
                tool_fq_arg = sample + ";" + ','.join(fq_list[0])
            if len(fq_list) >= 2:
                if not os.path.exists(''.join(fq_list[1])):
                    local_fq = self.work_dir + "/fastq/"
                    transfer.add_download(''.join(fq_list[1]), local_fq)
                    transfer.perform()
                    tool_fq_arg += ";" + ''.join(local_fq + os.path.basename(''.join(fq_list[1])))
                else:
                    tool_fq_arg += ";" + ','.join(fq_list[1])
            tool_opts = dict(
                transcriptome=self.option('transcriptome'),
                fastq=tool_fq_arg,
                method=self.option('method'),
                t2g=self.option('t2g'),
                pool=self.option("pool"),
                thread=self.option("thread"),
                read_len=self.option('read_len'),
                read_len_sd=self.option('read_len_sd'),
                map_tool=self.option('map_tool'),
            )
            tool = self.add_tool('denovo_rna_v2.single_quant')
            tool.set_options(tool_opts)
            self.tools.append(tool)
            self.samples.append(sample)
        if len(self.tools) == 1:
            self.tools[0].on("end", self.set_output)
            # self.tools[0].run()
        else:
            self.on_rely(self.tools, self.set_output, "quant")
        for tool in self.tools:
            tool.run()

    def set_output(self):
        self.logger.info("Set output of expression quantification")
        self.generate_exp_table()
        matrix_files = glob.glob(self.work_dir + '/*.matrix')
        if self.option("method").lower() == 'rsem':
            matrix_files += [self.work_dir + '/alignment_rate.txt']
            bamlist = self.work_dir + "/bam.list"
            self.option("bamlist", bamlist)
        for each in matrix_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)
            if each.endswith("transcript.count.matrix"):
                self.set_options(dict(transcript_count=link))
            elif each.endswith("transcript.tpm.matrix"):
                self.set_options(dict(transcript_tpm=link))
            elif each.endswith("gene.count.matrix"):
                self.set_options(dict(gene_count=link))
            elif each.endswith('gene.tpm.matrix'):
                self.set_options(dict(gene_tpm=link))
            elif each.endswith('gene.fpkm.matrix'):
                self.set_options(dict(gene_fpkm=link))
            elif each.endswith('transcript.fpkm.matrix'):
                self.set_options(dict(transcript_fpkm=link))
            elif each.endswith('alignment_rate.txt'):
                self.set_options(dict(align_rate=link))
            else:
                continue
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "表达定量分析结果目录"],
        ])
        super(QuantModule, self).end()

    def generate_exp_table(self):
        def merge_file(results, target_cols, new_col_names, out):
            # target_cols = ['TPM', 'count']
            # new_col_names = ['tpm', 'count']
            for each_col, new_name in zip(target_cols, new_col_names):
                column_list = list()
                for sample, quant in results:
                    tmp_col = pd.read_table(quant, index_col=0, header=0)[each_col]
                    # tmp_col.name = sample + '_' + new_name
                    tmp_col.name = sample
                    tmp_col.index.name = 'seq_id'
                    column_list.append(tmp_col)
                result_table = pd.concat(column_list, axis=1)
                result_table = result_table.loc[:, self.samples]
                result_table.to_csv(out+'.'+new_name+'.matrix', sep='\t')

        samples = self.samples
        join = os.path.join
        if self.option('method').lower() == 'salmon':
            samples = list()
            iso_results = list()
            gene_results = list()
            for each_tool in self.tools:
                target_dir = glob.glob(each_tool.work_dir + '/*_quant/quant.sf')
                iso_results.append(target_dir[0])
                target_dir = glob.glob(each_tool.work_dir + '/*_quant/quant.genes.sf')
                gene_results.append(target_dir[0])
                samples.append(each_tool.option('fastq').split(';')[0])
            merge_file(zip(samples, iso_results), ['TPM', 'NumReads'], ['tpm', 'count'],
                       join(self.work_dir, 'transcript'))
            merge_file(zip(samples, gene_results), ['TPM', 'NumReads'], ['tpm', 'count'],
                       join(self.work_dir, 'gene'))
        elif self.option('method').lower() == 'rsem':
            iso_results = list()
            gene_results = list()
            cnt_files = list()
            bam_list = list()
            samples = list()
            for each_tool in self.tools:
                target_dir = glob.glob(each_tool.work_dir + '/*_quant/*.isoforms.results')
                iso_results.append(target_dir[0])
                target_dir = glob.glob(each_tool.work_dir + '/*_quant/*.genes.results')
                gene_results.append(target_dir[0])
                target_dir = glob.glob(each_tool.work_dir + '/*_quant/*.stat/*.cnt')
                cnt_files.append(target_dir[0])
                target_dir = glob.glob(each_tool.work_dir + '/*_quant/*.sorted.bam')
                bam_list.append(target_dir[0])
                samples.append(target_dir[0].split('/')[-2].split('_quant')[0])
            self.logger.info(samples)
            self.logger.info(iso_results)

            def order_sample(data, ordered_samples):
                tmp_dict = dict(data)
                new_data = list()
                for each in ordered_samples:
                    new_data.append((each, tmp_dict[each]))
                return new_data

            merge_file(zip(samples, iso_results), ['TPM', 'expected_count', 'FPKM'], ['tpm', 'count', 'fpkm'],
                       join(self.work_dir, 'transcript'))
            merge_file(zip(samples, gene_results), ['TPM', 'expected_count', 'FPKM'], ['tpm', 'count', 'fpkm'],
                       join(self.work_dir, 'gene'))
            # get mapping rate
            with open(join(self.work_dir, 'alignment_rate.txt'), 'w') as f:
                f.write('sample\ttotal_reads\taligned_reads\taligned_rate\n')
                for s, each in order_sample(zip(samples, cnt_files), self.samples):
                    with open(each) as f2:
                        # ['un-alignable', 'alignable', 'too_many_align', 'total']
                        tmp_list = f2.readline().strip('\n').split()
                    map_num = int(tmp_list[1]) + int(tmp_list[2])
                    map_rate = float(map_num)/int(tmp_list[3])
                    f.write('{}\t{}\t{}\t{}\n'.format(s, tmp_list[3], map_num, map_rate))
            # write bam list file
            with open(join(self.work_dir, 'bam.list'), 'w') as f:
                for _x, bam in order_sample(zip(samples, bam_list), self.samples):
                    f.write(bam + '\n')
        elif self.option('method').lower() == 'kallisto':
            samples = list()
            isoform_results = list()
            for each_tool in self.tools:
                target_dir = glob.glob(each_tool.work_dir + '/*_quant/abundance.tsv')
                isoform_results.append(target_dir[0])
                samples.append(each_tool.option('fastq').split(';')[0])
            merge_file(zip(samples, isoform_results), ['tpm', 'est_counts'], ['tpm', 'count'],
                       join(self.work_dir, 'transcript'))
            # kallisto do not generate gene expression table, thus t2g file will be needed here.
            t2g_pd = pd.read_table(self.option('t2g'), index_col=0, header=None, usecols=[0, 1], names=['seq_id', 'gene'])
            gene_exp_list = list()
            gene_count_list = list()
            for sample, result in zip(samples, isoform_results):
                iso_table = pd.read_table(result, index_col=0, header=0)
                iso_table = pd.concat([iso_table, t2g_pd], axis=1)
                tmp_table = iso_table['est_counts'].groupby(iso_table['gene']).sum()
                tmp_table.name = sample
                tmp_table.index.name = 'seq_id'
                gene_count_list.append(tmp_table)
                # exp
                tmp_table = iso_table['tpm'].groupby(iso_table['gene']).sum()
                tmp_table.name = sample
                tmp_table.index.name = 'seq_id'
                gene_exp_list.append(tmp_table)
            gene_counts = pd.concat(gene_count_list, axis=1)
            gene_counts.to_csv(os.path.join(self.work_dir, 'gene.count.matrix'), sep='\t')
            gene_exp = pd.concat(gene_exp_list, axis=1)
            gene_exp.to_csv(os.path.join(self.work_dir, 'gene.tpm.matrix'), sep='\t')
        else:
            pass


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "Quant" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "denovo_rna_v2.quant",
            "instant": False,
            "options": dict(
                transcriptome="/mnt/ilustre/users/sanger-dev/workspace/20190709/Single_denovo_assemble21868/DenovoAssemble2/output/Trinity.filter.fasta",
                fastq="/mnt/ilustre/users/sanger-dev/workspace/20190708/Single_test_HiseqQc_5565-yyyyyyyy/HiseqQc/output/sickle_dir/fq_list.txt",
                method="rsem",
                t2g="/mnt/ilustre/users/sanger-dev/workspace/20190709/Single_denovo_assemble21868/DenovoAssemble2/output/Trinity.filter.gene_trans_map",
                pool=6,
                thread=6,
                output=None,
                read_len=149,
                read_len_sd=30,
                map_tool="bowtie2",
            )
        }
        # data['options']['method'] = 'rsem'
        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.run()
        # #
        # data['id'] += '1'
        # data['options']['method'] = 'salmon'
        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.run()
        #
        data['id'] += 'gdq'
        data['options']['method'] = 'rsem'
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
