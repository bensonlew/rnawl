#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from re import L
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import pandas as pd
import glob
import unittest
import shutil
from collections import OrderedDict
from biocluster.file import exists
from biocluster.file import download
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.config import Config


class QuantModule(Module):
    """
    对所有样本进行定釄1�7
    """
    def __init__(self, work_id):
        super(QuantModule, self).__init__(work_id)
        options = [
            dict(name="biotype", type="infile", format="prok_rna.common"),
            dict(name="transcriptome", type="infile", format="prok_rna.fasta"),
            dict(name="fastq", type="infile", format="prok_rna.fastq_list"),
            dict(name="bamlist", type="outfile", format="prok_rna.bamlist"),
            dict(name="method", type="string", default="rsem"),
            dict(name="libtype", type="string", default=None),
            dict(name="pool", type="int", default=6),
            dict(name="thread", type="int", default=6),
            dict(name="output", type="string", default=None),
            dict(name="read_len", type="int", default=149),
            dict(name="read_len_sd", type="int", default=30),
            dict(name="map_tool", type="string", default="bowtie2"),
            dict(name="gene_count", type="outfile", format="prok_rna.express_matrix"),  # for output
            dict(name="transcript_count", type="outfile", format="prok_rna.express_matrix"),  # for output
            dict(name="transcript_tpm", type="outfile", format="prok_rna.express_matrix"),  # for output
            dict(name="gene_tpm", type="outfile", format="prok_rna.express_matrix"),  # for output
            dict(name="gene_fpkm", type="outfile", format="prok_rna.express_matrix"),  # for output
            dict(name="transcript_fpkm", type="outfile", format="prok_rna.express_matrix"),  # for output
            dict(name="align_rate", type="string"),  # for output
            dict(name="id2name", type="string"),  # for convert gene_id to gene_name
            dict(name="task_id", type="string"),  # for getting upload dir
        ]
        self.add_option(options)
        self.tools = list()
        self.samples = list()

    def check_options(self):
        if self.option("method").lower() not in ["rsem", "salmon", "kallisto"]:
            raise OptionError("method is incorrect", code = "25000801")
        if not self.option("fastq").is_set:
            raise OptionError("fastq file not exist", code = "25000802")
        if not self.option("transcriptome").is_set:
            raise OptionError("transcriptome not exist", code = "25000803")
        if self.option("map_tool").lower() not in ["bowtie", "bowtie2", "star"]:
            raise OptionError("map_tool/aligner is incorrect", code = "25000804")
        if self.option("libtype") is not None:
            if self.option("libtype").lower() not in ["fr", "rf", "r", "f"]:
                raise OptionError("libtype argument is not in [None, 'fr', 'rf', 'r', 'f']", code = "25000805")

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
                if "workspace" in (''.join(fq_list[0])):
                    db = Config().get_mongo_client(mtype="prok_rna")[Config().get_mongo_dbname("prok_rna")]
                    col = db["sg_task"]
                    task_info = col.find_one({"task_id" : self.option("task_id")})
                    dir = os.path.dirname(task_info["fastq"])
                    base_name = os.path.basename(''.join(fq_list[0]))
                    fq_list[0] = os.path.join(dir, base_name)
                    self.logger.info(''.join(fq_list[0]))
                transfer.add_download(''.join(fq_list[0]), local_fq)
                transfer.perform()
                tool_fq_arg = sample + ";" + ''.join(local_fq + os.path.basename(''.join(fq_list[0])))
            else:
                tool_fq_arg = sample + ";" + ','.join(fq_list[0])
            if len(fq_list) >= 2:
                if not os.path.exists(''.join(fq_list[1])):
                    local_fq = self.work_dir + "/fastq/"
                    if "workspace" in (''.join(fq_list[1])):
                        db = Config().get_mongo_client(mtype="prok_rna")[Config().get_mongo_dbname("prok_rna")]
                        col = db["sg_task"]
                        task_info = col.find_one({"task_id" : self.option("task_id")})
                        dir = os.path.dirname(task_info["fastq"])
                        base_name = os.path.basename(''.join(fq_list[1]))
                        fq_list[1] = os.path.join(dir, base_name)
                        self.logger.info(''.join(fq_list[1]))
                    transfer.add_download(''.join(fq_list[1]), local_fq)
                    transfer.perform()
                    tool_fq_arg += ";" + ''.join(local_fq + os.path.basename(''.join(fq_list[1])))
                else:
                    tool_fq_arg += ";" + ','.join(fq_list[1])
        # for sample, fq_list in fastq_dict.items():
        #     tool_fq_arg = sample + ";" + ','.join(fq_list[0])
        #     if len(fq_list) >= 2:
        #         tool_fq_arg += ";" + ','.join(fq_list[1])
            tool_opts = dict(
                transcriptome=self.option('transcriptome'),
                fastq=tool_fq_arg,
                method=self.option('method'),
                pool=self.option("pool"),
                thread=self.option("thread"),
                read_len=self.option('read_len'),
                read_len_sd=self.option('read_len_sd'),
                map_tool=self.option('map_tool'),
                libtype=self.option('libtype'),
            )
            tool = self.add_tool('prok_rna.single_quant')
            tool.set_options(tool_opts)
            self.tools.append(tool)
            self.samples.append(sample)
        if len(self.tools) == 1:
            self.tools[0].on("end", self.set_output)
            # self.tools[0].run()
        else:
            self.on_rely(self.tools, self.set_output, "quant")
        for tool in self.tools:
            # self.logger.info(tool.option('fastq'))
            tool.run()
        # self.logger.info(self.samples)

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
            self.logger.info(results)
            self.logger.info(self.samples)
            # target_cols = ['TPM', 'count']
            # new_col_names = ['tpm', 'count']
            for each_col, new_name in zip(target_cols, new_col_names):
                column_list = list()
                for sample, quant in results:
                    self.logger.info(sample)
                    self.logger.info(quant)
                    tmp_col = pd.read_table(quant, index_col=0, header=0)[each_col]
                    # tmp_col.name = sample + '_' + new_name
                    tmp_col.name = sample
                    tmp_col.index.name = 'seq_id'
                    column_list.append(tmp_col)
                result_table = pd.concat(column_list, axis=1)
                result_table = result_table.loc[:, self.samples]
                result_table.to_csv(out+'.'+new_name+'.matrix', sep='\t')

        samples = self.samples

        # 下面为给文件添加丢�列来区分到底是srna，还是mrna, 上传的文件名字开头都是transcript

        rna_biotype = dict()
        #print self.option("biotype").prop["path"]
        with open(self.option("biotype").prop["path"], "r") as f1:
            for line in f1:
                items = line.strip().split("\t")
                if items[0] not in rna_biotype:
                    rna_biotype[items[0]] = items[1]

        def get_rna_biotype(line):
            return rna_biotype.get(line[0], "unknown")

        def get_rna_type(line):
            if line[0].startswith("sRNA"):
                return "sRNA"
            elif line[0].startswith("novel"):
                return "novel_gene"
            else:
                return "ref_gene"

        df_id2name = pd.read_table(self.option("id2name"), header=0, sep="\t")
        dict_id2name = OrderedDict(zip(df_id2name.iloc[:, 6], df_id2name.iloc[:, 5]))

        def convert_id_2_name(row):
            rna_type = row["seq_id"]
            if type(rna_type) == str and "srna" in rna_type.lower():
                return "-"
            else:
                if rna_type in dict_id2name.keys():
                    return dict_id2name[rna_type]

        join = os.path.join
        if self.option('method').lower() == 'salmon':
            iso_results = list()
            samples = list()
            # gene_results = list()  # 此处注释1衄1�7
            for each_tool in self.tools:
                target_dir = glob.glob(each_tool.work_dir + '/*_quant/quant.sf')
                iso_results.append(target_dir[0])
                samples.append(each_tool.option('fastq').split(';')[0])
                # target_dir = glob.glob(each_tool.work_dir + '/*_quant/quant.genes.sf')
                # gene_results.append(target_dir[0])  # 此处注释2衄1�7
            merge_file(zip(samples, iso_results), ['TPM', 'NumReads'], ['tpm', 'count'],
                       join(self.work_dir, 'transcript'))

            # 下划线之间为新添加分割这个结果文件的代码
            # salmon只包含tpm, count的结构1�7
            # -----------------------------------------------------------------------------------------------
            # if self.option("libtype"):
            list_type = ["transcript.tpm.matrix", "transcript.count.matrix"]
            for i in list_type:
                df_whole = pd.read_table(join(self.work_dir, i), header=0, dtype={0:str})
                df_whole['type'] = df_whole.apply(get_rna_type, axis=1)
                df_whole['biotype'] = df_whole.apply(get_rna_biotype, axis=1)
                df_whole['gene_name'] = df_whole.apply(convert_id_2_name, axis=1)
                new_columns = df_whole.columns[1:-1].tolist()  # 这个写成亄1�71＄1�7-2就把倒数第二列给漏掉了，要注意这个问预1�7
                new_columns.insert(0, df_whole.columns[0])
                new_columns.insert(1, df_whole.columns[-1])
                df_whole = df_whole[new_columns]
                os.remove(join(self.work_dir, i))
                df_whole.to_csv(join(self.work_dir, i), index=False, sep="\t")
                df_ref_gene = df_whole[df_whole["seq_id"].map(lambda x: not(x.startswith("sRNA") or x.startswith("novel")))]
                df_ref_gene.to_csv(join(self.work_dir, i.replace("transcript", "ref_gene")), index=False, sep="\t")
            # -----------------------------------------------------------------------------------------------

            # merge_file(zip(samples, gene_results), ['TPM', 'NumReads'], ['tpm', 'count'],
            #            join(self.work_dir, 'gene'))  # 此处注释2衄1�7
        elif self.option('method').lower() == 'rsem':
            iso_results = list()
            # gene_results = list()  # 此处注释1衄1�7
            cnt_files = list()
            bam_list = list()
            samples = list()
            for each_tool in self.tools:
                target_dir = glob.glob(each_tool.work_dir + '/*_quant/*.isoforms.results')
                iso_results.append(target_dir[0])
                # target_dir = glob.glob(each_tool.work_dir + '/*_quant/*.genes.results')
                # gene_results.append(target_dir[0])  # 此处注释2衄1�7
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
            # merge_file(zip(samples, gene_results), ['TPM', 'expected_count', 'FPKM'], ['tpm', 'count', 'fpkm'],
            #            join(self.work_dir, 'gene'))  # 此处注释2衄1�7

            # 下划线之间为新添加分割这个结果文件的代码
            # rsem朄1�73种类型的结果
            # -----------------------------------------------------------------------------------------------
            # if self.option("libtype"):
            list_type = ["transcript.tpm.matrix", "transcript.fpkm.matrix", "transcript.count.matrix"]
            for i in list_type:
                df_whole = pd.read_table(join(self.work_dir, i), header=0)
                df_whole['type'] = df_whole.apply(get_rna_type, axis=1)
                df_whole['biotype'] = df_whole.apply(get_rna_biotype, axis=1)
                df_whole['gene_name'] = df_whole.apply(convert_id_2_name, axis=1)
                new_columns = df_whole.columns[1:-1].tolist()  # 这个写成亄1�71＄1�7-2就把倒数第二列给漏掉了，要注意这个问预1�7
                new_columns.insert(0, df_whole.columns[0])
                new_columns.insert(1, df_whole.columns[-1])
                df_whole = df_whole[new_columns]
                os.remove(join(self.work_dir, i))
                df_whole.to_csv(join(self.work_dir, i), index=False, sep="\t")
                df_ref_gene = df_whole[df_whole["seq_id"].map(lambda x: not(x.startswith("sRNA") or x.startswith("novel")))]
                df_ref_gene.to_csv(join(self.work_dir, i.replace("transcript", "ref_gene")), index=False, sep="\t")
            # -----------------------------------------------------------------------------------------------

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

            # 下划线之间为新添加分割这个结果文件的代码
            # kallisto只包含tpm, count的结构1�7
            # -----------------------------------------------------------------------------------------------
            # if self.option("libtype"):
            list_type = ["transcript.tpm.matrix", "transcript.count.matrix"]
            for i in list_type:
                df_whole = pd.read_table(join(self.work_dir, i), header=0)
                df_whole['type'] = df_whole.apply(get_rna_type, axis=1)
                df_whole['biotype'] = df_whole.apply(get_rna_biotype, axis=1)
                df_whole['gene_name'] = df_whole.apply(convert_id_2_name, axis=1)
                new_columns = df_whole.columns[1:-1].tolist()  # 这个写成亄1�71＄1�7-2就把倒数第二列给漏掉了，要注意这个问预1�7
                new_columns.insert(0, df_whole.columns[0])
                new_columns.insert(1, df_whole.columns[-1])
                df_whole = df_whole[new_columns]
                os.remove(join(self.work_dir, i))
                df_whole.to_csv(join(self.work_dir, i), index=False, sep="\t")
                df_ref_gene = df_whole[df_whole["seq_id"].map(lambda x: not(x.startswith("sRNA") or x.startswith("novel")))]
                df_ref_gene.to_csv(join(self.work_dir, "ref_gene.tpm.matrix"), index=False, sep="\t")
            # -----------------------------------------------------------------------------------------------

            # kallisto do not generate gene expression table, thus t2g file will be needed here.
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
            "name": "prok_rna.quant",
            "instant": False,
            "options": dict(
                transcriptome="/mnt/ilustre/users/sanger-dev/workspace/20180801/Single_Srna_2589_fyt/Srna/Rockhopper/Rockhopper_Results/genome.feature.fa",
                fastq="/mnt/ilustre/users/sanger-dev/workspace/20180806/Single_HiseqQc_3546/HiseqQc/output/sickle_dir/fq_list.txt",
                method="rsem",
                pool=6,
                thread=6,
                output=None,
                read_len=149,
                read_len_sd=30,
                map_tool="bowtie2",
                # libtype=null,
                id2name="/mnt/ilustre/users/sanger-dev/sg-users/litangjian/prok_rna/gene2name",
            )
        }
        data['options']['method'] = 'rsem'
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()
        # #
        # data['id'] += '1'
        # data['options']['method'] = 'salmon'
        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.run()

        # data['id'] += 'gdq'
        # data['options']['method'] = 'kallisto'
        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.run()


if __name__ == '__main__':
    unittest.main()
