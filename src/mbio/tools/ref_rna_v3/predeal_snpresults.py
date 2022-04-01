# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180408

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import shutil
import pandas as pd
from collections import Counter
import unittest

class PredealSnpresultsAgent(Agent):
    """
    author = HONGDONG
    version = 1.0
    last modify = 20180408
    该tool用于对vcf进行处理，首先将分割后的vcf文件使用concat进行合并，然后annotate，最后在进行sort排序生成最后需要的vcf文件
    """
    def __init__(self, parent):
        super(PredealSnpresultsAgent, self).__init__(parent)
        options = [
            {"name": "snp_detail", "type": "string"},  # snp_anno.xls
            {"name": "snp_anno", "type": "string"},  # snp_annotation.xls
            {"name": "method", "type": "string", "default": "samtools"},  # snp_annotation.xls
        ]
        self.add_option(options)
        self._memory_increase_step = 16
        
    def check_options(self):
        if not os.path.exists(self.option("snp_detail")):
            OptionError("缺少snp_anno.xls文件", code = "33704111")
        if not os.path.exists(self.option("snp_anno")):
            OptionError("缺少snp_annotation.xls文件", code = "33704112")
        if not self.option("method") in ["samtools", "gatk"]:
            OptionError("只会有samtools和gatk这两种方法", code = "33704113")
        self.infile_size = os.path.getsize(self.option('snp_detail'))

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 1
        self._memory = '{}G'.format(int(float(self.infile_size) / 1024 ** 3 * 48 + 16))
        
    def end(self):
        super(PredealSnpresultsAgent, self).end()


class PredealSnpresultsTool(Tool):
    def __init__(self, config):
        super(PredealSnpresultsTool, self).__init__(config)

    def deal_snp_detail(self, snp_anno, new_output=None):
        """
        导入SNP详情表的函数
        :param snp_anno: snp模块的结果文件夹，即~/output_dir/
        :return:
        """
        snp_type_stat = {}
        snp_pos_stat = {}
        indel_pos_stat = {}
        all_depth_stat = {}
        all_freq_stat = {}
        depth_list = ["<=30", "31-100", "101-200", "201-300", "301-400", "401-500", ">501"]
        # graph_data_list = []
        chroms = set()
        distributions = set()
        data_list = []
        sample_old_index = []
        sample_old_gene = []

        if self.option('method') == 'samtools':
            df = pd.read_table(snp_anno, header=0, sep='\t', low_memory=False)
            os.remove(snp_anno)
            df.columns = df.columns.map(lambda x: x.split("bam/")[1] if "bam" in x else x)
            df.to_csv(snp_anno, header=True, sep="\t", index=False)

        with open(snp_anno, "r") as f:
            heads_info= f.readline().strip().split("\t")
            col_nums = len(heads_info)
            sam_nums = (len(heads_info) - 12) / 2
            sample_names = heads_info[11:11+sam_nums]
            sample_genotypes=heads_info[11+sam_nums:]
            for s in sample_names:
                snp_type_stat[s] = {}
                snp_pos_stat[s] = {}
                indel_pos_stat[s] = {}
                all_freq_stat[s] = [0]
                all_depth_stat[s] = [0, 0, 0, 0, 0, 0, 0]
                sample_old_index.append(-1)
                sample_old_gene.append('')
            # print f.next()
            for line in f:
                line = line.strip().split("\t")
                col_num=len(line)
                sam_num=(len(line)-12)/2
                sample_depthinfos = line[11:11+sam_num]
                sample_genotype_infos = line[11+sam_num:]
                new_gene_name = line[7]
                chroms.add(line[0])
                distributions.add(line[6])
                snp_type = line[3] + "/" + line[4]
                data = {
                    "type": line[-1],
                    "chrom": line[0],
                    "start": line[1],
                    "end": line[2],
                    "ref": line[3],
                    "alt": line[4],
                    "qual": float(line[5]),
                    "reads_num": int(line[8]),
                    "anno": line[6],
                    "gene": line[7],
                    "mut_type": line[9],
                    "mut_info": line[10],
                    "snp_type": snp_type
                }
                for n, s in enumerate(sample_names):
                    sad = sample_depthinfos[n]
                    mut_rate = 0
                    depth_num = -1
                    if sad != "-" and sad != ".":
                        sample_tdepth=sum([int(i) for i in sad.split(",")])
                        if sample_tdepth > 0:
                            # single_and_all = rate.split("/")[1]
                            # single_and_alt=rate.split("/")[0]
                            #mut_rate = round(float(rate.split("/")[0])/float(rate.split("/")[1]), 4)
                            # 统计各样本的突变频率数目
                            # if line[6] == "exonic":
                            #     sample_old_index[n], all_freq_stat[s] = self.gene_num_stat(new_gene_name, sample_old_gene[n], sample_old_index[n], all_freq_stat[s])
                            #     sample_old_gene[n] = new_gene_name
                            #  统计各样本的SNP类型数目
                            if not '-' in snp_type and len(snp_type) == 3:
                                snp_type_stat[s] = self.type_stat(snp_type, snp_type_stat[s])
                            depth_num = sample_tdepth

                            # 统计SNP/Indel位置信息
                            if data["type"] == "snp":
                                snp_pos_stat[s] = self.type_stat(data["anno"], snp_pos_stat[s])
                            else:
                                indel_pos_stat[s] = self.type_stat(data["anno"], indel_pos_stat[s])
                        #data[s + "_mut_rate"] = mut_rate
                        #data[s + "_reads_rate"] = single_and_all
                            #data[s]=str(int(single_and_all)-int(single_and_alt))+","+single_and_alt
                            data[s+"_sampledepth"] = sad
                        else:
                            data[s+"_sampledepth"] = sad
                        all_depth_stat[s] = self.get_depth_stat(depth_num, all_depth_stat[s])
                    else:
                        data[s + "_sampledepth"] = sad
                #data_list.append(data)
                for n, s in enumerate(sample_genotypes):
                    genotype = sample_genotype_infos[n]
                    # if type == './.' or type=='0/0':
                    #     genetype=line[3]+"/"+line[3]
                    # if type =="0/1":
                    #     genetype = line[3] + "/" + line[4]
                    # if type =="1/1":
                    #     genetype = line[4] + "/" + line[4]
                    data[s] = genotype
                data_list.append(data)
                # #统计reads深度
        if len(data_list) >= 1:
            df_data_list_pre = pd.DataFrame(data_list)
            df_data_list_pre.to_csv(new_output + "/data_anno_pre.xls", sep="\t", header=True, index=False)
        depth_data = self.get_stat_data(sample_names, depth_list, all_depth_stat, "depth_stat")
        df_depth = pd.DataFrame(depth_data)
        df_depth.to_csv(new_output + "/snp_depth_statistics.xls", sep="\t", header=True, index=False)

        #new_freq_stat = self.freq_stat(all_freq_stat)
        #freq_data = self.get_stat_data(sample_names, [1, 2, 3, 4, ">=5"], new_freq_stat, "freq_stat")
        #freq_data = pd.DataFrame(freq_data)
        #freq_data.to_csv(new_output + "/snp_freq_statistics.xls", sep="\t", header=True, index=False)

        snp_types = ["A/G", "A/C", "C/T", "G/A", "G/C", "C/A", "A/T", "C/G", "G/T", "T/C", "T/A", "T/G"]
        type_data = self.get_stat_data(sample_names, snp_types, snp_type_stat, "type_stat")
        type_data = pd.DataFrame(type_data)
        type_data.to_csv(new_output + "/snp_transition_tranversion_statistics.xls", sep="\t", header=True, index=False)

        snp_pos_data = self.get_stat_data(sample_names, list(distributions), snp_pos_stat, "snp_distribution")
        snp_pos_data = pd.DataFrame(snp_pos_data)
        snp_pos_data.to_csv(new_output + "/snp_position_distribution.xls", sep="\t", header=True, index=False)

        indel_pos_data = self.get_stat_data(sample_names, list(distributions), indel_pos_stat, "indel_distribution")
        indel_pos_data = pd.DataFrame(indel_pos_data)
        indel_pos_data.to_csv(new_output + "/indel_position_distribution.xls", sep="\t", header=True, index=False)
        with open(new_output + '/main_info', 'w') as main_w:
            main_w.write("snp_types" + '\t' + ';'.join(snp_types) + '\n')
            main_w.write("specimen" + '\t' + ';'.join(sample_names) + '\n')
            main_w.write("distributions" + '\t' + ';'.join(list(distributions)) + '\n')
            main_w.write("chroms" + '\t' + ';'.join(list(chroms)) + '\n')

    def get_depth_stat(self, depth_num, target_list):
        if depth_num == -1:
            pass
        else:
            if depth_num < 31:
                target_list[0] += 1
            elif 30 < depth_num < 101:
                target_list[1] += 1
            elif 100 < depth_num < 201:
                target_list[2] += 1
            elif 200 < depth_num < 301:
                target_list[3] += 1
            elif 300 < depth_num < 401:
                target_list[4] += 1
            elif 400 < depth_num < 501:
                target_list[5] += 1
            else:
                target_list[6] += 1
        # print target_dict
        return target_list

    def type_stat(self, dict_key, target_dict):
        if dict_key in target_dict:
            target_dict[dict_key] += 1
        else:
            target_dict[dict_key] = 1
        return target_dict

    def get_stat_data(self, sample_names, range_list, value_dict, stat_type):
        stat_data = []
        for n, ds in enumerate(range_list):
            data = {
                "range_key": ds,
                "stat_type": stat_type
            }
            for s in sample_names:
                if type(value_dict[s]) is list:
                    data["{}".format(s)] = value_dict[s][n]
                else:
                    if ds not in value_dict[s].keys():
                        value_dict[s][ds] = 0
                    data["{}".format(s)] = value_dict[s][ds]
            stat_data.append(data)
        # print stat_data
        return stat_data

    def gene_num_stat(self, new_gene_name, old_gene_name, old_index, new_list):
        if new_gene_name == old_gene_name:
            new_list[old_index] += 1
        else:
            new_list.append(1)
            old_index += 1
        return old_index, new_list

    def freq_stat(self, all_freq_stat):
        for a_s in all_freq_stat:
            c = Counter(all_freq_stat[a_s])
            values = c.values()
            keys = c.keys()
            new_d = {">=5": 0}
            for n, s in enumerate(keys):
                # print s
                if s < 5:
                    new_d[s] = values[n]
                else:
                    new_d[">=5"] += values[n]
            all_freq_stat[a_s] = new_d
        return all_freq_stat

    def merge(self, x, y):
        df_anno_pre1 = pd.read_table(x, header=0, sep="\t", low_memory=False)

        # tmp_list = df_anno_pre1.columns[:-13].tolist()
        # tmp_list.append(df_anno_pre1.columns[-1])
        tmp_list = df_anno_pre1.columns.drop(['alt','anno','chrom','end','gene','mut_info', 'mut_type', 'qual','reads_num','ref','snp_type','start', 'type']).tolist()
        tmp_list.append("type")
        df_anno_pre1_select = df_anno_pre1.loc[:, tmp_list]
        df_anno_pre1_select["index1"] = df_anno_pre1['alt'].apply(str) + df_anno_pre1['anno'].apply(str) + \
        df_anno_pre1['chrom'].apply(str) + df_anno_pre1["end"].apply(str) + df_anno_pre1["start"].apply(str) + df_anno_pre1["ref"].apply(str)
        df_anno_pre2 = pd.read_table(y, header=0, sep="\t", low_memory=False)
        df_anno_pre2_select_pre = df_anno_pre2.loc[:, df_anno_pre2.columns[:13]]
        df_anno_pre2_select_pre.rename(columns={"Depth": "Total depth", "CHROM": "Chrom", "ALT": "Alt", "ANNO": "Anno", "END": "End", "START": "Start", "REF": "Ref",
                                                "MUT_type": "MUT type", "MUT_info": "MUT info", "": ""}, inplace=True)
        df_anno_pre2_select = df_anno_pre2_select_pre[["GENE(in or nearby)", "Gene name", "Gene description", "Chrom", "Start", "End", "Ref", "Alt", "Total depth", "QUAL",
                                                      "Anno", "MUT type", "MUT info"]]
        df_anno_pre2_select["index1"] = df_anno_pre2['ALT'].apply(str) + df_anno_pre2['ANNO'].apply(str) + df_anno_pre2['CHROM'].apply(str) + \
                                        df_anno_pre2["END"].apply(str) + df_anno_pre2["START"].apply(str) + df_anno_pre2["REF"].apply(str)
        df_join_pre = pd.merge(df_anno_pre2_select, df_anno_pre1_select, on="index1", how="outer")
        df_join_pre.drop(columns=['index1'], inplace=True)
        df_join_pre.to_csv(self.output_dir + "/snp_annotation_statistics.xls", sep="\t", index=False)

    def run(self):
        """
        运行
        """
        super(PredealSnpresultsTool, self).run()
        self.logger.info("对snp生成的结果文件进行统计等预处理")
        self.deal_snp_detail(self.option('snp_detail'), new_output=self.output_dir)
        if os.path.exists(os.path.join(self.output_dir, 'data_anno_pre.xls')):
            self.merge(os.path.join(self.output_dir, 'data_anno_pre.xls'), self.option('snp_anno'))
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna'
        data = {
            "id": "snp" + str(random.randint(1, 10000))+"-xxx",
            "type": "tool",
            "name": "ref_rna_v3.predeal_snpresults",
            "instant": False,
            "options": dict(
                #ref_dict=test_dir + "/" + "Mus_musculus.GRCm38.dna_rm.toplevel.clean.dict",
                #bam_list="/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/bamlist_new",
                #call_type="sentieon",
                #ref_fasta=test_dir+"/"+"Mus_musculus.GRCm38.dna_rm.toplevel.clean.fa",
                # scm="complete",
                # scd="correlation",
                # corr_method='pearson',
                #output=None,
                #ref_gtf="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/gtf/Mus_musculus.GRCm38.89.gtf",
                #input_file="/mnt/ilustre/users/sanger-dev/workspace/20190304/Refrna_tsg_33538/SamRna/BcftoolVcf/output/variant.vcf",
                #ref_genome="customer_mode"
                #snp_detail="/mnt/ilustre/users/sanger-dev/workspace/20190617/Snp_tsg_33912_1820_7969/SnpRna/Annovar/snp_anno.xls",
                #snp_anno="/mnt/ilustre/users/sanger-dev/workspace/20190617/Snp_tsg_33912_1820_7969/SnpRna/Annovar/snp_annotation.xls",
                # snp_detail="/mnt/ilustre/users/sanger-dev/workspace/20190709/Snp_tsg_34785_2567_2561/CallSnpIndel/Annovar/snp_anno.xls",
                # snp_anno="/mnt/ilustre/users/sanger-dev/workspace/20190709/Snp_tsg_34785_2567_2561/CallSnpIndel/Annovar/snp_annotation.xls",
                snp_detail="/mnt/ilustre/users/sanger-dev/sg-users/xuxi/predeal_snpresults_tmp_file/snp_anno_1G.xls",
                snp_anno="/mnt/ilustre/users/sanger-dev/sg-users/xuxi/predeal_snpresults_tmp_file/snp_annotation_1G.xls",
                method="gatk"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()





