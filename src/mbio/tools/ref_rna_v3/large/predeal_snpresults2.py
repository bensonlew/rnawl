# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import shutil
import pandas as pd
from collections import Counter
import unittest
import dask
from itertools import chain

class PredealSnpresults2Agent(Agent):
    """
    author = SHICAIPING
    version = 1.0
    last modify = 20210413
    该tool用于对vcf进行处理，首先将分割后的vcf文件使用concat进行合并，然后annotate，最后在进行sort排序生成最后需要的vcf文件
    """
    def __init__(self, parent):
        super(PredealSnpresults2Agent, self).__init__(parent)
        options = [
            {"name": "snp_detail_list", "type": "string"},  # list of snp_anno.xls
            {"name": "snp_anno_list", "type": "string"},  # list of snp_annotation.xls
            {"name": "method", "type": "string", "default": "samtools"},
        ]
        self.add_option(options)
        self._memory_increase_step = 50
        
    def check_options(self):
        if not os.path.exists(self.option("snp_detail_list")):
            OptionError("缺少snp_detail_list文件")
        if not os.path.exists(self.option("snp_anno_list")):
            OptionError("缺少snp_anno_list文件")
        if not self.option("method") in ["samtools", "gatk"]:
            OptionError("只会有samtools和gatk这两种方法", code = "33704113")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 20
        self._memory = '100G'
        # self._memory = '{}G'.format('150')
        
    def end(self):
        super(PredealSnpresults2Agent, self).end()


class PredealSnpresults2Tool(Tool):
    def __init__(self, config):
        super(PredealSnpresults2Tool, self).__init__(config)

    def deal_snp_detail(self, snp_detail, snp_anno, new_output=None):
        if self.option('method') == 'samtools':
            with open(snp_detail,'r') as f:
                for line in f:
                    f1 = open(line.strip())
                    column_line = f1.readline().strip()
                    if "bam" in column_line:
                        cmd = "sed -i \'1s/bam\///\' {}".format(f1)
                        os.system(cmd)
    
    # 加一列REF/ALT并保存至新文件，设置列顺序，为了和旧代码结果一致
    # def rearrange(snp_detail,new_output + "/data_anno_pre.xls"):
        with open(new_output + "/data_anno_pre.xls",'w') as w:
            with open(snp_detail,'r') as f:
                file1 = f.readline().strip()
                f1 = open(file1)
                column_line = f1.readline().strip().split('\t')
                # sam_nums = (len(column_line) - 12) / 2
                # sample_names = column_line[11:11 + sam_nums]
                # sample_genotypes = column_line[11 + sam_nums:-1]
                all_len = len(column_line)
                ss = 11 + int((all_len - 12) / 2)
                ssss = all_len - 1
                column_name_new = []
                for sample_name in column_line[11:ss]:
                    column_name_new.append(sample_name + '_sampledepth')
                    column_name_new.append(sample_name + 'genotype')
                w.write('\t'.join(column_name_new)+'\talt\tanno\tchrom\tend\tgene\tmut_info\tmut_type\tqual\treads_num\tref\tsnp_type\tstart\ttype'+'\n')
                column_types = {}
                for column in column_line:
                    if column in ['START', 'END', 'Depth']:
                        column_types[column] = "int64"
                    elif column in ['QUAL']:
                        column_types[column] = "float64"
                    else:
                        column_types[column] = "category"
                snp_detail_df = pd.read_table(file1, sep='\t', dtype=column_types)
                for line_ in f1:
                    line = line_.strip().split('\t')
                    a = line[11:ss]
                    b = line[ss:ssss]
                    c = list(chain(*zip(a, b)))
                    c.extend([line[4],line[6],line[0],line[2],line[7],line[10],line[9],line[5],line[8],line[3],line[3]+'/'+line[4],line[1],line[-1]])
                    w.write('\t'.join(c)+'\n')
                for line in f:
                    file2 = line.strip()
                    f2 = open(file2)
                    head = f2.readline()
                    for line_ in f2:
                        line = line_.strip().split('\t')
                        a = line[11:ss]
                        b = line[ss:ssss]
                        c = list(chain(*zip(a, b)))
                        c.extend(
                            [line[4], line[6], line[0], line[2], line[7], line[10], line[9], line[5], line[8], line[3],
                             line[3] + '/' + line[4], line[1], line[-1]])
                        w.write('\t'.join(c) + '\n')
                    df1 = pd.read_table(file2, sep='\t', dtype=column_types)
                    snp_detail_df = pd.concat([snp_detail_df, df1], axis=0)

        # 读取anno和chrom的unique
        anno_uni = snp_detail_df['ANNO'].unique()
        chrom = snp_detail_df['CHROM'].unique()
        with open(new_output + '/main_info','w') as w:
            w.write('snp_types\tA/G;A/C;C/T;G/A;G/C;C/A;A/T;C/G;G/T;T/C;T/A;T/G\n')
            w.write('specimen\t' + ';'.join(sample_names) + '\n')
            w.write('distributions\t' + ';'.join(anno_uni) + '\n')
            w.write('chroms\t' + ';'.join(chrom) + '\n')

        # 新增REF/ALT列
        snp_detail_df['snp_type'] = snp_detail_df['REF'].str.cat(snp_detail_df['ALT'],sep='/')
        snp_detail_df['snp_type'] = snp_detail_df['snp_type'].astype('category')

        # 列名重命名
        new_columns = {'CHROM':'chrom', 'START':'start', 'END':'end', 'REF':'ref', 'ALT':'alt',\
        'QUAL':'qual', 'ANNO':'anno', 'GENE(in or nearby)':'gene', 'Depth':'reads_num','MUT_type':'mut_type', 'MUT_info':'mut_info'}
        new_sample_names = {}
        for sample_name in sample_names:
            new_sample_names[sample_name] = sample_name + '_sampledepth'
        new_columns.update(new_sample_names)
        snp_detail_df.rename(columns=new_columns,inplace=True)

        # 设置列顺序，为了和旧代码结果一致
        new_column_name = []
        for sample_name in sample_names:
            new_column_name.append(sample_name + '_sampledepth')
            new_column_name.append(sample_name + 'genotype')
        new_column_name.extend(["alt","anno","chrom","end","gene","mut_info","mut_type","qual","reads_num","ref","snp_type","start","type"])
        snp_detail_df = snp_detail_df[new_column_name]

        # 加一列REF/ALT并保存至新文件
        # snp_detail_df.to_csv(new_output + "/data_anno_pre.xls", sep="\t", header=True, index=False)
        # snp_detail_df.to_pickle(new_output + "/data_anno_pre.pkl")
        # rearrange(snp_detail,new_output + "/data_anno_pre.xls")
                
        with open(snp_anno,'r') as f:
            file1 = f.readline()
            f1 = open(file1)
            column_line_y = f1.readline().strip().split('\t')
            column_types_y = {}
            for column in column_line_y[:13]:
                if column in ['START','END','Depth']:
                    column_types_y[column] = "int64"
                elif column in ['QUAL']:
                    column_types_y[column] = "float64"
                else:
                    column_types_y[column] = "category"
            df_anno_pre2 = pd.read_table(file1, header=0, sep="\t", dtype=column_types_y, usecols=column_line_y[:13])
            for line_ in f:
                file2 = line_.strip()
                df2 = pd.read_table(file2, header=0, sep="\t", dtype=column_types_y, usecols=column_line_y[:13])
                df_anno_pre2 = pd.concat([df_anno_pre2, df2], axis=0)


        df_anno_needed_columns = [item for item in snp_detail_df.columns if item not in set(["gene", "mut_info", "mut_type", "reads_num", "snp_type", "qual"])]
        df_anno_pre1 = snp_detail_df.loc[:,df_anno_needed_columns]
        # df_anno_pre1 = pd.read_pickle(x)
        # df_anno_pre1.drop(["gene", "mut_info", "mut_type", "reads_num", "snp_type", "qual"], axis=1, inplace=True)
        df_anno_pre1.rename(columns={"chrom": "Chrom", "alt": "Alt", "anno": "Anno", "end": "End", "start": "Start",\
            "ref": "Ref"}, inplace=True)

        # df_anno_pre2.drop(df_anno_pre2.columns[13:], axis=1, inplace=True)
        df_anno_pre2.rename(columns={"Depth": "Total depth", "CHROM": "Chrom", "ALT": "Alt", "ANNO": "Anno",\
            "END": "End", "START": "Start", "REF": "Ref","MUT_type": "MUT type", "MUT_info": "MUT info", "": ""}, inplace=True)
        df_anno_pre2 = df_anno_pre2[["GENE(in or nearby)", "Gene name", "Gene description", "Chrom", "Start", \
            "End", "Ref", "Alt", "Total depth", "QUAL","Anno", "MUT type", "MUT info"]]
        df_join_pre = pd.merge(df_anno_pre2, df_anno_pre1, on=['Alt', 'Anno', 'Chrom', 'End', 'Start', 'Ref'], how="outer")
        df_join_pre.to_csv(self.output_dir + "/snp_annotation_statistics.xls", sep="\t", index=False, chunksize=500000)


        #####################
        #####################
        def statistic_snp(sample_name):
            sample_new_name = str(sample_name + '_sampledepth')
            needed_columns = [sample_new_name,'snp_type','type','anno']
            drop_columns = [i for i in new_column_name if i not in needed_columns]
            sample_detail_df = snp_detail_df.drop(drop_columns, axis=1, inplace=False)
            # sample_detail_df.drop(sample_detail_df[sample_detail_df[sample_new_name].isin(["-",'.'])].index, inplace=True)
            sample_detail_df = sample_detail_df[~sample_detail_df[sample_new_name].isin(["-",'.'])]
            sample_detail_df = sample_detail_df[[sample_new_name,'snp_type','type','anno']]
            # 样本值相加 再替换原值
            new_sample_detail_df = sample_detail_df[sample_new_name].str.split(',', expand=True)
            new_sample_detail_df.fillna(0, inplace=True)
            for column_num in range(len(new_sample_detail_df.columns)):
                new_sample_detail_df.iloc[:,column_num] = pd.to_numeric(new_sample_detail_df.iloc[:,column_num], downcast='unsigned')
            sample_detail_df[sample_new_name] = new_sample_detail_df.sum(axis=1)
            # 统计 样本值区间分布
            cuts = pd.cut(sample_detail_df[sample_new_name], bins=[0, 30, 100, 200, 300, 400, 500, float("inf")], right=True, \
                      labels=['<=30', '31-100', '101-200', '201-300', '301-400', '401-500', '>500'])
            sample_detail_df_depth_count = sample_detail_df[sample_new_name].groupby(cuts).count().to_frame()
            sample_detail_df_depth_count.index.name = None
            sample_detail_df_depth_count.columns = [sample_name]
            sample_detail_df_depth_count.index = sample_detail_df_depth_count.index.tolist()   # categoricalindex to index
            # 统计 REF/ALT
            sample_detail_df.drop(sample_detail_df[sample_detail_df[sample_new_name] == 0].index, inplace=True)
                                                                                 # 设定index：既能 增加没有的index，也能 删除多余的index
            snp_type_count = sample_detail_df['snp_type'].value_counts().reindex(index=\
                    ["A/G", "A/C", "C/T", "G/A", "G/C", "C/A", "A/T", "C/G", "G/T", "T/C", "T/A", "T/G"], fill_value=0).astype('int').to_frame()
            snp_type_count.columns = [sample_name]
            # 统计 snp和indel 对应的anno
            sample_detail_df.drop([sample_new_name,'snp_type'], axis=1, inplace=True)
            grouped = sample_detail_df.groupby('type')
            snp_df = pd.Series()
            indel_df = pd.Series()
            for name,group in grouped:
                if name == 'snp':
                    snp_df = group['anno']
                else:
                    indel_df = indel_df.append(group['anno'])
            snp_df_count = snp_df.value_counts().to_frame()
            snp_df_count.columns = [sample_name]
            indel_df_count = indel_df.value_counts().to_frame()
            indel_df_count.columns = [sample_name]
            indel_df_count = indel_df_count.reindex(index=list(snp_df_count.index), fill_value=0)
            ####
            return sample_detail_df_depth_count,snp_type_count,snp_df_count,indel_df_count
            
        
        sample_detail_df_depth_count_all = pd.DataFrame()
        snp_type_count_all = pd.DataFrame()
        snp_df_count_all = pd.DataFrame()
        indel_df_count_all = pd.DataFrame()
        ##
        detail_dfs = []
        for sample_name in sample_names:
            sample_result = dask.delayed(statistic_snp)(sample_name)
            detail_dfs.append(sample_result)
        detail_dfs = dask.compute(*detail_dfs)
        for sample_result in detail_dfs:
            sample_detail_df_depth_count_all = pd.concat([sample_detail_df_depth_count_all, sample_result[0]], axis=1)
            snp_type_count_all = pd.concat([snp_type_count_all, sample_result[1]], axis=1)
            snp_df_count_all = pd.concat([snp_df_count_all, sample_result[2]], axis=1)
            indel_df_count_all = pd.concat([indel_df_count_all,sample_result[3]], axis=1)
        #########################

        ## 保存 snp_depth_statistics
        sample_new_name_list_depth_count = sample_names + ['range_key','stat_type']
        sample_detail_df_depth_count_all.index.name = 'range_key'
        sample_detail_df_depth_count_all = sample_detail_df_depth_count_all.reset_index()
        sample_detail_df_depth_count_all['stat_type'] = 'depth_stat'
        sample_detail_df_depth_count_all = sample_detail_df_depth_count_all[sample_new_name_list_depth_count]
        sample_detail_df_depth_count_all.to_csv(new_output + "/snp_depth_statistics.xls", sep="\t", index=False)
        ## 保存 snp_transition_tranversion_statistics
        snp_type_count_all.index.name = 'range_key'
        snp_type_count_all = snp_type_count_all.reset_index()
        snp_type_count_all['stat_type'] = 'type_stat'
        snp_type_count_all = snp_type_count_all[sample_new_name_list_depth_count]
        snp_type_count_all.to_csv(new_output + "/snp_transition_tranversion_statistics.xls", sep="\t", index=False)
        ## 保存 snp_position_distribution
        snp_df_count_all.index.name = 'range_key'
        snp_df_count_all = snp_df_count_all.reset_index()
        snp_df_count_all['stat_type'] = 'snp_distribution'
        snp_df_count_all = snp_df_count_all[sample_new_name_list_depth_count]
        snp_df_count_all.to_csv(new_output + "/snp_position_distribution.xls", sep="\t", index=False)
        ## 保存 indel_position_distribution
        indel_df_count_all.index.name = 'range_key'
        indel_df_count_all = indel_df_count_all.reset_index()
        indel_df_count_all['stat_type'] = 'indel_distribution'
        indel_df_count_all = indel_df_count_all[sample_new_name_list_depth_count]
        indel_df_count_all.to_csv(new_output + "/indel_position_distribution.xls", sep="\t", index=False)

    def run(self):
        """
        运行
        """
        super(PredealSnpresults2Tool, self).run()
        self.logger.info("对snp生成的结果文件进行统计等预处理")
        self.deal_snp_detail(self.option('snp_detail_list'), self.option('snp_anno_list'), new_output=self.output_dir)
        # if os.path.exists(os.path.join(self.output_dir, 'data_anno_pre.xls')):
        #     self.merge(os.path.join(self.output_dir, 'data_anno_pre.xls'), self.option('snp_anno'))
        #     # os.remove(os.path.join(self.output_dir, 'data_anno_pre.pkl'))
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        # test_dir='/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna'
        data = {
            "id": "snp" + str(random.randint(1, 10000)) + "_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v3.large.predeal_snpresults2",
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
                snp_detail_list="/mnt/ilustre/users/isanger/workspace/20210315/Refrna_n34u_49qekvdhgok1ed33c1m5hf/CallSnpIndelSentieon/snp_anno_list.txt",
                snp_anno_list="/mnt/ilustre/users/isanger/workspace/20210315/Refrna_n34u_49qekvdhgok1ed33c1m5hf/CallSnpIndelSentieon/snp_annotation_list.txt",
                method="gatk"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()





