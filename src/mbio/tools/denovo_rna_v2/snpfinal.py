# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
import re
from biocluster.core.exceptions import OptionError
import pandas as pd
from collections import OrderedDict

class SnpfinalAgent(Agent):
    """
    对无参转录组序列进行SSR的预测，以及引物设计,序列位置分布预测。
    检查tool的详细的错误可以去tool的结果里面一个err的文件查看
    """
    def __init__(self, parent):
        super(SnpfinalAgent, self).__init__(parent)
        """
        unigene_fa这个有写好的文件格式检查，对于页面的具体参数传递过来的形式再讨论，再看怎么写到misa.ini里面,以及这个默认值是怎么设置
        """
        options = [
            {"name": "bamlist", "type": "int"},
            {"name": "call_vcf", "type": "infile", "format": "denovo_rna_v2.vcf"},
            {"name": "qual", "type": "float", "default": 20},
            {"name": "dp", "type": "int", "default": 1},
            # {"name": "isoform_unigene", "type": "string"},#由于call.vcf需要上传，而且为了避免在重运行里面转换，所以在snp分析的第一个tool，进行了unigene的转换
        ]
        self.add_option(options)
        self.step.add_steps("SSR_predict")
        self.on("start", self.step_start)
        self.on("end", self.step_end)


    def step_start(self):
        self.step.SSR_predict.start()
        self.step.update()

    def step_end(self):
        self.step.SSR_predict.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数

        :return:
        """
        if not self.option("call_vcf"):
            raise OptionError("必须设置输入文件:比对SSR位置的参考序列call_vcf文件", code = "32005601")


    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = "20G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_regexp_rules([
            [".", "", "结果输出目录"],
            ["snp_detail", " ", "snp详情表数据"],
            ["indel_detail", " ", "indel详情表数据"],
            ["snp_depth_statistics", " ", "SNP测序深度统计表  "],
            ["snp_homo_hete_statistics", " ", "SNP类型统计表  "],
            ["snp_transition_tranversion_statistics", " ", "SNP位点统计表"]
        ])
        super(SnpfinalAgent, self).end()


class SnpfinalTool(Tool):
    def __init__(self, config):
        super(SnpfinalTool, self).__init__(config)
        #self.bcftools_path = Config().SOFTWARE_DIR + "/bioinfo/seq/bcftools-1.4/bin/bcftools"
        self.bcftools_path = Config().SOFTWARE_DIR + "/bioinfo/align/samtools-1.6/bcftools-1.6/bcftools"
        self.script_path = 'bioinfo/rna/scripts'
        self.query_path = self.config.PACKAGE_DIR + "/denovo_rna_v2/bcftools_filter.sh"

    def bcffilter(self):
        # 写入的文件不要有空行，这里DP的取值是样本数的个数，默认为1，QUAL默认为20
        # 写在sh文件的引号要用双引号，这样才会正确解析命令，并且相应的参数都写在sh里面
        # with open(self.option("bamlist").prop["path"], "r") as len_bam:
        #     sample_num = len(len_bam.readlines())
        sample_num = self.option("bamlist")
        # cmd1 = "{} {} {} {} {} {}".format(self.script_path + '/bcftools_filter.sh', self.bcftools_path, self.option('qual'), self.option('dp')*sample_num, self.work_dir + '/filter_vcf', self.option("call_vcf").prop["path"])
        cmd1 = "{} {} {} {} {} {} {}".format("/sh", self.query_path, self.bcftools_path, self.option('qual'), self.option('dp')*sample_num, self.work_dir + '/filter_vcf', self.option("call_vcf").prop["path"])
        self.logger.info("运行bcffilter")
        self.logger.info(cmd1)
        cmd1_obj = self.add_command("cmd1", cmd1)
        cmd1_obj.software_dir = "/bin"
        cmd1_obj.run()
        self.wait(cmd1_obj)
        if cmd1_obj.return_code == 0:
            self.logger.info("运行bcffilter完成")
        else:
            cmd1_obj.rerun()
            self.wait(cmd1_obj)
            if cmd1_obj.return_code == 0:
                self.logger.info("运行bcffilter完成")
            else:
                self.set_error("运行bcffilter出错", code = "32005602")

    def splitvcf(self):
        # 将得到的filter.vcf分割为snp和indel，并且第二个是需要按照统计详情，进行写入的文件路径
        with open(self.work_dir + '/filter_vcf', "r") as filter, open(self.work_dir + "/indel", "w") as indel, open(self.work_dir + "/snp", "w") as snp:
            header_list = list()
            for line in filter:
                if len(line.split("#")) > 2 or "LowQual" in line:
                    continue
                if "CHROM" in line:
                    header_list = line.split("\t")
                    header_list = [x if "/" not in x else x.split("/")[-1].split(".")[0] for x in header_list]
                    indel.write("\t".join(tuple(header_list)) + "\n")
                    snp.write("\t".join(tuple(header_list)) + "\n")
                if "INDEL" in line:
                    if sum([int(x) for x in line.split("\t")[7].split("DP4=")[1].split(";MQ")[0].split(",")]) >= self.option('dp')*self.option("bamlist"):
                        indel.write(line)
                if "INDEL" not in line and "DP=" in line:
                    if sum([int(x) for x in line.split("\t")[7].split("DP4=")[1].split(";MQ")[0].split(",")]) >= self.option('dp')*self.option("bamlist"):
                        # self.logger.info("123456789")
                        snp.write(line)

    def snpindel_detail(self, r_path, w_path, new_w_path = None):
        # 第一个参数是从初步得到的vcf文件中拆分出来的snp或者indel路径，是需要读的文件，第二个是需要按照统计详情，进行写入的文件路径,会出现call不出来snp或者indel的情况，所以需要加一个文件行数的判断
        # 如果大于1行，才证明真正的有文件信息写入（除去表头）
        with open(r_path, "r") as r, open(w_path, "w") as w:
            head_indel = r.readline()
            head_indel_list = head_indel.strip("\n").split("\t")
            # print head_indel_list
            # 单独处理第二行，把列名按照规律添加进去
            _ = r.readline()
            line = r.readline()
            head = ['CHROM', 'POS', 'REF', 'ALT']
            # w.write("\t".join(['CHROM', 'POS', 'ID', 'REF', 'ALT']))
            count = len(open(r_path,'rU').readlines())
            if count > 2:
                geneid, pos, ref, alt = line.split("\t")[0], line.split("\t")[1], line.split("\t")[3], line.split("\t")[4]
                line_2_data = [geneid, pos, ref, alt]
                # w.write("\t".join([geneid, pos, ref, alt]) + "\t")
                # enumerate可以返回索引值和位置
                for ind, sample in enumerate(line.split("\t")[9:], start=9):
                    gt, dp, ad = sample.split(":")[0], sample.split(":")[2], sample.split(":")[3]
                    if gt.split("/") == ['0', '0']:
                        line_2_data += [ref, dp, ad, "Homo"]
                        #　这里如果是写成head += head_indel_list[ind]这就会把样品名给拆开了，因为这是字符串，字符串和列表相加会导致字符串被拆成一个个字母的列表元素
                        head.append(head_indel_list[ind])
                        head += [head_indel_list[ind] + '_' + x for x in ('depth', 'alledepth', 'type')]

                    elif gt.split("/") == ['1', '1']:
                        line_2_data += [alt, dp, ad, "Homo"]
                        head.append(head_indel_list[ind])
                        head += [head_indel_list[ind] + '_' + x for x in ('depth', 'alledepth', 'type')]

                    elif gt.split("/") == ['0', '1']:
                        line_2_data += [",".join([ref, alt]),dp, ad, "Hete"]
                        head.append(head_indel_list[ind])
                        head += [head_indel_list[ind] + '_' + x for x in ('depth', 'alledepth', 'type')]
                        # w.write("\t".join([",".join([ref, alt]), dp, ad, "Hete"])+ "\t")

                    elif gt.split("/") == ['1', '2']:
                        if ad.split(",").count('0') == 2:
                            line_2_data += [ref, dp, ad, "Homo"]
                            head.append(head_indel_list[ind])
                            head += [head_indel_list[ind] + '_' + x for x in ('depth', 'alledepth', 'type')]
                            # w.write("\t".join([alt, dp, ad, "Homo"])+ "\t")
                        else:
                            line_2_data += [ref, dp, ad, "Hete"]
                            head.append(head_indel_list[ind])
                            head += [head_indel_list[ind] + '_' + x for x in ('depth', 'alledepth', 'type')]
                            # w.write("\t".join([alt, dp, ad, "Hete"])+ "\t")

                    elif gt.split("/") == ['0', '2']:
                        line_2_data += [",".join([ref, alt]),dp, ad, "Hete"]
                        head.append(head_indel_list[ind])
                        head += [head_indel_list[ind] + '_' + x for x in ('depth', 'alledepth', 'type')]
                        # w.write("\t".join([",".join([ref, alt]), dp, ad, "Hete"])+ "\t")

                    else:
                        line_2_data += ["-", "-", "-" ,"-"]
                        head.append(head_indel_list[ind])
                        head += [head_indel_list[ind] + '_' + x for x in ('depth', 'alledepth', 'type')]
                        # w.write("\t".join(["-", "-", "-" ,"-"]) + "\t")
                w.write("\t".join(head) + "\n")
                w.write("\t".join(line_2_data) + "\n")

                for line in r:
                    line = line.strip("\n")
                    geneid, pos, ref, alt = line.split("\t")[0], line.split("\t")[1], line.split("\t")[3], line.split("\t")[4]
                    w.write("\t".join([geneid, pos, ref, alt]) + "\t")
                    for sample in line.split("\t")[9:]:
                        gt, dp, ad = sample.split(":")[0], sample.split(":")[2], sample.split(":")[3]
                        if gt.split("/") == ['0', '0']:
                            w.write("\t".join([ref, dp, ad, "Homo"]) + "\t")

                        elif gt.split("/") == ['1', '1']:
                            w.write("\t".join([alt, dp, ad, "Homo"]) + "\t")

                        elif gt.split("/") == ['0', '1']:
                            w.write("\t".join([",".join([ref, alt]), dp, ad, "Hete"]) + "\t")

                        elif gt.split("/") == ['1', '2']:
                            if ad.split(",").count('0') == 2:
                                w.write("\t".join([alt, dp, ad, "Homo"]) + "\t")
                            else:
                                w.write("\t".join([alt, dp, ad, "Hete"]) + "\t")

                        elif gt.split("/") == ['0', '2']:
                            w.write("\t".join([",".join([ref, alt]), dp, ad, "Hete"]) + "\t")

                        else:
                            w.write("\t".join(["-", "-", "-" ,"-"]) + "\t")
                    w.write("\n")

    def rewrite(self, file, new_w_path):
        #　去除可能文件末尾存在的\t，所以我们需要重新读取一次文件，要不然后面pandas读取文件，总会发现某些行的个数多出一列, 这个函数调用2次，分别计算snp和indel的详情情况
        with open(file, "r") as new_r_path, open(new_w_path, "w") as nwp:
            for line in new_r_path:
                line = line.strip()
                nwp.write(line + "\n")


    def snpdepth(self, s4):
        pd_snp = pd.read_table(s4, header = 0, sep = "\t")
        snp_col_len = pd_snp.shape[1]
        select_col = range(5, snp_col_len, 4)
        dp_selected = pd_snp.iloc[:,select_col].replace([r'-'], -1)
        column_names = dp_selected.columns
        num30, num31_100, num101_200, num201_300, num301_400, num401_500, num501 = 0, 0, 0, 0, 0, 0, 0
        dp_selected = dp_selected.astype('int64')
        df = pd.DataFrame()

        df['depth'] = ["<=30", "31-100", "101-200", "201-300", "301-400", "401-500", ">500"]
        for i in column_names:
            num30 = dp_selected[dp_selected[i].between(0, 30, inclusive=True)].shape[0]
            num31_100 = dp_selected[dp_selected[i].between(31, 1000, inclusive=True)].shape[0]
            num101_200 = dp_selected[dp_selected[i].between(101, 200, inclusive=True)].shape[0]
            num201_300 = dp_selected[dp_selected[i].between(201, 300, inclusive=True)].shape[0]
            num301_400 = dp_selected[dp_selected[i].between(301, 400, inclusive=True)].shape[0]
            num401_500 = dp_selected[dp_selected[i].between(401, 500, inclusive=True)].shape[0]
            num501 = dp_selected[dp_selected[i] >= 501].shape[0]
            df[i] = pd.DataFrame([num30, num31_100, num101_200, num201_300, num301_400, num401_500, num501])
            pd.concat([df, df[i]])
        df.to_csv(self.work_dir + "/depth",sep = "\t")
        tt_per = df.iloc[:,1:]/df.iloc[:,1:].sum()
        tt_per.columns = [x + "_per" for x in tt_per.columns]
        tt_new_per = pd.concat([df, tt_per], axis = 1)
        tt_new_per = tt_new_per.round(6)
        tt_new_per.to_csv(self.work_dir + "/depth_new_per",sep = "\t", index=False)

    # 统计snp的homo和hete的个数
    def snphh(self, s4):
        pd_snp = pd.read_table(s4, header = 0, sep = "\t")
        snp_col_len = pd_snp.shape[1]
        select_col = range(7, snp_col_len, 4)
        dp_selected = pd_snp.iloc[:,select_col].replace([r'-'], 'NA')
        column_names = dp_selected.columns
        columns = ["HomoSNP", "HeteSNP", "Total"]
        index = column_names
        df = pd.DataFrame(columns = columns, index = index)
        for i in column_names:
            num_homo = pd_snp[pd_snp[i].isin(['Homo'])].shape[0]
            num_hete = pd_snp[pd_snp[i].isin(['Hete'])].shape[0]
            total = num_hete + num_homo
            df.loc[i].HomoSNP = num_homo
            df.loc[i].HeteSNP = num_hete
            df.loc[i].Total = total
        df.index = [x.replace("_type", "") for x in column_names]
        df.index.name = 'sample'
        df.to_csv(self.work_dir + "/statis_hh", sep = "\t")

    # 统计snp的transition和tranversion
    def snptt(self, s4):
        with open(s4, "r") as s4:
            head_indel = s4.readline()
            head_indel_list = head_indel.strip("\n").split("\t")[4::4]
            # enumerate可以返回索引值和位置
            targets = [x.strip() for x in "AG,CT, GA, TC, AC, GC, CA, AT, CG, GT, TA, TG".split(',')]
            """如果写成下面这个样子，因为这样建造的字典的都是同一个，会得到所有的样本统计的个数都是一样，这是因为所有的字典都是同一个，会产生联动，下面正确的写法嵌套2层循环
            ，就会发现结果不变
            target_dict = {x:0 for x in targets}
            data_dict = dict(zip(head_indel_list, [target_dict for _ in head_indel_list]))
            """
            data_dict = OrderedDict(zip(head_indel_list, [{x:0 for x in targets} for _ in head_indel_list]))

            for line in s4:
                ref = line.split("\t")[2]
                tmp_list = line.split("\t")[4::4]

                for ind, sample in enumerate(tmp_list):
                    if ref == 'A':
                        if 'G' in sample.split(","):
                            data_dict[head_indel_list[ind]]['AG'] += 1
                        elif 'C' in sample.split(","):
                            data_dict[head_indel_list[ind]]['AC'] += 1

                        elif 'T' in sample.split(","):
                            data_dict[head_indel_list[ind]]['AT'] += 1

                    if ref == 'C':
                        if 'T' in sample.split(","):
                            data_dict[head_indel_list[ind]]['CT'] += 1

                        elif 'A' in sample.split(","):
                            data_dict[head_indel_list[ind]]['CA'] += 1

                        elif 'G' in sample.split(","):
                            data_dict[head_indel_list[ind]]['CG'] += 1

                    if ref == 'G':
                        if 'A' in sample.split(","):
                            data_dict[head_indel_list[ind]]['GA'] += 1

                        elif 'C' in sample.split(","):
                            data_dict[head_indel_list[ind]]['GC'] += 1

                        elif 'T' in sample.split(","):
                            data_dict[head_indel_list[ind]]['GT'] += 1

                    if ref == 'T':
                        if 'C' in sample.split(","):
                            data_dict[head_indel_list[ind]]['TC'] += 1

                        elif 'A' in sample.split(","):
                            data_dict[head_indel_list[ind]]['TA'] += 1

                        elif 'G' in sample.split(","):
                            data_dict[head_indel_list[ind]]['TG'] += 1

            result = pd.DataFrame(data_dict)
            result.index.name = 'type'
            result.to_csv(self.work_dir + "/transition_tranversion", sep = "\t")
            df_tt =  pd.read_table(self.work_dir + "/transition_tranversion", sep = "\t", index_col = 0)
            df_tt.index.name = 'type'
            # df_tt.index.name = 'Type'
            tt_per = df_tt/df_tt.sum()
            tt_per.columns = [x + "_per" for x in tt_per.columns]
            tt_new_per = pd.concat([df_tt, tt_per], axis = 1)
            tt_new_per = tt_new_per.round(6)
            tt_new_per.to_csv(self.work_dir + "/tt_new_per",sep = "\t")

    def run(self):
        super(SnpfinalTool, self).run()
        self.bcffilter()
        self.splitvcf()
        self.snpindel_detail(self.work_dir + "/indel" , self.work_dir + "/indel_new")
        if os.path.getsize(self.work_dir + "/indel_new") > 0:
            self.rewrite(self.work_dir + "/indel_new", self.work_dir + "/indel_rewrite")
            new_indel_rewrite =  pd.read_table(self.work_dir + "/indel_rewrite", sep = "\t", header = 0)
            new_indel_rewrite.insert(4, "REF_ALT", new_indel_rewrite["REF"] + new_indel_rewrite["ALT"])
            new_indel_rewrite.to_csv(self.work_dir + "/new_indel_rewrite",sep = "\t")

            # df_1 = pd.read_table('new_snp_rewrite', sep ="\t", header = 0)
            # df_2 = pd.read_table(self.option("isoform_unigene"), header=None, sep = '\t')
            # df_2.columns = ["CHROM", "unigene", "yes_no"]
            # merge_df12 = pd.merge(df_2, new_indel_rewrite, on="CHROM")
            # merge_df12_new = merge_df12[merge_df12.yes_no.isin(['yes'])]
            # merge_df12_new.drop(['CHROM', 'yes_no'], axis=1,inplace=True)
            # merge_df12_new.rename(columns={merge_df12_new.columns[0]: "CHROM"},inplace=True)
            # merge_df12_new.to_csv(self.work_dir + "/new_indel_rewrite",sep = "\t", index=False)

        self.snpindel_detail(self.work_dir + "/snp", self.work_dir + "/snp_new")
        if os.path.getsize(self.work_dir + "/snp_new") > 0:
            self.rewrite(self.work_dir + "/snp_new", self.work_dir + "/snp_rewrite")
            self.snpdepth(self.work_dir + "/snp_rewrite")
            self.snphh(self.work_dir + "/snp_rewrite")
            self.snptt(self.work_dir + "/snp_rewrite")
            new_snp_rewrite =  pd.read_table(self.work_dir + "/snp_rewrite", sep = "\t", header = 0)
            new_snp_rewrite.insert(4, "REF_ALT", new_snp_rewrite["REF"] + '/' + new_snp_rewrite["ALT"])
            new_snp_rewrite.to_csv(self.work_dir + "/new_snp_rewrite",sep = "\t", index = False)
            #
            # df_2 = pd.read_table(self.option("isoform_unigene"), header=None, sep = '\t')
            # df_2.columns = ["CHROM", "unigene", "yes_no"]
            # merge_df12 = pd.merge(df_2, new_snp_rewrite, on="CHROM")
            # merge_df12_new = merge_df12[merge_df12.yes_no.isin(['yes'])]
            # merge_df12_new.drop(['CHROM', 'yes_no'], axis=1,inplace=True)
            # merge_df12_new.rename(columns={merge_df12_new.columns[0]: "CHROM"},inplace=True)
            # merge_df12_new.to_csv(self.work_dir + "/new_snp_rewrite",sep = "\t", index=False)
        self.set_output()
        self.end()

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))

        filter_vcf = self.work_dir + '/filter_vcf'
        os.link(filter_vcf, os.path.join(self.output_dir, 'filter_vcf'))
        if os.path.getsize(self.work_dir + "/snp_new") > 0:
            snp_rewrite = self.work_dir + '/new_snp_rewrite'
            new_snp_rewrite = self.work_dir + '/new_snp_rewrite'
            depth = self.work_dir + '/depth'
            depth_new_per = self.work_dir + '/depth_new_per'
            tt_new_per = self.work_dir + '/tt_new_per'
            transition_tranversion = self.work_dir + '/transition_tranversion'
            statis_hh = self.work_dir + '/statis_hh'
            os.link(snp_rewrite, os.path.join(self.output_dir, 'snp_detail'))
            os.link(depth, os.path.join(self.output_dir, 'snp_depth_statistics'))
            os.link(statis_hh, os.path.join(self.output_dir, 'snp_homo_hete_statistics'))
            os.link(transition_tranversion, os.path.join(self.output_dir, 'snp_transition_tranversion_statistics'))
            self.logger.info("成功设置SNP结果目录")
        else:
            self.logger.info("此次分析没有call得到snp信息，设置snp_detail结果目录失败")

        if os.path.getsize(self.work_dir + "/indel_new") > 0:
            indle_rewrite = self.work_dir + '/new_indel_rewrite'
            os.link(indle_rewrite, os.path.join(self.output_dir, 'indel_detail'))
        else:
            self.logger.info("此次分析没有call得到indel信息，设置indel_detail结果目录失败")
