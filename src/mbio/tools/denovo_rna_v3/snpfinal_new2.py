# -*- coding: utf-8 -*-
# __author__ = 'litangjian ; xuxi improved pandas and dask in 20201026 '

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
import dask
from biocluster.core.exceptions import OptionError
import pandas as pd
import unittest


class SnpfinalNew2Agent(Agent):
    """
    对无参转录组序列进行SSR的预测，以及引物设计,序列位置分布预测。
    检查tool的详细的错误可以去tool的结果里面一个err的文件查看
    """

    def __init__(self, parent):
        super(SnpfinalNew2Agent, self).__init__(parent)
        """
        unigene_fa这个有写好的文件格式检查，对于页面的具体参数传递过来的形式再讨论，再看怎么写到misa.ini里面,以及这个默认值是怎么设置
        """
        options = [
            {"name": "bamlist", "type": "int"},
            {"name": "method", "type": "string", "default": "sentieon"},
            {"name": "call_vcf", "type": "infile", "format": "denovo_rna_v2.vcf"},
            {"name": "filted_snp_vcf", "type": "infile", "format": "denovo_rna_v2.vcf"},
            {"name": "filted_indel_vcf", "type": "infile", "format": "denovo_rna_v2.vcf"},
            {"name": "cds_bed", "type": "infile", "format": "denovo_rna_v2.common"},
            {"name": "allt2g", "type": "infile", "format": "denovo_rna_v2.common"},
            {"name": "anno", "type": "infile", "format": "denovo_rna_v2.common"},
        ]
        self.add_option(options)
        self.step.add_steps("SSR_predict")
        self.on("start", self.step_start)
        self.on("end", self.step_end)
        self._memory_increase_step = 100

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
            raise OptionError("必须设置输入文件:比对SSR位置的参考序列call_vcf文件", code="300101")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 20
        self._memory = "150G"

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
        super(SnpfinalNew2Agent, self).end()


class SnpfinalNew2Tool(Tool):
    def __init__(self, config):
        super(SnpfinalNew2Tool, self).__init__(config)
        self.bcftools_path = Config().SOFTWARE_DIR + "/bioinfo/align/samtools-1.6/bcftools-1.6/bcftools"
        self.script_path = 'bioinfo/rna/scripts'
        self.query_path = self.config.PACKAGE_DIR + "/denovo_rna_v2/bcftools_filter.sh"
        self.t2g_dict = {}
        self.g_pos = {}

    def splitvcf(self):
        if self.option("method") == "samtools":
            # 将得到的filter.vcf分割为snp和indel，并且第二个是需要按照统计详情，进行写入的文件路径
            with open(self.option("call_vcf").prop["path"], "r") as filter, \
                    open(self.work_dir + "/indel", "w") as indel, \
                    open(self.work_dir + "/snp", "w") as snp:
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
                        if sum([int(x) for x in
                                line.split("\t")[7].split("DP4=")[1].split(";MQ")[0].split(",")]) >= 1 * self.option(
                                "bamlist"):
                            indel.write(line)
                    if "INDEL" not in line and "DP=" in line:
                        if sum([int(x) for x in
                                line.split("\t")[7].split("DP4=")[1].split(";MQ")[0].split(",")]) >= 1 * self.option(
                                "bamlist"):
                            # self.logger.info("123456789")
                            snp.write(line)
        else:
            with open(self.option("filted_snp_vcf").prop["path"], "r") as filter_snp, \
                    open(self.option("filted_indel_vcf").prop["path"], "r") as filter_indel, \
                    open(self.work_dir + "/indel","w") as indel, open(self.work_dir + "/snp", "w") as snp:
                header_list = list()
                for line in filter_snp:
                    if len(line.split("#")) > 2 or "LowQual" in line:
                        continue
                    if "CHROM" in line:
                        header_list = line.split("\t")
                        header_list = [x if "/" not in x else x.split("/")[-1].split(".")[0] for x in header_list]
                        snp.write("\t".join(tuple(header_list)) + "\n")
                    else:
                        snp.write(line)
                for line in filter_indel:
                    if len(line.split("#")) > 2 or "LowQual" in line:
                        continue
                    if "CHROM" in line:
                        header_list = line.split("\t")
                        header_list = [x if "/" not in x else x.split("/")[-1].split(".")[0] for x in header_list]
                        indel.write("\t".join(tuple(header_list)) + "\n")
                    else:
                        indel.write(line)

    def bedpre(self):
        with open(self.option("allt2g").prop["path"]) as t2g:
            for line in t2g.readlines():
                line = line.strip().split("\t")
                if line[2] == "yes":
                    self.t2g_dict[line[0]] = line[1]
        with open(self.option("cds_bed").prop["path"], "r") as bd, open(os.path.join(self.work_dir, "cds_info"),
                                                                        "w") as ci:
            ci.write("CHROM\tstart\tend\n")
            for line in bd.readlines():
                line = line.strip().split("\t")
                if line[0] in self.t2g_dict:
                    ci.write(self.t2g_dict[line[0]] + "\t" + line[6] + "\t" + line[7] + "\n")
                    self.g_pos[self.t2g_dict[line[0]]] = [line[6], line[7]]

    def snpindel_detail(self, r_path, w_path, new_w_path=None):
        # 第一个参数是从初步得到的vcf文件中拆分出来的snp或者indel路径，是需要读的文件，第二个是需要按照统计详情，进行写入的文件路径,会出现call不出来snp或者indel的情况，所以需要加一个文件行数的判断
        # 如果大于1行，才证明真正的有文件信息写入（除去表头）
        with open(r_path, "r") as r, open(w_path, "w") as w:
            head_indel = r.readline()
            head_indel_list = head_indel.strip("\n").split("\t")
            # print head_indel_list
            # 单独处理第二行，把列名按照规律添加进去
            _ = r.readline()
            # line = r.readline()
            head = ['CHROM', 'POS', 'REF', 'ALT']
            for i in head_indel_list[9:]:
                head.append(i)
                head += [i + '_' + x for x in ('depth', 'alledepth', 'type')]
            w.write("\t".join(head) + "\tanno\n")
            count = len(open(r_path, 'rU').readlines())
            if count > 2:
                for line in r:
                    line = line.strip("\n")
                    geneid, pos, ref, alt = line.split("\t")[0], line.split("\t")[1], line.split("\t")[3], \
                                            line.split("\t")[4]
                    alt_infos = alt.split(",")
                    w.write("\t".join([geneid, pos, ref, alt]) + "\t")
                    for sample in line.split("\t")[9:]:
                        if len(sample.split(":")) == 6:
                            gt, dp, ad = sample.split(":")[0], sample.split(":")[2], sample.split(":")[3]
                        elif len(sample.split(":")) == 5:
                            gt, dp, ad = sample.split(":")[0], sample.split(":")[2], sample.split(":")[1]
                        elif len(sample.split(":")) == 3:
                            gt, dp, ad = sample.split(":")[0], sample.split(":")[1], "-"
                        elif len(sample.split(":")) == 7:
                            gt, dp, ad = sample.split(":")[0], sample.split(":")[2], sample.split(":")[1]
                        else:
                            gt, dp, ad = "-", "-", "-"
                        if gt.split("/") == ['0', '0']:
                            w.write("\t".join([ref, dp, ad, "Homo"]) + "\t")

                        elif gt.split("/") == ['1', '1']:
                            w.write("\t".join([alt_infos[0], dp, ad, "Homo"]) + "\t")

                        elif gt.split("/") == ['2', '2']:
                            w.write("\t".join([alt_infos[1], dp, ad, "Homo"]) + "\t")

                        elif gt.split("/") == ['0', '1']:
                            w.write("\t".join([",".join([ref, alt_infos[0]]), dp, ad, "Hete"]) + "\t")

                        elif gt.split("/") == ['1', '2']:
                            w.write("\t".join([",".join([alt_infos[0], alt_infos[1]]), dp, ad, "Hete"]) + "\t")

                        elif gt.split("/") == ['0', '2']:
                            w.write("\t".join([",".join([ref, alt_infos[1]]), dp, ad, "Hete"]) + "\t")

                        elif gt.split("/") == ['0', '3']:
                            w.write("\t".join([",".join([ref, alt_infos[2]]), dp, ad, "Hete"]) + "\t")

                        elif gt.split("/") == ['1', '3']:
                            w.write("\t".join([",".join([alt_infos[0], alt_infos[2]]), dp, ad, "Hete"]) + "\t")

                        elif gt.split("/") == ['2', '3']:
                            w.write("\t".join([",".join([alt_infos[1], alt_infos[2]]), dp, ad, "Hete"]) + "\t")

                        elif gt.split("/") == ['3', '3']:
                            w.write("\t".join([alt_infos[2], dp, ad, "Homo"]) + "\t")
                        else:
                            w.write("\t".join(["-", "-", "-", "-"]) + "\t")
                    if geneid in self.g_pos:
                        if int(self.g_pos[geneid][0]) <= int(pos) <= int(self.g_pos[geneid][1]):
                            w.write("CDS")
                        else:
                            w.write("non-CDS")
                    else:
                        w.write("non-CDS")
                    w.write("\n")

    def rewrite(self, file, new_w_path):
        # 　去除可能文件末尾存在的\t，所以我们需要重新读取一次文件，要不然后面pandas读取文件，总会发现某些行的个数多出一列, 这个函数调用2次，分别计算snp和indel的详情情况
        with open(file, "r") as new_r_path, open(new_w_path, "w") as nwp:
            for line in new_r_path:
                line = line.strip()
                nwp.write(line + "\n")

    def snp_analyse(self, s4):
        def snp_sample(i):
            column_name = list(pd_snp_columns)[i]
            column_depth_name = list(pd_snp_columns)[i + 1]
            ## main function area run for longest time
            snp_col_len_list = list(range(snp_col_len))
            t = [snp_col_len_list.remove(i) for i in [0, 2, i, i + 1, i + 3, snp_col_len - 1]]
            sample_df = pd_snp.drop(pd_snp_columns[snp_col_len_list], axis=1, inplace=False)
            sample_df = sample_df.where(
                (sample_df[column_name] != '-') & (sample_df[column_name] != sample_df.REF)).dropna()
            ## main function area run for longest time
            # snp depth
            sample_df_depth = sample_df.iloc[:, 3].astype('int')
            cuts = pd.cut(sample_df_depth, bins=[-1, 30, 100, 200, 300, 400, 500, float("inf")], right=True, \
                          labels=['<=30', '31-100', '101-200', '201-300', '301-400', '401-500', '>500'])
            sample_df_depth_count = sample_df_depth.groupby(cuts).count().to_frame()
            sample_df_depth_count.index.name = None
            sample_df_depth_count.columns = [column_depth_name]
            # snp type
            sample_df_type = sample_df.iloc[:, 4].astype('object')
            sample_df_type_statis = sample_df_type.value_counts().reindex(index=['Homo', 'Hete'], fill_value=0).astype(
                'int').to_frame()
            sample_df_type_statis.index.name = None
            sample_df_type_statis.columns = [column_name]
            # snp cds
            sample_df_cds = sample_df.iloc[:, 5].astype('object')  # 只展示这指定的两个索引的行
            sample_df_cds_statis = sample_df_cds.value_counts().reindex(index=['CDS', 'non-CDS'], fill_value=0).astype(
                'int').to_frame()
            sample_df_cds_statis.index.name = None
            sample_df_cds_statis.columns = [column_name]
            # snp transition
            sample_df_transition = sample_df.iloc[:, [1, 2]].astype('object')
            sample_df_transition_split_df = \
                sample_df_transition[column_name].str.split(',', expand=True).stack().reset_index(level=1,
                                                                                                  drop=True).rename(
                    column_name)
            sample_df_transition_split_df_merge = sample_df_transition.drop(column_name, axis=1).join(
                sample_df_transition_split_df)
            sample_df_transition_split_df_merge.rename(columns={column_name: 'sample'}, inplace=True)
            sample_df_transition_split_df_merge[column_name] = sample_df_transition_split_df_merge['REF'] + \
                                                               sample_df_transition_split_df_merge['sample']
            snp_transition = sample_df_transition_split_df_merge[column_name].value_counts().reindex(index \
                                                                                                         =["AC", "AG",
                                                                                                           "AT", "CA",
                                                                                                           "CG", "CT",
                                                                                                           "GA", "GC",
                                                                                                           "GT", "TA",
                                                                                                           "TC", "TG"],
                                                                                                     fill_value=0).astype(
                'int').to_frame()
            # snp annotation
            sample_df_anno = sample_df.iloc[:, 0].astype('object')
            sample_df_anno.drop_duplicates(inplace=True)
            pd_anno_sample = pd_anno.reindex(index=sample_df_anno.tolist())
            snp_anno = pd_anno_sample.loc[:, ["go", "KO_id", "cog", "swissprot", "nr"]].dropna(axis=0, how='all')
            snp_anno_count = snp_anno.count(axis=0).to_frame()
            snp_anno_count.columns = [column_name]
            snp_anno_count = snp_anno_count.T
            snp_anno_count['genes'] = str(sample_df_anno.shape[0])
            snp_anno_count['totalre'] = str(snp_anno.shape[0])
            snp_anno_count.rename(
                columns={'go': 'gore', 'KO_id': 'keggre', 'cog': 'cogre', 'swissprot': 'swissre', 'nr': 'nrre'},
                inplace=True)
            snp_anno_count = snp_anno_count[["genes", "totalre", "gore", "keggre", "cogre", "swissre", "nrre"]]
            return sample_df_depth_count, column_depth_name, sample_df_type_statis, column_name, sample_df_cds_statis, snp_transition, snp_anno_count

        df_depth = pd.DataFrame()
        df_type = pd.DataFrame()
        df_cds = pd.DataFrame()
        df_transition = pd.DataFrame()
        df_anno = pd.DataFrame()
        pd_snp = pd.read_table(s4, header=0, sep="\t")
        pd_anno = pd.read_table(self.option("anno").prop["path"], header=0)
        pd_anno = pd_anno[pd_anno["is_gene"] == "yes"].drop(columns=['transcript', 'is_gene']).set_index('gene_id')
        snp_col_len = pd_snp.shape[1]
        pd_snp_columns = pd_snp.columns
        detail_dfs = []
        for i in list(range(4, snp_col_len - 1, 4)):
            sample_result = dask.delayed(snp_sample)(i)
            detail_dfs.append(sample_result)
        detail_dfs = dask.compute(*detail_dfs)
        for sample_result in detail_dfs:
            df_depth = pd.concat([df_depth, sample_result[0][sample_result[1]]], axis=1)
            df_type = pd.concat([df_type, sample_result[2][sample_result[3]]], axis=1)
            df_cds = pd.concat([df_cds, sample_result[4][sample_result[3]]], axis=1)
            df_transition = pd.concat([df_transition, sample_result[5][sample_result[3]]], axis=1)
            df_anno = pd.concat([df_anno, sample_result[6]], axis=0)

        tuples = [('0', '<=30'), ('1', '31-100'), ('2', '101-200'), ('3', '201-300'), ('4', '301-400'),
                  ('5', '401-500'), ('6', '>500')]
        mul_index = pd.MultiIndex.from_tuples(tuples, names=('', 'depth'))
        df_depth.index = mul_index
        df_depth.to_csv(self.work_dir + "/depth", sep="\t")
        tt_per = df_depth.iloc[:, :] / df_depth.iloc[:, :].sum()
        tt_per.columns = [x + "_per" for x in tt_per.columns]
        tt_new_per = pd.concat([df_depth, tt_per], axis=1)
        tt_new_per.index = [indexx[1] for indexx in mul_index]
        tt_new_per.index.names = ['depth']
        tt_new_per = tt_new_per.round(6)
        tt_new_per.to_csv(self.work_dir + "/depth_new_per", sep="\t", index=True)

        df_type = df_type.T[['Homo', 'Hete']]
        df_type.index.name = 'sample'
        df_type.rename(columns={'Homo': 'HomoSNP', 'Hete': 'HeteSNP'}, inplace=True)
        df_type['Total'] = df_type.apply(lambda x: x.sum(), axis=1)
        df_type.to_csv(self.work_dir + "/statis_hh", sep="\t")

        df_cds = df_cds.T[['CDS', 'non-CDS']]  # 维持列的顺序
        df_cds.index.name = 'sample'
        df_cds.rename(columns={'CDS': 'cds', 'non-CDS': 'non_cds'}, inplace=True)
        df_cds['total'] = df_cds.apply(lambda x: x.sum(), axis=1)
        df_cds.to_csv(self.work_dir + "/statis_cds", sep="\t")

        df_transition.index.name = 'type'
        df_transition.to_csv(self.work_dir + "/transition_tranversion", sep="\t")
        df_tt = pd.read_table(self.work_dir + "/transition_tranversion", sep="\t", index_col=0)
        df_tt.index.name = 'type'
        tt_per = df_tt / df_tt.sum()
        tt_per.columns = [x + "_per" for x in tt_per.columns]
        tt_new_per = pd.concat([df_tt, tt_per], axis=1)
        tt_new_per = tt_new_per.round(6)
        tt_new_per.to_csv(self.work_dir + "/tt_new_per", sep="\t")

        df_anno.index.name = 'sample'
        df_anno.to_csv(self.work_dir + "/snp_anno_stat", sep="\t")

    def run(self):
        super(SnpfinalNew2Tool, self).run()
        self.splitvcf()
        self.bedpre()
        self.snpindel_detail(self.work_dir + "/indel", self.work_dir + "/indel_new")
        if os.path.getsize(self.work_dir + "/indel_new") > 0:
            self.rewrite(self.work_dir + "/indel_new", self.work_dir + "/indel_rewrite")
            new_indel_rewrite = pd.read_table(self.work_dir + "/indel_rewrite", sep="\t", header=0)
            new_indel_rewrite.insert(4, "REF_ALT", new_indel_rewrite["REF"] + new_indel_rewrite["ALT"])
            new_indel_rewrite.to_csv(self.work_dir + "/new_indel_rewrite", sep="\t")

        self.snpindel_detail(self.work_dir + "/snp", self.work_dir + "/snp_new")
        if os.path.getsize(self.work_dir + "/snp_new") > 0:
            self.rewrite(self.work_dir + "/snp_new", self.work_dir + "/snp_rewrite")
            self.snp_analyse(self.work_dir + "/snp_rewrite")
            new_snp_rewrite = pd.read_table(self.work_dir + "/snp_rewrite", sep="\t", header=0)
            new_snp_rewrite.insert(4, "REF_ALT", new_snp_rewrite["REF"] + '/' + new_snp_rewrite["ALT"])
            new_snp_rewrite.to_csv(self.work_dir + "/new_snp_rewrite", sep="\t", index=False)
        self.set_output()
        self.end()

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        if self.option("method") == "samtools":
            filter_vcf = self.option('call_vcf').prop["path"]
            os.link(filter_vcf, os.path.join(self.output_dir, 'filter_vcf'))
        else:
            filter_snp_vcf = self.option('filted_snp_vcf').prop["path"]
            filter_indel_vcf = self.option('filted_indel_vcf').prop["path"]
            os.link(filter_snp_vcf, os.path.join(self.output_dir, 'filter_snp_vcf'))
            os.link(filter_indel_vcf, os.path.join(self.output_dir, 'filter_indel_vcf'))
        if os.path.getsize(self.work_dir + "/snp_new") > 0:
            snp_rewrite = self.work_dir + '/new_snp_rewrite'
            new_snp_rewrite = self.work_dir + '/new_snp_rewrite'
            depth = self.work_dir + '/depth'
            depth_new_per = self.work_dir + '/depth_new_per'
            tt_new_per = self.work_dir + '/tt_new_per'
            transition_tranversion = self.work_dir + '/transition_tranversion'
            statis_hh = self.work_dir + '/statis_hh'
            statis_cds = self.work_dir + "/statis_cds"
            statis_anno = self.work_dir + "/snp_anno_stat"
            os.link(snp_rewrite, os.path.join(self.output_dir, 'snp_detail'))
            os.link(depth, os.path.join(self.output_dir, 'snp_depth_statistics'))
            os.link(statis_hh, os.path.join(self.output_dir, 'snp_homo_hete_statistics'))
            os.link(transition_tranversion, os.path.join(self.output_dir, 'snp_transition_tranversion_statistics'))
            os.link(statis_cds, os.path.join(self.output_dir, 'snp_cds_statistics'))
            os.link(statis_anno, os.path.join(self.output_dir, 'snp_anno_statistics'))
            self.logger.info("成功设置SNP结果目录")
        else:
            self.logger.info("此次分析没有call得到snp信息，设置snp_detail结果目录失败")

        if os.path.getsize(self.work_dir + "/indel_new") > 0:
            indle_rewrite = self.work_dir + '/new_indel_rewrite'
            os.link(indle_rewrite, os.path.join(self.output_dir, 'indel_detail'))
        else:
            self.logger.info("此次分析没有call得到indel信息，设置indel_detail结果目录失败")


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    test_dir = '/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/tools/denovo_rna_v2/test_files'

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "SnpfinalNew2" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v3.snpfinal_new2",
            "instant": False,
            "options": dict(
                filted_snp_vcf="/mnt/ilustre/users/isanger/workspace/20200913/Denovorna_majorbio_283598/Sentieon/VcfFilterGatk/output/pop.snp.filter.recode.vcf",
                filted_indel_vcf="/mnt/ilustre/users/isanger/workspace/20200913/Denovorna_majorbio_283598/Sentieon/VcfFilterGatk/output/pop.indel.filter.recode.vcf",
                cds_bed="/mnt/ilustre/users/isanger/workspace/20200913/Denovorna_majorbio_283598/AnnotOrfpfam/output/all_predicted.bed",
                allt2g="/mnt/ilustre/users/isanger/workspace/20200913/Denovorna_majorbio_283598/AnnotOrfpfam/output/all_tran2gen.txt",
                anno="/mnt/ilustre/users/isanger/workspace/20200913/Denovorna_majorbio_283598/AnnotClassBeta/output/all_annot.xls",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
