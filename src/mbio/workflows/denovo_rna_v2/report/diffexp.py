# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
import pandas as pd
import os
import glob
from biocluster.file import getsize, exists
from biocluster.file import download
from mbio.packages.ref_rna_v2.functions import tryforgood


class DiffexpWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffexpWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name="diff_main_id", type="string"),
            dict(name="count", type="string",),
            dict(name="group", type="string"),
            dict(name="cmp", type="string"),
            # pvalue 统计判断阈值
            dict(name="pvalue", type="float", default=0.05),
            dict(name="pvalue_padjust", type="string", default="padjust"),
            dict(name="fc", type="float", default=2),
            dict(name="filter_method", type="string", default=None),  # filter_method 20190801修改
            dict(name="tpm_filter_threshold", type="float", default="0"), #20190708 添加 by fwy
            dict(name="padjust_way", type='string', default="BH"),
            # method: DESeq2, edgeR, DEGseq
            dict(name="method", type="string", default="DESeq2"),
            # DESeq2 DE test method, Wald|LRT
            dict(name="deseq2_method", type="string", default="Wald"),
            # edger_method DE test method,exactTest|glmLRT|glmQLFTest
            dict(name="edger_method", type="string", default="glmQLFTest"),
            # degseq_method DE test method, LRT|CTR|FET|MARS|MATR|FC
            dict(name="degseq_method", type="string", default="MARS"),
            dict(name="task_id", type="string", default=""),
            dict(name="exp_level", type="string", default=""),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("denovo_rna_v2.diffexp")

    def run(self):
        self.tool.on("end", self.set_db)
        self.run_tool()
        super(DiffexpWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.all_exp = self.api.api("denovo_rna_v2.all_exp")
        # add result info
        self.all_exp.add_diffexp(self.tool.output_dir,
                            main_id=self.option('diff_main_id'),
                            diff_method=self.option('method'),
                            create_geneset=False,
                            pvalue_padjust=self.option('pvalue_padjust')
                            )
        self.paste_annotation()
        self.set_output()
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "差异分析结果目录", 0, '201158'],
        ])
        result_dir.add_regexp_rules([
            [r'.*_vs_.*\.xls', 'xls', '差异表达基因详情表',0,"201366"],
            [r'.*summary.*\.xls', 'xls', '差异表达基因统计表',0,"201367"],
            [r'.*total_diff_stat.*\.xls', 'xls', '差异表达基因详情总表',0,"201368"],
        ])
        super(DiffexpWorkflow, self).end()

    def set_output(self):
        diff_total=pd.read_table(os.path.join(self.tool.output_dir,"total_diff_stat.{}.xls".format(self.option("method").lower())),index_col="seq_id")
        diff_files_stat = glob.glob(os.path.join(self.tool.output_dir, '*.{}.xls'.format(self.option("method").lower())))
        df_anno = pd.read_table(diff_files_stat[0], index_col="seq_id")
        df1_genebase = pd.DataFrame(df_anno, columns=["length"])
        df1_geneanno=pd.DataFrame(df_anno,columns=["nr", "go", "KO_id", "KO_name", "paths", "cog", "cog_description", "pfam", "swissprot"])
        dff1 = pd.concat([df1_genebase, diff_total], axis=1)
        dff2 = pd.concat([dff1, df1_geneanno], axis=1)
        dff2.index.set_names("seq_id", inplace=True)
        dff2.to_csv(self.tool.output_dir + "/total_diff_stat.{}.xls".format(self.option("method").lower()), sep="\t")
        # if len(diff_files_stat) >1:
        #     with open(self.tool.output_dir + "/total_diff_stat.{}.xls".format(self.option("method").lower()),'r+') as f:
        #         content = f.read()
        #         f.seek(0, 0)
        #         f.write("seq_id" + content)

        # df_anno=pd.read_table(diff_files_stats[0],index_col="seq_id")
        # df1_genebase=pd.DataFrame(df_anno,columns=["length"])
        # df1=pd.concat([df1, df1_genebase], axis=1)
        # df1_geneanno=pd.DataFrame(df_anno,columns=["nr", "go", "KO_id", "KO_name", "paths", "cog", "cog_description", "pfam", "swissprot", "entrez"])
        # for diff_file_stat in diff_files_stats:
        #     ctrl, test = os.path.basename(diff_file_stat).split('.{}.xls'.format(self.option("method").lower()))[0].split('_vs_')
        #     fname = os.path.basename(diff_file_stat).split(".")[0]
        #     f_head = f_head + fname+"\t*\t*\t*\t*\t*\t*\t"
        #     df_t = pd.read_table(diff_file_stat, index_col="seq_id")
        #     df_core = pd.DataFrame(df_t, columns=["fc({}/{})".format(test,ctrl),"log2fc({}/{})".format(test,ctrl), "pvalue", "padjust", "significant", "regulate"])
        #     df1 = pd.concat([df1, df_core], axis=1)
        #     #df1 = df1.join(df_core, how="left")
        # df1 = pd.concat([df1,df1_geneanno], axis=1)
        # df1.to_csv(self.tool.output_dir + "/total_diff_stat.xls", sep="\t")
        # if len(diff_files_stats) >0:
        # with open(self.tool.output_dir + "/total_diff_stat.xls", 'r+') as f:
        #     content = f.read()
        #     f.seek(0, 0)
        #     f.write(f_head + "\n" + content)

    def run_tool(self):
        options = dict(
            exp=self.option('exp_matrix'),
            count=self.option('count'),
            group=self.option('group'),
            method=self.option('method'),
            cmp=self.option('cmp'),
            pvalue=self.option('pvalue'),
            pvalue_padjust=self.option('pvalue_padjust'),
            padjust_way=self.option('padjust_way'),
            fc=self.option('fc'),
        )
        if self.option("filter_method") is not None:
            options["filter_method"]=self.option("filter_method")
            options["tpm_filter_threshold"]=self.option("tpm_filter_threshold")
        if self.option("method").lower() == "edger":
            options["edger_method"]=self.option("edger_method")
        elif self.option("method").lower() == "degseq":
            options["degseq_method"]=self.option("degseq_method")
        elif self.option("method").lower() == "deseq2":
            options["deseq2_method"]=self.option("deseq2_method")
        self.tool.set_options(options)
        self.tool.run()

    @tryforgood
    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        elif os.path.exists(path):
            to_path = path
        elif exists(path):
            download(path, to_path)
        else:
            self.set_error('file can not find %s', variables=(path), code = '12003001')
        return to_path


    def paste_annotation(self):
        conn = self.all_exp.db['sg_annotation_stat']
        task_id = self.option('task_id')
        # try:
        #     find_result = conn.find_one({"task_id": task_id, "type": "latest", "status": "end"})
        # except:
        #     self.set_error('annotation result file cannot found', code = "13700501")
        # else:
        #     if find_result is None:
        find_result = conn.find_one({"task_id": task_id, "type": "origin"})

        exp_level = self.option('exp_level')
        if exp_level:
            self.logger.info("we find exp_level")
        else:
            self.logger.info("exp_level info miss")
        result_dir=find_result['result_dir']
        if exists(os.path.join(result_dir, 'all_annot.xls')):
            # remote_annot = os.path.join(result_dir, 'anno_stat/gene_anno_detail.xls')
            try:
                remote_annot = os.path.join(result_dir, 'all_annot.xls')
                annot = self.download_s3_file(remote_annot, "all_annot.xls")
            except:
                annot=os.path.join(self.work_dir,"all_annot.xls")
                os.link(os.path.join(result_dir, 'all_annot.xls'),annot)

            all_annot = pd.read_table(annot, header=0)
            if exp_level[0].upper() == 'G':
                annot_pd = all_annot[all_annot["is_gene"] == "yes"].drop(columns=['transcript', 'is_gene']).set_index(
                    'gene_id')
            else:
                annot_pd = all_annot.reset_index().drop(columns=['is_gene']).set_index('transcript')
        else:
            if exp_level[0].upper() == 'G':
                if exists(os.path.join(result_dir.replace("Annotation","AnnoQuery"), 'unigene_anno_detail.xls')):
                    remote_annot = os.path.join(os.path.join(result_dir.replace("Annotation","AnnoQuery"), 'unigene_anno_detail.xls'))
                    annot = self.download_s3_file(remote_annot, "gene_anno_detail.xls")
                    annot_pd = pd.read_table(annot, header=0,index_col="gene_id")
                else:
                    self.set_error('annotation result file cannot found', code="12003002")
            else:
                if exists(os.path.join(result_dir.replace("Annotation","AnnoQuery"), 'transcript_anno_detail.xls')):
                    remote_annot = os.path.join(os.path.join(result_dir.replace("Annotation","AnnoQuery"), 'transcript_anno_detail.xls'))
                    annot = self.download_s3_file(remote_annot, "trans_anno_detail.xls")
                    all_annot = pd.read_table(annot, header=0)
                    annot_pd = all_annot.reset_index().drop(columns=['gene_id']).set_index('transcript')
                else:
                    self.set_error('annotation result file cannot found', code="12003003")
        self.logger.info("annot dir is {}".format(annot))
        # if exp_level[0].upper() == 'G':
        #   if exists(os.path.join(find_result['result_dir'], 'anno_stat/gene_anno_detail.xls')):
        #     remote_annot = os.path.join(find_result['result_dir'], 'anno_stat/gene_anno_detail.xls')
        #     annot = self.download_s3_file(remote_annot, "annot.xls")
        #   else:
        #       self.set_error('annotation result file cannot found', code="13700501")
        # else:
        #   if exists(os.path.join(find_result['result_dir'], 'anno_stat/trans_anno_detail.xls')):
        #      remote_annot = os.path.join(find_result['result_dir'], 'anno_stat/trans_anno_detail.xls')
        #      annot = self.download_s3_file(remote_annot, "annot.xls")
        #   else:
        #       self.set_error('annotation result file cannot found', code="13700501")
        # self.logger.info("annot dir is {}".format(annot))


        for each in glob.glob(os.path.join(self.tool.output_dir, '*_vs_*.*.xls')):
            if each.endswith('.annot.xls'):
                continue
            if each.endswith('.normalize.xls'):
                continue
            if each.endswith('.sizeFactor.xls'):
                continue
            if each.endswith('.normFactor.xls'):
                continue
            diff_pd = pd.read_table(each, header=0, sep='\t', index_col="seq_id")
            ctrl, test = os.path.basename(each).split('.{}.xls'.format(self.option("method").lower()))[0].split('_vs_')
            b = pd.DataFrame(diff_pd, columns=["fc", "log2fc", "pvalue", "padjust", "significant", "regulate"])
            c = diff_pd.drop(["fc", "log2fc", "pvalue", "padjust", "significant", "regulate"], axis=1)
            result = pd.concat([b, c], axis=1)
            # result.rename(columns={"log2fc": "log2fc({}/{})".format(test, ctrl), "fc": "fc({}/{})".format(test, ctrl)},
            #               inplace=True)
            order = ["nr", "go", "KO_id", "KO_name", "paths", "cog", "cog_description", "pfam", "swissprot"]
            order2 = ["length"]
            pd1 = pd.DataFrame(annot_pd, columns=order)
            pd2 = pd.DataFrame(annot_pd, columns=order2)
            pd3 = pd.concat([pd2, result], axis=1, join_axes=[result.index])
            pd4 = pd.concat([pd3, pd1], axis=1, join_axes=[pd3.index])
            pd4.to_csv(each, sep='\t', header=True, index=True)


            #diff_pd = pd.read_table(each, header=0, sep='\t', index_col=0)
            #diff_pd = diff_pd.join(annot_pd, how='left')
            #out_diff = each[:-3] + 'annot.xls'
            #diff_pd.to_csv(out_diff, sep='\t', header=True, index=True)




