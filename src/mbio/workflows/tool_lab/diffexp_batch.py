# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import pandas as pd
import glob
import os
import json
import time
import re
from biocluster.file import getsize, exists
from biocluster.file import download
import pandas as pd


class DiffexpBatchWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffexpBatchWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_type', type='string', default="tpm"),
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name="diff_main_id", type="string"),
            dict(name="count", type="infile", format="ref_rna_v2.common"),
            dict(name="group", type="string"),
            dict(name="cmp", type="string"),
            dict(name="type", type="string"),
            dict(name="task_id", type="string"),
            dict(name="exp_level", type="string"),
            dict(name="is_batch", type="bool"),  #是否做批次效应处理
            dict(name="has_batch", type="bool"),  #是否有批次处理表
            dict(name="batch_matrix", type="infile", format="ref_rna_v2.common"),  #上传批次效应表
            # pvalue 统计判断阈值
            dict(name="pvalue", type="float", default=0.05),
            dict(name="pvalue_padjust", type="string", default="padjust"),
            dict(name="fc", type="float", default=2),
            dict(name="padjust_way", type='string', default="BH"),
            # method: DESeq2, edgeR, DEGseq
            dict(name="method", type="string", default="DESeq2"),
            dict(name="tpm_filter_threshold",type="float",default="0"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("tool_lab.diffexp_batch")
        self.all_exp = self.api.api("ref_rna_v2.all_exp")

    def run(self):
        self.tool.on("end", self.set_db)
        self.run_tool()
        super(DiffexpBatchWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        # add result info
        self.all_exp.add_diffexp(self.tool.output_dir,
                            main_id=self.option('diff_main_id'),
                            diff_method=self.option('method'),
                            create_geneset=False,
                            pvalue_padjust=self.option('pvalue_padjust'),
                            )
        self.paste_annotation()
        self.set_output()
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "差异分析结果目录", 0, "211063"],
        ])
        result_dir.add_regexp_rules([
            [r'.*_vs_.*\.xls', 'xls', '差异表达基因详情表',0,"211518"],
            [r'.*summary.*\.xls', 'xls', '差异表达基因统计表',0,"211519"],
            [r'.*total_diff_stat.*\.xls', 'xls', '差异表达基因详情总表',0,"211520"],
        ])
        super(DiffexpBatchWorkflow, self).end()

    def set_output(self):
        diff_total=pd.read_table(os.path.join(self.tool.output_dir,"total_diff_stat.{}.xls".format(self.option("method").lower())),index_col="seq_id")
        diff_files_stat = glob.glob(os.path.join(self.tool.output_dir, '*.{}.xls'.format(self.option("method").lower())))
        df_anno = pd.read_table(diff_files_stat[0], index_col="seq_id")
        df1_genebase = pd.DataFrame(df_anno, columns=["gene_name", "description", "length"])
        df1_geneanno=pd.DataFrame(df_anno,columns=["nr", "go", "KO_id", "KO_name", "paths", "cog", "cog_description", "pfam", "swissprot", "entrez"])
        dff1=pd.concat([df1_genebase,diff_total], axis=1)
        dff2=pd.concat([dff1,df1_geneanno],axis=1)
        dff2.index.set_names("seq_id", inplace=True)
        dff2.to_csv(self.tool.output_dir + "/total_diff_stat.{}.xls".format(self.option("method").lower()), sep="\t")
        # with open(self.tool.output_dir + "/total_diff_stat.{}.xls".format(self.option("method").lower()), 'r+') as f:
        #     content = f.read()
        #     f.seek(0, 0)
        #     f.write("seq_id" + content)


        # diff_files_stat = glob.glob(os.path.join(self.tool.output_dir, '*.{}.xls'.format(self.option("method").lower())))
        # df1 = pd.DataFrame()
        # #f_head = "seq_id\t \t \t"
        # df_anno=pd.read_table(diff_files_stat[0],index_col="seq_id")
        # df1_genebase=pd.DataFrame(df_anno,columns=["gene_name", "description", "length"])
        # df1=pd.concat([df1, df1_genebase], axis=1)
        # df1_geneanno=pd.DataFrame(df_anno,columns=["nr", "go", "KO_id", "KO_name", "paths", "cog", "cog_description", "pfam", "swissprot", "entrez"])
        # for diff_file_stat in diff_files_stat:
        #     ctrl, test = os.path.basename(diff_file_stat).split('.{}.xls'.format(self.option("method").lower()))[0].split('_vs_')
        #     fname = os.path.basename(diff_file_stat).split(".")[0]
        #     #f_head = f_head + fname+"\t*\t*\t*\t*\t*\t*\t"
        #     df_t = pd.read_table(diff_file_stat, index_col="seq_id")
        #     df_core = pd.DataFrame(df_t, columns=["fc({}/{})".format(test,ctrl),"log2fc({}/{})".format(test,ctrl), "pvalue", "padjust", "significant", "regulate"])
        #     df1 = pd.concat([df1, df_core], axis=1)
        #     #df1 = df1.join(df_core, how="left")
        # df1 = pd.concat([df1,df1_geneanno], axis=1)
        # df1.to_csv(self.tool.output_dir + "/total_diff_stat.xls", sep="\t")
        # with open(self.tool.output_dir + "/total_diff_stat.xls", 'r+') as f:
        #     content = f.read()
        #     f.seek(0, 0)
        #     f.write(f_head + "\n" + content)

    def run_tool(self):
        # filter count file
        set_total_fordiff=set()
        finalset_total_fordiff=set()
#        if self.option("type") == 'all':
#           count_file_t = self.option('count').prop['path']
#        else:
#            count_file_t = self.work_dir + '/total_known_seqs_count.matrix'
#            with open(count_file_t, 'w') as fw, open(self.option('count').prop['path']) as fr:
#                for line in fr:
#                    if not line.startswith('MSTRG') and not line.startswith('TCONS') and not line.startswith('XLOC'):
#                        fw.write(line)
#                        set_total_fordiff.add(line.strip().split("\t")[0])
#        count_file_final=self.work_dir + '/final_used_seqs_count.matrix'
#        exp_file=self.option('exp_matrix')
#        exp_file_final=self.work_dir + '/final_exp_matrix'
#        exp_df = pd.read_table(exp_file, index_col=0, header=0)
#        exp_final = exp_df[exp_df > 0].dropna(axis=0, how='any')
#        set_nosig_all = set(exp_df.index)-set(exp_final.index)
#        with open(count_file_final, 'w') as fw, open(count_file_t,"r") as fr:
#            head=fr.readline()
#            fw.write(head)
#            for line in fr.readlines():
#               if line.strip().split("\t")[0] not in set_nosig_all:
#                    fw.write(line)
#                    finalset_total_fordiff.add(line.strip().split("\t")[0])
#        f_list=list(finalset_total_fordiff)
#        set_nosig=set_total_fordiff-finalset_total_fordiff
#        exp_f=exp_final.loc[f_list]
#        exp_f.to_csv(exp_file_final, sep='\t')
#        nosig_list=self.work_dir+'/nosig_seq'
#        with open(nosig_list,"w") as nosig:
#            nosig.write("seq_id"+"\n")
#            nosig.write('\n'.join(set_nosig))

        if self.option("type") == 'all':
            count_file = self.option('count').prop['path']
        else:
            count_file = self.work_dir + '/known_seqs_count.matrix'
            with open(count_file, 'w') as fw, open(self.option('count').prop['path']) as fr:
                for line in fr:
                    if not line.startswith('MSTRG') and not line.startswith('TCONS') and not line.startswith('XLOC'):
                        fw.write(line)
        if self.option('is_batch') == False:

            options = dict(
                exp=self.option('exp_matrix'),
                exp_type=self.option('exp_type'),
                count=count_file,
                group=self.option('group'),
                method=self.option('method'),
                cmp=self.option('cmp'),
                pvalue=self.option('pvalue'),
                pvalue_padjust=self.option('pvalue_padjust'),
                padjust_way=self.option('padjust_way'),
                fc=self.option('fc'),
                tpm_filter_threshold=self.option("tpm_filter_threshold"),
                is_batch=self.option('is_batch'),
            )
        else:
            if self.option('has_batch') == True:
                options = dict(
                    exp=self.option('exp_matrix'),
                    exp_type=self.option('exp_type'),
                    count=count_file,
                    group=self.option('group'),
                    method=self.option('method'),
                    cmp=self.option('cmp'),
                    pvalue=self.option('pvalue'),
                    pvalue_padjust=self.option('pvalue_padjust'),
                    padjust_way=self.option('padjust_way'),
                    fc=self.option('fc'),
                    tpm_filter_threshold=self.option("tpm_filter_threshold"),
                    is_batch=self.option('is_batch'),
                    has_batch=self.option('has_batch'),
                    batch_matrix=self.option('batch_matrix')
                )
            else:
                options = dict(
                    exp=self.option('exp_matrix'),
                    exp_type=self.option('exp_type'),
                    count=count_file,
                    group=self.option('group'),
                    method=self.option('method'),
                    cmp=self.option('cmp'),
                    pvalue=self.option('pvalue'),
                    pvalue_padjust=self.option('pvalue_padjust'),
                    padjust_way=self.option('padjust_way'),
                    fc=self.option('fc'),
                    tpm_filter_threshold=self.option("tpm_filter_threshold"),
                    is_batch=self.option('is_batch'),
                    has_batch=self.option('has_batch'),
                )
        self.tool.set_options(options)
        self.tool.run()

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
            self.set_error('file can not find %s', variables=(path), code = '13700502')
        return to_path

    def paste_annotation(self):
        conn = self.all_exp.db['sg_annotation_stat']
        task_id = self.option('task_id')
        try:
            find_result = conn.find_one({"task_id": task_id, "type": "latest", "status": "end"})
        except:
            self.set_error('annotation result file cannot found', code = "13700501")
        else:
            if find_result is None:
                find_result = conn.find_one({"task_id": task_id, "type": "origin"})

        # 2019.01.16 bug 当不选择组装时，不产生allannot_class，仅有refannot_class，并且需要通过exists函数判断s3上的文件
        if exists(os.path.join(find_result['result_dir'], 'allannot_class/all_annot.xls')):
            remote_annot = os.path.join(find_result['result_dir'], 'allannot_class/all_annot.xls')
            annot = self.download_s3_file(remote_annot, "annot.xls")
        elif exists(os.path.join(find_result['result_dir'], 'refannot_class/all_annot.xls')):
            remote_annot = os.path.join(find_result['result_dir'], 'refannot_class/all_annot.xls')
            annot = self.download_s3_file(remote_annot, "annot.xls")
        else:
            self.set_error('annotation result file cannot found', code="13700501")
        self.logger.info("annot dir is {}".format(annot))

        exp_level = self.option('exp_level')
        all_annot = pd.read_table(annot, header=0, index_col=0, )
        if exp_level[0].upper() == 'G':
            annot_pd = all_annot[all_annot["is_gene"] == "yes"].drop(columns=['transcript_id', 'is_gene'])
        else:
            annot_pd = all_annot.reset_index().drop(columns=['is_gene']).set_index('transcript_id')

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
            b = pd.DataFrame(diff_pd,columns = ["fc","log2fc", "pvalue", "padjust", "significant", "regulate"])
            c=diff_pd.drop(["fc","log2fc", "pvalue", "padjust", "significant", "regulate"], axis=1)
            result = pd.concat([b, c], axis=1)
            result.rename(columns={"log2fc": "log2fc({}/{})".format(test, ctrl),"fc": "fc({}/{})".format(test, ctrl)}, inplace=True)
            order = ["nr", "go", "KO_id", "KO_name", "paths", "cog", "cog_description", "pfam", "swissprot", "entrez"]
            order2 = ["gene_name", "description", "length"]
            pd1 = pd.DataFrame(annot_pd, columns=order)
            pd2 = pd.DataFrame(annot_pd, columns=order2)
            pd3 = pd.concat([pd2, result], axis=1,join_axes=[result.index])
            pd4 = pd.concat([pd3, pd1], axis=1,join_axes=[pd3.index])
            pd4.to_csv(each, sep='\t', header=True, index=True)




