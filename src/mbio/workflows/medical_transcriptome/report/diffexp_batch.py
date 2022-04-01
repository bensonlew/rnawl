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
import re
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


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
            dict(name="filter_method", type="string" , default="no"),
            dict(name="count", type="infile", format="ref_rna_v2.common"),
            dict(name="group", type="string"),
            dict(name="cmp", type="string"),
            dict(name="type", type="string"),
            dict(name="task_id", type="string"),
            dict(name="exp_level", type="string"),
            dict(name="is_batch", type="bool", default=False),  #是否做批次效应处理
            dict(name="has_batch", type="bool"),  #是否有批次处理表
            dict(name="batch_matrix", type="infile", format="ref_rna_v2.common"),  #上传批次效应表
            # pvalue 统计判断阈值
            dict(name="pvalue", type="float", default=0.05),
            dict(name="pvalue_padjust", type="string", default="padjust"),
            dict(name="fc", type="float", default=2),
            dict(name="padjust_way", type='string', default="BH"),
            # method: DESeq2, edgeR, DEGseq, limma, NOIseq
            dict(name="method", type="string", default="DESeq2"),
            dict(name='prob', type='float', default=0.8),
            dict(name="tpm_filter_threshold",type="float",default="0"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("medical_transcriptome.batch.diffexp_batch")
        self.all_exp = self.api.api("medical_transcriptome.all_exp")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/01 Diff_Express')
        self.inter_dirs = []

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(DiffexpBatchWorkflow, self).send_log(data)

    def run(self):
        self.tool.on("end", self.set_db)
        self.get_run_log()
        self.run_tool()
        super(DiffexpBatchWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("medical_transcriptome", table="sg_diff", main_id=self.option('diff_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        # add result info
        if self.option("method").lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:
            self.all_exp.add_diffexp(self.tool.output_dir,
                                group_dict = self.option("group_dict"),
                                main_id=self.option('diff_main_id'),
                                diff_method=self.option('method'),
                                create_geneset=False,
                                pvalue_padjust=self.option('pvalue_padjust'),
                                )
        else:
            self.all_exp.add_diffexp_noiseq(self.tool.output_dir,
                main_id=self.option('diff_main_id'),
                diff_method=self.option('method'),
                create_geneset=False,
            )
        self.paste_annotation()
        self.set_output()
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.tool.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.tool.output_dir)
        self.inter_dirs = [
            ["01 Diff_Express", "", "差异基因数据挖掘结果目录", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "基因表达量差异分析文件", 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0]
        ])
        result_dir.add_regexp_rules([
            [r'.*summary\.bar\.pdf', 'pdf', '表达量统计柱状图',0],
            [r'.*summary\.bar2\.pdf', 'pdf', '表达量统计堆积图', 0],
            [r'.*_vs_.*\.scatter\.pdf', 'pdf', '表达量差异ma图', 0],
            [r'.*_vs_.*\.volcano\.pdf', 'pdf', '表达量差异火山图', 0],
            [r'.*_vs_.*\.xls', 'xls', '表达量差异结果表',0],
            [r'.*summary.*\.xls', 'xls', '表达量差异统计表',0],
            [r'.*total_diff_stat.*\.xls', 'xls', '表达量差异详情总表',0],
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

        opts = dict(
            exp=self.option('exp_matrix'),
            exp_type=self.option('exp_type'),
            count=count_file,
            group=self.option('group'),
            method=self.option('method'),
            cmp=self.option('cmp'),
            fc=self.option('fc'),
            tpm_filter_threshold=self.option("tpm_filter_threshold"),
            is_batch=self.option('is_batch'),
        )

        if self.option('is_batch') == False:
            if self.option("method").lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:
                opts.update(
                    pvalue=self.option('pvalue'),
                    pvalue_padjust=self.option('pvalue_padjust'),
                    padjust_way=self.option('padjust_way'),
                )
                # options = dict(
                #     exp=self.option('exp_matrix'),
                #     exp_type=self.option('exp_type'),
                #     count=count_file,
                #     group=self.option('group'),
                #     method=self.option('method'),
                #     cmp=self.option('cmp'),
                       #     pvalue=self.option('pvalue'),
                #     pvalue_padjust=self.option('pvalue_padjust'),
                #     padjust_way=self.option('padjust_way'),
                #     fc=self.option('fc'),
                #     tpm_filter_threshold=self.option("tpm_filter_threshold"),
                #     is_batch=self.option('is_batch'),
                # )
            else:
                opts.update(
                    prob=self.option('prob')
                )
        else:
            if self.option('has_batch') == True:
                opts = dict(
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
                opts = dict(
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
        self.tool.set_options(opts)
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
            find_result = conn.find_one({"task_id": task_id, "status": "end"})
        except:
            self.set_error('annotation result file cannot found', code = "13700501")
        else:
            if find_result is None:
                find_result = conn.find_one({"task_id": task_id})

        # 2019.01.16 bug 当不选择组装时，不产生allannot_class，仅有refannot_class，并且需要通过exists函数判断s3上的文件
        if exists(os.path.join(find_result['result_dir'], 'allannot_class/all_annot_tran.xls')):
            remote_annot = os.path.join(find_result['result_dir'], 'allannot_class/all_annot_tran.xls')
            annot = self.download_s3_file(remote_annot, "annot.xls")
        elif exists(os.path.join(find_result['result_dir'], 'refannot_class/all_annot_tran.xls')):
            remote_annot = os.path.join(find_result['result_dir'], 'refannot_class/all_annot_tran.xls')
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
            if each.endswith('.rawnormal.xls'):
                continue
            if each.endswith('.normalize.xls'):
                continue
            if each.endswith('.sizeFactor.xls'):
                continue
            if each.endswith('.normFactor.xls'):
                continue
            diff_pd = pd.read_table(each, header=0, sep='\t', index_col="seq_id")
            ctrl, test = os.path.basename(each).split('.{}.xls'.format(self.option("method").lower()))[0].split('_vs_')
            if self.option("method").lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:

                b = pd.DataFrame(diff_pd,columns = ["fc","log2fc", "pvalue", "padjust", "significant", "regulate"])
                c=diff_pd.drop(["fc","log2fc", "pvalue", "padjust", "significant", "regulate"], axis=1)
            else:
                try:
                    b = pd.DataFrame(diff_pd,columns = ['{}_mean'.format(ctrl), '{}_mean'.format(test),"fc","log2fc", "theta", "D", 'prob', "significant", "regulate"])
                    c = diff_pd.drop(['{}_mean'.format(ctrl), '{}_mean'.format(test),"fc","log2fc", "theta", "D", 'prob', "significant", "regulate"], axis=1)
                except:
                    b = pd.DataFrame(diff_pd,columns = ['{}_mean'.format(ctrl), '{}_mean'.format(test),"fc","log2fc", "D", 'prob', "significant", "regulate"])
                    c = diff_pd.drop(['{}_mean'.format(ctrl), '{}_mean'.format(test),"fc","log2fc", "D", 'prob', "significant", "regulate"], axis=1)
            result = pd.concat([b, c], axis=1)
            result.rename(columns={"log2fc": "log2fc({}/{})".format(test, ctrl),"fc": "fc({}/{})".format(test, ctrl)}, inplace=True)
            order = ["nr", "go", "KO_id", "KO_name", "paths", "cog", "cog_description", "pfam", "swissprot", "entrez"]
            order2 = ["gene_name", "description", "length"]
            pd1 = pd.DataFrame(annot_pd, columns=order)
            pd2 = pd.DataFrame(annot_pd, columns=order2)
            pd3 = pd.concat([pd2, result], axis=1,join_axes=[result.index])
            pd4 = pd.concat([pd3, pd1], axis=1,join_axes=[pd3.index])
            pd4.to_csv(each, sep='\t', header=True, index=True)




