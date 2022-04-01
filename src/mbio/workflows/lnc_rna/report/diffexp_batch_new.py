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
from mbio.packages.ref_rna_v2.functions import tryforgood
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import shutil

class DiffexpBatchNewWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffexpBatchNewWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_type', type='string', default="tpm"),
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name="diff_main_id", type="string"),
            dict(name="count", type="infile", format="lnc_rna.common"),
            dict(name="group", type="string"),
            dict(name="cmp", type="string"),
            dict(name="type", type="string"),
            dict(name="task_id", type="string"),
            dict(name="exp_level", type="string"),
            dict(name="seq_type_path",type="infile",format="lnc_rna.common"),#20190423 fwy新增，rna_seq_type文件路径
            # pvalue 统计判断阈值
            dict(name="pvalue", type="float", default=0.05),
            dict(name="pvalue_padjust", type="string", default="padjust"),
            dict(name="fc", type="float", default=2),
            dict(name="padjust_way", type='string', default="BH"),
            # method: DESeq2, edgeR, DEGseq
            dict(name="method", type="string", default="DESeq2"),
            dict(name="rna_type", type="string"),
            dict(name="filter_method", type="string", default=None),
            dict(name="tpm_filter_threshold", type="float", default="0"),
            # batch by zjx 20200917
            dict(name="is_batch", type="bool", default=False),
            dict(name="has_batch", type="bool"),
            dict(name="batch_matrix", type="infile", format="lnc_rna.common"),
            # noiseq
            dict(name='prob', type='float', default=0.8),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("lnc_rna.diffexp_batch")
        # self.all_exp = self.api.api("lnc_rna.all_exp")
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/02 Diff_Express')
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
        super(DiffexpBatchNewWorkflow, self).send_log(data)

    def run(self):
        self.tool.on("end", self.run_uniform)
        self.get_run_log()
        self.run_tool()
        super(DiffexpBatchNewWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("lnc_rna", table="sg_diff", main_id=self.option('diff_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        # add result info
        seq_type_path=self.option("seq_type_path").prop['path']
        self.stop_timeout_check()
        self.all_exp = self.api.api("lnc_rna.all_exp")
        # add result info
        self.all_exp.add_diffexp_all(self.uniform.output_dir, self.tool.output_dir,
                                     group_dict=self.option('group_dict'),
                                     main_id=self.option('diff_main_id'),
                                     diff_method=self.option('method'),
                                     create_geneset=False,
                                     pvalue_padjust=self.option('pvalue_padjust'),
                                     seq_type=seq_type_path,
                                     rna_type=self.option("rna_type")
                                     )
        self.paste_annotation()
        self.set_output()
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.tool.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))

        if os.path.exists(self.output_dir + '/upload'):
            os.remove(self.output_dir + '/upload')
        shutil.copytree(self.tool.output_dir, self.output_dir + '/upload')
        rm_files = glob.glob(self.output_dir + '/upload/*.DE.list')
        rm_files += glob.glob(self.output_dir + '/upload/*.normalize.xls')
        rm_files += glob.glob(self.output_dir + '/upload/*.sizeFactor.xls')
        rm_files += glob.glob(self.output_dir + '/upload/*.normFactor.xls')
        for each in rm_files:
            os.remove(each)

        result_dir = self.add_upload_dir(self.output_dir + '/upload')
        self.inter_dirs = [
            ["02 Diff_Express", "", "表达量差异分析结果目录", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "基因表达量差异分析文件", 0, "211063"],
            ['run_parameter.txt', 'txt', '运行参数日志', 0]
        ])
        result_dir.add_regexp_rules([
            [r'.*_diff_summary.xls', 'xls', '表达量差异统计表', 0],
            [r'.*_vs_.*\..*\.xls', 'xls', '表达量差异详情表', 0],
            [r'.*_vs_.*\..*\.annot.xls', 'xls', '表达量差异注释表', 0],
        ])

        super(DiffexpBatchNewWorkflow, self).end()

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

    def run_tool(self):
        count_file=os.path.join(self.work_dir,"seqs_count.matrix")
        with open(self.option('count').prop['path'], "r") as cf:
            with open(count_file,"w") as ff:
                        firstline = cf.readline()
                        ff.write(firstline)
                        exp_records = cf.readlines()
                        if self.option('type') == "ref":
                            for record in exp_records:
                                line_info = record.strip().split()
                                if line_info[-1] == "false":
                                    if self.option('rna_type') == "lncRNA":
                                        numinfo1 = len(line_info)
                                        if line_info[numinfo1-2] == "lncRNA":
                                            ff.write(record)
                                    elif self.option('rna_type') == "mRNA":
                                        numinfo1 = len(line_info)
                                        if line_info[numinfo1-2] == "mRNA":
                                            ff.write(record)
                                    else:
                                        ff.write(record)
                        else:
                            for record in exp_records:
                                line_info = record.strip().split()
                                if self.option('rna_type') == "lncRNA":
                                    numinfo1 = len(line_info)
                                    if line_info[numinfo1-2] == "lncRNA":
                                        ff.write(record)
                                elif self.option('rna_type') == "mRNA":
                                    numinfo1 = len(line_info)
                                    if line_info[numinfo1-2] == "mRNA":
                                        ff.write(record)
                                else:
                                    ff.write(record)
        ff.close()
        count_matrix = pd.read_table(count_file)
        count_matrix = count_matrix.set_index('seq_id')
        # exp_matrix.drop(exp_matrix.columns[len(exp_matrix.columns) - 1], axis=1, inplace=True)
        count_matrix.drop(['rna_type'], axis=1, inplace=True)
        count_matrix.drop(['is_new'], axis=1, inplace=True)
        count_output = os.path.join(self.work_dir,"count_matrix")
        count_matrix.to_csv(count_output, sep='\t', header=True, index=True)
        print('success to export count matrix')
        opts = dict(
            exp=self.option('exp_matrix'),
            exp_type=self.option('exp_type'),
            count=count_output,
            group=self.option('group'),
            method=self.option('method'),
            cmp=self.option('cmp'),
            tpm_filter_threshold=self.option("tpm_filter_threshold"),
            is_batch=self.option('is_batch'),
            filter_method=self.option("filter_method"),
            no_filter=True,
            analysis='split',
            fc=self.option('fc'),
        )
        if self.option('is_batch') == False:
            if self.option("method").lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:
                opts.update(
                    pvalue=self.option('pvalue'),
                    pvalue_padjust=self.option('pvalue_padjust'),
                    padjust_way=self.option('padjust_way'),
                )
            else:
                opts.update(
                    prob=self.option('prob')
                )
        else:
            if self.option('has_batch') == True:
                opts = dict(
                    exp=self.option('exp_matrix'),
                    exp_type=self.option('exp_type'),
                    count=count_output,
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
                    batch_matrix=self.option('batch_matrix'),
                    filter_method=self.option("filter_method"),
                    no_filter=True,
                    analysis='split'
                )
            else:
                opts = dict(
                    exp=self.option('exp_matrix'),
                    exp_type=self.option('exp_type'),
                    count=count_output,
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
                    filter_method=self.option("filter_method"),
                    no_filter=True,
                    analysis='split'
                )
        self.tool.set_options(opts)
        self.tool.run()

    def run_uniform(self):
        self.uniform = self.add_tool('denovo_rna_v3.batch.uniform')
        opts = {
            'method': self.option('method'),
            'input_dir': self.tool.output_dir,
        }
        if self.option('method').lower() in ['deseq2', 'degseq', 'edger', 'limma']:
            opts.update({'pvalue_padjust': self.option('pvalue_padjust')})
        self.uniform.set_options(opts)
        self.uniform.on('end', self.set_db)
        self.uniform.run()

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
        annot = os.path.join(find_result['result_dir'], 'allannot_class/all_annot.xls')
        annot = self.download_s3_file(annot, "annot.xls")
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
            # diff_pd = pd.read_table(each, header=0, sep='\t', index_col=0)
            # diff_pd = diff_pd.join(annot_pd, how='left')
            # out_diff = each[:-3] + 'annot.xls'
            # diff_pd.to_csv(out_diff, sep='\t', header=True, index=True)
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
            order = ["nr", "go", "KO_id", "KO_name", "paths", "cog", "cog_description", "pfam", "swissprot"]
            order2 = ["length"]
            pd1 = pd.DataFrame(annot_pd, columns=order)
            pd2 = pd.DataFrame(annot_pd, columns=order2)
            pd3 = pd.concat([pd2, result], axis=1, join_axes=[result.index])
            pd4 = pd.concat([pd3, pd1], axis=1, join_axes=[pd3.index])
            pd4.to_csv(each, sep='\t', header=True, index=True)
