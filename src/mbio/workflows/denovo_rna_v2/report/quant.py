# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
import os
from biocluster.file import getsize, exists
from biocluster.file import download
import time
import pandas as pd
from biocluster.config import Config
from bson.objectid import ObjectId
import types
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
import re
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import json



class QuantWorkflow(Workflow):
    """
    Snp结构功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(QuantWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name="transcriptome", type="infile", format="denovo_rna_v2.trinity_fasta"),
            dict(name="fastq", type="infile", format="denovo_rna_v2.fastq_list"),
            dict(name="method", type="string", default="rsem"),
            dict(name="libtype", type="string", default=None),
            dict(name="t2g", type="infile", format="denovo_rna_v2.common"),
            dict(name="read_len", type="int", default=149),
            dict(name="transcript_main_id", type="string"),
            dict(name="gene_main_id", type="string"),
            dict(name="task_type", type="string"),
            dict(name="raw_task_id", type="string"),
            dict(name="submit_location", type="string"),
            dict(name="exp_type", type="string"),
            dict(name="task_id", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        # self.tool = self.add_tool("denovo_rna_v2.quant")
        self.tool = self.add_module("denovo_rna_v2.quant")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/02 Express/01 Exp_Annalysis')
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
        super(QuantWorkflow, self).send_log(data)

    def run(self):
        self.get_run_log()
        # db = Config().get_mongo_client(mtype="denovo_rna_v2")[Config().get_mongo_dbname("denovo_rna_v2")]
        # col = db["sg_task"]
        # result = col.find_one({"task_id": self.option('raw_task_id')})
        # basepath = result["fastq"].split("/fq_list.txt")[0]
        # with open (self.option("fastq").prop["path"], "r") as f, open(self.work_dir + "/fq_list.txt", "w") as w:
        #     for line in f:
        #         items = line.strip().split("\t")
        #         if (r'workspace', line):
        #             if len(items) == 3:
        #                 w.write(items[0] + "\t" + basepath + "/" + os.path.basename(items[1]) + "\t" + basepath + "/" + os.path.basename(items[2]) + "\n")
        #             else:
        #                 w.write(items[0] + "\t" + basepath + "/" + os.path.basename(items[1]) + "\n")
        #         else:
        #             w.writeline(line)
        # self.option("fastq", self.work_dir + "/fq_list.txt")
        if self.option('exp_type').lower() == 'fpkm' and self.option('method').lower() == 'rsem':
            self.start_listener()
            self.fire("start")
            self.dump_fpkm(self.option('raw_task_id'),
                           self.option('transcript_main_id'),
                           self.option('gene_main_id'))
            self.logger.info('self.dump_fpkm has been done, sleep 20s then self.end')
            # time.sleep(20)
            self.end()
            # super(QuantWorkflow, self).end()
        elif self.option('exp_type').lower() == 'tpm' and self.option('method').lower() == 'rsem':
            self.start_listener()
            self.fire("start")
            self.dump_tpm(self.option('raw_task_id'),
                           self.option('transcript_main_id'),
                           self.option('gene_main_id'))
            self.logger.info('self.dump_fpkm has been done, sleep 20s then self.end')
            # time.sleep(20)
            self.end()
            # super(QuantWorkflow, self).end()
        else:
            self.tool.on("end", self.set_db)
            self.run_tool()
            super(QuantWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("denovo_rna_v2", table="sg_exp", main_id=self.option('gene_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        all_exp = self.api.api("denovo_rna_v2.all_exp")
        # add transcript info
        # exp_matrix = self.work_dir + '/transcript.tpm.matrix'
        exp_matrix = self.tool.output_dir + '/transcript.tpm.matrix'
        all_exp.add_exp(exp_matrix,
                        quant_method=self.option('method'),
                        exp_level='T',
                        main_id=self.option('transcript_main_id'),
                        exp_type='tpm',
                        add_distribution=False,
                        task_id=self.option("raw_task_id"),)
        # add gene info
        exp_matrix = self.tool.output_dir + '/gene.tpm.matrix'
        all_exp.add_exp(exp_matrix,
                        quant_method=self.option('method'),
                        exp_level='G',
                        main_id=self.option('gene_main_id'),
                        exp_type='tpm',
                        add_distribution=False,
                        task_id=self.option("raw_task_id"),)
        self.paste_annotation()
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        # for ff in ['gene.fpkm.matrix.xls','transcript.fpkm.matrix.xls','gene.tpm.matrix.annot.xls',\
        # 'transcript.tpm.matrix.annot.xls','gene.tpm.matrix.xls','transcript.tpm.matrix.xls']:
        #     if os.path.exists(os.path.join(self.output_dir, ff)):
        #         os.remove(os.path.join(self.output_dir, ff))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["02 Express", "", "表达量结果目录",0],
            ["02 Express/01 Exp_Annalysis", "", "表达量统计分析", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "表达量统计分析文件",0,"201437"],
            ['gene.fpkm.matrix.annot.xls', '', '基因fpkm表达定量注释结果表',0,"201439"],
            ['transcript.fpkm.matrix.annot.xls', '', '转录本fpkm表达定量注释结果表',0,"201441"],
            ['gene.tpm.matrix.annot.xls', '', '基因tpm表达定量注释结果表',0],
            ['transcript.tpm.matrix.annot.xls', '', '转录本tpm表达定量注释结果表',0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(QuantWorkflow, self).end()

    def run_tool(self):
        options = {
            "method": self.option('method'),
            "libtype": self.option('libtype'),
            "read_len": self.option('read_len'),
            "t2g": self.option('t2g').prop['path'],
            "transcriptome": self.option('transcriptome').prop['path'],
            "fastq": self.option('fastq').prop['path'],
        }
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
            self.set_error('file can not find', code="12001703")
        return to_path

    def dump_fpkm(self, task_id, transcript_main_id, gene_main_id):
        all_exp = self.api.api("denovo_rna_v2.all_exp")
        collection = all_exp.db['sg_exp']
        if transcript_main_id == "":
            results = collection.find({"task_id": task_id, "method": "RSEM", "exp_level": "G"})
        else:
            results = collection.find({"task_id": task_id, "method": "RSEM"})
        transcript_fpkm_file = None
        gene_fpkm_file = None
        transcript_fpkm_anno_file = None
        gene_fpkm_anno_file = None
        gene_count_file = None
        transcript_count_file = None
        for each in results:
            # 2019.04.10 bug each['count_file']也需要改名称，不然实际上只是把s3上的*.count.matrix.xls下载下来后，
            # 更改名称导表而已。正确的做法是将s3上的*.fpkm.matrix.xls下载下来，即指定正确的下载文件。
            if each['exp_level'] == 'T' and 'count_file' in each:
                local_t_base = os.path.basename(each['count_file']).replace('.count.matrix', '.fpkm.matrix')
                # transcript_fpkm_file = self.download_s3_file(each['count_file'], local_t_base)
                transcript_fpkm_file = self.download_s3_file(
                    each['count_file'].replace('.count.matrix', '.fpkm.matrix'), local_t_base
                )
                transcript_count_file=each['count_file']
                transcript_fpkm_anno_file = each['count_file'].replace('.count.matrix.xls', '.fpkm.matrix.annot.xls')
                if exists(transcript_fpkm_anno_file):
                    transcript_fpkm_anno_file = self.download_s3_file(transcript_fpkm_anno_file, os.path.basename(transcript_fpkm_anno_file))
                else:
                    transcript_fpkm_anno_file = None
            if each['exp_level'] == 'G' and 'count_file' in each:
                local_g_base = os.path.basename(each['count_file']).replace('.count.matrix', '.fpkm.matrix')
                # gene_fpkm_file =  self.download_s3_file(each['count_file'], local_g_base)
                gene_fpkm_file =  self.download_s3_file(
                    each['count_file'].replace('.count.matrix', '.fpkm.matrix'), local_g_base
                )
                gene_count_file=each['count_file']
                gene_fpkm_anno_file = each['count_file'].replace('.count.matrix.xls', '.fpkm.matrix.annot.xls')
                if exists(gene_fpkm_anno_file):
                    gene_fpkm_anno_file = self.download_s3_file(gene_fpkm_anno_file,os.path.basename(gene_fpkm_anno_file))
                else:
                    gene_fpkm_anno_file = None
            if transcript_fpkm_file and gene_fpkm_file:
                break
        if transcript_fpkm_file is not None:
            all_exp.add_exp(transcript_fpkm_file, quant_method="RSEM", exp_level='T',
                            main_id=transcript_main_id, exp_type='fpkm', add_distribution=False, task_id=self.option('raw_task_id'))
            # ftresults = collection.find_one({"task_id": task_id, "method": "RSEM", 'exp_level': 'T', 'exp_type': "FPKM"})
            # ftresults["count_file"] = transcript_count_file
            # collection.update({"task_id": task_id, "method": "RSEM", 'exp_level': 'T', 'exp_type': "FPKM"}, {
            #     '$set': {'count_file': transcript_count_file}}, upsert=True)
            if isinstance(transcript_main_id, types.StringTypes) or type(transcript_main_id) == unicode:
                transcript_main_id = ObjectId(transcript_main_id)
            collection.update({"task_id": task_id, 'main_id':transcript_main_id}, {
                '$set': {'count_file': transcript_count_file}}, upsert=True)
        elif transcript_main_id == "":
            self.logger.info('only Unigene level')
        else:
            self.set_error('transcript_fpkm/count_file is not found', code = "12001701")

        if gene_fpkm_file is not None:
            all_exp.add_exp(gene_fpkm_file, quant_method="RSEM", exp_level='G',
                            main_id=gene_main_id, exp_type='fpkm', add_distribution=False, task_id=self.option('raw_task_id'))
            # fgresults = collection.find_one({"task_id": task_id, "method": "RSEM", 'exp_level': 'G', 'exp_type': "FPKM"})
            # fgresults["count_file"] = gene_count_file
            if isinstance(gene_main_id, types.StringTypes) or type(gene_main_id) == unicode:
                gene_main_id = ObjectId(gene_main_id)
            # collection.update({"task_id": task_id, "method": "RSEM", 'exp_level': 'G', 'exp_type': "FPKM"}, {
            #     '$set': {'count_file': gene_count_file}}, upsert=True)
            collection.update({"task_id": task_id, 'main_id':gene_main_id}, {
                '$set': {'count_file': gene_count_file}}, upsert=True)
        else:
            self.set_error('gene_fpkm/count_file is not found', code = "12001702")

        ## set output
        if gene_fpkm_anno_file is not None:
            os.link(gene_fpkm_anno_file, os.path.join(self.output_dir, os.path.basename(gene_fpkm_anno_file)))
        else:
            os.link(gene_fpkm_file, os.path.join(self.output_dir, os.path.basename(gene_fpkm_file)))
        if transcript_fpkm_anno_file is not None:
            os.link(transcript_fpkm_anno_file,
                    os.path.join(self.output_dir, os.path.basename(transcript_fpkm_anno_file)))
        else:
            if transcript_main_id != "":
                os.link(transcript_fpkm_file, os.path.join(self.output_dir, os.path.basename(transcript_fpkm_file)))

    def dump_tpm(self, task_id, transcript_main_id, gene_main_id):
        all_exp = self.api.api("denovo_rna_v2.all_exp")
        collection = all_exp.db['sg_exp']
        if transcript_main_id == "":
            results = collection.find({"task_id": task_id, "method": "RSEM", "exp_level": "G"})
        else:
            results = collection.find({"task_id": task_id, "method": "RSEM"})
        transcript_fpkm_file = None
        gene_fpkm_file = None
        gene_fpkm_anno_file=None
        transcript_fpkm_anno_file=None
        gene_count_file=None
        transcript_count_file=None
        for each in results:
            # 2019.04.10 bug each['count_file']也需要改名称，不然实际上只是把s3上的*.count.matrix.xls下载下来后，
            # 更改名称导表而已。正确的做法是将s3上的*.fpkm.matrix.xls下载下来，即指定正确的下载文件。
            if each['exp_level'] == 'T' and 'count_file' in each:
                local_t_base = os.path.basename(each['count_file']).replace('.count.matrix', '.tpm.matrix')
                # transcript_fpkm_file = self.download_s3_file(each['count_file'], local_t_base)
                transcript_fpkm_file = self.download_s3_file(
                    each['count_file'].replace('.count.matrix', '.tpm.matrix'), local_t_base
                )
                transcript_count_file = each['count_file']
                transcript_fpkm_anno_file = each['count_file'].replace('.count.matrix.xls', '.tpm.matrix.annot.xls')
                #20191204 bug 需要修改anno文件路径
                if exists(transcript_fpkm_anno_file):
                    transcript_fpkm_anno_file = self.download_s3_file(transcript_fpkm_anno_file,os.path.basename(transcript_fpkm_anno_file))
                else:
                    transcript_fpkm_anno_file = None
            if each['exp_level'] == 'G' and 'count_file' in each:
                local_g_base = os.path.basename(each['count_file']).replace('.count.matrix', '.tpm.matrix')
                # gene_fpkm_file =  self.download_s3_file(each['count_file'], local_g_base)
                gene_fpkm_file =  self.download_s3_file(
                    each['count_file'].replace('.count.matrix', '.tpm.matrix'), local_g_base
                )
                gene_count_file = each['count_file']
                gene_fpkm_anno_file = each['count_file'].replace('.count.matrix.xls', '.tpm.matrix.annot.xls')
                if exists(gene_fpkm_anno_file):
                    gene_fpkm_anno_file = self.download_s3_file(gene_fpkm_anno_file,os.path.basename(gene_fpkm_anno_file))
                else:
                    gene_fpkm_anno_file = None
            if transcript_fpkm_file and gene_fpkm_file:
                break

        if transcript_fpkm_file is not None:
            all_exp.add_exp(transcript_fpkm_file, quant_method="RSEM", exp_level='T',
                            main_id=transcript_main_id, exp_type='tpm', add_distribution=False, task_id=self.option('raw_task_id'))

            # collection.update({"task_id": task_id, "main_id": "RSEM", 'exp_level': 'T', 'exp_type': "TPM"}, {
            #     '$set': {'count_file': transcript_count_file}}, upsert=True)
            # ftresults = collection.find_one({"task_id": task_id, "method": "RSEM", 'exp_level': 'T', 'exp_type': "TPM"})
            # ftresults["count_file"] = transcript_count_file
            if isinstance(transcript_main_id, types.StringTypes) or type(transcript_main_id) == unicode:
                transcript_main_id = ObjectId(transcript_main_id)
            collection.update({"task_id": task_id, 'main_id':transcript_main_id}, {
                '$set': {'count_file': transcript_count_file}}, upsert=True)
        elif transcript_main_id == "":
            self.logger.info('only Unigene level')
        else:
            self.set_error('transcript_fpkm/count_file is not found', code = "12001701")

        if gene_fpkm_file is not None:
            all_exp.add_exp(gene_fpkm_file, quant_method="RSEM", exp_level='G',
                            main_id=gene_main_id, exp_type='TPM', add_distribution=False, task_id=self.option('raw_task_id'))
            # collection.update({"task_id": task_id, "method": "RSEM", 'exp_level': 'G', 'exp_type': "TPM"}, {
            #     '$set': {'count_file': gene_count_file}}, upsert=True)
            if isinstance(gene_main_id, types.StringTypes) or type(gene_main_id) == unicode:
                gene_main_id = ObjectId(gene_main_id)
            collection.update({"task_id": task_id, 'main_id':gene_main_id}, {
                '$set': {'count_file': gene_count_file}}, upsert=True)
            # fgresults = collection.find_one({"task_id": task_id, "method": "RSEM", 'exp_level': 'G', 'exp_type': "TPM"})
            # fgresults["count_file"] = gene_count_file
        else:
            self.set_error('gene_tpm/count_file is not found', code = "12001702")

        ## set output
        if gene_fpkm_anno_file is not None:
            os.link(gene_fpkm_anno_file, os.path.join(self.output_dir, os.path.basename(gene_fpkm_anno_file)))
        else:
            os.link(gene_fpkm_file, os.path.join(self.output_dir, os.path.basename(gene_fpkm_file)))
        if transcript_fpkm_anno_file is not None:
            os.link(transcript_fpkm_anno_file,os.path.join(self.output_dir, os.path.basename(transcript_fpkm_anno_file)))
        else:
            if transcript_main_id != "":
                os.link(transcript_fpkm_file, os.path.join(self.output_dir, os.path.basename(transcript_fpkm_file)))


    def paste_annotation(self):
        all_exp = self.api.api("denovo_rna_v2.all_exp")
        conn = all_exp.db['sg_annotation_stat']
        task_id = self.option('task_id')
        try:
            find_result = conn.find_one({"task_id": task_id, "type": "latest", "status": "end"})
        except:
            self.set_error('annotation result file cannot found', code = "12001704")
        else:
            if find_result is None:
                find_result = conn.find_one({"task_id": task_id, "type": "origin"})
        gene_annot = os.path.dirname(find_result['result_dir']) + '/AnnoQuery/unigene_anno_detail.xls'
        trans_annot = os.path.dirname(find_result['result_dir']) + '/AnnoQuery/transcript_anno_detail.xls'
        gene_annot = self.download_s3_file(gene_annot, "unigene_anno_detail.xls")
        trans_annot = self.download_s3_file(trans_annot, "transcript_anno_detail.xls")
        gene_annot = pd.read_table(gene_annot, header=0, index_col=0)
        trans_annot = pd.read_table(trans_annot, header=0, index_col=0)
        # for gene
        exp_pd = pd.read_table(self.tool.output_dir + '/gene.tpm.matrix', header=0, sep='\t', index_col=0)
        exp_pd = exp_pd.join(gene_annot, how='left')
        out_exp = self.tool.output_dir + '/gene.tpm.matrix' + '.annot.xls'
        exp_pd.to_csv(out_exp, sep='\t', header=True, index=True)
        # for transcript
        exp_pd = pd.read_table(self.tool.output_dir + '/transcript.tpm.matrix', header=0, sep='\t', index_col=0)
        exp_pd = exp_pd.join(trans_annot, how='left')
        out_exp = self.tool.output_dir + '/transcript.tpm.matrix' + '.annot.xls'
        exp_pd.to_csv(out_exp, sep='\t', header=True, index=True)
