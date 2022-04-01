# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import shutil
import types
import json
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError
from mbio.packages.lnc_rna.copy_file import CopyFile
import gzip



class ProteinClassWorkflow(Workflow):
    """
    在目标DNA序列中搜寻开发阅读框，可以输出每个ORF所在的区域，并翻译成对应的蛋白序列
    此工具可以为新预测的DNA序列查找潜在的蛋白编码区
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ProteinClassWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "input", "type": "infile", "format": "sequence.fasta"},  # FASTA序列文件
            {"name": "pfam", "type": "string", "default": "yes"},
            {"name": "version", "type": "string", "default": "27"},
            {"name": "e_value", "type": "float", "default": 1e-3},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        print options
        self.add_option(options)
        self.hmm_tools = list()
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def check_options(self):
        if not self.option("input").is_set:
            raise OptionError('请输入FASTA文件')
        return True

    def run(self):
        self.run_splitfasta()
        super(ProteinClassWorkflow, self).run()

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run_splitfasta(self):
        self.splitfasta = self.add_tool("sequence.split_fasta")
        self.fasta_name = self.option("input").prop["path"]

        self.splitfasta.set_options({
            "fasta": self.fasta_name,
            "lines": 5000,
        })
        self.splitfasta.on('end', self.run_hmm)
        self.splitfasta.run()

    def run_hmm(self):
        opts = {
            "e_value": self.option("e_value"),
        }
        for f in os.listdir(self.splitfasta.output_dir):
            opts['pep'] = os.path.join(self.splitfasta.output_dir, f)
            hmm_tool = self.add_tool('denovo_rna_v2.hmm')
            hmm_tool.set_options(opts)
            # hmm_tool.run()
            self.hmm_tools.append(hmm_tool)
        if len(self.hmm_tools) == 1:
            self.hmm_tools[0].on("end", self.run_cathmmout)
        else:
            self.on_rely(self.hmm_tools, self.run_cathmmout)
        for tool in self.hmm_tools:
            tool.run()

    def run_cathmmout(self):
        gene2pfam = dict()

        def read_table_domtblout(target, gene2pfam):
            with open(target) as f:
                all_lines = f.readlines()
                for line in all_lines[1:]:
                    cols= line.strip().split("\t")
                    if cols[0] in gene2pfam:
                        pass
                    else:
                        gene2pfam[cols[0]] = cols[2].split(".")[0]
            return all_lines[1:]
        target_files = [x.work_dir + "/" + 'pfam_domain' for x in
                        self.hmm_tools]  # 处理所有的tool
        # 处理所有的pfam.domtblout文件
        # 这里任意选择一个tool的文件，不用使用enumerate来判断索引位置
        with open(target_files[0]) as f:
            tmp = f.readlines()
            head_lines = tmp[0]
        target_list = list()
        for x in target_files:
            target_list += read_table_domtblout(x, gene2pfam) # 这里不能使用append，
            # +=是把所有的变为一个列表，append加进来的列表还是列表
        target_list = [head_lines] + target_list
        output_path = self.work_dir + "/"
        with open(output_path + "pfam_domain", 'w') as f:
            for each in target_list:
                f.write(each)

        pfam_clan_file = self.config.SOFTWARE_DIR + "/database/Annotation/other2019/pfam32/Pfam-A.clans.tsv.gz"
        pfam2class = dict()
        pfam2clan = dict()
        with gzip.open(pfam_clan_file, 'rb') as f:
            for line in f:
                cols = line.strip().split("\t")
                if cols[3] != "":
                    pfam2class[cols[0]] = cols[3]
                if cols[2] != "":
                    pfam2clan[cols[0]] = cols[2]

        # self.logger.info("pfam_class{}".format(pfam2class))

        pfam2class_num = dict()
        pfam2clan_num = dict()
        # self.logger.info("gene2pfam{}".format(gene2pfam))
        for g,p in gene2pfam.items():
            if p in pfam2class:
                # print p
                if pfam2class[p] in pfam2class_num:
                    pfam2class_num[pfam2class[p]] += 1
                else:
                    pfam2class_num[pfam2class[p]] = 1
            if p in pfam2clan:
                if pfam2clan[p] in pfam2clan_num:
                    pfam2clan_num[pfam2clan[p]] += 1
                else:
                    pfam2clan_num[pfam2clan[p]] = 1

        with open(output_path + "pfam_domain.class_stat.xls", 'w') as f:
            f.write("class_name\tnum\n")
            for k in sorted(pfam2class_num.keys(), key=lambda x:pfam2class_num[x], reverse=True):
                f.write("{}\t{}\n".format(k, pfam2class_num[k]))

        with open(output_path + "pfam_domain.clan_stat.xls", 'w') as f:
            f.write("clan_name\tnum\n")
            for k in sorted(pfam2clan_num.keys(), key=lambda x:pfam2clan_num[x], reverse=True):
                f.write("{}\t{}\n".format(k, pfam2clan_num[k]))
        self.set_output()



    def set_db(self):
        protein_class_api = self.api.api("tool_lab.protein_class")
        protein_class_api.add_class_detail(self.option("main_id"), os.path.join(self.work_dir, "pfam_domain.class_stat.xls"))
        table_dict = {
            "column": [
                {"field": "class_name", "title": "class_name", "filter": "false", "sort": "false", "type": "string"},
                {"field": "num", "title": "num", "filter": "false", "sort": "false", "type": "int"}
            ],
            "condition": {}
        }
        table_info = json.dumps(table_dict, sort_keys=True, separators=(',', ':'))

        protein_class_api.update_db_record('protein_class',
                                           query_dict={"main_id": ObjectId(self.option("main_id"))},
                                           update_dict={'status': 'end',
                                                        'table_data': table_info})


    def set_output(self):
        self.set_db()
        for file in ['pfam_domain', 'pfam_domain.clan_stat.xls', 'pfam_domain.class_stat.xls']:
            CopyFile().linkfile(os.path.join(self.work_dir, file),  os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "ORF查找结果文件",0],
            [r'.*.fa', 'xls', 'ORF查找结果文件', 0],
        ])
        super(ProteinClassWorkflow, self).end()
