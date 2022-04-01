# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# last_modify:20200609

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import unittest
import gevent.subprocess as subprocess
import json
import glob
from collections import OrderedDict
from biocluster.config import Config
import pandas as pd


class DiffGenesetAllModule(Module):
    """
    该Module用于基因融合分析，默认使用方法
    """
    def __init__(self, work_id):
        super(DiffGenesetAllModule, self).__init__(work_id)
        options = [
            # {"name": "geneset_names", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "diff_path", "type": "string", "default": None},
            {"name": "annot_result", "type": "string", "default": None},
            {"name": "diff_method", "type": "string"},
            {"name": "transcript_exp_file", "type": "string"},
            {"name": "level", "type": "string", "default": "G"},
            {"name": "kegg_version", "type": "string", "default": None},
            {"name": "reactome_version", "type": "string", "default": None},
            {"name": "species", "type": "string", "default": "Homo_sapiens"}
        ]
        self.add_option(options)
        self.geneset_prepare_tools=[]
        self.genest_infos = OrderedDict()
        self.all_genesets = OrderedDict()
        self.common = self.add_tool("whole_transcriptome_v1_2.diffgt4work.annot_prepare")
        self.total_dict = OrderedDict()
        self.run_genesets_analysis = True
        self.diff_geneset_prepare = ""
        self.selecte_geneset = ""


    def check_options(self):
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def get_genesets_info(self):
        diff_dir = self.option("diff_path")
        map_dict = dict()
        self.all_genesets = map_dict
        category = 'mRNA'
        level = self.option("level")
        for fname in os.listdir(diff_dir):
            if fname.endswith('.detail.txt'):
                map_dict[fname[:-11]] = os.path.join(diff_dir, fname)
        for vs_pair, detail_table in map_dict.items():
            df = pd.read_table(detail_table)
            df = df[df['significant'] == 'yes']
            length = len(df)
            geneset_name = '{}_{}_{}'.format(level, category, vs_pair)
            regulate = "all"
            compare = list(df["compare"])[0]
            self.genest_infos[geneset_name] = {
                "geneset_name": geneset_name,
                "compare": compare,
                "gene_num" :length,
                "regulate": regulate,
                "compare_path": detail_table
            }
        self.prepare_common_file()

    def prepare_common_file(self):
        self.logger.info("开始准备差异数据挖掘公共文件")
        opts = {
            "annot_result": self.option("annot_result"),
            "level": self.option("level"),
            "kegg_version": self.option("kegg_version"),
            "transcript_exp_file": self.option("transcript_exp_file"),
            "species": self.option("species"),
        }
        self.common.set_options(opts)
        self.common.on("end", self.prepare_geneset_infos)
        self.common.run()

    def get_selected_genset_info(self):
        geneset_num = dict()
        for geneset in self.genest_infos:
            gene_num = self.genest_infos[geneset]["gene_num"]
            geneset_num[geneset] = gene_num
        def min_geneset_select(geneset_num):
            max_geneset_num = max([geneset_num[x] for x in geneset_num])
            if max_geneset_num <= 100:
                return False
            else:
                select_geneset = min(geneset_num, key=lambda x: geneset_num[x] if geneset_num[x] >= 100 else 10000000)
                return select_geneset
        select_geneset = min_geneset_select(geneset_num)
        if not select_geneset:
            self.run_genesets_analysis = False
        return select_geneset



    def prepare_geneset_infos(self):
        self.selecte_geneset= self.get_selected_genset_info()
        if not self.selecte_geneset:
            self.set_output()
        else:
            self.diff_geneset_prepare = self.add_tool("whole_transcriptome_v1_2.diffgt4work.diff_geneset_prepare")
            opts = {
                "annot_result": self.option("annot_result"),
                "compare": self.genest_infos[self.selecte_geneset]["compare"],
                "regulate": self.genest_infos[self.selecte_geneset]["regulate"],
                "geneset_name": self.genest_infos[self.selecte_geneset]["geneset_name"],
                "compare_path": self.genest_infos[self.selecte_geneset]["compare_path"],
                "level": self.option("level"),
                "species": self.option("species"),
            }
            self.diff_geneset_prepare.set_options(opts)
            self.diff_geneset_prepare.on('end', self.set_output)
            self.diff_geneset_prepare.run()

    def set_output(self):
        if not self.run_genesets_analysis :
            self.total_dict["has_geneset_pipline_results"] = False
        else:
            self.total_dict["all_genesets"] = self.all_genesets
            common_file_path = self.common.option("common_file_json").prop["path"]
            common_file_dict = json.load(open(self.common.option("common_file_json").prop["path"]))
            self.total_dict["common_file"] = common_file_dict
            self.total_dict["genesets"] =OrderedDict()
            geneset_json = self.diff_geneset_prepare.option("geneset_json").prop["path"]
            geneset_dict = json.load(open(geneset_json))
            geneset_name = geneset_dict["geneset_name"]
            self.total_dict["genesets"][geneset_name] = geneset_dict
        with open(os.path.join(self.output_dir,"prepare_json"),"w") as f:
            json.dump(self.total_dict,f,indent=2)
        self.end()


    def run(self):
        super(DiffGenesetAllModule, self).run()
        self.logger.info("开始运行差异基因集数据挖掘数据准备")
        self.logger.info("首先解析基因集列表文件")
        self.get_genesets_info()

    def end(self):
        super(DiffGenesetAllModule, self).end()



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
            "id": "geneset_prepare" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "medical_transcriptome.workflow_diffgt.diff_geneset_all",
            "instant": False,
            "options": dict(
                diff_path= "/mnt/ilustre/users/sanger-dev/workspace/20200923/MedicalTranscriptome_medical_transcriptome_workflow_20200923_153802009/DiffexpBatch/output",
                annot_result ="/mnt/ilustre/users/sanger-dev/workspace/20200923/MedicalTranscriptome_medical_transcriptome_workflow_20200923_153802009/AnnotMerge__1/output",
                diff_method="DESeq2",
                kegg_version ="202007",
                reactome_version="202007",
                level = "G",
                gene_count_file ="/mnt/ilustre/users/sanger-dev/workspace/20200923/MedicalTranscriptome_medical_transcriptome_workflow_20200923_153802009/Quant/output/gene.count.matrix",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()