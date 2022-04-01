# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modifiy = modified 2019.09.15

from biocluster.module import Module
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from mbio.workflows.denovo_rna_v2.denovo_test_api import DenovoTestApiWorkflow
from biocluster.wsheet import Sheet
import os
import re
import shutil
import json
import glob
import pandas as pd
from biocluster.config import Config
from mbio.packages.lnc_rna.copy_file import CopyFile
import unittest


class TargetMirnaModule(Module):
    """
    交互分析进行筛选blast参数重注释时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        #self.rpc = False
        super(TargetMirnaModule, self).__init__(wsheet_object)
        options = [
            {"name": "novol", "type": "infile", "format": "small_rna.fasta"},  # 输入文件
            {"name": "known", "type": "infile", "format": "small_rna.fasta"},  # 输入文件
            {"name": "m_novol", "type": "infile", "format": "small_rna.fasta"},  # 输入文件
            {"name": "m_known", "type": "infile", "format": "small_rna.fasta"},
            {"name": "l_novol", "type": "infile", "format": "small_rna.fasta"},  # 输入文件
            {"name": "l_known", "type": "infile", "format": "small_rna.fasta"},
            {"name": "c_novol", "type": "infile", "format": "small_rna.fasta"},  # 输入文件
            {"name": "c_known", "type": "infile", "format": "small_rna.fasta"},


            {"name": "ref", "type": "infile", "format": "small_rna.fasta"},  # 输入文件
            {"name": "type", "type": "string", "default": "animal"},
            {"name": "anno_detail", "type": "string", "default": None},
            {"name": "lnc_ref_detail", "type": "infile", "format": "small_rna.common"},
            {"name": "lnc_new_detail", "type": "infile", "format": "small_rna.common"},
            {"name": "circ_detail", "type": "string", "default": None},
            {"name": "outtable", "type": "outfile", "format": "small_rna.common"},
            {"name": "species_name", "type": "string", "default": None},
            {"name": "stat_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "origin_result", "type": "string"},
            {"name": "last_id", "type": "string"},
            {"name": "target_id", "type": "string"},
            {"name": "last_id_target", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "trans2gene", "type": "string"},
            {"name": "taxonomy", "type": "string", "default": "Animal"},
            {"name": "version", "type": "string", "default": "v1"},
            {"name": "anno_path", "type": "string", "default": None},

            {"name": "m_method", "type": "string", "default": "miranda"},
            {"name": "m_miranda_score", "type": "string", "default": "160"},
            {"name": "m_miranda_energy", "type": "string", "default": "-20"},
            {"name": "m_miranda_strict", "type": "string", "default": "on"},
            {"name": "m_rnahybird_num", "type": "string", "default": "100"},
            {"name": "m_rnahybird_energy", "type": "string", "default": "-20"},
            {"name": "m_rnahybird_pvalue", "type": "string", "default": "0.01"},
            {"name": "m_ps_robot_score", "type": "string", "default": "2.5"},
            {"name": "m_targetfinder_score", "type": "string", "default": "4"},
            {"name": "m_min_support", "type": "int", "default": 1},

            {"name": "l_method", "type": "string", "default": "miranda"},
            {"name": "l_miranda_score", "type": "string", "default": "160"},
            {"name": "l_miranda_energy", "type": "string", "default": "-20"},
            {"name": "l_miranda_strict", "type": "string", "default": "on"},
            {"name": "l_rnahybird_num", "type": "string", "default": "100"},
            {"name": "l_rnahybird_energy", "type": "string", "default": "-20"},
            {"name": "l_rnahybird_pvalue", "type": "string", "default": "0.01"},
            {"name": "l_ps_robot_score", "type": "string", "default": "2.5"},
            {"name": "l_targetfinder_score", "type": "string", "default": "4"},
            {"name": "l_min_support", "type": "int", "default": 1},

            {"name": "c_method", "type": "string", "default": "miranda"},
            {"name": "c_miranda_score", "type": "string", "default": "160"},
            {"name": "c_miranda_energy", "type": "string", "default": "-20"},
            {"name": "c_miranda_strict", "type": "string", "default": "on"},
            {"name": "c_rnahybird_num", "type": "string", "default": "100"},
            {"name": "c_rnahybird_energy", "type": "string", "default": "-20"},
            {"name": "c_rnahybird_pvalue", "type": "string", "default": "0.01"},
            {"name": "c_ps_robot_score", "type": "string", "default": "2.5"},
            {"name": "c_targetfinder_score", "type": "string", "default": "4"},
            {"name": "c_min_support", "type": "int", "default": 1},

        ]
        self.add_option(options)


        self.target_predict_m = self.add_module("small_rna.target_predict")
        self.target_predict_l = self.add_module("small_rna.target_predict")
        self.target_predict_c = self.add_module("small_rna.target_predict")
        self.target_predict_mn = self.add_module("small_rna.target_predict")
        self.target_predict_ln = self.add_module("small_rna.target_predict")
        self.target_predict_cn = self.add_module("small_rna.target_predict")

        # self.target_predict.on('end', self.set_db)

    def run_target_m(self):
        '''
        if os.path.exists(self.option("anno_path") + '/anno_stat/all_anno_detail.xls'):
            anno_detail = self.option("anno_path") + '/anno_stat/all_anno_detail.xls'
        else:
            anno_detail = self.option("anno_path") + '/anno_stat/trans_anno_detail.xls'
        '''

        options = {
            "novol": self.option("novol").prop['path'],
            "known": self.option("known").prop['path'],
            "ref": self.option("m_known").prop['path'],
            "method": self.option("m_method"),
            "anno_detail": self.option("anno_detail"),
            "version": "v1.2",
            "species": self.option("species_name"),
            "type": self.option("taxonomy").lower(),
            "min_support": self.option("m_min_support"),
            "miranda_score": self.option("m_miranda_score"),
            "miranda_energy": self.option("m_miranda_energy"),
            "miranda_strict": self.option("m_miranda_strict"),
            "rnahybird_num": self.option("m_rnahybird_num"),
            "rnahybird_energy": self.option("m_rnahybird_energy"),
            "rnahybird_pvalue": self.option("m_rnahybird_pvalue"),
            "ps_robot_score": self.option("m_ps_robot_score"),
            "targetfinder_score": self.option("m_targetfinder_score"),
        }

        self.target_predict_m.on("end", self.set_output, "m_known")
        self.target_predict_m.set_options(options)
        self.target_predict_m.run()


    def run_target_l(self):
        '''
        if os.path.exists(self.option("anno_path") + '/anno_stat/all_anno_detail.xls'):
            anno_detail = self.option("anno_path") + '/anno_stat/all_anno_detail.xls'
        else:
            anno_detail = self.option("anno_path") + '/anno_stat/trans_anno_detail.xls'
        '''


        options = {
            "novol": self.option("novol").prop['path'],
            "known": self.option("known").prop['path'],
            "ref": self.option("l_known").prop['path'],
            "method": self.option("l_method"),
            "anno_detail": self.work_dir + "/lnc_detail.xls",
            "version": "v1.2",
            "species": self.option("species_name"),
            "type": self.option("taxonomy").lower(),
            "min_support": self.option("l_min_support"),
            "miranda_score": self.option("l_miranda_score"),
            "miranda_energy": self.option("l_miranda_energy"),
            "miranda_strict": self.option("l_miranda_strict"),
            "rnahybird_num": self.option("l_rnahybird_num"),
            "rnahybird_energy": self.option("l_rnahybird_energy"),
            "rnahybird_pvalue": self.option("l_rnahybird_pvalue"),
            "ps_robot_score": self.option("l_ps_robot_score"),
            "targetfinder_score": self.option("l_targetfinder_score"),
        }

        self.target_predict_l.on("end", self.set_output, "l_known")
        self.target_predict_l.set_options(options)
        self.target_predict_l.run()

    def run_target_c(self):
        '''
        if os.path.exists(self.option("anno_path") + '/anno_stat/all_anno_detail.xls'):
            anno_detail = self.option("anno_path") + '/anno_stat/all_anno_detail.xls'
        else:
            anno_detail = self.option("anno_path") + '/anno_stat/trans_anno_detail.xls'
        '''

        options = {
            "novol": self.option("novol").prop['path'],
            "known": self.option("known").prop['path'],
            "ref": self.option("c_known").prop['path'],
            "method": self.option("c_method"),
            "anno_detail": self.option("anno_detail"),
            "version": "v1.2",
            "species": self.option("species_name"),
            "type": self.option("taxonomy").lower(),
            "min_support": self.option("c_min_support"),
            "miranda_score": self.option("c_miranda_score"),
            "miranda_energy": self.option("c_miranda_energy"),
            "miranda_strict": self.option("c_miranda_strict"),
            "rnahybird_num": self.option("c_rnahybird_num"),
            "rnahybird_energy": self.option("c_rnahybird_energy"),
            "rnahybird_pvalue": self.option("c_rnahybird_pvalue"),
            "ps_robot_score": self.option("c_ps_robot_score"),
            "targetfinder_score": self.option("c_targetfinder_score"),
            "circ_detail": self.option("circ_detail")
        }

        self.target_predict_c.on("end", self.set_output, "c_known")
        self.target_predict_c.set_options(options)
        self.target_predict_c.run()


    def run_target_mn(self):
        '''
        if os.path.exists(self.option("anno_path") + '/anno_stat/all_anno_detail.xls'):
            anno_detail = self.option("anno_path") + '/anno_stat/all_anno_detail.xls'
        else:
            anno_detail = self.option("anno_path") + '/anno_stat/trans_anno_detail.xls'
        '''

        options = {
            "novol": self.option("novol").prop['path'],
            "known": self.option("known").prop['path'],
            "ref": self.option("m_novol").prop['path'],
            "method": self.option("m_method"),
            "anno_detail": self.option("anno_detail"),
            "version": "v1.2",
            "species": self.option("species_name"),
            "type": self.option("taxonomy").lower(),
            "min_support": self.option("m_min_support"),
            "miranda_score": self.option("m_miranda_score"),
            "miranda_energy": self.option("m_miranda_energy"),
            "miranda_strict": self.option("m_miranda_strict"),
            "rnahybird_num": self.option("m_rnahybird_num"),
            "rnahybird_energy": self.option("m_rnahybird_energy"),
            "rnahybird_pvalue": self.option("m_rnahybird_pvalue"),
            "ps_robot_score": self.option("m_ps_robot_score"),
            "targetfinder_score": self.option("m_targetfinder_score"),
        }

        self.target_predict_mn.on("end", self.set_output, "m_novel")
        self.target_predict_mn.set_options(options)
        self.target_predict_mn.run()

    def run_target_ln(self):
        '''
        if os.path.exists(self.option("anno_path") + '/anno_stat/all_anno_detail.xls'):
            anno_detail = self.option("anno_path") + '/anno_stat/all_anno_detail.xls'
        else:
            anno_detail = self.option("anno_path") + '/anno_stat/trans_anno_detail.xls'
        '''

        options = {
            "novol": self.option("novol").prop['path'],
            "known": self.option("known").prop['path'],
            "ref": self.option("l_novol").prop['path'],
            "method": self.option("l_method"),
            "anno_detail": self.option("anno_detail"),
            "version": "v1.2",
            "species": self.option("species_name"),
            "type": self.option("taxonomy").lower(),
            "min_support": self.option("l_min_support"),
            "miranda_score": self.option("l_miranda_score"),
            "miranda_energy": self.option("l_miranda_energy"),
            "miranda_strict": self.option("l_miranda_strict"),
            "rnahybird_num": self.option("l_rnahybird_num"),
            "rnahybird_energy": self.option("l_rnahybird_energy"),
            "rnahybird_pvalue": self.option("l_rnahybird_pvalue"),
            "ps_robot_score": self.option("l_ps_robot_score"),
            "targetfinder_score": self.option("l_targetfinder_score"),
        }

        self.target_predict_ln.on("end", self.set_output, "l_novel")
        self.target_predict_ln.set_options(options)
        self.target_predict_ln.run()

    def run_target_cn(self):
        '''
        if os.path.exists(self.option("anno_path") + '/anno_stat/all_anno_detail.xls'):
            anno_detail = self.option("anno_path") + '/anno_stat/all_anno_detail.xls'
        else:
            anno_detail = self.option("anno_path") + '/anno_stat/trans_anno_detail.xls'
        '''

        options = {
            "novol": self.option("novol").prop['path'],
            "known": self.option("known").prop['path'],
            "ref": self.option("c_novol").prop['path'],
            "method": self.option("c_method"),
            "anno_detail": self.option("anno_detail"),
            "version": "v1.2",
            "species": self.option("species_name"),
            "type": self.option("taxonomy").lower(),
            "min_support": self.option("c_min_support"),
            "miranda_score": self.option("c_miranda_score"),
            "miranda_energy": self.option("c_miranda_energy"),
            "miranda_strict": self.option("c_miranda_strict"),
            "rnahybird_num": self.option("c_rnahybird_num"),
            "rnahybird_energy": self.option("c_rnahybird_energy"),
            "rnahybird_pvalue": self.option("c_rnahybird_pvalue"),
            "ps_robot_score": self.option("c_ps_robot_score"),
            "targetfinder_score": self.option("c_targetfinder_score"),
            "circ_detail": self.option("circ_detail")
        }

        self.target_predict_cn.on("end", self.set_output, "c_novel")
        self.target_predict_cn.set_options(options)
        self.target_predict_cn.run()


    def get_lnc_detail(self):
        all_df = pd.DataFrame()
        if self.option("lnc_ref_detail").is_set:
            ref = pd.read_table(self.option("lnc_ref_detail").prop['path'] , sep="\t", header=0)

            ref = ref.fillna("-").loc[:, ['lncrna_id', 'gene_id', 'gene_name', 'gene_description']].rename(columns={"lncrna_id": "transcript_id"})
            all_df = ref
            if self.option("lnc_new_detail").is_set:
                new = pd.read_table(self.option("lnc_new_detail").prop['path'], sep="\t", header=0)
                new = new.fillna("-").loc[:, ['transcript_id', 'gene_id', 'gene_name', 'gene_description']]

            all_df = ref.append(new)

        all_df.to_csv(self.work_dir + "/lnc_detail.xls", sep="\t", index=False)

    def run(self):
        super(TargetMirnaModule, self).run()
        self.get_lnc_detail()

        rely_list = []

        if self.option("m_known").is_set:
            rely_list.append(self.target_predict_m)
        if self.option("m_novol").is_set:
            rely_list.append(self.target_predict_mn)
        if self.option("l_known").is_set:
            rely_list.append(self.target_predict_l)
        if self.option("l_novol").is_set:
            rely_list.append(self.target_predict_ln)
        if self.option("c_known").is_set:
            rely_list.append(self.target_predict_c)
        if self.option("c_novol").is_set:
            rely_list.append(self.target_predict_cn)

        self.on_rely(rely_list, self.set_db)

        if self.option("m_known").is_set:
            self.run_target_m()
        if self.option("m_novol").is_set:
            self.run_target_mn()
        if self.option("l_known").is_set:
            self.run_target_l()
        if self.option("l_novol").is_set:
            self.run_target_ln()
        if self.option("c_known").is_set:
            self.run_target_c()
        if self.option("c_novol").is_set:
            self.run_target_cn()

    def set_output(self, event):
        obj = event["bind_object"]
        name = event['data']
        CopyFile().linkdir(obj.output_dir, os.path.join(self.output_dir, name))


    def set_db(self):
        self.end()

    def end(self):

        target_dir = self.output_dir


        repaths = [
            [".", "", "靶基因"],
        ]

        super(TargetMirnaModule, self).end()

if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        def test(self):
            import random
            from mbio.workflows.single import SingleWorkflow
            from biocluster.wsheet import Sheet

            test_dir  = "/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/whole_rna/"

            data = {
                'id': 'whole_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
                'type': "module",
                'name': 'whole_transcriptome.target_mirna',
                'options': {
                    "novol" : test_dir + "novel_mi.fa",
                    "known" : test_dir + "known_mi.fa",
                    "m_novol": test_dir + "novel_m.fa",
                    "m_known": test_dir + "known_m.fa",
                    "l_novol": test_dir + "novel_l.fa",
                    "l_known": test_dir + "known_l.fa",
                    "anno_detail": "/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/AnnotMerge/output/allannot_class/all_annot.xls",
                    "lnc_ref_detail": "/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/LargeGush/output/known_lnc_identify/known_lncrna_detail.xls",
                    "lnc_new_detail": "/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/LargeGush/output/new_lncrna_predict/novel_lncrna_predict_detail.xls"
                }
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()

    unittest.main()
