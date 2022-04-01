# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modifiy = modified 2017.12.23

from biocluster.workflow import Workflow
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
from biocluster.config import Config
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.whole_transcriptome.catalogue import mirna
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class TargetMiWorkflow(Workflow):
    """
    交互分析进行small rna靶基因预测
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        #self.rpc = False
        super(TargetMiWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "novol", "type": "infile", "format": "small_rna.fasta"},  # 输入文件
            {"name": "known", "type": "infile", "format": "small_rna.fasta"},  # 输入文件

            {"name": "m_novol", "type": "infile", "format": "small_rna.fasta"},  # 输入文件
            {"name": "m_known", "type": "infile", "format": "small_rna.fasta"},
            {"name": "l_novol", "type": "infile", "format": "small_rna.fasta"},  # 输入文件
            {"name": "l_known", "type": "infile", "format": "small_rna.fasta"},
            {"name": "c_novol", "type": "infile", "format": "small_rna.fasta"},  # 输入文件
            {"name": "c_known", "type": "infile", "format": "small_rna.fasta"},

            {"name": "type", "type": "string", "default": "animal"},
            {"name": "anno_detail", "type": "infile", "format": "small_rna.common"},
            {"name": "circ_detail", "type": "string", "default": None},
            {"name": "outtable", "type": "outfile", "format": "small_rna.common"},
            {"name": "species_name", "type": "string", "default": None},
            {"name": "task_id", "type": "string"},
            {"name": "origin_result", "type": "string"},
            {"name": "target_id", "type": "string"},
            {"name": "last_id_target", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "trans2gene", "type": "string"},
            {"name": "taxonomy", "type": "string", "default": None},
            {"name": "version", "type": "string", "default": "v1"},
            {"name": "anno_path", "type": "string", "default": None},
            {"name": "lnc_ref_detail", "type": "infile", "format": "small_rna.common"},
            {"name": "lnc_new_detail", "type": "infile", "format": "small_rna.common"},


            {"name": "m_miranda_score", "type": "string", "default": "160"},
            {"name": "m_miranda_energy", "type": "string", "default": "-20"},
            {"name": "m_miranda_strict", "type": "string", "default": "on"},
            {"name": "m_rnahybird_num", "type": "string", "default": "100"},
            {"name": "m_rnahybird_energy", "type": "string", "default": "-20"},
            {"name": "m_rnahybird_pvalue", "type": "string", "default": "0.01"},
            {"name": "m_ps_robot_score", "type": "string", "default": "2.5"},
            {"name": "m_targetfinder_score", "type": "string", "default": "4"},
            {"name": "m_min_support", "type": "int", "default": 1},

            {"name": "l_miranda_score", "type": "string", "default": "160"},
            {"name": "l_miranda_energy", "type": "string", "default": "-20"},
            {"name": "l_miranda_strict", "type": "string", "default": "on"},
            {"name": "l_rnahybird_num", "type": "string", "default": "100"},
            {"name": "l_rnahybird_energy", "type": "string", "default": "-20"},
            {"name": "l_rnahybird_pvalue", "type": "string", "default": "0.01"},
            {"name": "l_ps_robot_score", "type": "string", "default": "2.5"},
            {"name": "l_targetfinder_score", "type": "string", "default": "4"},
            {"name": "l_min_support", "type": "int", "default": 1},

            {"name": "c_miranda_score", "type": "string", "default": "160"},
            {"name": "c_miranda_energy", "type": "string", "default": "-20"},
            {"name": "c_miranda_strict", "type": "string", "default": "on"},
            {"name": "c_rnahybird_num", "type": "string", "default": "100"},
            {"name": "c_rnahybird_energy", "type": "string", "default": "-20"},
            {"name": "c_rnahybird_pvalue", "type": "string", "default": "0.01"},
            {"name": "c_ps_robot_score", "type": "string", "default": "2.5"},
            {"name": "c_targetfinder_score", "type": "string", "default": "4"},
            {"name": "c_min_support", "type": "int", "default": 1}

        ]
        for method in [
            'm_miranda', 'm_targetscan', 'm_psrobot', 'm_targetfinder', 'm_rnahybrid',
            'l_miranda', 'l_targetscan', 'l_psrobot', 'l_targetfinder', 'l_rnahybrid',
            'c_miranda', 'c_targetscan', 'c_psrobot', 'c_targetfinder', 'c_rnahybrid']:
            options.append({
                "name": method, "type": "string", "default": "no"
            })
        self.add_option(options)

        self.set_options(self._sheet.options())
        self.target_predict = self.add_module("whole_transcriptome.target_mirna")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/03 Target_Predict_Analysis/01 miRNA_Target')
        self.inter_dirs = []
        # self.target_predict.on('end', self.set_db)

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
        super(TargetMiWorkflow, self).send_log(data)

    def run_target(self):
        '''
        if os.path.exists(self.option("anno_path") + '/anno_stat/all_anno_detail.xls'):
            anno_detail = self.option("anno_path") + '/anno_stat/all_anno_detail.xls'
        else:
            anno_detail = self.option("anno_path") + '/anno_stat/trans_anno_detail.xls'
        '''

        m_method = list()
        l_method = list()
        c_method = list()
        for method in [
            'miranda', 'targetscan', 'psrobot', 'targetfinder', 'rnahybrid'
        ]:
            if self.option("m_" +method) == "yes":
                m_method.append(method)
            if self.option("l_" +method) == "yes":
                l_method.append(method)
            if self.option("c_" +method) == "yes":
                c_method.append(method)



        options = {
            "novol": self.option("novol").prop['path'],
            "known": self.option("known").prop['path'],
            "m_novol": self.option("m_novol"),
            "m_known": self.option("m_known"),
            "l_novol": self.option("l_novol"),
            "l_known": self.option("l_known"),
            "c_novol": self.option("c_novol"),
            "c_known": self.option("c_known"),

            "anno_detail":  self.option("anno_detail").prop['path'],
            "circ_detail":  self.option("circ_detail"),
            "lnc_ref_detail":  self.option("lnc_ref_detail"),
            "lnc_new_detail":  self.option("lnc_new_detail"),

            "taxonomy": self.option("taxonomy"),
            "species_name": self.option("species_name"),

            "m_method": ",".join(m_method),
            "m_min_support": self.option("m_min_support"),
            "m_miranda_score": self.option("m_miranda_score"),
            "m_miranda_energy": self.option("m_miranda_energy"),
            "m_miranda_strict": self.option("m_miranda_strict"),
            "m_rnahybird_num": self.option("m_rnahybird_num"),
            "m_rnahybird_energy": self.option("m_rnahybird_energy"),
            "m_rnahybird_pvalue": self.option("m_rnahybird_pvalue"),
            "m_ps_robot_score": self.option("m_ps_robot_score"),
            "m_targetfinder_score": self.option("m_targetfinder_score"),

            "l_method": ",".join(l_method),
            "l_min_support": self.option("l_min_support"),
            "l_miranda_score": self.option("l_miranda_score"),
            "l_miranda_energy": self.option("l_miranda_energy"),
            "l_miranda_strict": self.option("l_miranda_strict"),
            "l_rnahybird_num": self.option("l_rnahybird_num"),
            "l_rnahybird_energy": self.option("l_rnahybird_energy"),
            "l_rnahybird_pvalue": self.option("l_rnahybird_pvalue"),
            "l_ps_robot_score": self.option("l_ps_robot_score"),
            "l_targetfinder_score": self.option("l_targetfinder_score"),

            "c_method": ",".join(c_method),
            "c_min_support": self.option("c_min_support"),
            "c_miranda_score": self.option("c_miranda_score"),
            "c_miranda_energy": self.option("c_miranda_energy"),
            "c_miranda_strict": self.option("c_miranda_strict"),
            "c_rnahybird_num": self.option("c_rnahybird_num"),
            "c_rnahybird_energy": self.option("c_rnahybird_energy"),
            "c_rnahybird_pvalue": self.option("c_rnahybird_pvalue"),
            "c_ps_robot_score": self.option("c_ps_robot_score"),
            "c_targetfinder_score": self.option("c_targetfinder_score"),
        }

        self.option_dict = options

        self.target_predict.set_options(options)
        self.target_predict.on('end', self.set_output)
        self.target_predict.run()

    def run(self):
        self.get_run_log()
        self.run_target()
        super(TargetMiWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("whole_transcriptome", table="target_mi", main_id=self.option('target_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_output(self):
        self.logger.info("结果路径为{}".format(self.target_predict.output_dir))
        #output_dir = self.annotation.output_dir
        self.logger.info("结果路径为{}".format(self.target_predict.output_dir))
        self.set_db()

    def linkdir(self, olddir, newname, mode='link'):
        """
        移动目录下的输出文件/文件夹到输出文件夹下
        """
        if not os.path.isdir(olddir):
            raise Exception('需要移动到output目录的文件夹不存在。')
        newdir = os.path.join(self.output_dir, newname)
        if os.path.exists(newdir):
            shutil.rmtree(newdir)
        os.mkdir(newdir)
        allfiles = os.listdir(olddir)
        oldfiles = [os.path.join(olddir, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            else:
                new1 = os.path.join(newdir, os.path.basename(oldfiles[i]))
                os.system("mv {} {}".format(oldfiles[i], new1))



    def set_db(self):
        self.logger.info("保存结果到mongo")
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.test_api = self.api.api("whole_transcriptome.small_target")

        '''
        params_target = {
            "min_support": self.option("min_support"),
        }
        if "," in self.option("method"):
            for i in self.option("method").split(","):
                params_target.update({
                    i: "yes"
                })
        else:
            for i in self.option("method").split(";"):
                params_target.update({
                    i: "yes"
                })
        '''

        method_set = set()
        for method in [
            'm_miranda', 'm_targetscan', 'm_psrobot', 'm_targetfinder', 'm_rnahybrid',
            'l_miranda', 'l_targetscan', 'l_psrobot', 'l_targetfinder', 'l_rnahybrid',
            'c_miranda', 'c_targetscan', 'c_psrobot', 'c_targetfinder', 'c_rnahybrid']:
            if self.option(method) == "yes":
                method_set.add(method)

        params_target = dict()
        self.test_api.import_target_detail_web(
            self.option("target_id"),
            self.target_predict.output_dir,
            method_set,
            self.option("novol").prop['path'],
            self.option("known").prop['path'],
            anno_type="latest",
            species_name=self.option("species_name"),
            last_id_target=self.option("last_id_target"),
            target_dir = self.target_predict.output_dir,
            version=self.option("version")
        )
        self.set_file()
        self.end()

    def set_file(self):
        map_dict = {"target_dir": self.target_predict.output_dir}
        arg_dict = self.option_dict
        output_dir = self.output_dir
        mirna.set_mirna_target(map_dict, arg_dict, output_dir)

    def end(self):
        self.inter_dirs = [
            ["03 Target_Predict_Analysis", "", "靶基因预测分析结果目录",0],
            ["03 arget_Predict_Analysis/01 miRNA_Target", "", "miRNA靶基因预测", 0]
        ]
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        target_dir = self.output_dir

        repaths = [
            ['.', '', '靶基因预测文件', 0],
            ['tar_predict_detail.xls', 'xls', '靶基因预测详情表', 0],
            ['target_seqs', '', '靶基因序列', 0],
            ['target_align_detail', '', '靶基因比对详细信息结果目录', 0],
            ['target_seqs/mRNA_known_target.fa', 'fasta', '已知mRNA靶基因序列', 0],
            ['target_seqs/mRNA_novel_target.fa', 'fasta', '新mRNA靶基因序列', 0],
            ['target_seqs/lncRNA_known_target.fa', 'fasta', '已知lncRNA靶基因序列', 0],
            ['target_seqs/lncRNA_novel_target.fa', 'fasta', '新lncRNA靶基因序列', 0],
            ['target_seqs/circRNA_known_target.fa', 'fasta', '已知circRNA靶基因序列', 0],
            ['target_seqs/circRNA_novel_target.fa', 'fasta', '新circRNA靶基因序列', 0],
            ['target_align_detail/.*.txt.gz', 'gz', '靶基因比对详细信息', 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ]

        self.workflow_output_tmp = self._sheet.output
        '''
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        elif re.match(r'sanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
        elif re.match(r'^\w+://\S+/.+$', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp
        else:
            self.set_error("json output wrong")
        '''

        result_dir = self.add_upload_dir(target_dir)
        result_dir.add_regexp_rules(repaths)
        super(TargetMiWorkflow, self).end()
