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
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class TargetAnnotationWorkflow(Workflow):
    """
    交互分析进行筛选blast参数重注释时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        #self.rpc = False
        super(TargetAnnotationWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "novol", "type": "infile", "format": "small_rna.fasta"},  # 输入文件
            {"name": "known", "type": "infile", "format": "small_rna.fasta"},  # 输入文件
            {"name": "ref", "type": "infile", "format": "small_rna.fasta"},  # 输入文件
            {"name": "method", "type": "string", "default": "miranda"},
            {"name": "type", "type": "string", "default": "animal"},
            {"name": "min_support", "type": "int", "default": 1},
            {"name": "anno_detail", "type": "string", "default": None},
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
            {"name": "taxonomy", "type": "string", "default": None},
            {"name": "anno_path", "type": "string", "default": ""},
            {"name": "version", "type": "string", "default": "v1"},
            {"name": "miranda_score", "type": "string", "default": "160"},
            {"name": "miranda_energy", "type": "string", "default": "-20"},
            {"name": "miranda_strict", "type": "string", "default": "on"},
            {"name": "rnahybird_num", "type": "string", "default": "100"},
            {"name": "rnahybird_energy", "type": "string", "default": "-20"},
            {"name": "rnahybird_pvalue", "type": "string", "default": "0.01"},
            {"name": "ps_robot_score", "type": "string", "default": "2.5"},
            {"name": "targetfinder_score", "type": "string", "default": "4"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.target_predict = self.add_module("small_rna.target_predict")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/03 miRNA_Target')
        self.inter_dirs = []
        self.target_predict.on('end', self.set_db)

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
        super(TargetAnnotationWorkflow, self).send_log(data)

    def run_target(self):
        if os.path.exists(os.path.join(self.option("anno_path"), 'anno_stat/all_anno_detail.xls')):
            anno_detail = os.path.join(self.option("anno_path"), 'anno_stat/all_anno_detail.xls')
        else:
            anno_detail = os.path.join(self.option("anno_path"), 'anno_stat/trans_anno_detail.xls')
        options = {
            "novol": self.option("novol").prop['path'],
            "known": self.option("known").prop['path'],
            "ref": self.option("ref").prop['path'],
            "method": self.option("method"),
            "anno_detail": anno_detail,
            "version": self.option("version"),
            "species": self.option("species_name"),
            "type": self.option("taxonomy").lower(),
            "min_support": self.option("min_support"),
            "miranda_score": self.option("miranda_score"),
            "miranda_energy": self.option("miranda_energy"),
            "miranda_strict": self.option("miranda_strict"),
            "rnahybird_num": self.option("rnahybird_num"),
            "rnahybird_energy": self.option("rnahybird_energy"),
            "rnahybird_pvalue": self.option("rnahybird_pvalue"),
            "ps_robot_score": self.option("ps_robot_score"),
            "targetfinder_score": self.option("targetfinder_score"),
        }


        self.target_predict.set_options(options)
        self.target_predict.run()

    def run(self):
        self.logger.info(self.option('anno_path'))
        if not os.path.exists(self.option('anno_path')):
            anno_path = os.path.join(Config().SOFTWARE_DIR, self.option('anno_path').split("/app/")[1])
            self.option('anno_path', anno_path)
        self.logger.info(self.option('anno_path'))
        self.get_run_log()
        self.run_target()
        super(TargetAnnotationWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("small_rna", table="sg_target", main_id=self.option('target_id'),
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
        self.test_api = self.api.api("small_rna.target_annotation")

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


        self.test_api.import_target_detail_web(
            self.option("target_id"),
            self.target_predict.output_dir + '/novol_target.xls',
            self.target_predict.output_dir + '/known_target.xls',
            params_target,
            self.option("novol").prop['path'],
            self.option("known").prop['path'],
            anno_type="latest",
            species_name=self.option("species_name"),
            last_id_target=self.option("last_id_target"),
            target_dir = self.target_predict.output_dir,
            version=self.option("version")
        )

        params = {
            "nr_evalue": 1e-5,
            "nr_similarity": 0,
            "nr_identity": 0,
            "swissprot_evalue": 1e-5,
            "swissprot_similarity": 0,
            "swissprot_identity": 0,
            "cog_evalue": 1e-5,
            "cog_similarity": 0,
            "cog_identity": 0,
            "kegg_evalue": 1e-5,
            "kegg_similarity": 0,
            "kegg_identity": 0,
            "pfam_evalue": 1e-5,
        }
        self.test_api.anno_type = 'latest'
        self.test_api.run_web(
            self.target_predict.output_dir + '/known_target.xls',
            self.target_predict.output_dir + '/novol_target.xls',
            self.option("anno_path"),
            os.path.join(self.option("anno_path"), 'tran2gene.txt'),
            # os.path.dirname(self.option("anno_path").rstrip("/")) + '/g2t2p',
            params,
            self.option("task_id"),
            self.option("stat_id"),
            self.option("last_id"),
            version=self.option('version'),
        )
        self.end()

    '''
    def update_task_id(self, stat_id, blast_id, nr_id, swissprot_id):
        """更新主表task_id"""
        self.logger.info("更新主表task_id")
        db = Config().get_mongo_client(mtype="small_rna")[Config().get_mongo_dbname("small_rna")]
        #client = Config().mongo_client
        #db_name = Config().MONGODB + '_ref_rna'
        stat_coll = db['sg_annotation_stat']
        results = stat_coll.find_one({'_id': ObjectId(stat_id)})
        task_id = results['task_id']
        blast_coll = db['sg_annotation_blast']
        nr_coll = db['sg_annotation_nr']
        sw_coll = db['sg_annotation_swissprot']
        blast_coll.update({'_id': ObjectId(blast_id)}, {'$set': {'task_id': task_id}})
        nr_coll.update({'_id': ObjectId(nr_id)}, {'$set': {'task_id': task_id}})
        sw_coll.update({'_id': ObjectId(swissprot_id)}, {'$set': {'task_id': task_id}})
        self.logger.info("更新主表ID成功")
    '''

    def end(self):
        origin_dir = self.option("anno_path")
        for source in glob.glob(os.path.join(self.target_predict.output_dir, '*')):
            basename = os.path.basename(source)
            link_name = os.path.join(self.output_dir, basename)
            if os.path.exists(link_name):
                os.remove(link_name)
            os.link(source, link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))

        target_dir = self.output_dir
        trans_anno_detail = os.path.join(origin_dir + "/anno_stat/trans_anno_detail.xls")


        repaths = [
            [".", "", "miRNA靶基因预测文件"],
            ["All_annot_target.xls", "", "已知+新miRNA 对应的靶基因详情表"],
            ["known_target.xls", "", "已知miRNA 对应的靶基因详情表"],
            ["target.fa", "", "靶基因序列"],
            # [r".*detail.txt.gz", "", "靶基因比对详细信息"],
            ["novol_target.xls", "", "新miRNA 对应的靶基因详情表"],
            [r"known_*_detail.txt.gz", "", "已知miRNA靶基因比对详细信息"],
            [r"novol_*_detail.txt.gz", "", "新miRNA靶基因比对详细信息"],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ]

        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        elif re.match(r'sanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
        elif re.match(r'^\w+://\S+/.+$', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp
        else:
            self.set_error("json output wrong")

        if os.path.exists(os.path.join(target_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(target_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(target_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(target_dir)
        self.inter_dirs = [
            ["03 miRNA_Target", "", "miRNA靶基因预测结果目录",0],
        ]
        result_dir.add_regexp_rules(repaths)
        db = Config().get_mongo_client(mtype="small_rna")[Config().get_mongo_dbname("small_rna")]
        # col1 = db["sg_annotation_stat"]
        # col1.update({"_id": ObjectId(self.option("stat_id"))}, {"$set": {"result_dir": self.workflow_output + "/Annotation"}}, upsert=True)
        col2 = db["sg_target"]
        col2.update({"_id": ObjectId(self.option("target_id"))}, {"$set": {"result_dir": self.workflow_output}}, upsert=True)

        super(TargetAnnotationWorkflow, self).end()
