# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import os
import time
import sqlite3
from bson.objectid import ObjectId
from biocluster.config import Config
import datetime
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
import re
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder

class QuerySeqMultiWorkflow(Workflow):
    """
    表达量分布
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(QuerySeqMultiWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='seq_id', type='string'),
            dict(name="seq_type", type="string"),
            dict(name="main_id", type="string"),
            dict(name="task_id", type="string"),
            # to update sg_status
            {"name": "seq_db", "type": "infile", "format": "denovo_rna_v2.common"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/05 Structure_Analysis/01 CDS_predict')
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
        super(QuerySeqMultiWorkflow, self).send_log(data)

    def run(self):
        self.start_listener()
        self.fire("start")
        self.get_run_log()
        name = self.choose_seq()
        self.set_db(name)

    def get_run_log(self):
        get_run_log = GetRunLog("denovo_rna_v2", table="sg_query_seq_multi", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def choose_seq(self):
        conn = sqlite3.connect(self.option("seq_db").prop['path'])
        cursor = conn.cursor()
        if self.option("seq_type") == 'T':
            query_field = 't_id'
        else:
            query_field = 'g_id'

        time_now = datetime.datetime.now()
        time_str = time_now.strftime("%Y%m%d_%H%M%S")
        name = "cds_sequence_{}.fa".format(time_str)

        with open(self.output_dir + "/" + name, 'w') as f:
            for seq_id in set(self.option("seq_id").split(",")):
                cursor.execute("SELECT * FROM {} WHERE {}='{}'".format('seq_cds', query_field, seq_id))
                cds_info = cursor.fetchall()
                cds_data = list()
                if cds_info:
                    for each in cds_info:
                        tmp = dict(zip(['t_id', 'g_id', 'cds_pos', 'cds', 'pep', 'type'], each))
                        tmp['cds_length'] = len(tmp['cds'])
                        tmp['pep_length'] = len(tmp['pep'])
                        cds_data.append(tmp)
                else:
                    print("{} not found".format(seq_id))
                    tmp = dict(zip(['t_id', 'g_id', 'cds_pos', 'cds', 'pep', 'type'],
                                   [seq_id, seq_id, 'None', 'None', 'None', 'None']))
                    tmp['cds_length'] = 0
                    tmp['pep_length'] = 0
                    cds_data.append(tmp)
                #
                cursor.execute("SELECT sequence FROM {} WHERE {}='{}'".format('transcript_gene', query_field, seq_id))
                seq = cursor.fetchall()
                if seq:
                    sequence = seq[0][0]
                else:
                    sequence = 'None'

                f.write(">{}\n{}\n".format(seq_id, sequence))
                for cds in cds_data:
                    f.write(">{}\n{}\n".format(cds["t_id"]+ " " + cds["cds_pos"] + str(cds["cds_length"]), cds["cds"]))
                    f.write(">{}\n{}\n".format(cds["t_id"]+ " " + cds["cds_pos"] + str(cds["pep_length"]), cds["pep"]))
                f.write("\n")
        cursor.close()
        return name


    def set_db(self, name):
        """
        保存结果表到mongo数据库中

        """
        output_dir = self._sheet.output
        file_dir = output_dir + "/" + name
        # file_dir = self.output_dir + "/choose_seq.fa"
        api_base = self.api.api("denovo_rna_v2.api_base")
        api_base.update_db_record_2('sg_query_seq_multi', self.option("main_id"), status="end", main_id=ObjectId(self.option("main_id")), seq_dir=file_dir)

        db = Config().get_mongo_client(mtype="project")[Config().get_mongo_dbname("project")]
        time_now = datetime.datetime.now()

        download_col = db['sg_pipline_download_list']
        download_col.insert({
            "controller" : "denovo_rna_v2",
            "action" : "index",
            "collect" : "sg_cds_detail",
            "task_id" : self.option("task_id") ,
            "submit_location" : "sg_cds_detail",
            "field" : {
            },
            "condition" : "",
            "status" : 2,
            "created_ts": time_now.strftime('%Y-%m-%d %H:%M:%S'),
            "name": "sg_cds_detail_" + time_now.strftime('%Y-%m-%d %H:%M:%S').replace(" ","_").replace("-","").replace(":",""),
            "url": file_dir,
        })



        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["05 Structure_Analysis", "", "SNP分析结果目录",0],
            ["05 Structure_Analysis/01 CDS_predict", "", "CDS预测", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "序列下载文件",0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        result_dir.add_regexp_rules([
            [r'cds_sequence.*fa', 'fa', 'CDS序列',0]
        ])
        super(QuerySeqMultiWorkflow, self).end()
