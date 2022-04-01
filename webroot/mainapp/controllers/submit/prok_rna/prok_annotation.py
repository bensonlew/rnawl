# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import web
import json
import datetime
import unittest
import os
from mainapp.libs.signature import check_sig
from mainapp.libs.param_pack import *
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mbio.api.to_file.prok_rna import *


class ProkAnnotationAction(ProkRNAController):
    """
    蛋白注释重运行接口
    """
    def __init__(self):
        super(ProkAnnotationAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get("HTTP_CLIENT")
        print data
        return_result = self.check_options(data)

        if return_result:
            info = {"success": False, "info": return_result}
            return json.dumps(info)
        stat_info = self.prok_rna.get_main_info_by_record("sg_annotation_stat", task_id=str(data.task_id), type="origin")
        if not stat_info:
            info = {"success": False, "info": "stat_id不存在,请确认参数是否正确"}
            return json.dumps(info)
        print stat_info

        run_info = self.prok_rna.get_main_info_by_record("sg_annotation_stat", task_id=str(data.task_id), status="start", type="latest")
        if run_info:
            info = {"success": False, "info": "任务{}正在计算，请等待其完成后再运行".format(run_info['name'])}
            return json.dumps(info)

        latest_info = self.prok_rna.get_main_info_by_record("sg_annotation_stat", task_id=str(data.task_id), status="end",  type="latest")
        latest_main_id = ""

        task_info = self.prok_rna.get_task_info(data.task_id)

        if latest_info:
            latest_main_id = latest_info["_id"]

        params_json = json.loads(stat_info["params"])

        params_json.update({
            "task_id": str(data.task_id),
            "nr_evalue": str(data.nr_evalue),
            "cog_evalue": str(data.cog_evalue),
            "kegg_evalue": str(data.kegg_evalue),
            "pfam_evalue": str(data.pfam_evalue),
            "swissprot_evalue": str(data.swissprot_evalue),
            "submit_location": data.submit_location,
            "task_type": int(data.task_type)
        })
        main_table_name = "AnnotationStat_latest_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        mongo_data = [
            ("project_sn", stat_info["project_sn"]),
            ("task_id", stat_info["task_id"]),
            ("status", "start"),
            ("name", main_table_name),
            ("database", "nr,swissprot,go,cog,kegg,pfam"),
            ("desc", "注释统计主表"),
            ("result_dir", ""),
            ("type", "latest"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.prok_rna.insert_main_table("sg_annotation_stat", mongo_data)
        update_info = {str(main_table_id): "sg_annotation_stat"}
        origin_result = stat_info["result_dir"]
        if not origin_result.endswith('/'):
            origin_result += '/'
        origin_param = stat_info["params"]
        if "annot_group" in task_info:
            annot_group = task_info["annot_group"]
        else:
            annot_group = "REFRNA_GROUP_2019"
        options = {
            "nr_evalue": float(data.nr_evalue),
            "swissprot_evalue": float(data.swissprot_evalue),
            "cog_evalue": float(data.cog_evalue),
            "kegg_evalue": float(data.kegg_evalue),
            "pfam_evalue": float(data.pfam_evalue),
            "task_id": str(data.task_id),
            "stat_id": str(main_table_id),
            "last_id": str(latest_main_id),
            "origin_result": str(origin_result),
            "origin_param": origin_param,
            "annot_group": annot_group,
            "update_info": json.dumps(update_info)
        }
        to_file = []
        self.set_sheet_data(name="prok_rna.report.prok_annotation",
                            options=options,
                            #main_table_name="AnnotationStat/" + main_table_name,
                            main_table_name=main_table_name,
                            task_id=stat_info["task_id"],
                            project_sn=stat_info["project_sn"],
                            to_file=to_file, )
        task_info = super(ProkAnnotationAction, self).POST()
        task_info["content"] = {"ids": {"id": str(main_table_id), "name": main_table_name}}
        return json.dumps(task_info)

    def check_options(self, data):
        """
        检查网页端传来的参数是否正确
         ["task_id", "nr_evalue", "nr_score", "nr_similarity", "nr_identity",
          "swissprot_evalue", "swissprot_score", "swissprot_similarity",
           "swissprot_identity", "submit_location"]
        """
        warn_info = list()
        if not (hasattr(data, "task_id")):
            warn_info.append("缺少参数task_id")

        if not (hasattr(data, "task_type")):
            warn_info.append("缺少参数task_type")

        if not (hasattr(data, "nr_evalue")):
            data.nr_evalue = 1e-5
        else:
            try:
                if float(data.nr_evalue) <0:
                    warn_info.append("GO E-value值必须大于0")
            except:
                warn_info.append('输入的GO E-value不是数字 ')

        if not (hasattr(data, "swissprot_evalue")):
            data.swissprot_evalue = 1e-5
        else:
            try:
                if float(data.nr_evalue) <= 0:
                    warn_info.append("swissprot E-value值必需大于0")
            except:
                warn_info.append('输入的swissprot E-value不是数字 ')

        if not (hasattr(data, "cog_evalue")):
            data.cog_evalue = 1e-5
        else:
            try:
                if float(data.cog_evalue) < 0:
                    warn_info.append("COG E-value值必须大于0")
            except:
                warn_info.append('输入的COG E-value不是数字 ')

        if not (hasattr(data, "kegg_evalue")):
            data.kegg_evalue = 1e-5
        else:
            try:
                if float(data.kegg_evalue) < 0:
                    warn_info.append("KEGG E-value值必须大于0")
            except:
                warn_info.append('输入的KEGG E-value不是数字 ')

        return ";".join(warn_info)

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/prok_rna/prok_annotation "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="prok_rna_srna",
            nr_evalue=0.00001,
            cog_evalue=0.00001,
            kegg_evalue=1e-6,
            pfam_evalue=1e-4,
            submit_location="annotationstat",
            task_type=2,
        )
        arg_names, arg_values = args.keys(), args.values()
        arg_values = [str(x) for x in arg_values]
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
