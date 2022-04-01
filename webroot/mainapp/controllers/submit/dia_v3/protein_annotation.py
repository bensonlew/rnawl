# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import web
import json
import datetime
import unittest
import os
from mainapp.libs.signature import check_sig
from mainapp.libs.param_pack import *
from mainapp.controllers.project.dia_controller import DiaController
from mbio.api.to_file.dia import *


class ProteinAnnotationAction(DiaController):
    """
    蛋白注释重运行接口
    """
    def __init__(self):
        super(ProteinAnnotationAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get("HTTP_CLIENT")
        print data
        return_result = self.check_options(data)

        if return_result:
            info = {"success": False, "info": return_result}
            return json.dumps(info)
        stat_info = self.dia.get_main_info_by_record("sg_annotation_stat", task_id=str(data.task_id), type="origin")
        if not stat_info:
            info = {"success": False, "info": "stat_id不存在,请确认参数是否正确"}
            return json.dumps(info)
        print stat_info

        run_info = self.dia.get_main_info_by_record("sg_annotation_stat", task_id=str(data.task_id), status="start", type="latest")
        if run_info:
            info = {"success": False, "info": "任务{}正在计算，请等待其完成在下次运行".format(run_info['name'])}
            return json.dumps(info)

        latest_info = self.dia.get_main_info_by_record("sg_annotation_stat", task_id=str(data.task_id), status="end",  type="latest")
        latest_main_id = ""

        task_info = self.dia.get_task_info(data.task_id)

        if latest_info:
            latest_main_id = latest_info["_id"]

        taxon = "Animals"
        if stat_info.has_key("taxonomy"):
            taxon = stat_info["taxonomy"]
        else:
            pass

        params_json = json.loads(stat_info["params"])

        params_json.update({
            "task_id": str(data.task_id),
            "go_evalue": str(data.go_evalue),
            "go_identity": str(data.go_identity),
            "cog_evalue": str(data.cog_evalue),
            "cog_identity": str(data.cog_identity),
            "kegg_evalue": str(data.kegg_evalue),
            "kegg_identity": str(data.kegg_identity),
            "pfam_evalue": str(data.pfam_evalue),
            "submit_location": data.submit_location,
            "task_type": int(data.task_type)
        })
        main_table_name = "AnnoStat_latest_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        mongo_data = [
            ("project_sn", stat_info["project_sn"]),
            ("task_id", stat_info["task_id"]),
            ("status", "start"),
            ("name", main_table_name),
            ("database", "go,cog,kegg,pfam,subloc"),
            ("desc", "注释统计主表"),
            ("result_dir", ""),
            ("taxonomy", taxon),
            ("type", "latest"),
            ('version', 'v3'),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.dia.insert_main_table("sg_annotation_stat", mongo_data)
        update_info = {str(main_table_id): "sg_annotation_stat"}
        origin_result = stat_info["result_dir"]
        if not origin_result.endswith("/"):
            origin_result = origin_result + "/"
        origin_param = stat_info["params"]

        if "annot_group" in task_info:
            annot_group = task_info["annot_group"]
        else:
            annot_group = "GROUP_2017"
        options = {
            "go_evalue": float(data.go_evalue),
            "go_identity": float(data.go_identity),
            "cog_evalue": float(data.cog_evalue),
            "cog_identity": float(data.cog_identity),
            "kegg_evalue": float(data.kegg_evalue),
            "kegg_identity": float(data.kegg_identity),
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
        self.set_sheet_data(name="dia_v3.report.protein_annotation",
                            options=options,
                            #main_table_name="AnnotationStat/" + main_table_name,
                            main_table_name=main_table_name,
                            task_id=stat_info["task_id"],
                            project_sn=stat_info["project_sn"],
                            to_file=to_file, )
        task_info = super(ProteinAnnotationAction, self).POST()
        task_info["content"] = {"ids": {"id": str(main_table_id), "name": main_table_name}}
        ##
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.dia.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.dia.update_group_compare_is_use(data.task_id, data.control_id)
        ##
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

        if not (hasattr(data, "go_evalue")):
            data.go_evalue = 1e-3
        else:
            try:
                if float(data.go_evalue) > 1e-3:
                    warn_info.append("GO E-value值必须小于1e-3")
            except:
                warn_info.append('输入的GO E-value不是数字 ')

        if not (hasattr(data, "go_identity")):
            data.go_identity = 0
        else:
            try:
                if float(data.go_identity) < 0 or float(data.go_identity) > 100:
                    warn_info.append("GO Identity值需在0-100范围内")
            except:
                warn_info.append('GO Identity值必须是数字')


        if not (hasattr(data, "cog_evalue")):
            data.cog_evalue = 1e-3
        else:
            try:
                if float(data.cog_evalue) > 1e-3:
                    warn_info.append("COG E-value值必须小于1e-3")
            except:
                warn_info.append('输入的COG E-value不是数字 ')


        if not (hasattr(data, "cog_identity")):
            data.cog_identity = 0
        else:
            try:
                if float(data.cog_identity) < 0 or float(data.cog_identity) > 100:
                    warn_info.append("COG Identity值需在0-100范围内")
            except:
                warn_info.append('COG Identity值必须是数字')

        if not (hasattr(data, "kegg_evalue")):
            data.kegg_evalue = 1e-3
        else:
            try:
                if float(data.kegg_evalue) > 1e-3:
                    warn_info.append("KEGG E-value值必须小于1e-3")
            except:
                warn_info.append('输入的KEGG E-value不是数字 ')

        if not (hasattr(data, "kegg_identity")):
            data.kegg_identity = 0
        else:
            try:
                if float(data.kegg_identity) < 0 or float(data.kegg_identity) > 100:
                    warn_info.append("KEGG Identity值需在0-100范围内")
            except:
                warn_info.append('KEGG Identity值必须是数字')

        if not (hasattr(data, "pfam_evalue")):
            data.pfam_evalue = 1e-3
        else:
            try:
                if float(data.pfam_evalue) > 1e-3:
                    warn_info.append("PFAM E-value值必须小于1e-3")
            except:
                warn_info.append('输入的PFAM E-value不是数字 ')


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
        cmd += "s/dia_v3/protein_annotation "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="dia_test",
            go_evalue=0.00001,
            go_identity=80,
            cog_evalue=0.00001,
            cog_identity=80,
            kegg_evalue=1e-6,
            kegg_identity=90,
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
