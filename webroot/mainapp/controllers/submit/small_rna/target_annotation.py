# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import web
import json
import datetime
import unittest
import os
from mainapp.libs.signature import check_sig
from mainapp.libs.param_pack import *
from mainapp.controllers.project.small_rna_controller import SmallRnaController
from biocluster.config import Config


class TargetAnnotationAction(SmallRnaController):
    """
    small_rna 靶基因预测与注释重运行接口
    """
    def __init__(self):
        super(TargetAnnotationAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get("HTTP_CLIENT")
        print(data)
        return_result = self.check_options(data)

        if return_result:
            variables = []
            variables.append(return_result)
            info = {"success": False, "info": '%s' % return_result}
            return json.dumps(info)
        stat_info = self.small_rna.get_main_info_by_record("sg_annotation_stat", task_id=str(data.task_id), type="origin")
        if not stat_info:
            info = {"success": False, "info": "stat_id不存在,请确认参数是否正确"}
            return json.dumps(info)
        print(stat_info)

        task_info = self.small_rna.get_main_info_by_record("sg_task", task_id=str(data.task_id))
        if not stat_info:
            info = {"success": False, "info": "task_id不存在记录,请确认参数是否正确"}

        if "version" in task_info:
            version = task_info["version"]
        else:
            version = "v1.0"

        run_info = self.small_rna.get_main_info_by_record("sg_target", task_id=str(data.task_id), status="start", type="latest")
        if run_info:
            variables = []
            variables.append(run_info['name'])
            info = {"success": False, "info": "任务%s正在计算，请等待其完成在下次运行" % run_info['name']}
            return json.dumps(info)


        methods = []
        for method in ['miranda', 'targetscan', 'psrobot', 'targetfinder', 'rnahybrid']:
            if hasattr(data, method) and data[method] == "yes":
                methods.append(method)

        if len(methods) == 0:
            info = {"success": False, "info": "请至少选择一款软件"}
            return json.dumps(info)
        elif len(methods) < int(data.min_support):
            info = {"success": False, "info": "最少支持的数量多于选择的预测方法"}
            return json.dumps(info)

        ## 检查参数阈值范围
        for par in ['miranda_score', 'miranda_energy',
                    'rnahybird_num', 'rnahybird_energy', 'rnahybird_pvalue',
                    'ps_robot_score', 'targetfinder_score']:
            if hasattr(data, par):
                try:
                    par_value = float(getattr(data, par))
                except:
                    info = {"success": False, "info": "{} 参数设置错误".format(par)}
                    return json.dumps(info)
                '''
                if par == 'miranda_score':
                    if data.miranda_score
                if par == 'miranda_energy':

                if par == 'rnahybird_num':
                if par == 'rnahybird_energy':
                if par == 'rnahybird_pvalue':
                if par == 'ps_robot_score':
                if par == 'targetfinder_score':
                '''

        latest_info = self.small_rna.get_main_info_by_record("sg_annotation_stat", task_id=str(data.task_id), status="end",  type="latest")
        latest_main_id = ""

        spe_tax = {
                "Mus_musculus": 10090,
                "Rattus_norvegicus": 10116,
                "Monodelphis_domestica": 13616,
                "Xenopus_tropicalis": 8364,
                "Gallus_gallus": 9031,
                "Macaca_mulatta": 9544,
                "Pan_troglodytes": 9598,
                "Homo_sapiens":	9606,
                "Canis_lupus_familiaris": 9615,
                "Bos_taurus": 9913
            }

        if not task_info['organism_name'] in spe_tax and hasattr(data, 'targetscan') and data['targetscan'] == "yes":
            info = {"success": False, "info": "该项目物种 %s 不支持targetscan" % stat_info["species_name"]}
            return json.dumps(info)

        task_info = self.small_rna.get_task_info(data.task_id)
        if latest_info:
            latest_main_id = latest_info["_id"]

        latest_info_target = self.small_rna.get_main_info_by_record("sg_target", task_id=str(data.task_id), status="end",  type="latest")
        latest_main_id_target = ""

        target_info = self.small_rna.get_main_info_by_record("sg_target", task_id=str(data.task_id), status="end",  type="origin")
        result = target_info["result_dir"]

        task_info = self.small_rna.get_task_info(data.task_id)

        if 'target' in task_info and os.path.exists(task_info["target"]):
            target = task_info["target"]
        else:
            target = result + "/" + "target.fa"
        if latest_info_target:
            latest_main_id_target = latest_info_target["_id"]

        params_json = {
            "task_id": str(data.task_id),
            "min_support": str(data.min_support),
            "submit_location": data.submit_location,
            "task_type": int(data.task_type)
        }
        for method in ['miranda', 'targetscan', 'psrobot', 'targetfinder', 'rnahybrid']:
            if hasattr(data, method):
                params_json.update({
                    method: data[method]
                })

        for par in ['miranda_score', 'miranda_energy', 'miranda_strict',
                    'rnahybird_num', 'rnahybird_energy', 'rnahybird_pvalue',
                    'ps_robot_score', 'targetfinder_score']:
            if hasattr(data, par):
                params_json.update({
                    par: data[par]
                })

        main_table_name = "AnnotationStat_latest_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        mongo_data = [
            ("project_sn", stat_info["project_sn"]),
            ("task_id", stat_info["task_id"]),
            ("status", "start"),
            ("name", main_table_name),
            ("database", "nr,swissprot,cog,kegg,pfam"),
            ("desc", "注释统计主表"),
            ("result_dir", ""),
            ("taxonomy", stat_info["taxonomy"]),
            ("type", "latest"),
            ("species_name", stat_info["species_name"]),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("result_dir", stat_info["result_dir"]),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.small_rna.insert_main_table("sg_annotation_stat", mongo_data)

        main_table_name2 = "Target_latest_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        mongo_data2 = [
            ("project_sn", stat_info["project_sn"]),
            ("task_id", stat_info["task_id"]),
            ("status", "start"),
            ("name", main_table_name2),
            ("desc", "靶基因预测主表"),
            ("result_dir", ""),
            ("taxonomy", stat_info["taxonomy"]),
            ("version", version),
            ("type", "latest"),
            ("species_name", stat_info["species_name"]),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]

        main_table_id2 = self.small_rna.insert_main_table("sg_target", mongo_data2)

        update_info = {str(main_table_id2): "sg_target"}
        # origin_param = stat_info["params"]

        if os.path.exists(stat_info["result_dir"]):
            pass
        else:
            genome_dir = Config().SOFTWARE_DIR + "/database/Genome_DB_finish/"
            stat_info["result_dir"] = genome_dir + stat_info["result_dir"].split("/database/Genome_DB_finish/")[-1]
            if os.path.exists(stat_info["result_dir"]):
                pass
            else:
                info = {"success": False, "info": "找不到基因组文件"}

        options = {
            "task_id": str(data.task_id),
            "novol": task_info["novol"],
            "known": task_info["known"],
            "method": ",".join(methods),
            "anno_path": stat_info["result_dir"],
            "stat_id": str(main_table_id),
            "last_id": str(latest_main_id),
            "last_id_target": str(latest_main_id_target),
            "target_id": str(main_table_id2),
            "species_name": task_info["organism_name"],
            "taxonomy": stat_info["taxonomy"],
            "update_info": json.dumps(update_info),
            "min_support": data.min_support,
            "version": version,
            "ref":target,
        }

        for par in ['miranda_score', 'miranda_energy', 'miranda_strict',
            'rnahybird_num', 'rnahybird_energy', 'rnahybird_pvalue',
            'ps_robot_score', 'targetfinder_score']:
            if hasattr(data, par):
                options.update({
                    par: data[par]
                })


        to_file = []
        self.set_sheet_data(name="small_rna.report.target_annotation",
                            options=options,
                            #main_table_name="AnnotationStat/" + main_table_name,
                            main_table_name=main_table_name2,
                            task_id=stat_info["task_id"],
                            project_sn=stat_info["project_sn"],
                            to_file=to_file, )
        task_info = super(TargetAnnotationAction, self).POST()
        task_info["content"] = {"ids": {"id": str(main_table_id), "name": main_table_name2}}
        return json.dumps(task_info)

    def check_options(self, data):
        """
        检查网页端传来的参数是否正确
        """
        warn_info = list()
        method_num = 0
        for method in ['miranda', 'targetscan', 'psrobot', 'targetfinder', 'rnahybrid']:
            if hasattr(data, method):
                print "method is" + method
                method_num += 1

        if not (hasattr(data, "min_support")):
            warn_info.append("缺少参数min_support")
        elif int(data.min_support) > method_num:
            warn_info.append("最少支持的数量多于选择的预测方法")
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
        cmd += "s/small_rna/target_annotation "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict({
            "task_id": "small_rna",
            "rnahybrid": "yes",
            "miranda": "yes",
            "targetscan": "yes",
            "min_support": 2,
            "submit_location": "target_annotation",
            "task_type":2,
        })
        arg_names, arg_values = args.keys(), args.values()
        arg_values = [str(x) for x in arg_values]
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
