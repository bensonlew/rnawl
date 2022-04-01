# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180422

import os
import re
import unittest

import web
import json
import datetime

from bson import ObjectId
from biocluster.config import Config

from mainapp.libs.signature import check_sig
from mainapp.controllers.project.lnc_rna_controller import LncRnaController


class DrawCircosAction(LncRnaController):
    """
    画circos的接口
    """

    def __init__(self):
        super(DrawCircosAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        """
            {"name": "target_trans", "type": "infile", "format": "lnc_rna.lnc_common"},
            {"name": "target_cis", "type": "infile", "format": "lnc_rna.lnc_common"},
            {"name": "diff_exp", "type": "infile", "format": "lnc_rna.lnc_common"},
            {"name": "rna_type", "type": "infile", "format": "lnc_rna.lnc_common"},
            {"name": "top_ref_num", "type": "int", "default": 10},
            # {"name": "diff_group", "type": "string"},

            {"name": "gtf", "type": "infile", "format": "lnc_rna.lnc_gtf"},
            {"name": "ref_fa_fai", "type": "infile", "format": "lnc_rna.lnc_common"},

            {"name": "main_id", "type": "string"},
            {"name": "project_type", "type": "string"}
        :return:
        """
        data = web.input()
        params = ["diff_exp", "diff_group", "top_ref_num",
                  "task_id", "task_type", "submit_location"]
        add_params = ["target_trans", "target_cis", "target_cistrans"]

        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C3100501", "variables": var}
                return json.dumps(info)

        diff_exp_obj = ObjectId(data['diff_exp'])
        diff_info = \
            self.lnc_rna.get_main_info(main_id=diff_exp_obj, collection_name='sg_diff', task_id=data['task_id'])
        exp_level = diff_info['exp_level']
        if exp_level != 'T':
            info = {"success": False, "info": "表达水平必须为转录本水平：T", "code": "C3100501", "variables": ["diff_exp"]}
            return json.dumps(info)


        if hasattr(data, "target_trans") and hasattr(data, "target_cis"):
            params_json = {
                "target_trans": data.target_trans,
                "target_cis": data.target_cis,
                "diff_exp": data.diff_exp,
                "diff_group": data.diff_group,
                "top_ref_num": int(data.top_ref_num),
                "task_id": data.task_id,
                "task_type": int(data.task_type),
                "submit_location": data.submit_location
            }
        elif hasattr(data, "target_cistrans"):
            params_json = {
                "target_cistrans": data.target_cistrans,
                "diff_exp": data.diff_exp,
                "diff_group": data.diff_group,
                "top_ref_num": int(data.top_ref_num),
                "task_id": data.task_id,
                "task_type": int(data.task_type),
                "submit_location": data.submit_location
            }

        task_info = self.lnc_rna.get_task_info(data.task_id)
        # noinspection PyBroadException
        try:
            project_sn = task_info["project_sn"]
            lncrna_gtf = task_info["all_lnc_gtf"]
            mrna_gtf = task_info["mrna_gtf"]
            ref_fa_fai = task_info["ref_fa_fai"]
        except:
            info = {"success": False, "info": "sg_task表里没有project_sn, member_id, species_version_id信息，请检查!",
                    "code": "C3100502", "variables": ""}
            return json.dumps(info)

        if os.path.exists(ref_fa_fai):
            pass
        else:
            genome_dir = Config().SOFTWARE_DIR + "/database/Genome_DB_finish/"
            ref_fa_fai = genome_dir + ref_fa_fai.split("/database/Genome_DB_finish/")[-1]

            if os.path.exists(ref_fa_fai):
                pass
            else:
                info = {"success": False, "info": "找不到基因组文件"}
        params_json = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        main_table_name = 'Circos_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params_json),
            ("desc", "Circos图主表！"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = self.lnc_rna.insert_main_table(collection_name="sg_circos", data=mongo_data)
        update_info = {str(main_id): "sg_circos"}

        # prepare to file
        to_files = ["lnc_geneset.export_diff_exp_matrix(diff_exp)",
                    "lnc_geneset.export_rna_type_matrix(rna_type)"]
        if hasattr(data, "target_trans") and hasattr(data, "target_cis"):
            to_files.extend(["lnc_geneset.export_trans_matrix(target_trans)",
                             "lnc_geneset.export_cis_matrix(target_cis)"])
        elif hasattr(data, "target_cistrans"):
            to_files.extend(["lnc_geneset.export_trans_matrix2(target_trans)",
                             "lnc_geneset.export_cis_matrix2(target_cis)"])

        options = {
            "diff_exp": {'diff_group': data.diff_group, 'sg_diff_id': data.diff_exp},
            "rna_type": {'task_id': data.task_id, 'exp_level': 'T', 'diff_exp': data.diff_exp},
            "top_ref_num": int(data.top_ref_num),
            "update_info": json.dumps(update_info),
            "lncrna_gtf": lncrna_gtf,
            "mrna_gtf": mrna_gtf,
            "ref_fa_fai": ref_fa_fai,
            "main_id": str(main_id),
            'task_id': data.task_id
        }
        if hasattr(data, "target_trans") and hasattr(data, "target_cis"):
            options.update({
                "target_trans": data.target_trans,
                "target_cis": data.target_cis,
            })
        elif hasattr(data, "target_cistrans"):
            options.update({
                "target_trans": data.target_cistrans,
                "target_cis": data.target_cistrans,

            })


        self.set_sheet_data(name="lnc_rna.report.draw_circos",
                            project_sn=project_sn,
                            task_id=data.task_id,
                            main_table_name="circos/" + main_table_name,
                            options=options,
                            to_file=to_files,
                            params=params_json)

        task_info = super(DrawCircosAction, self).POST()
        task_info['id'] = str(main_id)

        return json.dumps(task_info)


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """

        def test_this(self):
            cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
            cmd += 'post '
            cmd += "-fr no "
            cmd += '-c {} '.format("client03")
            cmd += "s/lnc_rna/draw_circos "
            cmd += "-b http://bcl.tsg.com "
            '''
            args = dict(
                task_id="tsg_34051",
                # Conf配置
                task_type="2",
                # Conf配置
                submit_location="drawcircos",
                # "target_trans": trans作用靶基因预测结果表 主表ID
                target_trans='5cd01dd817b2bf25d2010234',
                # "target_cis": cis作用靶基因预测结果表 主表ID
                target_cis='5cd94be617b2bf0f9c3ec5a9',
                # "diff_exp": 表达量差异分析结果表 主表ID
                diff_exp='5cd3b5f317b2bf18dcbef37e',
                # "diff_group": 选择差异组别, 示例： Con|Vit
                diff_group='S_16C_16C_1|S_16C_25C_1',
                # "top_ref_num": 取染色体长度最长前 多少 个,
                top_ref_num="20",
            )
            '''
            args = dict(
                task_id="tsg_35612",
                # Conf配置
                task_type="2",
                # Conf配置
                submit_location="drawcircos",
                # "target_trans": trans作用靶基因预测结果表 主表ID
                target_cistrans='5d9c5a8817b2bf5c984873e3',
                diff_exp="5d9c57e817b2bf5c9843869f",
                # "diff_group": 选择差异组别, 示例： Con|Vit
                diff_group='S_16C_16C_1|S_16C_25C_1',
                # "top_ref_num": 取染色体长度最长前 多少 个,
                top_ref_num="20",
            )
            arg_names, arg_values = args.keys(), args.values()
            cmd += '-n "{}" -d "{}" '.format(";".join([str(x) for x in arg_names]),
                                             ";".join([str(x) for x in arg_values]))

            print(cmd)
            os.system(cmd)


    unittest.main()
