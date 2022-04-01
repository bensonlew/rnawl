# -*- coding: utf-8 -*-
# __author__ = "qinjincheng"

import web
import json
import datetime
from mainapp.controllers.project.small_rna_controller import SmallRnaController
from mbio.api.to_file.small_rna import *
from mainapp.libs.signature import check_sig
import unittest
import os

# 临时
# from webroot.mainapp.controllers.project.small_rna_controller import SmallRnaController

class GenesetVennAction(SmallRnaController):
    def __init__(self):
        super(GenesetVennAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()

        # 参数检测
        basic_args = ["task_id", "submit_location", "task_type", "geneset_ids", "gene_type"]
        # check whether arguments exist
        for arg in basic_args:
            if not hasattr(data, arg):
                return_info = {"success": False, "info": "Lack argument: {}".format(arg)}
                return json.dumps(return_info)

        # 主表创建相关信息
        tmp_time = datetime.datetime.now()

        created_ts = tmp_time.strftime("%Y-%m-%d %H:%M:%S")
        str_time = tmp_time.strftime("%Y%m%d_%H%M%S")

        desc = "Geneset_venn_main_table_{}".format(str_time)
        name = "Geneset_venn_{}_{}".format(data.gene_type, str_time)

        # 通过主表sg_task获取项目等信息
        task_id = data.task_id
        project_sn = self.small_rna.get_task_info(task_id=task_id)["project_sn"]

        # 创建参数json文档，用于检测是否重复投递任务
        geneset_id_list = str(data.geneset_ids).split(",")
        params = {
            "task_id": task_id,
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "gene_type": data.gene_type,
            "geneset_ids": geneset_id_list
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))

        # 创建venn任务表（用于前端查阅任务状态）
        main_info = {
            "task_id": task_id,
            "project_sn": project_sn,
            "name": name,
            "created_ts": created_ts,
            "desc": desc,
            "status": "end",
            "params": params
        }
        main_id = self.small_rna.insert_main_table("sg_geneset_venn", main_info)

        # work_flow import 路径
        task_name = "small_rna.report.geneset_venn"

        # 创建返回给前端的信息，以及work_flow所需参数
        #update_info = json.dumps({str(main_id): "sg_geneset_venn", "status_in": 1})
        update_info = json.dumps({str(main_id): "sg_geneset_venn"})
        options = {
            "main_id": str(main_id),
            "update_info": update_info  # to update sg_status
        }

        # 设置sheet参数对象
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,
                            task_id=task_id,
                            project_sn=project_sn,
                            module_type="workflow")

        # 运行work_flow跟新运行信息，task_info 为字典格式，包含sucess, 和 info 键
        task_info = super(GenesetVennAction, self).POST()

        # 跟新post请求返回信息
        success_info = "Parameters have been stored, geneset venn analysis is ready{}".format(data.get('geneset_ids'))
        fail_info = "Task fails {}".format(task_info.get("info", ""))
        task_info["info"] = success_info if (task_info.get('success', False) or task_info["info"]) else fail_info
        task_info["content"] = {"ids": {"id": str(main_id), "name": name}}

        # 1. sg_geneset中添加{"is_use": 1}信息，用来标记那些geneset被用过了
        # 2. 检测 "_id" 在sg_geneset_venn (在上面创建相应的表)
        # 3. 在sg_geneset_info表中关联 venn 表名 sg_geneset_venn 以及 上面创建的表_id , gene_set_id
        self.small_rna.insert_geneset_info(geneset_ids=str(data.geneset_ids),
                                           col_name="sg_geneset_venn",
                                           col_id=str(main_id))

        return json.dumps(task_info)


if __name__ == "__main__":
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """

        def test_this(self):
            cmd = "python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py "
            cmd += "post "
            cmd += "-fr no "
            cmd += "-c {} ".format("client03")
            cmd += "i/small_rna/geneset_venn "
            cmd += "-b http://192.168.12.102:9090 "
            args = dict(
                task_id="small_rna",
                task_type="1",  # 1 代表instance任务
                submit_location="geneset_venn",
                gene_type="M",  # mirna
                geneset_ids="5bc05d1da4e1af2a979fcd0c,5bc05d74a4e1af0473a02121",
            )

            arg_names, arg_values = args.keys(), args.values()
            cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
            print(cmd)
            os.system(cmd)

    # 调用
    unittest.main()
