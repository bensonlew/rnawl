# -*- coding: utf-8 -*-

import datetime
import json
import os
import unittest

import web

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig


class GenesetVennAction(WholeTranscriptomeController):
    def __init__(self):
        super(GenesetVennAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']  # 必须有的元素
        basic_args += ['geneset_id', 'level']  # 页面选择的参数：基因集，表达水平（基因，转录本） ps：全转录组需要添加新的参数category（）
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument:  %s" % arg, 'code': 'C2901901',
                        'variables': variables}
                return json.dumps(info)
        sg_task = self.whole_transcriptome.get_task_info(data.task_id)
        project_sn = sg_task["project_sn"]
        task_id = data.task_id
        # create main table record
        geneset_id_list = str(data.geneset_id).split(',')
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            level=data.level,
            # geneset_id=geneset_id_list,
            geneset_id=','.join(geneset_id_list),

        )
        if hasattr(data, "category"):
            params.update({
                "category": data.category
            })
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "GenesetVenn" + '_' + data.level + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='Geneset venn analysis main table',
            params=params,
            status="end",
            version="v1"
        )
        main_id = self.whole_transcriptome.insert_main_table('geneset_venn', main_info)

        # prepare option for workflow
        options = {
            "task_id": data.task_id,
            "main_id": str(main_id),
            "geneset_id": data.geneset_id,
            "submit_location": data.submit_location,
            "update_info": json.dumps({str(main_id): "geneset_venn"})  # to update sg_status
        }

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'whole_transcriptome.report.geneset_venn'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(GenesetVennAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
            },
            'info': "成功存储参数，可进行venn分析"
        }
        # 更新基因集的使用信息
        self.whole_transcriptome.insert_geneset_info(data.geneset_id, "geneset_venn", str(main_id))
        return json.dumps(task_info)


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "i/whole_transcriptome/geneset_venn "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_36088",
            task_type="1",
            submit_location="GenesetVenn",
            level="G",
            geneset_id="5da135c217b2bf14e012f744,5db0332a17b2bf3cdfb10639,5db0337e17b2bf3ce0b10634,5dbba38617b2bf6d13a3a818",
            category="mRNA"
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
