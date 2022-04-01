# -*- coding: utf-8 -*-
import os
import web
import json
import unittest
import datetime
import re
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.core.base import Base
from mainapp.controllers.core.basic import Basic
from biocluster.core.function import filter_error_info
from biocluster.wpm.client import *
# __author__ = 'shicaiping'


class CheckFileAction(Base):
    def __init__(self):
        super(CheckFileAction, self).__init__()

    @check_sig
    def POST(self):
        data = web.input()
        print("web input: {}".format(data))
        # check
        for arg in ['genome', 'gtf', 'in_type']:
            if not hasattr(data, arg):
                info = {'success': False, 'info': '缺少参数：{}'.format(arg)}
                return json.dumps(info)
            elif arg == "in_type" and data[arg] not in ["gtf", "gff"]:
                info = {'success': False, 'info': '参数in_type仅支持gtf/gff'}
                return json.dumps(info)

        workflow_id = "CheckFile_{}".format(datetime.datetime.now().strftime("%H%M%S%f")[:-3])
        data_json = {
            'id': workflow_id,
            'stage_id': 0,
            'name': "prok_rna.report.check_file",
            'type': 'workflow',
            'client': data.client,
            "IMPORT_REPORT_DATA": False,
            "IMPORT_REPORT_AFTER_END": False,
            'options': {
                "genome": data.genome,
                "gtf": data.gtf,
                "in_type": data.in_type,
            }
        }
        workflow_client = Basic(data=data_json, instant=True)
        try:
            run_info = workflow_client.run()
            run_info['info'] = filter_error_info(run_info['info'])
            return json.dumps(run_info)
        except Exception as e:
            return json.dumps({"success": False, "info": "running error: %s" % filter_error_info(str(e))})


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/lustre/users/sanger-dev/wpm2/sanger_bioinfo/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += '-dbversion {} '.format(1)
        cmd += "i/prok_rna/check_file "
        cmd += "-b http://wpm2.sanger.com "
        args = dict(
            genome="/mnt/lustre/users/sanger-dev/wpm2/workspace/20210722/Prokrna_2s9o_kmffc8ojen00hhqvbci5op/Download/output/GCF_000009345.1_ASM934v1/GCF_000009345.1_ASM934v1_genomic.fna",
            gtf="/mnt/lustre/users/sanger-dev/wpm2/workspace/20210722/Prokrna_2s9o_kmffc8ojen00hhqvbci5op/Download/output/GCF_000009345.1_ASM934v1/GCF_000009345.1_ASM934v1_genomic.gtf",
            in_type='gtf'
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
