# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import codecs
import datetime
import json
import os
import unittest

import web

from mainapp.controllers.project.tool_lab_controller import ToolLabController
from mainapp.libs.signature import check_sig


class ExtractGffFastaAction(ToolLabController):
    def __init__(self):
        super(ExtractGffFastaAction, self).__init__(instant=False)
        self.basic_args = ["task_id", "genome", "gff", "seq_type", "task_type"]
        self.expected_args = ["no_pseudo", "no_cds", "strand", "chr", "start", "end", "reverse"]
        self.input_data = web.input()

    @check_sig
    def POST(self):
        # check params
        try:
            check = self.check_params()
            if check is not True:
                return check
        except Exception, e:
            info = {'success': False, 'info': '{}'.format(e)}
            return json.dumps(info)

        # create main table
        try:
            params = self.pack_params()
            name = "ExtractGff" + '_' + self.input_data.seq_type + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=self.input_data.project_sn,
                task_id=self.input_data.task_id,
                name=name,
                version="v1",
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='Extract Gff sequence',
                params=params,
                status="start"
            )
            main_id = self.tool_lab.insert_main_table('extract_gff', main_info)
        except Exception, e:
            info = {'success': False, 'info': '{}'.format(e)}
            return json.dumps(info)

        # prepare option for workflow
        try:
            options = {
                "genome": self.input_data.genome,
                "gff": self.input_data.gff,
                "seq_type": self.input_data.seq_type,
                "no_pseudo": self.input_data.no_pseudo,
                "no_cds": self.input_data.no_cds,
                "strand": self.input_data.strand,
                "chr": self.input_data.chr,
                "start": self.input_data.start,
                "end": self.input_data.end,
                "reverse": self.input_data.reverse,
                "update_info": json.dumps({str(main_id): "extract_gff"})  # to update sg_status
            }

            # 把参数交给workflow运行相应的tool
            task_name = 'tool_lab.extract_gff_fasta'
            self.set_sheet_data(name=task_name,
                                options=options,
                                main_table_name=name,  # 设置交互分析结果目录名
                                module_type="workflow",
                                project_sn=self.input_data.project_sn,
                                task_id=self.input_data.task_id)
        except Exception, e:
            info = {'success': False, 'info': '{}'.format(e)}
            return json.dumps(info)

        # 运行workflow 并传回参数
        task_info = super(ExtractGffFastaAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
            }
        }
        return json.dumps(task_info)

    def check_params(self):
        for arg in self.basic_args:
            if not hasattr(self.input_data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2901801', 'variables': variables}
                return json.dumps(info)
            if arg.lower() == "null":
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "%s : is null or NULL" % arg, 'code': 'C2901802',
                        'variables': variables}
                return json.dumps(info)
        return True

    def pack_params(self):
        input_data_dict = dict(self.input_data)
        params_dict = dict()
        args = self.basic_args + self.expected_args
        for each in args:
            if each == "task_type":
                params_dict[each] = int(input_data_dict[each])
            else:
                params_dict[each] = input_data_dict[each]
        return json.dumps(params_dict, sort_keys=True, separators=(',', ':'))

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/tool_lab/extract_gff_fasta "
        cmd += "-b http://bcl.tsg.com "
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        args = dict(
            task_id="extract_gff_test1",
            project_sn="extract_gff_test1",
            task_type="2",
            submit_location="extract_gff",
            name="extract_gff_test1",
            genome=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
            gff=test_dir + "/" + "gtf/Saccharomyces_cerevisiae.R64-1-1.39.gtf",
            seq_type='all',
            reverse='False',
            no_pseudo='False',
            no_cds='False',
            chr="",
            strand="+",
            start='10',
            end='10000'
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
