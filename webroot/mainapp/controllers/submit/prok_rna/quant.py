# -*- coding: utf-8 -*-
import web
import json
import datetime
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mbio.api.to_file.prok_rna import *
from mainapp.libs.signature import check_sig
import unittest
import os
import copy
from mainapp.controllers.submit.prok_rna.apibase import ApiBase


class QuantAction(ProkRNAController):
    def __init__(self):
        super(QuantAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['method', 'exp_type', ]
        # check arg
        # exp_type是使用的软件
        '''
        参数实例，通过网址 http://www.tsg.com/data/temp_data/refrna.html 查询
                Array
        (
            [exp_type] => TPM
            [method] => kallisto
            [submit_location] => exp_detail
            [task_id] => tsanger_31035
            [task_type] => 2
        )
        '''
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        sg_task = self.prok_rna.get_task_info(data.task_id)
        project_sn = sg_task["project_sn"]
        id2name = sg_task["id2name"]  # 得到id和name之间的对应关系
        libtype = self.prok_rna.get_libtype(data.task_id)  # 随机选择一张主表获得libtype
        read_len = self.prok_rna.get_mean_read_len(data.task_id)
        result_dir = self.db['sg_exp'].find_one({'task_id': data.task_id, 'method': 'RSEM', 'exp_type': 'TPM'})['result_dir']
        # create main table record
        exp_type = data.exp_type.upper()
        ### 此处和邱于涵商量好存储字段大小写格式
        #工作流跑的默认存到mongo库params的method是RSEM，salmon，kallisto；exp_type是TPM，FPKM
        #交互分析运行后存到mongo库params的method是RSEM，salmon，kallisto；exp_type是TPM，FPKM
        if data.method.lower() == "rsem":
            quant_method = data.method.upper()
        else:
            quant_method = data.method.lower()
        time_now = datetime.datetime.now()
        params = dict(
            task_id=data.task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            method=quant_method,
            exp_type=data.exp_type,
        )
        main_info = dict(
            project_sn=project_sn,
            task_id=data.task_id,
            # name=name,
            exp_type=exp_type.upper(),
            method=quant_method,
            libtype=libtype,
            exp_level='mRNA+sRNA',
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            status="start",
            result_dir=result_dir
        )
        if libtype:
            if type(params) == dict:
                # params_ms = copy.deepcopy(params)
                params_m = copy.deepcopy(params)
                # params_s = copy.deepcopy(params)
                # if data.exp_level:
                #     params_m.update({"exp_level": data.exp_level})
                #     params_jm = json.dumps(params_m, sort_keys=True, separators=(',', ':'))
                # else:
                params_jm = json.dumps(params_m, sort_keys=True, separators=(',', ':'))
                # params_s.update({"exp_level": "sRNA"})
                # params_js = json.dumps(params_s, sort_keys=True, separators=(',', ':'))
                # params_ms.update({"exp_level": "mRNA_sRNA"})
                # params_jms = json.dumps(params_ms, sort_keys=True, separators=(',', ':'))
            # main_info_ms = copy.deepcopy(main_info)  # 不进行复制会导致重复的key, https://blog.csdn.net/Xiongtm/article/details/77650448
            main_info_m = copy.deepcopy(main_info)
            # main_info_s = copy.deepcopy(main_info)

            # exp_type, quant_method = data.exp_type.upper(), data.method
            # name = "Exp" + '_' + data.exp_level + '_' + quant_method + '_' + exp_type.upper() + '_'
            name = "Exp" + '_' + quant_method.upper() + '_' + exp_type.upper() + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            # 用两个 __ 分割这个字符串，分割得到的第一个就代表exp_level应该显示的字段
            # main_info_ms.update({"desc": '{} exp main table'.format("mRNA and sRNA"), "exp_level": "mRNA_sRNA", "params": params_jms, "name": "mRNA_sRNA__" + name})
            # main_id_ms = self.prok_rna.insert_main_table('sg_exp', main_info_ms)
            main_info_m.update({"desc": '{} exp main table'.format("mRNA"),  "params": params_jm, "name": name})
            main_id_m = self.prok_rna.insert_main_table('sg_exp', main_info_m)
            # main_info_s.update({"desc": '{} exp main table'.format("sRNA"), "exp_level": "sRNA", "params": params_js, "name": "sRNA__" + name})
            # main_id_s = self.prok_rna.insert_main_table('sg_exp', main_info_s)
            options = {
                "submit_location": data.submit_location,
                "task_type": data.task_type,
                "raw_task_id": data.task_id,
                "method": data.method,
                "libtype": libtype,
                "read_len": read_len,
                "transcriptome": self.use_s3(sg_task["assemble_fa"]),
                "biotype": data.task_id,
                "fastq": sg_task['fastq'],  # 返回一个fq.list， 包含样品名和路径
                # "main_id_ms": str(main_id_ms),
                "main_id_m": str(main_id_m),
                # "main_id_s": str(main_id_s),
                "exp_type": exp_type,
                "id2name": self.use_s3(id2name),
                # "exp_level": data.exp_level,
                # "update_info": json.dumps({str(main_id_ms): "sg_exp", str(main_id_m): "sg_exp", str(main_id_s): "sg_exp"})  # to update sg_status
                "update_info": json.dumps({str(main_id_m): "sg_exp"})  # to update sg_status
            }
            # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
            to_files = ["prok_rna.export_gene_type(biotype)"]
            task_name = 'prok_rna.report.quant'
            self.set_sheet_data(name=task_name,
                                options=options,
                                main_table_name=data.method + '_expression_quant',
                                module_type="workflow",
                                to_file=to_files,
                                project_sn=project_sn,
                                task_id=data.task_id)

            task_info = super(QuantAction, self).POST()
            task_info['content'] = {
                'ids': [{
                        # 'id': str(main_id_ms),
                        # 'name': 'mRNA_sRNA'
                        # }, {
                            'id': str(main_id_m),
                            'name': 'mRNA'
                        # }, {
                        #     'id': str(main_id_s),
                        #     'name': 'sRNA'
                        }]
            }
        else:
            if type(params) == dict:
                # params.update({"exp_level": data.exp_level})
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))

            #exp_type, quant_method = data.exp_type.upper(), data.method
            name = "Exp" + '_' + quant_method.upper() + '_' + exp_type.upper() + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")

            main_info.update({"desc": '{} exp main table'.format("mRNA"), "params": params, "name": name})
            main_id_m = self.prok_rna.insert_main_table('sg_exp', main_info)
            main_info.update({"desc": '{} exp main table'.format("mRNA")})
            # prepare option for tool
            options = {
                "submit_location": data.submit_location,
                "task_type": data.task_type,
                "raw_task_id": data.task_id,
                "method": data.method,
                "libtype": libtype,
                "read_len": read_len,
                "transcriptome": self.use_s3(sg_task["assemble_fa"]),
                "biotype": data.task_id,
                "fastq": sg_task['fastq'],  # 返回一个fq.list， 包含样品名和路径
                "main_id_m": str(main_id_m),
                "exp_type": exp_type,
                "id2name": self.use_s3(id2name),
                "update_info": json.dumps({str(main_id_m): "sg_exp"})  # to update sg_status
            }
            # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
            to_files = ["prok_rna.export_gene_type(biotype)"]
            task_name = 'prok_rna.report.quant'
            self.set_sheet_data(name=task_name,
                                options=options,
                                main_table_name=data.method + '_expression_quant',
                                module_type="workflow",
                                to_file=to_files,
                                project_sn=project_sn,
                                task_id=data.task_id)

            task_info = super(QuantAction, self).POST()
            task_info['content'] = {
                'ids': {
                    'id': str(main_id_m),
                    'name': "mRNA"
                    }
            }
        #  对sg_specimen_group_compare和sg_specimen_group增加字段
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.prok_rna.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.prok_rna.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/prok_rna/quant "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="prok_rna_srna",
            task_type="2",
            submit_location="exp_detail",
            method='RSEM',
            exp_type='fpkm',
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
