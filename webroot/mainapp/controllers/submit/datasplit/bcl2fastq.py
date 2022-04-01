# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

import os
import web
import json
from mainapp.libs.signature import check_sig
from mainapp.models.workflow import Workflow
from mainapp.models.mongo.submit.datasplit.datasplit import Datasplit
from mainapp.controllers.project.datasplit_controller import DatasplitController


class Bcl2fastqAction(DatasplitController):
    """
    数据拆分,一次拆分及自动拆分的接口
    """
    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["status_id", "run_type"]
        for name in params:
            if not hasattr(data, name):
                info = {"success": False, "info": "参数{}不存在".format(name)}
                return json.dumps(info)
        status_id = data.status_id
        run_type = data.run_type
        options = {}
        update_info = {str(status_id): 'seq_status'}
        options['update_info'] = json.dumps(update_info)
        board_data_path, board_seq_model, seq_number, barcode_mismatch = Datasplit().get_bcl2fastq_info(status_id)
        if run_type == 'auto' or run_type == 'split':
            task_name = 'datasplit.bcl2fastq'
            base_path = os.path.join(board_data_path, "Data/Intensities/BaseCalls/")
            if not os.path.exists(base_path):
                info = {"success": False, "info": "下机数据里缺少{}文件夹".format(base_path)}
                return json.dumps(info)
            options["data_path"] = board_data_path
            options['bases_mask'] = board_seq_model
            options['barcode_mismatch'] = barcode_mismatch
            options['sample_sheet'] = status_id
            options['run_type'] = run_type
            options['lib_id'] = status_id
            options['status_id'] = status_id
            to_file = ["datasplit.export_sample_sheet(sample_sheet)", "datasplit.export_lib_id(lib_id)"]
        elif run_type == "qc":
            task_name = 'datasplit.sample_split_qc'
            options['project_params'] = status_id
            options['status_id'] = status_id
            to_file = ["datasplit.export_split_qc_params(project_params)"]
        self.set_sheet_data(name=task_name, options=options, table_id=status_id, to_file=to_file, seq_number=seq_number)
        task_info = super(Bcl2fastqAction, self).POST()
        return json.dumps(task_info)
