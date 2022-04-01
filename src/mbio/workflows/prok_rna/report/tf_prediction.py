# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.workflow import Workflow
import unittest
import os
from biocluster.file import getsize, exists
from biocluster.file import download


class TfPredictionWorkflow(Workflow):
    """
    Prokrna v3.1 转录因子预测
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TfPredictionWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'cds_fa', 'type': 'string'},
            {"name": "evalue", "type": "float", "default": 1e-5},
            {"name": "update_info", "type": "string"},
            {"name": "name", "type": "string"},
            {"name": "main_table_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.predict_tf = self.add_tool("prok_rna.p2tf")
        # self.output_dir = self.preprocess.output_dir

    def run(self):
        if exists(os.path.join(self.option('cds_fa'), 'cds.fa')):
            cds_fa = os.path.join(self.option('cds_fa'), 'cds.fa')
        elif exists(os.path.join(self.option('cds_fa'), 'Sequence_database', 'cds.fa')):
            cds_fa = os.path.join(self.option('cds_fa'), 'Sequence_database', 'cds.fa')
        else:
            self.set_error('file can not find')
        path = self.download_s3_file(cds_fa, 'cds.fa')
        options = {
            "query": path,
            "evalue": self.option("evalue"),
        }
        self.predict_tf.set_options(options)
        self.predict_tf.on('end', self.set_db)
        self.predict_tf.run()
        super(TfPredictionWorkflow, self).run()

    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.predict_tf.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        if os.path.exists(path):
            to_path = path
        elif exists(path):
            download(path, to_path)
        else:
            self.set_error('file can not find')
        return to_path

    def set_db(self):
        """
        导表程序
        """
        tf_api = self.api.api('prok_rna.tf_prediction')
        result_path = os.path.join(self.predict_tf.output_dir, "p2tf_detail.xls")
        tf_api.add_tf_detail(tf_id=self.option("main_table_id"), result=result_path)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.predict_tf.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "转录因子预测文件夹", 0, "150001"],
            ["./*_vs_P2TF.xls", "XLS", "Blast比对结果表", 0, "150003"],
            ["./p2tf_detail.xls", "XLS", "P2TF注释详情表", 0, "150004"],
        ])
        super(TfPredictionWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.prok_rna_v3.report.tf_prediction import PreprocessWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'exp_preprocess_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'dia_v3.report.preprocess',
            'options': {
                # "raw_path": "/mnt/ilustre/users/sanger-dev/workspace/20201126/Diav3_202011261345/raw_treat_ref",
                # "group_table": "/mnt/ilustre/users/sanger-dev/sg-users/xuxi/dia_v3/test_main_workflow/remote_input_from_majorbio_297008/protein_group/group.txt",
                "raw_path": "/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/exp.txt",
                "group_table": "/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/group(1).txt",
                "all_eliminate": "all",
                "all_percent": 90,
                "group_specific": "any",
                "fillna": "min",
                "if_group": "yes",
                "fill_type": "group",
                "main_table_id": "5f7763e317b2bf57d7ed373b",
                # "group_dict": r' ({"C_3": ["C_3_1", "C_3_2", "C_3_3", "C_3_4", "C_3_5"], '
                #               r'"R_3": ["R_3_1", "R_3_2", "R_3_3", "R_3_4", "R_3_5"]}'.replace('"', '\\"'),
            }
        }
        wsheet_object = Sheet(data=data)
        wf = PreprocessWorkflow(wsheet_object)
        wf.sheet.id = 'dia_test'
        wf.sheet.project_sn = 'dia_test'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    unittest.main()