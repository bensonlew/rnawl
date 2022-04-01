# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError
import glob
import time

class CodingPredictWorkflow(Workflow):
    """
    Used for cds to protein code

    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CodingPredictWorkflow, self).__init__(wsheet_object)
        options = [
            # fasta文件
            {"name": "fasta_file", "type": "infile", "format": "ref_rna_v2.common"},
            #extra infos
            {'name': 'genome_name', 'type': 'string', 'default': 'Homo_sapiens'},
            # cpc params
            {'name': 'cpc', 'type': 'bool', 'default': True},
            {'name': 'cpc_score', 'type': 'float', 'default': 0.5},
            # cnci params
            {'name': 'cnci', 'type': 'bool', 'default': True},
            {'name': 'cnci_score', 'type': 'float', 'default': 0},
            {'name': 'taxonmy', 'type': 'string', 'default': 'Animal'},
            # cpat parmas
            {'name': 'pfamscan', 'type': 'bool', 'default': True},
            {'name': 'cpat', 'type': 'bool', 'default': True},
            {'name': 'cpat_score', 'type': 'float', 'default': 0.5},
            {'name': 'hexamer_dat', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'logit_model', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {"name": "main_id", "type": "string"},
            # {"name": "min_len", "type": "int", },
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.logger.info('物种是{}'.format(self.option("genome_name")))
        self.merge_tool = self.add_tool("tool_lab.coding_prediction.merge_predictions")
        self.cpc_tool = None
        self.cpat_tool = None
        self.cnci_tool = None
        self.pfam_tool = None
        self.predict_funcs = []
        self.predict_tools=[]
        self.tools_names = []
        self.set_options(self._sheet.options())
        if self.option('cpc') is True:
            self.cpc_tool = self.add_tool("tool_lab.coding_prediction.cpc_predictor")
            self.tools_names.append('cpc')
            self.step.add_steps('cpc_predict')
            self.predict_tools.append(self.cpc_tool)
            self.predict_funcs.append(self.cpc_predict)
        if self.option('cnci') is True:
            self.cnci_tool = self.add_tool("tool_lab.coding_prediction.cnci_predictor")
            self.tools_names.append('cnci')
            self.step.add_steps('cnci_predict')
            self.predict_tools.append(self.cnci_tool)
            self.predict_funcs.append(self.cnci_predict)
        if self.option('cpat') is True:
            # if not (self.option('hexamer_dat').is_set and self.option('logit_model').is_set):
            #     raise OptionError("使用cpat预测，hexamer_dat 和 logit_model 文件", code="23700801")
            if not self.option("genome_name") in ['Homo_sapiens', 'Mus_musculus', 'Danio_rerio', 'Drosophila_melanogaster']:
                raise OptionError("使用cpat预测，必须在Homo_sapiens，Mus_musculus，Danio_rerio，Drosophila_melanogaster中选择")
            # soft_num += 1
            if self.option("genome_name") in ['Homo_sapiens', 'Mus_musculus', 'Danio_rerio', 'Drosophila_melanogaster']:
                species2cpat = {
                    'Homo_sapiens': "Human",
                    'Mus_musculus': "Mouse",
                    'Danio_rerio': "Zebrafish",
                    'Drosophila_melanogaster': "Fly"
                }
                self.hexamer_dat = self.config.SOFTWARE_DIR + "/bioinfo/lnc_rna/CPAT-1.2.4/dat/{}_Hexamer.tsv".format(
                    species2cpat[self.option("genome_name")])
                self.logit_model = self.config.SOFTWARE_DIR + "/bioinfo/lnc_rna/CPAT-1.2.4/dat/{}_logitModel.RData".format(
                    species2cpat[self.option("genome_name")])
            self.cpat_tool = self.add_tool("tool_lab.coding_prediction.cpat_predictor")
            self.tools_names.append('cpat')
            self.step.add_steps('cpat_predict')
            self.predict_tools.append(self.cpat_tool)
            self.predict_funcs.append(self.cpat_predict)

        if self.option('pfamscan') is True:
            # soft_num += 1
            self.pfam_tool = self.add_tool("tool_lab.coding_prediction.pfam_scan_predictor")
            self.step.add_steps('pfamscan_predict')
            self.tools_names.append('pfam')
            self.predict_tools.append(self.pfam_tool)
            self.predict_funcs.append(self.pfamscan_predict)
        self.step.add_steps('predictions_merge')
        # self.set_options(self._sheet.options())

    # def check_options(self):
    #     if self.option('cpc') is True:
    #         self.cpc_tool = self.add_tool("tool_lab.coding_prediction.cpc_predictor")
    #         self.tools_names.append('cpc')
    #         self.step.add_steps('cpc_predict')
    #         self.predict_tools.append(self.cpc_tool)
    #         self.predict_funcs.append(self.cpc_predict)
    #     if self.option('cnci') is True:
    #         self.cnci_tool = self.add_tool("tool_lab.coding_prediction.cnci_predictor")
    #         self.tools_names.append('cnci')
    #         self.step.add_steps('cnci_predict')
    #         self.predict_tools.append(self.cnci_tool)
    #         self.predict_funcs.append(self.cnci_predict)
    #     if self.option('cpat') is True:
    #         # if not (self.option('hexamer_dat').is_set and self.option('logit_model').is_set):
    #         #     raise OptionError("使用cpat预测，hexamer_dat 和 logit_model 文件", code="23700801")
    #         if not self.option("genome_name") in ['Homo_sapiens', 'Mus_musculus', 'Danio_rerio', 'Drosophila_melanogaster']:
    #             raise OptionError("使用cpat预测，必须在Homo_sapiens，Mus_musculus，Danio_rerio，Drosophila_melanogaster中选择")
    #         # soft_num += 1
    #         if self.option("genome_name") in ['Homo_sapiens', 'Mus_musculus', 'Danio_rerio', 'Drosophila_melanogaster']:
    #             species2cpat = {
    #                 'Homo_sapiens': "Human",
    #                 'Mus_musculus': "Mouse",
    #                 'Danio_rerio': "Zebrafish",
    #                 'Drosophila_melanogaster': "Fly"
    #             }
    #             self.hexamer_dat = self.config.SOFTWARE_DIR + "/bioinfo/lnc_rna/CPAT-1.2.4/dat/{}_Hexamer.tsv".format(
    #                 species2cpat[self.option("genome_name")])
    #             self.logit_model = self.config.SOFTWARE_DIR + "/bioinfo/lnc_rna/CPAT-1.2.4/dat/{}_logitModel.RData".format(
    #                 species2cpat[self.option("genome_name")])
    #         self.cpat_tool = self.add_tool("tool_lab.coding_prediction.cpat_predictor")
    #         self.tools_names.append('cpat')
    #         self.step.add_steps('cpat_predict')
    #         self.predict_tools.append(self.cpat_tool)
    #         self.predict_funcs.append(self.cpat_predict)
    #
    #     if self.option('pfamscan') is True:
    #         # soft_num += 1
    #         self.pfam_tool = self.add_tool("tool_lab.coding_prediction.pfam_scan_predictor")
    #         self.step.add_steps('pfamscan_predict')
    #         self.tools_names.append('pfam')
    #         self.predict_tools.append(self.pfam_tool)
    #         self.predict_funcs.append(self.pfamscan_predict)
    #
    #     self.logger.info('物种是{}'.format(self.option("genome_name")))
    #
    #     return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def cpc_predict(self):
        """
            {'name': 'fasta_file', 'type': 'string'},
            {'name': 'cpc_score', 'type': 'float', 'default': 0.5},
        :return:
        """
        self.cpc_tool.set_options({
            "fasta_file": self.option("fasta_file").prop["path"],
            "cpc_score": self.option("cpc_score")
        })
        self.cpc_tool.on('start', self.set_step, {'start': self.step.cpc_predict})
        self.cpc_tool.on('end', self.set_output, {'file_regex': 'cpc_output.txt', 'tool_name': 'cpc'})
        self.cpc_tool.on('end', self.set_step, {'end': self.step.cpc_predict})
        self.cpc_tool.run()

    def cnci_predict(self):
        """
            {'name': 'taxonmy', 'type': 'string', 'default': 'Animal'},
            {'name': 'fasta_file', 'type': 'string'},
            {'name': 'cnci_score', 'type': 'float', 'default': 0},
        :return:
        """
        self.cnci_tool.set_options({
            "taxonmy": self.option('taxonmy'),
            "fasta_file": self.option("fasta_file").prop["path"],
            "cnci_score": self.option('cnci_score')
        })
        self.cnci_tool.on('start', self.set_step, {'start': self.step.cnci_predict})
        self.cnci_tool.on('end', self.set_output, {'file_regex': 'cnci_output.txt', 'tool_name': 'cnci'})
        self.cnci_tool.on('end', self.set_step, {'end': self.step.cnci_predict})
        self.cnci_tool.run()

    def cpat_predict(self):
        """
            {'name': 'fasta_file', 'type': 'string'},
            {'name': 'hexamer_dat', 'type': 'string', 'default': ''},
            {'name': 'logit_model', 'type': 'string', 'default': ''},
            {'name': 'cpat_score', 'type': 'float', 'default': 0.5}
        :return:
        """
        self.cpat_tool.set_options({
            "fasta_file": self.option("fasta_file").prop["path"],
            "hexamer_dat": self.hexamer_dat,
            "logit_model": self.logit_model,
            "cpat_score": self.option("cpat_score"),
        })
        self.cpat_tool.on('start', self.set_step, {'start': self.step.cpat_predict})
        self.cpat_tool.on('end', self.set_output, {'file_regex': 'cpat_output.txt', 'tool_name': 'cpat'})
        self.cpat_tool.on('end', self.set_step, {'end': self.step.cpat_predict})
        self.cpat_tool.run()

    def pfamscan_predict(self):
        """
            {'name': 'fasta_file', 'type': 'string'}
        :return:
        """
        self.pfam_tool.set_options({
            "fasta_file": self.option("fasta_file").prop["path"],
            "cpu": 12
        })
        self.pfam_tool.on('start', self.set_step, {'start': self.step.pfamscan_predict})
        self.pfam_tool.on('end', self.set_output, {'file_regex': 'pfam_output.txt', 'tool_name': 'pfam'})
        self.pfam_tool.on('end', self.set_step, {'end': self.step.pfamscan_predict})
        self.pfam_tool.run()

    def predictions_merge(self):
        self.merge_tool.set_options({
            'predictions_dir': self.output_dir,
            'tools': ','.join(self.tools_names),
        })
        self.merge_tool.on('start', self.set_step, {'start': self.step.predictions_merge})
        # self.merge_tool.on('end', self.set_output, {'file_regex': 'new_lncrna*', 'tool_name': 'merge'})
        self.merge_tool.on('end', self.set_output, {'file_regex': 'novel_*', 'tool_name': 'merge'})
        self.merge_tool.on('end', self.set_step, {'end': self.step.predictions_merge})
        self.merge_tool.run()

    def set_output(self, event):
        obj = event['bind_object']
        file_regex = os.path.join(obj.output_dir, event['data']['file_regex'])
        # self.logger.debug('%s ==== file name' % file_regex)
        tool_name = event['data']['tool_name']
        files = glob.glob(file_regex)
        for old_path in files:
            file_name = os.path.basename(old_path)
            new_path = os.path.join(self.output_dir, file_name)
            if os.path.exists(new_path):
                os.remove(new_path)
            os.system('ln {old} {new}'.format(old=old_path, new=new_path))
        self.logger.debug(tool_name + ' has completed prediction')
        if tool_name == 'merge':
            self.set_db()
            # self.end()

    def set_db(self):
        seq_stat = self.api.api("tool_lab.coding_predict")
        seq_stat.add_coding_stat_detail(os.path.join(self.merge_tool.output_dir,"predict_stat.xls"),main_id=self.option("main_id"))
        self.end()


    def run_function(self):
         for predict_func in self.predict_funcs:
             predict_func()
             time.sleep(5)
        # self.cpc_predict()


    def run(self):

        self.on_rely(self.predict_tools, self.predictions_merge)
        # self.cnci_predict()
        self.run_function()
        super(CodingPredictWorkflow, self).run()





    # def run_tool(self):
    #     opts = {
    #         'input_file': self.option('input_file'),
    #         'input_format': self.option('input_format'),
    #         'output_format' : self.option("output_format"),
    #
    #     }
    #     self.tool.set_options(opts)
    #     self.tool.on('end', self.set_output)
    #     self.tool.run()

    # def set_output(self):
    #     for file in os.listdir(self.tool.output_dir):
    #         os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
    #     self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "序列编码性预测结果文件",0],
        #     [r'.*\.cpm2tpm\.xls', 'xls', '定量指标cpm转tpm结果文件', 0],
        #     [r'.*\.count2tpm\.xls', 'xls', '定量指标count转tpm结果文件', 0],
        #     [r'.*\.count2cpm\.xls', 'xls', '定量指标count转cpm结果文件', 0],
        #     [r'.*\.count2fpkm\.xls', 'xls', '定量指标count转fpkm结果文件', 0],
        #     [r'.*\.fpkm2tpm\.xls', 'xls', '定量指标fpkm转tpm结果文件', 0],
        #     [r'.*\.cpm2fpkm\.xls', 'xls', '定量指标cpm转fpkm结果文件', 0],
        #     [r'.*\.fpkm2cpm\.xls', 'xls', '定量指标fpkm转cpm结果文件', 0],
         ])
        super(CodingPredictWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.coding_predict import CodingPredictWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "coding_predict" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.coding_predict",
            "options": dict(
                fasta_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/coding_predict/data/basic_filter.fa",
                genome_name="Homo_sapiens",
                # step=500
            )
        }
        wsheet = Sheet(data=data)
        wf =CodingPredictWorkflow(wsheet)
        wf.sheet.id = 'coding_predict'
        wf.sheet.project_sn = 'coding_predict'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)