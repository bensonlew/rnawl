#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/5/8 15:55
@file    : merge_results.py
"""
import csv
import os
import unittest

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
from mbio.packages.lnc_rna.lnc_identification.predict_common_tools import CommandFactory
from collections import OrderedDict
import pandas as pd

class MergePredictionsAgent(Agent):
    def __init__(self, parent):
        super(MergePredictionsAgent, self).__init__(parent)
        # 'fasta_file', 'gtf_file', 'hexamer_dat', 'logit_model',
        # , 'transcript_len', 'cnci_score',
        # 'orf_len', 'cpat_score', 'taxonmy', 'exon_num'
        options = [
            {'name': 'predictions_dir', 'type': 'string'},
            {'name': 'tools', 'type': 'string'},
        ]
        self.add_option(options)
        self.step.add_steps("merge_predictions")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.merge_predictions.start()
        self.step.update()

    def step_end(self):
        self.step.merge_predictions.finish()
        self.step.update()

    def check_options(self):
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'


    def end(self):
        # self.add_upload_dir(self.output_dir)
        super(MergePredictionsAgent, self).end()


class MergePredictionsTool(Tool):
    def __init__(self, config):
        super(MergePredictionsTool, self).__init__(config)
        self.python_path = "program/Python/bin/python"
        self.env_path = self.config.SOFTWARE_DIR + '/program/Python/bin'
        self.set_environ(PATH=self.env_path)


    def cmd_runner(self, cmd_name, cmd, check_stat=True):
        cmd_obj = self.add_command(cmd_name, cmd)
        cmd_obj.run()
        # self.wait(cmd_obj)
        if check_stat is True:
            self.check_stat(cmd_obj)
        return cmd_obj

    def check_stat(self, *cmd_objs):
        self.wait(*cmd_objs)
        for cmd_obj in cmd_objs:
            if cmd_obj.return_code == 0:
                self.logger.info('%s：运行完成' % cmd_obj.cmd)
            elif cmd_obj.return_code in (1, -9):
                self.add_state('memory_limit', 'memory is low!')
            else:
                self.set_error('%s: 运行错误%s' % cmd_obj.cmd)


    def run(self):
        super(MergePredictionsTool, self).run()
        tools = self.option('tools').strip().split(',')
        self.logger.info('%s：tool有这几个' % tools)
        preds_dir = self.option('predictions_dir')
        summary=OrderedDict()
        combine_infos = []
        for tool in tools:
            file = os.path.join(preds_dir, tool + '_output.txt')
            df = pd.read_table(file)
            coding_list= list(df[df["label"]=="coding"]["transcript_id"])
            summary[tool]=coding_list
            if not tool == "pfam":
                # df[tool+"_coding"]=df["label"].apply(lambda x : 0 if x=="coding" else 1)
                df[tool + "_coding"] = df["label"]
                df[tool + "_score"] = df["score"]
                df.set_index("transcript_id",inplace = True)
                combine_infos.append(df[tool+"_coding"])
                combine_infos.append(df[tool+"_score"])
                self.logger.info('%s：现在在做这个tool' % tool)
            else:
                # df[tool + "_coding"] = df["label"].apply(lambda x: 0 if x == "coding" else 1)
                df[tool + "_coding"] = df["label"]
                df.set_index("transcript_id", inplace=True)
                combine_infos.append(df[tool + "_coding"])
        combie_pd=pd.concat(combine_infos,axis=1,join="outer")
        combie_pd.to_csv(os.path.join(self.option("predictions_dir"),"combine_resluts"),sep="\t")
        with open(os.path.join(self.output_dir,"predict_stat.xls"),"w") as p:
            for i in summary:
                p.write(i+"\t")
                p.write(";".join(summary[i])+"\n")

        self.end()


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """

        def test(self):
            import random
            from mbio.workflows.single import SingleWorkflow
            from biocluster.wsheet import Sheet
            data = {
                "id": "MergePredictions_" + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "tool_lab.coding_prediction.merge_predictions",
                "instant": False,
                "options": dict(
                    predictions_dir="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/coding_predict/data",
                    tools="cpat,cpc,pfam",
                    # cnci_score=0,
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()

    unittest.main()
