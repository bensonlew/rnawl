#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/4/10 13:37
@file    : extr_known_lnc.py
"""
import csv
import json
import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool


class ExtrKnownLncAgent(Agent):
    def __init__(self, parent):
        super(ExtrKnownLncAgent, self).__init__(parent)
        options = [
            {'name': 'new_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'lnc_db_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
        ]
        self.add_option(options)
        self.step.add_steps("extr_known_lnc")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.extr_known_lnc.start()
        self.step.update()

    def step_end(self):
        self.step.extr_known_lnc.finish()
        self.step.update()

    def check_options(self):
        for name in ('new_gtf', 'lnc_db_gtf'):
            if not self.option(name).is_set:
                raise Exception('%s is not set' % name)
        return True

    def set_resource(self):
        self._cpu = 3
        self._memory = '6G'


class ExtrKnownLncTool(Tool):
    def __init__(self, config):
        super(ExtrKnownLncTool, self).__init__(config)
        env_path = os.path.join(self.config.SOFTWARE_DIR, '/program/Python/bin')
        self.set_environ(PATH=env_path)

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

    def cuffcompare(self):
        cuffcompare = 'bioinfo/rna/cufflinks-2.2.1/cuffcompare'
        cuffcompare_dir = os.path.join(self.work_dir, 'cuffcompare')
        if not os.path.exists(cuffcompare_dir):
            os.makedirs(cuffcompare_dir)
        else:
            os.system('rm -r %s' % cuffcompare_dir)
            os.makedirs(cuffcompare_dir)
        new_gtf = self.option('new_gtf').path
        gtf_basename = os.path.basename(new_gtf)

        ln_file = os.path.join(cuffcompare_dir, gtf_basename)
        return_code = os.system('ln -s {} {}'.format(new_gtf, ln_file))
        if return_code != 0:
            os.system('cp {} {}'.format(new_gtf, cuffcompare_dir))

        cmd = '{cuffcompare} -r {ref_gtf} {gtf_file}'
        cmd = cmd.format(cuffcompare=cuffcompare, ref_gtf=self.option('lnc_db_gtf').path, gtf_file=ln_file)
        self.cmd_runner('cuffcompare', cmd, check_stat=True)

        ref_map_file = os.path.join(cuffcompare_dir, 'cuffcmp.' + gtf_basename + '.refmap')
        # ref_gene_id     ref_id  class_code      cuff_id_list
        known_trans = set()
        with open(ref_map_file) as in_handler:
            for line_dic in csv.DictReader(in_handler, delimiter='\t'):
                class_code = line_dic['class_code']
                if class_code != '=':
                    continue
                for item in line_dic['cuff_id_list'].split(','):
                    known_trans.add(item.split('|')[1])

        return known_trans

    def output_file(self, ids_set):
        detail_json = os.path.join(self.output_dir, 'known_lnc_in_new.json')
        with open(detail_json, 'w') as json_handler:
            json.dump([i for i in ids_set], json_handler, indent=4)

    def run(self):
        super(ExtrKnownLncTool, self).run()
        ids_set = self.cuffcompare()
        self.output_file(ids_set=ids_set)
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
                "id": "extr_known_lnc_" + str(random.randint(1, 10000)) + '_' + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "lnc_rna.lncrna_identification.extr_known_lnc",
                "instant": False,
                "options": dict(
                    new_gtf='/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/assemble/method-cufflinks/'
                            'output/NewTranscripts/new_transcripts.gtf',
                    lnc_db_gtf='/mnt/ilustre/users/sanger-dev/workspace/20190409/LncDb_lnc_db_workflow_3956_9130/'
                               'output/lncrna.gtf',
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()


    unittest.main()
