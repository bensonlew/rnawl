#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/4/2 9:05
@file    : single_db.py
@author  : zhipeng.zhao
"""
import csv
import os
import time
import unittest

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool


class SingleDbAgent(Agent):
    def __init__(self, parent):
        super(SingleDbAgent, self).__init__(parent)
        options = [
            # {'name': 'lnc_ids', 'type': 'infile', 'format': 'lnc_rna.common'},
            # {'name': 'fasta', 'type': 'infile', 'format': 'lnc_rna.lnc_fasta'},
            {'name': 'ref_fa', 'type': 'infile', 'format': 'lnc_rna.lnc_fasta'},
            {'name': 'gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
            {'name': 'db_name', 'type': 'string', 'required': True},
            {'name': 'out_fa', 'type': 'outfile', 'format': 'lnc_rna.common'},
            {'name': 'out_gtf', 'type': 'outfile', 'format': 'lnc_rna.common'},
            {'name': 'out_ids_map', 'type': 'outfile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps("single_db")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.single_db.start()
        self.step.update()

    def step_end(self):
        self.step.single_db.finish()
        self.step.update()

    def check_options(self):
        for name in ('ref_fa', 'gtf'):
            if not self.option(name).is_set:
                raise OptionError(name + ' must be assigned')
        return True

    def set_resource(self):
        self._cpu = 3
        self._memory = "5G"

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(SingleDbAgent, self).end()


class SingleDbTool(Tool):

    def run_cmd(self, cmd_name, cmd, is_wait=True, shell=False):
        self.logger.debug(cmd_name + ': ' + cmd + '%s' % type(cmd_name))
        cmd_obj = self.add_command(str(cmd_name), cmd, shell=shell)
        if shell is True:
            cmd_obj.software_dir = ''
            cmd_obj._start_run_time = int(time.time())
        cmd_obj.run()
        if is_wait is True:
            self._check_stat(cmd_obj)
            return
        return cmd_obj

    def _check_stat(self, *cmd_objs):
        self.wait(*cmd_objs)
        for cmd_obj in cmd_objs:
            if cmd_obj.return_code == 0:
                self.logger.info("{} Finished successfully".format(cmd_obj.name))
            elif cmd_obj.return_code is None:
                self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_obj.name))
                cmd_obj.rerun()
                self.wait()
                if cmd_obj.return_code is 0:
                    self.logger.info("{} Finished successfully".format(cmd_obj.name))
                else:
                    self.set_error("{} Failed. >>>{}".format(cmd_obj.name, cmd_obj.cmd))
            else:
                self.set_error("{} Failed. >>>{}".format(cmd_obj.name, cmd_obj.cmd))

    def create_ids_map(self):
        out_file = os.path.join(self.output_dir, 'ids_map.xls')
        check_set = set()
        db_name = self.option('db_name')
        out_gtf = os.path.join(self.output_dir, 'lncrna.gtf')
        with self.option('gtf') as gtf_handler, open(out_file, 'w') as out_handler, \
                open(out_gtf, 'w') as gtf_out_handler:
            out_handler.write(
                'transcript_id\tsource\t{}_transcript_id\tgene_id\tgene_name\n'.format(self.option('db_name')))
            line_demo = '{transcript_id}\t{db_name}\t{transcript_id}\t{gene_id}\t{gene_name}\n'
            for line_splits, line in gtf_handler:
                attr_dict = line_splits[8]
                if line_splits[2] == "gene":
                    continue
                transcript_id = attr_dict['transcript_id']
                if transcript_id in check_set and line_splits[2] == "transcript":
                    continue
                check_set.add(transcript_id)
                if 'gene_name' not in attr_dict:
                    attr_dict['gene_name'] = ''

                # ids map file output
                out_handler.write(line_demo.format(db_name=db_name, **attr_dict))

                # gtf 输出
                line_splits[8] = ' '.join('{key} "{val}";'.format(key=k, val=v) for k, v in attr_dict.items())
                line_splits[1] = db_name
                line = '\t'.join(str(i) for i in line_splits) + '\n'
                gtf_out_handler.write(line)

        self.option('out_gtf').set_path(out_gtf)
        self.option('out_ids_map').set_path(out_file)

    def gtf2fa(self):
        gffread = 'bioinfo/lnc_rna/cufflinks-2.2.1.Linux_x86_64/gffread'
        new_fa = os.path.join(self.output_dir, 'lncrna.fa')
        cmd = '{} {} '.format(gffread, self.option('gtf').path)
        cmd += '-g {} '.format(self.option('ref_fa').path)
        cmd += '-w {}'.format(new_fa)

        self.run_cmd('gffread', cmd, is_wait=True)
        self.option('out_fa').set_path(new_fa)

    def run(self):
        super(SingleDbTool, self).run()
        self.create_ids_map()
        self.gtf2fa()
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
                "id": "SingDB_" + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "lnc_rna.lnc_db.single_db",
                "instant": False,
                "options": dict(
                    # {'name': 'gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
                    # {'name': 'db_name', 'type': 'string', 'required': True},
                    gtf='/mnt/ilustre/users/sanger-dev/workspace/20190318/LncDb_lnc_db_workflow_2785_2967/'
                        'GtfFilter/output/Homo_sapiens.GRCh38.95.gtf',
                    ref_fa='/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/'
                           'Homo_sapiens/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fa',
                    db_name='esembl'
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()


    unittest.main()
