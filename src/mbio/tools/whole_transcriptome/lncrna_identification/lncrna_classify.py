#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/3/11 10:42
@file    : lncrna_classify.py
"""

import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest


class LncrnaClassifyAgent(Agent):
    """
    lnc rna target prediction
    """

    def __init__(self, parent):
        super(LncrnaClassifyAgent, self).__init__(parent)
        options = [
            {'type': 'infile', 'name': 'mrna_gtf', 'format': 'lnc_rna.gtf'},
            {'type': 'infile', 'name': 'lncrna_gtf', 'format': 'lnc_rna.gtf'},
            {'type': 'string', 'name': 'out_file_name', 'default': 'lncRNA_classifications.xls'}
        ]
        self.add_option(options)

    def check_options(self):
        for name in ('mrna_gtf', 'lncrna_gtf'):
            if not self.option(name).is_set:
                raise OptionError("%s must be set")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"


class LncrnaClassifyTool(Tool):
    """
    lnc rna target prediction
        生成结果文件: lncRNA_classifications.xls
    """

    def __init__(self, config):
        super(LncrnaClassifyTool, self).__init__(config)
        self.python_path = 'miniconda2/bin/python'
        self.perl_path = 'miniconda2/bin/perl'
        self.bed2intersect_script = self.config.PACKAGE_DIR + '/lnc_rna/bed2intersect.py'
        self.gtf2bed_script = self.config.PACKAGE_DIR + '/lnc_rna/gtf2bed.pl'

    def run_command(self, cmd_name, cmd, is_wait=True):
        cmd_obj = self.add_command(cmd_name, cmd)
        cmd_obj.run()
        if is_wait is True:
            self._cmd_wait(cmd_obj)
        else:
            return cmd_obj

    def _cmd_wait(self, *cmd_objs):
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

    def gtf2bed(self):
        tmp_dir = os.path.join(self.work_dir, 'tmp_beds')
        if os.path.isdir(tmp_dir):
            pass
        else:
            os.makedirs(tmp_dir)
        mrna_bed = os.path.join(tmp_dir, 'mrna.bed')
        lncrna_bed = os.path.join(tmp_dir, 'lncrna.bed')

        data = ((self.option('lncrna_gtf').prop['path'], lncrna_bed, 'lncrna_gtf2bed'),
                (self.option('mrna_gtf').prop['path'], mrna_bed, 'mrna_gtf2bed'))
        cmd_objs = []
        for file, out_file, cmd_name in data:
            cmd = '{} {} {} {}'.format(self.perl_path, self.gtf2bed_script, file, out_file)
            cmd_objs.append(self.run_command(cmd_name, cmd, is_wait=False))
        self._cmd_wait(*cmd_objs)

        return mrna_bed, lncrna_bed

    def classify_lncrna(self):
        mrna_bed, lncrna_bed = self.gtf2bed()
        out_file = os.path.join(self.output_dir, self.option('out_file_name'))

        cmd = '{} {} '.format(self.python_path, self.bed2intersect_script)
        cmd += ' {} '.format(lncrna_bed)
        cmd += ' {} '.format(mrna_bed)
        cmd += ' {} '.format(out_file)
        self.run_command('lnc_b2d_intersect', cmd, is_wait=True)

        return out_file

    def run(self):
        super(LncrnaClassifyTool, self).run()
        self.classify_lncrna()
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
                "id": "lncrna_classify_" + str(random.randint(1, 10000)) + '_' + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "lnc_rna.lncrna_identification.lncrna_classify",
                "instant": False,
                "options": dict(
                    mrna_gtf=r'/mnt/ilustre/users/sanger-dev/workspace/20190402/LncDb_lnc_db_workflow_6536_6335/GtfFilter/output/mrna.gtf',
                    lncrna_gtf=r'/mnt/ilustre/users/sanger-dev/workspace/20190402/LncDb_lnc_db_workflow_6536_6335/GtfFilter/output/Homo_sapiens.GRCh38.95.gtf',
                    # lncrna_gtf=r'/mnt/ilustre/users/sanger-dev/workspace/20190312/Single_lncrna_identify_7888_1085/LncrnaIdentify/NewLncrnaPredict/MergePredictions/output/new_lncrna.gtf',
                    out_file_name='known_lncRNA_classifications.xls'
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()

    unittest.main()
