#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/4/11 14:27
@file    : known_lnc_identify.py
"""

import glob
import os
import unittest

from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class KnownLncIdentifyModule(Module):
    def __init__(self, work_id):
        """
        输出文件：
            known_lncrna_detail.xls
            known_lncrna.fa
            known_lncrna.gtf
            known_lncrna_ids.list
            known_mrna.fa
            known_mrna.gtf
            known_mrna_ids.list
        :param work_id:
        """
        super(KnownLncIdentifyModule, self).__init__(work_id)
        options = [
            {'name': 'biomart', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'biomart_type', 'type': 'string', 'default': 'type1'},
            {'name': 'lnc_db_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
            {'name': 'lnc_db_fa', 'type': 'infile', 'format': 'lnc_rna.lnc_fasta'},
            {'name': 'mrna_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
            {'name': 'mrna_fa', 'type': 'infile', 'format': 'lnc_rna.lnc_fasta'},
            {'name': 'exp_matrix', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'ids_mapping', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'database', 'type': 'string'}
        ]
        self.add_option(options)
        self.step.add_steps('biomart')
        self.biomart_tool = self.add_tool("lnc_rna.lncrna_identification.biomart")
        self.step.add_steps('lnc_classify')
        self.lnc_classify_tool = self.add_tool("lnc_rna.lncrna_identification.lncrna_classify")
        self.step.add_steps('extr_known')
        self.extr_known_tool = self.add_tool("lnc_rna.lncrna_identification.known_lncrna_identify")

    def check_options(self):
        self.step.add_steps('predictions_merge')
        for name in ('biomart', 'lnc_db_gtf', 'lnc_db_fa', 'mrna_gtf', 'mrna_fa', 'exp_matrix'):
            if not self.option(name).is_set:
                raise OptionError("必须输入 %s 文件" % name, code="23700801")
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def biomart(self):
        # 输出: biomart.xls, biomart.json
        self.biomart_tool.set_options({
            'biomart': self.option('biomart'),
            'biomart_type': self.option('biomart_type')
        })
        self.biomart_tool.on('start', self.set_step, {'start': self.step.biomart})
        self.biomart_tool.on('end', self.set_step, {'end': self.step.biomart})
        self.biomart_tool.run()

    def classify(self):
        """
            {'type': 'infile', 'name': 'mrna_gtf', 'format': 'lnc_rna.gtf'},
            {'type': 'infile', 'name': 'lncrna_gtf', 'format': 'lnc_rna.gtf'},
            {'type': 'string', 'name': 'out_file_name', 'default': 'known_lncrna_classifications.xls'}
        """
        self.lnc_classify_tool.set_options({
            'mrna_gtf': self.option('mrna_gtf').path,
            'lncrna_gtf': self.option('lnc_db_gtf').path,
            'out_file_name': 'known_lncrna_classifications.xls'
        })
        self.lnc_classify_tool.on('start', self.set_step, {'start': self.step.lnc_classify})
        self.lnc_classify_tool.on('end', self.set_step, {'end': self.step.lnc_classify})
        self.lnc_classify_tool.run()

    def extr_known(self):
        """
        {'name': 'lnc_db_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
        {'name': 'lnc_db_fa', 'type': 'infile', 'format': 'lnc_rna.lnc_fasta'},
        {'name': 'mrna_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
        {'name': 'mrna_fa', 'type': 'infile', 'format': 'lnc_rna.lnc_fasta'},
        {'name': 'exp_matrix', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
        {'name': 'biomart_json', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
        {'name': 'classify_info', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
        {'name': 'database', 'type': 'string'}
        :return: known_lnc_in_new.json
        """
        self.extr_known_tool.set_options({
            'lnc_db_gtf': self.option('lnc_db_gtf'),
            'lnc_db_fa': self.option('lnc_db_fa'),
            'mrna_gtf': self.option('mrna_gtf'),
            'mrna_fa': self.option('mrna_fa'),
            'exp_matrix': self.option('exp_matrix'),
            'biomart_json': os.path.join(self.biomart_tool.output_dir, 'biomart.json'),
            'classify_info': os.path.join(self.lnc_classify_tool.output_dir, 'known_lncrna_classifications.xls'),
            'database': self.option('database'),
            'ids_mapping': self.option('ids_mapping')
        })
        self.extr_known_tool.on('start', self.set_step, {'start': self.step.extr_known})
        self.extr_known_tool.on('end', self.set_output)
        self.extr_known_tool.run()

    def set_output(self):
        files = glob.glob(self.extr_known_tool.output_dir + '/*')
        for file in files:
            basename = os.path.basename(file)
            if os.path.exists(os.path.join(self.output_dir, basename)):
                os.remove(os.path.join(self.output_dir, basename))
            os.link(file, os.path.join(self.output_dir, basename))
            # os.system('ln {old} {new}'.format(old=file, new=os.path.join(self.output_dir, basename)))
        # os.system('rm {}'.format(' '.join(files)))
        self.end()

    def run(self):
        super(KnownLncIdentifyModule, self).run()
        # for func in self.predict_funcs:
        #     func()
        rely_tools = [self.biomart_tool, self.lnc_classify_tool]
        self.on_rely(rely_tools, self.extr_known)
        self.biomart()
        self.classify()


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
                "id": "known_lnc_identify_" + str(random.randint(1, 10000)),
                "type": "module",
                "name": "lnc_rna.known_lnc_identify",
                "instant": False,
                "options": dict(
                    # {'name': 'biomart', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
                    biomart='/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/'
                            'Homo_sapiens/Ensemble_release_89/biomart/Homo_sapiens.GRCh38.biomart_gene.txt',
                    # {'name': 'biomart_type', 'type': 'string', 'default': 'type1'},
                    biomart_type='type1',
                    # {'name': 'lnc_db_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
                    lnc_db_gtf="/mnt/ilustre/users/sanger-dev/workspace/20190402/LncDb_lnc_db_workflow_6536_6335/GtfFilter/output/Homo_sapiens.GRCh38.95.gtf",
                    # {'name': 'lnc_db_fa', 'type': 'infile', 'format': 'lnc_rna.lnc_fasta'},
                    lnc_db_fa="/mnt/ilustre/users/sanger-dev/workspace/20190402/LncDb_lnc_db_workflow_6536_6335/GtfFilter/output/lncrna_Homo_sapiens.GRCh38.95.fa",
                    # {'name': 'mrna_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
                    mrna_gtf='/mnt/ilustre/users/sanger-dev/workspace/20190402/LncDb_lnc_db_workflow_6536_6335/GtfFilter/output/mrna.gtf',
                    # {'name': 'mrna_fa', 'type': 'infile', 'format': 'lnc_rna.lnc_fasta'},
                    mrna_fa='/mnt/ilustre/users/sanger-dev/workspace/20190402/LncDb_lnc_db_workflow_6536_6335/GtfFilter/output/mrna.fa',
                    # {'name': 'exp_matrix', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
                    exp_matrix="/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncrna_tools/transcript.tpm.matrix",
                    # {'name': 'database', 'type': 'string'}
                    database='ensembl',
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()


    unittest.main()
