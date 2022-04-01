# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import shutil
import unittest
import glob
from biocluster.module import Module


class BuscoModule(Module):
    def __init__(self, work_id):
        super(BuscoModule, self).__init__(work_id)
        options = [
            {'name': 'fa_dir', 'type': 'infile', 'format': 'denovo_rna_v2.common_dir'},
            {'name': 'odb9', 'type': 'string'},
        ]
        self.add_option(options)
        self.tools = list()

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def run(self):
        super(BuscoModule, self).run()
        self.run_busco()

    def run_busco(self):
        for i in glob.glob(self.option('fa_dir').path + '/*'):
            self.busco = self.add_tool("tool_lab.busco")
            self.busco.set_options({
                'fa': i,
                'odb9': self.option('odb9')
            })
            self.tools.append(self.busco)
        else:
            self.on_rely(self.tools, self.set_output)
        for tool in self.tools:
            tool.run()

    def set_output(self):
        for tool in self.tools:
            busco_result = tool.option('busco_result').path
            summary_result = tool.option('summary_result').path
            sample_name = os.path.basename(tool.option('fa').path).split('.fa')[0]
            busco_result_new = os.path.join(self.output_dir, '{}_short_summary_busco_result.xls'.format(sample_name))
            summary_result_new = os.path.join(self.output_dir, '{}_summary_result.xls'.format(sample_name))
            if os.path.isfile(busco_result_new):
                os.remove(busco_result_new)
            os.link(busco_result, busco_result_new)
            if os.path.isfile(summary_result_new):
                os.remove(summary_result_new)
            os.link(summary_result, summary_result_new)
        self.end()

    def end(self):
        super(BuscoModule, self).end()




class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run this script to do test.
    '''

    def test_ath(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'Str_predict_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'tool_lab.str',
            'instant': False,
            'options': {
                'bamlist': '/mnt/ilustre/users/sanger-dev/workspace/20210114/Str_STR2926/RnaseqMapping/output/bamlist',
                'ref_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Solanum_lycopersicum/SL4.0_ITAG4.0/dna/S_lycopersicum_chromosomes.4.00.fa',
                'sample_list': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/STR/data/sample_list.txt'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_ath')])
    unittest.TextTestRunner(verbosity=2).run(suite)
