# -*- coding: utf-8 -*-
# __author__ = 'scp'

from biocluster.agent import Agent
import os
from biocluster.tool import Tool
import unittest

class GmapBuildAgent(Agent):
    def __init__(self, parent):
        super(GmapBuildAgent, self).__init__(parent)
        options = [
            {'name': 'reference_in', 'type': 'infile', 'format': 'ref_genome_db.fasta'},
            {'name': 'organism', 'type': 'string', 'default': ''},
            {'name': 'gmap_db_dir', 'type': 'outfile', 'format': 'lnc_rna.common_dir'},
            {'name': 'gmap_db_name', 'type': 'string'},
        ]
        self.add_option(options)

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        self.logger.debug('{} - {}'.format('reference_in', self.option('reference_in').prop['path']))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 20
        self._memory = '120G'

    def end(self):
        super(GmapBuildAgent, self).end()

class GmapBuildTool(Tool):
    def __init__(self, config):
        super(GmapBuildTool, self).__init__(config)
        self.gmap_build = 'bioinfo/lnc_rna/gmap/bin/gmap_build'

    def run(self):
        super(GmapBuildTool, self).run()
        self.run_gmap_build()
        self.set_output()
        self.end()

    def run_gmap_build(self):
        cmd = '{} -D {} -d {} {}'.format(
            self.gmap_build,
            self.work_dir,
            self.option('organism'),
            self.option('reference_in').prop['path']
        )
        cmd_name = 'gmap_build'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd):
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('succeed in running {}'.format(cmd_name))
        elif command.return_code is None:
            self.logger.warn('fail to run {}, try again'.format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info('succeed in rerunning {}'.format(cmd_name))
            else:
                self.set_error('fail to rerun {}, abord'.format(cmd_name))
        else:
            self.set_error('fail to run {}, abord'.format(cmd_name))

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        self.option("gmap_db_dir", self.work_dir + '/' + self.option("organism"))
        self.option("gmap_db_name", self.option("organism"))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'gmap_build_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db.gmap_build',
            'instant': False,
            'options': {
                'reference_in': '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Gallus_gallus/ensembl/Gallus_gallus.GRCg6a.dna.toplevel.fa',
                'organism': 'Gallus_gallus',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()