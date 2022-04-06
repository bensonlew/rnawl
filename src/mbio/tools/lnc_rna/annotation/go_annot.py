# -*- coding: utf-8 -*-
# __author__ = 'wangbixuan'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest
from mbio.packages.rna.annot_config import AnnotConfig

class GoAnnotAgent(Agent):
    '''
    last_modify: 2019.02.14
    '''
    def __init__(self, parent):
        super(GoAnnotAgent, self).__init__(parent)
        options = [
            {'name': 'blast2go_annot', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'go_list', 'type': 'outfile', 'format': 'lnc_rna.go_list'},
            {'name': 'go_level2', 'type': 'outfile', 'format': 'lnc_rna.go_level2'},
            {"name": "pir_version", "type": "string", "default": "2019"}, #pir database version
            {"name": "go_version", "type": "string", "default": "2019"}, #pir database version

        ]
        self.add_option(options)
        self.step.add_steps('go_annot')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self._memory_increase_step = 40
        # Delivery to specific queue (BLAST2GO)
        self.queue = 'BLAST2GO'

    def step_start(self):
        self.step.go_annot.start()
        self.step.update()

    def step_end(self):
        self.step.go_annot.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        if self.option('blast2go_annot').is_set:
            self.logger.debug('{} = {}'.format('blast2go_annot', self.option('blast2go_annot').prop['path']))
            self.infile_size = os.path.getsize(self.option('blast2go_annot').prop['path'])
        else:
            raise OptionError('annotation table of BLAST2GO must be provided')
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(int(float(self.infile_size) / 1024 ** 3 * 8 + 20))

    def end(self):
        super(GoAnnotAgent, self).end()

class GoAnnotTool(Tool):
    def __init__(self, config):
        super(GoAnnotTool, self).__init__(config)
        self.python = 'miniconda2/bin/python'
        self.go_merge_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/go_merge.py')
        self.go_annot_py = os.path.join(self.config.PACKAGE_DIR, 'rna/annotation/go_annotation2.py')
        self.go_split_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/go_split.py')
        self.query_gos_list = os.path.join(self.work_dir, 'query_gos.list')
        self.b2g_host = 'localhost'
        self.b2g_user = 'biocluster102'
        self.b2g_passwd = 'sanger-dev-123'
        self.b2g_db = 'b2gdb'
        self.go_detail_xls = os.path.join(self.work_dir, 'go_level_detail.xls')
        self.go_level2 = os.path.join(self.work_dir, 'go_level2.xls')
        self.idmapping_db = AnnotConfig().get_file_path(
            file ="idmapping.tb",
            db = "pir",
            version = self.option("pir_version"))
        self.go_obo = AnnotConfig().get_file_dict(db="go", version=self.option("go_version"))['go']

    def run(self):
        super(GoAnnotTool, self).run()
        self.run_go_merge()
        self.run_go_annot()
        self.run_go_split()
        self.set_output()
        self.end()

    def run_go_merge(self):
        cmd = '{} {}'.format(self.python, self.go_merge_py)
        cmd += ' -i {}'.format(self.option('blast2go_annot').prop['path'])
        cmd += ' -o {}'.format(self.query_gos_list)
        cmd_name = 'run_go_merge'
        self.run_code(cmd_name, cmd)

    def run_go_annot(self):
        cmd = '{} {} {} {} {}'.format(
            self.python,
            self.go_annot_py,
            self.query_gos_list,
            self.work_dir,
            self.go_obo
        )
        cmd_name = 'run_go_annot_t'
        self.run_code(cmd_name, cmd)
        '''
        cmd = '{} {} {} {} {} {} {} {}'.format(
            self.python,
            self.go_annot_py,
            self.query_gos_list,
            self.b2g_host,
            self.b2g_user,
            self.b2g_passwd,
            self.b2g_db,
            self.work_dir
        )
        cmd_name = 'run_go_annot'
        self.run_code(cmd_name, cmd)
        '''

    def run_go_split(self):
        cmd = '{} {} {} {}'.format(
            self.python,
            self.go_split_py,
            self.go_detail_xls,
            self.work_dir
        )
        cmd_name = 'run_go_split'
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
        go_list = os.path.join(self.output_dir, 'query_gos.list')
        if os.path.exists(go_list):
            os.remove(go_list)
        os.link(self.query_gos_list, go_list)
        self.logger.info('succeed in linking {} to {}'.format(self.query_gos_list, go_list))
        self.option('go_list').set_path(go_list)
        basenames = [
            'go_level2.xls', 'go_level3.xls', 'go_level4.xls'
        ]
        for basename in basenames:
            source = os.path.join(self.work_dir, basename)
            link_name = os.path.join(self.output_dir, basename)
            if os.path.exists(link_name):
                os.remove(link_name)
            os.link(source, link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        infiles = ["level4.stat.tsv", 'level3.stat.tsv', 'level2.stat.tsv']
        outfiles = ['go1234level_statistics.xls', 'go123level_statistics.xls', 'go12level_statistics.xls']
        for inf,item in zip(infiles, outfiles):
            linkfile = self.output_dir + '/' + item
            if os.path.exists(linkfile):
                os.remove(linkfile)
            os.link(self.work_dir + '/' + inf, linkfile)

        go_level2 = os.path.join(self.output_dir, 'go_level2.xls')
        if os.path.exists(go_level2):
            os.remove(go_level2)
        os.link(self.go_level2, go_level2)
        self.logger.info('succeed in linking {} to {}'.format(self.go_level2, go_level2))
        self.option('go_level2').set_path(go_level2)
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
            'id': 'go_annot_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.annotation.go_annot',
            'instant': False,
            'options': {
                'blast2go_annot': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_filter/output/go/blast2go.filter.xls',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
