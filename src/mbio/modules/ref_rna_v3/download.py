# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.module import Module
from mbio.packages.ref_rna_v3.functions import modlfuncdeco
import os
import unittest

class DownloadModule(Module):
    '''
    last_modify: 2019.06.13
    '''
    def __init__(self, work_id):
        super(DownloadModule, self).__init__(work_id)
        options = [
            {'name': 's3_file_list', 'type': 'infile', 'format': 'ref_rna_v3.common'},
            {'name': 'my_file_list', 'type': 'outfile', 'format': 'ref_rna_v3.common'},
        ]
        self.add_option(options)

    @modlfuncdeco
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    @modlfuncdeco
    def run(self):
        self.run_tool()
        super(DownloadModule, self).run()

    @modlfuncdeco
    def run_tool(self):
        self.tools = [self.set_tool(n, line) for n, line in enumerate(open(self.option('s3_file_list').path))]
        if len(self.tools) == 1:
            self.tools[0].on('end', self.set_output)
        else:
            self.on_rely(self.tools, self.set_output)
        for tool in self.tools:
            tool.run()

    @modlfuncdeco
    def set_tool(self, n, line):
        items = line.strip().split('\t')
        options = {'ifile': items[0]}
        if len(items) == 2:
            options.update({'basename': items[1]})
        self.step.add_steps('download_{}'.format(n))
        download = self.add_tool('ref_rna_v3.download')
        download.set_options(options)
        download.on('start', self.set_step, {'start': getattr(self.step, 'download_{}'.format(n))})
        download.on('end', self.set_step, {'end': getattr(self.step, 'download_{}'.format(n))})
        return download

    def set_output(self, event):
        lines = list()
        for tool in self.tools:
            source = tool.option('ofile').path
            link_name = os.path.join(self.output_dir, os.path.basename(source))
            if os.path.exists(link_name):
                os.remove(link_name)
            os.link(source, link_name)
            lines.append(link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        my_file_list = os.path.join(self.work_dir, os.path.basename(self.option('s3_file_list').path))
        open(my_file_list, 'w').writelines('{}\n'.format(line) for line in lines)
        self.option('my_file_list').set_path(my_file_list)
        self.end()

    @modlfuncdeco
    def end(self):
        super(DownloadModule, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'download_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'ref_rna_v3.download',
            'instant': False,
            'options': {
                's3_file_list': '/mnt/lustre/users/sanger/sg-users/qinjincheng/download/MJ20180726044/i-sanger_145213.txt'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
