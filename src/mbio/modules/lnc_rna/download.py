# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.module import Module
import os
import unittest

class DownloadModule(Module):
    '''
    last_modify: 2019.03.19
    '''
    def __init__(self, work_id):
        super(DownloadModule, self).__init__(work_id)
        options = [
            {'name': 's3_file_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'list_type', 'type': 'string', 'default': None},
            {'name': 'local_file_list', 'type': 'outfile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.tools = list()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        super(DownloadModule, self).run()
        self.run_tool()

    def run_tool(self):
        for n, line in enumerate(open(self.option('s3_file_list').path)):
            if self.option('list_type') == 'loc2name':
                s3_file, file_name = line.strip().split('\t')
            else:
                s3_file = line.strip()
            self.step.add_steps('download_{}'.format(n))
            download = self.add_tool('lnc_rna.download')
            options = {'file_in': s3_file}
            if self.option('list_type') == 'loc2name':
                options.update({'file_name': file_name})
            download.set_options(options)
            download.on('start', self.set_step, {'start': getattr(self.step, 'download_{}'.format(n))})
            download.on('end', self.set_step, {'end': getattr(self.step, 'download_{}'.format(n))})
            self.tools.append(download)
        if len(self.tools) == 1:
            self.tools[0].on('end', self.set_output)
        else:
            self.on_rely(self.tools, self.set_output)
        for tool in self.tools:
            tool.run()

    def set_output(self, event):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        lines = list()
        for tool in self.tools:
            source = tool.option('file_out').path
            link_name = os.path.join(self.output_dir, os.path.basename(source))
            if os.path.exists(link_name):
                os.remove(link_name)
            os.link(source, link_name)
            if self.option('list_type') == 'loc2name':
                lines.append('{}\t{}'.format(link_name, tool.option('file_name')))
            else:
                lines.append(link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        local_file_list = os.path.join(self.output_dir, os.path.basename(self.option('s3_file_list').path))
        open(local_file_list, 'w').writelines(['{}\n'.format(line) for line in lines])
        self.option('local_file_list').set_path(local_file_list)
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.end()

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
            'id': 'download_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'module',
            'name': 'lnc_rna.download',
            'instant': False,
            'options': {
                's3_file_list': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/download/bam.list'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_loc2name(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'download_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'module',
            'name': 'lnc_rna.download',
            'instant': False,
            'options': {
                's3_file_list': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/download/bam.loc2name.list',
                'list_type': 'loc2name'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()