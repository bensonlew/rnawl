# -*- coding: utf-8 -*-
# __author__ = 'jinlinfang, qinjincheng'

from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import os
import shutil
import re
import unittest

class RmatsModule(Module):
    '''
    last_modify: 2019.03.13
    '''
    def __init__(self, work_id):
        super(RmatsModule, self).__init__(work_id)
        options = [
            {'name': 'group_table', 'type': 'infile', 'format': 'sample.group_table'},
            {'name': 'control_table', 'type': 'infile', 'format': 'sample.control_table'},
            {'name': 'bam_loc', 'type': 'infile', 'format': 'lnc_rna.loc2name'},
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            # ['paired', 'single']
            {'name': 'read_type', 'type': 'string', 'default': 'paired'},
            # ['fr-unstranded', 'fr-firststrand', 'fr-secondstrand']
            {'name': 'lib_type', 'type': 'string', 'default': 'fr-unstranded'},
            {'name': 'read_length', 'type': 'int', 'default': 150},
            {'name': 'cstat', 'type': 'float', 'default': 0.0001},
        ]
        self.add_option(options)
        self.tools = list()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        if self.option('group_table').is_set:
            self.logger.debug('{} = {}'.format('group_table', self.option('group_table').prop['path']))
        else:
            raise OptionError('GROUP TABLE file must be provided')
        if self.option('control_table').is_set:
            self.logger.debug('{} = {}'.format('control_table', self.option('control_table').prop['path']))
        else:
            raise OptionError('CONTROL TABLE file must be provided')
        if self.option('bam_loc').is_set:
            self.logger.debug('{} = {}'.format('bam_loc', self.option('bam_loc').prop['path']))
        else:
            raise OptionError('BAM LOC file must be provided')
        if self.option('ref_gtf').is_set:
            self.logger.debug('{} = {}'.format('ref_gtf', self.option('ref_gtf').prop['path']))
        else:
            raise OptionError('input GTF file must be provided')
        self.logger.debug('{} = {}'.format('read_type', self.option('read_type')))
        self.logger.debug('{} = {}'.format('lib_type', self.option('lib_type')))
        self.logger.debug('{} = {}'.format('cstat', self.option('cstat')))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        super(RmatsModule, self).run()
        self.run_rmats()

    def run_rmats(self):
        num, vs_list = self.option('control_table').get_control_info()
        group_spname = self.option('group_table').get_group_spname()
        name2loc_dict = self.option('bam_loc').prop['locations']
        for n, (ctrl, test) in enumerate(vs_list):
            self.step.add_steps('rmats_{}'.format(n))
            rmats = self.add_tool('lnc_rna.structure.rmats')
            b1_config = os.path.join(rmats.work_dir, '{}.config'.format(test))
            open(b1_config, 'w').write(','.join([name2loc_dict[sp] for sp in group_spname[test]]))
            b2_config = os.path.join(rmats.work_dir, '{}.config'.format(ctrl))
            open(b2_config, 'w').write(','.join([name2loc_dict[sp] for sp in group_spname[ctrl]]))
            rmats.set_options({
                'ref_gtf': self.option('ref_gtf'),
                'B1_config': b1_config,
                'B2_config': b2_config,
                'read_type': self.option('read_type'),
                'lib_type': self.option('lib_type'),
                'read_length': self.option('read_length'),
                'cstat': self.option('cstat'),
            })
            rmats.vs_name = '{}_vs_{}'.format(test, ctrl)
            rmats.on('start', self.set_step, {'start': getattr(self.step, 'rmats_{}'.format(n))})
            rmats.on('end', self.set_step, {'end': getattr(self.step, 'rmats_{}'.format(n))})
            self.tools.append(rmats)
        if len(self.tools) == 1:
            self.tools[0].on('end', self.set_output)
        else:
            self.on_rely(self.tools, self.set_output)
        for tool in self.tools:
            tool.run()

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        result_files = (
            'all_events_detail_big_table.txt',
            'psi_stats.file.txt',
            'event_stats.file.txt',
            'event_type.file.txt'
        )
        for tool in self.tools:
            output_dir = os.path.join(self.output_dir, tool.vs_name)
            if not os.path.exists(output_dir):
                os.mkdir(output_dir)
            for output_file in os.listdir(tool.output_dir):
                source = os.path.join(tool.output_dir, output_file)
                link_name = os.path.join(output_dir, output_file)
                if os.path.exists(link_name):
                    os.remove(link_name)
                if os.path.split(output_file)[1] in result_files:
                    os.link(source, link_name)
                if re.search(r'^fromGTF\.(RI|A3SS|A5SS|SE|MXE)\.alter_id\.txt$', os.path.split(output_file)[1]):
                    os.link(source, link_name)
                if re.search(r'fromGTF\.novelEvents\.(RI|A3SS|A5SS|SE|MXE)\.alter_id\.txt$', os.path.split(output_file)[1]):
                    os.link(source, link_name)
                if re.search(r'(RI|A3SS|A5SS|SE|MXE)\.MATS\.JC\.alter_id\.txt$', os.path.split(output_file)[1]):
                    os.link(source, link_name)
                if re.search(r'(RI|A3SS|A5SS|SE|MXE)\.MATS\.JC\.alter_id\.psi_info\.txt$', os.path.split(output_file)[1]):
                    os.link(source, link_name)
                if re.search(r'(RI|A3SS|A5SS|SE|MXE)\.MATS\.JCEC\.alter_id\.txt$', os.path.split(output_file)[1]):
                    os.link(source, link_name)
                if re.search(r'(RI|A3SS|A5SS|SE|MXE)\.MATS\.JCEC\.alter_id\.psi_info\.txt$', os.path.split(output_file)[1]):
                    os.link(source, link_name)
                self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.end()

    def end(self):
        super(RmatsModule, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run script to do test.
    '''
    #
    # def test(self):
    #     import random
    #     from mbio.workflows.single import SingleWorkflow
    #     from biocluster.wsheet import Sheet
    #     data = {
    #         'id': 'rmats_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
    #         'type': 'module',
    #         'name': 'lnc_rna.rmats',
    #         'instant': False,
    #         'options': {
    #             'group_table': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/rmats/group.5.txt',
    #             'control_table': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/rmats/control.5.txt',
    #             'bam_loc': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/rmats/bam_loc.txt',
    #             'ref_gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf'
    #         }
    #     }
    #     wsheet = Sheet(data=data)
    #     wf = SingleWorkflow(wsheet)
    #     wf.run()

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'rmats_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'module',
            'name': 'lnc_rna.rmats',
            'instant': False,
            'options': {
                'group_table': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/rmats/group.txt',
                'control_table': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/rmats/control.txt',
                'bam_loc': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/rmats/bam_loc.txt',
                'ref_gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()