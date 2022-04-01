# -*- coding: utf-8 -*-
# __author__ = 'jinlinfang,shicaiping,qinjincheng'

import os
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.ref_rna_v2.rmats_process_func import *
from mbio.packages.ref_rna_v2.rmats_process_func import process_single_rmats_output_dir
import re
from mbio.files.sequence.file_sample import FileSampleFile
import unittest

class RmatsModule(Module):
    '''
    last_modify: 2019.06.10
    '''
    def __init__(self, work_id):
        super(RmatsModule, self).__init__(work_id)
        options = [
            {'name': 'sample_bam_dir', 'type': 'infile', 'format': 'align.bwa.bam_dir'},
            {'name': 'rmats_control', 'type': 'infile', 'format': 'sample.control_table'},
            {'name': 'group_table', 'type': 'infile', 'format': 'sample.group_table'},
            {'name': 'gname', 'type': 'string', 'default': 'group'},  # 分组方案名称
            {'name': 'seq_type', 'type': 'string', 'default': 'paired'},  # 两个选项：'paired'  or ’single‘
            {'name': 'read_length', 'type': 'int', 'default': 150},
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'gene_structure.gtf'},  # 一定要设置
            {'name': 'lib_type', 'type': 'string', 'default': 'fr-unstranded'},  # 建库类型
            {'name': 'as_diff', 'type': 'float', 'default': 0.0001},
        ]
        self.add_option(options)
        self.rmats_bam_tools = []
    
    def check_options(self):
        if not self.option('ref_gtf'):
            raise OptionError('必须设置参考基因组注释文件（ref_genome.gtf）', code = '23701102')
        if self.option('as_diff') < 0 or self.option('as_diff') >= 1:
            raise OptionError('差异剪接假设检验的置信度p应该: 0=< p <1', code = '23701103')
        if self.option('lib_type') not in ('fr-unstranded', 'fr-firststrand', 'fr-secondstrand'):
            raise OptionError('rMATS识别的建库类型只可设置为fr-unstranded或fr-firststrand或fr-secondstrand', code = '23701104')
        if self.option('seq_type') not in ('paired', 'single'):
            raise OptionError('rMATS识别的测序读长类型只可为paired（双端测序）或single（单端测序）', code = '23701105')
        if self.option('group_table').is_set and not self.option('gname'):
            raise OptionError('有分组文件时必须传入分组方案名字', code = '23701106')
        if self.option('group_table').is_set and self.option('gname') not in self.option('group_table').prop[
            'group_scheme']:
            raise OptionError('传入分组方案名字不在分组文件内', code = '23701107')
        if not self.option('rmats_control').is_set:
            raise OptionError('必须设置输入文件：上下调对照组参考文件', code = '23701108')
        return True
    
    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def run(self):
        super(RmatsModule, self).run()
        self.multi_rmats_bam_run()

    def multi_rmats_bam_run(self):
        vs_group_path_pair_lst, vs_list = self.get_group_str()
        n = 0
        for vs_pair in vs_group_path_pair_lst:
            rmats_bam = self.add_tool('ref_rna_v2.rmats_bam')
            vs_name = vs_list[n][0] + '_vs_' + vs_list[n][1]
            n = n + 1
            self.step.add_steps('rmats_bam_{}'.format(n))
            rmats_bam.set_options({
                'ref_gtf': self.option('ref_gtf').path,
                'seq_type': self.option('seq_type'),
                'read_length': self.option('read_length'),
                'A_group_bam': vs_pair[0],
                'B_group_bam': vs_pair[1],
                'lib_type': self.option('lib_type'),
                'cut_off': self.option('as_diff')
            })
            rmats_bam.vs_name = vs_name
            step = getattr(self.step, 'rmats_bam_{}'.format(n))
            step.start()
            rmats_bam.on('end', self.finish_update, 'rmats_bam_{}'.format(n))
            self.rmats_bam_tools.append(rmats_bam)
        self.on_rely(self.rmats_bam_tools, self.set_output)
        for tool in self.rmats_bam_tools:
            tool.run()

    def get_group_str(self):
        '''
        :return: [('a1_1.bam,a1_2.bam', 'b1_1.bam,b1_2.bam'), ('a2_1.bam,a2_2.bam', 'b2_1.bam,b2_2.bam'), ......]
        '''
        a_b_bam_path_tuple_lst = []
        sample_bams = [f for f in os.listdir(self.option('sample_bam_dir').path) if re.match(r'.*\.bam$', f)]
        sample_path_dic = {}
        for sample_bam in sample_bams:
            m = re.match(r'(\S+)\.bam$', sample_bam)
            sample_name = m.group(1)
            # 字典格式为：{sample_name: sample_bam_abs_path}
            sample_path_dic[sample_name] = os.path.join(self.option('sample_bam_dir').path, sample_bam)
        self.logger.info('sample_path_dic 为：{}'.format(sample_path_dic))
        # control_info：该比较方案下的比较对数量以及每个比较对的详情，如：
        # 3, [('A', 'B'), ('A', 'C'), ('B', 'C')]
        num, vs_list = self.option('rmats_control').get_control_info()
        self.logger.info('control_info为：{}，{}'.format(num, vs_list))
        # group_spname：该分组方案下分组类别对应的样本信息详细的字典，如：
        # {'A': ['A1', 'A2', 'A3'], 'B': ['B1', 'B2', 'B3']}
        group_sample_dic = self.option('group_table').get_group_spname()
        self.logger.info('group_sample_dic为：{}'.format(group_sample_dic))
        for vs_pair in vs_list:
            a_group_name = vs_pair[0]
            b_group_name = vs_pair[1]
            self.logger.debug('b_group_name（实验）为： {}'.format(b_group_name))
            self.logger.debug('a_group_name（对照）为： {}'.format(a_group_name))
            a_group_samples = group_sample_dic[a_group_name]
            a_group_path_lst = []
            b_group_path_lst = []
            b_group_samples = group_sample_dic[b_group_name]
            for a_group_sample in a_group_samples:
                a_group_path_lst.append(sample_path_dic[a_group_sample])
            for b_group_sample in b_group_samples:
                b_group_path_lst.append(sample_path_dic[b_group_sample])
            a_group_path_str = ','.join(a_group_path_lst)
            b_group_path_str = ','.join(b_group_path_lst)
            a_b_bam_path_tuple_lst.append((a_group_path_str, b_group_path_str))
        else:
            return a_b_bam_path_tuple_lst, vs_list

    def get_list(self):
        list_path = os.path.join(self.option('sample_bam_dir').path, 'list.txt')
        file_sample = FileSampleFile()
        file_sample.set_path(list_path)
        samples = file_sample.get_list()
        return samples

    def set_output(self):
        self.logger.info('set output')
        for rmats_bam_tool in self.rmats_bam_tools:
            output_dir = os.path.join(self.output_dir, rmats_bam_tool.vs_name)
            if not os.path.exists(output_dir):
                os.mkdir(output_dir)
            outfiles = os.listdir(rmats_bam_tool.output_dir)
            for f in outfiles:
                f_path = os.path.join(rmats_bam_tool.output_dir, f)
                target = os.path.join(output_dir, f)
                if os.path.exists(target):
                    os.remove(target)
                if os.path.split(f)[1] in ('all_events_detail_big_table.txt', 'psi_stats.file.txt', 'event_stats.file.txt', 'event_type.file.txt'):
                    os.link(f_path, target)
                if re.search(r'^fromGTF\.(RI|A3SS|A5SS|SE|MXE)\.alter_id\.txt$', os.path.split(f)[1]):
                    os.link(f_path, target)
                if re.search(r'fromGTF\.novelEvents\.(RI|A3SS|A5SS|SE|MXE)\.alter_id\.txt$', os.path.split(f)[1]):
                    os.link(f_path, target)
                if re.search(r'(RI|A3SS|A5SS|SE|MXE)\.MATS\.JC\.alter_id\.txt$', os.path.split(f)[1]):
                    os.link(f_path, target)
                if re.search(r'(RI|A3SS|A5SS|SE|MXE)\.MATS\.JC\.alter_id\.psi_info\.txt$', os.path.split(f)[1]):
                    os.link(f_path, target)
                if re.search(r'(RI|A3SS|A5SS|SE|MXE)\.MATS\.JCEC\.alter_id\.txt$', os.path.split(f)[1]):
                    os.link(f_path, target)
                if re.search(r'(RI|A3SS|A5SS|SE|MXE)\.MATS\.JCEC\.alter_id\.psi_info\.txt$', os.path.split(f)[1]):
                    os.link(f_path, target)
                if re.search(r'A_group_bam.txt', os.path.split(f)[1]):
                    os.link(f_path, target)
                if re.search(r'B_group_bam.txt', os.path.split(f)[1]):
                    os.link(f_path, target)
        self.logger.info('set output done')
        self.end()

    def end(self):
        super(RmatsModule, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run script to do test.
    '''
    def test(self):
        import random
        from biocluster.wsheet import Sheet
        from mbio.workflows.single import SingleWorkflow
        data = {
            'id': 'rmats_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'ref_rna_v2.rmats',
            'instant': False,
            'options': dict(
                sample_bam_dir='/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/test',
                rmats_control='/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/test_control.txt',
                group_table='/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/test_group.txt',
                lib_type='fr-unstranded',
                ref_gtf='/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/ref_and_new.gtf',
                seq_type='paired'
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
