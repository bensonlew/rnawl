# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import os

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
from biocluster.file import download
from mbio.packages.bac_comp_genome.gc_infos import DNAInfos


class MummerAgent(Agent):
    def __init__(self, parent):
        super(MummerAgent, self).__init__(parent)
        options = [
            {'name': 'mummer', 'type': 'string', 'default': 'nucmer'},  # nucmer promer
            {'name': 'ref', 'type': 'string', 'default':''},
            {'name': 'samples', 'type': 'string'},
            {'name': 'seq_dir', 'type': 'string'},
            {'name': 'super', 'type': 'bool', 'default': True},
            {'name': 'circle_mode', 'type': 'bool', 'default': False}
        ]
        self.add_option(options)
        self.step.add_steps('mummer')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)
        self.queue = 'BLAST'

    def stepstart(self):
        self.step.mummer.start()
        self.step.update()

    def stepfinish(self):
        self.step.mummer.finish()
        self.step.update()

    def check_options(self):
        if not self.option('samples'):
            raise OptionError('请设置待比较的样本名 samples')

    def set_resource(self):
        self._cpu = len(self.option('samples').split(',')) + 1
        self._memory = str(self._cpu * 3) + 'G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ['.', "", "结果输出目录"],
            ])
        super(MummerAgent, self).end()


class MummerTool(Tool):
    def __init__(self, config):
        super(MummerTool, self).__init__(config)
        self.python = '/program/Python/bin/python'
        self.package_path = self.config.PACKAGE_DIR + '/bac_comp_genome/mummer.py'
        self.mum_path = self.config.SOFTWARE_DIR + '/bioinfo/compare_genome/software/MUMmer3.23'

    def run(self):
        super(MummerTool, self).run()
        self.download_seq()
        self.mummer_run()
        self.end()

    def download_seq(self):
        sp_list = self.option('samples').split(',')
        if self.option('ref'):
            sp_list = [self.option('ref'), ] + sp_list
        seq_list = [sp + '.fna' for sp in sp_list]
        self.samples = [download(self.option('seq_dir') + '/' + seq) for seq in seq_list]

    def mummer_run(self):
        self.set_output()
        #samples = self.comb_genome()
        samples = self.samples
        [self.get_len(sp) for sp in samples]
        cmd = self.python + ' ' + self.package_path
        if self.option('ref'):
            cmd += ' -r ' + samples[0]
            samples = samples[1:]
        cmd +=' -m {} -l {} -p {} -o {}'.format(
                self.option('mummer'), ','.join(samples),
                self.mum_path, self.output_dir
            )
        if self.option('super'):
            cmd += ' -s'
        if self.option('circle_mode'):
            cmd += ' -circle_mode'
        command = self.add_command('mummer_run', cmd).run()
        self.wait(command)

        if command.return_code == 0:
            self.logger.info('package mummer_run.py运行成功')
        else:
            self.set_error('package mummer.py运行出错')

    def get_len(self, sample):
        self.logger.info('get {} seq length'.format(sample))
        seqs = DNAInfos(sample)
        seqs.parse_gc(1000000, 1000000)
        
        sp = os.path.basename(sample).split('.fna')[0]
        with open(sp + '.seq_len.xls', 'w') as w:
            [w.write('\t'.join(l) + '\n') for l in seqs.seq_len]

    def set_output(self):
        for root, dirs, files in os.walk(self.output_dir):
            for f in files:
                os.remove(os.path.join(root, f))

    def end(self):

        super(MummerTool, self).end()

    def comb_genome(self):
        sp_list = self.option('samples').split(',')
        samples = []
        i = 1
        for f in sp_list:
            fasta_dir = os.path.join(self.option('seq_dir'), f)
            cmd = '/bin/cat {}/* > {}.fna'.format(fasta_dir, f)
            samples.append(f + '.fna')
            self.add_command('mummer-comb_genome_' + str(i), cmd, shell=True).run()
            i += 1
            self.wait()

        return samples
