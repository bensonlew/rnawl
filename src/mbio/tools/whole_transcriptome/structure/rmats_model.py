# -*- coding: utf-8 -*-
# __author__ = 'jinlinfang, qinjincheng'

from biocluster.agent import Agent
from mbio.packages.whole_transcriptome.functions import toolfuncdeco, runcmd
from biocluster.tool import Tool
import os
import glob
import unittest
import re


class RmatsModelAgent(Agent):
    '''
    last_modify: 2019.06.18
    '''
    def __init__(self, parent):
        super(RmatsModelAgent, self).__init__(parent)
        options = [
            {'name': 'gene_id', 'type': 'string', 'default': None},
            {'name': 'event_file', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'event_type', 'type': 'string', 'default': None},
            {'name': 'l1', 'type': 'string', 'default': None}, # test
            {'name': 'l2', 'type': 'string', 'default': None}, # ctrl
            {'name': 'b1', 'type': 'string', 'default': None}, # test
            {'name': 'b2', 'type': 'string', 'default': None}, # ctrl
            {'name': 'exon_s', 'type': 'int', 'default': 1},
            {'name': 'intron_s', 'type': 'int', 'default': 1},
            {'name': 'background', 'type': 'string', 'default': 'white'},
            {'name': 'density', 'type': 'int', 'default': 600},
            {'name': 'quality', 'type': 'int', 'default': 100},
        ]
        self.add_option(options)
        self.step.add_steps('rmats_model')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    @toolfuncdeco
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def step_start(self):
        self.step.rmats_model.start()
        self.step.update()

    def step_end(self):
        self.step.rmats_model.finish()
        self.step.update()

    def set_resource(self):
        self._cpu = 1
        self._memory = '8G'

    @toolfuncdeco
    def end(self):
        super(RmatsModelAgent, self).end()

class RmatsModelTool(Tool):
    def __init__(self, config):
        super(RmatsModelTool, self).__init__(config)
        self._LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/program/ImageMagick/lib/:$LD_LIBRARY_PATH"
        self.samtools_path = os.path.join(self.config.SOFTWARE_DIR, "bioinfo/align/samtools-1.7/samtools")
        self.program = {
            'python': 'miniconda2/bin/python',
            'convert': 'program/ImageMagick/bin/convert'
        }
        self.script = {
            'rmats_model_prepare': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v3/structure/rmats_model_prepare.py'),
            'rmats2sashimiplot': os.path.join(
                self.config.SOFTWARE_DIR,
                'bioinfo/rna/rmats2sashimiplot-master/src/rmats2sashimiplot/rmats2sashimiplot.py'
            )
        }
        self.file = {
            'event_file': os.path.join(self.work_dir, 'event_file.txt'),
            'event_file_new': os.path.join(self.work_dir, 'event_file_new.txt')
        }
        self.dir = {
            'Sashimi_plot': os.path.join(self.work_dir, 'Sashimi_plot')
        }

    @toolfuncdeco
    def run(self):
        super(RmatsModelTool, self).run()
        self.run_samtools_chr()
        self.run_rmats_model_prepare()
        self.run_rmats2sashimiplot()
        self.run_convert()
        self.set_output()
        self.end()

    @toolfuncdeco
    def run_samtools_chr(self):
        bam = self.option('b1').strip().split(',')[0]
        cmd = '{} view -H {} > bam_head.txt'.format(self.samtools_path, bam)
        cmd_name = 'run_samtools'
        command = self.add_command(cmd_name, cmd, shell=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        else:
            self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "33704402")
        self.chr_list = list()
        with open('bam_head.txt', 'r') as b:
            for line in b.readlines():
                if line.startswith('@SQ'):
                    SN = line.strip().split('\t')[1]
                    chr = str(re.findall(r'SN:(.*?)$', SN)[0])
                    self.chr_list.append(chr)

    @toolfuncdeco
    def run_rmats_model_prepare(self):
        cmd = '{} {}'.format(self.program['python'], self.script['rmats_model_prepare'])
        cmd += ' -i {}'.format(self.option('event_file').path)
        cmd += ' -g {}'.format(self.option('gene_id'))
        cmd += ' -o {}'.format(self.file['event_file'])
        cmd_name = 'run_rmats_model_prepare'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        else:
            self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="33704402")
        with open(self.file['event_file'], 'r') as e, open(self.file['event_file_new'], 'w') as n:
            lines = e.readlines()
            first_line = lines[0]
            n.write(first_line)
            for line in lines[1:]:
                line_list = line.split('\t')
                chr = str(line_list[3])
                if chr not in self.chr_list:
                    chr_new = 'chr' + chr
                    if chr_new in self.chr_list:
                        line_list[3] = chr_new
                        n.write('\t'.join(line_list))
                else:
                    n.write(line)

    @toolfuncdeco
    def run_rmats2sashimiplot(self):
        cmd = '{} {}'.format(self.program['python'], self.script['rmats2sashimiplot'])
        cmd += ' -t {}'.format(self.option('event_type'))
        cmd += ' -e {}'.format(self.file['event_file_new'])
        cmd += ' --l1 {}'.format(self.option('l1'))
        cmd += ' --l2 {}'.format(self.option('l2'))
        cmd += ' -o {}'.format(self.work_dir)
        cmd += ' --b1 {}'.format(self.option('b1'))
        cmd += ' --b2 {}'.format(self.option('b2'))
        cmd += ' --exon_s {}'.format(self.option('exon_s'))
        cmd += ' --intron_s {}'.format(self.option('intron_s'))
        cmd_name = 'run_rmats2sashimiplot'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        else:
            self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="33704402")

    @toolfuncdeco
    def run_convert(self):
        for n, pdf in enumerate(glob.glob(os.path.join(self.dir['Sashimi_plot'], '*.pdf'))):
            cmd = '{}'.format(self.program['convert'])
            cmd += ' -background {}'.format(self.option('background'))
            cmd += ' -density {}'.format(self.option('density'))
            cmd += ' -quality {}'.format(self.option('quality'))
            cmd += ' -flatten {} {}'.format(pdf, '{}.png'.format(pdf[:-4]))
            cmd_name = 'run_convert_{}'.format(n)
            command = self.add_command(cmd_name, cmd)
            command.run()
            self.wait()
            if command.return_code == 0:
                self.logger.info("{} Finished successfully".format(cmd_name))
            else:
                self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="33704402")

    @toolfuncdeco
    def set_output(self):
        for basename in os.listdir(self.dir['Sashimi_plot']):
            source = os.path.join(self.dir['Sashimi_plot'], basename)
            link_name = os.path.join(self.output_dir, basename)
            if os.path.exists(link_name):
                os.remove(link_name)
            os.link(source, link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'rmats_model_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v3.structure.rmats_model',
            'instant': False,
            'options': {
                'gene_id': 'ENSG00000166579',
                'event_file': '/mnt/ilustre/users/sanger-dev/workspace/20190618/RmatsModel_tsg_33555_2602_4955/remote_input/result_dir/Con_12_vs_Vit_12/SE.MATS.JC.alter_id.txt',
                'event_type': 'SE',
                'l1': 'Vit1,Vit2',
                'l2': 'Con1,Con2',
                'b1': '/mnt/ilustre/users/sanger-dev/workspace/20190618/RmatsModel_tsg_33555_2602_4955/Download/output/Vit1.bam,/mnt/ilustre/users/sanger-dev/workspace/20190618/RmatsModel_tsg_33555_2602_4955/Download/output/Vit2.bam',
                'b2': '/mnt/ilustre/users/sanger-dev/workspace/20190618/RmatsModel_tsg_33555_2602_4955/Download/output/Con1.bam,/mnt/ilustre/users/sanger-dev/workspace/20190618/RmatsModel_tsg_33555_2602_4955/Download/output/Con2.bam'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
