# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from mbio.packages.ref_rna_v2.functions import toolfuncdeco, runcmd
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest

class MasigproAgent(Agent):
    '''
    last_modify: 2019.06.03
    '''
    def __init__(self, parent):
        super(MasigproAgent, self).__init__(parent)
        options = [
            {'name': 'matrix', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'geneset', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'design', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'cluster', 'type': 'int', 'default': 9},
            {'name': 'method', 'type': 'string', 'default': 'hclust'}, # ['hclust', 'kmeans', 'Mclust']
            {'name': 'exp_type', 'type': 'string', 'default': 'TPM'},  # ['TPM', 'FPKM']
            {'name': 'result', 'type': 'outfile', 'format': 'medical_transcriptome.common'},
            {'name': 'is_rmbe', 'type': 'string', 'default': 'false'}
        ]
        self.add_option(options)
        self.step.add_steps('masigpro')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.masigpro.start()
        self.step.update()

    def stepfinish(self):
        self.step.masigpro.finish()
        self.step.update()

    @toolfuncdeco
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(os.path.getsize(self.option('matrix').path) / 1024 ** 3 + 8)

    @toolfuncdeco
    def end(self):
        super(MasigproAgent, self).end()

class MasigproTool(Tool):
    def __init__(self, config):
        super(MasigproTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python',
            'rscript': 'bioinfo/rna/miniconda2/bin/Rscript'
        }
        self.script = {
            'pretreat': os.path.join(self.config.PACKAGE_DIR, 'medical_transcriptome/timeseries/pretreat.py'),
            'masigpro': os.path.join(self.config.PACKAGE_DIR, 'medical_transcriptome/timeseries/masigpro.r'),
            'heatmap.py': os.path.join(self.config.PACKAGE_DIR, 'medical_transcriptome/timeseries/heatmap.py'),
            'heatmap.r': os.path.join(self.config.PACKAGE_DIR, 'medical_transcriptome/timeseries/heatmap.r')
        }
        self.file = {
            'matrix': os.path.join(self.work_dir, 'matrix.delim.txt'),
            'desgin': os.path.join(self.work_dir, 'design.delim.txt'),
            # 'result': os.path.join(self.output_dir, 'result.tsv')
        }

    @toolfuncdeco
    def run(self):
        super(MasigproTool, self).run()
        self.run_pretreat()
        self.run_masigpro()
        self.run_heatmap()
        # self.set_output()
        self.end()

    @toolfuncdeco
    def run_pretreat(self):
        cmd = '{} {}'.format(self.program['python'], self.script['pretreat'])
        cmd += ' -m {}'.format(self.option('matrix').path)
        if self.option('geneset').is_set:
            cmd += ' -g {}'.format(self.option('geneset').path)
        cmd += ' -d {}'.format(self.option('design').path)
        cmd += ' -o {}'.format(self.work_dir)
        cmd += ' -rmbe {}'.format(self.option('is_rmbe'))
        cmd_name = 'run_pretreat'
        runcmd(self, cmd_name, cmd)

    @toolfuncdeco
    def run_masigpro(self):
        cmd = '{} {}'.format(self.program['rscript'], self.script['masigpro'])
        cmd += ' -i {}'.format(self.file['matrix'])
        cmd += ' -d {}'.format(self.file['desgin'])
        cmd += ' -k {}'.format(self.option('cluster'))
        cmd += ' -m {}'.format(self.option('method'))
        cmd += ' -o {}'.format(self.output_dir)
        cmd_name = 'run_masigpro'

        # try:
        #     command = self.add_command(cmd_name, cmd)
        #     command.run()
        #     self.wait(command)
        # except:
        #     with open("run_masigpro.o") as err:
        #         if "equal to 2" in err.read():
        #             self.set_error("number of sig.genes greater than or equal to 2")
        #         elif "equal to num of cluster" in err.read():
        #             self.set_error("number of sig.genes greater than or equal to num of cluster")
        #         else:
        #             self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd))

        command = self.add_command(cmd_name, cmd,ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
             self.logger.info("{} Finished successfully".format(cmd_name))
        else:
             with open("run_masigpro.o") as err:
                 if "no significant genes" in  err.read():
                     self.set_error("no significant genes")
                 if "equal to 2" in err.read():
                     self.set_error("number of sig.genes should greater than or equal to 2")
                 elif "equal to num of cluster" in err.read():
                     self.set_error("number of sig.genes should greater than or equal to num of cluster")
                 else:
                      self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd))


    @toolfuncdeco
    def run_heatmap(self):
        os.makedirs(os.path.join(self.work_dir,"heat_inter"))
        compare_infos = os.listdir(self.output_dir)
        for n,compare in enumerate(compare_infos):
            cmd = '{} {}'.format(self.program['python'], self.script['heatmap.py'])
            cmd += ' -m {}'.format(self.option('matrix').path)
            cmd += ' -d {}'.format(self.option('design').path)
            cmd += ' -r {}'.format(os.path.join(self.output_dir, compare, "result.tsv"))
            cmd += ' --interpreter {}'.format(os.path.join(self.config.SOFTWARE_DIR, self.program['rscript']))
            cmd += ' --script {}'.format(self.script['heatmap.r'])
            cmd += ' -t {}'.format(self.option('exp_type'))
            cmd += ' -o {}'.format(os.path.join(self.output_dir, compare))
            cmd_name = 'run_heatmap'+str(n)
            command = self.add_command(cmd_name, cmd)
            command.run()
            self.wait(command)
            # runcmd(self, cmd_name, cmd)

    # @toolfuncdeco
    # def set_output(self):
        # self.option('result').set_path(self.file['result'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test_hclust(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'masigpro_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'medical_transcriptome.timeseries.masigpro',
            'instant': False,
            'options': {
                'matrix': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/masigpro/matrix.tsv',
                'geneset': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/masigpro/geneset.list',
                # 'design': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/masigpro/design.txt',
                'design' :'/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/masigpro/test.design',
                'cluster': 9,
                'method': 'hclust'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_kmeans(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'masigpro_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v2.timeseries.masigpro',
            'instant': False,
            'options': {
                'matrix': '/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/test0721/masigpro/matrix.tsv',
                'design': '/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/test0721/masigpro/tt.txt',
                # 'matrix' :"/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/masigro/data/matrix.tsv",
                # 'design' : '/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/masigro/data/designtest.txt',
                'geneset':'/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/test0721/masigpro/geneset.list',
                'cluster': 2,
                'method': 'hclust'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_mclust(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'masigpro_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v2.timeseries.masigpro',
            'instant': False,
            'options': {
                # 'matrix': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/masigpro/vs.matrix.tsv',
                # 'design': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/masigpro/vs.design.tsv',
                'matrix':'/mnt/lustre/users/sanger/workspace/20200721/Masigpro_sanger_266906_1393_2128/matrix.tsv',
                'design': '/mnt/lustre/users/sanger/workspace/tmp/tmp_s3/sanger_266906/masigpro/tt.txt',
                'cluster': 4,
                'method': 'Mclust'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_kmeans')])
    unittest.TextTestRunner(verbosity=2).run(suite)
