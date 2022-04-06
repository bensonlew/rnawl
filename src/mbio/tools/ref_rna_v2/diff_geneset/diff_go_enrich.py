# -*- coding: utf-8 -*-
# __author__ = 'shenghe,qinjincheng'

from biocluster.agent import Agent
from mbio.packages.ref_rna_v2.functions import toolfuncdeco, runcmd
from biocluster.tool import Tool
import os
import shutil
import unittest
from mbio.packages.rna.annot_config import AnnotConfig

class DiffGoEnrichAgent(Agent):
    '''
    last_modify: 2019.06.12
    '''
    def __init__(self, parent):
        super(DiffGoEnrichAgent, self).__init__(parent)
        options = [
            {'name': 'diff_list', 'type': 'infile', 'format': 'ref_rna_v2.gene_list'},
            {'name': 'go_list', 'type': 'infile', 'format': 'ref_rna_v2.go_list'},
            {'name': 'go_version', 'type': 'string', 'default': '20200628'},
            {'name': 'alpha', 'type': 'float', 'default': 0.05},
            {'name': 'pval', 'type': 'float', 'default': 0.05},
            {'name': 'method', 'type': 'string', 'default': 'bh'},
            {'name': 'result', 'type': 'outfile', 'format': 'ref_rna_v2.common'}
        ]
        self.add_option(options)
        self.step.add_steps('go_enrich')
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.go_enrich.start()
        self.step.update()

    def step_finish(self):
        self.step.go_enrich.finish()
        self.step.update()

    @toolfuncdeco
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        else:
            method = self.option('method').lower()
            if method == 'bonferroni':
                self.option('method', 'sm_bonferroni')
            elif method == 'bh':
                self.option('method', 'fdr_bh')
            elif method == 'by':
                self.option('method', 'fdr_by')
            elif method == 'holm':
                self.option('method', 'sm_holm')
        self.logger.info('{} = {}'.format("diff_list", self.option("diff_list").prop["path"]))
        self.logger.info('{} = {}'.format("go_list", self.option("go_list").prop["path"]))

    def set_resource(self):
        self._cpu = 1
        self._memory = '8G'

    @toolfuncdeco
    def end(self):
        super(DiffGoEnrichAgent, self).end()

class DiffGoEnrichTool(Tool):
    def __init__(self, config):
        super(DiffGoEnrichTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python'
        }
        self.script = {
            'find_enrichment': os.path.join(
                self.config.SOFTWARE_DIR, 'bioinfo/annotation/goatools-0.6.5-shenghe/scripts/find_enrichment.py'
            ),
            'go_enrich_stats': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/geneset/go_enrich_stats.py')
        }
        self.file = {
            'obo': os.path.join(self.config.SOFTWARE_DIR, 'database/Annotation/other2019/go.obo'),
            'outfile': os.path.join(self.work_dir, 'outfile.txt'),
            'study': os.path.join(self.work_dir, 'diff.list'),
            'population': os.path.join(self.work_dir, 'all.list'),
            'association': os.path.join(self.work_dir, 'background.txt'),
            'result': os.path.join(self.output_dir, 'go_enrich_geneset_list_gene.xls')
        }

        if self.option("go_version") == "2018":
            self.file['obo'] = AnnotConfig().get_file_dict(db="go", version=self.option("go_version"))['go-basic']
        else:
            self.file['obo'] = AnnotConfig().get_file_dict(db="go", version=self.option("go_version"))['go']

        # self.goatools_path = '/bioinfo/annotation/goatools-0.6.5-shenghe'
        # self.go_enrich_path = self.goatools_path + '/scripts/find_enrichment.py'
        # self.obo = self.config.SOFTWARE_DIR + '/database/GO/go-basic.obo'
        # self.set_environ(PYTHONPATH=self.config.SOFTWARE_DIR + self.goatools_path)
        # self.python_path = 'miniconda2/bin/python'
        # self.out_enrich_fp = self.output_dir + '/go_enrich_' + os.path.splitext(os.path.basename(self.option('diff_list').path))[0] + '.xls'
        # self.out_go_graph = self.output_dir + '/go_lineage'
        # self.image_magick_path = self.config.SOFTWARE_DIR + '/program/ImageMagick/bin/'
        # self.out_adjust_graph = self.output_dir + '/adjust_lineage'
        # self.class_code_dict = {}
        # self.set_environ(FONTCONFIG_PATH=self.config.SOFTWARE_DIR + '/library/fontconfig-2.13.1/etc/fonts')

    @toolfuncdeco
    def run(self):
        super(DiffGoEnrichTool, self).run()
        self.pre_find_enrichment()
        self.run_find_enrichment()
        self.run_go_enrich_stats()
        self.set_output()
        self.end()

    @toolfuncdeco
    def pre_find_enrichment(self):
        shutil.copy(self.option('go_list').path, self.file['association'])
        poplst = {line.split('\t')[0] for line in open(self.file['association'])}
        open(self.file['population'], 'w').writelines('{}\n'.format(s) for s in poplst)
        open(self.file['study'], 'w').writelines(
            line for line in open(self.option('diff_list').path) if line.strip() in poplst
        )

    @toolfuncdeco
    def run_find_enrichment(self):
        cmd = '{} {}'.format(self.program['python'], self.script['find_enrichment'])
        cmd += ' --alpha {}'.format(self.option('alpha'))
        cmd += ' --pval {}'.format(self.option('pval'))
        cmd += ' --indent'
        cmd += ' --obo {}'.format(self.file['obo'])
        cmd += ' --outfile {}'.format(self.file['outfile'])
        cmd += ' --method {}'.format(self.option('method'))
        cmd += ' {study} {population} {association}'.format(**self.file)
        cmd_name = 'run_find_enrichment'
        runcmd(self, cmd_name, cmd)

    @toolfuncdeco
    def run_go_enrich_stats(self):
        cmd = '{} {}'.format(self.program['python'], self.script['go_enrich_stats'])
        cmd += ' -i {}'.format(self.file['outfile'])
        cmd += ' -o {}'.format(self.file['result'])
        cmd_name = 'run_go_enrich_stats'
        runcmd(self, cmd_name, cmd)

    # def run_enrich(self):
    #     back_ground = self.option('go_list').prop['path']
    #     # if self.get_geneset_type(self.option('diff_list').prop['path']) == 'known':
    #     #     back_ground = self.choose_known_background()
    #     cmd0 = 'less {}| cut -f1 > {}/all.list'.format(back_ground, self.work_dir)
    #     os.system(cmd0)
    #     new_file_name = self.check_list()  # edited by shijin 除去背景中不存在的基因
    #     cmd = self.python_path + ' ' + self.config.SOFTWARE_DIR + self.go_enrich_path + ' '
    #     cmd = cmd + new_file_name + ' ' + self.work_dir + '/all.list' + ' ' + back_ground
    #     cmd = cmd + ' --pval ' + self.option('pval') + ' --indent' + ' --method ' + self.option('method') + ' --outfile ' + self.out_enrich_fp
    #     cmd = cmd + ' --obo ' + self.obo
    #     command = self.add_command('go_enrich', cmd)
    #     command.run()
    #     self.wait()
    #     if command.return_code == 0:
    #         self.run_draw_go_graph()
    #     else:
    #         self.set_error('goatools计算错误', code = '33706103')

    # def check_list(self):
    #     '''
    #     去除diff_list中没有注释信息的数据
    #     new_file_name为在work_dir中生成的新diff_list文件的绝对路径
    #     :return:
    #     '''
    #     file1 = self.option('diff_list').prop['path']
    #     file2 = self.work_dir + '/all.list'
    #     f1 = open(file1, 'r')
    #     f2 = open(file2, 'r')
    #     lst_1 = f1.readlines()
    #     lst_2 = f2.readlines()
    #     f1.close()
    #     f2.close()
    #     new_file_name = self.work_dir + '/' + os.path.basename(file1)
    #     with open(new_file_name, 'w') as fw:
    #         for item in lst_1:
    #             if item in lst_2:
    #                 fw.write(item)
    #     return new_file_name

    # def get_geneset_type(self, diff_file):
    #     '''
    #     查看基因集是已知基因集还是已知基因+新基因 刘彬旭
    #     '''
    #     gene_set_f = open(diff_file, 'r')
    #     for line in gene_set_f.readlines():
    #         if line.startswith('MSTR') or line.startswith('TCON') or line.startswith('XLOC'):
    #             gene_set_f.close()
    #             return 'all'
    #     gene_set_f.close()
    #     return 'known'

    # def choose_known_background(self):
    #     '''
    #     如果基因集中仅含有已知基因，修改背景注释为已知基因注释 刘彬旭
    #     '''
    #     go_list_file = self.option('go_list').prop['path']
    #     new_go_list = self.work_dir + '/known_go.list'
    #     with open(new_go_list, 'w') as ngfw, open(go_list_file, 'r') as go_list_f:
    #         for line in go_list_f.readlines():
    #             if line.startswith('MSTR') or line.startswith('TCON') or line.startswith('XLOC'):
    #                 pass
    #             else:
    #                 ngfw.write(line)
    #         return new_go_list

    # def run_draw_go_graph(self):
    #     try:
    #         self.logger.info('run_draw_go_graph')
    #         go_pvalue,go_padjust = self.get_go_pvalue_dict()
    #         self.logger.info('rrrrrrrrrrrrrrrrun_draw_go_graph')
    #         self.logger.info(go_pvalue)
    #         self.logger.info('run_draw_go_graphhhhhhhhhhhhhhhhh')
    #         if go_pvalue:
    #             # cmd = self.image_magick_path + 'convert {} {}'.format(self.out_go_graph + '.png', self.out_go_graph + '.pdf')
    #             draw_GO(go_pvalue, out=self.out_go_graph, obo=self.obo)
    #             # subprocess.check_output(cmd, shell=True)
    #         if go_padjust:
    #             # cmd = self.image_magick_path + 'convert {} {}'.format(self.out_adjust_graph + '.png', self.out_adjust_graph + '.pdf')
    #             draw_GO(go_padjust, out=self.out_adjust_graph, obo=self.obo)
    #             # subprocess.check_output(cmd, shell=True)
    #         self.end()
    #     except Exception:
    #         self.set_error('绘图发生错误:\n%s', variables = (traceback.format_exc()), code = '33706104')

    # def get_go_pvalue_dict(self):
    #     go2pvalue = {}
    #     go2padjust = {}
    #     with open(self.out_enrich_fp) as f:
    #         f.readline()
    #         for line in f:
    #             line_sp = line.split('\t')
    #             p_bonferroni = float(line_sp[6])
    #             padjust = float(line_sp[9])
    #             go2padjust[line_sp[0]] = padjust
    #             go2pvalue[line_sp[0]] = p_bonferroni
    #     tar = sorted(go2pvalue.items(), key=lambda e:e[1], reverse=True)
    #     tar_adjust = sorted(go2padjust.items(), key=lambda e:e[1], reverse=True)
    #     new_go2padjust = dict(tar_adjust[-10:])
    #     new_go2pvalue = dict(tar[-10:])
    #     self.logger.info(new_go2pvalue)
    #     return new_go2pvalue, new_go2padjust

    @toolfuncdeco
    def set_output(self):
        self.option('result').set_path(self.file['result'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'go_enrich_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v2.geneset.go_enrich',
            'instant': False,
            'options': {
                'diff_list': '/mnt/ilustre/users/sanger-dev/workspace/20190611/GenesetEnrich_tsg_33555_9930_8830/geneset_list_gene.list',
                'go_list': '/mnt/ilustre/users/sanger-dev/workspace/20190611/GenesetEnrich_tsg_33555_9930_8830/GO.list',
                'method': 'BH'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
