# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.workflow import Workflow
from mbio.packages.ref_rna_v2.functions import workfuncdeco
import shutil
import pickle
import pandas as pd
import os
import unittest
import glob

class MasigproWorkflow(Workflow):
    '''
    last_modify: 2019.07.19
    '''
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MasigproWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'matrix', 'type': 'infile', 'format': 'ref_rna_v2.matrix'},
            {'name': 'target_gene_file', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'design', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'cluster', 'type': 'int', 'default': 9},
            {'name': 'method', 'type': 'string', 'default': 'hclust'},  # ['hclust', 'kmeans', 'Mclust']
            {'name': 'exp_type', 'type': 'string', 'default': 'TPM'},  # ['TPM', 'FPKM']
            {'name': 'result', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'main_id', 'type': 'string', 'default': None},
            {'name': 'update_info', 'type': 'string', 'default': None}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.target_genes = None

    @workfuncdeco
    def check_options(self):
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))
        if not self.option("target_gene_file").is_set:
            self.logger.info("没有指定的目标基因文件,使用表达量表的全部基因进行分析")


    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    @workfuncdeco
    def run(self):
        self.run_prepare_target_genes()
        self.run_masigpro()
        super(MasigproWorkflow, self).run()

    @workfuncdeco
    def run_prepare_target_genes(self):
        exp_df = pd.read_table(self.option("matrix").prop["path"], index_col=0)
        exp_df.index.name = "seq_id"
        exp_df.to_csv(self.option("matrix").prop["path"],sep="\t")
        if self.option("target_gene_file").is_set:
            self.target_genes = self.option("target_gene_file").prop["path"]
        else:
            exp_df = pd.read_table(self.option("matrix").prop["path"],index_col=0)
            target_genes = list(exp_df.index)
            with open(os.path.join(self.work_dir,"target_genes"),"w") as f:
                f.write("\n".join(target_genes))
            self.target_genes = os.path.join(self.work_dir,"target_genes")

    @workfuncdeco
    def run_masigpro(self):
        self.step.add_steps('masigpro')
        self.masigpro = self.add_tool('tool_lab.timeseries.masigpro')
        options = {
            'matrix': self.option('matrix'),
            'geneset': self.target_genes,
            'design': self.option('design'),
            'cluster': self.option('cluster'),
            'method': self.option('method'),
            'exp_type': self.option('exp_type')
        }
        self.masigpro.set_options(options)
        self.masigpro.on('start', self.set_step, {'start': self.step.masigpro})
        self.masigpro.on('end', self.set_step, {'end': self.step.masigpro})
        self.masigpro.on('end', self.set_output)
        self.masigpro.run()

    def set_output(self):
        shutil.rmtree(self.output_dir)
        self.modify_output_dir(self.masigpro.output_dir,self.output_dir)
        # shutil.copytree(self.masigpro.output_dir, self.output_dir)
        # os.rename(os.path.join(self.output_dir, 'result.tsv'), os.path.join(self.output_dir, 'result.xls'))
        # self.option('result').set_path(os.path.join(self.output_dir, 'result.xls'))

    def modify_output_dir(self,tool_out_dir,output_dir):
        shutil.copytree(tool_out_dir, output_dir)
        res_fs = glob.glob( output_dir + "/*/*result*")
        for res in res_fs:
            os.rename(res,os.path.join(os.path.dirname(res),"result.xls"))
        pkls = glob.glob( output_dir + "/*/*.pkl*")
        for pkl in pkls:
            seq_data = pickle.load(open(pkl))
            t = pd.DataFrame(seq_data)
            t["seq_ids"] = t["seq_ids"].apply(lambda x: ",".join(x))
            t= t[["time","seq_ids"]]
            name = os.path.basename(pkl)
            new_name = "heatmap" + os.path.splitext(name)[0] + ".gene_details.xls"
            t.to_csv(os.path.join(os.path.dirname(pkl),new_name),sep="\t",index=False)
            os.remove(pkl)


        self.end()
        # self.set_db()

    @workfuncdeco
    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [r'.', '', '时序差异分析结果目录',0,"211565"]
        ])
        result_dir.add_regexp_rules([
            [r'.', 'xls', '时序差异分析详情信息', 0],
            [r'.*/result.xls', 'xls', '时序差异分析详情表', 0],
            [r'.*/heatmap.*.xls', 'xls', '时序差异热图对应基因详情表', 0],
            [r'.*/heatmap.*.pdf', 'pdf', '时序差异热图pdf', 0],
            [r'.*/heatmap.*.png', 'png', '时序差异热图png', 0],
            [r'.*/heatmap.*.svg', 'svg', '时序差异热图svg', 0],
            [r'.*/groups.*.pdf', 'pdf', '不同处理(或同一处理)在不同时间点的表达模式图pdf', 0],
            [r'.*/groups.*.png', 'png', '不同处理(或同一处理)在不同时间点的表达模式图png', 0],
            [r'.*/groups.*.svg', 'svg', '不同处理(或同一处理)在不同时间点的表达模式图svg', 0],
            [r'.*/profile.*.pdf', 'pdf', '聚类指定基因在所有样本中的表达模式图pdf', 0],
            [r'.*/profile.*.png', 'png', '聚类指定基因在所有样本中的表达模式图png', 0],
            [r'.*/profile.*.svg', 'svg', '聚类指定基因在所有样本中的表达模式图svg', 0],
        ])
        super(MasigproWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''
    def test_kmeans(self):
        import random
        from mbio.workflows.tool_lab.masigpro import MasigproWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'tool_lab.masigpro',
            'instant': False,
            'options': {
                'matrix': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/masigpro/vs.matrix.tsv',
                'design': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/masigpro/vs.design.tsv',
                'cluster': 4,
                'method': 'kmeans'
            }
        }
        wsheet = Sheet(data=data)
        wf = MasigproWorkflow(wsheet)
        wf.sheet.id = 'ref_rna_v2_upgrade'
        wf.sheet.project_sn = 'ref_rna_v2_upgrade'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_kmeans')])
    unittest.TextTestRunner(verbosity=2).run(suite)
