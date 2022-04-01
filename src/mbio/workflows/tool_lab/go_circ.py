# -*- coding: utf-8 -*-
# __author__ = 'shicaiping,zoujiaxun'

import os
from biocluster.workflow import Workflow
from biocluster.config import Config
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError
import pandas as pd
import math
import time
from biocluster.file import getsize, exists
from biocluster.file import download


class GoCircWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GoCircWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'source', 'type': 'string', 'default': 'tool_lab'},
            {'name': 'diff_exp', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'go_file', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'term_num', 'type': 'int', 'default': 10},
            {'name': 'relate_id', 'type': 'string'},
            {'name': 'task_id', 'type': 'string'},
            {'name': 'project_task_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.go_circ")
        self.set_options(self._sheet.options())

    def run(self):
        if self.option('source') == 'project':
            self.diff_exp, self.go_path = self.check_file_path()
            self.diff_file = self.download_s3_file(self.diff_exp, 'diff_exp.txt')
            self.go_file = self.download_s3_file(self.go_path, 'go_file.txt')
        elif self.option('source') == 'ref':
            self.diff_file = self.option('diff_exp').path
            self.go_file = self.option('go_file').path

        if self.option('source') == 'tool_lab':
            self.go_file_new = self.option('go_file').path
        elif self.option('source') in ['ref', 'project']:
            diff_dict = dict()
            diff_exp = pd.read_table(self.diff_file, sep='\t', index_col=0)
            diff_exp.fillna('-', inplace=True)
            columns = diff_exp.columns
            for col in columns:
                if "log2fc" in col:
                    log2fc_col = col
                if col.startswith('fc'):
                    fc_col = col
            for i in diff_exp.index.tolist():
                diff_dict[i] = {'gene_name': diff_exp.loc[i]['gene_name'], 'log2fc': diff_exp.loc[i][log2fc_col], 'fc':  diff_exp.loc[i][fc_col]}
            go_file = pd.read_table(self.go_file, header=0, index_col=0, sep='\t')
            self.go_file_new = os.path.join(self.work_dir, 'go_file.txt')
            with open(self.go_file_new, 'w') as g:
                g.write('category' + '\t' + 'ID' + '\t' + 'term' + '\t' + 'count' + '\t' + 'genes' + '\t' + 'logFC' + '\t' + 'adj_pval' + '\t' + 'zscore' + '\n')
                for i in go_file.index.tolist():
                    go_type = go_file.loc[i]['go_type']
                    term = go_file.loc[i]['discription']
                    seq_list = go_file.loc[i]['seq_list'].split(';')
                    count_seq = len(seq_list)
                    up_down = list()
                    for c in seq_list:
                        if diff_dict[c]['log2fc'] > 0:
                            up_down.append(1)
                        else:
                            up_down.append(-1)
                    zscore = sum(up_down)/math.sqrt(count_seq)
                    for j in seq_list:
                        g.write(go_type + '\t' + i + '\t' + term + '\t' + str(count_seq) + '\t' + diff_dict[j]['gene_name']
                                + '\t' + str(diff_dict[j]['log2fc']) + '\t' + str(go_file.loc[i]['p_corrected']) + '\t' + str(zscore) + '\n')





        self.run_tool()
        super(GoCircWorkflow, self).run()

    def check_file_path(self):
        collection_name = 'tool_thurl'
        project_type = 'tool_lab'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        conn_upset = db[collection_name]
        status = 'start'
        count_time = 0
        while status == 'start':
            if count_time > 600:
                self.set_error('超过十分钟还没有结果文件生成，请检查是否生成文件时报错')
                break
            time.sleep(10)
            print 'sleep 10s'
            try:
                upset = conn_upset.find_one(
                    {'task_id': self.option('project_task_id'), 'relate_id': ObjectId(self.option('relate_id'))})
                status = upset['status']
            except:
                pass
            count_time += 10
        upset = conn_upset.find_one(
            {'task_id': self.option('project_task_id'), 'relate_id': ObjectId(self.option('relate_id'))})
        diff_path = upset['diff_path']
        enrich_path = upset['enrich_path']
        return diff_path, enrich_path

    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        if os.path.exists(path):
            to_path = path
        elif exists(path):
            download(path, to_path)
        else:
            self.set_error('file can not find %s', variables=(path,), code='13700502')
        return to_path

    def run_tool(self):
        opts = {
            'go_file': self.go_file_new,
            'term_num': self.option('term_num')
        }
        self.tool.set_options(opts)
        self.tool.on('end', self.set_db)
        self.tool.run()

    def set_db(self):

        """
         保存结果标准化数据到mongo数据库中
        """
        # picedit = self.api.api('tool_lab.picedit')
        # picedit.add_translation(self.option('main_id'), self.tool.option('translation_file').path)
        s3_path = '{}/{}'.format(self._sheet.output, 'go_circ.pdf')
        go_circ = self.api.api('tool_lab.go_circ')
        go_circ.add_go_circ(self.option('main_id'), s3_path)
        self.set_output()


    def set_output(self):
        for file in os.listdir(self.tool.output_dir):
            os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "GO富集圈图", 0],
        ])
        super(GoCircWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.exp_norm import ExpNormWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "exp_norm" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.exp_norm",
            "options": dict(
                exp_matrix="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/de_tools/known_seqs_count.matrix",
                convert_type="DESeq2",
                # float_num=4,
            )
        }
        wsheet = Sheet(data=data)
        wf =ExpNormWorkflow(wsheet)
        wf.sheet.id = 'exp_norm'
        wf.sheet.project_sn = 'exp_norm'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
