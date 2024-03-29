# -*- coding: utf-8 -*-
# __author__ = 'shicaiping,zoujiaxun'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError
from Bio import SeqIO
from biocluster.file import getsize, exists
from biocluster.file import download
from biocluster.config import Config
import time

class GenomeCircWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenomeCircWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'genome_type', 'type': 'string', 'default': 'finish'},
            {'name': 'genome_fa', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'genome_struction', 'type': 'infile', 'format':'medical_transcriptome.common'},
            {'name': 'update_info', 'type': 'string'},
            {"name": 'source', 'type': 'string', 'default': 'tool_lab'},  # ['tool_lab', 'project'],
            {"name": 'relate_id', 'type': 'string'},
            {'name': 'task_id', 'type': 'string'},
            {'name': 'project_task_id', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.genome_circ")
        self.set_options(self._sheet.options())

    def run(self):
        if self.option('source') == 'project':
            self.fa_path, self.file_path = self.check_file_path()
            self.genome_path = self.download_s3_file(self.fa_path, 'genome.fa')
            self.struction_path = self.download_s3_file(self.file_path, 'table.txt')
        self.run_tool()
        super(GenomeCircWorkflow, self).run()

    def run_tool(self):
        if self.option('source') == 'tool_lab':
            genome_fa = self.option('genome_fa').path
            genome_struction = self.option('genome_struction').path
        elif self.option('source') == 'project':
            genome_fa = self.genome_path
            genome_struction = self.struction_path
        opts = {
            'genome_type': self.option('genome_type'),
            'genome_fa': genome_fa,
            'genome_struction': genome_struction
        }
        self.tool.set_options(opts)
        self.tool.on('end', self.set_db)
        self.tool.run()

    def set_db(self):

        """
         保存结果标准化数据到mongo数据库中
        """
        genome_circ_file = os.path.join(self.tool.output_dir, 'genome_circ.txt')
        genome_circ = self.api.api('tool_lab.genome_circ')
        if self.option('source') == 'tool_lab':
            genome_fa = self.option('genome_fa').path
            genome_struction = self.option('genome_struction').path
        elif self.option('source') == 'project':
            genome_fa = self.genome_path
            genome_struction = self.struction_path
        genome_circ.add_genome_circ(self.option('main_id'), genome_circ_file, self.option('genome_type'), genome_fa, genome_struction)
        self.set_output()

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
        fa_path = upset['path']
        file_path = upset['path1']
        return fa_path,file_path

    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        if exists(path):
            download(path, to_path)
        else:
            self.set_error('file can not find %s', variables=(path), code = '13700502')
        return to_path


    def set_output(self):
        self.end()

    def end(self):
        super(GenomeCircWorkflow, self).end()


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
