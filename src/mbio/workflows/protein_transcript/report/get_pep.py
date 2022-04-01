# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

from biocluster.workflow import Workflow
from mbio.packages.ref_rna_v3.functions import workfuncdeco
import os
import shutil
import datetime
from collections import OrderedDict
import json
from biocluster.config import Config


class GetPepWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GetPepWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'project_type', 'type': 'string', 'default': 'ref_rna_v2'},
            {'name': 'task_id', 'type': 'string', 'default': None},
            {'name': 'dir_suffix', 'type': 'string', 'default': '12345'}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    @workfuncdeco
    def check_options(self):
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))
        if not self.option('project_type'):
            self.set_error("project_type必须提供")
        if not self.option("task_id"):
            self.set_error('必须提供task_id')

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    @workfuncdeco
    def run(self):
        pep_path = ''
        self.genome = Config().get_mongo_client(mtype='ref_rna_v2', db_version=1)[
                    Config().get_mongo_dbname(mtype='ref_rna_v2', db_version=1)]
        try:
            self.db = Config().get_mongo_client(mtype=self.option('project_type'), db_version=1)[
                Config().get_mongo_dbname(mtype=self.option('project_type'), db_version=1)]
            self.logger.info(self.db)
            if self.option('project_type') == 'whole_transcriptome':
                self.task_info =  self.db['task'].find_one({'task_id': self.option('task_id')})
            else:
                self.task_info = self.db['sg_task'].find_one({'task_id': self.option('task_id')})
            if not self.task_info:
                try:
                    self.db = \
                    Config().get_mongo_client(mtype=self.option('project_type'), db_version=0)[
                        Config().get_mongo_dbname(mtype=self.option('project_type'), db_version=0)]
                    self.logger.info(self.db)
                    if self.project_type == 'whole_transcriptome':
                        self.task_info = self.db['task'].find_one({'task_id': self.option('task_id')})
                    else:
                        self.task_info = self.db['sg_task'].find_one({'task_id': self.option('task_id')})
                    if not self.task_info:
                        try:
                            self.db = Config().get_mongo_client(mtype=self.option('project_type'), db_version=1,
                                                                task_id=self.option('task_id'))[
                                Config().get_mongo_dbname(mtype=self.option('project_type'), db_version=1,
                                                          task_id=self.option('task_id'))]
                            self.logger.info(self.db)
                            if self.project_type == 'whole_transcriptome':
                                self.task_info = self.db['task'].find_one({'task_id': self.option('task_id')})
                            else:
                                self.task_info = self.db['sg_task'].find_one({'task_id': self.option('task_id')})
                            if not self.task_info:
                                self.set_error("数据库中找不到记录")
                        except:
                            pass
                except:
                    pass
        except:
            self.set_error("数据库中找不到记录")
        if self.option('project_type') in ['ref_rna', 'ref_rna_v2', 'whole_transcriptome', 'medical_transcriptome']:
            self.start_listener()
            genome_id = self.task_info['genome_id']
            pep_path = os.path.join(Config().SOFTWARE_DIR, 'database/Genome_DB_finish',self.genome['sg_genome_db'].find_one({'genome_id': genome_id})['pep'])
            self.transfer(pep_path=pep_path)
        elif self.option('project_type') in ['denovo_rna', 'denovo_rna_v2']:
            pep_s3 = self.task_info['bedpath'].replace('.bed', '.pep')
            self.download_from_s3(pep_s3)
            super(GetPepWorkflow, self).run()
        elif self.option('project_type') == 'prok_rna':
            pep_s3 = self.task_info['rock_index'] + '/cds.faa'
            self.download_from_s3(pep_s3)
            super(GetPepWorkflow, self).run()
        elif self.option('project_type') in ['metagenomic']:
            ppm = self.db['geneset'].find_one({"task_id": self.option('task_id')})["ppm"]
            pep_s3 = ppm.split('gene_profile')[0] + 'uniGeneset/gene.uniGeneset.faa'
            self.download_from_s3(pep_s3)
            super(GetPepWorkflow, self).run()

    @workfuncdeco
    def download_from_s3(self, infile):
        self.download = self.add_tool('ref_rna_v3.download')
        options = {
            'ifile': infile
        }
        self.download.set_options(options)
        self.download.on('end', self.transfer)
        self.download.run()

    def transfer(self, pep_path=None):
        if not pep_path:
            pep_path = self.download.option('ofile').path
        self.logger.info("开始向线下服务器传递结果文件，请耐心等待")
        target_dir = "{}_{}".format(self.option("task_id"), self.option("dir_suffix"))
        self.offline_results = os.path.join(self.work_dir, target_dir)
        if os.path.isdir(self.offline_results):
            shutil.rmtree(self.offline_results)
        os.mkdir(self.offline_results)
        os.link(pep_path, os.path.join(self.offline_results, os.path.basename(pep_path)))
        cmd = "scp -r -i ~/.ssh/id_rsa {} dongmei.fu@10.2.4.236:/mnt/ilustre/users/yitong.feng/scripts/sanger_webapi/sanger_db_dir".format(
            self.offline_results)
        try:
            code = os.system(cmd)
            if code == 0:
                self.logger.info("命令{}执行成功！".format(cmd))
                os.system("touch {}ok".format(self.work_dir + "/"))
                cmd1 = "scp -i ~/.ssh/id_rsa {} dongmei.fu@10.2.4.236:/mnt/ilustre/users/yitong.feng/scripts/sanger_webapi/sanger_db_dir/{}".format(
                    self.work_dir + "/ok", target_dir)
                os.system(cmd1)
                self.logger.info("已创建ok文件！")
            else:
                self.logger.info("命令{}执行失败！".format(cmd))
                self.set_error("向线下服务器传递数据失败")
        except:
            self.set_error("向线下服务器传递数据失败")
        self.end()

    @workfuncdeco
    def end(self):
        super(GetPepWorkflow, self).end()