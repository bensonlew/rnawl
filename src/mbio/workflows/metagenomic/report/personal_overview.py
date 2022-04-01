# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'

from biocluster.workflow import Workflow
from biocluster.api.file.remote import RemoteFileManager
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from mainapp.models.mongo.metagenomic import Metagenomic
import os
import re
import gevent
import shutil
import json
import types
from mainapp.libs.param_pack import group_detail_sort
import copy
from biocluster.core.exceptions import OptionError
from mbio.packages.metagenomic.id_convert import id2name


anno_col_info = {
    'go': ['anno_go', 'GO_Origin'],  # [col_name, table_name]
    'phi': ['anno_phi', 'PHI_Origin'],
    'mvirdb': ['anno_mvir', 'AnnoMvirdb_Origin'],
    'tcdb': ['anno_tcdb', 'AnnoTcdb_Origin'],
    'qs': ['anno_qs', 'QS_Origin'],
    'pfam': ['anno_pfam', 'PFAM_Origin'],
    'sec': ['anno_sec', 'AnnoSec_Origin'],
    'sec_de': ['anno_sec', 'AnnoSec_Origin_Deunclassified'],
    'sec_lca': ['anno_sec', 'AnnoSec_Origin_LCA'],
    't3ss': ['anno_ttss', 'TTSS_Origin'],
    't3ss_de': ['anno_ttss', 'TTSS_Origin_Deunclassified'],
    't3ss_lca': ['anno_ttss', 'TTSS_Origin_LCA'],
    'probio': ['anno_probio', 'Probio_Origin'],
    'probio_de': ['anno_probio', 'Probio_Origin_Deunclassified'],
    'probio_lca': ['anno_probio', 'Probio_Origin_LCA'],
    'p450': ['anno_cyps', 'CYPS_Origin'],
    'nr_de': ['anno_nr', 'NR_Origin_Deunclassified'],
    'nr_lca': ['anno_nr', 'NR_Origin_LCA'],
}


class PersonalOverviewWorkflow(Workflow):
    """
    宏基因组个性化注释
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PersonalOverviewWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "anno_list", "type": "string", "default": ""},  # 个性化注释类型
            {"name": "task_id", "type": "string"},  # 所属项目的task_id
            {"name": "update_info", "type": "string", "default": ""},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.common_api = self.api.api('metagenomic.common_api')
        self.col = self.common_api.db['personal_anno']
        self.overview = self.add_tool('annotation.database_merge')
        self.overview_api = self.api.api('metagenomic.mg_anno_overview')
        self.running = set()

    def run(self):
        self.logger.info("Start Run Workflow")
        self.add_event('wait_anno')
        self.on('wait_anno', self.run_overview)
        self.stop_timeout_check()
        self.on('start', self.wait_anno)
        self.overview.on('end', self.set_db)
        super(PersonalOverviewWorkflow, self).run()

    def run_overview(self):
        geneset = self.get_one('geneset', 'GENESET_Origin')
        self.geneset_id = geneset['_id']
        self.ori_overview = self.get_one('anno_overview', 'Overview_Origin')
        gene_len = geneset['gene_list_length']
        if 'new_sum' in self.ori_overview:
            pre_sum = self.ori_overview['new_sum']
        else:
            pre_sum = gene_len.split('gene_profile')[0] + 'anno_overview.xls'
        types = []
        files = []
        for anno in self.running:
            if anno.endswith('_de') or anno.endswith('_lca'):
                continue
            anno_main = self.get_one(*anno_col_info[anno])
            self.logger.info("anno_main:{} ; anno_col: {}".format(anno_main, anno_col_info[anno]))
            if anno == 'sec':
                for j in ['sec_Gram_neg', 'sec_Gram_pos', 'sec_Euk']:
                    types.append(j)
                for j in ['signalp_Gram-_SignalP.txt',
                          'signalp_Gram+_SignalP.txt',
                          'signalp_Euk_SignalP.txt']:
                    files.append(anno_main['anno_file'] + '/' + j)
            elif anno in ['probio_lca', 'probio_de_unclassified']:
                continue
            else:
                types.append(anno)
                files.append(anno_main['anno_file'])
        opts = {
            'gene_length_table': self.download_s3(geneset['gene_list_length']),
            'other': ','.join(types),
            'otherf': self.download_files(files),
            'pre_sum': self.download_s3(pre_sum)
        }
        self.overview.set_options(opts)
        self.overview.run()

    def get_one(self, col, name):
        col_info = self.common_api.find_one(col,
                                            {'task_id': self.option('task_id'),
                                             'name': name})
        return col_info

    def ov_delete(self, filters, main=True):
        if main:
            col = self.common_api.db['anno_overview']
            col.delete_one(filters)
        else:
            col = self.common_api.db['anno_overview_detail']
            col.delete_many(filters)

    def wait_anno(self):
        self.step.update()
        self.stop_timeout_check()
        done = set()
        while True:
            self.logger.info('等待完成后，进行后续注释')
            p_ov = self.col.find_one({'task_id': self.option('task_id'),
                                      'type': 'overview'})
            if not p_ov:
                self.set_error('表personal_anno中无overview信息，请确定正常发起运行')
            self.running = set(p_ov['running_list'].split(','))
            for one in self.running:
                one_info = self.col.find_one({'task_id': self.option('task_id'),
                                              'type': one})
                if one_info['status'] == 'failed':
                    self.set_error(one + '注释运行出错，结束总览合并运行')
                if one_info['status'] == 'end':
                    done.add(one)
            self.logger.info('等待{}注释完成后，进行overview步骤'.format(
                self.running - done
            ))
            if done == self.running:
                break
            gevent.sleep(5)
        self.fire('wait_anno')

    def download_files(self, file_list):
        '''
        下载远程注释文件
        '''
        file_path = []
        if not os.path.exists(os.path.join(self.work_dir, 'remote_input')):
            os.mkdir(os.path.join(self.work_dir, 'remote_input'))
        for f in file_list:
            file_path.append(self.download_s3(f))
        return ','.join(file_path)

    def download_s3(self, file_path):
        remote_file = RemoteFileManager(file_path)
        remote_file.download(os.path.join(self.work_dir, 'remote_input'))
        return remote_file.local_path

    def set_db(self):
        """
        保存结果output，导mongo数据库
        """
        # 设置结果文件
        file_path = os.path.join(self.overview.output_dir,
                                 'gene_overview_anno.xls')
        self.link(file_path)
        overview_id = self.ori_overview['_id']
        if 'last_sum' in self.ori_overview:
            last_sum = self.ori_overview['last_sum']
        else:
            last_sum = []
        task_id = self.option('task_id')
        # 导入新的overview主表及详情表
        new_overview_id = self.overview_api.add_anno_overview(self.geneset_id,
                                                              task_id=task_id)
        anno_path = os.path.join(self.sheet.output, 'gene_overview_anno.xls')
        last_sum.append(anno_path)
        self.overview_api.add_anno_overview_detail(new_overview_id, file_path)
        self.logger.info('总览表导完')
        overview_col = self.common_api.db['anno_overview']
        overview_col.update_one({'_id': ObjectId(new_overview_id)},
                                {'$set': {'new_sum': anno_path,
                                          'last_sum': last_sum}})
        for anno_type in self.running:
            self.col.update_one({'task_id': self.option('task_id'),
                                 'type': anno_type},
                                {'$set': {'in_overview': 1}})
        # 删除之前的主表和对应详情表
        self.logger.info("开始删除原详情表")
        self.logger.info("开始删除原详情表overview_detail表数据！")
        self.ov_delete({'anno_overview_id': overview_id}, main=False)
        self.ov_delete({'_id': overview_id}, main=True)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ['gene_overview_anno.xls', 'xls', '个性化注释总览表', 0, "120346"]
        ])
        super(PersonalOverviewWorkflow, self).end()
