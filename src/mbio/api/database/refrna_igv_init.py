# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/5/22 14:55

import re, os, Bio, argparse, sys, fileinput

import re, os, Bio, argparse, sys, fileinput, urllib2

from pymongo import MongoClient
from bson.objectid import ObjectId
import types
from types import StringTypes
import re, subprocess
import json, time
import pandas as pd
import numpy as np
import datetime, os
from bson.son import SON
from collections import Counter
import glob
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
import gridfs
from mainapp.models.mongo.ref_rna import RefRna
'''
rmats model 导表函数
'''


class RefrnaIgvInit(Base):
    def __init__(self, bind_object):
        super(RefrnaIgvInit, self).__init__(bind_object)
        self._db_name = Config().MONGODB + '_ref_rna'
        
    
    
    @report_check
    def add_igv_init(self, name=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        task_info = RefRna().get_main_info(task_id, 'sg_task')
        reference_fasta_url = task_info["params"]["ref_fasta_path"]
        reference_fai_url = task_info["params"]["ref_fasta_fai_path"]
        anno_bed_url = task_info["params"]["ref_genes_bed_path"]
        anno_bed_index_url = task_info["params"]["ref_genes_bed_tbi_path"]
        locus = task_info["params"]["init_locus"]
        params = '''{
            locus: "%s",
            showKaryo: false,
            showIdeogram: false,
            showNavigation: true,
            showRuler: true,
            showCursorTrackingGuide: true,
            showCenterGuide: true,
            trackLabelsVisible: true,
            showControls: true,
            showCommandBar: true,
            showSequence: true,
            trackDefaults: {
                bam: {coverageThreshold: 0.2, coverageQualityWeight: true}
            },
            reference: {fastaURL: "%s", indexURL: "%s", indexed: true
            },
            tracks: [{name: "Genes", url: "%s", indexURL: "%s", displayMode: "EXPANDED", format: "bed", sourceType: "file"}]}''' % (
            locus, reference_fasta_url, reference_fai_url, anno_bed_url, anno_bed_index_url)
        
        data = {
            'name': name if name else 'igv_init_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'desc': 'igv初始化主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            # 'params':
            #     json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params,
            'task_id': task_id,
            'project_sn': project_sn
        }
        self.bind_object.logger.info('要插入的数据为：%s' % data)
        collection_obj = self.db['sg_igv_init']
        try:
            igv_init_id = collection_obj.insert_one(data).inserted_id
            self.bind_object.logger.info('导入igv_init主表成功')
        except Exception as e:
            raise Exception('导入igv_init主表失败:%s' % e)
        return igv_init_id
