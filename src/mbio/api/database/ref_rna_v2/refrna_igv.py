# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/5/11 16:44

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
from mbio.files.sequence.fasta import FastaFile


class RefrnaIgv(base):
    def __init__(self, bind_object):
        super(RefrnaIgv, self).__init__(bind_object)
        self._db_name = Config().MONGODB + '_ref_rna'
    
    @report_check
    def add_init_option(self, params=None, major=True, ref_gtf=None, name=None, gtf_gz=None, gtf_tbi=None, fai=None,
                        ref_fa=None, host=None):
        
        if not (ref_gtf or gtf_gz):
            raise Exception('没有设置ref gtf路径')
        if not ref_fa:
            raise Exception('没有设置ref gtf路径')
        if not host:
            raise Exception('没有设置host路径')
        if not fai:
            fa_obj = FastaFile(ref_fa)
            fa_obj.check()
            fai = fa_obj.index()
        fai_url = self.get_url(path=fai, host=host)
        fa_url = self.ger_url(path=ref_fa, host=host)
        self.get_contig_len_dic(ref_fa)
