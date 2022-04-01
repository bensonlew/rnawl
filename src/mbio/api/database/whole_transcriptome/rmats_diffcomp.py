# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mbio.api.database.whole_transcriptome.api_base import ApiBase
from biocluster.api.database.base import report_check
import pandas as pd
from bson.objectid import ObjectId
import unittest

class RmatsDiffcomp(ApiBase):
    '''
    last_modify: 2019.06.14
    '''
    def __init__(self, bind_object):
        super(RmatsDiffcomp, self).__init__(bind_object)
        self._e2t = {'SE': ('InclusionTranscripts', 'SkippingTranscripts'),
                     'RI': ('RetainTranscripts', 'AbandonTranscripts'),
                     'MXE': ('1stExonTranscripts', '2ndExonTranscripts'),
                     'A5SS': ('LongExonTranscripts', 'ShortExonTranscripts'),
                     'A3SS': ('LongExonTranscripts', 'ShortExonTranscripts')}

    @report_check
    def add_rmats_diffcomp_detail(self, diffcomp_txt, main_id):
        documents = list()
        df = pd.read_table(diffcomp_txt)
        df['rmats_diffcomp_id'] = ObjectId(main_id)
        for d in df.to_dict('r'):
            for k, v in d.items():
                if k.endswith('_JC') or k.endswith('_JCEC'):
                    try:
                        d[k] = float(v)
                    except:
                        pass
            else:
                k1, k2 = self._e2t[d['Type']]
                d[k1] = d['Transcripts_1']
                d[k2] = d['Transcripts_2']
                d.pop('Transcripts_1')
                d.pop('Transcripts_2')
                documents.append(d)
        self.create_db_table('splicing_rmats_diffcomp_detail', documents)
        insert_dict = {
            'main_id': ObjectId(main_id),
            'status': 'end'
        }
        self.update_db_record('splicing_rmats_diffcomp', main_id, insert_dict=insert_dict)
