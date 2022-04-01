# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import glob
import os
import re
import shutil

import numpy as np
import pandas as pd
from biocluster.config import Config
from bson.objectid import ObjectId


class Upload(object):
    def __init__(self,project_type="whole_transcriptome"):
        self.database = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]

# self.database = Config().get_mongo_client(mtype='whole_transcriptome')[Config().get_mongo_dbname('whole_transcriptome')]


    def set_assembly(self,idir, odir):
        if os.path.isdir(odir):
            shutil.rmtree(odir)
        else:
            os.mkdir(odir)
        # transcript length distribution table
        for file in glob.glob(os.path.join(idir, 'Statistics/trans_count_stat_*.txt')):
            step = re.search(r'trans_count_stat_(\d+?).txt', file).group(1)
            df0 = pd.read_table(file)
            df0.index.name = 'Length'
            df0 = df0.rename({'all_transcripts': 'Number'}, axis=1)
            df0.to_csv(os.path.join(odir, 'length_distribution.{}.xls'.format(step)), sep='\t')
        # new transcript type statistics table
        df1 = pd.DataFrame([
            ('=', 'Complete match of intron chain'),
            ('c', 'Contained'),
            ('e',
             'Single exon transfrag overlapping a reference exon and at least 10 bp of a reference intron,indiction a possible pre-mRNA fragment'),
            ('i', 'Transfrag falling entirely within a reference intron'),
            (
                'j',
                'Potentially novel isoform (fragment):at least one splice junction is shared with a reference transcript'),
            ('o', 'Generic exonic overlap with a referfence transcript'),
            ('p', 'Possible polymerase run-on fragment(within 2Kbases of a reference transcript'),
            ('s',
             'An intron of the transfrag overlaps a reference intron on the opposite strand(likely due to read mapping errors'),
            ('u', 'Unknown,intergenic transcript'),
            ('x', 'Exonic overlap with reference on the opposite strand'),
            ('r',
             'Repeat. Currently determined by looking at the soft-masked reference sequence and applied to transcripts where at least 50% of the bases are lower case'),
            ('.', 'Tracking file only,indicates multiple classifications')], columns=['Class code', 'Description'])
        df2 = pd.read_table(os.path.join(idir, 'Statistics/code_num.txt'), names=['Class code', 'Number'], usecols=[0, 2])
        df3 = pd.merge(df1, df2, how='outer')
        df3 = df3.fillna(0)
        df3['Number'] = df3['Number'].astype(np.int64)
        df3.to_csv(os.path.join(odir, 'classcode_statistics.xls'), sep='\t', index=False)
        # sequence files
        sdir = os.path.join(odir, 'Sequence')
        if os.path.isdir(sdir):
            shutil.rmtree(sdir)
        os.mkdir(sdir)
        os.link(os.path.join(idir, 'NewTranscripts/all_transcripts.fa'), os.path.join(sdir, 'all_transcripts.fa'))
        os.link(os.path.join(idir, 'NewTranscripts/new_transcripts.fa'), os.path.join(sdir, 'new_transcripts.fa'))
        os.link(os.path.join(idir, 'NewTranscripts/ref_and_new.gtf'), os.path.join(sdir, 'all_transcripts.gtf'))
        os.link(os.path.join(idir, 'NewTranscripts/new_transcripts.gtf'), os.path.join(sdir, 'new_transcripts.gtf'))
        os.link(os.path.join(idir, 'NewTranscripts/trans2gene'), os.path.join(sdir, 'trans2gene.txt'))


    def set_rmats(self,idir, odir, task_id):
        if os.path.isdir(odir):
            shutil.rmtree(odir)
        os.mkdir(odir)
        # alternative splicing event statistics table
        os.link(os.path.join(idir, 'sample.event.count.JC.txt'), os.path.join(odir, 'event_type.JC.xls'))
        os.link(os.path.join(idir, 'sample.event.count.JCEC.txt'), os.path.join(odir, 'event_type.JCEC.xls'))
        # alternative splicing event detail table
        gid2des_dct,gid2name_dct = self.get_gid2des_dct(task_id)
        for document in self.database['splicing_rmats'].find({'task_id': task_id}):
            splicing_id = document['main_id']
            s1, s2 = document['compare_plan'].split('|')
            output_dir = os.path.join(odir, '{}_vs_{}'.format(s1, s2))
            # modify 20210830 调整cmp_order顺序
            self.export_rmats_detail(splicing_id, gid2des_dct, gid2name_dct,s2=s2, s1=s1, output_dir=output_dir)
            # self.export_rmats_detail(splicing_id, gid2des_dct, s2=s1, s1=s2, output_dir=output_dir)
        # alternative splicing intra-group event statistics table and pattern statistics table
        for document in self.database['splicing_rmats_stats'].find({'task_id': task_id}):
            stat_id = document['main_id']
            print 'stat_id -> ({})'.format(stat_id)
            print 'group -> ({})'.format(document['group'])
            for k, v in document['group'].items():
                if v == 's1':
                    s1 = k
                if v == 's2':
                    s2 = k
            else:
                print 's1 -> ({})'.format(s1)
                print 's2 -> ({})'.format(s2)
            output_dir = os.path.join(odir, '{}_vs_{}'.format(s1, s2))
            self.export_rmats_diff_stats(stat_id, output_dir)
            self.export_rmats_psi(stat_id, output_dir)


    def get_gid2des_dct(self,task_id):
        query_id = self.database['exp'].find_one({'task_id': task_id})['main_id']
        lst = list({(document['gene_id'], document.get('description', '')) for document in
                    self.database['exp_detail'].find({'exp_id': query_id})})
        lst2 = list({(document['gene_id'], document.get('gene_name', '')) for document in
                    self.database['exp_detail'].find({'exp_id': query_id})})
        return dict(lst),dict(lst2)


    def export_rmats_detail(self,splicing_id, gid2des_dct,gid2name_dct, s2, s1, output_dir):
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        dct = {k: {'JC': dict(), 'JCEC': dict()} for k in ('SE', 'MXE', 'A3SS', 'A5SS', 'RI')}
        for d in self.database['splicing_rmats_detail'].find({'splicing_id': ObjectId(splicing_id)}):
            eid = d['event_id']
            # gname = d['gname'] if d['gname'] != 'nan' else '-'
            gname = gid2name_dct[d['gid']]
            for jt in ('JC', 'JCEC'):
                dct[d['type']][jt][eid] = {'AS ID': d['event_id'], 'Gene ID': d['gid'], 'Gene name': gname,
                                           'Gene description': gid2des_dct[d['gid']],
                                           'Novel AS': d['novel_as'], 'Chr': d['chr'], 'Strand': d['strand']}
                if d['type'] == 'SE':
                    dct['SE'][jt][eid].update({
                        'InclusionTranscripts': d['InclusionTranscripts'],
                        'SkippingTranscripts': d['SkippingTranscripts'],
                        'ExonStart': d['es'],
                        'ExonEnd': d['ee'],
                        'UpstreamES': d['up_es'],
                        'UpstreamEE': d['up_ee'],
                        'DownstreamES': d['down_es'],
                        'DownstreamEE': d['down_ee']
                    })
                elif d['type'] == 'MXE':
                    dct['MXE'][jt][eid].update({
                        '1stExonTranscripts': d['1stExonTranscripts'],
                        '2ndExonTranscripts': d['2ndExonTranscripts'],
                        '1stExonStart': d['firstes'],
                        '1stExonEnd': d['firstee'],
                        '2ndExonStart': d['secondes'],
                        '2ndExonEnd': d['secondee']
                    })
                elif d['type'] == 'A3SS':
                    dct['A3SS'][jt][eid].update({
                        'LongExonTranscripts': d['LongExonTranscripts'],
                        'ShortExonTranscripts': d['ShortExonTranscripts'],
                        'LongExonStart': d['les'],
                        'LongExonEnd': d['lee'],
                        'ShortES': d['ses'],
                        'ShortEE': d['see'],
                        'FlankingES': d['fes'],
                        'FlankingEE': d['fee']
                    })
                elif d['type'] == 'A5SS':
                    dct['A5SS'][jt][eid].update({
                        'LongExonTranscripts': d['LongExonTranscripts'],
                        'ShortExonTranscripts': d['ShortExonTranscripts'],
                        'LongExonStart': d['les'],
                        'LongExonEnd': d['lee'],
                        'ShortES': d['ses'],
                        'ShortEE': d['see'],
                        'FlankingES': d['fes'],
                        'FlankingEE': d['fee']
                    })
                elif d['type'] == 'RI':
                    dct['RI'][jt][eid].update({
                        'RetainTranscripts': d['RetainTranscripts'],
                        'AbandonTranscripts': d['AbandonTranscripts'],
                        'RiExonStart': d['ries'],
                        'RiExonEnd': d['riee']
                    })
            dct[d['type']]['JC'][eid].update({
                'Diff significant': d['diff_jc'],
                'IncLevelDiff ({}/{},ΔPSI)'.format(s2, s1): d['inc_diff_jc'],
                'P Value JunctionCountOnly': d['pvalue_jc'],
                'FDR JunctionCountOnly': d['fdr_jc'],
                'IJC {}'.format(s2): d['ijc_s1'],
                'SJC {}'.format(s2): d['sjc_s1'],
                'IJC {}'.format(s1): d['ijc_s2'],
                'SJC {}'.format(s1): d['sjc_s2'],
                'IncFormLen': d['inclen_jc'],
                'SkipFormLen': d['skiplen_jc'],
                'IncLevel1 (PSI {})'.format(s2): d['inc1_jc'],
                'IncLevel2 (PSI {})'.format(s1): d['inc2_jc'],
                'Average IncLevel1 (Average PSI {})'.format(s2): d['aver_inc1_jc'],
                'Average IncLevel2 (Average PSI {})'.format(s1): d['aver_inc2_jc'],
                'Increase Inclusion {}'.format(s2): d['upinc_s1_jc'],
                'Increase Exclusion {}'.format(s2): d['upexc_s1_jc'],
                'Increase Inclusion {}'.format(s1): d['upinc_s2_jc'],
                'Increase Exclusion {}'.format(s1): d['upexc_s2_jc']
            })
            dct[d['type']]['JCEC'][eid].update({
                'Diff significant': d['diff_all'],
                'IncLevelDiff ({}/{},ΔPSI)'.format(s2, s1): d['inc_diff_all'],
                'P Value ReadsOnTargetAndJunctionCounts': d['pvalue_all'],
                'FDR ReadsOnTargetAndJunctionCounts': d['fdr_all'],
                'IC {}'.format(s2): d['ic_s1'],
                'SC {}'.format(s2): d['sc_s1'],
                'IC {}'.format(s1): d['ic_s2'],
                'SC {}'.format(s1): d['sc_s2'],
                'IncFormLen': d['inclen_all'],
                'SkipFormLen': d['skiplen_all'],
                'IncLevel1 (PSI {})'.format(s2): d['inc1_all'],
                'IncLevel2 (PSI {})'.format(s1): d['inc2_all'],
                'Average IncLevel1 (Average PSI {})'.format(s2): d['aver_inc1_all'],
                'Average IncLevel2 (Average PSI {})'.format(s1): d['aver_inc2_all'],
                'Increase Inclusion {}'.format(s2): d['upinc_s1_all'],
                'Increase Exclusion {}'.format(s2): d['upexc_s1_all'],
                'Increase Inclusion {}'.format(s1): d['upinc_s2_all'],
                'Increase Exclusion {}'.format(s1): d['upexc_s2_all']
            })

        def get_columns(event_type, junction_type):
            lst = ['AS ID', 'Gene ID', 'Gene name', 'Gene description', 'Novel AS', 'Chr', 'Strand',
                   'Diff significant', 'IncLevelDiff ({}/{},ΔPSI)'.format(s2, s1)]
            lst.extend({
                           'JC': ['P Value JunctionCountOnly', 'FDR JunctionCountOnly'],
                           'JCEC': ['P Value ReadsOnTargetAndJunctionCounts', 'FDR ReadsOnTargetAndJunctionCounts']
                       }[junction_type])
            lst.extend({
                           'SE': ['InclusionTranscripts', 'SkippingTranscripts', 'ExonStart', 'ExonEnd', 'UpstreamES',
                                  'UpstreamEE', 'DownstreamES', 'DownstreamEE'],
                           'MXE': ['1stExonTranscripts', '2ndExonTranscripts', '1stExonStart', '1stExonEnd', '2ndExonStart',
                                   '2ndExonEnd'],
                           'A3SS': ['LongExonTranscripts', 'ShortExonTranscripts', 'LongExonStart', 'LongExonEnd',
                                    'ShortES', 'ShortEE', 'FlankingES', 'FlankingEE'],
                           'A5SS': ['LongExonTranscripts', 'ShortExonTranscripts', 'LongExonStart', 'LongExonEnd',
                                    'ShortES', 'ShortEE', 'FlankingES', 'FlankingEE'],
                           'RI': ['RetainTranscripts', 'AbandonTranscripts', 'RiExonStart', 'RiExonEnd']
                       }[event_type])
            lst.extend({
                           'JC': ['IJC {}'.format(s2), 'SJC {}'.format(s2), 'IJC {}'.format(s1), 'SJC {}'.format(s1)],
                           'JCEC': ['IC {}'.format(s2), 'SC {}'.format(s2), 'IC {}'.format(s1), 'SC {}'.format(s1)]
                       }[junction_type])
            lst.extend([
                'IncFormLen', 'SkipFormLen', 'IncLevel1 (PSI {})'.format(s2), 'IncLevel2 (PSI {})'.format(s1),
                'Average IncLevel1 (Average PSI {})'.format(s2), 'Average IncLevel2 (Average PSI {})'.format(s1),
                'Increase Inclusion {}'.format(s2), 'Increase Exclusion {}'.format(s2),
                'Increase Inclusion {}'.format(s1), 'Increase Exclusion {}'.format(s1)
            ])
            return lst

        if not os.path.isdir(os.path.join(output_dir, 'JC')):
            os.mkdir(os.path.join(output_dir, 'JC'))
        if not os.path.isdir(os.path.join(output_dir, 'JCEC')):
            os.mkdir(os.path.join(output_dir, 'JCEC'))
        for et, jt2dcts in dct.items():
            print 'event_type -> ({})'.format(et)
            for jt, dcts in jt2dcts.items():
                print 'junction_type -> ({})'.format(jt)
                df = pd.DataFrame(dcts.values())
                df = df.reindex(get_columns(et, jt), axis=1)
                df.to_csv(os.path.join(output_dir, '{}/{}.detail.xls'.format(jt, et)), sep='\t', index=False)
        else:
            print 'succeed in exporting files to {}'.format(output_dir)


    def export_rmats_diff_stats(self,stat_id, output_dir):
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        document = self.database['splicing_rmats_diff_stats'].find_one({'stat_id': ObjectId(stat_id)})
        index = ['JunctionCountOnly(JC)', 'ReadsOnTargetAndJunctionCounts(JCEC)', 'JC&JCEC', 'JC|JCEC']
        df = pd.DataFrame(document['diff_stats'], index=index).T
        df = df.rename({i: i.upper() for i in df.index})
        df = df.reindex(['SE', 'MXE', 'A3SS', 'A5SS', 'RI', 'TOTAL'])
        df.index.name = 'AS type'
        df.to_csv(os.path.join(output_dir, 'diff_event_stats.xls'), sep='\t')
        print 'succeed in exporting {}'.format(os.path.join(output_dir, 'diff_event_stats.txt'))


    def export_rmats_psi(self,stat_id, output_dir):
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        document = self.database['splicing_rmats_psi'].find_one({'stat_id': ObjectId(stat_id)})
        index = ['Exclusion (Increase exclusion in case，ΔPSI<0)',
                 'Inclusion (Increase inclusion in case, ΔPSI>0)',
                 'Total events']
        df = pd.DataFrame(document['s1_jc'], index=index).T
        df = df.rename({i: i.upper() for i in df.index})
        df = df.reindex(['SE', 'MXE', 'A3SS', 'A5SS', 'RI', 'TOTAL'])
        df.index.name = 'AS type'
        jc_ofile = os.path.join(output_dir, 'diff_pattern_stats.JC.xls')
        df.to_csv(jc_ofile, sep='\t')
        print 'succeed in exporting {}'.format(jc_ofile)
        df = pd.DataFrame(document['s1_all'], index=index).T
        df = df.rename({i: i.upper() for i in df.index})
        df = df.reindex(['SE', 'MXE', 'A3SS', 'A5SS', 'RI', 'TOTAL'])
        df.index.name = 'AS type'
        jcec_ofile = os.path.join(output_dir, 'diff_pattern_stats.JCEC.xls')
        df.to_csv(jcec_ofile, sep='\t')
        print 'succeed in exporting {}'.format(jcec_ofile)
