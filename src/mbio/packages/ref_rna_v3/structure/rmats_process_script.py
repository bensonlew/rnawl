# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import pickle
import re

import pandas as pd

from mbio.packages.ref_rna_v3.functions import pkgsfuncdeco


class RmatsToolkit(object):
    def __init__(self, gtf, lastod, b1, b2):
        if os.path.isfile(gtf):
            self.gtf = gtf
        if os.path.isdir(lastod):
            self.lastod = lastod
        if os.path.isfile(b1):
            self.b1conf = b1
        if os.path.isfile(b2):
            self.b2conf = b2
        self.process_func_dict = {
            'SE': self.process_se_row,
            'RI': self.process_ri_row,
            'MXE': self.process_mxe_row,
            'A5SS': self.process_a5ss_row,
            'A3SS': self.process_a3ss_row,
        }
        self.event_header_dict = {
            'SE': ['InclusionTranscripts', 'SkippingTranscripts'],
            'RI': ['RetainTranscripts', 'AbandonTranscripts'],
            'MXE': ['1stExonTranscripts', '2ndExonTranscripts'],
            'A5SS': ['LongExonTranscripts', 'ShortExonTranscripts'],
            'A3SS': ['LongExonTranscripts', 'ShortExonTranscripts']
        }
        self.g2t2e = dict()
        self.etdat = dict()
        self.jcdct = dict()
        self.jcecdct = dict()

    @pkgsfuncdeco
    def run(self):
        self.get_sample_event_coordinates()
        self.dump_sample_event_coordinates()
        self.get_exon_coordinates()
        self.get_event_transcripts()
        self.export_transcrips_data()

    @pkgsfuncdeco
    def get_exon_coordinates(self):
        g2t2e = dict()
        pg = re.compile(r'gene_id "(\S+)";')
        pt = re.compile(r'transcript_id "(\S+)";')
        for line in open(self.gtf):
            items = line.strip().split('\t')
            if len(items) == 9 and items[0][0] != '#' and items[2] == 'exon':
                es = int(items[3]) - 1
                ee = int(items[4])
                mg = re.search(pg, items[8])
                mt = re.search(pt, items[8])
                if mg and mt:
                    g = mg.group(1)
                    t = mt.group(1)
                    if g in g2t2e:
                        if t in g2t2e[g]:
                            g2t2e[g][t].append((es, ee))
                        else:
                            g2t2e[g][t] = [(es, ee)]
                    else:
                        g2t2e[g] = {t: [(es, ee)]}
        else:
            for g in g2t2e:
                for t in g2t2e[g]:
                    g2t2e[g][t] = sorted(g2t2e[g][t])
            else:
                self.g2t2e = g2t2e

    @pkgsfuncdeco
    def get_event_transcripts(self):
        all_event_txpt_info = dict()
        for event_type in self.process_func_dict:
            event_source_table = os.path.join(self.lastod, 'fromGTF.{}.alter_id.txt'.format(event_type))
            df = pd.read_table(event_source_table)
            process_func = self.process_func_dict[event_type]
            event_txpt_data = dict()
            for i, row in df.iterrows():
                event_id = row['ID']
                txpt_info = process_func(row)
                event_txpt_data.update({event_id: txpt_info})
            else:
                all_event_txpt_info.update({event_type: event_txpt_data})
        else:
            self.etdat = all_event_txpt_info

    @pkgsfuncdeco
    def export_transcrips_data(self):
        for event_type, event_txpt_data in self.etdat.items():
            df = pd.DataFrame([pd.Series(txpt_info, index=self.event_header_dict[event_type], name=event_id)
                               for event_id, txpt_info in event_txpt_data.items()])
            df.index.name = 'event_id'
            df.to_csv(os.path.join(self.lastod, '{}.transcripts.txt'.format(event_type)), sep='\t')

    def process_se_row(self, series):
        """
        ese ~ tuple for exon (start, end)
        start -> 0 base offset
        """
        gene_id = series['GeneID']
        s_ese = (int(series['exonStart_0base']), int(series['exonEnd']))
        u_ese = (int(series['upstreamES']), int(series['upstreamEE']))
        d_ese = (int(series['downstreamES']), int(series['downstreamEE']))
        inclusion_transcripts = list()
        skipping_transcripts = list()
        for transcript_id, ese_list in self.g2t2e[gene_id].items():
            if u_ese in ese_list and d_ese in ese_list:
                if s_ese in ese_list:
                    inclusion_transcripts.append(transcript_id)
                else:
                    skipping_transcripts.append(transcript_id)
        else:
            ret1 = ','.join(inclusion_transcripts) if inclusion_transcripts else '-'
            ret2 = ','.join(skipping_transcripts) if skipping_transcripts else '-'
            return ret1, ret2

    def process_ri_row(self, series):
        """
        ese ~ tuple for exon (start, end)
        start -> 0 base offset
        """
        gene_id = series['GeneID']
        r_ese = (int(series['riExonStart_0base']), int(series['riExonEnd']))
        u_ese = (int(series['upstreamES']), int(series['upstreamEE']))
        d_ese = (int(series['downstreamES']), int(series['downstreamEE']))
        retain_transcripts = list()
        abandon_transcripts = list()
        for transcript_id, ese_list in self.g2t2e[gene_id].items():
            if r_ese in ese_list:
                retain_transcripts.append(transcript_id)
            elif u_ese in ese_list and d_ese in ese_list:
                abandon_transcripts.append(transcript_id)
        else:
            ret1 = ','.join(retain_transcripts) if retain_transcripts else '-'
            ret2 = ','.join(abandon_transcripts) if abandon_transcripts else '-'
            return ret1, ret2

    def process_mxe_row(self, series):
        """
        ese ~ tuple for exon (start, end)
        start -> 0 base offset
        """
        gene_id = series['GeneID']
        m1_ese = (int(series['1stExonStart_0base']), int(series['1stExonEnd']))
        m2_ese = (int(series['2ndExonStart_0base']), int(series['2ndExonEnd']))
        u_ese = (int(series['upstreamES']), int(series['upstreamEE']))
        d_ese = (int(series['downstreamES']), int(series['downstreamEE']))
        first_exon_transcripts = list()
        second_exon_transcripts = list()
        for transcript_id, ese_list in self.g2t2e[gene_id].items():
            if u_ese in ese_list and d_ese in ese_list:
                if m1_ese in ese_list:
                    first_exon_transcripts.append(transcript_id)
                if m2_ese in ese_list:
                    second_exon_transcripts.append(transcript_id)
        else:
            ret1 = ','.join(first_exon_transcripts) if first_exon_transcripts else '-'
            ret2 = ','.join(second_exon_transcripts) if second_exon_transcripts else '-'
            return ret1, ret2

    def process_a5ss_row(self, series):
        """
        ese ~ tuple for exon (start, end)
        start -> 0 base offset
        """
        gene_id = series['GeneID']
        l_ese = (int(series['longExonStart_0base']), int(series['longExonEnd']))
        s_ese = (int(series['shortES']), int(series['shortEE']))
        f_ese = (int(series['flankingES']), int(series['flankingEE']))
        long_exon_transcripts = list()
        short_exon_transcripts = list()
        for transcript_id, ese_list in self.g2t2e[gene_id].items():
            if f_ese in ese_list:
                if l_ese in ese_list:
                    long_exon_transcripts.append(transcript_id)
                if s_ese in ese_list:
                    short_exon_transcripts.append(transcript_id)
        else:
            ret1 = ','.join(long_exon_transcripts) if long_exon_transcripts else '-'
            ret2 = ','.join(short_exon_transcripts) if short_exon_transcripts else '-'
            return ret1, ret2

    def process_a3ss_row(self, series):
        """
        ese ~ tuple for exon (start, end)
        start -> 0 base offset
        """
        gene_id = series['GeneID']
        l_ese = (int(series['longExonStart_0base']), int(series['longExonEnd']))
        s_ese = (int(series['shortES']), int(series['shortEE']))
        f_ese = (int(series['flankingES']), int(series['flankingEE']))
        long_exon_transcripts = list()
        short_exon_transcripts = list()
        for transcript_id, ese_list in self.g2t2e[gene_id].items():
            if f_ese in ese_list:
                if l_ese in ese_list:
                    long_exon_transcripts.append(transcript_id)
                if s_ese in ese_list:
                    short_exon_transcripts.append(transcript_id)
        else:
            ret1 = ','.join(long_exon_transcripts) if long_exon_transcripts else '-'
            ret2 = ','.join(short_exon_transcripts) if short_exon_transcripts else '-'
            return ret1, ret2

    @pkgsfuncdeco
    def get_sample_event_coordinates(self):
        test_conf = self.b1conf
        ctrl_conf = self.b2conf
        big_table = os.path.join(self.lastod, 'all_events_detail_big_table.txt')
        test_samples = [os.path.basename(p)[:-4] for p in open(test_conf).read().split(',')]
        ctrl_samples = [os.path.basename(p)[:-4] for p in open(ctrl_conf).read().split(',')]
        self.jcdct = self.initialize_sample2event2set_dict(test_samples + ctrl_samples)
        self.jcecdct = self.initialize_sample2event2set_dict(test_samples + ctrl_samples)
        df = pd.read_table(big_table)
        func = lambda x: int(float(x))
        for n, row in df.iterrows():
            ijc1s = map(func, str(row['IJC_SAMPLE_1']).split(',')) if not pd.isnull(row['IJC_SAMPLE_1']) \
                else [0] * len(test_samples)
            sjc1s = map(func, str(row['SJC_SAMPLE_1']).split(',')) if not pd.isnull(row['SJC_SAMPLE_1']) \
                else [0] * len(test_samples)
            ijc2s = map(func, str(row['IJC_SAMPLE_2']).split(',')) if not pd.isnull(row['IJC_SAMPLE_2']) \
                else [0] * len(ctrl_samples)
            sjc2s = map(func, str(row['SJC_SAMPLE_2']).split(',')) if not pd.isnull(row['SJC_SAMPLE_2']) \
                else [0] * len(ctrl_samples)
            ic1s = map(func, str(row['IC_SAMPLE_1']).split(',')) if not pd.isnull(row['IC_SAMPLE_1']) \
                else [0] * len(test_samples)
            sc1s = map(func, str(row['SC_SAMPLE_1']).split(',')) if not pd.isnull(row['SC_SAMPLE_1']) \
                else [0] * len(test_samples)
            ic2s = map(func, str(row['IC_SAMPLE_2']).split(',')) if not pd.isnull(row['IC_SAMPLE_2']) \
                else [0] * len(ctrl_samples)
            sc2s = map(func, str(row['IC_SAMPLE_2']).split(',')) if not pd.isnull(row['IC_SAMPLE_2']) \
                else [0] * len(ctrl_samples)
            for i, s in enumerate(test_samples):
                jc_flag, jcec_flag = 0, 0
                if ijc1s[i] + sjc1s[i]:
                    jc_flag = 1
                if ic1s[i] + sc1s[i]:
                    jcec_flag = 1
                self.add_coordinates(row, jc_flag, jcec_flag, s)
            for i, s in enumerate(ctrl_samples):
                jc_flag, jcec_flag = 0, 0
                if ijc2s[i] + sjc2s[i]:
                    jc_flag = 1
                if ic2s[i] + sc2s[i]:
                    jcec_flag = 1
                self.add_coordinates(row, jc_flag, jcec_flag, s)

    @staticmethod
    def initialize_sample2event2set_dict(samples):
        return {s: {'SE': set(), 'RI': set(), 'MXE': set(), 'A3SS': set(), 'A5SS': set()} for s in samples}

    def add_coordinates(self, series, jc_flag, jcec_flag, sample):
        if series['type'] == 'SE':
            if jc_flag:
                self.jcdct[sample]['SE'].add((series['chr'], (series['upstreamES'], series['upstreamEE'],
                                                              series['exonStart_0base'], series['exonEnd'],
                                                              series['downstreamES'], series['downstreamEE'])))
            if jcec_flag:
                self.jcecdct[sample]['SE'].add((series['chr'], (series['upstreamES'], series['upstreamEE'],
                                                                series['exonStart_0base'], series['exonEnd'],
                                                                series['downstreamES'], series['downstreamEE'])))
        elif series['type'] == 'RI':
            if jc_flag:
                self.jcdct[sample]['RI'].add((series['chr'], (series['upstreamES'], series['upstreamEE'],
                                                              series['riExonStart_0base'], series['riExonEnd'],
                                                              series['downstreamES'], series['downstreamEE'])))
            if jcec_flag:
                self.jcecdct[sample]['RI'].add((series['chr'], (series['upstreamES'], series['upstreamEE'],
                                                                series['riExonStart_0base'], series['riExonEnd'],
                                                                series['downstreamES'], series['downstreamEE'])))
        elif series['type'] == 'MXE':
            if jc_flag:
                self.jcdct[sample]['MXE'].add((series['chr'], (series['upstreamES'], series['upstreamEE'],
                                                               series['1stExonStart_0base'], series['1stExonEnd'],
                                                               series['2ndExonStart_0base'], series['2ndExonEnd'],
                                                               series['downstreamES'], series['downstreamEE'])))
            if jcec_flag:
                self.jcecdct[sample]['MXE'].add((series['chr'], (series['upstreamES'], series['upstreamEE'],
                                                                 series['1stExonStart_0base'], series['1stExonEnd'],
                                                                 series['2ndExonStart_0base'], series['2ndExonEnd'],
                                                                 series['downstreamES'], series['downstreamEE'])))
        elif series['type'] == 'A3SS':
            if jc_flag:
                self.jcdct[sample]['A3SS'].add((series['chr'], (series['longExonStart_0base'], series['longExonEnd'],
                                                                series['shortES'], series['shortEE'],
                                                                series['flankingES'], series['flankingEE'])))
            if jcec_flag:
                self.jcecdct[sample]['A3SS'].add((series['chr'], (series['longExonStart_0base'], series['longExonEnd'],
                                                                  series['shortES'], series['shortEE'],
                                                                  series['flankingES'], series['flankingEE'])))
        elif series['type'] == 'A5SS':
            if jc_flag:
                self.jcdct[sample]['A5SS'].add((series['chr'], (series['longExonStart_0base'], series['longExonEnd'],
                                                                series['shortES'], series['shortEE'],
                                                                series['flankingES'], series['flankingEE'])))
            if jcec_flag:
                self.jcecdct[sample]['A5SS'].add((series['chr'], (series['longExonStart_0base'], series['longExonEnd'],
                                                                  series['shortES'], series['shortEE'],
                                                                  series['flankingES'], series['flankingEE'])))

    @pkgsfuncdeco
    def dump_sample_event_coordinates(self):
        pickle.dump(self.jcdct, open(os.path.join(self.lastod, 'sample.event.coordinates.JC.pk'), 'w'))
        pickle.dump(self.jcecdct, open(os.path.join(self.lastod, 'sample.event.coordinates.JCEC.pk'), 'w'))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Script for digging data from rMATS last output')
    parser.add_argument('--gtf', action='store', required=True, dest='gtf',
                        help='an annotation of genes and transcripts in GTF format')
    parser.add_argument('--od', action='store', required=True, dest='od',
                        help='output folder of last step')
    parser.add_argument('--b1', action='store', required=True, dest='b1',
                        help='test group BAM configuration file')
    parser.add_argument('--b2', action='store', required=True, dest='b2',
                        help='ctrl group BAM configuration file')
    args = parser.parse_args()

    inst = RmatsToolkit(gtf=args.gtf, lastod=args.od, b1=args.b1, b2=args.b2)
    inst.run()
