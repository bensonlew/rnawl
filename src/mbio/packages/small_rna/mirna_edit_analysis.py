# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import os
import sys
import json
from scipy.stats import binom_test
import numpy as np
import pandas as pd

parser = OptionParser(description='Conduct miRNA edit analysis by bowtie results')
parser.add_option('-p', '--hairpin', dest='hairpin' ,help='bowtie result of hairpin.fa against to genome')
parser.add_option('-q', '--sample', dest='sample', help='bowtie result of sample.fq against to genome')
parser.add_option('-s', '--site', dest='site', help='json file of site relationship between hairpin and mature')
parser.add_option('-n', '--name', dest='name', help='json file of name relationship between hairpin and mature')
parser.add_option('-o', '--output', dest='output', help='path of result file')
(opts, args) = parser.parse_args()

class MirnaEditAnalysis(object):
    '''
    Process miRNA bases editing anlysis by (1) and (2) and (3) and (4) and (5)
    (1) bowtie result of hairpin.fa against to genome
    (2) bowtie result of sample.fq against to genome
    (3) json file of site relationship between hairpin and mature
    (4) json file of name relationship between hairpin and mature
    (5) path of result file
    '''

    # ----- initialization and check -----------------------------------------------------------------------------------

    def __init__(self, hairpin_hit, sample_hit, mirna_site_json, mirna_name_json, output_file):
        self.phred_plus = 33
        self.output_file = output_file
        if self.check(hairpin_hit):
            self.hairpin_hit = hairpin_hit
        if self.check(sample_hit, detect_phred=True):
            self.sample_hit = sample_hit
        with open(mirna_site_json) as f:
            self.hairpin_mature_site = json.load(f)
        with open(mirna_name_json) as f:
            self.hairpin_relate_name = json.load(f)

    def check(self, hit, detect_phred=False):
        if not os.path.isfile(hit):
            raise Exception('Error: {} is not a file'.format(hit))
        if os.path.getsize(hit) == 0:
            print 'Warning: {} is empty, abord'.format(hit)
            self.result_list = list()
            self.generate_result_file()
            sys.exit(0)
        try:
            first_line = open(hit).readline()
            line_list = first_line.split('\t')
            if line_list[1] not in ['+', '-'] and not line_list[3].isdigit():
                raise Exception('Error: find format error in {}'.format(hit))
        except:
            raise Exception('Error: {} is not readable'.format(hit))
        if detect_phred:
            with open(hit) as f:
                for line in f:
                    quality = line.strip().split('\t')[5]
                    for q in quality:
                        if q in list('!"#$%&()*+,-./0123456789:;<=>?FI'):
                            self.phred_plus = 33
                            return True
                        elif q in list('KLMNOPQRSTUVWXYZ[\]^_`abcdefghi'):
                            self.phred_plus = 64
                            return True
        else:
            return True

    # ----- main flow --------------------------------------------------------------------------------------------------

    def main(self):
        self.process_hairpin_hit()
        self.process_sample_hit()
        self.binomial_analysis()
        self.adjust_p_value()
        self.match_hairpin_name()
        self.thin_to_ura()
        self.generate_result_file()

    # ----- sub process ------------------------------------------------------------------------------------------------

    def process_hairpin_hit(self):
        self.genome_location_base = dict()
        self.genome_location_hairpin = dict()
        self.hairpin_location_relative_position = dict()
        with open(self.hairpin_hit) as f:
            for line in f:
                self.line_to_genome_location_base(line)
                self.line_to_genome_location_hairpin(line)
                self.line_to_hairpin_location_relative_position(line)

    def process_sample_hit(self):
        self.genome_location_base_distribution = dict()
        self.genome_location_mismatch = dict()
        self.genome_location_mismatch_quality = dict()
        self.genome_location_mismatch_average_error_rate = dict()
        with open(self.sample_hit) as f:
            for line in f:
                self.line_to_genome_location_base_distribution(line)
                self.line_to_genome_location_mismatch(line)
                self.line_to_genome_location_mismatch_quality(line)
        self.make_genome_location_mismatch_average_error_rate()

    def binomial_analysis(self):
        self.result_list = list()
        for chromosome in self.genome_location_mismatch:
            for location in self.genome_location_mismatch[chromosome]:
                for hairpin_id in self.genome_location_hairpin[chromosome][location]:
                    if hairpin_id not in self.hairpin_mature_site:
                        continue
                    strand = self.genome_location_mismatch[chromosome][location]
                    w = self.genome_location_base[chromosome][location]
                    distribution_dict = self.genome_location_base_distribution[chromosome][location]
                    n = sum(distribution_dict.values())
                    p = self.genome_location_mismatch_average_error_rate[chromosome][location]

                    for e, x in distribution_dict.iteritems():
                        if e != w and x > 0:
                            p_value = binom_test(x, n, p, alternative='greater')
                            if strand == '+':
                                _w = w
                                _e = e
                            elif strand == '-':
                                _w = self.complementary(w)
                                _e = self.complementary(e)
                            position_in_hairpin = self.hairpin_location_relative_position[hairpin_id][location]
                            for mature_id in self.hairpin_mature_site[hairpin_id]:
                                if position_in_hairpin in self.hairpin_mature_site[hairpin_id][mature_id]:
                                    _mature_id = mature_id
                                    position_in_mature = self.hairpin_mature_site[hairpin_id][mature_id].index(position_in_hairpin)
                                    self.result_list.append({
                                        'mature_id': _mature_id,
                                        'hairpin_id': hairpin_id,
                                        'w': _w,
                                        'e': _e,
                                        'mp': position_in_mature,
                                        'hp': position_in_hairpin,
                                        'mis_num': x,
                                        'total_num': n,
                                        'p_value': p_value
                                    })

    def adjust_p_value(self):
        if not self.result_list:
            return
        self.p_value_dict = dict()
        for i, dct in enumerate(self.result_list):
            self.p_value_dict[str(i)] = dct['p_value']
        self.p_adjust_dict = self.benjaminiand_hochberg_method(self.p_value_dict)
        for i, dct in enumerate(self.result_list):
            self.result_list[i]['p_adjust'] = self.p_adjust_dict[str(i)]

    def match_hairpin_name(self):
        if not self.result_list:
            return
        for i, dct in enumerate(self.result_list):
            hairpin_id = self.result_list[i]['hairpin_id']
            _hairpin_id = self.hairpin_relate_name[hairpin_id]['harpin_id']
            self.result_list[i]['hairpin_id'] = _hairpin_id

    def thin_to_ura(self):
        if not self.result_list:
            return
        for i, dct in enumerate(self.result_list):
            if dct['w'] == 'T':
                self.result_list[i]['w'] = 'U'
            if dct['e'] == 'T':
                self.result_list[i]['e'] = 'U'

    def generate_result_file(self):
        columns = ['mature_id', 'hairpin_id', 'w', 'e', 'mp', 'hp', 'mis_num', 'total_num', 'p_value', 'p_adjust']
        if not self.result_list:
            with open(self.output_file, 'w') as f:
                f.write('{}\n'.format('\t'.join(columns)))
        else:
            df = pd.DataFrame(self.result_list, columns=columns)
            df.to_csv(self.output_file, sep='\t', index=None)

    # ----- fucntion ---------------------------------------------------------------------------------------------------

    @staticmethod
    def complementary(base):
        return 'TGCA'['ACGT'.index(base.upper())]

    def line_to_genome_location_base(self, line):
        line_list = line.strip().split('\t')
        hairpin_id = line_list[0]
        strand = line_list[1]
        chromosome = line_list[2]
        location = int(line_list[3])
        sequence = line_list[4]
        for base in sequence:
            if chromosome in self.genome_location_base:
                self.genome_location_base[chromosome][str(location)] = base
            else:
                self.genome_location_base[chromosome] = {str(location): base}
            location += 1

    def line_to_genome_location_hairpin(self, line):
        line_list = line.strip().split('\t')
        hairpin_id = line_list[0]
        strand = line_list[1]
        chromosome = line_list[2]
        location = int(line_list[3])
        sequence = line_list[4]
        for base in sequence:
            if chromosome in self.genome_location_hairpin:
                if str(location) in self.genome_location_hairpin[chromosome]:
                    self.genome_location_hairpin[chromosome][str(location)].append(hairpin_id)
                else:
                    self.genome_location_hairpin[chromosome][str(location)] = [hairpin_id]
            else:
                self.genome_location_hairpin[chromosome] = {str(location): [hairpin_id]}
            location += 1

    def line_to_hairpin_location_relative_position(self, line):
        line_list = line.strip().split('\t')
        hairpin_id = line_list[0]
        strand = line_list[1]
        chromosome = line_list[2]
        location = int(line_list[3])
        sequence = line_list[4]
        if strand == '+':
            for position, base in enumerate(sequence):
                if hairpin_id in self.hairpin_location_relative_position:
                    self.hairpin_location_relative_position[hairpin_id][str(location)] = position
                else:
                    self.hairpin_location_relative_position[hairpin_id] = {str(location): position}
                location += 1
        if strand == '-':
            location = location + len(sequence) - 1
            for position, base in enumerate(sequence):
                if hairpin_id in self.hairpin_location_relative_position:
                    self.hairpin_location_relative_position[hairpin_id][str(location)] = position
                else:
                    self.hairpin_location_relative_position[hairpin_id] = {str(location): position}
                location -= 1

    def line_to_genome_location_base_distribution(self, line):
        line_list = line.strip().split('\t')
        hairpin_id = line_list[0]
        strand = line_list[1]
        chromosome = line_list[2]
        location = int(line_list[3])
        sequence = line_list[4]
        if chromosome not in self.genome_location_base:
            return
        for base in sequence:
            if str(location) in self.genome_location_base[chromosome]:
                if chromosome in self.genome_location_base_distribution:
                    if str(location) in self.genome_location_base_distribution[chromosome]:
                        self.genome_location_base_distribution[chromosome][str(location)][base.upper()] += 1
                    else:
                        self.genome_location_base_distribution[chromosome][str(location)] = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
                        self.genome_location_base_distribution[chromosome][str(location)][base.upper()] += 1
                else:
                    self.genome_location_base_distribution[chromosome] = {str(location): {'A': 0, 'C': 0, 'G': 0, 'T': 0}}
                    self.genome_location_base_distribution[chromosome][str(location)][base.upper()] += 1
            location += 1

    def line_to_genome_location_mismatch(self, line):
        line_list = line.strip().split('\t')
        if len(line_list) == 8:
            hairpin_id = line_list[0]
            strand = line_list[1]
            chromosome = line_list[2]
            location = int(line_list[3])
            sequence = line_list[4]
            quality = line_list[5]
            mismatch = line_list[7]
        else:
            return
        if chromosome not in self.genome_location_base:
            return
        if strand == '+':
            mismatch_location = location + int(mismatch.split(':')[0])
        elif strand == '-':
            mismatch_location = location + len(sequence) - 1 -int(mismatch.split(':')[0])
        if str(mismatch_location) in self.genome_location_base[chromosome]:
            if chromosome in self.genome_location_mismatch:
                self.genome_location_mismatch[chromosome][str(mismatch_location)] = strand
            else:
                self.genome_location_mismatch[chromosome] = {str(mismatch_location): strand}

    def line_to_genome_location_mismatch_quality(self, line):
        line_list = line.strip().split('\t')
        if len(line_list) == 8:
            hairpin_id = line_list[0]
            strand = line_list[1]
            chromosome = line_list[2]
            location = int(line_list[3])
            sequence = line_list[4]
            quality = line_list[5]
            mismatch = line_list[7]
        else:
            return
        if chromosome not in self.genome_location_base:
            return
        if strand == '+':
            mismatch_location = location + int(mismatch.split(':')[0])
            mismatch_index = int(mismatch.split(':')[0])
        elif strand == '-':
            mismatch_location = location + len(sequence) - 1 - int(mismatch.split(':')[0])
            mismatch_index = len(quality) - 1 - int(mismatch.split(':')[0])
        if str(mismatch_location) in self.genome_location_base[chromosome]:
            if chromosome in self.genome_location_mismatch_quality:
                if str(mismatch_location) in self.genome_location_mismatch_quality[chromosome]:
                    self.genome_location_mismatch_quality[chromosome][str(mismatch_location)].append(quality[mismatch_index])
                else:
                    self.genome_location_mismatch_quality[chromosome][str(mismatch_location)] = [quality[mismatch_index]]
            else:
                self.genome_location_mismatch_quality[chromosome] = {str(mismatch_location): [quality[mismatch_index]]}

    @staticmethod
    def quality_to_error_rate(quality, plus=33):
        return pow(10, ((ord(quality)-plus) / -10.0))

    def make_genome_location_mismatch_average_error_rate(self):
        for chromosome in self.genome_location_mismatch_quality:
            self.genome_location_mismatch_average_error_rate[chromosome] = dict()
        for chromosome in self.genome_location_mismatch_quality:
            for location in self.genome_location_mismatch_quality[chromosome]:
                error_rate_list = list()
                for quality in self.genome_location_mismatch_quality[chromosome][location]:
                    error_rate_list.append(self.quality_to_error_rate(quality, self.phred_plus))
                average_error_rate = np.mean(error_rate_list)
                if average_error_rate > 1.0:
                    average_error_rate = 1.0
                self.genome_location_mismatch_average_error_rate[chromosome][location] = average_error_rate

    @staticmethod
    def benjaminiand_hochberg_method(p_value_dict):
        sort_lst = list()
        for k, v in p_value_dict.iteritems():
            if sort_lst == list():
                sort_lst.append([k, v])
                continue
            for n, i in enumerate(sort_lst):
                if v <= i[1]:
                    sort_lst.insert(n, [k, v])
                    break
            else:
                sort_lst.append([k, v])
        rank_lst = list()
        tmp_p = sort_lst[0][1]
        rank = 1
        for n, i in enumerate(sort_lst):
            if i[1] != tmp_p:
                rank_lst.append([i[0], i[1], n + 1])
                tmp_p = i[1]
                rank = n + 1
            else:
                rank_lst.append([i[0], i[1], rank])
        corrected_lst = list()
        for i in rank_lst:
            corrected_lst.append([i[0], i[1] * len(rank_lst) / (i[2])])
        ret = dict()
        for i in corrected_lst:
            if i[1] > 1:
                ret[i[0]] = 1.0
            else:
                ret[i[0]] = i[1]
        return ret

if __name__ == '__main__':
    if opts.hairpin and opts.sample and opts.site and opts.name and opts.output:
        inst = MirnaEditAnalysis(
            hairpin_hit=opts.hairpin,
            sample_hit=opts.sample,
            mirna_site_json=opts.site,
            mirna_name_json=opts.name,
            output_file=opts.output
        )
        inst.main()
    else:
        parser.print_help()
