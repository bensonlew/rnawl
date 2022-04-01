#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/3/13 11:23
@file    : gtf_filter.py
"""

import csv
import json
import os
import unittest
from itertools import chain

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
from mbio.packages.lnc_rna.lnc_identification.predict_common_tools import CommandFactory


class GtfFilterAgent(Agent):
    def __init__(self, parent):
        """接受参数说明：
            biotype_key: gtf或gff属性部分biotype对应键的名称: ensembl中: transcript_biotype,
                该参数可选, 如果没传参数则通过retained_types过滤 gtf第三列
            lnc_biotypes: lncRNA对应的biotype类型, 以逗号分隔
            chrom_mapping: chromosome列如果不是数字或者chr1格式, 则要转换[文件两列, gtf_chrom_field<TAB>chrom_id]
            retained_attrs: 如果是gff文件是要提供, 格式: json string of list
                '[["gene", "gene_id"], ["transcript_id", "transcript_id"], ["Dbxref", "Dbxref"]]'
        :param parent:
        """
        super(GtfFilterAgent, self).__init__(parent)
        # lnc_biotypes = 'sense_overlapping,antisense,retained_intron,sense_intronic,macro_lncRNA,' \
        #                'bidirectional_promoter_lncRNA,3prime_overlapping_ncRNA,lincRNA,non_coding'
        options = [
            {'name': 'gtf', 'type': 'infile', 'format': 'lnc_rna.comm_gtf'},
            {'name': 'biotype_key', 'type': 'string', 'default': ''},
            {'name': 'lnc_biotypes', 'type': 'string', 'default': ''},
            {'name': 'anno_types', 'type': 'string'},
            {'name': 'chrom_mapping', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'retained_attrs', 'type': 'string', 'default': ''},
            {'name': 'out_gtf', 'type': 'outfile', 'format': 'lnc_rna.common'},
            {'name': 'out_ids', 'type': 'outfile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps("gtf_filter")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.gtf_filter.start()
        self.step.update()

    def step_end(self):
        self.step.gtf_filter.finish()
        self.step.update()

    def check_options(self):
        if not self.option('gtf').is_set:
            raise OptionError('缺少gtf or gff文件')
        if not self.option('biotype_key'):
            if not self.option('retained_attrs'):
                raise Exception('when file is gff file, retained_attrs must be set')
            self.logger.warning('如果想通过"注释信息的类型"提取一定要传入相应的类型')
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(GtfFilterAgent, self).end()


class GtfFilterTool(Tool):
    def __init__(self, config):
        super(GtfFilterTool, self).__init__(config)
        self.__filter_func = None
        self.__convert_chrom = lambda key, line,: line.lstrip('chr')

    def _confirm_filter_func(self):
        key = self.option('biotype_key')
        anno_types = {i.strip() for i in self.option('anno_types').strip().split(',')}
        self.logger.debug(', '.join(anno_types))
        if key:
            self.logger.debug('use __biotype_filter function filter')
            biotypes = {i.strip() for i in self.option('lnc_biotypes').strip().split(',')}
            self.__filter_func = self.__biotype_filter(key=key,
                                                       biotypes=biotypes,
                                                       anno_types=anno_types)
        else:
            self.logger.debug('use __anno_type_filter function filter')
            if not self.option('retained_attrs'):
                raise Exception('when file is gff file, retained_attrs must be set')
            self.__filter_func = self.__anno_type_filter(anno_types=anno_types)
        if self.option('chrom_mapping').is_set:
            self.__convert_chrom = self.__chrom_convert(self.option('chrom_mapping'))

    def __chrom_convert(self, file_obj):
        """序列的编号转换, 比如cnbi中该位置为chr_accessions, 需要转换为染色体编号

        :param file_obj:
        :return:
        """
        cv_dict = None
        with file_obj.get_reader() as in_handler:
            cv_dict = {k: v for k, v in (line.strip().split('\t') for line in in_handler if line.strip)}
        chroms_set = {str(i) for i in chain(range(1, 23), ('X', 'Y', 'MT'))}

        def _convert(key, line):
            chrom = cv_dict.get(key)
            chr_id, other = line.split('\t', 1)
            if chrom is None:
                line = None
            else:
                chrom = chrom.lower().lstrip('chr').upper()
                if chrom in chroms_set:
                    line = chrom + '\t' + other
                else:
                    line = None
            return line

        return _convert

    def __biotype_filter(self, key, biotypes, anno_types):
        """此方法过滤GTF文件
            通过biotype过滤lncRNAs

        :param key: gtf文件中属性列中biotype对应的键名, 例如在ensembl中transcript_biotype
        :param biotypes: lncRNA对应的所有biotype类型
        :param anno_types: gtf文件第3列, 注释类型过滤
        :return:
        """
        self.logger.debug(' =====  USE biotype_filter function ===== ')
        if self.option('retained_attrs'):
            try:
                attrs_keys = json.loads(self.option('retained_attrs'))
            except json.decoder.JSONDecodeError:
                attrs_keys = []
            except ValueError:
                attrs_keys = []
            if len(attrs_keys) < 2:
                t_id_key = [i for i, k in attrs_keys if k == 'transcript_id']
                t_id_key = t_id_key[0] if len(t_id_key) == 1 else 'transcript_id'
                attrs_keys = None
            else:
                has_gene_name = False
                old_gene = None
                for o, n in attrs_keys:
                    if n == 'gene_name':
                        has_gene_name = True
                    if n == 'gene_id':
                        old_gene = o
                if has_gene_name is not True:
                    attrs_keys.append((old_gene, 'gene_name'))

                attr_demo = ' '.join('%s {%s};' % (k, i) for i, k in attrs_keys)
        else:
            attrs_keys = None
            t_id_key = 'transcript_id'
        transcript_dic = {}

        def _filter(split_items, line):
            attr_dict = split_items[8]
            bio_type = attr_dict.get(key)
            if bio_type and bio_type in biotypes and split_items[2].strip() in anno_types:
                t_id = attr_dict[t_id_key]
                # if 'ID' in attr_dict:
                #     attr_dict['ID'] = '-'.join(attr_dict['ID'].split('-')[1: -1])
                recored_type = None
                if split_items[2] == 'exon':
                    recored_type = 'exon'
                    if 'ID' in attr_dict:
                        attr_dict['ID'] = transcript_dic['trans_id']
                else:
                    recored_type = 'transcript'
                    if 'ID' in attr_dict:
                        transcript_dic['trans_id'] = attr_dict['ID'].replace('rna-', '')
                        attr_dict['ID'] = transcript_dic['trans_id']
                if attrs_keys is not None:
                    line = '\t'.join(str(i) for i in split_items[: 8]) + '\t' + attr_demo.format(**attr_dict)
                else:
                    if t_id_key != 'transcript_id':
                        line = line.replace(t_id_key + ' ', 'transcript_id ', 1)
                line = line.strip()
                if 'gene_name' not in attr_dict:
                    gene = attr_dict.get('gene_id') or attr_dict.get('gene')
                    if gene is None:
                        self.logger.debug('this gtf file has no gene_id or gene key in attribute [9th column]')
                    else:
                        line += ' gene_name "{}"'.format(gene)

                return recored_type, t_id, self.__convert_chrom(split_items[0], line)
            return None, None, None

        return _filter

    def __anno_type_filter(self, anno_types):
        """此方法是用于过滤 GFF 文件
            通过注释类型列[例如: gene, transcript, exon, ...]过滤lncRNAs

        :param anno_types: 建议转录本对应的注释类型, 例如: lnc_RNA,antisense_RNA
        :return:
        """
        self.logger.debug(' =====  USE anno_type_filter function ===== ')
        attrs_keys = json.loads(self.option('retained_attrs'))
        has_gene_name = False
        old_gene = None
        for o, n in attrs_keys:
            if n == 'gene_name':
                has_gene_name = True
            if n == 'gene_id':
                old_gene = o
        if has_gene_name is not True:
            attrs_keys.append((old_gene, 'gene_name'))
        attr_demo = ' '.join('%s "{%s}";' % (k, i) for i, k in attrs_keys)
        t_id_key = [i for i, k in attrs_keys if k == 'transcript_id'][0]
        st_dict = {'is_w': False, 'exon_num': 0}
        transcript_dic = {}

        if 'exon' in anno_types:
            anno_types.discard('exon')
        gene_list = [i for i in anno_types if 'gene' in i]
        if len(gene_list) > 0:
            for i in gene_list:
                anno_types.discard(i)

        def _filter(split_items, line):

            attr_dict = split_items[8]
            # if 'ID' in attr_dict:
            #     attr_dict['ID'] = '-'.join(attr_dict['ID'].split('-')[1: -1])
            if split_items[2].strip() in anno_types:
                st_dict['is_w'] = True
                split_items[2] = 'transcript'
                if 'ID' in attr_dict:
                    transcript_dic['trans_id'] = attr_dict['ID'].replace('rna-', '')
                    attr_dict['ID'] = transcript_dic['trans_id']
                for i, k in attrs_keys:
                    if i not in attr_dict:
                        attrs_keys.remove([i, k])
                attr_demo = ' '.join('%s "{%s}";' % (k, i) for i, k in attrs_keys)
                new_line = '\t'.join(str(i) for i in split_items[: 8]) + '\t' + attr_demo.format(**attr_dict) + '\n'
                return 'transcript', attr_dict[t_id_key], self.__convert_chrom(split_items[0], new_line)
            elif split_items[2] == 'exon' and st_dict['is_w']:
                if 'ID' in attr_dict:
                    attr_dict['ID'] = transcript_dic['trans_id']
                for i, k in attrs_keys:
                    if i not in attr_dict:
                        attrs_keys.remove([i, k])
                attr_demo = ' '.join('%s "{%s}";' % (k, i) for i, k in attrs_keys)
                new_line = '\t'.join(str(i) for i in split_items[: 8]) + '\t' + attr_demo.format(
                    **attr_dict) + '\n'
                return 'exon', attr_dict[t_id_key], self.__convert_chrom(split_items[0], new_line)
            else:
                st_dict['is_w'] = False
                return None, None, None

        return _filter

    def run_filter(self):
        lnc_set = set()
        recored_dict = {}
        out_exon_num = None
        with self.option('gtf') as handler:
            basename = os.path.basename(handler.path)
            out_file = os.path.join(self.output_dir, basename[: -3] + 'gtf')
            with open(out_file, 'w') as out_handler:
                for items in handler:
                    if out_exon_num is None:
                        out_exon_num = 'exon_number' not in items[0][8]
                    recored_type, t_id, line = self.__filter_func(*items)
                    if line is None:
                        continue
                    lnc_set.add(t_id)
                    # out_handler.write(line)
                    if recored_type == 'transcript':
                        self.logger.debug(line)
                        if recored_dict:
                            trans_line = recored_dict['transcript'].strip() + '\n'
                            for i, (_, w_line) in enumerate(sorted(recored_dict['exon'], key=lambda item: item[0]), 1):
                                if out_exon_num:
                                    if w_line.strip().endswith(";"):
                                        trans_line += w_line.strip() + ' exon_number "%s";\n' % i
                                    else:
                                        trans_line += w_line.strip() + '; exon_number "%s";\n' % i
                                else:
                                    if w_line.strip().endswith(";"):
                                        trans_line += w_line + '\n'
                                    else:
                                        trans_line += w_line + ';\n'
                            out_handler.write(trans_line)
                        recored_dict[recored_type] = line
                        recored_dict['exon'] = []
                    else:
                        start, end = items[0][3: 5]
                        recored_dict['exon'].append(((start, end), line))

                if recored_dict:
                    trans_line = recored_dict['transcript'].strip() + '\n'
                    for i, (_, w_line) in enumerate(sorted(recored_dict['exon'], key=lambda item: item[0]), 1):
                        if out_exon_num:
                            if w_line.strip().endswith(";"):
                                trans_line += w_line.strip() + ' exon_number "%s";\n' % i
                            else:
                                trans_line += w_line.strip() + '; exon_number "%s";\n' % i
                        else:
                            if w_line.strip().endswith(";"):
                                trans_line += w_line + '\n'
                            else:
                                trans_line += w_line + ';\n'
                    out_handler.write(trans_line)

        lnc_ids_path = os.path.join(self.output_dir, basename[: -3] + 'lnc_ids.list')
        with open(lnc_ids_path, 'w') as out_handler:
            out_handler.write('\n'.join(lnc_set))

        # {'name': 'out_gtf', 'type': 'outfile', 'format': 'lnc_rna.common'},
        # {'name': 'out_ids', 'type': 'outfile', 'format': 'lnc_rna.common'},
        self.option('out_gtf').set_path(out_file)
        self.option('out_ids').set_path(lnc_ids_path)

    def run(self):
        super(GtfFilterTool, self).run()
        self._confirm_filter_func()
        self.run_filter()
        self.end()


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """

        def test(self):
            import random
            from mbio.workflows.single import SingleWorkflow
            from biocluster.wsheet import Sheet
            '''
            {'name': 'gtf', 'type': 'infile', 'format': 'lnc_rna.comm_gtf'},
            {'name': 'biotype_key', 'type': 'string', 'default': None},
            {'name': 'lnc_biotypes', 'type': 'string', 'default': lnc_biotypes},
            {'name': 'anno_types', 'type': 'string'},
            {'name': 'chrom_mapping', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'retained_attrs', 'type': 'string', 'default': None},
            '''
            data = {
                "id": "gtf_filter_" + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "lnc_rna.lnc_db.gtf_filter",
                "instant": False,
                "options": dict(
                    # gtf="/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/human/38/2019_02_27/data/GCF_000001405.38_GRCh38.p12_genomic.gff",
                    # gtf="/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/human/38/2019_02_27/data/Homo_sapiens.GRCh38.95.gtf",
                    # anno_types="lnc_RNA,antisense_RNA",
                    gtf='/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Saccharomyces_cerevisiae/ensembl/Saccharomyces_cerevisiae.R64-1-1.43.gtf',
                    biotype_key='transcript_biotype',
                    anno_types="transcript,exon",
                    # lnc_biotypes='sense_overlapping,antisense,retained_intron,sense_intronic,macro_lncRNA,' \
                    #              'bidirectional_promoter_lncRNA,3prime_overlapping_ncRNA,lincRNA,non_coding',
                    lnc_biotypes='ncRNA',
                    # retained_attrs='[["gene", "gene_id"], ["transcript_id", "transcript_id"], ["Dbxref", "Dbxref"]]',
                    retained_attrs='',
                    # chrom_mapping='/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/human/38/2019_02_27/data/chr_accessions_Human_38_p12'
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()


    unittest.main()
