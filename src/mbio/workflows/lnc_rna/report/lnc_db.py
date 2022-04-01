#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/3/16 23:21
@file    : lnc_db.py
"""
import json
import os
import unittest
from collections import OrderedDict

from biocluster.core.exceptions import OptionError
from biocluster.workflow import Workflow


class LncDbWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """option中dbs_info参数格式说明
            gtf文件过滤:
                biotype_key: 属性列中biotype关键字
                lnc_biotypes: 属于lncRNA的biotype的值, 例如ensembl中: sense_overlapping,antisense,等等
                anno_types: 需要留下的注释类型列[第3列], 该列可以作为过滤条件[ 作为过滤条件是上面两个为空]
                retained_attrs: 用于提取属性中值, 当anno_types作为过滤时必须提供用于过滤，当利用biotype过
                    滤时会修改原关键字名, 如：gene,gene_id; 就是把属性中gene名的改为gene_id
                chrom_mapping: 第一列名称转换表, 例如NCBI中第一列为location ID 需要转换为染色体编号
                    [ 文件两列: location_ID <TAB> chrom_number ]
            params = [
                {  # 列表中第一个元素为主要参考数据库, 其余会并入它
                    'gtf': '/.../Homo_sapiens.GRCh38.genome.gtf',
                    'biotype_key': 'transcript_biotype',
                    'lnc_biotypes': 'sense_overlapping,antisense,retained_intron,...',
                    'anno_types': 'transcript,exon',
                    'retained_attrs': '',
                    'chrom_mapping': '',
                    'db_name': 'ensembl',
                    'fasta': '/.../Homo_sapiens.GRCh38.ncrna.fa',
                },
                {
                    'gtf': '/.../GCF_000001405.38_GRCh38.p12_genomic.gff',
                    'biotype_key': '',
                    'lnc_biotypes': '',
                    'anno_types': 'lnc_RNA,antisense_RNA',
                    # retained_attrs 格式: 原属性名,新属性名;原属性名,新属性名;...
                    'retained_attrs': 'gene,gene_id;transcript_id,transcript_id;Dbxref,Dbxref',
                    'chrom_mapping': '/.../chr_accessions_Human_38_p12',
                    'db_name': 'ncbi',
                    'fasta': '/.../Homo_sapiens.GRCh38.ncrna.fa',
                }
            ]

            dbs_info = json.dumps(params)
        :param wsheet_object:
        """
        self._sheet = wsheet_object
        super(LncDbWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'ref_fa', 'type': 'infile', 'format': 'lnc_rna.lnc_fasta'},
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
            {'name': 'gmap_db_dir', 'type': 'infile', 'format': 'lnc_rna.common_dir'},
            {'name': 'gmap_db_name', 'type': 'string'},
            {'name': 'dbs_info', 'type': 'string'},  # json string of list<dict> - docstring
            {'name': 'db_order', 'type': 'string', 'default': 'ensembl,ncbi,noncode'}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.ref_lnc_db = None
        self.gtf_filter_tools = OrderedDict()
        self.params_dic = None
        self.tool_upload_list = None
        self.fa_filter_tool = None
        self.seq_align_tool = None
        self.merge_info_tool = None
        self.single_db_tool = None
        # self.dump_tool = self.api.api("lnc_rna.tf_predict")

    def check_options(self):
        for name in ('ref_fa', 'gmap_db_dir'):
            if not self.option(name).is_set:
                raise OptionError(name + ' must be assigned')

    def params_parser(self):
        demo_dict = {'gtf': '', 'biotype_key': '', 'lnc_biotypes': '', 'anno_types': 'transcript,exon',
                     'retained_attrs': '', 'chrom_mapping': '', 'db_name': '', 'fasta': ''}
        info_dic = OrderedDict()
        index_flag = 0
        for dic in json.loads(self.option('dbs_info')):
            tmp_dic = demo_dict.copy()
            tmp_dic.update(dic)
            db_name = tmp_dic['db_name']
            if not db_name:
                raise OptionError('db_name must be assigned, and can not be empty string or None')
            if index_flag == 0:
                index_flag += 1
                self.ref_lnc_db = db_name
            if tmp_dic['retained_attrs']:
                tmp_dic['retained_attrs'] = json.dumps(
                    [(k.strip(), v.strip()) for k, v in
                     (item.split(',') for item in tmp_dic['retained_attrs'].strip('; \n\r').split(';'))])
            if tmp_dic['gtf']:
                self.gtf_filter_tools[db_name] = self.add_tool('lnc_rna.lnc_db.gtf_filter')
            info_dic[db_name] = tmp_dic

        if len(info_dic) == 1:
            self.single_db_tool = self.add_tool('lnc_rna.lnc_db.single_db')
        else:
            self.fa_filter_tool = self.add_tool('lnc_rna.lnc_db.ids_mapping')
            self.seq_align_tool = self.add_tool('lnc_rna.lnc_db.seq_align')
            self.merge_info_tool = self.add_tool('lnc_rna.lnc_db.merge_info')

        self.params_dic = info_dic

    def gtf_filter(self):
        """
        {'name': 'gtf', 'type': 'infile', 'format': 'lnc_rna.comm_gtf'},
        {'name': 'biotype_key', 'type': 'string', 'default': ''},
        {'name': 'lnc_biotypes', 'type': 'string', 'default': lnc_biotypes},
        {'name': 'anno_types', 'type': 'string'},
        {'name': 'chrom_mapping', 'type': 'infile', 'format': 'lnc_rna.common'},
        {'name': 'retained_attrs', 'type': 'string', 'default': ''},
        {'name': 'out_gtf', 'type': 'outfile', 'format': 'lnc_rna.common'},
        {'name': 'out_ids', 'type': 'outfile', 'format': 'lnc_rna.common'},
        :return:
        """
        index_flag = 0
        for db_name, tool in self.gtf_filter_tools.items():
            if index_flag == 0:
                t_param = self.params_dic[db_name].copy()
                t_param.pop('db_name')
                t_param.pop('fasta')
                if not t_param['chrom_mapping']:
                    t_param.pop('chrom_mapping')
                if 'ids_mapping' in t_param:
                    t_param.pop('ids_mapping')
                self.logger.debug(tool.name + ' params: ' + json.dumps(t_param))
                tool.set_options(t_param)
                tool.run()

    def fa_filter(self):
        """
            [
                {
                    'db_name': '该fa文件所属数据库名',
                    'fa': 'fasta path',
                    'ids_mapping': '自己transcript id 和 参考id 间映射, 文件格式: self_id<TAB>ref_id',
                    'lnc_ids': 'lncrna id list file',
                },
            ]
        {'name': 'params', 'type': 'string', 'required': True},
        {'name': 'ref_db_name', 'type': 'string', 'required': True}
        {'name': 'out_fa', 'type': 'outfile', 'format': 'lnc_rna.common'},
        {'name': 'out_ids_map', 'type': 'outfile', 'format': 'lnc_rna.common'}
        :return:
        """
        index_flag = 0
        params_list = []
        for db_name, dic in self.params_dic.items():
            if index_flag == 0:
                index_flag += 1
                continue
            gtf_tool = self.gtf_filter_tools.get(db_name)
            lnc_ids = gtf_tool.option('out_ids').path if gtf_tool else ''
            tool_params = {
                'db_name': dic['db_name'],
                'fasta': dic['fasta'],
                'lnc_ids': lnc_ids
            }
            if 'ids_mapping'in dic and os.path.exists(dic["ids_mapping"]):
                tool_params["ids_mapping"] = dic["ids_mapping"]
            params_list.append(tool_params)
        if params_list:
            params_str = json.dumps(params_list)
            self.logger.debug(self.fa_filter_tool.name + ' params: ' + params_str)
            options = {
                'params': params_str,
                'ref_db_name': self.ref_lnc_db,
                'ref_lnc_list': self.gtf_filter_tools[self.ref_lnc_db].option('out_ids').path
            }
            self.fa_filter_tool.set_options(options)
            self.fa_filter_tool.on('end', self.seq_align)
            self.fa_filter_tool.run()
        else:
            self.single_db_tool.set_options(dict(
                ref_fa=self.option('ref_fa'),
                gtf=self.gtf_filter_tools[self.ref_lnc_db].option('out_gtf').path,
                db_name=self.ref_lnc_db
            ))
            self.tool_upload_list = (
                ('out_gtf', self.single_db_tool, 'lncrna.gtf'),
                ('out_fa', self.single_db_tool, 'lncrna.fa'),
                ('out_ids_map', self.single_db_tool, 'ids_matrix.xls')
            )
            self.single_db_tool.on('end', self.upload2workflow)
            self.single_db_tool.run()

    def seq_align(self):
        """
        {'name': 'ref_fa', 'type': 'infile', 'format': 'lnc_rna.common'},
        {'name': 'ref_db_gtf', 'type': 'infile', 'format': 'lnc_rna.common'},
        {'name': 'fasta', 'type': 'infile', 'format': 'lnc_rna.common'},
        {'name': 'gmap_db_dir', 'type': 'string', 'required': True},
        {'name': 'gmap_db_name', 'type': 'string', 'required': True},
        {'name': 'out_gtf', 'type': 'outfile', 'format': 'lnc_rna.common'},
        {'name': 'out_ids_map', 'type': 'outfile', 'format': 'lnc_rna.common'},
        :return:
        """
        params_dic = {
            'ref_fa': self.option('ref_fa').path,
            'ref_db_gtf': self.gtf_filter_tools[self.ref_lnc_db].option('out_gtf').path,
            'fasta': self.fa_filter_tool.option('out_fa').path,
            'gmap_db_dir': self.option('gmap_db_dir').path,
            'gmap_db_name': self.option('gmap_db_name')
        }
        self.seq_align_tool.set_options(params_dic)
        self.logger.debug(self.seq_align_tool.name + ' params: ' + json.dumps(params_dic))
        self.seq_align_tool.on('end', self.merge_info)
        self.seq_align_tool.run()

    def merge_info(self):
        """
        {'name': 'merged_gtf', 'type': 'infile', 'format': 'lnc_rna.comm_gtf'},
        {'name': 'ref_db_gtf', 'type': 'infile', 'format': 'lnc_rna.comm_gtf'},
        {'name': 'ref_fa', 'type': 'infile', 'format': 'lnc_rna.common'},
        {'name': 'known_ids', 'type': 'infile', 'format': 'lnc_rna.common'},
        {'name': 'novel_ids', 'type': 'infile', 'format': 'lnc_rna.common'},
        # {'name': 'ids_map', 'type': 'outfile', 'format': 'lnc_rna.common'},
        # {'name': 'out_fa', 'type': 'outfile', 'format': 'lnc_rna.common'},
        # {'name': 'ref_db_name', 'type': 'string'},
        :return:
        """
        params_dic = {
            'known_ids': self.fa_filter_tool.option('out_ids_map').path,
            'novel_ids': self.seq_align_tool.option('out_ids_map').path,
            'ref_db_gtf': self.gtf_filter_tools[self.ref_lnc_db].option('out_gtf').path,
            'merged_gtf': self.seq_align_tool.option('out_gtf').path,
            'ref_fa': self.option('ref_fa').path,
            'ref_db_name': self.ref_lnc_db,
            'db_order': self.option('db_order'),
            'ref_gtf': self.option('ref_gtf').path
        }
        self.logger.debug(self.merge_info_tool.name + ' params: ' + json.dumps(params_dic))
        self.merge_info_tool.set_options(params_dic)
        self.merge_info_tool.on('end', self.upload2workflow)
        self.merge_info_tool.run()

    def upload2workflow(self):
        if self.tool_upload_list is None:
            tool_list = (('out_gtf', self.merge_info_tool, 'lncrna.gtf'),
                         ('out_fa', self.merge_info_tool, 'lncrna.fa'),
                         ('ids_map', self.merge_info_tool, 'ids_matrix.xls'))
        else:
            tool_list = self.tool_upload_list
        for name, tool, new_name in tool_list:
            file_obj = tool.option(name)
            new_path = os.path.join(self.output_dir, new_name)
            file_obj.hard_link(new_path=new_path)
        self.end()

    def run(self):
        self.params_parser()
        self.on_rely(list(self.gtf_filter_tools.values()), self.fa_filter)
        self.gtf_filter()
        super(LncDbWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ["./lncrna.gtf", "", "lncRNA gtf"],
            ["./lncrna.fa", "", "lncRNA fasta"],
            ["./ids_matrix.xls", "", "lncRNA ids matrix in databases"],
        ])
        super(LncDbWorkflow, self).end()


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        '''
        This is test for the workflow. Just run this script to do test.
        '''

        # def test(self):
        #     from biocluster.wsheet import Sheet
        #     import random
        #     params = [
        #         {
        #             'gtf': '/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/human/38/'
        #                    '2019_02_27/data/Homo_sapiens.GRCh38.95.gtf',
        #             'biotype_key': 'transcript_biotype',
        #             'lnc_biotypes': 'sense_overlapping,antisense,retained_intron,sense_intronic,macro_lncRNA,'
        #                             'bidirectional_promoter_lncRNA,3prime_overlapping_ncRNA,lincRNA,non_coding',
        #             'anno_types': 'transcript,exon',
        #             'retained_attrs': '',
        #             'chrom_mapping': '',
        #             'db_name': 'ensembl',
        #             'fasta': '/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/human/'
        #                      '38/ensembl/Homo_sapiens.GRCh38.ncrna.fa',
        #         },
        #         # {
        #         #     'gtf': '/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncrna_tools/'
        #         #            'GCF_000001405.38_GRCh38.p12_genomic.gff',
        #         #     'biotype_key': '',
        #         #     'lnc_biotypes': '',
        #         #     'anno_types': 'lnc_RNA,antisense_RNA',
        #         #     # retained_attrs 格式: 原属性名,新属性名;原属性名,新属性名;...
        #         #     'retained_attrs': 'gene,gene_id;ID,transcript_id;Dbxref,Dbxref',
        #         #     'chrom_mapping': '/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/'
        #         #                      'human/38/2019_02_27/data/chr_accessions_Human_38_p12',
        #         #     'db_name': 'ncbi',
        #         #     'fasta': '/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/'
        #         #              'human/38/2019_02_27/data/GCF_000001405.38_GRCh38.p12_rna.fna',
        #         #     'ids_mapping': '/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/'
        #         #                    'human/38/2019_02_27/data/ncbi_transcript2ensembl'
        #         # }
        #     ]
        #     params = json.dumps(params)
        #     data = {
        #         'id': 'lnc_db_workflow_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
        #         'type': 'workflow',
        #         'name': 'lnc_rna.report.lnc_db',
        #         'options': {
        #             # {'name': 'ref_fa', 'type': 'infile', 'format': 'lnc_rna.lnc_fasta'},
        #             # {'name': 'gmap_db_dir', 'type': 'infile', 'format': 'lnc_rna.common_dir'},
        #             # {'name': 'gmap_db_name', 'type': 'string'},
        #             # {'name': 'dbs_info', 'type': 'string'}
        #             'ref_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/'
        #                       'Homo_sapiens/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fa',
        #             'gmap_db_dir': '/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/gmap/GRCh38_DB',
        #             'gmap_db_name': 'GRCh38',
        #             'dbs_info': params
        #         }
        #     }
        #     wsheet = Sheet(data=data)
        #     wf = LncDbWorkflow(wsheet)
        #     wf.sheet.id = 'lnc_rna'
        #     wf.sheet.project_sn = 'lnc_rna'
        #     wf.IMPORT_REPORT_DATA = True
        #     wf.IMPORT_REPORT_AFTER_DATA = False
        #     wf.run()

        def test(self):
            from biocluster.wsheet import Sheet
            import random
            params = [
                {
                    'gtf': '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Arabidopsis_thaliana/ensembl/Arabidopsis_thaliana.TAIR10.43.gtf',
                    'biotype_key': 'transcript_biotype',
                    'lnc_biotypes': 'lncRNA,ncRNA',
                    'anno_types': 'transcript,exon',
                    'retained_attrs': '',
                    'chrom_mapping': '',
                    'db_name': 'ensembl',
                    'fasta': '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Arabidopsis_thaliana/ensembl/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa',
                },
                {
                    'db_name': 'ncbi',
                    'fasta': '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Arabidopsis_thaliana/ncbi/ncbi_lncRNA.fa',
                    'ids_mapping': '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Arabidopsis_thaliana/ncbi/ncbi_2_ensmebl.txt'
                },
                {
                    'db_name': 'noncode',
                    'fasta': '/mnt/ilustre/users/sanger-dev/app/database/lnc_rna/noncode/NONCODEv5_arabidopsis.fa',
                    'ids_mapping': '/mnt/ilustre/users/sanger-dev/app/database/lnc_rna/noncode/noncode_id_2_ensembl_id.txt'
                },
                {
                    'db_name': 'greenc',
                    'fasta': '/mnt/ilustre/users/sanger-dev/app/database/lnc_rna/greenc/Arabidopsis_thaliana.fasta',
                    'ids_mapping': '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Arabidopsis_thaliana/greenc/greenc_2_ensembl.txt'
                }
            ]
            params = json.dumps(params)
            data = {
                'id': 'lnc_db_ath_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
                'type': 'workflow',
                'name': 'lnc_rna.report.lnc_db',
                'options': {
                    # {'name': 'ref_fa', 'type': 'infile', 'format': 'lnc_rna.lnc_fasta'},
                    # {'name': 'gmap_db_dir', 'type': 'infile', 'format': 'lnc_rna.common_dir'},
                    # {'name': 'gmap_db_name', 'type': 'string'},
                    # {'name': 'dbs_info', 'type': 'string'}
                    'ref_fa': '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Arabidopsis_thaliana/ensembl/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa',
                    'ref_gtf': '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Arabidopsis_thaliana/ensembl/Arabidopsis_thaliana.TAIR10.43.gtf',
                    'gmap_db_dir': '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Arabidopsis_thaliana/ensembl/Arabidopsis_thaliana',
                    'gmap_db_name': 'Arabidopsis_thaliana',
                    'db_order': 'ensembl,noncode,greenc',
                    'dbs_info': params
                }
            }
            wsheet = Sheet(data=data)
            wf = LncDbWorkflow(wsheet)
            wf.sheet.id = 'lnc_rna'
            wf.sheet.project_sn = 'lnc_rna'
            wf.IMPORT_REPORT_DATA = True
            wf.IMPORT_REPORT_AFTER_DATA = False
            wf.run()


    unittest.main()
