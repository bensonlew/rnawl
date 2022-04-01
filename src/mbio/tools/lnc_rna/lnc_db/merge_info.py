#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/3/16 22:17
@file    : merge_info.py
"""

import csv
import json
import os
import time
import unittest
from collections import defaultdict
from itertools import chain

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool

from mbio.files.lnc_rna.lnc_gtf import LncGtfFile


class MergeInfoAgent(Agent):
    def __init__(self, parent):
        super(MergeInfoAgent, self).__init__(parent)
        options = [
            {'name': 'merged_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
            {'name': 'ref_db_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
            {'name': 'ref_fa', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'known_ids', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'novel_ids', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'ids_map', 'type': 'outfile', 'format': 'lnc_rna.common'},
            {'name': 'out_fa', 'type': 'outfile', 'format': 'lnc_rna.common'},
            {'name': 'out_gtf', 'type': 'outfile', 'format': 'lnc_rna.common'},
            {'name': 'ref_db_name', 'type': 'string'},
            {'name': 'cpu', 'type': 'int', 'default': 10},
            {'name': 'db_order', 'type': 'string', 'default': 'ensembl,ncbi,noncode,greenc'}
        ]
        self.__num = 2
        self.__memory = '4G'
        self.add_option(options)
        self.option('cpu', self.__num)
        self.step.add_steps("merge_info")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.merge_info.start()
        self.step.update()

    def step_end(self):
        self.step.merge_info.finish()
        self.step.update()

    def check_options(self):
        for name in ('known_ids', 'novel_ids', 'merged_gtf', 'ref_db_gtf', 'ref_fa', 'ref_gtf'):
            if not self.option(name).is_set:
                raise OptionError(name + ' must be assigned')
        return True

    def set_resource(self):
        self.logger.debug('cpu: %s, memory: %s' % (self.__num, self.__memory))
        self._cpu = self.__num
        self._memory = self.__memory

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(MergeInfoAgent, self).end()


class MergeInfoTool(Tool):
    def __init__(self, config):
        super(MergeInfoTool, self).__init__(config)
        self.tmp_work_dir = os.path.join(self.work_dir, 'compare_dir')
        if os.path.exists(self.tmp_work_dir):
            os.system('rm {}/*'.format(self.tmp_work_dir))
        else:
            os.mkdir(self.tmp_work_dir)

    def run_cmd(self, cmd_name, cmd, is_wait=True, shell=False):
        self.logger.debug(cmd_name + ': ' + cmd + '%s' % type(cmd_name))
        cmd_obj = self.add_command(str(cmd_name), cmd, shell=shell)
        if shell is True:
            cmd_obj.software_dir = ''
            cmd_obj._start_run_time = int(time.time())
        cmd_obj.run()
        if is_wait is True:
            self._check_stat(cmd_obj)
            return
        return cmd_obj

    def _check_stat(self, *cmd_objs):
        self.wait(*cmd_objs)
        for cmd_obj in cmd_objs:
            if cmd_obj.return_code == 0:
                self.logger.info("{} Finished successfully".format(cmd_obj.name))
            elif cmd_obj.return_code is None:
                self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_obj.name))
                cmd_obj.rerun()
                self.wait()
                if cmd_obj.return_code is 0:
                    self.logger.info("{} Finished successfully".format(cmd_obj.name))
                else:
                    self.set_error("{} Failed. >>>{}".format(cmd_obj.name, cmd_obj.cmd))
            else:
                self.set_error("{} Failed. >>>{}".format(cmd_obj.name, cmd_obj.cmd))

    def ids_map_reader(self, new_gtf, new_gtf_two):
        yield_set = set()
        with LncGtfFile(file_path=new_gtf) as new_gtf_handler, \
                LncGtfFile(file_path=new_gtf_two) as new_gtf_two_handler, \
                LncGtfFile(file_path=self.option('ref_db_gtf').path) as ref_gtf_handler:
            for split_list, _ in chain(ref_gtf_handler, new_gtf_handler, new_gtf_two_handler):
                attr_dict = split_list[8]
                transcript_id = attr_dict.get('transcript_id')
                if transcript_id in yield_set:
                    continue
                gene_name = attr_dict.get('gene_name', '')
                gene_id = attr_dict.get('gene_id')
                yield transcript_id, gene_id, gene_name
                yield_set.add(transcript_id)

    def merge(self, new_lnc_set, new_gtf_file, new_gtf_two, known2new_ids_map, new_trans2gene):
        # (new_trans_set, new_gtf, new_gtf_two, known_ids_map, new_trans2gene)
        db_order = [(i.strip(), i.strip() + '_transcript_id') for i in
                    self.option('db_order').strip(' \r\n,\t').split(',')]
        out_file = os.path.join(self.output_dir, 'ids_matrix.xls')
        ids_map = defaultdict(dict)

        with self.option('known_ids').get_reader() as known_handler, \
                self.option('novel_ids').get_reader() as novel_handler, \
                open(out_file, 'w') as out_handler:
            novel_fields = [i for i in novel_handler.readline().strip().split('\t') if i]
            known_fields = [i for i in known_handler.readline().strip().split('\t') if i]
            # 目的处理异常情况[用户没有提供所有数据库优先级列表]
            all_db_fields = {(i.strip().replace('_transcript_id', ''), i.strip()) for i in
                             chain(novel_fields, known_fields) if i != 'transcript_id'}
            db_order.extend(all_db_fields - set(db_order))

            novel_handler.seek(0)
            known_handler.seek(0)
            tmp_known_fields = known_fields[:]
            for i in novel_fields:
                if i not in tmp_known_fields:
                    known_fields.append(i)

            db_t_field = self.option('ref_db_name') + '_transcript_id'
            header_fields = known_fields[:]
            header_fields.insert(1, db_t_field)
            header_fields.append('gene_id')
            header_fields.append('gene_name')
            header_fields.insert(1, 'source')

            header = '\t'.join(header_fields) + '\n'
            out_handler.write(header)  # 输出文件字段名
            line_demo = '\t'.join('{%s}' % i for i in header_fields if i) + '\n'

            novel_dic = {dic.pop('transcript_id'): dic for dic in csv.DictReader(novel_handler, delimiter='\t')}
            known_dic = {dic.pop('transcript_id'): dic for dic in csv.DictReader(known_handler, delimiter='\t')}
            tmp_dic = {}  # 临时使用，目的减少创建字典对象个数
            # for key in all_ids:
            db_trans_fields = known_fields[1:]
            for transcript_id, gene_id, gene_name in self.ids_map_reader(new_gtf_file, new_gtf_two):
                is_new = False
                if transcript_id in new_lnc_set:
                    if transcript_id not in new_trans2gene:
                        continue
                    is_new = True

                novel = novel_dic.get(known2new_ids_map.get(transcript_id, transcript_id), tmp_dic)
                known = known_dic.get(transcript_id, tmp_dic)
                # self.logger.debug(str(novel) + ' === ' + str(known))
                if is_new:
                    dic = {
                        k: ','.join({i for i in novel.get(k, '').split(',') if i})
                        for k in db_trans_fields}
                    dic[db_t_field] = ''
                else:
                    def get_ids(key):
                        return [i for i in known.get(key, '').split(',') if i] or novel.get(key, '').split(',')

                    dic = {
                        k: ','.join({i for i in get_ids(k) if i}) for
                        k in db_trans_fields}
                    dic[db_t_field] = transcript_id
                # self.logger.debug(new_t_id + '  ' + str(dic))
                if 'gene_id' not in dic or not dic['gene_id']:
                    dic['gene_id'] = gene_id

                values_len = len([i for i in dic.values() if i])
                if values_len <= 1:
                    continue

                source = self.option('ref_db_name')
                if is_new:
                    for db_name, db_key in db_order:
                        new_trans_id = dic.get(db_key)
                        if new_trans_id:
                            sub_ids_map = ids_map[transcript_id]
                            source = db_name
                            transcript_id = new_trans_id.split(',', 1)[0]
                            sub_ids_map['source'] = source
                            sub_ids_map['new_id'] = transcript_id
                            break

                # dic['transcript_id'] = transcript_id
                # dic['gene_name'] = gene_name
                # dic['source'] = source
                self.logger.info(dic)
                self.logger.info(line_demo)
                new_line = line_demo.format(transcript_id=transcript_id, gene_name=gene_name, source=source, **dic)
                out_handler.write(new_line)

        self.option('ids_map').set_path(out_file)
        return ids_map

    def known_genes_info(self):
        res_dict = {}
        with self.option('ref_gtf') as gtf_handler:
            for split_list, _ in gtf_handler:
                attr_dict = split_list[8]
                gene_id = attr_dict['gene_id']
                gene_name = attr_dict.get('gene_name')
                if gene_name and gene_name not in res_dict:
                    res_dict[gene_name] = gene_id
                if gene_id and gene_id not in res_dict:
                    res_dict[gene_id] = gene_id
        return res_dict

    def out_gtf(self, gene_dic):
        # f_gtf = os.path.join(self.tmp_work_dir, 'filtered_merged.gtf')
        # cmd = '''/bin/grep -v 'class_code \\"=\\"' {gtf} > {out_gtf}'''.format(gtf=self.option('merged_gtf').path,
        #                                                                            out_gtf=f_gtf)
        # self.run_cmd('filter_merged_gtf', cmd, is_wait=True, shell=True)

        f_2_gtf = os.path.join(self.tmp_work_dir, 'filtered_1_merged.gtf')
        new_lnc_gtf = os.path.join(self.tmp_work_dir, 'new_lnc_merged.gtf')
        new_trans_set = set()
        trans_ids_map = {}
        new_trans2gene = {}
        with self.option('merged_gtf') as gtf_handler, \
                open(f_2_gtf, 'w') as out_handler, \
                open(new_lnc_gtf, 'w') as new_lnc_out:
            for split_list, _ in gtf_handler:
                attr_dict = split_list[8]
                trans_id = attr_dict['transcript_id']
                if 'class_code' in attr_dict:
                    class_code = attr_dict['class_code'].strip()
                else:
                    class_code = "none"
                if class_code == '=':
                    old_id = attr_dict['oId'].strip()
                    trans_ids_map[old_id] = trans_id
                    continue
                trans_ids_map[trans_id] = trans_id
                new_trans_set.add(trans_id)
                if 'gene_name' in attr_dict:
                    gene_name = attr_dict['gene_name']
                    new_name = gene_dic.get(gene_name)
                    if not new_name:
                        continue
                    attr_dict['gene_id'] = new_name
                    new_trans2gene[trans_id] = new_name
                    split_list[8] = ' '.join('{key} "{val}";'.format(key=k, val=v) for k, v in attr_dict.items())
                    line = '\t'.join(str(i) for i in split_list) + '\n'
                    out_handler.write(line)
                else:
                    split_list[8] = ' '.join('{key} "{val}";'.format(key=k, val=v) for k, v in attr_dict.items())
                    line = '\t'.join(str(i) for i in split_list) + '\n'
                    new_lnc_out.write(line)

        return f_2_gtf, new_trans_set, trans_ids_map, new_trans2gene, new_lnc_gtf

    def compare_gtf(self, new_gtf):

        ref_gtf = self.option('ref_gtf').path
        new_ref_gtf = os.path.join(self.tmp_work_dir, os.path.basename(ref_gtf))
        os.system('ln -s {} {}'.format(ref_gtf, new_ref_gtf))

        outprefix = os.path.join(self.tmp_work_dir, 'cuffcmp')
        outfile = os.path.join(self.tmp_work_dir, 'cuffcmp.' + os.path.basename(new_gtf) + '.tmap')
        cuffcmp = 'bioinfo/lnc_rna/cufflinks-2.2.1.Linux_x86_64/cuffcompare'

        cmd = '{cuffcmp} '.format(cuffcmp=cuffcmp)
        cmd += '-r {} '.format(new_ref_gtf)
        cmd += '-o {} '.format(outprefix)
        cmd += new_gtf

        self.run_cmd('gffcompare', cmd, is_wait=True, shell=False)
        return outfile

    def filter_gtf(self, tmap_file, new_gtf, ids_map_dict, new_trans2gene):
        res_dic = {}
        with open(tmap_file) as in_handler:
            for dic in csv.DictReader(in_handler, delimiter='\t'):
                ref_g_id = dic['ref_gene_id']
                cuff_id = dic['cuff_id']
                class_code = dic['class_code']
                # 只取内含子完美匹配的结果
                if class_code != 'c':
                    continue
                res_dic[cuff_id] = ref_g_id

        out_file = os.path.join(self.tmp_work_dir, 'filtered_' + os.path.basename(new_gtf))
        with LncGtfFile(file_path=new_gtf) as in_handler, \
                open(out_file, 'w') as out_handler:
            for split_list, _ in in_handler:
                attr_dict = split_list[8]
                trans_id = attr_dict['transcript_id']
                gene_name = res_dic.get(trans_id)
                if gene_name is None:
                    continue
                gene_id = ids_map_dict[gene_name]
                attr_dict['gene_id'] = gene_id
                attr_dict['gene_name'] = gene_name
                new_trans2gene[trans_id] = gene_id
                split_list[8] = ' '.join('{key} "{val}";'.format(key=k, val=v) for k, v in attr_dict.items())
                line = '\t'.join(str(i) for i in split_list) + '\n'
                out_handler.write(line)

        return out_file

    def merge_gtf(self, new_gtf, new_gtf_two, ids_map, new_trans_set):
        out_file = os.path.join(self.output_dir, 'lncrna.gtf')

        with LncGtfFile(file_path=new_gtf) as new_gtf_handler, \
                LncGtfFile(new_gtf_two) as new_two_handler, \
                LncGtfFile(file_path=self.option('ref_db_gtf').path) as ref_gtf_handler, \
                open(out_file, 'w') as out_handler:

            for split_list, _ in chain(ref_gtf_handler, new_gtf_handler, new_two_handler):
                attr_dict = split_list[8]
                if "transcript_id" not in attr_dict:
                    continue
                transcript_id = attr_dict['transcript_id']
                source = self.option('ref_db_name')
                if transcript_id in new_trans_set:
                    if transcript_id not in ids_map:
                        continue
                    ids_dic = ids_map[transcript_id]
                    attr_dict['transcript_id'] = ids_dic['new_id']
                    source = ids_dic['source']
                    # self.logger.debug(' ==========================================  ')
                    # self.logger.debug(json.dumps(ids_map))
                split_list[1] = source
                out_handler.write('\t'.join(str(i) for i in split_list[: 8]) + '\t' + ' '.join(
                    '{key} "{value}";'.format(key=key, value=value) for key, value in attr_dict.items()) + '\n')

        return out_file

    def gtf2fa(self, gtf):
        gffread = 'bioinfo/lnc_rna/cufflinks-2.2.1.Linux_x86_64/gffread'
        new_fa = os.path.join(self.output_dir, 'lncrna.fa')
        cmd = '{} {} '.format(gffread, gtf)
        cmd += '-g {} '.format(self.option('ref_fa').path)
        cmd += '-w {}'.format(new_fa)

        cmd_obj = self.run_cmd('gffread', cmd, is_wait=False)

        return cmd_obj, new_fa

    def out_stat(self, known_ids_map):
        with open(os.path.join(self.work_dir, 'known2new_ids.list'), 'w') as out_handler:
            out_handler.write('known_id\tnew_id\n')
            out_handler.write(
                ''.join('{known}\t{new}\n'.format(known=k, new=n) for k, n in known_ids_map.items() if k != n))

    def run(self):
        super(MergeInfoTool, self).run()
        genes_dic = self.known_genes_info()
        new_gtf, new_trans_set, known_ids_map, new_trans2gene, new_lnc_gtf = self.out_gtf(genes_dic)
        self.out_stat(known_ids_map)
        tmap_file = self.compare_gtf(new_lnc_gtf)
        new_gtf_two = self.filter_gtf(tmap_file, new_lnc_gtf, genes_dic, new_trans2gene)
        ids_map = self.merge(new_trans_set, new_gtf, new_gtf_two, known_ids_map, new_trans2gene)
        out_gtf = self.merge_gtf(new_gtf, new_gtf_two, ids_map, new_trans_set)
        cmd_obj, new_fa = self.gtf2fa(out_gtf)
        self.wait(cmd_obj)
        self.option('out_gtf').set_path(out_gtf)
        self.option('out_fa').set_path(new_fa)
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
            data = {
                "id": "MergeInfo_" + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "lnc_rna.lnc_db.merge_info",
                "instant": False,
                "options": dict(
                    # {'name': 'merged_gtf', 'type': 'infile', 'format': 'lnc_rna.comm_gtf'},
                    # {'name': 'ref_db_gtf', 'type': 'infile', 'format': 'lnc_rna.comm_gtf'},
                    # {'name': 'ref_fa', 'type': 'infile', 'format': 'lnc_rna.common'},
                    # {'name': 'known_ids', 'type': 'infile', 'format': 'lnc_rna.common'},
                    # {'name': 'novel_ids', 'type': 'infile', 'format': 'lnc_rna.common'},
                    # {'name': 'ids_map', 'type': 'outfile', 'format': 'lnc_rna.common'},
                    # {'name': 'out_fa', 'type': 'outfile', 'format': 'lnc_rna.common'},
                    # {'name': 'ref_db_name', 'type': 'string'},
                    merged_gtf='/mnt/ilustre/users/sanger-dev/workspace/20190425/LncDb_lnc_db_sac_1751_4036/SeqAlign/'
                               'tmp_work_dir/merged_asm/merged.gtf',
                    ref_db_gtf='/mnt/ilustre/users/sanger-dev/workspace/20190425/Single_gtf_filter_8086/GtfFilter/'
                               'output/Saccharomyces_cerevisiae.R64-1-1.43.gtf',
                    ref_gtf='/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Saccharomyces_cerevisiae/'
                            'ensembl/Saccharomyces_cerevisiae.R64-1-1.43.gtf',
                    ref_fa='/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Saccharomyces_cerevisiae/'
                           'ensembl/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa',
                    known_ids='/mnt/ilustre/users/sanger-dev/workspace/20190425/LncDb_lnc_db_sac_1751_4036/IdsMapping/'
                              'output/20190425_165825644489_merged_ids_maping.xls',
                    novel_ids='/mnt/ilustre/users/sanger-dev/workspace/20190425/LncDb_lnc_db_sac_1751_4036/SeqAlign/'
                              'output/2019_04_25_170002_ids_mapping.xls',
                    ref_db_name='ensembl'
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()


    unittest.main()
