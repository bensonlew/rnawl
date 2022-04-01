#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/3/16 14:45
@file    : seq_align.py
"""
import csv
import datetime
import os
import re
import time
import unittest
from collections import defaultdict

import pandas as pd

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
from mbio.files.lnc_rna.comm_gtf import CommGtfFile


class SeqAlignAgent(Agent):
    def __init__(self, parent):
        super(SeqAlignAgent, self).__init__(parent)
        options = [
            {'name': 'ref_fa', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'ref_db_gtf', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'fasta', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'gmap_db_dir', 'type': 'string', 'required': True},
            {'name': 'gmap_db_name', 'type': 'string', 'required': True},
            {'name': 'out_gtf', 'type': 'outfile', 'format': 'lnc_rna.common'},
            {'name': 'out_ids_map', 'type': 'outfile', 'format': 'lnc_rna.common'},
            {'name': 'cpu', 'type': 'int', 'default': 20},
        ]
        self.__num = 21
        self.__memory = '50G'
        self.add_option(options)
        self.option('cpu', self.__num)
        self.step.add_steps("ids_mapping")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.ids_mapping.start()
        self.step.update()

    def step_end(self):
        self.step.ids_mapping.finish()
        self.step.update()

    def check_options(self):
        for name in ('ref_fa', 'ref_db_gtf', 'fasta'):
            if not self.option(name).is_set:
                raise OptionError(name + ' must be assigned')
        return True

    def set_resource(self):
        self.logger.debug('cpu: %s, memory: %s' % (self.__num, self.__memory))
        self._cpu = self.__num - 1
        self._memory = self.__memory

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(SeqAlignAgent, self).end()


class SeqAlignTool(Tool):
    def __init__(self, config):
        super(SeqAlignTool, self).__init__(config)
        self.tmp_work_dir = os.path.join(self.work_dir, 'tmp_work_dir')
        env_path = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/lnc_rna/cufflinks-2.2.1.Linux_x86_64/')
        self.set_environ(PATH=env_path)

    def run_cmd(self, cmd_name, cmd, is_wait=True, shell=False):
        cmd_obj = self.add_command(str(cmd_name), cmd, shell=shell, ignore_error=True)
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

    def gmap(self):
        """
        /mnt/ilustre/users/sanger-dev/app/bioinfo/align/gmap-2016-08-24/src/gmap
        -k 15
        --no-chimeras
        -f gff3_gene
        -n 1
        --gff3-add-separators=0
        -t 10
        -D /mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/gmap/GRCh38_DB
        -d GRCh38 20190317_162116246496_merged.fa > map.gff
        :return:
        """
        gmap = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/lnc_rna/gmap/bin/gmap')
        gmapl = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/lnc_rna/gmap/bin/gmapl')
        date = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")
        outfile = os.path.join(self.tmp_work_dir, date + '_gmap.gff3')
        cmd = '{gmap_soft} '
        cmd += ' -k 15 '
        cmd += ' --no-chimeras '
        cmd += ' -f gff3_gene '
        cmd += ' -n 1 '
        cmd += ' --gff3-add-separators=0 '
        cmd += ' -t {thread_num} '.format(thread_num=self.option('cpu'))
        cmd += ' -D {gmap_db_dir} '.format(gmap_db_dir=self.option('gmap_db_dir'))
        cmd += ' -d {gmap_db} '.format(gmap_db=self.option('gmap_db_name'))
        cmd += ' {fa_path} '.format(fa_path=self.option('fasta').path)
        cmd += ' > {outfile}'.format(outfile=outfile)

        cmd_obj = self.run_cmd('gmap', cmd.format(gmap_soft=gmap), is_wait=False, shell=True)
        self.wait(cmd_obj)
        if cmd_obj.return_code != 0:
            self.run_cmd('gmapl', cmd.format(gmap_soft=gmapl), is_wait=True, shell=True)
        return outfile

    def gff2gtf(self, outfile):
        gtf_obj = CommGtfFile()
        gtf_obj.set_path(outfile).set_type('gff')
        anno_type_set = {'mRNA', 'exon'}
        # 外显子变换匹配pattern
        exon_num_pattern = re.compile(r'.+\.exon(?P<num>\d+)$')
        # 文件输出格式化模板
        attr_demo = 'gene_id "{gene_id}"; transcript_id "{transcript_id}";'
        exon_demo = ' exon_number "{exon_number}"; target "{target}";'
        out_gtf = os.path.join(self.tmp_work_dir, os.path.basename(outfile)[: -4] + 'gtf')
        with gtf_obj as gtf_handler, open(out_gtf, 'w') as out_handler:
            for split_list, _ in gtf_handler:
                anno_type = split_list[2]
                if anno_type not in anno_type_set:
                    continue
                attr_dict = split_list[8]
                transcript_id = attr_dict['Name']
                gene_id = transcript_id + '.gene'
                split_list[2] = 'transcript' if anno_type == 'mRNA' else anno_type
                line = '\t'.join(str(i) for i in split_list[: 8])
                line += '\t' + attr_demo.format(gene_id=gene_id, transcript_id=transcript_id)
                if anno_type == 'exon':
                    exon_number = exon_num_pattern.match(attr_dict['ID']).group('num')
                    line += exon_demo.format(exon_number=exon_number, target=attr_dict['Target'])

                out_handler.write(line + '\n')
        gtf_list = os.path.join(self.tmp_work_dir, 'gtf_list')
        with open(gtf_list, 'w') as out_handler:
            out_handler.write(out_gtf)
        return gtf_list, out_gtf

    def cuffmerge(self, gtf_list_file):
        outdir = os.path.join(self.tmp_work_dir, 'merged_asm')
        outfile = os.path.join(outdir, 'merged.gtf')
        cuffmerge_soft = 'bioinfo/lnc_rna/cufflinks-2.2.1.Linux_x86_64/cuffmerge'

        cmd = '{} '.format(cuffmerge_soft)
        cmd += ' -o {} '.format(outdir)
        cmd += ' -g {} '.format(self.option('ref_db_gtf').path)
        cmd += ' -s {} '.format(self.option('ref_fa').path)
        cmd += ' -p {} '.format(self.option('cpu'))
        cmd += ' {} '.format(gtf_list_file)

        self.run_cmd('cuffmerge', cmd, is_wait=True, shell=False)
        # 设置gtf输出路径
        self.option('out_gtf').set_path(outfile)
        return outfile

    def cuffcompare(self, merge_gtf, *gtf_files):
        gtfs = ' '.join(gtf_files)
        outdir = os.path.join(self.tmp_work_dir, 'cuffcompare')
        os.makedirs(outdir)
        outprefix = os.path.join(outdir, 'cuffcmp')
        outfile = os.path.join(self.tmp_work_dir, 'cuffcmp.' + os.path.basename(gtf_files[-1]) + '.tmap')
        cuffcmp = 'bioinfo/lnc_rna/cufflinks-2.2.1.Linux_x86_64/cuffcompare'

        cmd = '{cuffcmp} '.format(cuffcmp=cuffcmp)
        cmd += '-r {} '.format(merge_gtf)
        cmd += '-o {} '.format(outprefix)
        cmd += gtfs

        self.run_cmd('gffcompare', cmd, is_wait=True, shell=False)
        return outfile

    def ids_mapping(self, tmap_file):
        res_dic = defaultdict(dict)
        with open(tmap_file) as in_handler:
            for dic in csv.DictReader(in_handler, delimiter='\t'):
                ref_id = dic['ref_id']
                cuff_id, db_name = dic['cuff_id'].rsplit('_', 1)
                class_code = dic['class_code']
                # 只取内含子完美匹配的结果
                if class_code != '=':
                    continue
                if ref_id not in res_dic[db_name]:
                    res_dic[db_name][ref_id] = []
                res_dic[db_name][ref_id].append(cuff_id)

        # ids间关系统计，用于查看一对多映射 [ 差错用 ]
        stat = []
        stat_file = os.path.join(self.output_dir, 'ids.stat')
        for key, val in res_dic.items():
            for sub_k, sub_v in val.items():
                if len(sub_v) > 1:
                    duff_ids = ','.join(sub_v)
                    stat.append('\t'.join((key, sub_k, duff_ids)))
                    val[sub_k] = duff_ids
                else:
                    val[sub_k] = sub_v[0]

        with open(stat_file, 'w') as out_handler:
            out_handler.write('db_name\tref_id\tids\n')
            out_handler.write('\n'.join(stat))

        # 输出不同数据库间id映射关系
        res_dic = {k + '_transcript_id': pd.Series(v) for k, v in res_dic.items()}
        map_df = pd.DataFrame(res_dic)
        map_df.fillna(value='')
        map_df.index.name = 'transcript_id'
        name = datetime.datetime.now().strftime('%Y_%m_%d_%H%M%S') + '_ids_mapping.xls'
        ids_mapping_file = os.path.join(self.output_dir, name)
        map_df.to_csv(ids_mapping_file, sep='\t')
        # 设置ids_map输出路径
        self.option('out_ids_map').set_path(ids_mapping_file)

    def run(self):
        super(SeqAlignTool, self).run()
        if not os.path.isdir(self.tmp_work_dir):
            os.makedirs(self.tmp_work_dir)
        gff_file = self.gmap()
        gtf_list, map_gtf = self.gff2gtf(gff_file)
        merged_gtf = self.cuffmerge(gtf_list)
        tmap_file = self.cuffcompare(merged_gtf, map_gtf)
        self.ids_mapping(tmap_file)
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
                "id": "SeqAlignTool_" + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "lnc_rna.lnc_db.seq_align",
                "instant": False,
                "options": dict(
                    # {'name': 'out_gtf', 'type': 'outfile', 'format': 'lnc_rna.common'},
                    # {'name': 'out_ids_map', 'type': 'outfile', 'format': 'lnc_rna.common'},

                    # {'name': 'ref_fa', 'type': 'infile', 'format': 'lnc_rna.common'},
                    ref_fa='/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/'
                           'Homo_sapiens/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fa',
                    # {'name': 'ref_db_gtf', 'type': 'infile', 'format': 'lnc_rna.common'},
                    ref_db_gtf='/mnt/ilustre/users/sanger-dev/workspace/20190315/Single_gtf_filter_143/'
                               'GtfFilter/output/Homo_sapiens.GRCh38.95.gtf',
                    # {'name': 'fasta', 'type': 'infile', 'format': 'lnc_rna.common'},
                    fasta='/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/human/38'
                          '/2019_02_27/data/20190317_162116246496_merged_part.fa',
                    # {'name': 'gmap_db_dir', 'type': 'string', 'required': True},
                    gmap_db_dir='/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/gmap/GRCh38_DB',
                    # {'name': 'gmap_db_name', 'type': 'string', 'required': True},
                    gmap_db_name='GRCh38'
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()


    unittest.main()
