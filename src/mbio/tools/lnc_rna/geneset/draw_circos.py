#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/3/28 12:44
@file    : draw_circos_bak.py
"""
import json
import os
import time
import unittest
from collections import OrderedDict, defaultdict
from itertools import chain

import pandas as pd

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
from mbio.tools.lnc_rna.geneset import circos_config
from mbio.packages.lnc_rna.geneset.chrom_dict import ChromDict


class DrawCircosAgent(Agent):
    """
    该tool用于整理出画circos的数据
    version 1.0
    """

    def __init__(self, parent):
        super(DrawCircosAgent, self).__init__(parent)
        options = [
            {"name": "target_trans", "type": "infile", "format": "lnc_rna.lnc_common"},
            {"name": "target_cis", "type": "infile", "format": "lnc_rna.lnc_common"},
            {"name": "lncrna_gtf", "type": "infile", "format": "lnc_rna.lnc_gtf"},
            {"name": "mrna_gtf", "type": "infile", "format": "lnc_rna.lnc_gtf"},
            {"name": "diff_exp", "type": "infile", "format": "lnc_rna.lnc_common"},
            {"name": "rna_type", "type": "infile", "format": "lnc_rna.lnc_common"},
            {"name": "ref_fa_fai", "type": "infile", "format": "lnc_rna.lnc_common"},
            {"name": "top_ref_num", "type": "int", "default": 10},
            {"name": "geneset_type", "type": "string", "default": 'T'},
            # lncRNA和mRNA圈是否只展示差异基因
            {"name": "only_diff", "type": "bool", "default": True},
        ]
        self.add_option(options)
        self.step.add_steps('DrawCircos')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.DrawCircos.start()
        self.step.update()

    def step_end(self):
        self.step.DrawCircos.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        """
        for name in ('target_trans', 'target_cis', 'lncrna_gtf', 'mrna_gtf', 'diff_exp', 'rna_type', 'ref_fa_fai'):
            if not self.option(name).is_set:
                raise OptionError('必须提供%s结果表' % name, code="34502201")
        if not isinstance(self.option('top_ref_num'), int):
            raise OptionError('必须提供top_ref_num结果表', code="34502201")
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(DrawCircosAgent, self).end()


class DrawCircosTool(Tool):

    def __init__(self, config):
        super(DrawCircosTool, self).__init__(config)
        self._version = "1.0.1"
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/WGS/circos-0.69-6/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/miniconda2/bin')
        self.script_path = self.config.PACKAGE_DIR + '/wgs/draw.circos.pl'
        self.perl_path = 'miniconda2/bin/perl '
        self.image_magick = "program/ImageMagick/bin/convert"

        self.tmp_dir = os.path.join(self.work_dir, 'conf')
        if not os.path.isdir(self.tmp_dir):
            os.mkdir(self.tmp_dir)

        self.geneset_type = self.option('geneset_type')

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

    def get_chr_detail(self, abud_or_len='len', chroms2ids_detail=None):
        chrom_detail = []
        chrom_len_dic = {}

        for items in self.option('ref_fa_fai').csv_reader():
            chrom = items[0].strip()
            chrom_len = int(items[1])
            if abud_or_len == 'len':
                chrom_detail.append((chrom, chrom_len))
            chrom_len_dic[chrom] = chrom_len
        if not (abud_or_len == 'len'):
            assert chroms2ids_detail is not None, 'must provide chroms2ids_detail param ' \
                                                  'when abud_or_len is abud[abundance]'
            exp_df = pd.read_table(self.option('diff_exp').path, header=0)
            exp_dic = {k: v for k, v in zip(exp_df['seq_id'], exp_df['avg_exp'])}
            for items in self.option('ref_fa_fai').csv_reader():
                chrom = items[0]
                total_exp = sum(exp_dic.get(i, 0) for i in chroms2ids_detail[chrom])
                chrom_detail.append((chrom, total_exp))

        tmp_chroms = sorted(chrom_detail, key=lambda item: item[1], reverse=True)
        tmp_chroms = {i for i, _ in tmp_chroms[: self.option('top_ref_num')]}
        self.logger.debug(str(tmp_chroms))
        colors = circos_config.chr_colors
        colors_len = len(colors)

        res_dict = OrderedDict()
        for i, (chr, chr_len) in enumerate(item for item in chrom_detail if item[0] in tmp_chroms):
            index = i % colors_len
            self.logger.debug('== {} color: {} == '.format(chr, colors[index]))
            res_dict['chr_' + chr] = {'chr_len': chrom_len_dic[chr], 'color': colors[index]}

        return res_dict

    def get_gene_detail(self):
        # seq_id_name = None
        if self.geneset_type == 'T':
            seq_id_name = 'transcript_id'
        elif self.geneset_type == 'G':
            seq_id_name = 'gene_id'
        else:
            raise Exception('geneset type must be T/G, but your type is "%s"' % self.geneset_type)

        gene_detail = defaultdict(dict)
        res_dic = {}
        chroms_detail = defaultdict(set)
        with self.option('lncrna_gtf') as lnc_handler, \
                self.option('mrna_gtf') as mrna_handler:
            for line_splits, line in chain(lnc_handler, mrna_handler):
                anno_type = line_splits[2]
                if anno_type != 'exon':
                    continue

                attr_dict = line_splits[8]
                seq_id = attr_dict[seq_id_name]
                seq_dic = res_dic.get(seq_id)
                start, end = line_splits[3], line_splits[4]
                chrom = line_splits[0]
                chroms_detail[chrom].add(seq_id)
                if seq_dic is None:
                    seq_dic = {'start': start, 'end': end, 'chr': 'chr_' + chrom}
                    res_dic[seq_id] = seq_dic
                    continue
                gene_sub = gene_detail['gene_id']
                if len(gene_sub) == 0:
                    gene_sub['start'] = start
                    gene_sub['end'] = end
                    gene_sub['chr'] = 'chr_' + chrom

                gene_sub['start'] = min(gene_sub['start'], start)
                gene_sub['end'] = max(gene_sub['end'], end)
                seq_dic['start'] = min(seq_dic['start'], start)
                seq_dic['end'] = max(seq_dic['end'], end)

        return res_dic, chroms_detail, gene_detail

    def get_diff(self):
        res_dict = {}
        for dic in self.option('diff_exp').csv_dict_reader():
            res_dict[dic.pop('seq_id')] = dic
        return res_dict

    def solve_map(self, chrom_dict, gtf_dict):
        # chrom, start, end, seq_id
        for seq_id, dic in gtf_dict.items():
            chrom_dict.add_map(dic['chr'], dic['start'], dic['end'], seq_id)

    def solve_exp(self, gtf_info, chrom_dict, exp_dict, is_log2fc=True):
        needed_chroms = chrom_dict.chrom_detail_dict
        self.logger.debug('solve_exp' + str(needed_chroms))
        v_key = 'log2fc' if is_log2fc else 'avg_exp'
        significant_set = set()
        only_diff = self.option('only_diff')
        for dic in self.option('rna_type').csv_dict_reader():
            # seq_id, rna_type
            try:
                seq_id = dic['seq_id']
                rna_type = dic['rna_type']
            except KeyError:
                self.logger.warning('rna_type文件中没有 seq_id和rna_type字段请检查避免错误')
                continue

            seq_dic = gtf_info.get(seq_id)
            # self.logger.debug('seq_dic: ' + str(seq_dic))
            if seq_dic is None:  # 过滤掉没有信息的，理论不存在，此处为了增加可用性
                continue
            chrom = seq_dic['chr']
            if chrom not in needed_chroms:  # 跳过不需要的染色体
                print(chrom)
                continue

            sub_exp_dic = exp_dict.get(seq_id, None)
            if sub_exp_dic is None:
                continue

            significant = sub_exp_dic['significant']
            if significant == 'yes':
                significant_set.add(seq_id)
            elif not only_diff:
                continue

            status = sub_exp_dic.get('regulate', 'no change')
            exp = sub_exp_dic[v_key]
            try:
                exp = float(exp)
            except ValueError:
                self.logger.warning('exp value error, please check out file')
                continue

            assert exp >= 0, 'solve_exp: exp must be positive real number'
            chrom_dict.add_exp(chrom, int(seq_dic['start']), int(seq_dic['end']), exp, status, seq_id,
                               is_lnc=True if rna_type == 'lncRNA' else False)
        return significant_set

    def solve_targets(self, chrom_dict, significant_set=None):
        l_id_name = 'lncrna_id' if self.geneset_type == 'T' else 'gene_id'
        m_id_name = 'mrna_id' if self.geneset_type == 'T' else 'mgene_id'
        significant_set = significant_set or set()
        for dic in self.option('target_trans').csv_dict_reader():
            lnc_id = dic[l_id_name]
            if lnc_id not in significant_set:  # 过滤掉非差异基因
                continue
            rna_id = dic[m_id_name]
            chrom_dict.add_targets(rna_id, target_is_cis=False, lnc_seq_id=lnc_id)
        for dic in self.option('target_cis').csv_dict_reader():
            lnc_id = dic[l_id_name]
            if lnc_id not in significant_set:  # 过滤掉非差异基因
                continue
            rna_id = dic[m_id_name]
            chrom_dict.add_targets(rna_id, target_is_cis=True, lnc_seq_id=lnc_id)

    @staticmethod
    def histogram_plot(file, r0, r1, min, max, thickness, fill_color, orientation):
        """

        :param file: # chromosome  start  end  data
        :param r0: 内半径[最好相对的, 例如: 05r]
        :param r1: 外半径
        :param min:
        :param max:
        :param thickness:
        :param fill_color: out/in
        :return:
        """
        return circos_config.histogram.format(file=file,
                                              r0=r0,
                                              r1=r1,
                                              min=min,
                                              max=max,
                                              thickness=thickness,
                                              fill_color=fill_color,
                                              orientation=orientation)

    @staticmethod
    def links_plot(file, radius, color, thickness):
        """
            <link>
            file = {file}
            radius = {radius}
            bezier_radius = 0r
            color = {color}
            thickness = {thickness}
            </link>

        :param file: # chromosome  start  end  data
        :param r0: 内半径[最好相对的, 例如: 05r]
        :param r1: 外半径
        :param min:
        :param max:
        :param thickness:
        :param fill_color: out/in
        :return:
        """
        return circos_config.links.format(file=file, radius=radius, color=color, thickness=thickness)

    def out_conf(self, chr_dict):
        """
        :param chr_dict:
        :return:
        """
        # chrom_colors.conf 染色体对应颜色
        chr_colors_conf_path = os.path.join(self.tmp_dir, 'chrom_colors.conf')
        with open(chr_colors_conf_path, 'w') as out_handler:
            color_str = '\n'.join(
                '{chrom} = {color}'.format(chrom='chr_color_' + chrom, color=dic['color']) for chrom, dic in
                chr_dict.items())
            chroms_str = ''.join('%s;' % i for i in chr_dict)
            out_handler.write(circos_config.chromosomes_and_color.format(chromosomes=chroms_str, colors=color_str))

        # ideogram.conf 配置文件
        ideogram_conf = circos_config.ideogram
        ideogram_conf_path = os.path.join(self.tmp_dir, 'ideogram.conf')
        with open(ideogram_conf_path, 'w') as out_handler:
            out_handler.write(ideogram_conf)

        # ticks.conf 配置文件
        tick_conf = circos_config.ticks
        tick_conf_path = os.path.join(self.tmp_dir, 'ticks.conf')
        with open(tick_conf_path, 'w') as out_handler:
            out_handler.write(tick_conf)

        # colors_conf
        color_conf = circos_config.colors_string
        color_conf_path = os.path.join(self.tmp_dir, 'm_colors.conf')
        with open(color_conf_path, 'w') as out_handler:
            out_handler.write(color_conf)

        # housekeeping_conf
        housekeeping_conf = circos_config.housekeeping
        housekeeping_conf_path = os.path.join(self.tmp_dir, 'housekeeping.conf')
        with open(housekeeping_conf_path, 'w') as out_handler:
            out_handler.write(housekeeping_conf)

        return {'chrom_colors_conf': chr_colors_conf_path,
                'ideogram_conf': ideogram_conf_path,
                'tick_conf': tick_conf_path,
                'color_conf': color_conf_path,
                'housekeeping_conf': housekeeping_conf_path}

    def out_data(self, chrom_dict, is_log2fc):
        """
        'chr': 'chr	-	{chrom}	chr1	0	{end}	{chrom}',
        'gene': 'band	{chrom}	{band}	{band}	{start}	{end}	{color}'

        :param chrom_dict:
        :return:
        """
        data_dir = os.path.join(self.tmp_dir, 'data')
        if not os.path.isdir(data_dir):
            os.mkdir(data_dir)
        file_dict = {}
        # karyotype 染色体外环
        karyotype = os.path.join(data_dir, 'chr_band.txt')
        file_dict['karyotype'] = karyotype
        # mRNA 柱状图
        mrna_histogram_out = os.path.join(data_dir, 'mrna_histogram_out.txt')
        file_dict['mrna_histogram_out'] = mrna_histogram_out
        mrna_histogram_in = os.path.join(data_dir, 'mrna_histogram_in.txt')
        file_dict['mrna_histogram_in'] = mrna_histogram_in
        # lRNA 柱状图
        lncrna_histogram_out = os.path.join(data_dir, 'lncrna_histogram_out.txt')
        file_dict['lncrna_histogram_out'] = lncrna_histogram_out
        lncrna_histogram_in = os.path.join(data_dir, 'lncrna_histogram_in.txt')
        file_dict['lncrna_histogram_in'] = lncrna_histogram_in
        # cis 柱状图数据
        cis_histogram_out = os.path.join(data_dir, 'cis_histogram_out.txt')
        cis_histogram_in = os.path.join(data_dir, 'cis_histogram_in.txt')
        file_dict['cis_histogram_out'] = cis_histogram_out
        file_dict['cis_histogram_in'] = cis_histogram_in
        # trans 柱状图数据
        comm_chr_trans_cis = os.path.join(data_dir, 'comm_trans_histogram.txt')
        diff_chr_trans_cis = os.path.join(data_dir, 'diff_trans_histogram.txt')
        file_dict['comm_chr_trans_cis'] = comm_chr_trans_cis
        file_dict['diff_chr_trans_cis'] = diff_chr_trans_cis

        chr_demo = circos_config.chr_data['chr']
        histogram_line_demo = circos_config.mrna_histo_demo
        cis_line_demo = circos_config.cir_line_demo

        max_min_stat = defaultdict(dict)
        # 核型数据，取消了
        # gene_demo = circos_config.chr_data['gene']
        # band_colors = circos_config.band_colors
        # band_color_num = len(band_colors)
        # columns_num = chrom_dict.columns_num
        with open(karyotype, 'w') as out_handler, \
                open(mrna_histogram_out, 'w') as mout_handler, \
                open(mrna_histogram_in, 'w') as min_handler, \
                open(lncrna_histogram_out, 'w') as lout_handler, \
                open(lncrna_histogram_in, 'w') as lin_handler, \
                open(cis_histogram_out, 'w') as cis_out_handler, \
                open(cis_histogram_in, 'w') as cis_in_handler, \
                open(diff_chr_trans_cis, 'w') as diff_trans_handler, \
                open(comm_chr_trans_cis, 'w') as comm_trans_handler:
            sum_basic = 0
            chroms = []
            # 染色体 条带 [最外环]
            for key, dic in chrom_dict.chrom_detail_dict.items():
                chroms.append(key)
                chr_len = dic['chr_len']
                sum_basic += chr_len
                out_handler.write(
                    chr_demo.format(chrom=key, end=chr_len, color='chr_color_' + key, chrom_lable=key[4:].replace("_", "-")))

            mexp_up_num = 0
            lexp_up_num = 0
            item_grep = 'item_grep'  # '128,128,128'
            item_red = 'item_red'  # '255,0,0'
            item_blue = 'item_blue'  # 0,0,255
            item_chen = 'item_chen'  # '255,102,51'
            item_green = 'item_green'  # '0,255,0'

            mrna_out_info = []
            mrna_in_info = []
            lrna_out_info = []
            lrna_in_info = []
            cis_out_info = []
            cis_in_info = []
            comm_chr_trans = []
            diff_chr_trans = []

            for chrom in chroms:
                sub_chrom_info = chrom_dict.index_dict[chrom]
                # self.logger.debug(chrom + ' == ' + chrom_dict[sub_chrom_info[0]].chrom + '===' + str(
                #     chrom_dict[sub_chrom_info[0]].start))
                # self.logger.debug(chrom + ' == ' + chrom_dict[sub_chrom_info[-1]].chrom + '===' + str(
                #     chrom_dict[sub_chrom_info[-1]].start))
                mrna_out_one = False
                mrna_in_one = False
                lrna_out_one = False
                lrna_in_one = False
                cis_out_one = False

                for index in sub_chrom_info:
                    # circos图外圈数据
                    item = chrom_dict[index]
                    start, end = int(item.start) + 1, int(item.end) - 1

                    # color = band_colors[flag % band_color_num]
                    # band_name = '{}_{}_{}'.format(chrom, 'band', flag)
                    # 核型数据 [染色体上条带]
                    # out_handler.write(
                    #     gene_demo.format(chrom=chrom, band=band_name, start=start, end=end, color=color))

                    # mid_idx = (start + end) // 2
                    # MRNA柱状图数据, 方向朝外的柱子 [上调]
                    m_dic = max_min_stat['h_mrna']
                    m_up_num, m_up_exp = item.exp(is_lnc=False, up_down='up')
                    if m_up_num > 0 and (is_log2fc is False and m_up_exp > 0):
                        m_dic['up'] = m_dic.get('up', 0) + m_up_exp
                        mexp_up_num += 1
                        mrna_out_info.append(histogram_line_demo.format(
                            chrom=chrom, start=start, end=end, value=m_up_exp, color=item_red))
                    elif mrna_out_one is False:
                        mrna_out_one = True
                        mrna_out_info.append(histogram_line_demo.format(
                            chrom=chrom, start=start, end=end, value=0.01, color=item_red))
                    if len(mrna_out_info) >= 500:
                        mout_handler.write(''.join(mrna_out_info))
                        mrna_out_info = []

                    # MRNA柱状图数据, 方向朝内的柱子 [下调]
                    m_down_num, m_down_exp = item.exp(is_lnc=False, up_down='down')
                    if m_down_num > 0 and (is_log2fc is False and m_down_exp > 0) or is_log2fc is True:
                        m_dic['up'] = m_dic.get('up', 0) + m_down_exp
                        mexp_up_num += 1
                        m_down_exp = -m_down_exp if m_down_exp > 0 else 0  # 防止浮点数造成影响
                        mrna_in_info.append(histogram_line_demo.format(
                            chrom=chrom, start=start, end=end, value=m_down_exp, color=item_blue))
                    elif mrna_in_one is False:
                        mrna_in_one = True
                        mrna_in_info.append(histogram_line_demo.format(
                            chrom=chrom, start=start, end=end, value=-0.01, color=item_blue))
                    if len(mrna_in_info) >= 500:
                        min_handler.write(''.join(mrna_in_info))
                        mrna_in_info = []

                    # lncRNA柱状图数据, 方向朝外的柱子 [上调]
                    l_dic = max_min_stat['h_lncrna']
                    l_up_num, l_up_exp = item.exp(is_lnc=True, up_down='up')
                    if l_up_num == 0 and lrna_out_one is False:
                        lrna_out_one = True
                        lrna_out_info.append(histogram_line_demo.format(
                            chrom=chrom, start=start, end=end, value=0.01, color=item_chen))
                    elif l_up_exp > 0:
                        l_dic['up'] = l_dic.get('up', 0) + l_up_exp
                        lexp_up_num += 1
                        lrna_out_info.append(histogram_line_demo.format(
                            chrom=chrom, start=start, end=end, value=l_up_exp, color=item_chen))
                    if len(lrna_out_info) >= 500:
                        lout_handler.write(''.join(lrna_out_info))
                        lrna_out_info = []

                    # lncRNA柱状图数据, 方向朝内的柱子 [下调]
                    l_down_num, l_down_exp = item.exp(is_lnc=True, up_down='down')
                    if l_down_num > 0 and (is_log2fc is False and l_down_exp > 0) or is_log2fc is True:
                        l_dic['up'] = l_dic.get('up', 0) + l_down_exp
                        lexp_up_num += 1
                        l_down_exp = -l_down_exp if l_down_exp > 0 else 0.01  # 防止浮点数造成影响
                        lrna_in_info.append(histogram_line_demo.format(
                            chrom=chrom, start=start, end=end, value=l_down_exp, color=item_green))
                    elif lrna_in_one is False:
                        lrna_in_one = True
                        lrna_in_info.append(histogram_line_demo.format(
                            chrom=chrom, start=start, end=end, value=0.01, color=item_green))
                    if len(lrna_in_info) >= 500:
                        lin_handler.write(''.join(lrna_in_info))
                        lrna_in_info = []

                    # cis 柱状图数据
                    targets = item.cis_set
                    if len(targets) > 0:
                        cis_out_info.append(histogram_line_demo.format(
                            chrom=chrom, start=start, end=end, value=0.25, color='item_blue'))
                    elif not cis_out_one:
                        cis_out_one = True
                        cis_out_info.append(histogram_line_demo.format(
                            chrom=chrom, start=start, end=end, value=0.00001, color='item_blue'))
                        cis_in_info.append(
                            histogram_line_demo.format(chrom=chrom, start=start, end=end, value=-0.00001,
                                                       color='item_red'))
                    if len(cis_out_info) >= 500:
                        cis_out_handler.write(''.join(cis_out_info))
                        cis_out_info = []

                    for cis_item in targets:
                        cis_in_info.append(
                            histogram_line_demo.format(chrom=cis_item.chrom, start=cis_item.start, end=cis_item.end,
                                                       value=-0.25, color='item_red'))

                    if len(cis_in_info) >= 500:
                        cis_in_handler.write(''.join(cis_in_info))
                        cis_in_info = []

                    # trans links 柱状图数据
                    targets = sorted(item.trans_set, key=lambda itm: (itm.chrom_index, itm.start))
                    for target_item in targets:
                        t_chrom, t_start, t_end = target_item.chrom, target_item.start, target_item.end
                        if item.in_comm_chrom(target_item):
                            # '{cir_name}	{chrom}	{start}	{end}	{t_name}	{t_chrom}	{t_start}	{t_end}	color={color}'
                            comm_chr_trans.append(cis_line_demo.format(
                                chrom=chrom, start=start, end=end,
                                t_chrom=t_chrom, t_start=t_start, t_end=t_end, color=item_chen
                            ))
                        else:
                            diff_chr_trans.append(cis_line_demo.format(
                                chrom=chrom, start=start, end=end,
                                t_chrom=t_chrom, t_start=t_start, t_end=t_end, color=item_blue
                            ))
                    if len(comm_chr_trans) >= 500:
                        comm_trans_handler.write(''.join(comm_chr_trans))
                        comm_chr_trans = []
                    if len(diff_chr_trans) >= 500:
                        diff_trans_handler.write(''.join(diff_chr_trans))
                        diff_chr_trans = []
            # mrna_out_info = []
            #             mrna_in_info = []
            #             lrna_out_info = []
            #             lrna_in_info = []
            #             cis_out_info = []
            #             cis_in_info = []
            #             comm_chr_trans = []
            #             diff_chr_trans = []
            if mrna_out_info:
                mout_handler.write(''.join(mrna_out_info))
            if mrna_in_info:
                min_handler.write(''.join(mrna_in_info))
            if lrna_out_info:
                lout_handler.write(''.join(lrna_out_info))
            if lrna_in_info:
                lin_handler.write(''.join(lrna_in_info))
            if cis_out_info:
                cis_out_handler.write(''.join(cis_out_info))
            if cis_in_info:
                cis_in_handler.write(''.join(cis_in_info))
            if comm_chr_trans:
                comm_trans_handler.write(''.join(comm_chr_trans))
            if diff_chr_trans:
                diff_trans_handler.write(''.join(diff_chr_trans))

        if mexp_up_num > 0:
            max_min_stat['h_mrna']['up'] = 2 * round((max_min_stat['h_mrna'].get('up') / mexp_up_num), 5)
        if lexp_up_num > 0:
            max_min_stat['h_lncrna']['up'] = 2 * round((max_min_stat['h_lncrna'].get('up') / mexp_up_num), 5)

        return file_dict, max_min_stat

    def out_main_conf(self, conf_dic, data_dic, max_min_stat):
        """
            karyotype = {karyotype}
            <<include {ideogram_conf}>>
            <<include {ticks_conf}>>
            <<include {chromosomes_and_color_conf}>>
            {plots}
            {links}
            <<include {colors_conf}>>
            <<include {housekeeping_conf}>>
        :param conf_dic:
        :param data_dic:
        :return:
        """
        h_mrna = max_min_stat['h_mrna']
        h_lncrna = max_min_stat['h_lncrna']
        l_max_val = max(abs(h_lncrna.get('down', 0)), h_lncrna.get('up', 0))
        m_max_val = max(abs(h_mrna.get('down', 0)), h_mrna.get('up', 0))

        if l_max_val == 0 and m_max_val == 0:
            l_max_val = m_max_val = 0.5
        elif l_max_val == 0:
            l_max_val = m_max_val
        elif m_max_val == 0:
            m_max_val = l_max_val

        mrna_histogram_out = data_dic['mrna_histogram_out']
        m_out_histo = self.histogram_plot(mrna_histogram_out, '0.8501r', '0.9r', min=0, max=m_max_val,
                                          thickness=1, fill_color='', orientation='out')
        mrna_histogram_in = data_dic['mrna_histogram_in']
        m_in_histo = self.histogram_plot(mrna_histogram_in, '0.8r', '0.8499r', min=-m_max_val, max=0,
                                         thickness=1, fill_color='', orientation='out')

        lncrna_histogram_out = data_dic['lncrna_histogram_out']
        l_out_histo = self.histogram_plot(lncrna_histogram_out, '0.7001r', '0.75r', min=0, max=l_max_val,
                                          thickness=1, fill_color='', orientation='out')
        lncrna_histogram_in = data_dic['lncrna_histogram_in']
        l_in_histo = self.histogram_plot(lncrna_histogram_in, '0.65r', '0.6999r', min=-l_max_val, max=0,
                                         thickness=1, fill_color='', orientation='out')

        cis_histogram_out = data_dic['cis_histogram_out']
        c_out_histo = self.histogram_plot(cis_histogram_out, '0.55r', '0.6r', min='0', max='1', thickness=1,
                                          fill_color='', orientation='out')
        cis_histogram_in = data_dic['cis_histogram_in']
        c_in_histo = self.histogram_plot(cis_histogram_in, '0.5r', '0.55r', min='-1', max='0', thickness=1,
                                         fill_color='', orientation='out')

        comm_chr_trans_cis = data_dic['comm_chr_trans_cis']
        comm_cis = self.links_plot(comm_chr_trans_cis, radius='0.45r', color='item_chen', thickness=4)
        diff_chr_trans_cis = data_dic['diff_chr_trans_cis']
        diff_cis = self.links_plot(diff_chr_trans_cis, radius='0.45r', color='blue', thickness=4)

        # 配置文件
        main_conf_str = circos_config.main_conf.format(
            karyotype=data_dic['karyotype'],
            ideogram_conf=conf_dic['ideogram_conf'],
            ticks_conf=conf_dic['tick_conf'],
            plots='\n\n'.join((m_out_histo, m_in_histo, l_out_histo, l_in_histo, c_out_histo, c_in_histo)),
            links='\n\n'.join((comm_cis, diff_cis)),
            colors_conf=conf_dic['color_conf'],
            housekeeping_conf=conf_dic['housekeeping_conf'],
            chromosomes_and_color_conf=conf_dic['chrom_colors_conf'],
        )
        main_conf_file = os.path.join(self.tmp_dir, 'circos_main.conf')
        with open(main_conf_file, 'w') as out_handler:
            out_handler.write(main_conf_str)

        return main_conf_file

    def draw_picture(self, main_conf):
        """
        $ ./app/bioinfo/WGS/circos-0.69-6/bin/circos --help
            Usage:
                  # without -conf Circos will search for configuration
                  circos

                  # use specific configuration file
                  circos -conf circos.conf

                  # detailed debugging for code components
                  # see http://www.circos.ca/documentation/tutorials/configuration/debugging
            circos -conf $outdir/draw.circos/draw.conf -outputfile circos -outputdir $outdir/
        :param main_conf:
        :return:
        """
        cmd = self.config.SOFTWARE_DIR + '/bioinfo/WGS/circos-0.69-6/bin/circos '
        cmd += ' -conf {} '.format(main_conf)
        cmd += ' -outputfile circos '
        cmd += '-outputdir {} '.format(self.output_dir)

        self.run_cmd('draw_picture', cmd, is_wait=True, shell=True)

    def convert_png_pdf(self):
        cmd = self.image_magick
        cmd += ' -flatten '
        cmd += ' -quality 100 '
        cmd += ' -density 130 '
        cmd += ' -background white '
        cmd += os.path.join(self.output_dir, 'circos.png')
        cmd += ' '
        cmd += os.path.join(self.output_dir, 'circos.pdf')

        self.run_cmd('png2pdf', cmd, is_wait=True)

    def convert_png_pdf_new(self):
        try:
            import cairosvg
            png = os.path.join(self.output_dir, 'circos.svg')
            pdf = os.path.join(self.output_dir, 'circos.pdf')
            cairosvg.svg2pdf(url=png, write_to=pdf)
        except:
            self.convert_png_pdf()

    def run(self):
        """
        运行
        """
        super(DrawCircosTool, self).run()
        # 获取基因详情, gene或者transcript长度、染色体编号、start、end
        gtf_dict, chroms_detail, gene_detail = self.get_gene_detail()
        # 获取染色体长度信息
        # abud_or_len='len' / 'abud'
        chr_info_dict = self.get_chr_detail(abud_or_len='abud', chroms2ids_detail=chroms_detail)
        # 差异基因信息
        diff_dict = self.get_diff()
        # 染色体分类对象
        chrom_dict = ChromDict(chr_info_dict, bind_obj=self, gene_detail=gene_detail, col_num=500)
        # self.logger.debug('chrom_dict: ' + str(dict(chrom_dict)))
        # self.logger.debug('chrom_dict.chroms: ' + str(dict(chrom_dict.chrom_detail_dict)))
        # 处理seq_id和chrom映射
        self.solve_map(chrom_dict, gtf_dict)
        # 添加表达信息各个seq到划分的区域
        self.logger.debug('FUNCTION: solve_exp start =====')
        significant_set = self.solve_exp(gtf_info=gtf_dict, chrom_dict=chrom_dict, exp_dict=diff_dict, is_log2fc=False)
        # 添加靶点到划分区域
        self.logger.debug('FUNCTION: solve_targets start =====')
        self.solve_targets(chrom_dict=chrom_dict, significant_set=significant_set)
        # 输出配置文件
        self.logger.debug('FUNCTION: out_conf start =====')
        conf_files_dic = self.out_conf(chr_dict=chr_info_dict)
        # 输出数据文件
        self.logger.debug('FUNCTION: out_data start =====')
        data_files_dic, max_min_stat = self.out_data(chrom_dict=chrom_dict, is_log2fc=False)
        # 输出主配置文件
        self.logger.debug('FUNCTION: out_main_conf start =====')
        main_conf_file = self.out_main_conf(conf_dic=conf_files_dic, data_dic=data_files_dic, max_min_stat=max_min_stat)
        # 输出图片
        self.draw_picture(main_conf_file)
        self.convert_png_pdf_new()

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
            {"name": "target_trans", "type": "infile", "format": "lnc_rna.lnc_common"},
            {"name": "target_cis", "type": "infile", "format": "lnc_rna.lnc_common"},
            {"name": "lncrna_gtf", "type": "infile", "format": "lnc_rna.lnc_gtf"},
            {"name": "mrna_gtf", "type": "infile", "format": "lnc_rna.lnc_gtf"},
            {"name": "diff_exp", "type": "infile", "format": "lnc_rna.lnc_common"},
            {"name": "rna_type", "type": "infile", "format": "lnc_rna.lnc_common"},
            {"name": "ref_fa_fai", "type": "infile", "format": "lnc_rna.lnc_common"},
            {"name": "top_ref_num", "type": "int", "default": 10},
            {"name": "geneset_type", "type": "string", "default": 'T'}
            '''
            data = {
                "id": "draw_circos_" + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "lnc_rna.geneset.draw_circos",
                "instant": False,

                "options": dict(
                    target_trans='/mnt/ilustre/users/sanger-dev/workspace/20190401/DrawCircos_lnc_rna_8841_1244/'
                                 'trans_targets.txt',
                    target_cis='/mnt/ilustre/users/sanger-dev/workspace/20190401/DrawCircos_lnc_rna_8841_1244/'
                               'cis_targets.txt',
                    gtf='/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/assemble/method-cufflinks/'
                        'output/NewTranscripts/ref_and_new.gtf',
                    diff_exp='/mnt/ilustre/users/sanger-dev/workspace/20190425/DrawCircos_lnc_rna_7972_151/'
                             'diff_exp.txt',
                    rna_type='/mnt/ilustre/users/sanger-dev/workspace/20190401/DrawCircos_lnc_rna_8841_1244/'
                             'seq_id_type.txt',
                    ref_fa_fai='/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens'
                               '/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fa.fai',
                    top_ref_num=10,
                    geneset_type='T'
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()


    unittest.main()
