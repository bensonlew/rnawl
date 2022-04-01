#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/3/29 8:43
@file    : chrom_detail_dict.py
@author  : zhipeng.zhao
"""
import json
from collections import defaultdict

import math


class ExpItem(object):
    def __init__(self, item_type):
        self.item_type = item_type
        self.up_num = 0
        self.up_exp = 0
        self.max_up_exp = 0
        self.down_num = 0
        self.down_exp = 0
        self.min_down_exp = 0
        self.neutral_num = 0
        self.neutral_exp = 0
        self.__status = None

    def add_exp(self, exp, exp_type):
        if exp_type == 'up':
            # assert exp >= 0, 'when the type is up, exp must be greater than 0 == ' + str(exp)
            self.up_num += 1
            print('up: ' + str(exp))
            self.up_exp = (self.up_exp * self.up_num + exp) / (self.up_num + 1)
            self.max_up_exp = max(self.max_up_exp, exp)
        elif exp_type == 'down':
            # assert exp >= 0, 'when the type is down, exp must be greater than 0 == ' + str(exp)
            self.down_num += 1
            print('down: ' + str(exp))
            # assert exp < 0, 'when the type is down, exp must be less 0'
            self.down_exp = (self.down_exp * self.down_num + exp) / (self.down_num + 1)
            self.min_down_exp = min(self.min_down_exp, exp)
        else:
            # assert exp >= 0, 'when the type is neutral, exp must be greater than 0 == ' + str(exp)
            self.neutral_num += 1
            self.neutral_exp = (self.neutral_exp * self.neutral_num) / (self.neutral_num + 1)

    @property
    def status(self):
        if self.__status:
            return self.__status
        if self.up_num and self.up_num > self.down_num:
            return 'up'
        elif self.down_num and self.down_num > self.up_num:
            return 'down'
        else:
            return 'no change'

    def exp(self):
        if self.status == 'up':
            return self.up_exp
        elif self.status == 'down':
            return self.down_exp
        else:
            return self.neutral_exp


class Item(object):
    def __init__(self, start=None, end=None, chrom=None, chrom_index=None):
        self.__start = start
        self.__end = end
        self.__chrom = chrom
        self.__chrom_index = chrom_index

        self.trans_set = set()
        self.cis_set = set()

        self.lnc_exp_item = ExpItem('lncrna')
        self.m_exp_item = ExpItem('mrna')

        self.__mrna_num = 0
        self.__lrna_num = 0

        self.__status = None

    def in_item(self, chrom, start, end):
        if chrom == self.__chrom:
            return False
        if self.start <= start <= self.end:
            return True
        else:
            return False

    def is_in(self, start, end):
        if start < self.start:
            return -1
        elif start <= self.end:
            return 0
        else:
            return 1

    def add_exp(self, exp, up_or_down, is_lnc=False):
        if is_lnc:
            self.__lrna_num += 1
            exp_item = self.lnc_exp_item
        else:
            self.__mrna_num += 1
            exp_item = self.m_exp_item
        exp_item.add_exp(exp, up_or_down)

    def add_targets(self, *items, **kwargs):
        is_cis = kwargs.get('is_cis', False)
        tmp_set = self.cis_set if is_cis is True else self.trans_set
        for item in items:
            tmp_set.add(item)

    def in_comm_chrom(self, item):
        return self.chrom == item.chrom

    @property
    def is_lnc(self):
        # return True if self.cis_set or self.trans_set else False
        return True if self.__lrna_num >= self.__mrna_num else False

    def status(self, is_lnc):
        if is_lnc:
            return self.lnc_exp_item.status
        else:
            return self.m_exp_item.status

    def exp(self, is_lnc, up_down, ndigits=2):
        if is_lnc:
            if up_down == 'up':
                return self.lnc_exp_item.up_num, round(self.lnc_exp_item.up_exp, ndigits)
            elif up_down == 'down':
                return self.lnc_exp_item.down_num, round(self.lnc_exp_item.down_exp, ndigits)
            else:
                return self.lnc_exp_item.neutral_num, round(self.lnc_exp_item.neutral_exp, ndigits)
        else:
            if up_down == 'up':
                return self.m_exp_item.up_num, round(self.m_exp_item.up_exp, ndigits)
            elif up_down == 'down':
                return self.m_exp_item.down_num, round(self.m_exp_item.down_exp, ndigits)
            else:
                return self.m_exp_item.neutral_num, round(self.m_exp_item.neutral_exp, ndigits)

    def ceil_exp(self, is_lnc, up_down):
        if is_lnc:
            if up_down == 'up':
                return self.lnc_exp_item.up_num, math.ceil(self.lnc_exp_item.up_exp)
            elif up_down == 'down':
                return self.lnc_exp_item.down_num, math.ceil(self.lnc_exp_item.down_exp)
            else:
                return self.lnc_exp_item.neutral_num, math.ceil(self.lnc_exp_item.neutral_exp)
        else:
            if up_down == 'up':
                return self.m_exp_item.up_num, math.ceil(self.m_exp_item.up_exp)
            elif up_down == 'down':
                return self.m_exp_item.down_num, math.ceil(self.m_exp_item.down_exp)
            else:
                return self.m_exp_item.neutral_num, math.ceil(self.m_exp_item.neutral_exp)

    @property
    def start(self):
        return self.__start

    @start.setter
    def start(self, value):
        self.__start = value

    @property
    def end(self):
        return self.__end

    @end.setter
    def end(self, value):
        self.__end = value

    @property
    def chrom(self):
        return self.__chrom

    @chrom.setter
    def chrom(self, value):
        self.__chrom = value

    @property
    def chrom_index(self):
        return self.__chrom_index

    @chrom_index.setter
    def chrom_index(self, value):
        self.__chrom_index = value


class ChromDict(dict):
    def __init__(self, chr_dic, *args, **kwargs):
        self.bind_obj = None
        if 'bind_obj' in kwargs:
            self.bind_obj = kwargs.pop('bind_obj')
            self.log_flag = 0
        super(ChromDict, self).__init__()
        self.chrom_detail_dict = chr_dic
        self.col_num = 500 if 'col_num' not in kwargs else kwargs['col_num']
        self.step = None
        self.index_dict = defaultdict(list)
        self.seq2index_dict = {}
        self.columns_num = 0
        self._init_items()

    def get_index(self, chrom, start, end):
        index_list = self.index_dict[chrom]
        index_len = len(index_list)
        s_index = 0
        e_index = index_len - 1
        while True:
            if s_index > e_index:
                return None

            mid = (s_index + e_index) // 2
            chrom_obj_key = index_list[mid]
            chrom_obj = self[chrom_obj_key]

            is_in = chrom_obj.is_in(start, end)

            if is_in == 0:
                return chrom_obj
            elif is_in == -1:
                e_index = mid - 1
            else:
                s_index = mid + 1

    def add_map(self, chrom, start, end, seq_id):
        if chrom not in self.index_dict:
            return
        item = self.get_index(chrom, start, end)
        if item is None:
            return
        self.seq2index_dict[seq_id] = item.chrom_index

    def add_exp(self, chrom, start, end, exp, up_neutral_dwon, seq_id=None, is_lnc=False):
        item = self.get_index(chrom, start, end)
        if item is None:
            return
        item.add_exp(exp=exp, up_or_down=up_neutral_dwon, is_lnc=is_lnc)

    def add_targets(self, *target_items, **kwargs):
        target_is_cis = kwargs.get('target_is_cis', False)
        lnc_seq_id = kwargs.get('lnc_seq_id')
        item = self.get_item(lnc_seq_id)
        # self, *items, target_is_cis=True
        if item is None:
            return
        target_items = [j for j in (self.get_item(i) for i in target_items) if j]
        item.add_targets(*target_items, is_cis=target_is_cis)

    def get_item(self, seq_id):
        index = self.seq2index_dict.get(seq_id)
        if index is None:
            return None
        return self[index]

    def _init_items(self):
        """
        {
            'chr1': {
                'chr_len': 11111, 'color': '12,223,34'
            }
        }
        :return:
        """
        sum_len = sum(i['chr_len'] for _, i in self.chrom_detail_dict.items())
        self.step = math.floor(sum_len / self.col_num)
        index = 0
        for chrom, chr_dic in self.chrom_detail_dict.items():
            chr_len = chr_dic['chr_len']
            start = 1
            end = self.step
            self[index] = Item(start=start, end=end, chrom=chrom, chrom_index=index)
            self.index_dict[chrom].append(index)
            index += 1
            while end < chr_len:
                self.columns_num += 1
                start = start + self.step
                end = end + self.step
                self[index] = Item(start=start, end=end, chrom=chrom, chrom_index=index)
                self.index_dict[chrom].append(index)
                index += 1
            chr_dic['chr_len'] = end


if __name__ == '__main__':
    chr_obj = ChromDict({'chr1': {'chr_len': 44400000}, 'chr2': {'chr_len': 37800000}, 'chr3': {'chr_len': 39700000}})
    chroms_info = {'1': '248956422', '2': '242193529', '3': '198295559', '4': '190214555', '5': '181538259',
                   '6': '170805979', '7': '159345973', '8': '145138636', '9': '138394717', '10': '133797422',
                   '11': '135086622', '12': '133275309', '13': '114364328', '14': '107043718', '15': '101991189',
                   '16': '90338345', '17': '83257441', '18': '80373285', '19': '58617616', '20': '64444167',
                   '21': '46709983', '22': '50818468', 'X': '156040895', 'Y': '57227415'}

    with open(r"C:\Users\zhipeng.zhao\Desktop\chroms_gene_detail__1.json") as in_handler:
        dic = json.load(in_handler)

    # for i, d in chr_obj.items():
    #     print(i, d.__dict__)
    # q_val, m_val = divmod(2462380, 24380)
    # chr_index = int(q_val if m_val == 0 else q_val + 1)
    # print(chr_obj.step)