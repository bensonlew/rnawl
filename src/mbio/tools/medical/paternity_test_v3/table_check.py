# -*- coding: utf-8 -*-
# __author__ = 'chenhongyu'

import shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import os
import re
import codecs
import copy
from chardet.universaldetector import UniversalDetector
import datetime

class TableCheckAgent(Agent):
    def __init__(self, parent=None):
        super(TableCheckAgent, self).__init__(parent)
        options = [
            {'name' : 'online_table', 'type' : 'infile', 'format' : 'paternity_test_V2.sample_information'},
            {'name' : 'offline_table', 'type' : 'infile', 'format' : 'paternity_test_V2.sample_information'},
            {'name' : 'samples', 'type' : 'infile' , 'format' : 'paternity_test_V2.samples_ref'},
            {'name': 'experiment', 'type': 'infile', 'format': 'paternity_test_V2.samples_ref'},
            {'name': 'customer', 'type': 'infile', 'format': 'paternity_test_V2.samples_ref'}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option('online_table').is_set:
            raise OptionError('线上上机表必须提供！')
        return True

    def set_resourse(self):
        '''
        定义所需要资源
        :return:
        '''
        self._cpu = 1
        self._memory = '100M'

    def end(self):
        super(TableCheckAgent, self).end()

class TableCheckTool(Tool):
    def __init__(self, config):
        super(TableCheckTool, self).__init__(config)
        self.header = list()
        self.content = list()
        self.names = dict()                                     # 存储上机表中的样品名称及索引，去重
        self.names_revised = list()                             # 去除p表中的生产返工样品
        self.barcode_dic = dict()                               # 存储上机表中的barcode及索引，去重
        ############################################错误捕捉#######################################################
        self.name_not_filled = list()
        self.barcode_not_filled = list()
        self.name_dup = list()
        self.barcode_dup = list()
        self.samples_dup = list()
        self.name_not_match_1 = list()
        self.name_not_match_2 = list()
        self.name_not_match_3 = list()

    def merge_and_de_dup(self):
        '''
        将两个上机表合并，将样本名称与本地数据库作比较，捕捉重复。
        捕捉内部订单编号和index重复的错误
        在线上上机表中删除生产返工样品对应的样品
        :return:
        '''
        online_table = self.option('online_table').prop['table_reformed']
        offline_table = self.option('offline_table').prop['table_reformed'] if self.option('offline_table').is_set else []
        total_table = online_table + offline_table[1:]


        self.header = total_table[0]
        self.content = total_table[1:]
        sn_index = self.get_index(u'内部订单编号', self.header)
        parent_index = self.get_index(u'亲本', self.header)
        barcode_index = self.get_index(u'index序列', self.header)
        FXLX_index = self.get_index(u'分析类型', self.header)

        temp = copy.copy(self.content)                # 复制一份self.content，当遇到返工样品时，在temp中删除先前失败的样品，
                                                      # 最后将temp赋给self.content

        for i in range(0, len(self.content)):
            sn = self.content[i][sn_index]
            parent = self.content[i][parent_index]
            barcode = self.content[i][barcode_index]
            name = sn + '-' + parent if parent else sn
            FXLX = self.content[i][FXLX_index]

            # 检测样本名称是否重复
            if not sn:
                self.name_not_filled.append(self.content[i])
            else:
                if name not in self.names.keys():
                    self.names[name] = i
                else:
                    self.name_dup.append(self.content[i])

            # 检测barcode是否重复
            if not barcode:
                self.barcode_not_filled.append(self.content[i])
            else:
                if barcode not in self.barcode_dic.keys():
                    self.barcode_dic[barcode] = i
                elif sn == self.content[self.barcode_dic[barcode]][sn_index]:
                    # barcode重复：检查订单内部编号是否重复，若重复，即为生产返工样本，检查通过
                    if FXLX == u'生产返工':                   # 删除生产返工样品所对应的失败样品，先判断下失败样品的位置
                        temp.remove(self.content[self.barcode_dic[barcode]])
                    else:
                        temp.remove(self.content[i])
                else:
                    # 若不重复，即不同样本同一barcode，报错
                    self.barcode_dup.append(self.content[i])
        self.content = temp                                   # 删除p表中的返工样品


        for item in self.content:                            # 重新获得一份样品名列表
            sn = item[sn_index]
            parent = item[parent_index]
            name = sn + '-' + parent if parent else sn
            self.names_revised.append(name)



        ################检测上机表中的样本是否与samples表重复###############################
        names_of_samples_table = self.get_specific_information(self.option('samples').prop['path'], '\"name\"', ',')
        for item in self.names_revised:
            if item in names_of_samples_table:
                self.samples_dup.append(item)

    def check_sample_name(self):
        '''
        检测上机表、批次表、客户信息表中样本名称的一致性
        :return:
        '''
        names_of_experiment_table = self.get_specific_information(self.option('experiment').prop['path'], '样本名称', ',')
        names_of_customer_table  = self.get_specific_information(self.option('customer').prop['path'], '样本名称', ',')
        for item in self.names_revised:
            if item not in names_of_experiment_table:
                self.name_not_match_1.append(item)
            if item not in names_of_customer_table:
                self.name_not_match_2.append(item)

        for item in names_of_experiment_table:
            if item not in names_of_customer_table:
                self.name_not_match_3.append(item)

        return len(names_of_experiment_table)

    def parser_all_mistake(self):
        library_sn_not_match = self.merge_mistake_info('library_sn_not_match')
        library_illegal = self.merge_mistake_info('library_illegal')
        ctDNA_sample_type_error = self.merge_mistake_info('ctDNA_sample_type_error')
        test_sample_without_test_prefix = self.merge_mistake_info('test_sample_without_test_prefix')
        FXLX_error = self.merge_mistake_info('FXLX_error')

        self.merge_and_de_dup()
        samples_count_of_experiment_table = self.check_sample_name()

        self.logger.info('上机表样品个数：{}'.format(len(self.names_revised)))
        self.logger.info('实验批次表样品个数：{}'.format(samples_count_of_experiment_table))

        checker = 0

        if library_illegal:
            self.logger.info('建库类型填写有误：')
            checker += 1
            for item in library_illegal:
                self.logger.info(item)

        if library_sn_not_match:
            self.logger.info('建库类型与内部编号格式不匹配：')
            checker += 1
            for item in library_sn_not_match:
                self.logger.info(item)

        if ctDNA_sample_type_error:
            self.logger.info('ctDNA样品类型没有填写或填写有误：')
            checker += 1
            for item in ctDNA_sample_type_error:
                self.logger.info(item)

        if test_sample_without_test_prefix:
            self.logger.info('测试样本的内部编号必须以TEST或test起始：')
            checker += 1
            for item in test_sample_without_test_prefix:
                self.logger.info(item)

        if FXLX_error:
            self.logger.info('线下表分析类型列为空或填写错误：')
            checker += 1
            for item in FXLX_error:
                self.logger.info(item)

        if self.name_not_filled:
            self.logger.info('内部编号没有填写：')
            checker += 1
            for item in self.name_not_filled:
                self.logger.info(item)

        if self.barcode_not_filled:
            self.logger.info('barcode没有填写：')
            checker += 1
            for item in self.barcode_not_filled:
                self.logger.info(item)

        if self.name_dup:
            self.logger.info('样品名称重复：')
            checker += 1
            for item in self.name_dup:
                self.logger.info(item)

        if self.barcode_dup:
            self.logger.info('barcode重复：')
            checker += 1
            for item in self.barcode_dup:
                self.logger.info(item)

        if self.samples_dup:
            self.logger.info('样品在线下数据库中存在：')
            checker += 1
            for item in self.samples_dup:
                self.logger.info(item)

        if self.name_not_match_1:
            self.logger.info('上机表样品名称与实验批次表中不一致：')
            checker += 1
            for item in self.name_not_match_1:
                self.logger.info(item)

        if self.name_not_match_2:
            self.logger.info('上机表样品名称与客户信息表中不一致：')
            checker += 1
            for item in self.name_not_match_2:
                self.logger.info(item)

        if self.name_not_match_3:
            self.logger.info('实验批次表样品名称与客户信息表中不一致：')
            checker += 1
            for item in self.name_not_match_3:
                self.logger.info(item)

        if checker == 0:
            self.logger.info('该批样品正常')
            self.make_sg_sample_info()
            self.end()
        else:
            self.logger.info('上机表内容有误，详情见以上log')
            self.set_error('上机表内容有误')

    def make_sg_sample_info(self):
        number = self.get_index(u'序号',self.header)
        sample_accept_time = self.get_index(u'收样日期', self.header)
        sample_name = self.get_index(u'样本姓名', self.header)
        library_type = self.get_index(u'建库类型', self.header)
        sample_number = self.get_index(u'样本编号', self.header)
        case_name = self.get_index(u'内部订单编号', self.header)
        sample_id = self.get_index(u'亲本', self.header)
        sample_type = self.get_index(u'样本类型', self.header)
        library_name = self.get_index(u'文库名称', self.header)
        s_sequence_num = self.get_index(u'上机测序量', self.header)
        index = self.get_index(u'index序列', self.header)
        sequence_size = self.get_index(u'片段大小', self.header)
        experimental_result = self.get_index(u'系统实验结果', self.header)
        use_type = self.get_index(u'分析类型',self.header)
        remarks = self.get_index(u'备注',self.header)

        data_list_update = list()
        for data in self.content:
            for i in range(0, len(data)):
                data[i] = data[i].encode('utf-8')

            update_data = {
                'number' : data[number],
                'sample_accept_time' : data[sample_accept_time],
                'sample_name' : data[sample_name],
                'library_type' : data[library_type],
                'sample_number' : data[sample_number],
                'case_name' : data[case_name],
                'sample_id' : data[sample_id],
                'sample_type' : data[sample_type],
                'library_name' : data[library_name],
                's_sequence_num' : data[s_sequence_num],
                'index' : data[index],
                'sequence_size' : data[sequence_size],
                'experimental_result' : data[experimental_result],
                'use_type' : data[use_type],
                'remarks' : data[remarks],
                'insert_time': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                'name' : data[case_name] + data[sample_id]
            }
            data_list_update.append(update_data)


        query_tuple = ('name', 'case_name', 'sample_id')
        self.api.paternity_test_v2.add_to_mongo(data_list_update, 'sg_sample_info', 'update', query_tuple)



    def get_specific_information(self, path, title, sep):
        column_list = list()

        encoding = self.detect_encoding(path)

        with codecs.open(path, 'r', encoding) as FIN:
            header_split = FIN.readline().strip('\n').split(sep)             # 待完善
            index = self.get_index(title, header_split)
            for line in FIN:
                line_split = line.strip('\n').split(sep)
                column_list.append(line_split[index])
        return column_list

    def merge_mistake_info(self, error_name, option_1='online_table', option_2='offline_table'):
        '''
        合并file中捕获的错误
        :return:
        '''
        error_1 = self.option(option_1).prop[error_name]
        error_2 = self.option(option_2).prop[error_name] if self.option(option_2).is_set else []
        error_list = error_1.extend(error_2)
        return error_list

    def get_index(self, title, header_split):
        '''
        根据列名返回对应的索引
        :param title:
        :return:
        '''
        for i in range(0, len(header_split)):
            if header_split[i] == title:
                return i

    def detect_encoding(self, path):
        detector = UniversalDetector()
        with open(path, 'r') as FIN:
            for line in FIN:
                detector.feed(line)
                if detector.done:
                    break
            detector.close()
        return detector.result['encoding']

    def run(self):
        super(TableCheckTool, self).run()
        self.parser_all_mistake()



