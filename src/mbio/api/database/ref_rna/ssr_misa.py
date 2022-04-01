# -*- coding: utf-8 -*-
from collections import OrderedDict
import os
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import types
# from biocluster.api..base import Base, report_check
from biocluster.config import Config


# 这个脚本是用来测试导表的，需要先准备好需要使用的文件
#　测试导表的时候，继承不再是Base, 而是object，然后截停函数也不能运行

# class Ssr(Base):
#     def __init__(self, bind_object):
#         super(Ssr, self).__init__(bind_object)
#         self._db_name = Config().MONGODB + 'ref_rna'
class Ssr(object):
    def __init__(self):
        super(Ssr, self).__init__()
        # self.db = pymongo.MongoClient(host="192.168.10.189", port=27017).tsanger_ref_rna
        self.db = Config().MONGODB
        self.db = Config().mongo_client[Config().MONGODB + "_denovo_rna_v2"]

    # @report_check
    def add_ssr(self, ssr_statistic_path, name=None, params=None, stat_id=None):
        # task_id = self.bind_object.sheet.id
        # project_sn = self.bind_object.sheet.project_sn
        # 这个task_id，task_id也是自己随便写
        task_id = "tsg_2000"
        project_sn = "20000000"
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                raise Exception('stat_id必须为ObjectId对象或其对应的字符串！')
        ssr_statistic_dict = OrderedDict()
        with open(ssr_statistic_path, 'r') as f2:
            header2 = f2.readline()
            for line in f2:
                line = line.strip().split('\t')
                ssr_statistic_dict[line[0]] = int(line[1])
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'Ssr_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': params,
            'status': 'end',
            'desc': 'ssr结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'stat_id': stat_id,
            'ssr_statistic_dict':ssr_statistic_dict,
        }
        collection = self.db['sg_ssr']
        print collection
        ssr_id = collection.insert_one(insert_data).inserted_id
        print ("denovo_rna_v2!")
        return ssr_id
        # self.bind_object.logger.info("add denovo_rna_v2_ssr!")
        # print "add ref_ssr!"

    # @report_check
    def add_ssr_analyze_detail(self, ssr_id, ssr_analyze_detail_path):
        if not isinstance(ssr_id, ObjectId):
            if isinstance(ssr_id, types.StringTypes):
                ssr_id = ObjectId(ssr_id)
            else:
                raise Exception('ssr_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(ssr_analyze_detail_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(ssr_analyze_detail_path))
        data_list = []
        with open(ssr_analyze_detail_path, 'r') as f1:
            header1 = f1.readline()
            for line in f1:
                line = line.strip().split('\t')
                data = [
                    ('ssr_id', ssr_id),
                    ('seq_id', line[0]),
                    ('ssr_nr', int(line[1])),
                    ('ssr_type', line[2]),
                    ('ssr_motif', line[3]),
                    ('ssr_motif_length', int(line[4])),
                    ('start', int(line[5])),
                    ('end', int(line[6])),

                    ('fp1', line[7] if line[7] else " "),
                    ('tm1_1', line[8] if line[8] else " "),
                    ('size_1_1', (line[9]) if (line[9]) else " "),
                    ('rp1', line[10] if line[10] else " "),
                    ('tm1_2', line[11] if line[11] else " "),
                    ('size_1_2', line[12] if line[12] else " "),
                    ('product1_size', line[13] if line[13] else " "),
                    ('start_1', line[14] if line[14] else " "),
                    ('end_1', line[15] if line[15] else " "),

                    ('fp2', line[16] if line[16] else " "),
                    ('tm2_1', line[17] if line[17] else " "),
                    ('size_2_1', line[18] if line[18] else " "),
                    ('rp2', line[19] if line[19] else " "),
                    ('tm2_2', line[20] if line[20] else " "),
                    ('size_2_2', line[21] if line[21] else " "),
                    ('product2_size', line[22] if line[22] else " "),
                    ('start_2', line[23] if line[23] else " "),
                    ('end_2', line[24] if line[24] else " "),

                    ('fp3', line[25] if line[25] else " "),
                    ('tm3_1', line[26] if line[26] else " "),
                    ('size_3_1', line[27] if line[27] else " "),
                    ('rp3', line[28] if line[28] else " "),
                    ('tm3_2', line[29] if line[29] else " "),
                    ('size_3_2', line[30] if line[30] else " "),
                    ('product3_size', line[31] if line[31] else " "),
                    ('start_3', line[32] if line[32] else " "),
                    ('end_3', line[33] if line[33] else " "),

                    ('position', line[34])
                ]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_ssr_analyze_detail']
            collection.insert_many(data_list)
        # except Exception, e:
        #     self.bind_object.logger.error("导入注释统计信息：%s出错!" % (ssr_analyze_detail_path, e))
        # else:
        #     self.bind_object.logger.info("导入注释统计信息：%s成功!" % (ssr_analyze_detail_path))
        except Exception, e:
            print("导入注释ssr_analyze_detail信息：%s出错!" % ssr_analyze_detail_path)
        else:
            print("导入注释ssr_analyze_detail信息：%s成功!" % ssr_analyze_detail_path)
    #
    # def add_ssr_statistic(self, ssr_id, ssr_statistic_path):
    #     if not isinstance(ssr_id, ObjectId):
    #         if isinstance(ssr_id, types.StringTypes):
    #             ssr_id = ObjectId(ssr_id)
    #         else:
    #             raise Exception('ssr_id必须为ObjectId对象或其对应的字符串！')
    #     if not os.path.exists(ssr_statistic_path):
    #         raise Exception('{}所指定的路径不存在，请检查！'.format(ssr_statistic_path))
    #     data_list = []
    #     with open(ssr_statistic_path, 'r') as f2:
    #         header2 = f2.readline()
    #         for line in f2:
    #             line = line.strip().split('\t')
    #             # if 'total_number_of_identified_ssrs' == line[0]:
    #             data = [
    #                 ('ssr_id', ssr_id),
    #                 (line[0], int(line[1]))
    #                ]
    #             data = SON(data)
    #             data_list.append(data)
    # 之前上面的这段函数定义是专门写一个细节表,后面觉得这个表格的字段比较少，所以直接插入到主表中
                # elif 'number_of_unigenes_containing_ssrs' == line[0]:
                #     data = [
                #        ('ssr_id', ssr_id),
                #        (line[0], int(line[1]))
                #       ]
                #     data = SON(data)
                #     data_list.append(data)
                #
                # elif 'number_of_unigenes_containing_more_than_1_ssr' == line[0]:
                #     data = [
                #         ('ssr_id', ssr_id),
                #         (line[0], int(line[1]))
                #        ]
                #     data = SON(data)
                #     data_list.append(data)
                #
                # elif 'mono_nucleotide_repeats' == line[0]:
                #     data = [
                #         ('ssr_id', ssr_id),
                #         (line[0], int(line[1]))
                #        ]
                #     data = SON(data)
                #     data_list.append(data)
                #
                # elif 'di_nucleotide_repeats' == line[0]:
                #     data = [
                #         ('ssr_id', ssr_id),
                #         (line[0], int(line[1]))
                #        ]
                #     data = SON(data)
                #     data_list.append(data)
                #
                # elif 'tri_nucleotide_repeats' == line[0]:
                #     data = [
                #         ('ssr_id', ssr_id),
                #         (line[0], int(line[1]))
                #        ]
                #     data = SON(data)
                #     data_list.append(data)
                #
                # elif 'tetra_nucleotide_repeats' == line[0]:
                #     data = [
                #         ('ssr_id', ssr_id),
                #         (line[0], int(line[1]))
                #        ]
                #     data = SON(data)
                #     data_list.append(data)
                #
                #
                # elif 'penta_nucleotide_repeats' == line[0]:
                #     data = [
                #         ('ssr_id', ssr_id),
                #         (line[0], int(line[1]))
                #        ]
                #     data = SON(data)
                #     data_list.append(data)
                #
                # elif 'hexa_nucleotide_repeats' == line[0]:
                #     data = [
                #          ('ssr_id', ssr_id),
                #          (line[0], int(line[1]))
                #         ]
                    # data = SON(data)
                    # data_list.append(data)
        # try:
        #     collection = self.db['sg_ssr_statistic']
        #     collection.insert_many(data_list)
        # # except Exception, e:
        # #     self.bind_object.logger.error("导入注释统计信息：%s出错!" % (ssr_analyze_detail_path, e))
        # # else:
        # #     self.bind_object.logger.info("导入注释统计信息：%s成功!" % (ssr_analyze_detail_path))
        # except Exception, e:
        #     print("导入注释ssr_statistic信息：%s出错!" % ssr_statistic_path)
        # else:
        #     print("导入注释ssr_statistic信息：%s成功!" % ssr_statistic_path)

    def add_ssr_repeats_class(self, ssr_id, ssr_repeats_class_path):
        if not isinstance(ssr_id, ObjectId):
            if isinstance(ssr_id, types.StringTypes):
                ssr_id = ObjectId(ssr_id)
            else:
                raise Exception('ssr_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(ssr_repeats_class_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(ssr_repeats_class_path))
        data_list = []
        with open(ssr_repeats_class_path, 'r') as f3:
            # strip会去除第二行多余的\t，所以第二行会出现超出索引范围
            header3_1 = f3.readline()
            header3_2= f3.readline()
            for line in f3:
                line = line.strip().split('\t')
                # if len(line[0].replace('/', "")) == 2:
                data = [
                    ('ssr_id', ssr_id),
                    ('ssr_type', line[0].replace('/', "")),
                    ('repeat_number_1_5', line[1]),
                    ('repeat_number_6_10', line[2]),
                    ('repeat_number_11_15', line[3]),
                    ('repeat_number_more_than_15', line[4])
                   ]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_ssr_repeats_class']
            collection.insert_many(data_list)
        # except Exception, e:
        #     self.bind_object.logger.error("导入注释统计信息：%s出错!" % (ssr_analyze_detail_path, e))
        # else:
        #     self.bind_object.logger.info("导入注释统计信息：%s成功!" % (ssr_analyze_detail_path))
        except Exception, e:
            print("导入注释ssr_repeats_class信息：%s出错!" % ssr_repeats_class_path)
        else:
            print("导入注释ssr_repeats_class信息：%s成功!" % ssr_repeats_class_path)





if __name__ == '__main__':
    anno = Ssr()
    stat_id = "5915060aa4e1af022575712b"
    ssr_statistic_path = '/mnt/ilustre/users/sanger-dev/workspace/20171030/Single_test_ssr72976/Ssr/ssr_type.txt'
    ssr_id = anno.add_ssr(ssr_statistic_path, name=None, params=None, stat_id=stat_id)
    ssr_analyze_detail_path = '/mnt/ilustre/users/sanger-dev/workspace/20171030/Single_test_ssr72976/Ssr/tmp.txt'
    ssr_repeats_class_path = '/mnt/ilustre/users/sanger-dev/workspace/20171030/Single_test_ssr72976/Ssr/ssr_repeats_class.txt'
    anno.add_ssr_analyze_detail(ssr_id, ssr_analyze_detail_path)
    # anno.add_ssr_statistic(ssr_id, ssr_statistic_path)
    anno.add_ssr_repeats_class(ssr_id, ssr_repeats_class_path)
