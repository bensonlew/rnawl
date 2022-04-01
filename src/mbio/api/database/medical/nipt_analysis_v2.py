# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
from biocluster.api.database.base import Base, report_check
from mainapp.libs.param_pack import param_pack
from biocluster.config import Config
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
import datetime
import xlrd
import os
import re


class NiptAnalysisV2(Base):
    def __init__(self, bind_object):
        """
        无创产前筛查的所有的导表函数
        __author__ = HONGDONG
        __last_modify__ = 20171210
        :param bind_object:
        """
        super(NiptAnalysisV2, self).__init__(bind_object)
        self._project_type = "nipt_v2"

    def add_zz_result(self, file_path, table_id=None):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                raise Exception("table_id必须为ObjectId对象或者其对应的字符串！")
        data_list = []
        with open(file_path, "rb") as r:
            data1 = r.readlines()[1:]
            for line1 in data1:
                temp1 = line1.rstrip().split("\t")
                if temp1[1] == 'NA':
                    zz_data = 'NA'
                else:
                    zz_data = '%.3f' % eval(temp1[1])
                data = [("interaction_id", table_id), ("sample_id", str(temp1[0])), ("zz", zz_data)]
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["sg_interaction_zzresult"]
            collection.insert_many(data_list)
        except Exception, e:
            raise Exception("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list, table_id

    def update_main(self, main_id, file_path):
        insert = {}
        with open(file_path, "rb") as r:
            data1 = r.readlines()
            for line1 in data1:
                temp1 = line1.rstrip().split("\t")
                if temp1[1] == 'total_zz':
                    insert = {
                        "zz": '%.3f' % eval(temp1[2]),
                        'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    }
                else:
                    pass
        try:
            collection = self.db["sg_main"]
            collection.update({"_id": main_id}, {'$set': insert})
        except Exception, e:
            raise Exception("更新主表出错：{}".format(e))
        else:
            self.bind_object.logger.info("更新主表成功!")

    def get_id(self, name, batch_id):
        """
        查找样本的main_id,一个样本就对应一个main_id
        :param name:
        :param batch_id:
        :return:
        """
        collection = self.db["sg_main"]
        temp = collection.find_one({"sample_id": name, "batch_id": ObjectId(batch_id)})['_id']
        self.bind_object.logger.info("sg_main_id:%s" % temp)
        return temp

    def add_z_result(self, file_path, table_id=None):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                raise Exception("table_id必须为ObjectId对象或其对应的字符串!")
        data_list = []
        with open(file_path, 'rb') as r:
            data1 = r.readlines()[1:]
            for line in data1:
                temp = line.rstrip().split("\t")
                data = [("interaction_id", table_id), ("sample_id", str(temp[0])), ("chr", str(temp[1])),
                        ("cn", str(temp[2])), ("bin", str(temp[3])), ("n", str(temp[4])), ("sd", str(temp[5])),
                        ("mean", str(temp[6])), ("z", str(temp[7]))]
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["sg_interaction_zresult"]
            collection.insert_many(data_list)
        except Exception, e:
            raise Exception("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list

    def add_bed_file(self, file_path):
        """
        用于添加nipt分析shell部分生成的bed文件到mongo表中
        """
        data_list = []
        with open(file_path, 'rb') as r:
            data1 = r.readlines()
            for line in data1:
                temp = line.rstrip().split("\t")
                data = [("chr", str(temp[0])),
                        ("start", str(temp[1])), ("end", str(temp[2])), ("gc", str(temp[3])), ("map", str(temp[4])),
                        ("pn", str(temp[5])), ("reads", str(temp[6])), ("sample_id", str(temp[7]))]
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["sg_bed"]
            collection.insert_many(data_list)
        except Exception, e:
            raise Exception("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息到sg_nipt_bed成功!" % file_path)
        return True

    def check_exist_bed(self, sample_id):
        collection = self.db["sg_bed"]
        if collection.find_one({"sample_id": sample_id}):
            return True
        else:
            return False

    def nipt_customer(self, files):
        try:
            bk = xlrd.open_workbook(files)
            sh = bk.sheet_by_name(u'Report')
        except:
            self.bind_object.logger.info("nipt家系表的表格不符合格式要求")
            raise Exception('nipt家系表的表格不符合格式要求')
        nrows = sh.nrows
        insert = list()
        # 获取各行数据
        report_num_index = ''
        sample_date_index = ''
        patient_name_index = ''
        accpeted_date_index = ''
        gestation_index = ''
        number_index = ''
        pregnancy_index = ''
        ivfet_index = ''
        hospital_index = ''
        doctor_index = ''
        tel_index = ''
        status_index = ''
        final_period_index = ''
        gestation_week_index = ''
        age_index = ''
        type_index = ''
        ws_number_index = ''
        for i in range(0, nrows):
            row_data = sh.row_values(i)
            if i == 0:
                try:
                    report_num_index = row_data.index(u'\u8ba2\u5355\u5185\u90e8\u7f16\u53f7')  # (内部订单编号)报告编号
                    sample_date_index = row_data.index(u'\u91c7\u6837\u65e5\u671f')  # (采样日期)采样日期
                    patient_name_index = row_data.index(u'\u68c0\u6d4b\u4eba\u59d3\u540d')  # (检测人姓名)患者姓名
                    accpeted_date_index = row_data.index(u'\u63a5\u6536\u65e5\u671f')  # (接收日期)收样日期
                    number_index = row_data.index(u'\u8ba2\u5355\u7f16\u53f7')  # (订单编号)样本编号
                    # register_number_index = row_data.index(u'\u4f4f\u9662/\u95e8\u8bca\u53f7')  # 住院号
                    gestation_index = row_data.index(u'\u5b55\u4ea7\u53f2')  # 孕产史
                    final_period_index = row_data.index(u'\u672b\u6b21\u6708\u7ecf')  # 末次月经
                    gestation_week_index = row_data.index(u'\u5b55\u5468')  # 孕周
                    pregnancy_index = row_data.index(u'\u5355\u80ce/\u53cc\u80ce')  # 单双胎
                    ivfet_index = row_data.index(u'IVF-ET\u598a\u5a20')  # IVF-ET妊娠
                    hospital_index = row_data.index(u'\u9001\u68c0\u5355\u4f4d/\u533b\u9662')  # 送检单位、医院
                    doctor_index = row_data.index(u'\u9001\u68c0\u533b\u751f')  # 送检医生
                    tel_index = row_data.index(u'\u60a3\u8005\u8054\u7cfb\u7535\u8bdd')  # (患者联系电话)患者联系方式
                    status_index = row_data.index(u'\u6837\u672c\u72b6\u6001')  # (样本状态)标本状态异常
                    age_index = row_data.index(u'\u5e74\u9f84')  # 年龄
                    type_index = row_data.index(u'\u6807\u672c\u7c7b\u578b')  # (标本类型)样本类型
                    ws_number_index = row_data.index(u'\u8ba2\u5355\u5185\u90e8\u7f16\u53f7')  # (内部订单编号)样本缩写
                except:
                    self.bind_object.logger.info("nipt家系表的表头信息不全")
                    raise Exception('nipt家系表的表头信息不全')
            else:
                para_list = []
                report_num = row_data[report_num_index]
                para_list.append(report_num)
                sample_date = row_data[sample_date_index]
                # if sample_date == 'Nf':
                #     para_list.append(sample_date)
                # elif sample_date == '/':
                #     para_list.append(sample_date)
                if type(sample_date) == float:
                    sample_date_tuple = xlrd.xldate_as_tuple(sample_date, 0)
                    para_list.append(
                        str(sample_date_tuple[0]) + '-' + str(sample_date_tuple[1]) + '-' + str(sample_date_tuple[2]))
                else:
                    para_list.append(sample_date)

                patient_name = row_data[patient_name_index]
                para_list.append(patient_name)
                accpeted_date = row_data[accpeted_date_index]

                # if accpeted_date == 'Nf':
                #     para_list.append(accpeted_date)
                # elif accpeted_date == '/':
                #     para_list.append(accpeted_date)
                if type(accpeted_date) == float:
                    accpeted_date_tuple = xlrd.xldate_as_tuple(accpeted_date, 0)
                    para_list.append(
                        str(accpeted_date_tuple[0]) + '-' + str(accpeted_date_tuple[1]) + '-' + str(
                            accpeted_date_tuple[2]))
                else:
                    para_list.append(accpeted_date)
                number = row_data[number_index]
                para_list.append(number)
                # register_number = row_data[register_number_index]
                # para_list.append(register_number)
                gestation = row_data[gestation_index]
                para_list.append(gestation)
                pregnancy = row_data[pregnancy_index]
                para_list.append(pregnancy)
                ivfet = row_data[ivfet_index]
                para_list.append(ivfet)
                hospital = row_data[hospital_index]
                para_list.append(hospital)
                doctor = row_data[doctor_index]
                para_list.append(doctor)
                tel = row_data[tel_index]
                para_list.append(tel)
                status = row_data[status_index]
                para_list.append(status)

                final_period = row_data[final_period_index]
                if type(final_period) == float:
                    final_period_tuple = xlrd.xldate_as_tuple(final_period, 0)
                    para_list.append(str(final_period_tuple[0]) + '-' + str(final_period_tuple[1]) + '-' + str(
                                final_period_tuple[2]))
                else:
                    para_list.append(final_period)

                gestation_week = row_data[gestation_week_index]
                para_list.append(gestation_week)

                age = row_data[age_index]
                para_list.append(str(age))
                sample_type = row_data[type_index]
                para_list.append(sample_type)
                ws_number = row_data[ws_number_index]
                para_list.append(ws_number)
                collection = self.db['sg_customer']
                if type(para_list[10]) == float:
                    tel = int(para_list[10])
                else:
                    tel = para_list[10]

                if collection.find_one({"report_num": para_list[0]}):
                    continue
                else:
                    insert_data = {
                        "report_num": para_list[0],
                        "sample_date": para_list[1],
                        "patient_name": patient_name,
                        "accpeted_date": para_list[3],
                        "number": para_list[4],
                        # "register_number": para_list[5], # 去除住院门诊号
                        "gestation": '/' if str(para_list[5]) == '孕/产/' else para_list[5],
                        "pregnancy": para_list[6],
                        "IVFET": para_list[7],
                        "hospital": para_list[8],
                        "doctor": para_list[9],
                        "tel": tel,
                        "status": para_list[11],
                        "final_period": para_list[12],
                        "gestation_week": self.set_week_date(para_list[13]),
                        "age": para_list[14],
                        "sample_type": para_list[15],
                        "sample_id": para_list[16]
                    }
                    insert.append(insert_data)
        if len(insert) == 0:
            self.bind_object.logger.info('可能无新的客户信息')
        else:
            try:
                collection = self.db['sg_customer']
                collection.insert_many(insert)
            except Exception as e:
                raise Exception('插入客户信息表出错：{}'.format(e))
            else:
                self.bind_object.logger.info("插入客户信息表成功")

    def set_week_date(self, gestation_week):
        """
        备注提示，xls读进来的数据都是unicode，要注意处理
        :param gestation_week:
        :return:
        """
        if re.search("\+", gestation_week):
            temp = gestation_week.strip().split('+')
            if len(temp[0]) > 2:
                m = temp[0][0:len(temp[0]) - 1]
            else:
                m = temp[0]
            if len(temp[1]) > 1:
                n = temp[1][0:len(temp[1]) - 1]
            else:
                n = temp[1]
            data = m + '+' + n
        else:
            if len(gestation_week) >= 3:
                data = gestation_week[0:2]
            else:
                data = gestation_week
        return data

    def export_bed_file(self, sample, dir_path):
        """
        用于导出bed文件，用于后面计算
        """
        collection = self.db['sg_bed']
        sample_bed = str(sample) + ".bed"
        files = os.path.join(dir_path, sample_bed)
        if os.path.exists(files):
            self.bind_object.logger.info("work_dir中已经存在bed文件！")
            pass
        else:
            search_result = collection.find({"sample_id": str(sample)})
            if search_result.count() != 0:
                self.bind_object.logger.info("mongo表中存在样本的bed文件！")
                final_result = search_result
                files = os.path.join(dir_path, sample + '.bed')
            else:
                raise Exception("没有在数据库中搜到%s" % sample)
            with open(files, "w+") as f:
                for n in final_result:
                    f.write(str(n['chr']) + '\t' + str(n['start']) + '\t' + str(n['end']) +
                            '\t' + str(n['gc']) + '\t' + str(n['map']) + '\t' + str(n['pn']) + '\t' +
                            str(n['reads']) + '\t' + str(n['sample_id']) + '\n')
            if os.path.getsize(files):
                return files
            else:
                raise Exception("样本 %s 的bed文件为空！" % sample)

    def add_main(self, member_id, sample_id, batch_id):
        collection = self.db['sg_main']
        insert_data = {
            'member_id': member_id,
            'sample_id': sample_id,
            'batch_id': ObjectId(batch_id),
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        try:
            main_id = collection.insert_one(insert_data).inserted_id
        except Exception as e:
            raise Exception('插入主表出错：{}'.format(e))
        else:
            self.bind_object.logger.info("插入主表成功")
        return main_id

    def add_main_(self, member_id, sample_id, batch_id):
        collection = self.db['sg_main']
        insert_data = {
            'member_id': member_id,
            'sample_id': sample_id,
            'batch_id': ObjectId(batch_id),
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'zz': '',
            'result': "error"
        }
        try:
            main_id = collection.insert_one(insert_data).inserted_id
        except Exception as e:
            raise Exception('插入异常信息主表出错：{}'.format(e))
        else:
            self.bind_object.logger.info("插入异常信息主表成功")
        return main_id

    def add_interaction(self, main_id, bw, bs, ref_group, sample_id):
        collection = self.db['sg_interaction']
        params = dict()
        params['bw'] = bw
        params['bs'] = bs
        params['ref_group'] = ref_group
        params['main_id'] = str(main_id)
        new_params = param_pack(params)
        name = 'bw-' + str(bw) + '_bs-' + str(bs)+'_ref-'+str(ref_group)
        insert_data = {
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'main_id': main_id,
            'params': new_params,
            'name': name,
            "sample_id": sample_id,
            'status': 'start'
        }
        try:
            interaction_id = collection.insert_one(insert_data).inserted_id
        except Exception as e:
            raise Exception('插入交互表出错：{}'.format(e))
        else:
            self.bind_object.logger.info("插入交互表成功{}-{}!".format(bw, bs))
        return interaction_id

    def update_interaction(self, name, interaction_id):
        collection = self.db['sg_interaction']
        try:
            collection.update({"sample_id": name, "_id": interaction_id}, {'$set': {'status': 'end'}})
        except Exception as e:
            raise Exception('更新交互表出错：{}'.format(e))
        else:
            self.bind_object.logger.info("更新交互表成功")

    def add_fastqc(self, files, batch_id):
        collection = self.db['sg_sample_fastqc']
        insert_data = {}
        file_name = files
        insert_data["batch_id"] = ObjectId(batch_id)
        insert_data['created_ts'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        if re.search(r'.*.map.valid_fastqc.html$', file_name):
            insert_data['fastqc_link'] = file_name
            sample_id = os.path.basename(file_name).split('.')[0]
            insert_data['sample_id'] = sample_id
            insert_data['display_name'] = "map.valid"
        else:
            insert_data['fastqc_link'] = file_name
            sample_id = os.path.basename(file_name).split('.')[0]
            temp = sample_id.split('_')[0]
            insert_data['sample_id'] = temp
            insert_data['display_name'] = "fastq.gz"
        try:
            collection.insert_one(insert_data)
        except Exception as e:
            raise Exception('插入fastqc表出错：{}'.format(e))
        else:
            self.bind_object.logger.info("插入fastqc表成功")

    def add_qc(self, files, file_gc, file_new_gc, batch_id):
        collection = self.db['sg_sample_qc']
        insert = {}
        file_name = files.split('/')[-1]
        sample_name = file_name.split('.')[0]
        insert["sample_id"] = sample_name
        insert["batch_id"] = ObjectId(batch_id)
        insert['created_ts'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open(files, 'rb') as f:
            for line in f:
                line = line.strip()
                line = line.split(': ')
                if re.search('.*num$', line[0]):
                    insert["num"] = line[1]
                elif re.search('.*n_map$', line[0]):
                    insert["n_map"] = line[1]
                elif re.search('.*n_dedup$', line[0]):
                    insert["n_dedup"] = line[1]
                elif re.search('.*valid_reads$', line[0]):
                    insert["valid_reads"] = line[1]
                elif re.search('.*properly_paired$', line[0]):
                    insert["properly_paired"] = line[1]
        insert['mapping_ratio'] = '%.3f' % (float(insert["n_map"]) / float(insert["num"]))
        with open(file_gc, 'rb') as gc:
            for line in gc:
                line = line.strip()
                line = line.split('\t')
                if line[0] == "%GC":
                    insert['GC_old'] = '%.3f' % float(line[1])
                elif line[0] == 'adapt_percent':
                    insert['adapt_percent'] = '%.3f' % float(line[1])
        with open(file_new_gc, 'rb') as new_gc:
            line = new_gc.readlines()[0]
            line = line.strip().split('\t')
            insert['GC'] = line[1]
        try:
            collection.insert_one(insert)
        except Exception as e:
            raise Exception('插入qc表出错：{}'.format(e))
        else:
            self.bind_object.logger.info("插入qc表成功")

    def report_result(self, main_id, files, batch_id):
        insert = {}
        collection = self.db['sg_interaction_result']
        # collection_interaction = self.db['sg_interaction']
        collection_main = self.db['sg_main']
        insert['main_id'] = main_id
        insert['batch_id'] = ObjectId(batch_id)
        insert['created_ts'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        # temp = collection_interaction.find_one({"_id": interaction_id})
        # main_id = temp[u'main_id']
        with open(files, 'r') as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                if line[1] == 'total_zz':
                    insert['sample_id'] = line[0]
                    # if -3 < float(line[2]) < 3:
                    #     insert['result'] = 'normal'
                    # else:
                    #     insert['result'] = 'abnormal'
                elif line[1] == '13':
                    insert['chr_13'] = '%.3f' % float(line[2])
                    if -3 < float(line[2]) < 3:
                        insert['13_result'] = 'low'
                        insert['13_desc'] = ''
                    elif float(line[2]) <= -3:
                        insert['13_result'] = 'high'
                        insert['13_desc'] = '13染色体异常'  # 染色体异常
                    elif float(line[2]) >= 3:
                        insert['13_result'] = 'high'
                        insert['13_desc'] = '13染色体三体综合征'
                elif line[1] == '18':
                    insert['chr_18'] = '%.3f' % float(line[2])
                    if -3 < float(line[2]) < 3:
                        insert['18_result'] = 'low'
                        insert['18_desc'] = ''
                    elif float(line[2]) <= -3:
                        insert['18_result'] = 'high'
                        insert['18_desc'] = '18染色体异常'  # 染色体异常
                    elif float(line[2]) >= 3:
                        insert['18_result'] = 'high'
                        insert['18_desc'] = '18染色体三体综合征'
                elif line[1] == '21':
                    insert['chr_21'] = '%.3f' % float(line[2])
                    if -3 < float(line[2]) < 3:
                        insert['21_result'] = 'low'
                        insert['21_desc'] = ''
                    elif float(line[2]) <= -3:
                        insert['21_result'] = 'high'
                        insert['21_desc'] = '21染色体异常'  # 染色体异常
                    elif float(line[2]) >= 3:
                        insert['21_result'] = 'high'
                        insert['21_desc'] = '21染色体三体综合征'
                elif line[1] == '23':
                    insert['X'] = '%.3f' % float(line[2])
                elif line[1] == '24':
                    insert['Y'] = '%.3f' % float(line[2])
        try:
            collection.insert_one(insert)
            if insert['21_result'] == 'low' and insert['18_result'] == 'low' and insert['13_result'] == 'low':
                collection_main.update({'_id': main_id}, {'$set': {'result': 'normal'}})
            else:
                collection_main.update({'_id': main_id}, {'$set': {'result': 'abnormal'}})
        except Exception as e:
            raise Exception('插入结果表和更新主表出错：{}'.format(e))
        else:
            self.bind_object.logger.info("插入结果表和更新主表成功")

    def add_sample_summary(self, samples):
        """
        用于整合产筛样本的信息数据
        :param samples:
        :return:
        """
        sample_info = []
        for sample_id in samples:
            result = self.db['sg_interaction_result'].find_one({"sample_id": sample_id})
            result1 = self.db['sg_sample_qc'].find_one({"sample_id": sample_id})
            result2 = self.db['sg_interaction_zzresult'].find_one({"sample_id": sample_id})
            if result1 and "created_ts" in result1.keys():
                date = result1['created_ts']
            else:
                date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            insert_data = {
                "sample_name": sample_id,
                "num": str(result1['num']) if result1 and "num" in result1.keys() else "--",
                "mapping_ratio": str(result1['mapping_ratio']) if result1 and "mapping"
                                                                              "_ratio" in result1.keys() else "--",
                "num_valid_reads": str(result1['valid_reads']) if result1 and "valid"
                                                                              "_reads" in result1.keys() else "--",
                "GC_mapped": str(result1['GC']) if result1 and "GC" in result1.keys() else "--",
                "adapt_percent": str(result1['adapt_percent']) if result1 and "adapt_p"
                                                                              "ercent" in result1.keys() else "--",
                "ZZ": str(result2['zz']) if result1 and "zz" in result2.keys() else "--",
                "z_21": result['chr_21'] if result1 and "chr_21" in result.keys() else "--",
                "z_18": result['chr_18'] if result1 and "chr_18" in result.keys() else "--",
                "z_13": result['chr_13'] if result1 and "chr_13" in result.keys() else "--",
                "X": result['X'] if result1 and 'X' in result.keys() else "--",
                "Y": result['Y'] if result1 and 'Y' in result.keys() else "--",
                "date": date
            }
            sample_info.append(insert_data)
        try:
            self.db['sg_sample_summary'].insert_many(sample_info)
        except Exception as e:
            raise Exception("导入样本的汇总信息出错！{}".format(e))
        else:
            self.bind_object.logger.info("导入样本的汇总信息成功！")

    def add_report_summary(self, samples, batch_id):
        """
        用于整理报告的相关数据, 这里面的下机时间是板子的下机时间，要去获取到板号然后进行处理
        :param samples:
        :param batch_id:这里的batch_id是sg_analysis_status的主表id，因为现在的拆分的主表的id是唯一的了，所以这里就进行针对处理
        :return:
        """
        sample_info = []
        db = Config().get_mongo_client(mtype="pt_v2")[Config().get_mongo_dbname("pt_v2")]
        datasplit_id = db['sg_analysis_status'].find_one({"_id": ObjectId(batch_id)})
        if datasplit_id and "batch_id" in datasplit_id.keys():
            batch_id_ = datasplit_id['batch_id']
        else:
            raise Exception("sg_analysis_status中没有查找到batch_id字段！")
        batch_info = db['sg_datasplit'].find_one({"_id": batch_id_})
        if batch_info and 'board_batch' in batch_info.keys():
            down_time_ = str(batch_info['board_batch']).split('_')[0]
            down_time = datetime.datetime(int('20' + down_time_[0:2]), int(down_time_[2:4]),
                                          int(down_time_[4:6])).strftime("%Y-%m-%d")
        else:
            down_time = '--'
        for sample_id in samples:
            sample_full_name = sample_id                                         # 报告信息汇总页面保留-2样本的记录
            if '-' in sample_id:
                sample_id = sample_id.strip().split('-')[0]
            result = self.db['sg_customer'].find_one({"sample_id": sample_id})
            insert_data = {
                "down_time": down_time,
                "number": result['number'] if result and 'number' in result.keys() else '--',
                "sample_id": sample_full_name,
                "patient_name": result['patient_name'] if result and 'patient_name' in result.keys() else '--',
                "sample_date": result['sample_date'] if result and 'sample_date' in result.keys() else '--',
                "accepted_date": result['accepted_date'] if result and 'accepted_date' in result.keys() else '--',
                'sex': "女",
                'age': result['age'] if result and 'age' in result.keys() else '--',
                'gestation': result['gestation'] if result and 'gestation' in result.keys() else '--',
                'final_period': result['final_period'] if result and 'final_period' in result.keys() else '--',
                'gestation_week': result['gestation_week'] if result and 'gestation_week' in result.keys() else '--',
                'pregnancy': result['pregnancy'] if result and 'pregnancy' in result.keys() else '--',
                'IVFET': result['IVFET'] if result and 'IVFET' in result.keys() else '--',
                'tel': result['tel'] if result and 'tel' in result.keys() else '--',
                'hospital': result['hospital'] if result and 'hospital' in result.keys() else '--',
                'doctor': result['doctor'] if result and 'doctor' in result.keys() else '--',
                'type': result['sample_type'] if result and 'sample_type' in result.keys() else '--',
                'status': result['status'] if result and 'status' in result.keys() else '--',
                'others': '--',
                'company_name': result['company_name'] if result and 'company_name' in result else '--',
                'contacts_people': result['contacts_people'] if result and 'contacts_people' in result else '--',
                'report_finish_time': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                'report_send_time': '--'
            }
            sample_info.append(insert_data)
        try:
            self.db['sg_report_summary'].insert_many(sample_info)
        except Exception as e:
            raise Exception("导入报告的汇总信息出错！{}".format(e))
        else:
            self.bind_object.logger.info("导入报告的汇总信息成功！")

    def update_analysis_status(self, status_id, update_type, samples=None):
        """
        设置产筛的进度条，status_id是sg_analysis_status的主表id
        :param status_id:
        :param samples:
        :param update_type: 更新的类型主要有3个，1,2,3   1的时候是正常的更新end_count个数， 2的时候是没有样本进行分析的时候直
        接更新end_count与all_count值一致， 3的时候，是每次运行的时候，初始化下all_count，因为第二次运行的时候all_count 中有一部
        分已经计算完成了
        :return:
        """
        db = Config().get_mongo_client(mtype="pt_v2")[Config().get_mongo_dbname("pt_v2")]
        collection = db['sg_analysis_status']
        try:
            if update_type == "2":
                result = collection.find_one({"_id": ObjectId(status_id)})
                if result:
                    all_counts = result['all_counts']
                    collection.update({'_id': ObjectId(status_id)}, {'$set': {'end_counts': all_counts}}, upsert=True)
            elif update_type == "3":
                if samples:
                    all_counts = int(samples)
                    collection.update({'_id': ObjectId(status_id)}, {'$set': {'all_counts': all_counts}}, upsert=True)
                else:
                    raise Exception("更新类型为3时，必须要有样本的个数！")
            elif update_type == "1":
                collection.update({'_id': ObjectId(status_id)}, {'$inc': {'end_counts': 1}}, upsert=True)
            else:
                raise Exception("更新类型{}不合法！".format(update_type))
        except Exception as e:
            raise Exception("更新产筛进度条失败--类型{}！{}".format(update_type, e))
        else:
            self.bind_object.logger.info("更新产筛进度条成功--类型{}！".format(update_type))
