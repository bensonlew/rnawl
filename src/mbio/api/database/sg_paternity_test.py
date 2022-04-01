# -*- coding: utf-8 -*-
# __author__ = 'moli.zhou'
import datetime
import random
from biocluster.api.database.base import Base, report_check
import re
from biocluster.config import Config
from pymongo import MongoClient
import gridfs
from mainapp.libs.param_pack import param_pack
from bson import ObjectId


class SgPaternityTest(Base):
    '''
    将亲子鉴定的结果内容存入数据库中
    '''
    def __init__(self, bind_object):
        super(SgPaternityTest, self).__init__(bind_object)
        # self.mongo_client = Config().mongo_client
        # self.database = self.mongo_client[Config().MONGODB+'_paternity_test']
        # self.mongo_client_ref = Config().biodb_mongo_client
        # self.database_ref = self.mongo_client_ref['sanger_paternity_test_ref']  # 正式机的参考库
        self._project_type = "pt"

    @report_check
    def add_sg_father(self,dad,mom,preg,batch_id,member_id,type=None):
        '''
        添加father主表，一个批次有多少个样本就有多少个主表
        :param dad:父本id
        :param mom:母本id
        :param preg:胎儿id
        :param batch_id:批次表的_id
        :param member_id:前端传入的用户id
        :return:返回值为主表的_id
        '''
        temp_d = re.search("WQ([0-9]*)-F.*", dad)
        temp_m = re.search(".*-(M.*)", mom)
        temp_s = re.search(".*-(S.*)", preg)
        if type == 'free':  # 自由交互的主表需要完整记录样本名
            name = 'check_' + dad + '_' + mom + '_' + preg
        else:
            name = dad + "-" + temp_m.group(1) + "-" + temp_s.group(1)
        # 信息增加modify by zhouxuan 20170705

        pt_collection = self.db["sg_pt_customer"]
        # -T 表示重上机信息不变 # modify by zhouxuan 20170728
        # 根据现在的样本名绑定家系信息中的name方便查找家系信息
        if temp_d and temp_m and temp_s:
            if re.match('(.*-T)([0-9])', dad):
                dad_ = ('-').join(dad.split('-')[:-1])
            else:
                dad_ = dad
            if re.match('(.*-T)([0-9])', mom):
                mom_ = ('-').join(mom.split('-')[:-1])
                temp_m_ = re.search(".*-(M.*)", mom_)
            else:
                temp_m_ = temp_m
            message_id = dad_ + "-" + temp_m_.group(1)  # 只有父本和母本的名字
        else:
            message_id = 'free_WQ'
        # 对于自由组合而言，只有正确的家系生成的message_id才会有用

        result = pt_collection.find_one({"name": message_id})  # 自由组合可能没有这个
        if result:
            report_status = result['report_status']
            accept = result['accept_time']
            if report_status == '是':
                report_status = '1'
            else:
                report_status = '0'
            time = accept.split('-')
            accept_time = datetime.datetime(int(time[0]), int(time[1]), int(time[2]), 0, 0)
            report_time = accept_time + datetime.timedelta(days=5)
        else:  # 当自由组合找不到的时候，
            report_status = '0'
            report_time = datetime.datetime.now()
            self.bind_object.logger.info('该家系信息不全，请查看：{}'.format(message_id))
            # 这边不需要raise 可以忍受有错误，一方面前面已经判断过一次了，所以不用再次进行判断，另一方面自由组合不需要一定符合家系信息
        if len(str(report_time.month)) == 1:
            ti = str(report_time.year) + '0' + str(report_time.month)
        else:
            ti = str(report_time.year) + str(report_time.month)
        if len(str(report_time.day)) == 1:
            ti = ti + '0' + str(report_time.day)
        else:
            ti += str(report_time.day)
        insert_data = {
            "dad_id": dad,
            "mom_id": mom,
            "preg_id": preg,
            "family_id": temp_d.group(1),
            "name": name,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "batch_id": ObjectId(batch_id),
            "member_id": member_id,
            "message_id": message_id,
            "report_time": ti,
            "report_status": report_status,
        }
        try:
            collection = self.db['sg_father']
            father_id = collection.insert_one(insert_data).inserted_id
        except Exception as e:
            self.bind_object.logger.error('导入家系主表出错：{}'.format(e))
        else:
            self.bind_object.logger.info("导入家系主表成功")
        return father_id

    def add_father_result(self,father_id,pt_father_id, dad_id):
        '''
        将最终的分析匹配结果和交互表id添加到主表当中去，方便网页展示时取数据以及筛选数据
        :param father_id:father主表的_id
        :param pt_father_id:pt_father交互表的_id
        :param dad_id:父本id
        '''
        self.bind_object.logger.info("father主表更新1")
        collection_result = self.db['sg_pt_father_analysis']
        collection = self.db['sg_father']
        # case = collection.find_one({"_id":father_id})
        # dad_id = case['dad_id']
        if collection_result.find_one({'dad_id':dad_id}):
            self.bind_object.logger.info("dad_id存在")
            result_case = collection_result.find_one({'pt_father_id':pt_father_id, "dad_id":dad_id})
        else:
            result_case = collection_result.find_one({'pt_father_id': pt_father_id, "dad_id": 'NA'})
            self.bind_object.logger.info("dad_id为NA")
        result = result_case['result']

        try:
            collection.update({"_id": father_id}, {'$set': {"pt_father_id": pt_father_id,'result':result}}, multi=True)
        except Exception as e:
            self.bind_object.logger.error('更新father主表结果出错：{}'.format(e))
        else:
            self.bind_object.logger.info("更新father主表结果成功")

    def update_infoshow(self, pt_father_id,mom,preg):
        '''
        如果分析结果有问题，如样本深度不够等等。即在结果表中标记qc字段为不合格
        '''
        collection_result = self.db['sg_pt_father_result_info']
        insert={
            "pt_father_id":pt_father_id,
            "qc":"unqualified",
            "mom_id":mom,
            "preg_id":preg
        }

        try:
            collection_result.insert_one(insert)
        except Exception as e:
            self.bind_object.logger.error('更新有问题的母子信息表出错：{}'.format(e))
        else:
            self.bind_object.logger.info("更新有问题的母子信息表成功")

    def add_father_qc(self, father_id, pt_father_id):
        '''
        将分析出的样本qc值加入到主表中，方便页面展示
        '''
        collection_result = self.db['sg_pt_father_result_info']
        collection = self.db['sg_father']

        result_case = collection_result.find_one({'pt_father_id': pt_father_id})
        qc = result_case['qc']

        try:
            collection.update({"_id": father_id}, {'$set': {'qc': qc}}, multi=True)
        except Exception as e:
            self.bind_object.logger.error('更新father主表家系质控出错：{}'.format(e))
        else:
            self.bind_object.logger.info("更新father主表家系质控成功")

    def update_sg_pt_father(self, pt_father_id):
        '''
        流程结束时更新交互主表的状态

        '''
        try:
            collection = self.db['sg_pt_father']
            # collection.update({"_id": pt_father_id}, {'$set': {"status": "end"}})
            collection.update({"_id": pt_father_id}, {'$set': {"status": "end"}}, multi=True)
        except Exception as e:
            self.bind_object.logger.error('更新pt_father主表状态出错：{}'.format(e))
        else:
            self.bind_object.logger.info("更新pt_father主表状态成功")



    # "status": "end",
    # @report_check
    def add_pt_father(self, father_id, err_min, dedup):
        '''
        增加交互主表。第一次运行时自动添加一个（即主表生成时，交互主表也生成）。后在交互页面投递任务时，
        每一个任务对应一个交互主表。每一个主表可能对应不同的交互主表（视交互次数而定，但至少对应一个）
        '''
        params = dict()
        params['err_min'] = err_min
        params['dedup'] = dedup
        name = 'err-' + str(err_min) + '_dedup-'+ str(dedup)
        insert_data = {
            "father_id": father_id,
            "name": name,
            "status": "start"
        }

        collection = self.db['sg_pt_father']
        new_params = param_pack(params)
        insert_data["params"] = new_params
        # collection.insert_data["params"] = params
        try:
            pt_father_id = collection.insert_one(insert_data).inserted_id
            # collection.insert_one(insert_data)
        except Exception as e:
            self.bind_object.logger.error('导入交互主表出错：{}'.format(e))
        else:
            self.bind_object.logger.info("导入交互主表成功")
        return pt_father_id

    @report_check
    def add_sg_ref_file(self,father_id, ref_fasta,targets_bedfile,ref_point,fastq_path):
        '''
        参考文件的记录

        '''
        insert_data={
            "father_id": father_id,
            "ref_fasta": ref_fasta,
            "targets_bedfile": targets_bedfile,
            "ref_point": ref_point,
            "fastq_path":fastq_path,
        }
        try:
            collection = self.db['sg_pt_ref_file']
            collection.insert_one(insert_data)
        # collection.insert_one(insert_data)
        except Exception as e:
            self.bind_object.logger.error('导入参考文件表出错：{}'.format(e))
        else:
            self.bind_object.logger.info("导入参考文件表成功")

    @report_check
    def add_sg_pt_father_detail(self,file_path,pt_father_id):
        '''
        调试表的导入
        '''
        sg_pt_family_detail = list()
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                if line[0] == "chrom":
                    continue
                # if line[44] == 'Mis':
                #     Mis = '错配'
                # else:
                #     Mis = '-'
                if line[8] == 'NA':
                    dad_rf = 'NA'
                else:
                    dad_rf = round(float(line[8]),8)
                if line[17] == 'NA':
                    preg_rf = 'NA'
                else:
                    preg_rf = round(float(line[17]),8)
                if line[26] == 'NA':
                    mom_rf = 'NA'
                else:
                    mom_rf = round(float(line[26]),8)

                insert_data = {
                    # "task_id": self.bind_object.id,
                    "pt_father_id": pt_father_id,
                    "chrom": line[0],
                    "pos":line[1],
                    "dad_id": line[2],
                    "dad_ref": line[3],
                    "dad_alt": line[4],
                    "dad_dp": line[5],
                    "dad_ref_dp": line[6],
                    "dad_alt_dp": line[7],
                    "dad_rf": dad_rf,
                    "dad_geno": line[9],
                    "dad_geno_bases": line[10],
                    "preg_id": line[11],
                    "preg_ref": line[12],
                    "preg_alt": line[13],
                    "preg_dp": line[14],
                    "preg_ref_dp": line[15],
                    "preg_alt_dp": line[16],
                    "preg_rf": preg_rf,
                    "preg_geno": line[18],
                    "preg_geno_bases": line[19],
                    "mom_id": line[20],
                    "mom_ref": line[21],
                    "mom_alt": line[22],
                    "mom_dp": line[23],
                    "mom_ref_dp": line[24],
                    "mom_alt_dp": line[25],
                    "mom_rf": mom_rf,
                    "mom_geno": line[27],
                    "mom_geno_bases": line[28],
                    "reg": line[29],
                    "from": line[30],
                    "to": line[31],
                    "rs": line[32],
                    "hapmap_rf": line[33],
                    "hapmap_geno": line[34],
                    "n": line[35],
                    "mj_ref": line[36],
                    "pA": line[37],
                    "pG": line[38],
                    "pC": line[39],
                    "pT": line[40],
                    "mj_dp": line[41],
                    "mj_gene": line[42],
                    "is_test": line[43],
                    "is_mis": line[44],
                    "mustbe": line[45],
                    "mustnotbe": line[46],
                    "good": line[47],
                    "pi": line[48]
                }
                sg_pt_family_detail.append(insert_data)
            try:
                collection = self.db['sg_pt_father_detail']
                collection.insert_many(sg_pt_family_detail)
            except Exception as e:
                self.bind_object.logger.error('导入调试页面表格出错：{}'.format(e))
            else:
                self.bind_object.logger.info("导入调试页面表格成功")

    @report_check
    def add_pt_father_figure(self, file_dir,pt_father_id):
        '''
        导入结果图片
        '''
        fs = gridfs.GridFS(self.db)
        family_fig = fs.put(open(file_dir + '_family.png', 'r'))
        figure1 = fs.put(open(file_dir + '_fig1.png', 'r'))
        figure2 = fs.put(open(file_dir + '_fig2.png', 'r'))
        preg_percent = fs.put(open(file_dir + '_preg_percent.png', 'r'))
        update_data = {
            # "task_id": self.bind_object.id,
            "pt_father_id": pt_father_id,
            'family_fig': family_fig,
            'figure1': figure1,
            'figure2': figure2,
            'preg_percent': preg_percent
        }
        try:
            collection = self.db['sg_pt_father_figure']
            figure_id = collection.insert_one(update_data).inserted_id
        except Exception as e:
            self.bind_object.logger.error('导入图片表格出错：{}'.format(e))
        else:
            self.bind_object.logger.info("导入图片表格成功")
        return figure_id

    @report_check
    def add_analysis_tab(self, file_path,pt_father_id):
        '''
        结果信息存入表格，包括测试位点数，有效率无效率等等
        '''
        sg_pt_family_detail = list()
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                if line[0] == "dad.id":
                    continue
                temp_fp = eval(line[4])
                RCP = temp_fp / (temp_fp + 1)
                if RCP > 0.5:
                    rcp_result = "> 99.99%"
                else:
                    rcp_result = "< 0.01%"
                insert_data = {
                    # "task_id": self.bind_object.id,
                    "pt_father_id": pt_father_id,
                    "dad_id": line[0],
                    "test_pos_n": line[1],
                    "err_pos_n": line[2],
                    "err_rate": line[3],
                    "fq": line[4],
                    "dp": line[5],
                    "eff_rate": line[6],
                    "ineff_rate": line[7],
                    "result": line[8],
                    "rcp": rcp_result,
                }
                sg_pt_family_detail.append(insert_data)
            try:
                collection = self.db['sg_pt_father_analysis']
                collection.insert_many(sg_pt_family_detail)
            except Exception as e:
                self.bind_object.logger.error('导入是否匹配表格出错：{}'.format(e))
            else:
                self.bind_object.logger.info("导入是否匹配表格成功")

    @report_check
    def add_info_detail(self, file_path,pt_father_id):
        '''
        基本信息存入数据库，包括母本胎儿是否匹配，胎儿信号比例等等
        '''
        sg_pt_family_detail = list()
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                if line[0] == "bed.preg.id":
                    continue
                if line[1] >= 30 and line[0] >= 4 and line[7] >= 95:
                    qc = 'qualified'
                else:
                    qc = 'unqualified'
                if str(line[7]) == 'NA':
                    mom_preg = line[7]
                else:
                    if eval(line[7]) >= 95:
                        mom_preg = '{}  Yes'.format(line[7])
                    else:
                        mom_preg = '{}  No'.format(line[7])
                insert_data = {
                    # "task_id": self.bind_object.id,
                    "pt_father_id": pt_father_id,
                    "preg_id": line[0],
                    "dp_preg": line[1],
                    "percent": line[2],
                    "error": line[3],
                    "s_signal": line[4],
                    "mom_id": line[5],
                    "dp_mom": line[6],
                    "mom_preg": mom_preg,
                    "qc": qc
                }
                sg_pt_family_detail.append(insert_data)
            try:
                collection = self.db['sg_pt_father_result_info']
                collection.insert_many(sg_pt_family_detail)
            except Exception as e:
                self.bind_object.logger.error('导入基本信息表格出错：{}'.format(e))
            else:
                self.bind_object.logger.info("导入基本信息表格成功")

    # @report_check
    def add_test_pos(self, file_path, pt_father_id):
        '''
        测试位点信息导入数据库
        '''
        sg_pt_family_detail = list()
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                if line[0] == "检测位点编号":
                    continue
                if line[5] == 'Mis':
                    Mis = '错配'
                else:
                    Mis = '-'
                insert_data = {
                    # "task_id": self.bind_object.id,
                    "pt_father_id": pt_father_id,
                    "test_no": line[0],
                    "chrom": line[1],
                    "dad_geno": line[2],
                    "mom_geno": line[3],
                    "preg_geno": line[4],
                    "is_mis": Mis
                }
                sg_pt_family_detail.append(insert_data)
            try:
                collection = self.db['sg_pt_father_test_pos']
                collection.insert_many(sg_pt_family_detail)
            except Exception as e:
                self.bind_object.logger.error('导入位点信息表格出错：{}'.format(e))
            else:
                self.bind_object.logger.info("导入位点信息表格成功")

    def has_problem(self,pt_father_id,dad_id):
        '''
        如果在分析家系时，有样本质检不过关，此时不绘制结果图，匹配结果字段做异常标记
        '''
        collection = self.db['sg_pt_father_analysis']
        if collection.find_one({'dad_id':dad_id}):
            collection.update({"pt_father_id":pt_father_id,'dad_id':dad_id},{"$set":{"result":'MARK'}}, multi=True)
        else:
            collection.update({"pt_father_id": pt_father_id, 'dad_id': 'NA'}, {"$set": {"result": 'MARK'}}, multi=True)

    def check_pt_message(self, family_id_, member_id_, type):
        collection = self.db["sg_pt_customer"]
        if type == 'mom':
            m = collection.find_one({"pt_serial_number": family_id_, 'mom_id_': member_id_})
        else:
            m = collection.find_one({"pt_serial_number": family_id_, 'dad_id_': member_id_})
        if m:
            return 'True'
        else:
            return 'False'

    def import_dedup_data(self, file_path, pt_father_id):
        sg_pt_family_detail = list()
        with open(file_path, 'r') as f:
            data = f.readlines()[1:]
            for line in data:
                line = line.strip().split('\t')
                temp_fp = eval(line[4])
                RCP = float(temp_fp) / (float(temp_fp) + 1)
                if RCP > 0.5:
                    rcp_result = ">99.99%"
                else:
                    rcp_result = "<0.01%"
                result = self.ref_db['sg_pt_ref_main'].find_one({"sample_id": line[0],
                                                                      "storage_time": {"$exists": True}})  # 正式机
                # result = self.database['sg_pt_ref_main'].find_one({"sample_id": line[0],
                #                                                    "storage_time": {"$exists": True}})   # 测试机
                if result:
                    dad_time = result['storage_time']  # 改样本的入库时间，之前的都没有(在sanger上跑过的才会有)
                else:
                    dad_time = ''
                insert_data = {
                    "pt_father_id": pt_father_id,
                    "dad_id": line[0],
                    "test_pos_n": line[1],
                    "err_pos_n": line[2],
                    "err_rate": line[3],
                    "fq": format(eval(line[4]), '.2e'),  # 科学计数法保留两位
                    "dp": line[5],
                    "eff_rate": line[6],
                    "ineff_rate": line[7],
                    "result": line[8],
                    "rcp": rcp_result,
                    "dad_time": dad_time
                }
                sg_pt_family_detail.append(insert_data)
            try:
                collection = self.db['sg_pt_father_analysis']
                collection.insert_many(sg_pt_family_detail)
            except Exception as e:
                self.bind_object.logger.error('导入查重表格出错：{}'.format(e))
            else:
                self.bind_object.logger.info("导入查重表格成功")

    def sample_size(self, sample_id, batch_id, tab_none=None):
        collection = self.db['sg_pt_problem_sample']
        # self.mongo_client_ref = Config().biodb_mongo_client   # 线上
        # self.database_ref = self.mongo_client_ref['sanger_paternity_test_ref']
        ref_data = self.ref_db['sg_pt_ref_main'].find_one({"sample_id": sample_id})  # 线上
        # ref_data = self.database['sg_pt_ref_main'].find_one({"sample_id": sample_id})  # 线下
        split_data_name = ref_data["split_data_name"]
        try:
            collection.insert_one({'sample_id': sample_id,
                                   'split_data_name': split_data_name,
                                   'batch_id': batch_id,
                                   'time': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                                   'tab_none': tab_none})
        except Exception as e:
            self.bind_object.logger.error('导入问题样本出错：{}'.format(e))
            raise Exception('导入问题样本出错：{}'.format(e))
        else:
            self.bind_object.logger.info("导入问题样本成功")
