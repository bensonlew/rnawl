# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
# last_modify:20170922
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
# from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from mbio.packages.metagenomic.id_convert import name2id


class AssembleGene(Base):
    # 对应表格链接：http://git.majorbio.com/liu.linmeng/metagenomic/wikis/collection/assemble_gene/assemble
    def __init__(self, bind_object):
        super(AssembleGene, self).__init__(bind_object)
        self._project_type = "metagenomic"
        self.name_convert = name2id(self.bind_object.sheet.id, type='task')
        # self._db_name = Config().MONGODB + '_metagenomic'

    @report_check
    def add_assemble_stat(self, show_method, min_contig, software_ver={}):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        assem_method = "multiple" if "Multiple" in show_method else show_method
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': '组装拼接主表',
            'created_ts': created_ts,
            'name': 'assemble_stat',
            'params': 'null',
            'status': 'end',
            'show_method' : show_method,
            'assem_method': assem_method,
            'min': min_contig,
        }
        if software_ver:
            insert_data.update(software_ver)
        collection = self.db['assemble_stat']
        sequence_id = collection.insert_one(insert_data).inserted_id
        return sequence_id

    @report_check
    def add_assemble_stat_detail(self, sequence_id, stat_path):
        if not isinstance(sequence_id, ObjectId):
            if isinstance(sequence_id, types.StringTypes):
                sequence_id = ObjectId(sequence_id)
            else:
                self.bind_object.set_error('sequence_id必须为ObjectId对象或其对应的字符串！', code="52800201")
        if not os.path.exists(stat_path):  # 检查要上传的数据表路径是否存在
            self.bind_object.set_error('stat_path所指定的路径不存在，请检查！', code="52800202")
        data_list = list()  # 存入表格中的信息，然后用insert_many批量导入
        with open(stat_path, 'rb') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                if line[0] in self.name_convert.keys():
                    specimen = self.name_convert[line[0]]
                else:
                    specimen = line[0]
                data = [
                    ('assem_id', sequence_id),
                    ('specimen_name', specimen),
                    ('contigs', int(line[1])),
                    ('contigs_bases', int(line[2])),
                    ('n50', int(line[3])),
                    ('n90', int(line[4])),
                    ('max', int(line[5])),
                    ('min', int(line[6])),
                ]
                if line[0] not in ['Idba_Mix', 'Megahit_Mix', 'Newbler_Mix']:
                    data.append(('method', 'danpin'))
                else:
                    data.append(('method', 'hunpin'))
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['assemble_stat_detail']
            # 将detail表名称写在这里
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (stat_path, e))
            self.bind_object.set_error("导入拼接信息出错", code="52800203")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % stat_path)

    @report_check
    def add_assemble_stat_bar(self, sequence_id, length_path):  # 拼接预测模块长度分布图detail表
        if not isinstance(sequence_id, ObjectId):
            if isinstance(sequence_id, types.StringTypes):
                sequence_id = ObjectId(sequence_id)
            else:
                self.bind_object.set_error('sequence_id必须为ObjectId对象或其对应的字符串！', code="52800201")
        if not os.path.isdir(length_path):
            self.bind_object.set_error('length_path所指定的路径不存在，请检查！', code="52800204")
        # step_data = dict()
        # data_list = list()
        status = 1
        file_list = os.listdir(length_path)
        for files in file_list:
            step_data = dict()
            spe_step = files.strip().split('.')[0]
            spe = spe_step.strip().split('_step_')[0]
            step = spe_step.strip().split('_step_')[1]
            with open(os.path.join(length_path, files), 'rb') as f:
                lines = f.readlines()
                for line in lines[1:-1]:
                    line = line.strip().split('\t')
                    step_data[line[0]] = int(line[1])
            if spe in self.name_convert.keys():
                    spe = self.name_convert[spe]
            data = [
                ('assem_id', sequence_id),
                ('specimen_name', spe),
                ('step', int(step)),
                ('step_data', step_data),
            ]
            data = SON(data)
            # data_list.append(data)
            try:
                collection = self.db["assemble_stat_bar"]
                collection.insert_one(data)
            except Exception, e:
                status = 0
                self.bind_object.logger.error('导入%s信息出错：%s' % (os.path.join(length_path, files), e))
                self.bind_object.set_error("导入拼接信息出错", code="52800203")
            status *= status
        if status == 1:
            self.bind_object.logger.info('导入%s信息成功！' % length_path)

    @report_check
    def add_predict_gene(self, show_method, min_gene, software_ver={}):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        assem_method = "multiple" if "Multiple" in show_method else show_method
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': '基因预测主表',
            'created_ts': created_ts,
            'name': 'predict_gene',
            'assem_method': assem_method,
            'params': 'null',
            'status': 'end',
            'min': min_gene,
        }
        if software_ver:
            insert_data.update(software_ver)
        collection = self.db['predict_gene']
        sequence_id = collection.insert_one(insert_data).inserted_id
        return sequence_id

    @report_check
    def add_predict_gene_detail(self, sequence_id, stat_path):
        if not isinstance(sequence_id, ObjectId):
            if isinstance(sequence_id, types.StringTypes):
                sequence_id = ObjectId(sequence_id)
            else:
                self.bind_object.set_error('sequence_id必须为ObjectId对象或其对应的字符串！', code="52800201")
        if not os.path.exists(stat_path):
            self.bind_object.set_error('stat_path所指定的路径不存在，请检查！', code="52800202")
        data_list = list()
        with open(stat_path, 'rb') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                if line[0] in self.name_convert.keys():
                    specimen = self.name_convert[line[0]]
                else:
                    specimen = line[0]
                data = [
                    ('predict_gene_id', sequence_id),
                    ('specimen_name', specimen),
                    ('orfs', int(line[1])),
                    ('total_length', int(line[2])),
                    ('average_length', round(float(line[3]), 2)),
                    ('max', int(line[4])),
                    ('min', int(line[5])),
                ]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['predict_gene_detail']
            # 将detail表名称写在这里
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (stat_path, e))
            self.bind_object.set_error("导入基因预测信息出错", code="52800205")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % stat_path)

    @report_check
    def add_predict_gene_bar(self, sequence_id, length_path):  # 拼接预测模块长度分布图detail表
        if not isinstance(sequence_id, ObjectId):
            if isinstance(sequence_id, types.StringTypes):
                sequence_id = ObjectId(sequence_id)
            else:
                self.bind_object.set_error('sequence_id必须为ObjectId对象或其对应的字符串！', code="52800201")
        if not os.path.isdir(length_path):
            self.bind_object.set_error('length_path所指定的路径不存在，请检查！', code="52800204")
        # step_data = dict()
        # data_list = list()
        status = 1
        file_list = os.listdir(length_path)
        for files in file_list:
            step_data = dict()
            # data_list = list()
            spe_step = files.strip().split('.')[0]
            spe = spe_step.strip().split('_step_')[0]
            step = spe_step.strip().split('_step_')[1]
            if spe in self.name_convert.keys():
                    spe = self.name_convert[spe]
            with open(os.path.join(length_path, files), 'rb') as f:
                lines = f.readlines()
                for line in lines[1:-1]:
                    line = line.strip().split('\t')
                    step_data[line[0]] = int(line[1])
            data = [
                ('predict_gene_id', sequence_id),
                ('specimen_name', spe),
                ('step', int(step)),
                ('step_data', step_data),
            ]
            data = SON(data)
            try:
                collection = self.db["predict_gene_bar"]
                collection.insert_one(data)
            except Exception, e:
                status = 0
                self.bind_object.logger.error('导入%s信息出错：%s' % (os.path.join(length_path, files), e))
                self.bind_object.logger.set_error("导入基因预测信息出错", code="52800205")
            status *= status
        if status == 1:
            self.bind_object.logger.info('导入%s信息成功！' % length_path)
        '''
            data_list.append(data)
        try:
            collection = self.db["predict_gene_bar"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error('导入%s信息出错：%s' % (length_path, e))
        else:
            self.bind_object.logger.info('导入%s信息成功！' % length_path)
        '''
    @report_check
    def add_predict_gene_total(self, stat_path):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        if not os.path.exists(stat_path):
            self.bind_object.set_error('stat_path所指定的路径不存在，请检查！', code="52800202")
        data_list = list()
        with open(stat_path, 'rb') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [
                    ('orfs', int(line[1])),
                    ('total_length', int(line[2])),
                    ('average_length', round(float(line[3]),2)),
                ]
                data = SON(data)
                data['project_sn'] = project_sn
                data['task_id'] = task_id
                data['desc'] = '基因预测Total统计'
                data['created_ts'] = created_ts
                data['name'] = 'predict_gene_total'
                data['params'] = 'null'
                data['status'] = 'end'
                if line[0] in ['Total', 'total']:
                    break
        collection = self.db['predict_gene_total']
        collection.insert_one(data)
