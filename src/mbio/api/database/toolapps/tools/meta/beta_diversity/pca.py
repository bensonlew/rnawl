# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
import json
from biocluster.api.database.base import Base, report_check
import re
import datetime
from bson import SON
from biocluster.config import Config
import os


class Pca(Base):
    """
    last_modify by zhouxuan 20170825
    last_modified by binbinzhao 20180821
    """
    def __init__(self, bind_object):
        super(Pca, self).__init__(bind_object)
        self.output_dir = self.bind_object.output_dir
        if Config().MONGODB == 'sanger':
            self._db_name = 'toolapps'
        else:
            self._db_name = 'ttoolapps'
        self._project_type = 'toolapps'
        self.specimen_ids_dict = {}
        self.check()

    @report_check
    def run(self):
        """
        运行函数
        """
        if self.bind_object._task.option('group_table'):
            group_name = os.listdir(self.output_dir)
            for i in group_name:
                self.main_id = self.scatter_in(i)
                self.table_ids = self.table_in(i)
        else:
            self.main_id = self.scatter_in()  # 插入主表
            self.table_ids = self.table_in()  # 要展示在页面的table表的插入
        return self.main_id
        pass

    def table_in(self, group=None):  # 调用函数 insert_table 进行table表的插入
        """
        导入表格相关信息
        """
        if group:
            output_dir = os.path.join(self.output_dir, group)
        else:
            output_dir = self.output_dir
        ratation = self.insert_table(output_dir + '/pca_rotation.xls', 'PCA_特征主成分贡献度表',
                                     'PCA分析中样本特征的贡献度统计结果，例如在OTU表的PCA分析中代表物种/OTU的贡献度', group)
        importance = self.insert_table(output_dir + '/pca_importance.xls', 'PCA_解释度表', 'PCA结果坐标轴的解释度值', group)
        return [ratation, importance]

    def insert_table(self, fp, name, desc, group_name=None):  # 实现table表格插入
        with open(fp) as f:
            columns = f.readline().strip().split('\t')
            insert_data = []
            if 'pca_importance.xls' in fp:
                table_id = self.db['table'].insert_one(SON(
                    project_sn=self.bind_object.sheet.project_sn,
                    task_id=self.bind_object.id,
                    name=name,
                    attrs=columns,
                    desc=desc,
                    status='end',
                    created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    type="importance",
                    group_name=group_name,
                )).inserted_id
            else:
                table_id = self.db['table'].insert_one(SON(
                    project_sn=self.bind_object.sheet.project_sn,
                    task_id=self.bind_object.id,
                    name=name,
                    attrs=columns,
                    desc=desc,
                    status='end',
                    created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    group_name=group_name,
                )).inserted_id
            for line in f:
                line_split = line.strip().split('\t')
                data = dict(zip(columns, line_split))
                data['table_id'] = table_id
                insert_data.append(data)
            self.db['table_detail'].insert_many(insert_data)
        return table_id

    def scatter_in(self, group=None):  # 主表写入
        """
        导入散点图相关信息
        """
        if group:
            output_dir = os.path.join(self.output_dir, group)
        else:
            output_dir = self.output_dir
        with open(output_dir + '/pca_sites.xls') as f:
            attrs = f.readline().strip().split('\t')[1:]  # 获取该二维表的列名
            samples = []
            insert_data = []
            scatter_id = self.db['scatter'].insert_one(SON(
                project_sn=self.bind_object.sheet.project_sn,  # 项目id
                task_id=self.bind_object.id,  # 任务id
                name='pca',  # 分析名称
                attrs=attrs,  # pca_sites.xls文件的列名
                desc='PCA主成分分析',
                status='faild',  # 默认表状态是失败
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),  # 表的创建时间
                group_name=group  # 记录分组名称
            )).inserted_id
            for line in f:
                line_split = line.strip().split('\t')
                samples.append(line_split[0])  # 一行行取样本名
                data = SON(specimen_name=line_split[0],
                           scatter_id=scatter_id)
                data.update(zip(attrs, line_split[1:]))
                insert_data.append(data)
            if self.specimen_ids_dict == {}:
                self.specimen_ids_dict = self.insert_specimens(samples)  # 导入分组原始文件，也就是样本信息
            if group:
                self.insert_group(self.specimen_ids_dict, group, output_dir + "/" + group + '_group.xls')  # 把分组导入到表中去
            self.db['scatter_detail'].insert_many(insert_data)  # 导入散点信息，画图使用
            self.db['scatter'].update_one({'_id': scatter_id}, {'$set': {'status': 'end', 'specimen_ids': self.specimen_ids_dict.values()}})
            # 更新主表为成功
            return scatter_id

    def insert_specimens(self, specimen_names):  # 此处的specimen_names为样本名
        """
        """
        task_id = self.bind_object.id
        project_sn = self.bind_object.sheet.project_sn
        datas = [SON(project_sn=project_sn, task_id=task_id, name=i) for i in specimen_names]
        ids = self.db['specimen'].insert_many(datas).inserted_ids  # 记录样本信息
        return SON(zip(specimen_names, ids))

    def insert_group(self, group_id, group_name, group_path):
        group_dict = {}
        group_cat = []
        with open(group_path, 'r') as g:
            r = g.readline()
            for line in g:
                line = line.strip().split("\t")
                group_dict[line[0]] = line[1]
                if line[1] not in group_cat:
                    group_cat.append(line[1])
        group_list = []
        for i in group_cat:
            dict_ = {}
            for key in group_id:
                if group_dict[key] == i:
                    dict_[str(group_id[key])] = key
            group_list.append(dict_)
        self.db['specimen_group'].insert_one(SON(
            task_id=self.bind_object.id,
            category_names=group_cat,
            specimen_names=group_list,
            group_name=group_name,
        ))

    def check(self):
        """
        检查文件格式是否正确
        """
        pass
