# -*- coding: utf-8 -*-
# __author__ = 'zhangpeng'
from biocluster.api.database.base import Base, report_check
import re
import datetime
from bson import SON
from biocluster.config import Config
import os


class Hcluster(Base):
    def __init__(self, bind_object):
        super(Hcluster, self).__init__(bind_object)
        self.output_dir = self.bind_object.output_dir
        self.work_dir = self.bind_object._task.work_dir
        if Config().MONGODB == 'sanger':
            self._db_name = 'toolapps'
        else:
            self._db_name = 'ttoolapps'
        self._project_type = 'toolapps'
        self.check()

    @report_check
    def run(self):
        """
        运行函数
        """
        if self.bind_object._task.option('group_table').is_set:
            group_name_list = os.listdir(self.output_dir)
            for i in group_name_list:
                sample_list = self.hcluster_in(i)
                self.table_ids = self.table_in(i)
        else:
            sample_list = self.hcluster_in()
            self.table_ids = self.table_in()
        specimen_ids_dict = self.insert_specimens(sample_list)
        if self.bind_object._task.option('group_table').is_set:
            for group_name in group_name_list:
                group_path = os.path.join(self.work_dir, group_name + '_group.xls')
                if os.path.exists(group_path):
                    self.insert_group(specimen_ids_dict, group_name, group_path)

    def table_in(self, group_name=None):
        '''
        导入二维表
        :param group_name: 分组方案名称
        :return:
        '''
        if group_name:
            output_dir = os.path.join(self.output_dir, group_name)
        else:
            output_dir = self.output_dir
        ratation = self.insert_table(output_dir + '/data_table', '聚类树原始数据表',
                                     '生成聚类树的原始数据表', group_name)
        return [ratation]

    def insert_table(self, fp, name, desc, group_name):
        with open(fp) as f:
            columns = f.readline().strip().split('\t')
            insert_data = []
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

    def hcluster_in(self, group_name=None):
        """
        导入venn图相关信息
        """
        if group_name:
            output_dir = os.path.join(self.output_dir, group_name)
        else:
            output_dir = self.output_dir
        with open(output_dir + '/hcluster.tre') as f:
            hcluster_id = self.db['tree'].insert_one(SON(
                project_sn=self.bind_object.sheet.project_sn,
                task_id=self.bind_object.id,
                name='hcluster',
                desc='层次聚类树图',
                status='end',
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                group_name=group_name,
            )).inserted_id
            line = f.readline()
            tree = line.strip()
            raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', tree)
            sample_list = [i[1] for i in raw_samp]
            self.bind_object.logger.info(tree)
            try:
                collection = self.db["tree"]
                collection.update_one({"_id": hcluster_id}, {"$set": {"value": line}})
            except Exception, e:
                self.bind_object.logger.error("导入tree信息出错:%s" % e)
                raise Exception("导入tree信息出错:%s" % e)
            else:
                self.bind_object.logger.info("导入tree信息成功!")
        return sample_list

    def insert_specimens(self, specimen_names):  # add by zhouxuan 20170508
        task_id = self.bind_object.id
        project_sn = self.bind_object.sheet.project_sn
        datas = [SON(project_sn=project_sn, task_id=task_id, name=i) for i in specimen_names]
        ids = self.db['specimen'].insert_many(datas).inserted_ids
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
