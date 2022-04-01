# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
import json
from biocluster.api.database.base import Base, report_check
import re
import datetime
from bson import SON
from biocluster.config import Config
import os


class PhyTree(Base):
    """
    version 1.0
    author zouxuan
    last_modified:20180123
    """

    def __init__(self, bind_object):
        super(PhyTree, self).__init__(bind_object)
        self.output_dir = self.bind_object.output_dir
        self.work_dir = self.bind_object.work_dir
        if Config().MONGODB == 'sanger':
            self._db_name = 'toolapps'
        else:
            self._db_name = 'ttoolapps'
        self._project_type = 'toolapps'
        self.tree_id = ''
        self.check()

    @report_check
    def run(self):
        """
        运行函数
        """
        self.tree_file = self.output_dir + "/phylo_tree.nwk"
        self.tree_id = self.tree_in(self.tree_file)
        if self.bind_object._task.option("abundance_table").is_set:
            self.bar_table = self.output_dir + "/abundance_table.xls"
            self.bar_in(self.bar_table)
            if self.bind_object._task.option("sample_group").is_set:
                self.group_detail = self.bind_object._task.option("sample_group").get_group_detail()
                self.bind_object.logger.info(self.group_detail)
                self.bind_object._task.option("sample_group").get_info()
                sample = self.bind_object._task.option("sample_group").prop['sample_name']
                self.specimen_ids_dict = self.insert_specimens(sample, type='sample')
                for group in self.group_detail:
                    self.insert_group(group, type="sample")
            else:
                sample = self.bind_object._task.option('abundance_table').get_sample_info()
                self.bind_object.logger.info(sample)
                self.specimen_ids_dict = self.insert_specimens(sample, type='sample',group=True)
        if self.bind_object._task.option("otu_group").is_set:
            self.group_detail = self.bind_object._task.option("otu_group").get_group_detail()
            self.bind_object.logger.info(self.group_detail)
            for group in self.group_detail:
                self.insert_group(group, type="sequence")


    @report_check
    def bar_in(self, file_path):
        species = []
        data_list = []
        with open(file_path) as f:
            lines = f.readlines()
            head = lines[0].strip().split("\t")
            head = head[1:]
            bar_id = self.db['bar'].insert_one(SON(
                _id=self.tree_id,
                project_sn=self.bind_object.sheet.project_sn,
                task_id=self.bind_object.id,
                name='bar',
                desc='一个样本一个值的柱形图表格',
                status='faild',
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                categories=head
            )).inserted_id
            for line in lines[1:]:
                line_split = line.strip().split("\t")
                self.bind_object.logger.info(line_split)
                data = [
                    ("bar_id", bar_id),
                    ("sample_name", line_split[0]),
                    ("value", [float(x) for x in line_split[1:]])
                ]
                species.append(line_split[0])
                data_son = SON(data)
                data_list.append(data_son)
            self.db['bar_detail'].insert_many(data_list)
            self.db['bar'].update_one({'_id': bar_id},
                                      {'$set': {'status': 'end', 'attrs': species}})

    def insert_specimens(self, specimen_names, type=None, group=False):
        """
        """
        task_id = self.bind_object.id
        project_sn = self.bind_object.sheet.project_sn
        datas = [SON(project_sn=project_sn, task_id=task_id, name=i, type=type) for i in specimen_names]
        ids = self.db['specimen'].insert_many(datas).inserted_ids
        id = [str(i) for i in ids]
        if group:
            self.db['specimen_group'].insert_one(
                SON(task_id=task_id, category_names=specimen_names, group_name='all', type=type,
                    specimen_names=SON(zip(id,specimen_names))))
        self.bind_object.logger.info("样本id导入结束")
        return SON(zip(specimen_names, ids))

    def insert_group(self, group, type=None):
        category_names = []
        specimen_names = []
        for s1 in self.group_detail[group]:
            category_names.append(s1)
            group_specimen_ids = {}
            for s2 in self.group_detail[group][s1]:
                try:
                    if type == 'sample':
                        group_specimen_ids[str(self.specimen_ids_dict[s2])] = s2
                    if type == 'sequence':
                        group_specimen_ids[str(self.seq_ids_dict[s2])] = s2
                except:
                    raise Exception("分组方案和结果文件里的样本不一致，请检查特征值是否错误")
            specimen_names.append(group_specimen_ids)
        self.db['specimen_group'].insert_one(SON(
            task_id=self.bind_object.id,
            category_names=category_names,
            specimen_names=specimen_names,
            group_name=group,
            type=type
        ))

    def tree_in(self, tree_file, group_name=None):
        """
        导入树文件
        """
        with open(tree_file) as f:
            line = f.readline()
            tree = line.strip()
            raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\\| 0-9a-zA-Z_-]+?):[0-9])', tree)
            # line = re.sub('[)]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):', "):", line)  #去除bootstrap值
            sample_list = [i[1] for i in raw_samp]
            self.seq_ids_dict = self.insert_specimens(sample_list, 'sequence')
            self.bind_object.logger.info(tree)
        try:
            tree_id = self.db['bar_tree'].insert_one(SON(
                project_sn=self.bind_object.sheet.project_sn,
                task_id=self.bind_object.id,
                name='phylogenetic tree',
                desc='进化树',
                status='end',
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                group_name=group_name,
                value=line
            )).inserted_id

            self.db['circle_tree'].insert_one(SON(
                _id=tree_id,
                project_sn=self.bind_object.sheet.project_sn,
                task_id=self.bind_object.id,
                name='phylogenetic tree',
                desc='环形进化树图',
                status='end',
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                group_name=group_name,
                value=line
            ))
        except Exception, e:
            self.bind_object.logger.error("导入tree信息出错:%s" % e)
            raise Exception("导入tree信息出错:%s" % e)
        else:
            self.bind_object.logger.info("导入tree信息成功!")
        return tree_id

    def check(self):
        """
        检查文件格式是否正确
        """
        pass
