# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from bson import ObjectId
from biocluster.api.database.base import Base, report_check
from collections import OrderedDict

# from biocluster.config import Config
"""
分组方案格式如下:
{
    "group_name" : "G6",
    "category_names" : [
        "Vad_s",
        "Vad_d",
    ],
    "specimen_names" : [
        {
            "58b4eca103eeb1542dcae5e1" : "2",
            "58b4eca103eeb1542dcae5e3" : "4",
            "58b4eca103eeb1542dcae5e9" : "12",
            "58b4eca103eeb1542dcae5ef" : "18",
            "58b4eca103eeb1542dcae5db" : "20"
        },
        {
            "58b4eca103eeb1542dcae5e0" : "3",
            "58b4eca103eeb1542dcae5e2" : "5",
            "58b4eca103eeb1542dcae5e8" : "11",
            "58b4eca103eeb1542dcae5ee" : "19",
            "58b4eca103eeb1542dcae5dc" : "21"
        },
    ],
    "task_id" : "sanger_10389",
}
"""


class Group(Base):
    def __init__(self, bind_object=None):
        super(Group, self).__init__(bind_object)
        self._project_type = 'metaasv'
        task_id = "_".join(self.bind_object.sheet.id.split('_')[0:2])
        self.otu_id = []  # 存放otu文件
        collection = self.db['asv']
        result = collection.find({"task_id": task_id})
        for one_result in result:
            one_otu_id = one_result["_id"]
            self.otu_id.append(one_otu_id)

    @report_check
    def add_ini_group_table(self, file_path, spname_spid=None, task_id=None, group_name=None, add=None,
                            sort_samples=False):
        self.group_id_list = []
        self.collection = self.db['specimen_group']
        # 解析文件
        with open(file_path) as f:
            names = f.readline().strip()
            if group_name is None:
                names = names.split('\t')[1:]  # 分组方案名称
            else:
                names = group_name
            groups_dict = [OrderedDict() for i in names]  # 各分组方案具体内容，元素为每个分组方案的有序字典
            samples = []  # 样本名列表
            for i in f:
                groups = i.strip().split('\t')
                sample = groups[0]
                samples.append(sample)
                for index, v in enumerate(groups[1:]):
                    if add is None:
                        v = v
                    else:
                        v = str(add) + str(v)
                    if v in groups_dict[index]:
                        groups_dict[index][v].append(sample)
                    else:
                        groups_dict[index][v] = [sample]
        self.bind_object.logger.info('groups_dict: {}'.format(groups_dict))
        # 组合样本id
        if not spname_spid:
            self.species = self.db['specimen_detail']
            self.specimen = self.db['specimen']
            if task_id:
                self.bind_object.logger.info('task_id: {}'.format(task_id))
                result = self.specimen.find_one({'task_id': "_".join(task_id.split('_')[0:2])})
            else:
                result = self.specimen.find_one({'task_id': "_".join(self.bind_object.sheet.id.split('_')[0:2])})
            if result:
                main_id = result["_id"]
            self.bind_object.logger.info('main_id: {}'.format(main_id))
            spname_spid = {}
            for one in self.species.find({'specimen_id': main_id}):
                spname_spid[one['specimen']] = str(one['_id'])
        # else:
        #     self.species = self.db['specimen_detail']
        #     new_spname_spid = {}
        #     for _id in spname_spid:
        #         result = self.species.find_one({'_id': _id})
        #         new_specimen = result["new_specimen"]
        #         new_spname_spid[new_specimen] = str(_id)
        self.bind_object.logger.info('spname_spid: {}'.format(spname_spid))
        for group in groups_dict:
            for category in group:
                id_name = []
                self.bind_object.logger.info('category: {}'.format(category))
                self.bind_object.logger.info('group: {}'.format(group))
                for i in group[category]:
                    if i not in spname_spid:
                        self.bind_object.logger.info('i: {}'.format(i))
                        self.bind_object.set_error('分组文件中的样本在提供的序列/样本中没有找到')
                    else:
                        id_name.append((str(spname_spid[i]), i))
                group[category] = OrderedDict(id_name)
        # 导入数据
        for index, name in enumerate(names):
            main_group_id = self.insert_one_group(name, groups_dict[index])
            if main_group_id not in self.group_id_list:
                self.group_id_list.append(main_group_id)
        self.bind_object.logger.info('分组文件中的所有分组方案导入完成。')
        if sort_samples:
            self.bind_object.logger.info("添加分组的样本顺序到样本表中")
            self.update_num_sgotuspecimen(spname_spid, samples)
        return self.group_id_list

    def update_num_sgotuspecimen(self, spname_spid, samples):
        """
        更新sg_otu_specimen表加入一组编号字段，用于按照本分组文件排序样本。
        :params spname_spid: 样本名对样本id的字典
        :params samples: 样本列表
        """
        sg_otu_specimen = self.db.asv_specimen
        for index, name in enumerate(samples):
            for otu_id in self.otu_id:
                result = sg_otu_specimen.update_many({'specimen_id': ObjectId(spname_spid[name]), "asv_id": otu_id},
                                                     {"$set": {'order_num': index}})
                if result.matched_count < 1:
                    self.bind_object.set_error('没有正确将样本分组中的样本顺序更新到mongo数据表中')
        self.bind_object.logger.info('样本顺序信息更新到数据库中完成。')

    def insert_one_group(self, group_name, group):
        data = {
            'project_sn': self.bind_object.sheet.project_sn,
            'task_id': "_".join(self.bind_object.sheet.id.split('_')[0:2]),
            'group_name': group_name,
            'category_names': group.keys(),
            'specimen_names': group.values(),
            "is_use": "1"
        }
        inserted_id = self.collection.insert_one(data).inserted_id
        self.bind_object.logger.info('导入样本分组方案{}成功'.format(group_name))
        return inserted_id
