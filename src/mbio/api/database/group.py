# -*- coding: utf-8 -*-
# __author__ = 'xuting'
# lastmodied = 'shenghe'  # 重构的导入方式

from bson import ObjectId
from biocluster.api.database.base import Base, report_check
from collections import OrderedDict
import json

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
        self._project_type = 'meta'
        task_id = "_".join(self.bind_object.sheet.id.split('_')[0:2])  # get task_id by guhaidong 20171115
        self.otu_id = []  # 存放otu文件
        collection = self.db['sg_otu']
        result = collection.find({"task_id": task_id})
        for one_result in result:
            one_otu_id = one_result["_id"]
            self.otu_id.append(one_otu_id)  # get end
            # self._db_name = Config().MONGODB

    @report_check
    def add_ini_group_table(self, file_path, spname_spid=None, task_id=None, group_name=None, add=None,
                            sort_samples=False):
        # group_name /add /  add by zouxuan 20180426
        self.collection = self.db['sg_specimen_group']
        # 解析文件
        group_detail_dict = {}
        group_id_dict = OrderedDict()
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
                    if v:
                        if add is None:
                            v = v
                        else:
                            v = str(add) + str(v)
                        if v in groups_dict[index]:
                            groups_dict[index][v].append(sample)
                        else:
                            groups_dict[index][v] = [sample]
        # 组合样本id
        if not spname_spid:
            self.species = self.db['sg_specimen']
            spname_spid = {}
            for one in self.species.find({'task_id': "_".join(self.bind_object.sheet.id.split('_')[0:2])}):
                spname_spid[one['specimen_name']] = str(one['_id'])
        for group in groups_dict:
            for category in group:
                id_name = []
                for i in group[category]:
                    if i not in spname_spid:
                        self.bind_object.set_error('分组文件中的样本在提供的序列/样本中没有找到', code="51003401")
                    else:
                        id_name.append((str(spname_spid[i]), i))
                group[category] = OrderedDict(id_name)
        for index, name in enumerate(names):
            group_detail_dict[name] = {}
            for gname in groups_dict[index]:
                group_detail_dict[name][gname] = []
                for xx in groups_dict[index][gname]:
                    group_detail_dict[name][gname].append(str(xx))
        # 导入数据
        for index, name in enumerate(names):
            group_id = self.insert_one_group(name, groups_dict[index])
            group_id_dict[name] = str(group_id)
        self.bind_object.logger.info('分组文件中的所有分组方案导入完成。')
        if sort_samples:
            self.bind_object.logger.info("添加分组的样本顺序到样本表中")
            self.update_num_sgotuspecimen(spname_spid, samples)
        return group_id_dict,group_detail_dict

    def update_num_sgotuspecimen(self, spname_spid, samples):
        """
        更新sg_otu_specimen表加入一组编号字段，用于按照本分组文件排序样本。
        :params spname_spid: 样本名对样本id的字典
        :params samples: 样本列表
        """
        sg_otu_specimen = self.db.sg_otu_specimen
        for index, name in enumerate(samples):
            for otu_id in self.otu_id:  # add otu_id loop by guhaidong 20171115
                result = sg_otu_specimen.update_many({'specimen_id': ObjectId(spname_spid[name]), "otu_id": otu_id},
                                                     {"$set": {'order_num': index}})  # add otu_id by guhaidong 20171115
                if result.matched_count < 1:
                    self.bind_object.set_error('没有正确将样本分组中的样本顺序更新到mongo数据表中', code="51003402")
        self.bind_object.logger.info('样本顺序信息更新到数据库中完成。')

    def insert_one_group(self, group_name, group):
        data = {
            'project_sn': self.bind_object.sheet.project_sn,
            'task_id': "_".join(self.bind_object.sheet.id.split('_')[0:2]),
            'group_name': group_name,
            'category_names': group.keys(),
            'specimen_names': group.values()
        }
        group_id = self.collection.insert_one(data).inserted_id
        self.bind_object.logger.info('导入样本分组方案{}成功'.format(group_name))
        return group_id
