# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modified = 20180704


from bson import ObjectId
from biocluster.api.database.base import Base, report_check
from collections import OrderedDict


# from biocluster.config import Config
class SpecimenGroup(Base):
    def __init__(self, bind_object=None):
        super(SpecimenGroup, self).__init__(bind_object)
        self._project_type = "metabolome"

    @report_check
    def add_group_table(self, file_path, task_id=None, is_used="0", origin_name=None, group_name=None):
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
            group_info = []
            for line in f:
                line = line.strip().split('\t')
                sample = line[0]
                samples.append(sample)
                for index, v in enumerate(line[1:]):
                    if v in groups_dict[index]:
                        groups_dict[index][v].append(sample)
                    else:
                        groups_dict[index][v] = [sample]
                        # 组合样本id
        for group in groups_dict:
            for category in group:
                id_name = []
                for i in group[category]:
                    id_name.append(i)
                group[category] = id_name
        # 导入数据
        for index, name in enumerate(names):
            specimen_group_id = self.insert_one_group(origin_name, groups_dict[index], is_used=is_used, task_id=task_id)
            group_info.append({
                "names": name,
                "group_id": specimen_group_id,
                "group_detail": groups_dict[index],
            })
        self.bind_object.logger.info('分组文件中的所有分组方案导入完成。')
        return group_info

    def insert_one_group(self, group_name, group, is_used="0", task_id=None):  #
        data = {
            'project_sn': self.bind_object.sheet.project_sn,
            'task_id': task_id if task_id else self.bind_object.sheet.id,
            'group_name': group_name,
            'category_names': group.keys(),
            'specimen_names': group.values(),
            'is_used': is_used,
        }
        inserted_id = self.collection.insert_one(data).inserted_id
        self.bind_object.logger.info('导入样本分组方案{}成功'.format(group_name))
        return inserted_id

    @report_check
    def add_sample_table(self, file_path, task_id=None):
        self.collection = self.db['specimen']
        project_sn = self.bind_object.sheet.project_sn
        task_id = task_id if task_id else self.bind_object.sheet.id
        with open(file_path) as f:
            f.next()
            for line in f:
                if not line.strip():
                    continue
                line = line.strip().split('\t')
                sample = line[0]
                group = line[1]
                data = {
                    'project_sn': project_sn,
                    'task_id': task_id,
                    "name": sample,
                    "new_name": sample,
                    "group": group,
                    "desc": ""
                }
                result = self.collection.find_one(
                    {"project_sn": project_sn, "task_id": task_id, "name": sample, "group": group})
                if not result:
                    inserted_id = self.collection.insert_one(data).inserted_id
                else:
                    pass
        self.bind_object.logger.info('导入样本详情成功')

    @report_check
    def add_compare_group(self, diff_file_path, group_id, diff_group_name=None, task_id=None):
        self.collection = self.db['specimen_group_compare']
        compare_names = []
        with open(diff_file_path) as f:
            f.next()
            for line in f:
                if not line.strip():
                    continue
                line = line.strip().split('\t')
                control = line[0]
                other = line[1]
                diff_group = other + "_vs_" + control
                compare_names.append(diff_group)
            data = {
                'project_sn': self.bind_object.sheet.project_sn,
                'task_id': task_id if task_id else self.bind_object.sheet.id,
                "compare_names": compare_names,
                "compare_category_name": "all",
                "compare_group_name": diff_group_name,
                "group_id": ObjectId(group_id)
            }
        inserted_id = self.collection.insert_one(data).inserted_id
        self.bind_object.logger.info('导入差异对照组成功')
        #compare_names_str = ";".join(compare_names)
        return inserted_id, compare_names
