# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

from biocluster.api.database.base import Base, report_check
import re
from collections import defaultdict
from biocluster.config import Config
from collections import OrderedDict


class DenovoGroup(Base):
    def __init__(self, bind_object=None):
        super(DenovoGroup, self).__init__(bind_object)
        self._db_name = Config().MONGODB + '_rna'
        self.scheme = list()
        self.info_dict = dict()
        self.group_id = list()

    @report_check
    def add_ini_group_table(self, file_path, spname_spid, task_id=None):
        self.collection = self.db['sg_denovo_specimen_group']
        if task_id is None:
            task_id = self.bind_object.sheet.id

        (self.info_dict, self.scheme) = self._get_table_info(file_path, spname_spid)
        for s in self.scheme:
            inserted_id = self._insert_one_table(s, task_id, file_path, spname_spid)
            self.group_id.append(inserted_id)
        return self.group_id

    def _insert_one_table(self, s, task_id, file_path, spname_spid):
        category_names = list()
        specimen_names = list()
        for k in self.info_dict:
            if k[0] == s:
                category_names.append(k[1])
        for i in range(len(category_names)):
            for k in self.info_dict:
                if k[0] == s and k[1] == category_names[i]:
                    tmp_dic = dict()
                    for sp in self.info_dict[k]:
                        tmp_dic[str(spname_spid[sp])] = sp
                    specimen_names.append(tmp_dic)
        insert_date = {
            "project_sn": self.bind_object.sheet.project_sn,
            # "project_sn": "test_project_sn",
            "task_id": task_id,
            "group_name": s,
            "category_names": category_names,
            "specimen_names": specimen_names
        }
        try:
            inserted_id = self.collection.insert_one(insert_date).inserted_id
        except Exception as e:
                self.bind_object.logger.error("导入sg_denovo_specimen_group表格{}失败：{}".format(file_path, e))
        else:
            self.bind_object.logger.info("导入sg_denovo_specimen_group表格{}成功".format(file_path))
        return inserted_id

    def _get_table_info(self, file_path, spname_spid):
        info_dic = defaultdict(list)  # info_dict[(分组方案名, 组名)] = [样本名1,  样本名2,...]
        scheme = list()  # 分组方案
        index_gpname = dict()
        with open(file_path, 'rb') as r:
            line = r.next().rstrip("\r\n")
            line = re.split('\t', line)
            for i in range(1, len(line)):
                index_gpname[i] = line[i]
                scheme.append(line[i])
            for line in r:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                for i in range(1, len(line)):
                    if line[0] not in spname_spid:
                        raise Exception("意外错误,样本名{}在以导入的样本当中未找到".format(line[0]))
                    info_dic[(index_gpname[i], line[i])].append(str(line[0]))
        return (info_dic, scheme)

    @report_check
    def get_group_detail(self, file_path, spname_spid, scheme_name):
        (info_dic, scheme) = self._get_table_info(file_path, spname_spid)
        group_detail = {}
        for s in info_dic:
            if s[0] == scheme_name:
                sp_id = list()
                for name in info_dic[s]:
                    sp_id.append(str(spname_spid[name]))
                sp_id.sort()
                group_detail[s[1]] = sp_id
        group_detail = OrderedDict(sorted(group_detail.items(), key=lambda t: t[0]))
        return group_detail
