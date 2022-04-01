# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
import json
from biocluster.api.database.base import Base, report_check
import re
import datetime
from bson import SON
from biocluster.config import Config


class Nmds(Base):
    def __init__(self, bind_object):
        super(Nmds, self).__init__(bind_object)
        self.output_dir = self.bind_object.output_dir
        self.work_dir = self.bind_object.work_dir
        if Config().MONGODB == 'sanger':
            self._db_name = 'toolapps'
        else:
            self._db_name = 'ttoolapps'
        self._project_type = 'toolapps'
        self.check()
        self.group_detail = {}
        self.specimen_ids_dict = {}

    @report_check
    def run(self):
        """
        运行函数
        """
        if self.bind_object._task.option("group_table").is_set:
            self.group_detail = self.bind_object._task.option("group_table").get_group_detail()
            self.bind_object.logger.info(self.group_detail)
            for group in self.group_detail:
                sample_list = self.table_in(group)
                if self.specimen_ids_dict == {}:
                    self.specimen_ids_dict = self.insert_specimens(sample_list)
                self.insert_group(group)
                self.main_id = self.nmds_in(group)
        else:
            sample_list = self.table_in()
            self.specimen_ids_dict = self.insert_specimens(sample_list)
            self.main_id = self.nmds_in()

    def table_in(self, group=None):
        """
		导入表格相关信息
		"""
        sample_list = self.insert_table(self.output_dir + '/Nmds/nmds_sites.xls', 'nmds分析结果表',
                                     '画图的原始数据表', group)
        return sample_list

    def insert_table(self, fp, name, desc, group=None):
        self.bind_object.logger.info('开始导入table表')
        with open(fp) as f:
            columns = f.readline().strip().split('\t')
            if columns[1] == "MDS1":
                columns[1] = "NMDS1"
            if columns[2] == "MDS2":
                columns[2] = "NMDS2"
            insert_data = []
            table_id = self.db['table'].insert_one(SON(
                project_sn=self.bind_object.sheet.project_sn,
                task_id=self.bind_object.id,
                name=name,
                attrs=columns,
                desc=desc,
                status='end',
                group_name=group,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )).inserted_id
            sample_list = []
            for line in f:
                line_split = line.strip().split('\t')
                sample_list.append(line_split[0])
                data = dict(zip(columns, line_split))
                data['table_id'] = table_id
                insert_data.append(data)
            self.db['table_detail'].insert_many(insert_data)
            self.bind_object.logger.info('table表导入结束')
            self.bind_object.logger.info('开始导入样本id')
        return sample_list

    def nmds_in(self, group=None):
        """
        导入nmds相关信息
        """
        with open(self.output_dir + '/Nmds/nmds_stress.xls') as f:
            self.bind_object.logger.info('主表导入')
            for line in f:
                stress = line.strip('\n')
            scatter_id = self.db['scatter'].insert_one(SON(
                project_sn=self.bind_object.sheet.project_sn,
                task_id=self.bind_object.id,
                name='nmds',
                desc='nmds分析',
                status='failed',
                specimen_ids=self.specimen_ids_dict.values(),
                group_name=group,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                stress=stress,
            )).inserted_id
            self.bind_object.logger.info('主表导入结束')
        with open(self.output_dir + '/Nmds/nmds_sites.xls', 'r') as s:
            self.bind_object.logger.info("nmds画图数据表导入开始")
            insert_data = []
            head = s.next().strip('\r\n')  # windows换行符
            head = re.split('\t', head)
            new_head = head[1:]
            if new_head == ["MDS1", "MDS2"]:
                new_head = ["NMDS1", "NMDS2"]
            r_list = []
            for line in s:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                sample_num = line[1:]
                otu_detail = dict()
                otu_detail['scatter_id'] = scatter_id
                otu_detail['specimen_name'] = line[0]
                otu_detail['specimen_id'] = self.specimen_ids_dict[line[0]]
                r_list.append(line[0])  # 后续画图做准备
                for i in range(0, len(sample_num)):
                    otu_detail[new_head[i]] = float(sample_num[i])  # 保证画图时取到的数据是数值型
                insert_data.append(otu_detail)
            try:
                self.db['scatter_detail'].insert_many(insert_data)
                self.db['scatter'].update_one({"_id": scatter_id}, {"$set": {"attrs": new_head, "status": "end"}})
            except Exception as e:
                self.bind_object.logger.info("nmds画图数据导入出错{}".format(e))
                raise Exception("nmds画图数据导入出错{}".format(e))
            else:
                self.bind_object.logger.info("nmds画图数据表导入完成")
        return scatter_id

    def insert_specimens(self, specimen_names):
        """
		"""
        task_id = self.bind_object.id
        project_sn = self.bind_object.sheet.project_sn
        datas = [SON(project_sn=project_sn, task_id=task_id, name=i) for i in specimen_names]
        ids = self.db['specimen'].insert_many(datas).inserted_ids
        self.bind_object.logger.info("样本id导入结束")
        return SON(zip(specimen_names, ids))

    def insert_group(self, group):
        category_names = []
        specimen_names = []
        for s1 in self.group_detail[group]:
            category_names.append(s1)
            group_specimen_ids = {}
            for s2 in self.group_detail[group][s1]:
                try:
                    group_specimen_ids[str(self.specimen_ids_dict[s2])] = s2
                except:
                    raise Exception("分组方案和结果文件里的样本不一致，请检查特征值是否错误")
            specimen_names.append(group_specimen_ids)
        self.db['specimen_group'].insert_one(SON(
            task_id=self.bind_object.id,
            category_names=category_names,
            specimen_names=specimen_names,
            group_name=group
        ))

    def check(self):
        """
        检查文件格式是否正确
        """
        pass
