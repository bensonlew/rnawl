# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'
import json
from biocluster.api.database.base import Base, report_check
import re
import os
import datetime
from bson import SON
from biocluster.config import Config


class VennTable(Base):
    def __init__(self, bind_object):
        super(VennTable, self).__init__(bind_object)
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
        all_file = os.listdir(self.output_dir)
        if len(all_file) == 2:
            self.main_id = self.venn_in()
            self.table_ids = self.table_in()
        else:
            for i in self.bind_object._task.option("group_table").prop['group_scheme']:
                self.main_id = self.venn_in(i)
                self.table_ids = self.table_in(i)
        return self.main_id
        pass

    def table_in(self, group=None):
        """
        导入表格相关信息
        """
        all_file = os.listdir(self.output_dir)
        if len(all_file) == 2:
            venn_otu = self.insert_table(self.output_dir + '/venn_table.xls', 'Venn数据表', '分组间共有和分组中特有的物种的数量统计', group)
        else:
            venn_otu = self.insert_table(self.output_dir + '/' + group + '_venn_table.xls', 'Venn数据表', '分组间共有和分组中特有的物种的数量统计', group)
        return [venn_otu]

    def insert_table(self, fp, name, desc, group_name=None):
        columns = ['Group_label', 'Coincidence_num']
        insert_data = []
        table_id = self.db['table'].insert_one(SON(
            project_sn=self.bind_object.sheet.project_sn,
            task_id=self.bind_object.id,
            name=name,
            attrs=columns,
            desc=desc,
            status='end',
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            group_name=group_name
        )).inserted_id
        with open(fp) as f:
            for line2 in f:
                line_split = line2.strip().split('\t')
                data = SON(table_id=table_id)
                data['Group_label'] = line_split[0]
                data['Coincidence_num'] = line_split[1]
                insert_data.append(data)
            self.db['table_detail'].insert_many(insert_data)
        return table_id

    def venn_in(self, group=None):
        """
        导入venn图相关信息
        """
        all_file = os.listdir(self.output_dir)
        if len(all_file) == 2:
            venn_otu = self.insert_venn(self.output_dir + '/venn_graph.xls', 'Venn', 'venn图', group)
        else:
            venn_otu = self.insert_venn(self.output_dir + '/' + group + '_venn_graph.xls', 'Venn', 'venn图', group)
        return [venn_otu]

    def insert_venn(self, fp, name, desc, group_name =None):
        with open(fp) as f:
            venn_id = self.db['venn'].insert_one(SON(
                project_sn=self.bind_object.sheet.project_sn,
                task_id=self.bind_object.id,
                name=name,
                desc=desc,
                status='faild',
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                group_name=group_name,
            )).inserted_id
            samples = []
            insert_data = []
            # f = f.next()
            for line in f:
                if re.search("#", line):
                    pass
                else:
                    line_list = []
                    line_split = line.strip().split('\t')
                    samples.append(line_split[0])
                    member = line_split[1].split(",")
                    for i in member:
                        line_list.append(i)
                    data = SON(category_name=line_split[0], venn_id=venn_id)
                    data['venn_list'] = line_list
                    insert_data.append(data)
            self.db['venn_detail'].insert_many(insert_data)
            self.db['venn'].update_one({'_id': venn_id}, {'$set': {'status': 'end', 'attrs': samples}})
            return venn_id

    def check(self):
        """
        检查文件格式是否正确
        """
        pass
