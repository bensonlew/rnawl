# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
# __last_author__ = 'zengjing'
import json
from biocluster.api.database.base import Base, report_check
import re
import datetime
from bson import SON
from biocluster.config import Config
import os


class HcHeatmap(Base):
    def __init__(self, bind_object):
        super(HcHeatmap, self).__init__(bind_object)
        self.output_dir = self.bind_object.output_dir
        self.work_dir = self.bind_object.work_dir
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
        if self.bind_object._task.option("group_table").is_set and self.bind_object._task.option("group_method") in ["average", "sum", "middle"]:
            group_name = self.bind_object._task.option("group_table").prop["group_scheme"]
            self.bind_object.logger.info(group_name)
            for g in group_name:
                fp = self.output_dir + "/{}_result_data".format(g)
                row_path = self.output_dir + "/{}_row_tre".format(g)
                col_path = self.output_dir + "/{}_col_tre".format(g)
                self.main_id = self.heatmap_in(fp, row_path, col_path, g)
                self.table_ids = self.table_in(g)
        else:
            fp = self.output_dir + "/new_result_data"
            row_path = self.output_dir + "/new_row_tre"
            col_path = self.output_dir + "/new_col_tre"
            self.main_id = self.heatmap_in(fp, row_path, col_path)
            self.table_ids = self.table_in()

    def table_in(self, group=None):
        """
        导入表格相关信息
        """
        if group:
            self.insert_table(fp=self.output_dir + '/{}_result_data'.format(group), name='heatmap图结果表', desc='画heatmap图时使用的数据', group_name=group)
        else:
            self.insert_table(self.output_dir + '/new_result_data', 'heatmap图结果表', '画heatmap图时使用的数据')

    def insert_table(self, fp, name, desc, group_name=None):
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
                group_name=group_name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )).inserted_id
            for line in f:
                line_split = line.strip().split('\t')
                data = dict(zip(columns, line_split))
                data['table_id'] = table_id
                insert_data.append(data)
            self.db['table_detail'].insert_many(insert_data)
        return table_id

    def heatmap_in(self, fp, row_path, col_path, group_name=None,):
        """
        导入heatmap图相关信息
        """
        with open(fp) as f:
            self.bind_object.logger.info("heatmap主表写入")
            heatmap_id = self.db['heatmap'].insert_one(SON(
                project_sn=self.bind_object.sheet.project_sn,
                task_id=self.bind_object.id,
                name='heatmap',
                desc='聚类热图',
                status='failed',
                group_name=group_name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            )).inserted_id
            self.bind_object.logger.info("heatmap主表写入结束")
            insert_data = []
            head = f.next().strip('\r\n')  # windows换行符
            head = re.split('\t', head)
            new_head = head[1:]
            r_list = []
            for line in f:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                sample_num = line[1:]
                otu_detail = dict()
                otu_detail['heatmap_id'] = heatmap_id
                otu_detail['row_name'] = line[0]
                r_list.append(line[0])  # 后续画图做准备
                for i in range(0, len(sample_num)):
                    otu_detail[new_head[i]] = float(sample_num[i])  # 保证画图时取到的数据是数值型
                insert_data.append(otu_detail)
            self.bind_object.logger.info("heatmap表格数据导入开始")
            try:
                self.db['heatmap_detail'].insert_many(insert_data)
            except Exception as e:
                self.bind_object.logger.info("heatmap数据表导入出错{}".format(e))
                raise Exception("heatmap数据表导入出错{}".format(e))
            else:
                self.bind_object.logger.info("heatmap数据表导入完成")
            if os.path.exists(row_path):
                self.bind_object.logger.info("拥有行聚类树")
                with open(row_path, "r") as m:
                    row_tree = m.readline().strip()
                    raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', row_tree)
                    row_list = [i[1] for i in raw_samp]
            else:
                row_tree = ""
                row_list = r_list
            if os.path.exists(col_path):
                self.bind_object.logger.info("拥有列聚类树")
                with open(col_path, "r") as n:
                    col_tree = n.readline().strip()
                    raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', col_tree)
                    col_list = [i[1] for i in raw_samp]
            else:
                col_tree = ""
                col_list = new_head
            self.bind_object.logger.info("heatmap_id：{}".format(heatmap_id))
            try:
                self.db['heatmap'].update_one({'_id': heatmap_id}, {'$set':
                                                            {'status': 'end',
                                                            'row_tree': row_tree,
                                                            'row_list': row_list,
                                                            'col_tree': col_tree,
                                                            'col_list': col_list}})
            except Exception as e:
                self.bind_object.logger.info("heatmap主表更新出错{}".format(e))
                raise Exception("heatmap主表更新出错{}".format(e))
            else:
                self.bind_object.logger.info("heatmap主表更新完成")
            return heatmap_id

    def check(self):
        """
        检查文件格式是否正确
        """
        pass
