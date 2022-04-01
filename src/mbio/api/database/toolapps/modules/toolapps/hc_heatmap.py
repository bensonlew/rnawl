# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.api.database.base import Base, report_check
import re,os
import datetime
from bson import SON
from biocluster.config import Config


class HcHeatmap(Base):
    def __init__(self, bind_object):
        super(HcHeatmap, self).__init__(bind_object)
        self.output_dir = self.bind_object.output_dir
        self.work_dir = self.bind_object.work_dir
        self._project_type = 'toolapps'
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
        self.main_id = self.heatmap_in(self.output_dir)
        self.table_ids = self.table_in()
        return self.main_id

    @report_check
    def run(self):
        """
        运行函数
        """
        self.main_id = self.heatmap_in(self.output_dir)
        self.table_ids = self.table_in()
        return self.main_id
        pass

    def table_in(self):
        """
		导入表格相关信息
		"""
        correlation = self.insert_table(os.path.join(self.output_dir, "heatmap.taxa.table.xls"), '聚类heatmap热图结果表',
                                     '聚类heatmap热图结果表相关数据')
        return correlation

    def insert_table(self, fp, name, desc):
        """
        导入主表和详情表
        """
        with open(fp) as f:
            lines = f.readlines()
            columns = lines[0].strip().split('\t')
            envs = columns[1:]
            insert_data = []
            table_id = self.db['table'].insert_one(SON(
                project_sn=self.bind_object.sheet.project_sn,
                task_id=self.bind_object.id,
                name=name,
                attrs=columns,
                desc=desc,
                status='end',
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )).inserted_id
            for line in lines[1:]:
                line = line.strip().split("\t")
                sp_name = line[0]
                data = {
                    "table_id": table_id,
                    "name": sp_name,
                }
                for n, e in enumerate(envs):
                    data[e] = line[n + 1]
                insert_data.append(data)
            self.db['table_detail'].insert_many(insert_data)
        return table_id

    def heatmap_in(self,dir):
        """
        导入heatmap图相关信息
        """
        self.bind_object.logger.info("聚类热图开始写入主表")
        options = self.bind_object.sheet.options
        self.bind_object.logger.info("options:{}".format(options))
        heatmap_id = self.db['heatmap'].insert_one(SON(
            project_sn=self.bind_object.sheet.project_sn,
            task_id=self.bind_object.id,
            name='heatmap',
            desc='聚类热图',
            status='start',
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        )).inserted_id
        self.bind_object.logger.info("hc_heatmap主表写入主表结束")

        heatmap_table = os.path.join(dir, "heatmap.taxa.table.xls")
        with open (heatmap_table,'r') as f:
            columns = f.readline().strip().split('\t')
            insert_data = []
            species_list = []
            specimen_list = columns[1:]
            for line in f:
                line = line.strip().split("\t")
                data = {
                   "heatmap_id": heatmap_id,
                   "row_name": line[0],
                   "value_type": "pearsons_cor"
                }
                for n, e in enumerate(specimen_list):
                   data[e] = line[n + 1]
                if line[0] not in species_list:
                    species_list.append(line[0])
                insert_data.append(data)
            self.bind_object.logger.info("hc_heatmap_detail表格数据导入开始")
        try:
            self.db['heatmap_detail'].insert_many(insert_data)
        except Exception as e:
            self.bind_object.logger.info("hc_heatmap_detail数据表导入出错{}".format(e))
            self.bind_object.set_error("hc_heatmap_detail数据表导入出错%s", variables=(e), code="54402601")
        else:
            self.bind_object.logger.info("hc_heatmap_detail数据表导入完成")

        row_tree_path = self.output_dir + '/specimen_hcluster.tre'
        if os.path.exists(row_tree_path):
            self.bind_object.logger.info("拥有行聚类树")
            with open(row_tree_path, "r") as m:
                col_tree = m.readline().strip()
                raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', col_tree)
                col_list = [i[1] for i in raw_samp]
        else:
            col_tree = ""
            col_list = specimen_list
        col_tree_path = self.output_dir + '/species_hcluster.tre'
        if os.path.exists(col_tree_path):
            self.bind_object.logger.info("拥有列聚类树")
            with open(col_tree_path, "r") as n:
                row_tree = n.readline().strip()
                raw_samp2 = re.findall(r'([(,]([\[\]\|\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', row_tree)
                row_list = [i[1] for i in raw_samp2]
        else:
            row_tree = ""
            row_list = species_list
        self.bind_object.logger.info("heatmap_id：{}".format(heatmap_id))
        try:
            self.db['heatmap'].update_one({'_id': heatmap_id}, {'$set':
                                                            {'status': 'end',
                                                            'row_tree': row_tree,
                                                            'row_list': row_list,
                                                            'col_tree': col_tree,
                                                            'col_list': col_list}})
        except Exception as e:
            self.bind_object.logger.info("cor_heatmap主表更新出错{}".format(e))
            self.bind_object.set_error("cor_heatmap主表更新出错%s", variables=(e), code="54402602")
        else:
            self.bind_object.logger.info("cor_heatmap主表更新完成")
        return heatmap_id

    def check(self):
        """
        检查文件格式是否正确
        """
        pass
