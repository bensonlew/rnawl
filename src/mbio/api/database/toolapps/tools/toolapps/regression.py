# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import json
from biocluster.api.database.base import Base, report_check
import re
import datetime
from bson import SON
from biocluster.config import Config
from bson.objectid import ObjectId
import os
from types import StringTypes


class Regression(Base):
    def __init__(self, bind_object):
        super(Regression, self).__init__(bind_object)
        self.output_dir = self.bind_object.output_dir ##结果目录
        self.work_dir = self.bind_object.work_dir ##项目目录
        if Config().MONGODB == 'sanger':   ##定义小工具的mongo数据库名称
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
        if self.bind_object._task.option("group").is_set:
            self.group_detail = self.bind_object._task.option("group").get_group_detail()
            self.bind_object.logger.info(self.group_detail)
        sample_list = self.bind_object._task.option("group").prop['sample_name']
        for group in self.group_detail:
            if self.specimen_ids_dict == {}:
                self.specimen_ids_dict = self.insert_specimens(sample_list)
            self.insert_group(group)
        self.main_id = self.regression_in(self.output_dir)
        return self.main_id

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

    def regression_in(self,dir):
        """
        导入regression图相关信息
        """
        files =os.listdir(dir)
        regression_id = ''
        R_2  =float()
        regression_equation =''
        xmin =float()
        ymin =float()
        xmax =float()
        ymax =float()
        adj_R_2=float()
        p_value=float()
        for file in files:
            if file == 'Regression.data.xls':
                with open (dir+"/"+file,'r') as f:
                    columns = f.readline().strip().split('\t')[1:]
                    insert_data = []
                    self.bind_object.logger.info("regression主表的regression_curve写入")
                    regression_id= self.db['regression'].insert_one(SON(
                        project_sn=self.bind_object.sheet.project_sn,
                        task_id=self.bind_object.id,
                        name='regression_curve',
						group_name=None,
                        desc='回归分析点数据',
                        status='failed',
                        sample_name=columns,
                        created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    )).inserted_id
                    self.bind_object.logger.info("regression主表的regression_curve写入结束")
                    for line in f:
                        line = line.strip().split("\t")
                        data = [("regression_id", regression_id), ("sample_name", line[
                            0]), ("X_value", line[1]), ("Y_value", line[2])]
                        data_son = SON(data)
                        insert_data.append(data_son)
                    self.bind_object.logger.info("regression_curve表格数据导入开始")
                try:
                    self.db['regression_detail'].insert_many(insert_data)
                except Exception as e:
                    self.bind_object.logger.info("regression_curve数据表导入出错{}".format(e))
                    raise Exception("regression_curve数据表导入出错{}".format(e))
                else:
                    self.bind_object.logger.info("regression_curvep数据表导入完成")
            if file == 'Regression.message.xls':
                insert_data = []
                with open(dir + "/"+file, 'rb') as r:
                    i = 0
                    for line in r:
                        if i == 0:
                            i = 1
                        else:
                            line = line.strip('\n')
                            line_data = line.split('\t')
                            line_d = ''
                            k = float(line_data[5])
                            if float('%.3f' % float(line_data[6])) > 0:
                                line_d = 'y=' + str('%.3f' % float(line_data[5])) + 'x+' + str(
                                    '%.3f' % float(line_data[6]))
                            elif float('%.3f' % float(line_data[6])) < 0:
                                line_d = 'y=' + str('%.3f' % float(line_data[5])) + 'x' + str(
                                    '%.3f' % float(line_data[6]))
                            elif float('%.3f' % float(line_data[6])) == 0.000:
                                line_d = 'y=' + str('%.3f' % float(line_data[5])) + 'x'
                            R_2= '%.3f' % float(line_data[7])
                            adj_R_2 = '%.3f' % float(line_data[0])
                            p_value = '%.3f' % float(line_data[9])
                            xmin=line_data[1]
                            xmax=line_data[3]
                            if k > 0:
                                yn = line_data[2]
                                yx = line_data[4]
                                if yn < yx:
                                    ymin = yn
                                    ymax = yx
                            else:
                                yn = line_data[4]
                                yx = line_data[2]
                                if yn < yx:
                                    ymin = yx
                                    ymax = yn
                                else:
                                    ymin = yn
                                    ymax = yx
                            regression_equation=line_d
        try:
            self.db['regression'].update_one({'_id': regression_id}, {'$set':{'status': 'end',
                         'R_2':R_2, 'xmin':xmin,'ymin': ymin, 'xmax': xmax,'ymax': ymax,'regression_equation': regression_equation,"adj_R_2": adj_R_2, "p_value": p_value}})
        except Exception as e:
            self.bind_object.logger.info("regression主表更新出错{}".format(e))
            raise Exception("regression主表更新出错{}".format(e))
        else:
            self.bind_object.logger.info("regression主表更新完成")
        return regression_id


    def check(self):
        """
        检查文件格式是否正确
        """
        pass
