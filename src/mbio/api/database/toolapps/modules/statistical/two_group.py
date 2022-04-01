# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
import json
from biocluster.api.database.base import Base, report_check
import re
import datetime
from bson import SON
from biocluster.config import Config
import os
from bson.objectid import ObjectId
from types import StringTypes


class TwoGroup(Base):
    def __init__(self, bind_object):
        super(TwoGroup, self).__init__(bind_object)
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
        self.stat_path = self.bind_object._task.two_group.output_dir + '/' + self.bind_object._task.option("test") + '_result.xls'
        self.boxfile_path = self.bind_object._task.two_group.output_dir + '/' + self.bind_object._task.option("test") + '_boxfile.xls'
        self.ci_path = self.bind_object._task.two_group.output_dir + '/' + self.bind_object._task.option("test") + '_CI.xls'
        self.bar_path = self.bind_object._task.two_group.work_dir + '/' + self.bind_object._task.option(
            "test") + '_plot_group_bar.xls'

        self.table_id = ''
        self.group_name = 'group'

        if not os.path.isfile(self.stat_path):
            raise Exception("找不到报告文件:{}".format(self.stat_path))
        if not os.path.isfile(self.boxfile_path):
            raise Exception("找不到报告文件:{}".format(self.boxfile_path))
        if not os.path.isfile(self.ci_path):
            raise Exception("找不到报告文件:{}".format(self.ci_path))


    @report_check
    def run(self):
        """
        运行函数
        """
        self.group_detail = self.bind_object._task.option("group_file").get_group_detail()
        sample_list = self.bind_object._task.option("group_file").prop['sample_name']
        for group in self.group_detail:
            if self.specimen_ids_dict == {}:
                self.specimen_ids_dict = self.insert_specimens(sample_list)
            self.insert_group(group)
        self.table_id = self.table_in()
        # self.difference_bar_in()
        self.bar_in(self.bar_path)
        self.box_in(self.boxfile_path)


    def table_in(self):
        """
        导入表格相关信息
        """
        table_id = self.insert_table(self.stat_path, [self.ci_path], '组间差异显著性结果表',
                                        '两组间差异性检验')
        return table_id

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
            group_name=self.group_name
        ))

    @report_check
    def insert_table(self, statfile, cifiles, name, desc, posthoc=None, group=None):
        stat_info, sort_list = self.cat_files(statfile, cifiles)
        table_id = self.db['table'].insert_one(SON(
            project_sn=self.bind_object.sheet.project_sn,
            task_id=self.bind_object.id,
            name=name,
            desc=desc,
            group_name=self.group_name,
            status='faild',
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            graphic='bar,box'
        )).inserted_id

        self.db['difference_bar'].insert_one(SON(
            _id=table_id,
            project_sn=self.bind_object.sheet.project_sn,
            task_id=self.bind_object.id,
            name=name,
            desc=desc,
            group_name=self.group_name,
            status='faild',
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            #graphic='bar,box',
            test=self.bind_object._task.option("test"),
            coverage=float(self.bind_object._task.option("coverage")),
            ci=float(self.bind_object._task.option("ci")),
        ))

        data_list = []
        data_list_bar = []
        attrs = []
        alias = []
        for name in sort_list:
            data = [
                ("sample_name", name),
                ("table_id", table_id),
                ("graphic", "bar,box")
            ]

            data_bar = [
                ("sample_name", name),
                ("differencebar_id", table_id),
                ("graphic", "bar,box")
            ]

            # self.bind_object.logger.info("run second for loop")
            attrs = ["graphic", "sample_name"]
            alias = ["graphic", "appellation"]
            for i in stat_info[name].keys():
                if i in ['pvalue', 'qvalue', 'corrected_pvalue'] and stat_info[name][i] == 'NA':
                    if i == 'pvalue':
                        data.append((i, stat_info[name][i]))
                        data_bar.append((i, stat_info[name][i]))
                        # attrs.append(i)
                        # alias.append(i)
                    else:
                        data.append(('corrected_pvalue', stat_info[name][i]))
                        data_bar.append(('corrected_pvalue', stat_info[name][i]))
                        # attrs.append('corrected_pvalue')
                        # alias.append('corrected_pvalue')
                elif re.search(r'_pvalue$', i) and posthoc in ['tukeykramer', 'gameshowell']:
                    data.append((i, stat_info[name][i]))
                    data_bar.append((i, stat_info[name][i]))
                    attrs.append(i)
                    alias.append(i)
                elif i in ['qvalue','corrected_pvalue']:
                    data.append(('corrected_pvalue', float(stat_info[name][i])))
                    data_bar.append(('corrected_pvalue', float(stat_info[name][i])))
                    # attrs.append('corrected_pvalue')
                    # alias.append('corrected_pvalue')
                elif i == 'pvalue':
                    data.append(('pvalue', float(stat_info[name][i])))
                    data_bar.append(('pvalue', float(stat_info[name][i])))
                else:
                    data.append((i, float(stat_info[name][i])))
                    data_bar.append((i, float(stat_info[name][i])))
                    attrs.append(i)
                    alias.append(i)
            alias.append("pvalue")
            attrs.append("pvalue")
            attrs.append('corrected_pvalue')
            alias.append('corrected_pvalue')

            # self.bind_object.logger.info("second loop end")
            data_son = SON(data)
            data_bar_son = SON(data_bar)
            data_list.append(data_son)
            data_list_bar.append(data_bar_son)
        # self.bind_object.logger.info("first loop end")
        try:
            self.db['table_detail'].insert_many(data_list)
            self.db['difference_bar_detail'].insert_many(data_list_bar)
        except Exception, e:
            self.bind_object.logger.error("导入species_difference_check_detail信息出错:%s" % e)
        else:
            self.bind_object.logger.info("导入species_difference_check_detail信息成功!")
        self.bind_object.logger.info(table_id)
        self.db['table'].update({"_id": table_id}, {'$set': {"attrs": attrs,"alias": alias,'status': 'end'}})
        self.db['difference_bar'].update({"_id": table_id}, {'$set': {"attrs": attrs, 'status': 'end'}})
        return table_id

    @report_check
    def bar_in(self, file_path):
        species = []
        data_list = []
        with open(file_path) as f:
            lines = f.readlines()
            head = lines[0].strip().split("\t")
            bar_id = self.db['bar'].insert_one(SON(
                _id=self.table_id,
                project_sn=self.bind_object.sheet.project_sn,
                task_id=self.bind_object.id,
                name='bar',
                desc='一个样本一个值的柱形图表格',
                status='faild',
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                categories=head,
            )).inserted_id
            for line in lines[1:]:
                line_split = line.strip().split("\t")
                self.bind_object.logger.info(line_split)
                data = [
                    ("bar_id", self.table_id),
                    ("sample_name", line_split[0]),
                    ("value", line_split[1:])
                ]
                species.append(line_split[0])
                data_son = SON(data)
                data_list.append(data_son)
            self.db['bar_detail'].insert_many(data_list)
            self.db['bar'].update_one({'_id': self.table_id},
                                      {'$set': {'status': 'end', 'attrs': species}})

    @report_check
    def box_in(self, file_path):
        box_plot_id = self.db['box_plot'].insert_one(SON(
            _id=self.table_id,
            project_sn=self.bind_object.sheet.project_sn,
            task_id=self.bind_object.id,
            name='box_plot',
            desc='箱线图的数据',
            status='faild',
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        )).inserted_id
        data_list = []
        with open(file_path, 'rb') as r:
            l = r.readline().strip("\n")
            group_list = re.findall(r'min\((.*?)\)', l)
            while True:
                line = r.readline().strip('\n')
                if not line:
                    break
                line_data = line.split("\t")
                i = 1
                for name in group_list:
                    sample_name = line_data[0]
                    data = [("box_id", self.table_id), ("sample_name", sample_name)]
                    data.append(("category_name", name))
                    data.append(("box_data", [float(line_data[i]),float(line_data[i + 1]),float(line_data[i + 2]),float(line_data[i + 3]),float(line_data[i + 4])]))
                    i += 5
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            collection = self.db["box_plot_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功！" % file_path)
        self.db['box_plot'].update_one({'_id': self.table_id},
                                      {'$set': {'status': 'end', 'attrs': group_list}})
        return box_plot_id

    def cat_files(self, statfile, cifiles=None):
        """
        将组间差异比较的统计结果文件与posthoc结果文件cat成一个文件
        cifiles为posthoc文件路径的列表
        """
        data = {}
        sort_list = []
        with open(statfile, 'rb') as s:
            head = s.readline().strip('\n').split()
            for line in s:
                line = line.strip('\n').split('\t')
                info = {}
                for i in head:
                    info[i] = line[head.index(i) + 1]
                data[line[0]] = info
                sort_list.append(line[0])
        # self.bind_object.logger.info("check if cifiles or not")
        # self.bind_object.logger.info(cifiles)
        if cifiles:
            for f in cifiles:
                with open(f, 'rb') as r:
                    head = r.readline().strip('\n').split('\t')
                    for line in r:
                        line = line.strip('\n').split('\t')
                        if line[0] in sort_list:
                            info = {}
                            for i in head[1:]:
                                info[i] = line[head.index(i)]
                            data[line[0]].update(info)
                        else:
                            pass
        # self.bind_object.logger.info("len cifiles is : %s" % len(cifiles))
        if cifiles and len(cifiles) > 1:
            for name in data:
                keys = data[name].keys()
                mean = []
                low_ci = []
                ci = []
                for i in keys:
                    if re.search(r'mean$', i):
                        mean.append(float(data[name][i]))
                    elif re.search(r'lowerCI$', i):
                        group = i.split('_lowerCI')[0]
                        up = group + '_upperCI'
                        low_ci.append(float(data[name][i]))
                        ci.append(float(data[name][up]) - float(data[name][i]))
                max_mean = max(mean)
                max_ci = max(ci)
                min_ci = min(low_ci)
                if max_ci <= max_mean:
                    n = 1
                else:
                    n = round(max_ci / max_mean)
                l = round(abs(min_ci / n) + max_mean + 3)
                data[name].update({'n': n, 'l': l})
        return data, sort_list

    def check(self):
        """
        检查文件格式是否正确
        """
        pass
