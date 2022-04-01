# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
# last_modify: 20190906
import datetime
from types import StringTypes

from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from bson.son import SON
import pandas as pd


class DiffComp(Base):
    def __init__(self, bind_object):
        super(DiffComp, self).__init__(bind_object)
        self._project_type = 'bac_comparative'
        self._project_sn = ''
        self._task_id = self.bind_object._sheet.task_id
        self._name = self.bind_object.name
        self._work_dir = self.bind_object.work_dir
        self._output_dir = self.bind_object.output_dir

        self._main_table_id = None

    @report_check
    def add_main(self, main_table, params=None):
        '''
        :param main_table: 主表名称
        :param spe_list: 排序的物种名，字符串分隔
        :param params: 前段传过来的参数列表
        :param graph_dir: 如果为 KEGG 比较分析，需要传入KEGG地图文件夹的路径
        :return : 返回主表id
        '''
        self._main_table_id = main_table + '_id'
        created_ts = datetime.datetime.now()
        data = {
            'project_sn': self._project_sn,
            'task_id': self._task_id,
            'name': self._task_id + created_ts.strftime("%Y%m%d_%H%M%S"),
            'params': params,
            'desc': self._name + '差异分析',
            'status': 'end',
            'created_ts': created_ts.strftime("%Y-%m-%d %H:%M:%S")
        }
        return self._run_add(main_table, data, main=True)

    @report_check
    def add_detail(self, filepath, table, main_id, main_table,
                   columns=None, mongo_keys=None, header=True, header_lower=False,
                   tag_key=[], tag_value=[]):
        '''
        :param table: 表的名称
        :param main_id: 对应的主表id
        :param filepath: 用于导表的文件路径，其包括表头和数据行，表头为对应mongo
                         表的字段
        :param main_table: mongo表中main_id对应的字段名
        :param columns: list 选择用于导表的列，若为None，即全选 columns=header
        :param mongo_keys: dict 新旧名字的字典{oldname: newname}
        :param header: 导表文件中是否带有表头，即列名，没有则为数字索引，从0起始
                       header=False时，columns为所选的列号
        :param header_lower: 是否将表头改为小写
        :param tag_key: list 详情表中额外添加的字段 tag_key 和 tag_value 等长
        :param tag_key: list 详情表中额外添加的字段对应的值，
        '''
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error("{}必须为ObjectID对象或其字符串形式".format(main_id))
        h = 0 if header else None
        df = pd.read_csv(filepath, header=h, index_col=False, sep='\t')
        if columns:
            df = df[columns]
        if mongo_keys:
            df = df.rename(columns=mongo_keys)
        df[main_table] = main_id
        if tag_key:
            if len(tag_key) != len(tag_value):
                self.set_error('新加的tag key列表必须和value列表等长')
            for i in range(len(tag_key)):
                df[tag_key[i]] = tag_value[i]
        if header_lower:
            df.columns = map(lambda x: x.lower(), df.columns)
        data = self.df_to_mongo(df)
        self._run_add(table, data)

    def df_to_mongo(self, df):
        keys = map(lambda x: x.replace('.', '_'), df.columns)
        mongo_data = []
        df.apply(
            lambda x: mongo_data.append(SON(dict(zip(keys, x)))), axis=1
            )
        return mongo_data

    def update_sg_status(self, table_id, data):
        self.update_table('sg_status', table_id, data, search_id='table_id')

    def update_table(self, table, table_id, data, search_id='_id'):
        tb = self.db[table]
        tb.update({search_id: ObjectId(table_id)},
                  {'$set': data})

    def _run_add(self, name, data, main=False):
        try:
            if main:
                collection = self.db[name]
                main_id = collection.insert_one(SON(data)).insert_id
                collection.update_one({'_id': main_id},
                                      {'$set': {'main_id': main_id}},
                                      )
                return main_id
            else:
                collection = self.db[name]
                collection.insert_many(data)
        except Exception, e:
            self.bind_object.set_error('导入表%s出错：%s') % (name, e)
        else:
            self.bind_object.logger.info('导入表%s成功！' % name)
