"""
table io for the random forest tool
version 0.0.1
author: yingnn
date: 2017.10.16

"""
from __future__ import print_function
import os
import datetime
from bson import SON
import pandas as pd
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
import rf_const as CONSTS
from rf_const import read, get_header


class Rf(Base):
    """
    
    """
    def __init__(self, bind_object):
        super(Rf, self).__init__(bind_object)
        self.output_dir = self.bind_object.output_dir
        self.work_dir = self.bind_object.work_dir
        self._project_type = 'toolapps'
        self._db_name = 'ttoolapps'
        self.tables = CONSTS.tables
        self.bind_object.logger.info('%r' % self.db['table'])
        self.bind_object.logger.info('%s' % os.getcwd())
        os.chdir(self.output_dir)  # !
        self.bind_object.logger.info('%s' % os.getcwd())

    @report_check
    def run(self):
        # table_ids = []
        table_id_ = None
        for table in self.tables:
            print(table)
            self.bind_object.logger.info('%r' % table)
            if os.path.exists(table['data_path']):
                if table['db'] == CONSTS.db_bar:
                    table_id_ = self.save2table(fp=table['data_path'],
                                                table_name=table['db'],
                                                table_columns=table['columns'],
                                                db=table['db'],
                                                attr_col=0,
                                                fmt_col=1,
                                                info=table, )
                elif table['db'] == CONSTS.db_curve and table_id_ is not None:
                    self.save2table(fp=table['data_path'],
                                    table_name=table['db'],
                                    table_columns=table['columns'],
                                    db=table['db'],
                                    id_table=table_id_,
                                    attr_col=0,
                                    # fmt_col=1,
                                    columns_table=['feature_number', 'error', ],
                                    info=table, )
                else:
                    self.save2table(fp=table['data_path'],
                                    table_name=table['name'],
                                    table_columns=table['columns'],
                                    db=table['db'], )
                # table_ids.append(table_id)
        # return table_ids

    def save2table(self, fp, table_name, table_desc='', table_columns=None,
                   db=CONSTS.db_public, id_table=None,
                   columns_table=None,
                   category_label=None, attr_col=None, fmt_col=None, info=None):
        if table_columns is None:
            table_columns = get_header(fp)
        if columns_table is None:
            columns_table = table_columns
        if category_label is None:
            category = [table_columns[1]]
        if attr_col is None:
            attr = table_columns
        else:
            df = pd.read_table(fp, index_col=attr_col)
            attr = df.index.tolist()

        son = SON(
            project_sn=self.bind_object.sheet.project_sn,
            task_id=self.bind_object.id,
            attrs=attr,
            name=table_name,
            desc=table_desc,
            categories=category,
            status=CONSTS.status_start,
            created_ts=datetime.datetime.now().strftime(
                "%Y-%m-%d %H:%M:%S"),
            table_name=table_name,
            table_columns=table_columns,
            table_desc=table_desc,
            table_status=CONSTS.status_start,
            table_datetime=datetime.datetime.now().strftime(
                "%Y-%m-%d %H:%M:%S"), )
        if info is not None:
            for key in ['title', 'x_label', 'y_label', 'show_legend']:
                if info.has_key(key):
                    if info[key] is None:
                        # table best_feature_numbers
                        son[key] = table_columns[1]
                    else:
                        son[key] = info[key]
        if id_table is not None:
            self.bind_object.logger.info('id_table: %r' % id_table)
            son['_id'] = id_table
        table_id = self.db[db].insert_one(son).inserted_id
        self.bind_object.logger.info('table_id: %r' % table_id)
        
        data = read(fp)
        if data is None:  # has no data body, only data header
            self.db[db].update_one({CONSTS.table_id: table_id},
                                   {CONSTS.table_set:
                                   {CONSTS.table_status: CONSTS.status_end,
                                    'status': CONSTS.status_end}})
            return table_id

        data_insert = []
        id_ = {'_'.join([table_name, 'id']): table_id,
               'table_id': table_id}
        for dat in data:
            data_s = SON(id_)
            for i in range(len(columns_table)):
                if fmt_col == i:
                    data_s[columns_table[i]] = [float(dat[i]), ]
                else:
                    data_s[columns_table[i]] = dat[i]
            data_insert.append(data_s)

        self.db[CONSTS.db_name_sep.join(
            [db, CONSTS.details])].insert_many(data_insert)

        self.db[db].update_one({CONSTS.table_id: table_id},
                               {CONSTS.table_set:
                               {CONSTS.table_status: CONSTS.status_end,
                                'status': CONSTS.status_end}})
        return table_id
