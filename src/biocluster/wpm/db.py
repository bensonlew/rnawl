# -*- coding: utf-8 -*-
# __author__ = 'HD'
import json
import datetime
import sys
import psycopg2
from ..config import Config
import traceback
from DBUtils.PooledDB import PooledDB


class Postgresql(object):
    """
    进行postgresql操作, 该脚本不用了
    """
    def __init__(self):
        self._pool = None
        self._config = Config()
        self._conn = self.__get_conn()
        self._cursor = self._conn.cursor()

    def __del__(self):
        self.close()

    def conn(self):
        if not self._conn:
            self._conn = self.__get_conn()
        return self._conn

    def __get_conn(self):
        if self._pool is None:
            self.init_pgsql_pool()
        return self._pool.connection()

    def cursor(self):
        if not self._cursor:
            self._cursor = self.conn.cursor()
        return self._cursor

    def close(self):
        try:
            if self._cursor:
                self._cursor.close()
            if self._conn:
                self._conn.close()
            if self._pool:
                self._pool.close()
            self._cursor = None
            self._conn = None
        except:
            pass

    def init_pgsql_pool(self):
        try:
            pool = PooledDB(
                creator=psycopg2,  # 使用连接数据库的模块 psycopg2
                maxconnections=6,  # 连接池允许的最大连接数，0 和 None 表示不限制连接数
                mincached=1,  # 初始化时，链接池中至少创建的空闲的链接，0 表示不创建
                maxcached=4,  # 链接池中最多闲置的链接，0 和 None 不限制
                blocking=True,  # 连接池中如果没有可用连接后，是否阻塞等待。True，等待；False，不等待然后报错
                maxusage=None,  # 一个链接最多被重复使用的次数，None 表示无限制
                setsession=[],  # 开始会话前执行的命令列表
                host=self._config.DB_HOST,
                port=int(self._config.DB_PORT),
                user=self._config.DB_USER,
                password=self._config.DB_PASSWD,
                database=self._config.DB_NAME)
            self._pool = pool
        except Exception as e:
            print 'ERROR: create postgresql pool failed on {}:{}.\n'.format(datetime.datetime.now(), e)
            exstr = traceback.format_exc()
            print exstr
            print e
            sys.stdout.flush()

    def findone(self, sql):
        """
        查询一条记录
        :return:
        """
        count = self.cursor.execute(sql)
        if count > 0:
            return self.cursor.fetchone()
        else:
            return False

    def insert(self, sql, value):
        self.cursor.execute(sql, value)
        pass
