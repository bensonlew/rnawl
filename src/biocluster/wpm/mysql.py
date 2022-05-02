# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

import MySQLdb
from ..config import Config
from DBUtils.PooledDB import PooledDB
from MySQLdb.cursors import DictCursor


class Mysql(object):
    """
        MYSQL数据库对象，负责产生数据库连接 , 此类中的连接采用连接池实现
        获取连接对象：conn = Mysql.getConn()
        释放连接对象;conn.close()或del conn
    """
    # 连接池对象
    pool = None

    def __init__(self):
        """
        数据库构造函数，从连接池中取出连接，并生成操作游标
        """
        self._config = Config()
        self._conn = self.__get_conn()
        self._cursor = self._conn.cursor()

    def __del__(self):
        self.close()

    def close(self):
        try:
            if self._cursor:
                self._cursor.close()
            if self._conn:
                self._conn.close()
            self._cursor = None
            self._conn = None
        except:
            pass

    @property
    def conn(self):
        if not self._conn:
            self._conn = self.__get_conn()
        return self._conn

    @property
    def cursor(self):
        if not self._cursor:
            self._cursor = self.conn.cursor()
        return self._cursor

    def __get_conn(self):
        """
        @summary: 静态方法，从连接池中取出连接
        @return MySQLdb.connection
        """
        if Mysql.pool is None:
            Mysql.pool = PooledDB(creator=MySQLdb, mincached=1, maxcached=5, maxconnections=500,
                                  blocking=False, host=self._config.DB_HOST, user=self._config.DB_USER, maxusage=10000,
                                  passwd=self._config.DB_PASSWD, db=self._config.DB_NAME,
                                  port=int(self._config.DB_PORT), charset='utf8', cursorclass=DictCursor,
                                  setsession=['SET AUTOCOMMIT = 1'])
        return Mysql.pool.connection()

    def get_all(self, sql, param=None):
        """
        @summary: 执行查询，并取出所有结果集
        @param sql:查询ＳＱＬ，如果有查询条件，请只指定条件列表，并将条件值使用参数[param]传递进来
        @param param: 可选参数，条件列表值（元组/列表）
        @return: result list/boolean 查询到的结果集
        """
        if param is None:
            count = self.cursor.execute(sql)
        else:
            count = self.cursor.execute(sql, param)
        if count > 0:
            result = self.cursor.fetchall()
        else:
            result = False
        # print self.cursor._last_executed
        return result

    def get_one(self, sql, param=None):
        """
        @summary: 执行查询，并取出第一条
        @param sql:查询ＳＱＬ，如果有查询条件，请只指定条件列表，并将条件值使用参数[param]传递进来
        @param param: 可选参数，条件列表值（元组/列表）
        @return: result list/boolean 查询到的结果集
        """
        if param is None:
            count = self.cursor.execute(sql)
        else:
            count = self.cursor.execute(sql, param)
        if count > 0:
            result = self.cursor.fetchone()
        else:
            result = False
        # self._conn.commit()
        # print self.cursor._last_executed
        return result

    def get_many(self, sql, num, param=None):
        """
        @summary: 执行查询，并取出num条结果
        @param sql:查询ＳＱＬ，如果有查询条件，请只指定条件列表，并将条件值使用参数[param]传递进来
        @param num:取得的结果条数
        @param param: 可选参数，条件列表值（元组/列表）
        @return: result list/boolean 查询到的结果集
        """
        if param is None:
            count = self.cursor.execute(sql)
        else:
            count = self.cursor.execute(sql, param)
        if count > 0:
            result = self.cursor.fetchmany(num)
        else:
            result = False
        # self._conn.commit()
        # print self.cursor._last_executed
        return result

    def insert_one(self, sql, value):
        """
        @summary: 向数据表插入一条记录
        @param sql:要插入的ＳＱＬ格式
        @param value:要插入的记录数据tuple/list
        @return: insertId 受影响的行数
        """
        self.cursor.execute(sql, value)
        # self._conn.commit()
        # print self.cursor._last_executed
        return self.__get_insert_id()

    def insert_many(self, sql, values):
        """
        @summary: 向数据表插入多条记录
        @param sql:要插入的ＳＱＬ格式
        @param values:要插入的记录数据tuple(tuple)/list[list]
        @return: count 受影响的行数
        """
        count = self.cursor.executemany(sql, values)
        # self._conn.commit()
        # print self.cursor._last_executed
        return count

    def __get_insert_id(self):
        """
        获取当前连接最后一次插入操作生成的id,如果没有则为０
        """
        self.cursor.execute("SELECT @@IDENTITY AS id")
        result = self.cursor.fetchall()
        if result:
            return result[0]['id']
        else:
            return None

    def query(self, sql, param=None):
        if param is None:
            count = self.cursor.execute(sql)
        else:
            count = self.cursor.execute(sql, param)
        return count

    def update(self, sql, param=None):
        """
        @summary: 更新数据表记录
        @param sql: ＳＱＬ格式及条件，使用(%s,%s)
        @param param: 要更新的  值 tuple/list
        @return: count 受影响的行数
        """
        num = self.query(sql, param)
        # self._conn.commit()
        # print self.cursor._last_executed
        return num

    def delete(self, sql, param=None):
        """
        @summary: 删除数据表记录
        @param sql: ＳＱＬ格式及条件，使用(%s,%s)
        @param param: 要删除的条件 值 tuple/list
        @return: count 受影响的行数
        """
        return self.query(sql, param)

    def end(self, option='commit'):
        """
        @summary: 结束事务
        """
        if option == 'commit':
            self.conn.commit()
        else:
            self.conn.rollback()
