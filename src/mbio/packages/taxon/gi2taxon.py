# -*- coding: utf-8 -*-
# __author__ = 'sheng.he'
# last_modified:20160621
# modified by konghualei 20161116

import subprocess
from biocluster.config import Config


db = Config().get_mongo_client(mtype="metagenomic", ref=True)[Config().get_mongo_dbname("metagenomic", ref=True)]
# db = Config().biodb_mongo_client.sanger_biodb
nr_collection = db.NR_sequence
taxon_collection = db.species_taxon


def create_gi2taxid(name, db, dbfile):
    """用于在sqlite3的数据库中创建表，导入gi2taxid的数据"""
    sql = 'CREATE TABLE {} (gi INT PRIMARY KEY NOT NULL, taxid INT NOT NULL);'.format(
        name)
    pop = subprocess.Popen('sqlite3 {} -cmd \"{}\"'.format(db, sql),
                           stdin=subprocess.PIPE,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    # print '.separator "\\t" "\\n"\n.q'
    out, err = pop.communicate('.q')
    pop.wait()
    if 'Error' in out or err:
        raise Exception(err)
    pop = subprocess.Popen('sqlite3 {}'.format(db),
                           stdin=subprocess.PIPE,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    out, err = pop.communicate(
        '.separator \\t\n.import {} {}\n.q\n'.format(dbfile, name))
    pop.wait()
    if 'Error' in out or err:
        raise Exception(err)


# 默认的数据库位置，待定
# gi_db = sqlite3.connect(Config().SOFTWARE_DIR + '/database/ncbi_taxon/2taxid_sqlite3.db')
# taxon_db = sqlite3.connect(Config().SOFTWARE_DIR + '/database/ncbi_taxon/taxa.sqlite')


class taxon(object):
    """ncbi物种分类中一个节点分类的对象，提供一些查询方法和信息存储"""

    def __init__(self):
        self.taxid = None  # 分类id
        self.gi_prot = None  # 对应的gi_prot_id
        self.gi_nucl = None  # 对应的gi_nucl_id
        self.track = []  # id父分类链
        self.rank = None  # 等级
        self.track_taxon = None  # 父分类(spname,rank)链
        self.parent = None  # 父id
        self.spname = None  # 分类名称

    def get_taxon(self):
        if self.taxid:
            # result = taxon_db.execute('select spname, rank from species where taxid={}'.format(self.taxid))
            result = taxon_collection.find({"_id": self.taxid})
            count = 0
            for one in result:
                # self.spname = one[0]
                self.spname = one['species_name']
                # self.rank = one[1]
                self.rank = one['rank']
                count += 1
            if count == 1:
                return (self.spname, self.rank)
            elif count == 0:
                print 'WARNNING:taxid:{}在数据库中没有信息'.format(self.taxid)
            else:
                raise Exception('taxid:{}在数据库中存在多条记录'.format(self.taxid))

    def get_taxid(self, gi, table='prot'):
        gi = int(gi)
        if table == 'prot':
            self.gi_prot = gi
            table = 'prot_gi2taxid'
        elif table == 'nucl':
            self.gi_nucl = gi
            table = 'nucl_gi2taxid'
        else:
            raise TypeError('错误的table格式，必须为prot或者nucl')
        # result = gi_db.execute('select taxid from {} where gi={}'.format(table, gi))
        result = nr_collection.find({"_id": gi})
        count = 0
        for one in result:
            _taxid = str(one["taxid"])
            if _taxid != "None":
                self.taxid = int(one["taxid"])
                count += 1
            if count == 1:
                return self
            elif count == 0:
                print 'WARNNING:gi:{}在{}没有物种分类注释信息'.format(gi, table)
            else:
                raise Exception(
                    '数据库中存在相同gi的信息table:{},gi:{}'.format(table, gi))

    def get_track(self):
        if self.taxid:
            # result = taxon_db.execute('select track from species where taxid={}'.format(self.taxid))
            result = taxon_collection.find({"_id": self.taxid})
            count = 0
            for one in result:
                count += 1
                # for i in one[0].split(','):
                for i in one['tree_path'].split(','):
                    track_one = taxon()
                    track_one.taxid = int(i)
                    track_one.get_taxon()
                    self.track.append(track_one)
            if count == 1:
                return self.track
            elif count == 0:
                print 'WARNNING:taxid:{}在表中不存在'
            else:
                raise Exception('数据库中存在多条taxid:{}的信息'.format(self.taxid))

    def get_track_id(self):
        if self.taxid:
            # result = taxon_db.execute('select track from species where taxid={}'.format(self.taxid))
            result = taxon_collection.find({"_id": self.taxid})
            count = 0
            for one in result:
                # self.trackid = one[0]
                self.track_taxon = one["tree_path"]
                count += 1
            if count == 1:
                # return self.trackid
                return self.track_taxon
            elif count == 0:
                print 'WARNNING:taxid:{}在数据库中没有信息'.format(self.taxid)
            else:
                raise Exception('taxid:{}在数据库中存在多条记录'.format(self.taxid))


def gi_taxon(gilist, db):
    if not isinstance(gilist, list) and not isinstance(gilist, set):
        raise TypeError('传入参数必须是gi的列表')
    taxons = [taxon().get_taxid(i, db) for i in gilist]
    gi_taxons = dict(zip(gilist, taxons))
    for item in gi_taxons.iteritems():
        if item[1]:
            item[1].get_track()
            item[1].track.reverse()
            track = ['{}{{{}}}'.format(i.spname, i.rank)
                     for i in item[1].track]
            track = ';'.join(track)
            gi_taxons[item[0]] = track
    return gi_taxons

if __name__ == '__main__':  # for test
    a = taxon()
    a.get_taxid('148678297')
    a.get_track()
    print([(i.spname, i.rank) for i in a.track])
