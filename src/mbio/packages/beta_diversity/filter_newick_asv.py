# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
# __version__ = 'v1.0'
# __last_modified__ = '20171101'
"""

"""

from Bio import Phylo
from bson.objectid import ObjectId
from bson.errors import InvalidId as bosn_InvalidID
from collections import defaultdict
from biocluster.config import Config
import json
import types
import re

# db_name = Config().MONGODB
# mongo_client = Config().mongo_client
def get_mongo():
    mongo_client = Config().get_mongo_client(mtype="metaasv")
    db_name = Config().get_mongo_dbname("metaasv")
    return mongo_client, db_name


class Asv(object):
    def __init__(self, asv_content, samples):
        self.__level_dict = {1: 'd__', 2: 'k__', 3: 'p__', 4: 'c__', 5: 'o__', 6: 'f__', 7: 'g__', 8: 's__', 9: 'asv'}
        if isinstance(asv_content, types.DictType):
            self.__dict = asv_content
        else:
            try:
                self.__dict = json.loads(asv_content)
            except ValueError as e:
                raise Exception('初始化参数不正确！info:%s' % e)
        self.total_reads = self.count_samples(samples)

    def __getattr__(self, name):
        if name in self.__dict:
            return self.__dict[name]
        else:
            raise Exception('错误的属性名称：可能是asv表detail中缺失数据！')

    @property
    def dict_value(self):
        return self.__dict

    def get_level_name(self, level=9):
        if isinstance(level, types.IntType):
            pass
        else:
            raise Exception('错误的分类水平类型,必须是数字(1-9)')
        if level == 9:
            return getattr(self, self.num_to_level(level))
        name_list = []
        i = 1
        while i <= level:
            name_list.append(getattr(self, self.num_to_level(i)))
            i += 1
        name = '; '.join(name_list)
        return name

    def count_samples(self, samples):
        """
        pass
        """
        total_value = 0
        for i in samples:
            total_value = total_value + self.dict_value[i]
        return total_value

    def num_to_level(self, num):
        if isinstance(num, types.IntType):
            if 10 > num > 0:
                return self.__level_dict[num]
            else:
                raise Exception('错误的分类水平大小(1-9):%s' % num)
        else:
            raise Exception('错误的分类水平类型,必须是数字(1-9)')


class asv_table(object):
    def __init__(self, table_content,samples_in_group=None):
        if isinstance(table_content, types.DictType):
            self.__dict = table_content
        else:
            try:
                self.__dict = json.loads(table_content)
            except ValueError, e:
                raise Exception('初始化参数不正确！info:%s' % e)
        self.__mongodb, self.db_name = get_mongo()
        if samples_in_group == None:
            self._samples = self._get_samples_name()
        else:
            self._samples = self._get_samples_name(samples_in_group)
        self.asvs = self._get_all_asvs()

    def __getattr__(self, name):
        if name in self.__dict:
            return self.__dict[name]
        else:
            raise Exception('错误的属性名称！')

    @property
    def dict_value(self):
        return self.__dict

    @property
    def samples(self):
        return self._samples

    def _get_samples_name(self,groups=None):      #增加groups参数，根据分组中的样品名称挑选样品
        collection = self.__mongodb[self.db_name]['asv_specimen']
        specimen_collection = self.__mongodb[self.db_name]['specimen_detail']
        results = collection.find({'asv_id': self._id})
        samples = []
        if results.count():
            for i in results:
                specimen = specimen_collection.find_one({'_id': i['specimen_id']})
                if specimen:
                    try:
                        if groups == None:        #如果不提供分组，则挑选全部样品
                            samples.append(specimen['specimen'])
                        elif specimen['specimen'] in groups:        #如果此样品在选择的分组中，则挑选此样品
                            samples.append(specimen['specimen'])
                    except KeyError:
                        raise Exception('样本collection中不存在specimen_name字段：%s' % specimen)
                else:
                    raise Exception('asv_specimen中的specimen_id无法在specimen_detail中找到数据')
            return samples
        else:
            raise Exception('asv表没有样本信息')

    def _get_all_asvs(self):
        collection = self.__mongodb[self.db_name]['asv_detail']
        results = collection.find({'asv_id': self._id})
        asvs = []
        if results.count():
            for i in results:
                asvs.append(Asv(i, self.samples))
            return asvs
        else:
            raise Exception('没有找到asv表对应的detail信息')

    def get_level_asv(self, level=8):
        u"""在某一个水平上的代表asv，除去其他asv."""
        if isinstance(level, types.IntType):
            if 9 > level > 0:
                level_asv = defaultdict(lambda: ('', 0))  # 生成字典，值为一个字符串（level名）和一个数字（数量）的元组
                remove_asv = []
                for asv in self.asvs:
                    level_name = asv.get_level_name(level)
                    if level_asv[level_name][1] <= asv.total_reads:
                        if level_asv[level_name][0]:
                            remove_asv.append(level_asv[level_name][0])
                        level_asv[level_name] = (asv.asv, asv.total_reads)
                    else:
                        remove_asv.append(asv.asv)
                level_asv = {i[1][0]: i[0] for i in level_asv.iteritems()}
                remove_asv = list(set(remove_asv))
                return level_asv, remove_asv
            else:
                raise Exception('错误的分类水平大小(1-8):%s' % level)
        else:
            raise Exception('错误的分类水平类型,必须是数字(1-8)')

    def level_rank(self, level):
        rank_level = defaultdict(int)
        for asv in self.asvs:
            level_name = asv.get_level_name(level)
            rank_level[level_name] += asv.total_reads
        rank_species = sorted(rank_level, key=lambda x: rank_level[x], reverse=True)
        return rank_species


def get_origin_asv(asv_id, connecter=None, database=None, collection='asv',bind_obj=None):
    u"""从一个asv ID找到这个asv的原始asv ID."""
    mongo_client, db_name = get_mongo()
    bind_obj.logger.info('ok1')
    if not connecter:
        connecter = mongo_client
    if not database:
        database = db_name
    collect = connecter[database][collection]
    if isinstance(asv_id, types.StringTypes):
        try:
            ObjectId(asv_id)
        except bosn_InvalidID as e:
            return False, e
    elif isinstance(asv_id, ObjectId):
        asv_id = str(asv_id)
    else:
        return False, '输入asv_id参数必须为可ObjectId字符串或者ObjectId类型'
    origin_id = None
    while True:
        result = collect.find_one({'_id': ObjectId(asv_id)})
        if not result:
            return False, 'asv_id无法找到对应的数据表'
        else:
            ObjectId(result['asv_id'])
            if result['asv_id'] == ObjectId(asv_id):  # asv_id初始化ID为自己
                origin_id = asv_id
                if isinstance(origin_id, ObjectId):
                    pass
                else:
                    origin_id = ObjectId(origin_id)
                break
            else:
                origin_id = result['asv_id']
                break
    return True, origin_id


def get_asv_phylo_newick(asv_id, connecter=None, database=None,
                         collection='phylo_tree'):  #tree_type: phylo  表NJ。phylo_ml 表ML， phylo_mp 表 MP
    """
    根据原始表的asv ID找到其对应的进化树文件
    返回一个两个元素的元组,第一个bool代表找到或者没找到，第二个是错误信息或者结果进化树字符串
    """
    mongo_client, db_name = get_mongo()
    if not connecter:
        connecter = mongo_client
    if not database:
        database = db_name
    if isinstance(asv_id, types.StringTypes):
        try:
            asv_id = ObjectId(asv_id)
        except bosn_InvalidID as e:
            return False, e
    elif not isinstance(asv_id, ObjectId):
        return False, '输入id参数必须为字符串或者ObjectId类型'
    collect = connecter[database][collection]
    #result = collect.find_one({'table_id': asv_id, "table_type": "asv", "tree_type": "phylo"})
    result = collect.find_one({'asv_id': asv_id, "name": "Tree_Origin"})
    if not result:
        return False, '没有找到id对应的newick数据'
    else:
        #if result['table_type'] == table_type and result['tree_type'] == tree_type:
        return True, result['newicktree'].rstrip()
        #else:
        #    return False, '找到id对应数据，但是table或者tree类型不正确'


#def get_level_newicktree(asv_id, level=9, tempdir='./', return_file=False, topN=None, bind_obj=None):
# asv使用
def get_level_newicktree(asv_id, group=None, level=9, tempdir='./', return_file=False, topN=None, bind_obj=None,tree_type='phylo'):  #增加group参数
    mongo_client, db_name = get_mongo()
    collection = mongo_client[db_name]['asv']
    tempdir = tempdir.rstrip('/') + '/'
    temptre = tempdir + 'temp_newick.tre'
    filter_tre = tempdir + 'temp_filter_newick.tre'
    origin_id = get_origin_asv(asv_id,bind_obj=bind_obj)
    if bind_obj:
        bind_obj.logger.info('origin_id:' + str(origin_id))
    if origin_id[0]:
        level_newick = get_asv_phylo_newick(origin_id[1])  # get_asv_phylo_newick返回一个二元 元组，第一个代表是否找到  ##add tree_type
        if level_newick[0]:
            pass
        else:
            raise Exception('原始asv表找不到对应进化树文件')
    else:
        raise Exception('asv表没有找到对应的原始表')
    if bind_obj:
        bind_obj.logger.info('origin_newick:' + str(level_newick)[:200])
    if isinstance(asv_id, ObjectId):
        filter_find = collection.find_one({'_id': asv_id})
    else:
        filter_find = collection.find_one({'_id': ObjectId(asv_id)})
    if group == None :
        filter_asv_table = asv_table(filter_find)  # 生成查询asv表对象
    else:
        filter_asv_table = asv_table(filter_find,group) #生成查询asv表对象，只包含分组中的样品
    result = collection.find_one({'_id': origin_id[1]})
    this_asv_table = asv_table(result)  # 生成原始表的asv表对象

    '''#
    此处将filter_asv_table按分组进行挑选
    '''#
#    if group ！= 'All':
 #       with open(group) as g:
  #          g.readline()
   #         sample_group = {}


    remove_asvs = _get_remove_asvs(filter_asv_table, this_asv_table)
    wbf = open(temptre, 'wb')
    wbf.write(level_newick[1])
    wbf.close()
    phylo_newick = Phylo.read(temptre, 'newick')  # 生成原始进化树对象，没有找到直接从内存字符串解析的方式
    for i in remove_asvs:  # 移除查询表和原始表的差异asv
        phylo_newick.prune(i)
    if isinstance(level, types.IntType):
        if level == 9:
            pass
        elif 9 > level > 0:
            try:
                level_asv = filter_asv_table.get_level_asv(level)  # 返回元组，0为asv对应level名的字典，1为需要移除的asv
                for remove_one in level_asv[1]:  # 移除选定特定level而移除的asv
                    phylo_newick.prune(remove_one)
                terminals = phylo_newick.get_terminals()
                for i in terminals:  # 将枝的名字改为该水平名称
                    i.name = level_asv[0][i.name]
            except IOError as e:
                raise Exception('进化树依据水平过滤出错,info:%s' % e)
        else:
            raise Exception('错误的分类水平大小(1-9):%s' % level)
        if topN:
            rank_species = filter_asv_table.level_rank(level)
            if len(rank_species) > topN:
                for one in rank_species[topN:]:
                    phylo_newick.prune(one)

        Phylo.write(phylo_newick, filter_tre, 'newick')
        
        ####zouguanqing.zou
        fr = open(filter_tre)
        line = fr.readline()
        def format_change(matched):
            instr = matched.group("group")
            return instr[1:]

        replacedStr1 = re.sub("(?P<group>:[\d\.]+:[\d\.]+)", format_change, line)
        fr.close()

        def simple_name(name):
            name = name.group()
            return name.split(';')[-1].strip().replace(':', '-').strip('\'')  # replace用于去掉名称中带有冒号

        replacedStr = re.sub(r'\'(.+?)\'', simple_name, replacedStr1)

        fw = open(filter_tre+'.bootstrap','w')
        fw.write(replacedStr)
        fw.close()
        
        ####

        phylo_newick = Phylo.read(filter_tre, 'newick')  # 一次写入和一次读取可以对进化树的格式简略化，主要是由于发现经过prune修剪的的树会出现空枝
        format_phylo_newick = phylo_newick.format('newick')  # 格式化返回字符串newick树
        if return_file:  # 如果需要返回文件则，生成文件在提供的文件夹中，并返回文件路径
            tempfile = open(filter_tre, 'w')
            tempfile.write(format_phylo_newick)
            return filter_tre
        else:
            return format_phylo_newick
    else:
        raise Exception('错误的分类水平类型,必须是数字(1-9)')


def _get_remove_asvs(small_asvs, total_asvs):
    """
    从两个asvtable对象中找到，小的哪些asvs不存在，用于从原始进化树种删除该部分不存在的asv
    """
    small = set([asv.asv for asv in small_asvs.asvs])
    total = set([asv.asv for asv in total_asvs.asvs])
    if small & total != small:
        raise Exception('子asv表并不完全在父asv表中')
    remove_asvs = list(total - small)
    return remove_asvs


# asv_id = '589817e7a4e1af69eccae310'
# get_level_newicktree(asv_id=asv_id, level=7, tempdir='./', return_file=True, topN=10)
