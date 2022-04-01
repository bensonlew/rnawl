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


class Otu(object):
    def __init__(self, otu_content, samples):
        self.__level_dict = {1: 'd__', 2: 'k__', 3: 'p__', 4: 'c__', 5: 'o__', 6: 'f__', 7: 'g__', 8: 's__', 9: 'asv'}
        if isinstance(otu_content, types.DictType):
            self.__dict = otu_content
        else:
            try:
                self.__dict = json.loads(otu_content)
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


class otu_table(object):
    def __init__(self, table_content,samples_in_group=None):
        if isinstance(table_content, types.DictType):
            self.__dict = table_content
        else:
            try:
                self.__dict = json.loads(table_content)
            except ValueError, e:
                raise Exception('初始化参数不正确！info:%s' % e)
        self.__mongodb = Config().get_mongo_client(mtype="metaasv")
        self.db_name = Config().get_mongo_dbname("metaasv")
        if samples_in_group == None:
            self._samples = self._get_samples_name()
        else:
            self._samples = self._get_samples_name(samples_in_group)
        self.otus = self._get_all_otus()

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
                        raise Exception('样本collection中不存在specimen字段：%s' % specimen)
                else:
                    raise Exception('asv_specimen中的specimen_id无法在specimen_detail中找到数据')
            return samples
        else:
            raise Exception('ASV表没有样本信息')

    def _get_all_otus(self):
        collection = self.__mongodb[self.db_name]['asv_detail']
        results = collection.find({'asv_id': self._id})
        otus = []
        if results.count():
            for i in results:
                otus.append(Otu(i, self.samples))
            return otus
        else:
            raise Exception('没有找到ASV表对应的detail信息')

    def get_level_otu(self, level=8):
        u"""在某一个水平上的代表OTU，除去其他OTU."""
        if isinstance(level, types.IntType):
            if 9 > level > 0:
                level_otu = defaultdict(lambda: ('', 0))  # 生成字典，值为一个字符串（level名）和一个数字（数量）的元组
                remove_otu = []
                for otu in self.otus:
                    level_name = otu.get_level_name(level)
                    if level_otu[level_name][1] <= otu.total_reads:
                        if level_otu[level_name][0]:
                            remove_otu.append(level_otu[level_name][0])
                        level_otu[level_name] = (otu.asv, otu.total_reads)
                    else:
                        remove_otu.append(otu.asv)
                level_otu = {i[1][0]: i[0] for i in level_otu.iteritems()}
                remove_otu = list(set(remove_otu))
                return level_otu, remove_otu
            else:
                raise Exception('错误的分类水平大小(1-8):%s' % level)
        else:
            raise Exception('错误的分类水平类型,必须是数字(1-8)')

    def level_rank(self, level):
        rank_level = defaultdict(int)
        for otu in self.otus:
            level_name = otu.get_level_name(level)
            rank_level[level_name] += otu.total_reads
        rank_species = sorted(rank_level, key=lambda x: rank_level[x], reverse=True)
        return rank_species


def get_origin_otu(otu_id, connecter=None, database=None, collection='asv'):
    u"""从一个OTU ID找到这个OTU的原始OTU ID."""
    if not connecter:
        connecter = Config().get_mongo_client(mtype="metaasv")
    if not database:
        database = Config().get_mongo_dbname("metaasv")
    collect = connecter[database][collection]
    if isinstance(otu_id, types.StringTypes):
        try:
            ObjectId(otu_id)
        except bosn_InvalidID as e:
            return False, e
    elif isinstance(otu_id, ObjectId):
        otu_id = str(otu_id)
    else:
        return False, '输入otu_id参数必须为可ObjectId字符串或者ObjectId类型'
    origin_id = None
    print("+++++++++__________+++++++++:::{}".format(otu_id))
    while True:
        result = collect.find_one({'_id': ObjectId(otu_id)})
        if not result:
            return False, 'otu_id无法找到对应的数据表'
        else:
            ObjectId(result['asv_id'])
            if result['asv_id'] == otu_id:  # otu_id初始化ID为自己
                origin_id = otu_id
                if isinstance(origin_id, ObjectId):
                    pass
                else:
                    origin_id = ObjectId(origin_id)
                break
            else:
                otu_id = result['asv_id']
    return True, origin_id


def get_otu_phylo_newick(otu_id, connecter=None, database=None,
                         collection='phylo_tree', table_type='otu', tree_type='phylo'):  #tree_type: phylo  表NJ。phylo_ml 表ML， phylo_mp 表 MP
    """
    根据原始表的OTU ID找到其对应的进化树文件
    返回一个两个元素的元组,第一个bool代表找到或者没找到，第二个是错误信息或者结果进化树字符串
    """
    if not connecter:
        connecter = Config().get_mongo_client(mtype="metaasv")
    if not database:
        database = Config().get_mongo_dbname("metaasv")
    if isinstance(otu_id, types.StringTypes):
        try:
            otu_id = ObjectId(otu_id)
        except bosn_InvalidID as e:
            return False, e
    elif not isinstance(otu_id, ObjectId):
        return False, '输入id参数必须为字符串或者ObjectId类型'
    collect = connecter[database][collection]
    result = collect.find_one({'asv_id': otu_id, "name": "Tree_Origin"})
    if not result:
        return False, '没有找到id对应的newick数据'
    else:
        #if result['table_type'] == table_type and result['tree_type'] == tree_type:
        return True, result['newicktree'].rstrip()


#def get_level_newicktree(otu_id, level=9, tempdir='./', return_file=False, topN=None, bind_obj=None):
def get_level_newicktree(asv_id, group=None, level=9, tempdir='./', return_file=False, topN=None, bind_obj=None,tree_type='phylo'):  #增加group参数
    mongo_client = Config().get_mongo_client(mtype="metaasv")
    db_name = Config().get_mongo_dbname("metaasv")
    collection = mongo_client[db_name]['asv']
    tempdir = tempdir.rstrip('/') + '/'
    temptre = tempdir + 'temp_newick.tre'
    filter_tre = tempdir + 'temp_filter_newick.tre'
    origin_id = get_origin_otu(asv_id)
    if bind_obj:
        bind_obj.logger.info('origin_id:' + str(origin_id))
        bind_obj.logger.info("+++++++++__________+++++++++:::{}\n{}".format(origin_id, asv_id))
    if origin_id[0]:
        level_newick = get_otu_phylo_newick(origin_id[1])  # get_otu_phylo_newick返回一个二元 元组，第一个代表是否找到  ##add tree_type
        if level_newick[0]:
            pass
        else:
            raise Exception('原始asv表找不到对应进化树文件')
    else:
        raise Exception('asv表没有找到对应的原始表')
    if bind_obj:
        bind_obj.logger.info('origin_newick:' + str(level_newick)[:200])
        bind_obj.logger.info("+++++++++__________+++++++++:::{}".format(origin_id))
    if isinstance(asv_id, ObjectId):
        filter_find = collection.find_one({'_id': asv_id})
    else:
        filter_find = collection.find_one({'_id': ObjectId(asv_id)})
    if group == None :
        filter_otu_table = otu_table(filter_find)  # 生成查询OTU表对象
    else:
        filter_otu_table = otu_table(filter_find,group) #生成查询OTU表对象，只包含分组中的样品
    result = collection.find_one({'_id': origin_id[1]})
    this_otu_table = otu_table(result)  # 生成原始表的OTU表对象

    '''#
    此处将filter_otu_table按分组进行挑选
    '''#
#    if group ！= 'All':
 #       with open(group) as g:
  #          g.readline()
   #         sample_group = {}


    remove_otus = _get_remove_otus(filter_otu_table, this_otu_table)
    wbf = open(temptre, 'wb')
    wbf.write(level_newick[1])
    wbf.close()
    phylo_newick = Phylo.read(temptre, 'newick')  # 生成原始进化树对象，没有找到直接从内存字符串解析的方式
    for i in remove_otus:  # 移除查询表和原始表的差异otu
        phylo_newick.prune(i)
    if isinstance(level, types.IntType):
        if level == 9:
            pass
        elif 9 > level > 0:
            try:
                level_otu = filter_otu_table.get_level_otu(level)  # 返回元组，0为otu对应level名的字典，1为需要移除的otu
                for remove_one in level_otu[1]:  # 移除选定特定level而移除的otu
                    phylo_newick.prune(remove_one)
                terminals = phylo_newick.get_terminals()
                for i in terminals:  # 将枝的名字改为该水平名称
                    i.name = level_otu[0][i.name]
            except IOError as e:
                raise Exception('进化树依据水平过滤出错,info:%s' % e)
        else:
            raise Exception('错误的分类水平大小(1-9):%s' % level)
        if topN:
            rank_species = filter_otu_table.level_rank(level)
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


def _get_remove_otus(small_otus, total_otus):
    """
    从两个OTUtable对象中找到，小的哪些OTUs不存在，用于从原始进化树种删除该部分不存在的OTU
    """
    small = set([otu.asv for otu in small_otus.otus])
    total = set([otu.asv for otu in total_otus.otus])
    if small & total != small:
        raise Exception('子ASV表并不完全在父ASV表中')
    remove_otus = list(total - small)
    return remove_otus