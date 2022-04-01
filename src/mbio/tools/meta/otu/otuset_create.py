# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
from bson.objectid import ObjectId
import pandas as pd
import re
import os
import json


class OtusetCreateAgent(Agent):
    def __init__(self, parent):
        super(OtusetCreateAgent, self).__init__(parent)
        options = [
            {"name": "table_name", "type": "string"},
            {"name": "table_id", "type": "string"},
            {"name": "pvalue", "type": "float", "default": 1.0},
            {"name": "qvalue", "type": "float", "default": 1.0},
            {"name": "species_name", "type": "string", "default": ""},
            {"name": "label", "type": "string"},
            {"name": "lda", "type": "float", "default": 0.0},
            {"name": "top", "type": "int", "default": 100},
            {"name": "main_id", "type": "string"},
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数二次检查
        """
        if not self.option('table_name'):
            self.set_error('必须设置参数 table_name')
        if self.option("table_name") == 'sg_otu_venn':
            if not self.option('label'):
                self.set_error('在venn图分析页面必须配置参数 "label"')
        if self.option("table_name") in ["sg_species_difference_lefse"]:
            if (not self.option("pvalue")) and (not self.option("lda")):
                raise OptionError("{},{}不存在请检查".format(self.option("pvalue"), self.option("lda")))
        if self.option("table_name") in ["sg_randomforest"]:
            if not self.option("top"):
                raise OptionError("{}不存在请检查".format(self.option("top")))
        if self.option("table_name") in ["sg_species_difference_check"]:
            if (not self.option("pvalue")) and (not self.option("qvalue")) and (not self.option("species_name")):
                raise OptionError("{},{},{}不存在请检查".format(self.option("pvalue"), self.option("qvalue"), self.option("species_name")))

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'


class OtusetCreateTool(Tool):
    def __init__(self, config):
        super(OtusetCreateTool, self).__init__(config)
        Config().DBVersion = self.config.DBVersion
        self.client = Config().get_mongo_client(mtype='meta')
        self.db = self.client[Config().get_mongo_dbname(mtype='meta')]
        self.level_dict = {
            1: "d__", 2: "k__", 3: "p__", 4: "c__", 5: "o__",
            6: "f__", 7: "g__", 8: "s__", 9: "otu"
        }

    def run(self):
        super(OtusetCreateTool, self).run()
        self.get_main_infos()
        names, field = self.get_filtered_data()
        self.create_sets(names, field)
        self.client.close()
        self.end()

    def get_main_infos(self):
        main_table = self.db[self.option('table_name')]
        self.main_info = main_table.find_one({'_id': ObjectId(self.option('table_id'))})
        self.params = json.loads(self.main_info['params'])

    def create_sets(self, names, field):
        self.logger.info('names {}\nfield: {}'.format(names, field))
        if field == 'otu':
            self._update_otuset(names)
        else:
            otu_main = self.db['sg_otu'].find_one({'task_id': self.main_info['task_id'],
                                                   'name': 'OTU_Taxon_Origin'})
            if not otu_main:
                otu_main = self.db['sg_otu'].find_one({'_id': self.main_info['otu_id']})
            if 'main_id' in otu_main:
                otu_id = otu_main['main_id']
            else:
                otu_id = otu_main['_id']
            otu_table = self.get_table_data(
                'sg_otu_detail',
                {'otu_id': otu_id},
                ['otu', ] + field
            )
            otu_data = []
            for i in otu_table:
                otu_data.append(i)
            otu_table = pd.DataFrame(otu_data)
            print otu_table.head()
            selection = [False for _ in range(len(otu_data))]
            for one in field:
                if self.option('table_name') == 'sg_species_difference_lefse':
                    otu_table[one] = otu_table[one].apply(
                        lambda x: x.replace('.', '_')
                    )
                selection |= otu_table[one].isin(names)

            self._update_otuset(otu_table[selection]['otu'])

    def _update_otuset(self, sets):
        li = list(set(sets))
        li_len = len(li)
        self.db['sg_otuset'].update_one(
            {'_id': ObjectId(self.option('main_id'))},
            {'$set': {'otuset_list': li, 'otuset_size': li_len}}
        )

    def get_filtered_data(self):
        '''
        统一的mongo表筛选入口，
        返回值`names`: 用于在otu表中筛选内容，
        返回值`field`: 用于指定 `names` otu表中所在的项
        `query_hander`中所有键值对应的方法必须返回 name 和 field 两个值
        '''
        query_hander = {
            'sg_otu_venn': self.for_venn,
            'sg_randomforest': self.for_randomforset,
            'sg_species_difference_check': self.for_diff,
            'sg_species_difference_lefse': self.for_lefse,
        }
        names, field = query_hander[self.option('table_name')]()
        return names, field

    def get_table_data(self, table, searcher, include):
        '''
        筛选获得所有满足条件建的数据
        :params talbe: collection名称
        :params searcher: 筛选条件
        :params include: 指定返回的内容包含的字段
        '''
        tb = self.db[table]
        data = []
        for one in tb.find(searcher):
            data.append({k: one[k] for k in include})
        if not data:
            self.set_error('no data after filter the table '
                           '{} by {} and {}'.format(table, searcher, include))
        print ">>top 10 table data for {}\n".format(table)
        print data[:10]
        print "\n<<top 10 table data for {}\n".format(table)
        return data

    def for_venn(self):
        '''
        venn图分析页面发起创建otu集
        '''
        searcher = {
            'venn_id': ObjectId(self.option('table_id')),
        }
        include = ['otu_names', 'category_name']
        data = self.get_table_data('sg_otu_venn_graph', searcher, include)
        target = []
        others = []
        labels = self.option('label').replace("&amp;", "&")
        labels = re.split(r'\s*&\s*|\s*only', labels)
        for one in data:
            if one['category_name'] in labels:
                #target.extend(one['otu_names'].split(','))
                if target:                                   # bu zhaozhigang 20210406
                    target = list(set(target).intersection(set(one['otu_names'].split(','))))
                else:
                    target = one['otu_names'].split(',')
            else:
                others.extend(one['otu_names'].split(','))
        uniq_items = set(target) - set(others)
        if self.option('species_name'):
            uniq_items = filter(
                lambda x: self.option('species_name').lower() in x.lower(),
                uniq_items)
        return uniq_items, [self.level_dict[self.main_info['level_id']], ]

    def for_diff(self):
        '''
        两样本，两组，多组比较分析页面发起创建otu集
        '''
        searcher = {
            'species_check_id': ObjectId(self.option('table_id')),
            'pvalue': {'$lt': self.option('pvalue')},
            'corrected_pvalue': {'$lt': self.option('qvalue')},
        }
        include = ['species_name', ]
        diff_data = self.get_table_data('sg_species_difference_check_detail', searcher, include)
        names = map(lambda x: x['species_name'], diff_data)
        if self.option('species_name'):
            names = filter(
                lambda x: self.option('species_name').lower() in x.lower(),
                names)
        return names, [self.level_dict[self.main_info['level_id']], ]

    def for_lefse(self):
        '''
        lefse比较分析页面发起创建otu集
        '''
        searcher = {
            'species_lefse_id': ObjectId(self.option('table_id')),
            'lda': {'$ne': ''},
            'pvalue': {'$ne': '-'}
        }
        include = ['species_name', 'lda', 'pvalue']
        lefse_data = self.get_table_data('sg_species_difference_lefse_detail',
                                         searcher, include)
        names = [x['species_name'].split('.')[-1] for x in lefse_data
                 if float(x['lda']) > self.option('lda')and float(x['pvalue']) < self.option('pvalue')]
        levels = [l for k, l in self.level_dict.items()
                  if self.params['start_level'] <= k <= self.params['end_level']]
        return names, levels

    def for_randomforset(self):
        '''
        randomforset分析页面发起创建otu集
        '''
        searcher = {
            'randomforest_id': ObjectId(self.option('table_id')),
        }
        include = ['species_name', 'accuracy']
        rf_data = self.get_table_data('sg_randomforest_species_bar',
                                      searcher, include)
        rf_data = sorted(rf_data, key=lambda x: x['accuracy'], reverse=True)
        names = map(lambda x: x['species_name'], rf_data[:self.option('top')])
        return names, [self.level_dict[self.main_info['level_id']], ]
