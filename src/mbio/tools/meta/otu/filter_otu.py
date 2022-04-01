# -*- coding: utf-8 -*-
# __author__ = 'xuting'
# modified: hesheng  # 对物种筛选进行重构
# lastmodified: liulinmeng # 对筛选条件中的物种名的特殊字符进行转义，以完成相应筛选功能
import json
import re
import os
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
from bson.objectid import ObjectId


class FilterOtuAgent(Agent):
    """
    根据传入的json，对一张OTU表进行过滤，过滤的条件有三种:
    species_filter: 物种过滤，用于保留或者滤去特定的物种
    sample_filter: 用于滤去在x个样本中序列数小于y的OTU
    reads_filter: 用于滤去序列数小于x的OTU
    """
    def __init__(self, parent):
        super(FilterOtuAgent, self).__init__(parent)
        options = [
            {"name": "in_otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入的OTU文件
            {"name": "filter_json", "type": "string", "default": ""},  # 过滤规则
            {"name": "out_otu_table", "type": "outfile", "format": "meta.otu.otu_table"}  # 输出的结果OTU表
        ]
        self.add_option(options)
        self.step.add_steps("filter_otu")

    def start_filter_otu(self):
        self.step.filter_otu.start()
        self.step.update()

    def end_filter_otu(self):
        self.step.filter_otu.end()
        self.step.update()

    def check_options(self):
        """
        参数检测
        """
        if not self.option("in_otu_table").is_set:
            raise OptionError("输入的OTU文件不能为空", code="32704701")
        if self.option("filter_json") == "":
            raise OptionError("输入的筛选条件不能为空", code="32704702")

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["out_otu.xls", "xls", "结果OTU表格"]
        ])
        super(FilterOtuAgent, self).end()

    def set_resource(self):
        """
        设置所需的资源
        """
        self._cpu = 2
        self._memory = "10G"


class FilterOtuTool(Tool):
    def __init__(self, config):
        super(FilterOtuTool, self).__init__(config)
        self.otu_sample_dict = self.option("in_otu_table").extract_info()
        self.json = list()
        self.otu_json = list()# 将整个OTU表读入，生成一个列表，方便后面进行过滤操作
        with open(self.option("in_otu_table").prop["path"], 'rb') as r:
            self.otu_json = r.readlines()
        self.otu_head = self.otu_json.pop(0)  # OTU表的表头
        self.LEVEL = {
            1: "d__", 2: "k__", 3: "p__", 4: "c__", 5: "o__",
            6: "f__", 7: "g__", 8: "s__", 9: "otu"
        }
        self.keep_list = list()  # 处理物种筛选中的保留的逻辑
        Config().DBVersion = self.config.DBVersion
        self.client = Config().get_mongo_client(mtype='meta')  # 用于获取otuset主表
        self.db = self.client[Config().get_mongo_dbname(mtype='meta')]

    def flat_cond(self):
        filter_json = json.loads(self.option("filter_json"))
        self.logger.info(filter_json)
        for cond in sorted(filter_json):
            if cond['name'] == "species_filter":
                for va in cond['value'].split(','):
                    c = {k: v for k, v in cond.items()}
                    c['value'] = va
                    filter_json.append(c)
                filter_json.remove(cond)
        self.logger.info(filter_json)
        return filter_json

    def filter_otu_table(self):
        """
        筛选 OTU表中的物种，样本和reads

        params:
        """
        filter_json = self.flat_cond()
        sp_condition_keep = []
        sp_condition_remove = []
        sam_condition = []
        reads_condition = []
        otuset_condition = []
        for cond in filter_json:
            if cond['name'] == "species_filter" and cond["type"] == "keep":
                cond['value'] = self.escapeExprSpecialWord(cond['value'])
                sp_condition_keep.append(cond)
            elif cond['name'] == "species_filter" and cond["type"] == "remove":
                cond['value'] = self.escapeExprSpecialWord(cond['value'])
                sp_condition_remove.append(cond)
            elif cond["name"] == "sample_filter":
                sam_condition.append(cond)
            elif cond['name'] == "reads_filter":
                reads_condition.append(cond)
            elif cond['name'] == "otuset_filter":
                otuset_condition.append(cond)
            else:
                pass
        if sp_condition_keep:
            self.otu_json = self.keep_species(sp_condition_keep)
        if sp_condition_remove:
            self.otu_json = self.remove_species(sp_condition_remove)
        for i in sam_condition:
            self.filter_samples(i)
        for i in reads_condition:
            self.filter_reads(i)
        #for i in otuset_condition:
        self.filter_otuset(otuset_condition)
        self.client.close()

    def escapeExprSpecialWord(self, words):
        '''
        转义正则表达式中的特殊字符
        :params words: 需要转义的字符串
        :return : 转义之后的字符串
        '''
    
        fbsArr = [ "\\", "$", "(", ")", "*", "+", ".", "[", "]", "?", "^", "{", "}", "|"]  
        #fbsArr = ["(", ")", "*", "+", ".", "[", "]", "?"]  
        for key in fbsArr:    
            if words.find(key)>=0:  
                words = words.replace(key, "\\"+key); 
        return words

    def keep_species(self, conditions, fuzzy=False):
        """
        保留特定物种的OTU, 多条件时选交集

        :params conditions: 条件列表
        :params fuzzy: 是否模糊匹配(不区分大小写)
        :return : 筛选出来的OTU列表
        """
        temp_otus = []
        for cond in conditions:
            if not fuzzy:
                cond['pattern'] = re.compile('^{}$'.format(cond['value'].strip()))
            else:
                cond['pattern'] = re.compile('{}'.format(cond['value'].strip()), flags=re.IGNORECASE)
        for i in self.otu_json:
            for cond in conditions:
                level = int(cond["level_id"]) - 1
                sp_name = re.split(r';', (re.split(r'\t', i, maxsplit=1)[0]))
                sp_name = sp_name[level].strip()
                if cond['pattern'].search(sp_name):
                    temp_otus.append(i)
                    break
        return temp_otus

    def remove_species(self, conditions, fuzzy=False):
        """
        去除特定物种的OTU, 多条件时，任何条件去除即去除

        :params conditions: 条件列表
        :params fuzzy: 是否模糊匹配(不区分大小写)
        :return : 筛选出来的OTU列表
        """
        temp_otus = []
        for cond in conditions:
            if not fuzzy:
                cond['pattern'] = re.compile('^{}$'.format(cond['value'].strip()))
            else:
                cond['pattern'] = re.compile('{}'.format(cond['value'].strip()), flags=re.IGNORECASE)
        for i in self.otu_json:
            for cond in conditions:
                level = int(cond["level_id"]) - 1
                sp_name = re.split(r';', (re.split(r'\t', i, maxsplit=1)[0]))
                sp_name = sp_name[level].strip()
                if cond['pattern'].search(sp_name):
                    break
            else:
                temp_otus.append(i)
        return temp_otus

    def filter_samples(self, my_json):
        """
        保留或去除至少在x个样本中序列数大于y的物种(OTU)
        """
        tmp_list = self.otu_json[:]
        my_c = 0
        for line in self.otu_json:
            otu = line.split("\t")[0]  # otu即是第一列，也是self.otu_sample_dict的第一个key值
            count = 0
            for sp in self.otu_sample_dict[otu]:
                if self.otu_sample_dict[otu][sp] >= int(my_json["reads_num"]):
                    count += 1
            if count < int(my_json["sample_num"]) and my_json['type'] == 'keep':
                tmp_list.remove(line)
            elif count >= int(my_json["sample_num"]) and my_json['type'] == 'remove':
                tmp_list.remove(line)
            else:
                my_c += 1
        self.otu_json = tmp_list[:]

    def filter_reads(self, my_json):
        """
        保留或去除序列数总和大于x的物种(OTU)
        """
        tmp_list = self.otu_json[:]
        for line in self.otu_json:
            otu = line.split("\t")[0]  # otu即是第一列，也是self.otu_sample_dict的第一个key值
            summary = 0
            for sp in self.otu_sample_dict[otu]:
                summary += self.otu_sample_dict[otu][sp]
            if summary < int(my_json["reads_num"]) and my_json['type'] == 'keep':
                tmp_list.remove(line)
            elif summary >= int(my_json["reads_num"]) and my_json['type'] == 'remove':
                tmp_list.remove(line)
        self.otu_json = tmp_list[:]

    def filter_otuset(self, my_json):
        '''
        保留或去除包含当前otu集的
        '''
        keep_otuset = []
        remove_otuset = []
        for i in my_json:
            otuset_id = ObjectId(i['otuset'])
            otuset_table = self.db['sg_otuset'].find_one({'_id': otuset_id, 'status': 'end'})
            otuset = otuset_table['otuset_list']
            for x in otuset:
                if i['type'] == 'keep':
                    if x not in keep_otuset:
                        keep_otuset.append(x)
                else:
                    if x not in keep_otuset:
                        remove_otuset.append(x)
        for xx in remove_otuset:
            if xx in keep_otuset:
                remove_otuset.remove(xx)

        tmp_list = self.otu_json[:]
        if keep_otuset:
            [tmp_list.remove(i) for i in self.otu_json
             if i.split('\t')[0].rpartition(';')[-1] not in keep_otuset]
        if remove_otuset:
            [tmp_list.remove(i) for i in self.otu_json
             if i.split('\t')[0].rpartition(';')[-1] in remove_otuset]
        self.otu_json = tmp_list[:]

    def run(self):
        super(FilterOtuTool, self).run()
        self.filter_otu_table()
        # if len(self.otu_json) == 0:
        #    raise Exception("过滤之后的结果OTU是空的, 请查看过滤的条件是否正确！")
        with open(os.path.join(self.output_dir, "filter_otu.xls"), "wb") as w:
            w.write(self.otu_head)
            for line in self.otu_json:
                w.write(line)
        self.option("out_otu_table").set_path(os.path.join(self.output_dir, "filter_otu.xls"))
        self.logger.info("OTU过滤完成")
        self.end()
