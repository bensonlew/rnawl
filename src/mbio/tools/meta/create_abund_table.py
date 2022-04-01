# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last_modifiy:2018.06.27
#last_modify:2018.11.05 qc.zhang

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import pandas as pd


class CreateAbundTableAgent(Agent):
    """
    生成各项分析的不同数据库的标准丰度表格
    """
    def __init__(self, parent):
        super(CreateAbundTableAgent, self).__init__(parent)
        options = [
            {"name": "anno_table", "type": "infile", "format": "meta.profile"},  # 各数据库的注释表格
            {"name": "geneset_table", "type": "infile", "format": "meta.otu.otu_table,sequence.profile_table"},
            {"name": "gene_list", "type": "infile", "format": "meta.profile"},
            {"name": "level_type", "type": "string", "default": ""},
            # 注释表的字段，eg：Pathway，Level1，Level2
            {"name": "level_type_name", "type": "string", "default": ""},
            # 注释表字段的具体levelname，eg：Level1下Metabolism(对应的KO)
            {"name": "lowest_level", "type": "string", "default": ""},  # 注释表数据库对应的最低分类，eg：KEGG的ko
            {"name": "out_table", "type": "outfile", "format": "meta.otu.otu_table"},
            {"name": "anno_type", "type": "string","default": ""}, #当anno_type=nr时，用于增加物种以上的分类；
            {"name": "graphic_type", "type": "string", "default": ""},  # bar,heatmap,circos，bubble
            {"name": "level_color", "type": "string", "default": ""},   #add bu qingchen.zhang@20181202主要是用于增加颜色水平的判断
            #{"name": "project_name", "type": "string", "default": ""},  #add by qingchen.zhang@20181204 用于区分与meta的项目
        ]
        self.add_option(options)
        self._memory_increase_step = 50  # 每次重运行增加50G内存 add by GHD @ 20190215

    def check_options(self):
        if not self.option("anno_table").is_set and not self.option("gene_list").is_set:
            raise OptionError("请传入注释表格或者基因list文件！", code="32703801")
        if not self.option("geneset_table").is_set:
            raise OptionError("请传入基因丰度文件！", code="32703802")
        if self.option("level_type") == "" and self.option("anno_table").is_set:
            raise OptionError("请提供筛选的level水平", code="32703803")
        if self.option("level_type_name") != "":
            if self.option("lowest_level") == "":
                raise OptionError("请提供数据库对应的最低分类", code="32703804")

    def set_resource(self):
        self._cpu = 5
        self._memory = '20G'  # 内存5G增加到10G增加到15G by GHD @20180309 改回 by GHD @ 20180428
        # tmp_mem = 10 * (self._rerun_time + 1)  # 每次重运行的内存增加10G by GHD @ 20180320
        # self._memory = '%sG' % tmp_mem
        # self.logger.info('abund_table use memory : ' + self._memory)

    def end(self):
        super(CreateAbundTableAgent, self).end()


class CreateAbundTableTool(Tool):
    def __init__(self, config):
        super(CreateAbundTableTool, self).__init__(config)

    def create_abund_table(self):
        geneset_table_path = self.option("geneset_table").prop["path"]
        self.logger.info("geneset_table is : " + geneset_table_path)
        geneset_table = pd.read_csv(geneset_table_path, sep='\t', header=0, chunksize = 2000000)
        new_otu_file_path = os.path.join(self.output_dir, "new_abund_table.xls")
        rename = {}
        if self.option("gene_list").is_set:
            gene_list_path = self.option("gene_list").prop["path"]
            self.logger.info("gene_list is : " + gene_list_path)
            gene_list = pd.read_csv(gene_list_path, sep='\t', header=0)
            gene_list = list(gene_list['GeneID'])
            out = []
            for ck in geneset_table:
                tmp = ck[ck['GeneID'].isin(gene_list)]
                out.append(tmp)
            gene_table = pd.concat(out)
            del gene_table["Total"]
            gene_table = gene_table[list((gene_table > 0).any(axis=1))]
            if len(gene_table) < 1:
                self.set_error('在所选gene集参数下数据为空，请重新设置该参数!', code="32703801")
                self.set_error('在所选gene集参数下数据为空，请重新设置该参数!', code="32703804")
            gene_table.to_csv(new_otu_file_path, sep="\t", index=False)
        else:
            anno_table_path = self.option("anno_table").prop["path"]
            self.logger.info("anno table is : " + anno_table_path)
            anno_table = pd.read_table(anno_table_path, sep='\t', header=0)
            if 'Class_description' in anno_table.columns and self.option("anno_type") not in ["card", "ardb"]:  # 注释更新后，传参中的class 对应的Class_description内容
                if self.option("level_type") == 'Class':
                    self.option("level_type", 'Class_description')
                    rename["Class_description"] = "Class"
                if self.option("lowest_level") == 'Class':
                    self.option("lowest_level", 'Class_description')
                    rename["Class_description"] = "Class"
            a = pd.DataFrame(anno_table)
            a = a.fillna('-')
            if self.option('graphic_type')=="heatmap" or self.option('graphic_type') == "bubble": #add by qingchen.zhang@20181205 用于增加颜色水平的筛选
                if self.option('level_color') != "":
                    if (self.option("level_type")=='Function' and self.option('level_color')=='Category') or ((self.option("level_type")=='Level2' and self.option('level_color')=='Level1') or (self.option("level_type")=='Level3' and self.option('level_color')=='Level1') or (self.option("level_type")=='Level3' and self.option('level_color')=='Level2')):
                        a[self.option('level_type')] = a[self.option('level_type')].apply(lambda x: x.replace('; ', ';'))
                        a[self.option('level_color')] = a[self.option('level_color')].apply(lambda x: x.replace('; ', ';'))
                        type = a[self.option("level_type")].str.split(';',expand=True).stack().reset_index(level=1, drop=True)
                        color = a[self.option('level_color')].str.split(';',expand=True).stack().reset_index(level=1,drop=True)
                        table = pd.concat([type,color], axis=1)
                        table.columns = [self.option("level_type"), self.option('level_color')]
                        table2 = table[self.option("level_type")].str.cat(table[self.option('level_color')],sep='|')
                        result = table2.groupby(table2.index).apply(';'.join)
                        self.logger.info("result表的head：%s," % result.ix[0,:])
                        del a[self.option("level_type")]
                        a = pd.concat([a,result],axis=1)
                        if self.option("level_type_name") != "":
                            a = a.ix[:, ["#Query", self.option("lowest_level"), self.option("level_type")]]
                            a.columns = ["GeneID", self.option("lowest_level"), self.option("level_type")]
                        else:
                            a = a.ix[:, ["#Query", self.option("level_type")]]
                            a.columns = ["GeneID", self.option("level_type")]
                    else:
                        if self.option("level_type_name") != "":
                            a = a.ix[:, ["#Query", self.option("lowest_level"), self.option("level_type"), self.option('level_color')]]
                            a.columns = ["GeneID", self.option("lowest_level"), self.option("level_type"), self.option('level_color')]
                        else:
                            a = a.ix[:, ["#Query", self.option("level_type"), self.option('level_color')]]
                            a.columns = ["GeneID", self.option("level_type"), self.option('level_color')]
                else:
                    if self.option("level_type_name") != "":
                        a = a.ix[:, ["#Query", self.option("lowest_level"), self.option("level_type")]]
                        a.columns = ["GeneID", self.option("lowest_level"), self.option("level_type")]
                    else:
                        a = a.ix[:, ["#Query", self.option("level_type")]]
                        a.columns = ["GeneID", self.option("level_type")]
            else:
                if self.option("level_type_name") != "":
                    a = a.ix[:, ["#Query", self.option("lowest_level"), self.option("level_type")]]
                    a.columns = ["GeneID", self.option("lowest_level"), self.option("level_type")]
                else:
                    a = a.ix[:, ["#Query", self.option("level_type")]]
                    a.columns = ["GeneID", self.option("level_type")]
            # 针对大数据优化
            out = []
            for ck in geneset_table:
                out.append(self.by_anno_table(a, ck))

            if self.option("level_type_name") != "":
                level_type_abund_table = pd.concat(out).groupby(self.option("lowest_level")).sum()
            else:
                level_type_abund_table = pd.concat(out).groupby(self.option("level_type")).sum()
            if len(level_type_abund_table)< 1:  # add by zhujuan 当基因集和注释表下数据为空空时，报错并提示 20180124
                self.set_error('在所选基因集和注释表下数据为空，请重新设置基因集或注释表参数!', code="32703802")
                self.set_error('在所选基因集和注释表下数据为空，请重新设置基因集或注释表参数!', code="32703805")

            if rename:
                level_type_abund_table = level_type_abund_table.rename(columns=rename)
            level_type_abund_table.to_csv(new_otu_file_path, sep="\t")

    def by_anno_table(self, anno_table, geneset_table):
        a = anno_table.copy()
        b = geneset_table.copy()
        del b["Total"]
        abund = a.merge(b, on='GeneID', how='inner')
        self.logger.info(">>>>>>>>>>>>>>>>>>>>>>>>")
        self.logger.info(self.option("level_type"))
        self.logger.info(self.option("level_type_name"))
        self.logger.info(abund.head())
        if len(abund)< 1:
            return None  # 当前 dataframe为空的话直接返回
        abund[self.option('level_type')] = abund[self.option('level_type')].apply(lambda x: x.replace('; ', ';'))
        level_type_abund = abund.drop(self.option("level_type"), axis=1).join(
            abund[self.option("level_type")].astype('str').str.split(';', expand=True).stack().reset_index(
                level=1, drop=True).rename(self.option("level_type")))
        if self.option("level_type_name") != "":
            if self.option("lowest_level") in ["Pathway"]:
                abund[self.option('lowest_level')] = abund[self.option('lowest_level')].apply(lambda x: x.replace('; ', ';'))
                pathway = abund[self.option("lowest_level")].str.split(";", expand=True).stack().reset_index(
                    level=1, drop=True).rename(self.option("lowest_level"))
                level_type_abund[self.option("lowest_level")] = list(pathway)
                level_type_abund = level_type_abund[level_type_abund[self.option("level_type")].astype('str') == self.option(
                    "level_type_name")]
                level_type_abund = level_type_abund[level_type_abund[self.option("lowest_level")].astype('str') != "-"]
                level_type_abund = level_type_abund.drop_duplicates()
            else:
                level_type_abund = level_type_abund[level_type_abund[self.option("level_type")].astype('str') == self.option(
                    "level_type_name")]
                level_type_abund = level_type_abund[level_type_abund[self.option("lowest_level")].astype('str') != "-"]
            level_type_abund_table = level_type_abund.groupby(self.option("lowest_level")).sum()
        else:
            if self.option('level_color') != "":#分为一对一、一对多和多对多关系
                if (self.option("level_type")=='Function' and self.option('level_color')=='Category') or ((self.option("level_type")=='Level2' and self.option('level_color')=='Level1') or (self.option("level_type")=='Level3' and self.option('level_color')=='Level1') or (self.option("level_type")=='Level3' and self.option('level_color')=='Level2')):
                    level_type_abund = level_type_abund[level_type_abund[self.option("level_type")].astype('str') != "-|-"]
                    level_type_abund_table = level_type_abund.groupby(self.option("level_type")).sum()
                else:
                    level_type_abund = level_type_abund[level_type_abund[self.option("level_type")].astype('str') != "-"]
                    self.logger.info("level_type_abund表的head：%s" % level_type_abund.ix[0,:])
                    new_level_type = level_type_abund[self.option("level_type")] + "|" + level_type_abund[self.option('level_color')]
                    del level_type_abund[self.option("level_type")]
                    level_type_abund = pd.concat([level_type_abund,new_level_type],axis=1)
                    level_type_abund.rename(columns = {0: self.option("level_type")}, inplace=True)
                    self.logger.info("level_type_abund表的head：%s" % level_type_abund.ix[0,:])
                    #level_type_abund = abund.drop(self.option("level_type"), axis=1).join(level_type_abund['type_color'].rename(self.option("level_type")))
                    level_type_abund_table = level_type_abund.groupby(self.option("level_type")).sum()
            else:
                level_type_abund = level_type_abund[level_type_abund[self.option("level_type")].astype('str') != "-"]
                level_type_abund_table = level_type_abund.groupby(self.option("level_type")).sum()
        level_type_abund_table = level_type_abund_table.ix[list((level_type_abund_table > 0).any(axis=1))]
        if len(level_type_abund_table) < 1:  # add by zhujuan 当分类水平下数据都空时，报错并提示 20171124
            self.set_error('在所选物种/功能的分类参数下数据为空，请重新设置该参数!', code="32703803")
            self.set_error('在所选物种/功能的分类参数下数据为空，请重新设置该参数!', code="32703806")
        return level_type_abund_table

    def set_output(self):
        self.logger.info('开始设置输出结果文件')
        try:
            self.option('out_table', os.path.join(self.output_dir, "new_abund_table.xls"))
            self.logger.info(self.option('out_table').prop['path'])
        except Exception as e:
            self.set_error("输出结果文件异常——%s", variables=(e), code="32703807")

    def run(self):
        super(CreateAbundTableTool, self).run()
        self.create_abund_table()
        '''
        try:
            self.create_abund_table()
        except:
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180322
        '''
        self.set_output()
        self.end()


