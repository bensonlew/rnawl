# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os, shutil
from mbio.packages.metabolome.common import Relation
from  mainapp.models.mongo.metabolome import Metabolome
import glob
import re
import pandas as pd
from biocluster.config import Config
from pymongo.errors import ServerSelectionTimeoutError
from pymongo.errors import NetworkTimeout
import time

class KeggpAgent(Agent):
    """
    代谢集kegg代谢物注释
    last_modify: 2018.6.19
    """

    def __init__(self, parent):
        super(KeggpAgent, self).__init__(parent)
        options = [
            {"name": "metabset", "type": "infile", "format": "metabolome.mul_metabset"},  # 代谢集文件
            {"name": "anno_overview", "type": "infile", "format": "sequence.profile_table"},  # 代谢集总览表ko表
            {"name": "level_out", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "trans_file", "type": "infile", "format": "sequence.profile_table"}, #  转化id使用
            {"name": "task_id", "type": "string", "default": ""}, # 用于区分新老工作流
            {"name": "database_version", "type": "string", "default": ""}, # 用于区分新老工作流
        ]
        self.add_option(options)
        self.step.add_steps("keggp")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.keggp.start()
        self.step.update()

    def stepfinish(self):
        self.step.keggp.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        # if not self.option("metabset").is_intersection:
        #     raise OptionError("两个基因集中没有交集")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(KeggpAgent, self).end()


class KeggpTool(Tool):
    def __init__(self, config):
        super(KeggpTool, self).__init__(config)
        self.image_magick = self.config.SOFTWARE_DIR + "/program/ImageMagick/bin/convert"
        self.metablome = Metabolome()
        self.metablome._config = Config()
        self.version_value = ''

        if self.option("database_version"):## 工作流
            self.version_value = self.option("database_version")
        elif self.option("task_id"):## 工作流
            self.search_mongo_num = 0
            self._search_mongo()

        if self.version_value == 'v2021.09.18':
            self.html_path = self.config.SOFTWARE_DIR + "/database/Annotation/all/KEGG/v2021.09.18/html/"
        elif self.version_value == 'v94.2':
            self.html_path = self.config.SOFTWARE_DIR + "/database/Annotation/all/KEGG/version_202007_meta/html/"
        elif self.version_value in [""]:
            self.html_path = self.config.SOFTWARE_DIR + "/database/KEGG/map_html/"
        else:
            self.html_path = self.config.SOFTWARE_DIR + "/database/Annotation/all/KEGG/{}/html/".format(self.version_value)

    def _search_mongo(self):
        self.search_mongo_num += 1
        try:
            self.logger.info('开始查询mongo库：%s次'%(self.search_mongo_num))
            self.version_value = self.metablome.find_version_from_task(self.option("task_id"))
        except (ServerSelectionTimeoutError, NetworkTimeout):
            if self.search_mongo_num < 11:
                time.sleep(5)
                self._search_mongo()
            else:
                self.set_error('数据库重查10次仍失败')


    def run(self):
        """
        运行
        :return:
        """
        super(KeggpTool, self).run()
        #try:
        self.logger.info('value：')
        self.logger.info(self.image_magick)
        self.logger.info(self.html_path)
        self.logger.info(self.option('metabset').path)
        self.logger.info(self.option('anno_overview').path)

        if self.version_value not in [""]:
            from mbio.packages.metabolome.metabsetp import Metabsetp
            obj = Metabsetp(self.image_magick, self.html_path)
        else:
            from mbio.packages.metabolome.metabsetp2 import Metabsetp
            obj = Metabsetp(self.image_magick, self.html_path)
        obj.run(self.option('metabset').path, self.option('anno_overview').path, self.output_dir)
        #except Exception,e:
        #    self.logger.error(e)
        #    raise Exception(e)
        self.logger.info("keggp运行完毕")
        self.id_to_name()
        self.set_output()
        self.end()

    def set_output(self):
        self.logger.info("设置结果目录")
        self.option('level_out',self.output_dir + '/level.xls')
        self.logger.info("修改mark的title")
        # try:
        self.change_mark_title(self.output_dir + '/pathway_img')
        # except Exception:
        #     self.logger.info("修改mark的title 失败")

    def id_to_name(self):
        # 增加一列代谢物名称
        # 影响pathway_img 目录位置，故改用move的方式, old_path为移动位置，new_path为最终位置
        annofile = self.option("anno_overview").prop["path"].replace("ko.xls","anno.xls")
        if self.option("trans_file").is_set:
            table_path = self.option("trans_file").prop["path"]
        elif os.path.exists(annofile):
            table_path = annofile
        else:
            raise Exception("please input trans id file")
        self.desc_table_path = table_path
        self.metab_trans = Relation()
        self.logger.info(table_path)
        map_table, id_name_dict = self.metab_trans.get_dataframe_and_dict(table_path, metab_name="metab")
        oldfiles = os.listdir(self.output_dir)
        for eachfile in oldfiles:
            oldfile = os.path.join(self.work_dir, eachfile)
            newfile = os.path.join(self.output_dir, eachfile)
            if os.path.isfile(newfile):
                if os.path.exists(oldfile):
                    os.remove(oldfile)
                shutil.move(newfile, oldfile)
                self.metab_trans.add_metabolites_column(id_name_dict, oldfile, newfile)


    def get_metabset(self):
        self.name_map_metabs = {}
        with open(self.option('metabset').path) as fr:
            for line in fr:
                sp = line.strip().split('\t')
                if len(sp) > 1:
                    self.name_map_metabs[sp[0]] = sp[1].split(',')


    def change_mark_title(self,out_dir):  #{cid:{'a':['metaname','metaname2'],'b':[]}}
        self.get_metabset()
        self.com_metab_map = self.get_compound_metab_map()
        files = glob.glob(out_dir+'/*html.mark')
        for f in files:
            fd = pd.read_table(f,sep='\t')
            for i in range(len(fd)):
                ori_title = fd['title'][i]
                if fd['shape'][i] != 'circle':
                    continue
                res = re.findall('C\d{5}',ori_title)
                if res:
                    tmp_list = []
                    for r in res:
                        if r in self.com_metab_map.keys():
                            setname_res = []
                            for set_name in sorted(self.name_map_metabs.keys()):
                                if len(self.com_metab_map[r][set_name]) !=0:
                                    set_name_value = set_name+'('+';'.join(self.com_metab_map[r][set_name])+')'
                                    setname_res.append(set_name_value)
                            if len(setname_res) > 0:
                                tmp_c = r + ' (' +','.join(setname_res)+')'
                                tmp_list.append(tmp_c)
                    if len(tmp_list) != 0 :
                        fd['title'][i] = ' '.join(tmp_list)   # C1111 (setname(AA;BB),setname(CC;DD)) C2222 []
            fd.to_csv(f,sep='\t',index=False)


    def get_compound_metab_map(self):
        metab_table = self.desc_table_path  #anno.xls
        self.logger.info(metab_table)
        data = pd.read_table(metab_table,sep='\t')
        compound = data['compound_id'].tolist()
        metab = data['metab'].tolist()
        metab_id = data['metab_id'].tolist()

        com_metab_id = zip(compound,metab,metab_id)
        com_metab_map = {}
        for c,m,mid in com_metab_id:
            if re.match('^metab_\d*$',m):
                continue

            spc = c.split(';')
            for each_c in spc:
                if each_c not in com_metab_map:
                    com_metab_map[each_c] ={}
                    for set_name in self.name_map_metabs:
                        com_metab_map[each_c][set_name] = []

                for set_name in self.name_map_metabs:
                    if mid in self.name_map_metabs[set_name]:
                        if m not in com_metab_map[each_c][set_name]:
                            com_metab_map[each_c][set_name].append(m)

        return com_metab_map

