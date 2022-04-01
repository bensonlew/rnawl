# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import os, shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.metabolome.common import Relation
from  mainapp.models.mongo.metabolome import  Metabolome
import pandas as pd
import re
import glob
from biocluster.config import Config

class AnnoKeggpAgent(Agent):
    """
    代谢组kegg通路注释
    version: v1.0
    author: guhaidong
    last_modify: 2018.06.15
    """

    def __init__(self, parent):
        super(AnnoKeggpAgent, self).__init__(parent)
        options = [
            {"name": "metab_table", "type": "infile", "format": "metabolome.metab_desc"},  # 预处理的metab_desc结果
            {"name": "organism", "type": "string", "default": "False"},  # 参考物种
            {"name": "level_out", "type": "outfile", "format": "sequence.profile_table"},  # 注释层级结果
            {"name": "stat_out", "type": "outfile", "format": "sequence.profile_table"},  # 注释统计结果
            {"name": "task_id", "type": "string", "default": ""}, # 用于区分新老工作流
            {"name": "database_version", "type": "string", "default": ""},
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
        检查参数
        :return:
        """
        if not self.option('metab_table'):
            raise OptionError('必须输入代谢物详情', code="34700601")
        return True

    def set_resource(self):
        """
        设置所需资源，需在子类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "8G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", ""],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(AnnoKeggpAgent, self).end()


class AnnoKeggpTool(Tool):
    def __init__(self, config):
        super(AnnoKeggpTool, self).__init__(config)
        self.level_out = os.path.join(self.output_dir, 'level.xls')
        self.stat_out = os.path.join(self.output_dir, 'stat.xls')
        self.metablome = Metabolome()
        self.metablome._config = Config()
        self.version_value = ''
        if self.option("database_version"):
            self.version_value = self.option("database_version")
        elif self.option("task_id") in ['']:## 工作流
            self.version_value = 'v2021.09.18'
        else:# 交互分析
            self.version_value = self.metablome.find_version_from_task(self.option("task_id"))
        if self.version_value in ["kegg"]:
            self.html_db = self.config.SOFTWARE_DIR + "/database/KEGG/map_html/"
        elif self.version_value in ['v94.2', '94.2']: ## 94.2
            self.html_db = self.config.SOFTWARE_DIR + "/database/Annotation/all/KEGG/version_202007_meta/html/"
        else:
            self.html_db = self.config.SOFTWARE_DIR + "/database/Annotation/all/KEGG/{}/html/".format(self.version_value)

    def run(self):
        """
        运行
        :return:
        """
        super(AnnoKeggpTool, self).run()
        if self.version_value not in [""]:
            from mbio.packages.metabolome.anno_keggp import AnnoKeggp
            self.orginisms = self.option('organism').strip("_")
            obj = AnnoKeggp()
            self.logger.info("from mbio.packages.metabolome.anno_keggp import AnnoKeggp")
        else:## 老的库脚本
            self.orginisms = self.option('organism')
            from mbio.packages.metabolome.anno_keggp2 import AnnoKeggp
            obj = AnnoKeggp()
            self.logger.info("from mbio.packages.metabolome.anno_keggp2 import AnnoKeggp")
        try:
            self.logger.info(self.option('metab_table').path)
            self.logger.info(self.option('organism'))
            self.logger.info(self.level_out)
            self.logger.info(self.stat_out)
            self.logger.info(self.html_db)

            obj.run(self.option('metab_table').path, self.orginisms, self.level_out, self.stat_out,html_db=self.html_db,color='red',db_v=self.version_value)
        except Exception,e:
            self.logger.error(e)
            self.set_error(e, code="34700601")
        self.id_to_name()
        self.set_output()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        self.option('level_out').set_path(self.level_out)
        self.option('stat_out').set_path(self.stat_out)
        self.logger.info("设置anno_keggp分析结果目录成功")
        self.logger.info("替换mark文件的title值")
        self.change_mark_title(self.output_dir+'/pathway_img')
        self.end()

    def id_to_name(self):
        # 增加一列代谢物名称
        # 影响pathway_img 目录位置，故改用move的方式
        old_level = os.path.join(self.work_dir, 'level.xls')
        old_stat = os.path.join(self.work_dir, 'stat.xls')
        if os.path.exists(old_level):
            os.remove(old_level)
        if os.path.exists(old_stat):
            os.remove(old_stat)
        shutil.move(self.level_out, old_level)
        shutil.move(self.stat_out, old_stat)
        table_path = self.option("metab_table").prop["path"]
        self.metab_trans = Relation()
        self.logger.info(table_path)
        map_table, id_name_dict = self.metab_trans.get_dataframe_and_dict(table_path)
        self.metab_trans.add_metabolites_column(id_name_dict, old_level, self.level_out)
        self.metab_trans.add_metabolites_column(id_name_dict, old_stat, self.stat_out)

    def change_mark_title(self,out_dir):
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
                            if len(self.com_metab_map[r]) > 0:
                                tmp = r + ' (' +';'.join(self.com_metab_map[r])+')'
                                tmp_list.append(tmp)
                    if len(tmp_list) != 0 :
                        fd['title'][i] = ' '.join(tmp_list)
            fd.to_csv(f,sep='\t',index=False)


    def get_compound_metab_map(self):
        metab_table = self.option('metab_table').path
        data = pd.read_table(metab_table,sep='\t')
        compound = data['KEGG Compound ID'].tolist()
        metab = data['Metabolite'].tolist()
        com_metab = zip(compound,metab)
        com_metab_map = {}
        for c,m in com_metab:
            if re.match('^metab_\d*$',m):
                continue
            spc = c.split(';')
            for each_c in spc:
                if each_c in com_metab_map:
                    if m not in com_metab_map[each_c]:
                        com_metab_map[each_c].append(m)
                else:
                    com_metab_map[each_c] = [m]
        return com_metab_map






