# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
from  mainapp.models.mongo.metabolome import  Metabolome
import pandas as pd
from biocluster.config import Config


class EnrichTopoModule(Module):
    def __init__(self, work_id):
        super(EnrichTopoModule, self).__init__(work_id)
        options = [
            {"name": "metabset", "type": "infile", "format": "metabolome.metabset"},  # 代谢集文件
            {"name": "anno_overview", "type": "infile", "format": "metabolome.overview_table"},  # 代谢集总览表
            {"name": "ko_overview", "type": "infile", "format": "sequence.profile_table"},  # 代谢集总览表ko表 by ghd @20191015
            {"name": "correct", "type": "string", "default": "BH"},  # 多重检验校正方法
            {"name": "bg", "type": "string", "default": "project"},
            # 背景，project:本项目鉴定到的代谢物合集; species:本物种全部代谢物合集; kegg:KEGG数据库全部代谢物合集
            {"name": "species", "type": "string", "default": "all"},
            {"name": "method", "type": "string", "default": "rbc"}, # rbc or  rod or none
            {"name": "task_id", "type": "string", "default": ""}, # 用于区分新老工作流
            {"name": "database_version", "type": "string", "default": ""}, # 用于区分新老工作流
        ]
        self.add_option(options)
        self.enrich_tool = self.add_tool("metabolome.metabset.enrich")
        self.topo_tool = self.add_tool("metabolome.metabset.topology")



    def check_options(self):

        if not self.option('metabset').is_set:
            raise OptionError("必须设置代谢集", code="34701001")
        if not self.option('anno_overview').is_set:
            raise OptionError("必须设置代谢总览表", code="34701002")
        if self.option("correct") not in ["BH", "BY", "bonferroni", "holm"]:
            raise OptionError("矫正参数不在范围内，错误参数值：%s", variables=(self.option('correct')), code="34701003")

        return True

    def run_enrich(self):
        opts = {
            'anno_overview': self.option("anno_overview"),
            'ko_overview': self.option("ko_overview"),  # by ghd @20191015
            'metabset': self.option("metabset"),
            'correct': self.option('correct'),
            'bg': self.option('bg'),
            'species' : self.option("species"),
            "version": self.version_value
        }
        self.enrich_tool.set_options(opts)
        self.enrich_tool.run()

    def run_topology(self):
        set_list = pd.read_table(self.option("metabset").path,sep='\t',header=-1)[0].tolist()
        data = pd.read_table(self.enrich_tool.work_dir + '/2compound.info',sep='\t',header=0)
        data2= data[data['metab'].apply(lambda x:x in set_list)]
        new_compound = self.work_dir+'/select_2compound.info'
        data2.to_csv(new_compound,sep='\t',index=False)
        opts = {
            'compound_table': new_compound,
            'method' : self.option('method'),
            'backgroup' : os.path.join(self.enrich_tool.output_dir,'DE.list.check.kegg_enrichment.xls'),
            "version": self.version_value
        }
        self.topo_tool.set_options(opts)
        self.topo_tool.run()


    def set_output(self):

        self.move_dir(self.enrich_tool.output_dir,self.output_dir)
        if self.option('method') != "none":
            self.move_dir(self.topo_tool.output_dir,self.output_dir+'/topology_png')
        self.end()

    def run(self):
        super(EnrichTopoModule, self).run()
        self.metablome = Metabolome()
        self.metablome._config = Config()
        self.logger.info("###database_version: {}".format(self.option("database_version")))
        self.version_value = self.option("database_version") or self.metablome.find_version_from_task(self.option("task_id"))
        if self.option('method') == "none":
            self.enrich_tool.on("end",self.set_output)
        else:
            self.enrich_tool.on("end",self.run_topology)
            self.topo_tool.on("end",self.set_output)
        self.run_enrich()


    def end(self):
        super(EnrichTopoModule, self).end()

    def move_dir(self,ori_dir,t_dir):
        if not os.path.exists(t_dir):
            os.mkdir(t_dir)
        if not os.path.exists(ori_dir):
            self.logger.set_error('%s 不存在'%ori_dir)
        files = os.listdir(ori_dir)
        for f in files:
            o_file = os.path.join(ori_dir,f)
            t_file = os.path.join(t_dir,f)
            if os.path.exists(t_file):
                os.remove(t_file)
            os.link(o_file, t_file)

