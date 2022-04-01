# -*- coding: utf-8 -*-

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
from  mainapp.models.mongo.metabolome import  Metabolome
from biocluster.config import Config
import pandas as pd
import numpy as np

class EnrichClusterModule(Module):
    def __init__(self, work_id):
        super(EnrichClusterModule, self).__init__(work_id)
        options = [
            {"name": "metabset", "type": "infile", "format": "metabolome.mul_metabset"},  # 代谢集文件
            {"name": "anno_overview", "type": "infile", "format": "metabolome.overview_table"},  # 代谢集总览表
            {"name": "ko_overview", "type": "infile", "format": "sequence.profile_table"},  # 代谢集总览表ko表 by ghd @20191015
            {"name": "correct", "type": "string", "default": "BH"},  # 多重检验校正方法
            {"name": "bg", "type": "string", "default": "project"},
            # 背景，project:本项目鉴定到的代谢物合集; species:本物种全部代谢物合集; kegg:KEGG数据库全部代谢物合集
            {"name": "species", "type": "string", "default": "all"},
            ##clustr 的参数
            {"name": "select", "type": "string", "default": "pvalue"},
            {'name': 'sct', 'type': 'string', 'default': 'hierarchy'},  # 样本或代谢集聚类算法，hierarchy、无, kmeans
            {'name': 'scd', 'type': 'string', 'default': 'euclidean'},  # 样本或代谢集距离计算方式
            {'name': 'scm', 'type': 'string', 'default': 'single'},  # 样本或代谢集聚类方式,"complete","average","single"
            {'name': 'mct', 'type': 'string', 'default': 'hierarchy'},  # 代谢物或通路聚类算法，hierarchy、kmeans、无
            {'name': 'mcd', 'type': 'string', 'default': 'euclidean'},  # 代谢物或通路距离计算方式
            {'name': 'n_cluster', 'type': 'int', 'default': 1},  # 代谢物或通路聚类数目，kmeans时使用
            {'name': 'mcm', 'type': 'string', 'default': 'single'}, # 代谢物或通路聚类方式, hierarchy时使用，"complete","average","single"
            {"name": "scale", "type": "string", "default" :"T"},
            {"name": "task_id", "type": "string", "default": ""}, # 用于区分新老工作流
        ]
        self.add_option(options)
        self.cluster_tool = self.add_tool("metabolome.metabset.pathway_cluster")
        self.metabset_map_enrich = {}



    def check_options(self):

        if not self.option('metabset').is_set:
            raise OptionError("必须设置代谢集", code="34701001")
        if not self.option('anno_overview').is_set:
            raise OptionError("必须设置代谢总览表", code="34701002")
        if self.option("correct") not in ["BH", "BY", "bonferroni", "holm"]:
            raise OptionError("矫正参数不在范围内，错误参数值：%s", variables=(self.option('correct')), code="34701003")

        return True

    def split_metabset_file(self):
        self.metabset_names = []
        with open(self.option('metabset').path) as f:
            ##f.readline()
            for line in f:
                spline = line.strip().split("\t")
                metabset_name = spline[0]
                metabset_list = spline[1].split(',')
                self.metabset_names.append(metabset_name)
                with open(self.work_dir+'/'+metabset_name +'_metabset.list','w') as fw:
                    fw.write('\n'.join(metabset_list))

    def run_enrich(self):
        self.split_metabset_file()
        for metabset_name in self.metabset_names:
            metabset_file = self.work_dir+'/'+metabset_name +'_metabset.list'
            enrich_tool = self.add_tool("metabolome.metabset.enrich")
            opts = {
                'anno_overview': self.option("anno_overview"),
                'ko_overview': self.option("ko_overview"),
                'metabset': metabset_file,
                'correct': self.option('correct'),
                'bg': self.option('bg'),
                'species' : self.option("species"),
                "version": self.version_value
            }
            enrich_tool.set_options(opts)
            self.metabset_map_enrich[metabset_name] = enrich_tool

        self.on_rely(self.metabset_map_enrich.values(), self.run_cluster)
        for tool in self.metabset_map_enrich.values():
            tool.run()

    def run_cluster(self):
        params_name_map = {
            'pvalue' : 'P-Value',
            'qvaule' : 'Corrected P-Value',
            'in_pop' : 'Ratio_in_pop',
            'in_study' : 'Ratio_in_study'
        }
        select = params_name_map[self.option('select')]

        data_list = []
        pvalue_list = []
        desc_list = []
        metab_list = []
        for name in self.metabset_map_enrich:
            enrich_tool = self.metabset_map_enrich[name]
            enrich_result = os.path.join(enrich_tool.output_dir ,'DE.list.check.kegg_enrichment.xls')
            data = pd.read_table(enrich_result, sep='\t',index_col='ID')
            pick_data = data[select]
            data_list.append(pd.DataFrame({name: pick_data}))
            #if self.option('select') in ['pvalue','qvalue']:
            #    pick_data = data[select]
            #else:
            #    pick_data = map(eval,data[select]+'.00')  # eval(1/2.0) = 0.5
            #    data_list.append(pd.DataFrame({name: pick_data}))

            pick2_data = data['P-Value']
            desc_data = data[['Term', 'Database']]
            desc_list.append(pd.DataFrame(desc_data))

            pvalue_list.append(pd.DataFrame({name: pick2_data}))

            metab_list.append(pd.DataFrame({name+'_Study_num': data['#Study_num'],name+'_Metab_ids':data['Metab_ids']}))


        desc_all = pd.concat(desc_list, axis=0)
        desc_all = desc_all.reset_index()
        desc_all = desc_all.drop_duplicates(subset=['ID'])
        self.desc_xls = self.work_dir + '/pathway_desc.xls'
        desc_all.to_csv(self.desc_xls, sep='\t',index=False)
        metab_ids_xls = self.work_dir + '/metab_ids.xls'

        pick_all = pd.concat(data_list,axis=1, join='outer')
        pvalue_all = pd.concat(pvalue_list, axis=1, join='outer')
        metab_ids = pd.concat(metab_list, axis=1, join='outer')
        metab_ids.to_csv(metab_ids_xls, sep='\t')

        pvalue_all.fillna(1, inplace=True)
        if self.option('select') in ['pvalue','qvaule']:
            pick_all.fillna(1, inplace=True)
        else:
            pick_all.fillna(0, inplace=True)
            for col in pick_all.columns:
                pick_all[col] = map(eval, pick_all[col].astype('string')+'.00')  # eval(1/2.0) = 0.5


        self.ori_exp = '{}/cluster_input.xls'.format(self.work_dir)
        self.sort_exp = ''
        self.pvalue_xls = '{}/enrich_pvalue.xls'.format(self.work_dir)
        pvalue_all.to_csv(self.pvalue_xls, sep='\t')
        pick_all.to_csv(self.ori_exp, sep='\t')
        if self.option("scale") == 'T':
            scale_exp =self.scale_data(self.ori_exp,self.work_dir, method='UV')
        else:
            scale_exp = self.ori_exp

        opts = {
            "exp" : scale_exp,
            "sct" : self.option("sct"),
            "scd" : self.option("scd"),
            "scm" : self.option("scm"),
            "mct" : self.option("mct"),
            "mcd" : self.option("mcd"),
            "mcm" : self.option("mcm"),
            "n_cluster" : self.option("n_cluster")
        }
        self.cluster_tool.set_options(opts)
        self.cluster_tool.on('end', self.set_output)
        self.cluster_tool.run()


    def scale_data(self, ori_table, outDir, method='UV',sample=None):
        table = pd.read_table(ori_table, sep='\t', index_col=0)
        if sample:
            table = table[sample]   ##只用所选的样本做标准化

        all = np.array(table).flatten()
        m = np.mean(all)
        s = np.std(all)

        if method=='Par':
            scaled_data = table.apply(lambda x: (x - m) / np.sqrt(s), axis=1)
        elif method == 'Ctr':
            scaled_data = table.apply(lambda x: (x - m), axis=1)
        else:
            scaled_data = table.apply(lambda x: (x - m) / s, axis=1)
        scaled_data =scaled_data.fillna(0)  # 防止处理后出现空值
        exp_profile = os.path.join(outDir, "scale_data.xls")
        scaled_data.to_csv(exp_profile, index=True, header=True, sep="\t")
        return exp_profile


    def merge_pathway_desc(self, in_file, out_file):
        desc_data = pd.read_table(self.desc_xls, sep='\t',index_col=0)
        data = pd.read_table(in_file, sep='\t', index_col=0)
        data_add_desc = desc_data.join(data)
        data_add_desc.to_csv(out_file, sep='\t')

    def set_output(self):
        self.move_dir(self.cluster_tool.output_dir,self.output_dir)

        self.merge_pathway_desc(self.ori_exp, self.output_dir+'/cluster_exp.xls')

        ## 跟树的顺序对原始表达量表排序
        order_file = self.cluster_tool.work_dir + '/cluster_exp.xls'
        self.order_origin(self.output_dir+'/cluster_exp.xls', order_file,self.output_dir+'/cluster_exp.xls')

        if os.path.exists(self.output_dir+'/enrich_pvalue.xls'):
            os.remove(self.output_dir+'/enrich_pvalue.xls')
        os.link(self.pvalue_xls, self.output_dir+'/enrich_pvalue.xls')
        self.end()

    def run(self):
        super(EnrichClusterModule, self).run()
        self.metablome = Metabolome()
        self.metablome._config = Config()
        if self.option("task_id") in ['']:## 工作流
            self.version_value = '94.2'
        else:
            self.version_value = self.metablome.find_version_from_task(self.option("task_id"))
        self.run_enrich()


    def end(self):
        super(EnrichClusterModule, self).end()

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

    def order_origin(self, origin_file, order_file, outfile):
        order_file = os.path.join(self.work_dir, order_file)
        outfile = os.path.join(self.work_dir, outfile)
        table = pd.read_table(order_file, sep="\t", index_col=0)
        select_names = table.index
        origin_table = pd.read_table(origin_file, sep="\t", index_col=0)
        selecl_origin = origin_table.loc[select_names,]
        selecl_origin.to_csv(outfile, index=True, header=True,sep="\t")