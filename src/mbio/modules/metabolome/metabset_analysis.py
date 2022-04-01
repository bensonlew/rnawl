# -*- coding: utf-8 -*-
# __author__ = zouguanqing
# last_modifiy = modified 201908

import os
import glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import unittest
import pandas as pd

class MetabsetAnalysisModule(Module):
    """
    单个代谢集分析
    """
    
    def __init__(self,work_id):
        super(MetabsetAnalysisModule, self).__init__(work_id)
        options = [
            {"name": "anno_overview", "type": "infile", "format": "metabolome.overview_table"},
            {"name": "metabset", "type": "string"}, #1列，无表头  "format": "metabolome.metabset,sequence.profile_table"
            {"name": "metabset_row", "type": "string"},  # 2列，分别为集名，meta id（,分割） , "format": "metabolome.mul_metabset,sequence.profile_table"
            ##enrich
            {"name": "correct", "type": "string","default": "BH"},
            {"name": "bg", "type": "string", "default": "species"},
             # 背景，project:本项目鉴定到的代谢物合集; species:本物种全部代谢物合集; kegg:KEGG数据库全部代谢物合集
            {"name": "species", "type": "string", "default": "all"},
            {"name": "method", "type": "string", "default": "rbc"}, ##topo方法
            ##keggp
            {"name": "ko_overview","type": "infile", "format":"sequence.profile_table" },
            ##keggc
            {"name": "database", "type":"string","default":"CBR"},
            #hmdb
            {"name": "hmdb_overview","type": "infile", "format":"sequence.profile_table" },
            # vip
            {'name': 'diff_dir', 'type': 'infile', 'format': 'annotation.mg_anno_dir'},
            ## 公用
            {"name": "metab_trans", "type": "infile", "format": "sequence.profile_table"},
            {"name": "metab_table", "type": "infile", "format": "sequence.profile_table"}, # 标准化时用，用原始丰度表
            ##cluster
            # {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            # {"name": "group_detail", "type": "string"},
            {"name": "scale_exp","type": "infile", "format": "sequence.profile_table"},  #!!!
            {"name": "metab_cluster_method", "type": "string", "default": "hierarchy"},
            {"name": "sam_cluster_method", "type": "string", "default": "hierarchy"},
            {"name": "group_method", "type": "string", "default": "average"},
            {"name": "top_meta", "type": "int", "default": 50},
            {"name": "metab_dist", "type": "string", "default": "euclidean"},
            {"name": "metab_cluster", "type": "string", "default": "complete"},
            {"name": "n_cluster", "type": "int", "default": 10},
            {"name": "sam_dist", "type": "string", "default": "euclidean"},
            {"name": "sam_cluster", "type": "string", "default": "complete"},
            {"name": "skip_analysis","type": "string", "default":""},
            {"name": "corr_clust_top", "type": "int", "default": 50},
            {"name": "task_id", "type": "string", "default": ""}, # 用于区分新老工作流
            {"name": "database_version", "type": "string", "default": ""}, # 用于区分新老工作流


        ]
        self.add_option(options)
        self.metab_cluster = self.add_tool("metabolome.metabset.metab_cluster")
        self.corr_cluster = self.add_tool("metabolome.compare.corr_tree")
        self.metab_vip = self.add_tool("metabolome.metabset.metab_vip")
        self.keggc = self.add_tool("metabolome.annotation.anno_keggc")
        self.keggp = self.add_tool("metabolome.metabset.keggp")
        self.enrich_module = self.add_module("metabolome.enrich_topo")
        self.anno = self.add_tool("metabolome.annotation.anno_hmdb")
        self.ipath = self.add_tool("metabolome.metabset.ipath3")
        self.run_list = []

    def check_options(self):

        if not self.option("anno_overview").is_set:
            raise OptionError("请传入anno_overview！")
        if not self.option("ko_overview").is_set:
            raise OptionError("请传入ko_overview！")
        # if not self.option("hmdb_overview").is_set:
        #     raise OptionError("请传入hmdb_overview！")


    def run(self):
        super(MetabsetAnalysisModule, self).run()
        self.skips = self.option('skip_analysis').split(',')
        self.sub_trans_table()
        if 'cluster' not in self.skips:
            self.run_cluster()
        if 'vip' not in self.skips:
            self.run_vip()
        if 'keggc' not in self.skips:
            self.run_keggc()
        if 'keggp' not in self.skips:
            self.run_keggp()
        if 'enrich' not in self.skips:
            self.run_enrich()
        if 'hmdb' not in self.skips:
            self.run_hmdb_anno()
        if 'ipath' not in self.skips:
            self.run_ipath()
        if 'corr_cluster' not in self.skips:
            self.run_corr_cluster()
        self.on_rely(self.run_list,self.set_output)
        for t in self.run_list:
            t.run()

    def sub_trans_table(self):
        trans_table = self.option('metab_table').path
        data = pd.read_table(trans_table,sep='\t',header=0)
        metab_table = pd.read_table(self.option('metabset'),sep='\t',header=-1)
        meta_list = metab_table[0].tolist()
        sub_data = data[data['metab_id'].apply(lambda x: x in meta_list)]
        self.sub_exp_table = self.work_dir+'/sub_trans_table'
        sub_data.to_csv(self.sub_exp_table,sep='\t',index=False)
        ##202009 top corr
        if len(sub_data) <= self.option("corr_clust_top"):
            self.corr_sub_exp_table = self.sub_exp_table
        else:
            self.corr_sub_exp_table = self.work_dir + '/corr_sub_tarans_table'

            sub_data.set_index('metab_id',inplace=True)
            sub_data['TotalTotalTotal'] = sub_data.apply(lambda x: x.sum(),axis=1)
            corr_sub_data = sub_data.sort_values(by='TotalTotalTotal', ascending=False)[0:self.option("corr_clust_top")]
            corr_sub_data.drop('TotalTotalTotal',1,inplace=True)
            corr_sub_data.to_csv(self.corr_sub_exp_table, sep='\t')



    def sub_scal_trans_table(self):
        if self.option("scale_exp").is_set:
            trans_table = self.option('scale_exp').path
            data = pd.read_table(trans_table,sep='\t',header=0)
            metab_table = pd.read_table(self.option('metabset'),sep='\t',header=-1)
            meta_list = metab_table[0].tolist()
            sub_data = data[data['metab_id'].apply(lambda x: x in meta_list)]
            self.sub_scale_exp_table = self.work_dir+'/sub_scale_trans_table'
            sub_data.to_csv(self.sub_scale_exp_table,sep='\t',index=False)


    def run_cluster(self):

        self.logger.info("start run metabset_cluster !")

        options = {
            'exp': self.sub_exp_table,
            'sct': self.option("sam_cluster_method"),
            'mct': self.option("metab_cluster_method"),
            'metab_trans': self.option('metab_trans'),

        }
        if  self.option("metab_cluster_method") == "hierarchy":
            options['mcm'] = self.option("metab_cluster")
            options['mcd'] = self.option("metab_dist")
            options['n_cluster'] = self.option("n_cluster")
        if  self.option("sam_cluster_method") == "hierarchy":
            options['scm'] = self.option("sam_cluster")
            options['scd'] = self.option("sam_dist")
        if  self.option("metab_cluster_method") == "kmeans":
            options['n_cluster'] = self.option("n_cluster")
            options['mcd'] = self.option("metab_dist")
        if self.option("scale_exp").is_set:  #!!!   #标准化转到tool中运行
            self.sub_scal_trans_table()
            options['exp'] = self.sub_scale_exp_table
            options["before_scale"] = self.sub_exp_table  ###

        self.logger.info(options)
        self.metab_cluster.set_options(options)
        self.run_list.append(self.metab_cluster)


    def run_vip(self):
        self.logger.info("start run_vip !")
        opts = {
            'diff_dir': self.option('diff_dir'),
            "metab_set_table": self.option('metabset'),
            'vip_type': 'oplsda',
            'vip_cut': 1.0,
            'vip_top': 30,
            'mct': 'hierarchy',
            'mcd': 'euclidean',
            'mcm': 'complete',
            'scale': True,
            'metab_trans': self.option('metab_trans'),
            'metab_table': self.option('metab_table')
        }

        self.metab_vip.set_options(opts)
        self.run_list.append(self.metab_vip)
        # self.metab_vip.on('end', self.set_db)
        # self.metab_vip.run()

    def run_corr_cluster(self):
        self.logger.info("start corr cluster !")
        opts = {
            'exp':  self.corr_sub_exp_table,
            'scm': "complete",
            'scd': "euclidean",
            'sct': "hierarchy",
            'corr_method': "pearson",
            'file_tran': True,
            'metab_trans' : self.option("metab_trans").path
        }
        self.corr_cluster.set_options(opts)
        self.run_list.append(self.corr_cluster)


    def run_keggc(self):
        overview = pd.read_table(self.option('anno_overview').path,sep='\t')
        metabs = pd.read_table(self.option('metabset'),sep='\t',header=-1)  #1列，无表头
        metabs_list = metabs[0].tolist()
        new_names = {}
        if 'metab' in overview.columns:
            new_names['metab'] = 'Metabolite'
        if 'compound_id' in overview.columns:
            new_names['compound_id'] = 'KEGG Compound ID'
        overview.rename(columns=new_names,inplace=True)
        sub_overview = overview[overview['metab_id'].apply(lambda x: x in metabs_list)]
        self.new_overview = self.work_dir+'/sub_overview.xls'
        sub_overview.to_csv(self.new_overview,sep='\t',index=False)

        options = {
            'metab_table': self.new_overview,
            'database_name':self.option('database'),
            'database_version': self.option("database_version")
        }
        self.keggc.set_options(options)
        self.run_list.append(self.keggc)
        # self.keggc.on('end', self.set_db)
        # self.keggc.run()

    def run_keggp(self):
        options = {
            "anno_overview": self.option("ko_overview"),  #
            "metabset": self.option("metabset_row"),   ##行，2列，无表头，集名+ 元素
            "trans_file": self.option("anno_overview").path,  #总览 ，转id
            'database_version': self.option("database_version")
        }
        self.keggp.set_options(options)
        self.run_list.append(self.keggp)
        # self.keggp.on('end', self.set_db)
        # self.keggp.run()

    def run_enrich(self):
        # super(MetabsetEnrichWorkflow, self).run()

        options = {
            "anno_overview": self.option('anno_overview'),  #总览
            "ko_overview": self.option('ko_overview'), # by ghd @20191015
            "metabset": self.option("metabset"),  #1列，无表头
            "correct": self.option("correct"),
            "bg": self.option("bg"),
            "species":  self.option("species"),
            "method" :   self.option("method"),
            'database_version': self.option("database_version")
        }
        self.enrich_module.set_options(options)
        self.run_list.append(self.enrich_module)
        # self.enrich_module.on('end', self.set_db)
        # self.enrich_module.run()

    def run_hmdb_anno(self):  #hmdb
        options = {
            "anno_overview": self.option('hmdb_overview'), #HmdbLevel_Origin
            "metabset": self.option('metabset'),  ##metabtable #1列，无表头
            "type" : 'metabsethmbd',
        }
        self.anno.set_options(options)
        self.run_list.append(self.anno)
        # self.anno.on('end', self.set_db)
        # self.anno.run()


    def run_ipath(self):
        opts = {
            "metabset": self.option("metabset_row"),  ##行，2列，无表头，集名+ 元素
            "anno_overview": self.option("anno_overview"),   #总览
        }
        self.ipath.set_options(opts)
        self.run_list.append(self.ipath)


    def set_output(self):

        map = {
            'MetabsetCluster':self.metab_cluster,
            'MetabsetVip':self.metab_vip,
            'MetabsetCorr': self.corr_cluster,
            'MetabsetKeggc':self.keggc,
            'MetabsetKeggp': self.keggp,
            'MetabsetEnrich': self.enrich_module,
            'MetabsetHmdb' :  self.anno,
            'MetabsetIpath': self.ipath
        }
        for d in map:
            out_dir = self.output_dir +'/'+d
            if os.path.exists(out_dir):
                shutil.rmtree(out_dir)
            if os.path.exists(map[d].output_dir):
                shutil.copytree(map[d].output_dir,out_dir)
        self.end()


    def end(self):
        super(MetabsetAnalysisModule, self).end()



class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/liulinmeng/metabolome/package/procrustes/module'
        method = "pcoa"
        data = {
            "id": "procrustes" + str(random.randint(1,10000)),
            "type": "module",
            "name": "metabolome.procrustes",
            "instant": False,
            "options": dict(
                metab_table = test_dir + "/metab.select_table.xls",
                asso_table = test_dir + "/corr.select_table.xls",
                #group_file = test_dir + "/group.txt",
                metab_dist = "bray_curtis",
                asso_dist = "bray_curtis",
                method = method,
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()

