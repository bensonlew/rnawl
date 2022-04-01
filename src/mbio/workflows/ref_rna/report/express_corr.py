# -*- coding: utf-8 -*-
# __author__ = 'konghualei, 20170421'
import web
import json
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import group_detail_sort
from bson import ObjectId
from biocluster.config import Config
import datetime
import pandas as pd
import shutil
import re
import os

from biocluster.workflow import Workflow

class ExpressCorrWorkflow(Workflow):
    """
    计算表达量相关性
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(ExpressCorrWorkflow, self).__init__(wsheet_object)
        options = [
            {"name":"express_file","type":"string"},
            {"name":"group_id","type":"string"},
            {"name":"group_detail","type":"string"},
            {"name":"correlation_id","type":"string"},
            {"name":"update_info","type":"string"},
            {"name":"type","type":"string","default":"gene"}, #传给to_file 参数
            {"name":"method","type":"string","default":"pearson"}, #聚类方法
            # {"name":"hclust_method",'type':"string",'default':'complete'}, #层次聚类方式
            {"name":"express_level","type":"string"}, #对应fpkm/tpm
            {"name":"corr_pca","type":"string"} #pca 或 correlation分析
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        if self.option("corr_pca") == "corr":
            self.corr = self.add_tool('statistical.correlation')
        if self.option("corr_pca") == "pca":
            self.pca = self.add_tool('rna.pca')
        with open(self.option('group_id'), 'r+') as f1:
            f1.readline()
            if not f1.readline():
                self.group_id = 'all'
            else:
                self.group_id = self.option('group_id')
        self.logger.info("开始打印group_id")
        self.logger.info(self.group_id)

    def run_corr(self):
        """样本间相关性分析"""
        if self.group_id in ['all','All','ALL']:
            new_fpkm = self.option("express_file").split(",")[0]
        else:
            specimen = self.get_samples()
            new_fpkm = self.fpkm(specimen)
        opts = {
            "fpkm":new_fpkm,
            "method":self.option('method')
        }
        self.corr.set_options(opts)
        self.corr.on('end', self.set_db)
        self.corr.run()
    
    def get_samples(self):
        edger_group_path = self.option("group_id")
        self.logger.info(edger_group_path)
        samples=[]
        with open(edger_group_path,'r+') as f1:
            f1.readline()
            for lines in f1:
                line=lines.strip().split("\t")
                samples.append(line[0])
        return samples
    
    def run_pca(self):
        """样本间pca分析"""
        if self.group_id in ['all','All','ALL']:
            new_fpkm = self.option("express_file").split(",")[0]
        else:
            new_specimen = self.get_samples()
            print 'haha1'
            print new_specimen
            new_fpkm = self.fpkm(new_specimen)
        self.logger.info(new_fpkm)
        opts = {
            "otutable": new_fpkm
        }
        self.pca.set_options(opts)
        self.pca.on('end', self.set_db)
        self.pca.run()
        
    def fpkm(self,samples):
        fpkm_path = self.option("express_file").split(",")[0]
        fpkm = pd.read_table(fpkm_path, sep="\t")
        print "heihei1"
        print fpkm.columns[1:]
        print 'heihei2'
        print samples
        no_samp = []
        sample_total = fpkm.columns[1:]

        for sam in sample_total:
            if sam not in samples:
                no_samp.append(sam)
        self.logger.info("heihei3")
        self.logger.info(no_samp)
        if no_samp:
            new_fpkm = fpkm.drop(no_samp, axis=1)
            print new_fpkm.columns
            if self.option("corr_pca") == 'corr':
                self.new_fpkm = self.corr.work_dir + "/fpkm"
            if self.option("corr_pca") == "pca":
                self.new_fpkm = self.pca.work_dir + "/fpkm"
            header=['']
            header.extend(samples)
            new_fpkm.columns=header
            new_fpkm.to_csv(self.new_fpkm, sep="\t",index=False)
            return self.new_fpkm
        else:
            return fpkm_path
    
    def set_db(self):
        api_corr = self.api.refrna_corr_express
        if self.option('corr_pca') == 'corr':
            self.logger.info(self.corr.output_dir)
            # _id = api_corr.add_correlation_table(self.corr.output_dir,express_id="58ef0bcba4e1af740ec4c14c",\
                # detail=False, seq_type="gene")
            api_corr.add_correlation_detail(self.corr.output_dir, self.option("correlation_id"),updata_tree=True)

        if self.option('corr_pca') == 'pca':
            self.logger.info(self.pca.output_dir)
            pca_path = self.pca.output_dir
            inserted_id = self.option("correlation_id")
            pca_file = os.path.join(pca_path, 'pca_importance.xls')
            pca_rotation = os.path.join(pca_path, 'pca_rotation.xls')
            site_file = os.path.join(pca_path, 'pca_sites.xls')
            #api_corr.add_pca(pca_file=pca_file, correlation_id=inserted_id)
            #api_corr.add_pca_rotation(input_file=pca_rotation, db_name='sg_express_pca_rotation',
            #                         correlation_id=inserted_id)
            api_corr.add_pca_rotation(input_file=site_file, db_name='sg_express_pca_rotation', correlation_id=inserted_id)
            self._update_pca(pca_file, self.option("correlation_id"))
        self.end()

    def _update_pca(self, pca_file, pca_id):
        #db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        db =  Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
        con = db['sg_express_pca']
        pc_num = []
        pc = {}
        if not os.path.exists(pca_file):
            raise Exception("pca分析错误，没有生成{}文件".format(self.pca.output_dir + "/pca_importance.xls"))
        with open(pca_file, 'r+') as f1:
            f1.readline()
            for lines in f1:
                line = lines.strip().split("\t")
                pc_num.append(line[0])
                pc[line[0]] = line[1]
        con.update({"_id": ObjectId(pca_id)}, {"$set": {"pc_num": pc_num}})
        self.logger.info("更新pc_num成功")
        for keys, values in pc.items():
            con.update({"_id": ObjectId(pca_id)}, {"$set": {keys: values}})
        self.logger.info("更新pc成功！")

    def run(self):
        if self.option("corr_pca") == "corr":
            self.run_corr()
        if self.option("corr_pca") == "pca":
            self.run_pca()
        super(ExpressCorrWorkflow, self).run()

    def end(self):
        if self.option('corr_pca') == 'corr':
            output1_dir = self.corr.output_dir
            result = self.add_upload_dir(output1_dir)
            result.add_relpath_rules([
                [".", "", "表达量相关性分析结果文件"],
            ])
        if self.option('corr_pca') == 'pca':
            output2_dir = self.pca.output_dir
            result2 = self.add_upload_dir(output2_dir)
            result2.add_relpath_rules([[".", "", "表达量PCA分析结果文件"], ])
        super(ExpressCorrWorkflow, self).end()



