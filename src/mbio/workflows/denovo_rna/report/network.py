# -*- coding: utf-8 -*-
# __author__: konghualei
# last modify: 20161208

"""无参转录组网络共表达分析"""

from biocluster.workflow import Workflow
from biocluster.config import Config
import os
import re
from bson.objectid import ObjectId
from cStringIO import StringIO
import json
from bson.son import SON
from bson.objectid import ObjectId
import bson.binary
import datetime
import types
from gridfs import GridFS
class NetworkWorkflow(Workflow):
    """
    报告中调用网络共表达分析时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(NetworkWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "express_file", "type": "string", "default": "none"},
            {"name": "update_info", "type": "string"},
            {"name": "softpower", "type": "int", "default": 9},
            {"name": "dissimilarity", "type": "float", "default": 0.25},
            {"name": "module", "type": "float", "default": 0.1},
            #{"name": "network", "type": "float", "default": 0.2},
            {"name": "network_express_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.network = self.add_tool("denovo_rna.express.network")
        self.output_dir = self.network.output_dir

    def run_network(self):
        exp_files = self.option('express_file').split(',')[0]
        """生成gene_list文件"""
        fpkm_path = os.path.split(exp_files)[0]
        gene_list_path = get_gene_list(exp_files, fpkm_path)
        options = {
            'diff_fpkm': exp_files,
            'softpower': self.option('softpower'),
            'module': self.option('module'),
            #'network': self.option('network'),
            'gene_file': gene_list_path
        }
        self.network.set_options(options)
        self.network.on('end', self.set_db)
        self.network.run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        relpath = [
           ['.', '', '结果输出目录'],
           ['all_nodes.txt', 'txt', 'nodes结果信息'],
           ['all_edges.txt', 'txt', 'edges结果信息'],
           ["ModuleTree.pdf", "pdf", "ModuleTree pdf图"],
           ["softPower.pdf", "pdf", "softpower pdf图"],
           ["ModuleTree.png", "png", "ModuleTree png图"],
           ["softPower.png", "png", "softpower png图"],
           ['removeSample.xls', 'xls', 'removeSample.xls'],
           ['networkHeatmap.pdf', "pdf", "networkHeatmap.pdf"],
           ['sampleClustering.pdf', 'pdf', 'sampleClustering.pdf'],
           ['eigenGeneHeatmap.pdf', 'pdf', 'eigenGeneHeatmap.pdf'],
           ['eigengeneClustering.pdf', 'pdf', 'eigengeneClustering.pdf'],
           ['netcolor2gene.xls', 'xls', 'netcolor2gene.xls'],
           ['removeGene.xls', 'xls', 'removeGene.xls']
        ]
        result_dir.add_regexp_rules([
            [r"CytoscapeInput.*", 'txt', 'Cytoscape作图数据'],
            #[r'*.png', 'png', 'module和softpower图']
        ])
        result_dir.add_relpath_rules(relpath)
        super(NetworkWorkflow, self).end()

    def set_db(self):
        """保存结果表保存到mongo数据库中"""
        api_network = self.api.denovo_network
        network_files = os.listdir(self.output_dir)
        nodes_edges=list()
        all_color = list()
        for f in network_files:
            if re.search(r'CytoscapeInput-edges*', f):
                module_color = re.search(r'\w+-\w+-(\w+).\w+',f).group(1)
                if module_color:
                    all_color.append(module_color)
                else:
                   raise Exception('没有获取到module颜色信息!')
                api_network.add_network_module(network_id = self.option('network_express_id'),
                    module_path = self.output_dir + '/' + f, module_color = module_color)
            if re.search(r'all_nodes.txt', f):
                nodes_path =  self.output_dir + '/' + f
            if re.search(r'all_edges.txt', f):
                edges_path = self.output_dir + '/' + f
            if re.search(r'softPower.png', f):
                softpower_path = self.output_dir + '/' + f
            if re.search(r'ModuleTree.png', f):
                module_path = self.output_dir + '/' + f
        api_network.add_network_detail(network_id = self.option('network_express_id'),\
                node_path = nodes_path, edge_path = edges_path)
        self.update_network(table_id = self.option('network_express_id'), module_path = module_path,
            softpower_path = softpower_path, module= all_color)
        self.end()

    def update_network(self, table_id, module_path, softpower_path, module):
        db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
        #client = Config().mongo_client
        #db_name = Config().MONGODB + '_rna'
        print db_name
        collection = db['sg_denovo_network']
        self.logger.info(module_path)
        self.logger.info(softpower_path)
        self.logger.info(module)
        #with open(softpower_path, "rb") as s, open(module_path, 'rb') as m:
            #softpower_id = StringIO(s.read())
            #softpower_id = bson.binary.Binary(softpower_id.getvalue())
            #module_id = StringIO(m.read())
            #module_id = bson.binary.Binary(module_id.getvalue())
        fs_module=GridFS(client[db_name])
        module_id=fs_module.put(open(module_path,'r'))
        fs_softpower=GridFS(client[db_name])
        softpower_id=fs_softpower.put(open(softpower_path,'r'))
        self.logger.info(softpower_id)
        self.logger.info(module_id)
        self.logger.info(module)
        if not module:
           raise Exception('module 信息不存在！')
        collection.update({'_id': ObjectId(table_id)}, {'$set': {'module_png': ObjectId(module_id), 'softpower_png': ObjectId(softpower_id), 'module': module}})

    def run(self):
        self.run_network()
        print "network run"
        super(NetworkWorkflow, self).run()

def get_gene_list(diff_fpkm, fpkm_path):
    gene_list_path = os.path.join(fpkm_path, 'gene_list.txt')
    gene_list = open(gene_list_path, 'w+')
    i=0
    #gene_list.write("tracking_id" + "\t" + "tracking_id" + "\n")
    with open(diff_fpkm, "r+") as files:
        for f in files:
            i += 1
            if i == 1:
                next
            else:
                file = f.strip().split("\t")
                gene_name = file[0]
                gene_list.write(gene_name+"\n")
    gene_list.close()
    return gene_list_path
