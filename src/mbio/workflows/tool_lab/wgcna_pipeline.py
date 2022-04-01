# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
from mainapp.models.mongo.ref_rna_v2 import RefRnaV2
from collections import OrderedDict
from bson.objectid import ObjectId
import pandas as pd
import glob
import datetime
import os
from biocluster.config import Config
import json
import time
from biocluster.file import getsize, exists
from biocluster.file import download
import re
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.ref_rna_v2.wgcna.wgcna import Wgcna
from mbio.packages.ref_rna_v2.chart_advance import ChartAdvance
import unittest


class WgcnaPipelineWorkflow(Workflow):
    """
    666
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(WgcnaPipelineWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'exp_matrix', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'gene_list', 'type':'infile', 'format':'ref_rna_v2.common'},
            {'name': 'trait_path', 'type': 'infile', 'format':'ref_rna_v2.common'},
            {'name': 'trait_type', 'type': 'string', 'default': None},
            {'name': 'me', 'type': 'float', 'default': None},
            {'name': 'cv', 'type': 'float', 'default': None},
            {'name': 'networkType', 'type': 'string'},
            {'name': 'power', 'type': 'string'},
            {'name': 'minModuleSize', 'type': 'int'},
            {'name': 'minKMEtoStay', 'type': 'float'},
            {'name': 'mergeCutHeight', 'type': 'float'},
            {'name': 'corr_method', 'type': 'string'},
            {'name': 'top', 'type': 'string'},
            {'name': 'threshold', 'type': 'string'},
            {'name': 'source', 'type': 'string', 'default': 'tool_lab'},
            {'name': 'relate_id', 'type': 'string'},
            {'name': 'task_id', 'type': 'string'},
            {'name': 'project_task_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.wgcna_prepare = self.add_tool("tool_lab.wgcna.wgcna_prepare")
        self.wgcna_module = self.add_tool("tool_lab.wgcna.wgcna_module")
        self.wgcna_relate = self.add_tool("tool_lab.wgcna.wgcna_relate")
        self.wgcna_network = list()
        self.dump_tool = self.api.api("tool_lab.wgcna")
        self.ref_rna = Wgcna()  # 用于创建主表
        # self.ref_rna.DBVersion =  self.config.DBVersion
        # self.ref_rna.db = Config().get_mongo_client(mtype='ref_rna_v2')[Config().get_mongo_dbname('ref_rna_v2')]
        # Config().get_mongo_client(mtype='ref_rna_v2')[Config().get_mongo_dbname('ref_rna_v2')]

    def check_trait(self, trait_path, exp_path):

        with open(trait_path, 'rb') as f, open(exp_path, 'r') as e:
            exp_sample_list = e.readline().strip().split('\t')[1:]


            def is_number(s):
                try:
                    float(s)
                    return True
                except ValueError:
                    pass

            phenoypes = dict()
            phenoype_list = list()
            # first_sample, test_one = f.readline().strip().split()[0:2]
            for phenoype in f.readline().strip().split()[1:]:
                phenoypes[phenoype] = list()
                phenoype_list.append(phenoype)
            if self.option('trait_type') == "discrete":
                if len(phenoype_list) >= 2:
                    self.set_error('选择非连续类型时，表型数据只可以为1列')
            trait_samples = list()
            for line in f:
                cols = line.strip().split()
                if cols[0] not in exp_sample_list:
                    self.set_error("表型表中的样本不在表达量中，请核实")
                if len(cols) < len(phenoype_list) + 1:
                    self.set_error("样本表型数据缺失", )
                for n, p in enumerate(phenoype_list):
                    phe = cols[n+1]
                    if self.option('trait_type') == "discrete":
                        if is_number(phe):
                            self.set_error('选择非连续类型时，表型数据不应该是数字')
                    else:
                        if not is_number(phe):
                            self.set_error("选择连续类型时，表型数据必须是数字")
                    phenoypes[p].append(phe)
                trait_samples.append(cols[0])
            for p, phe in phenoypes.items():
                if len(set(phe)) == 1:
                    self.set_error('表型只有单一的值请删除该列')
                    return json.dumps({'success': False, 'info': "表型%s只有单一的值请删除该列", "variables":[p], "code" : "C2903414"})

    def run(self):
        if self.option('source') == 'project':
            self.exp_path, self.group_path = self.check_file_path()
            self.exp_file = self.download_s3_file(self.exp_path, 'exp_file.txt')
            self.trait_file = self.download_s3_file(self.group_path, 'group_table.txt')
        if self.option('source') == 'tool_lab':
            self.exp_file = self.option('exp_matrix').prop['path']
            self.trait_file = self.option('trait_path').prop['path']
        self.check_trait(self.trait_file, self.exp_file)
        self.wgcna_prepare.on("end", self.run_wgcna_module)
        self.wgcna_module.on("end", self.run_wgcna_relate)
        self.wgcna_relate.on("end", self.run_wgcna_network)
        if self.option('gene_list').is_set:
            self.run_exp()
        self.run_wgcna_prepare()
        super(WgcnaPipelineWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        workflow_output = self._sheet.output
        # # add result info
        self.add_module_result()
        self.add_relate_result()
        self.add_network_result()
        self.dump_tool.update_db_record('sg_wgcna_pipeline', self.option('main_id'), main_id=ObjectId(self.option('main_id')),
                                        output_dir=workflow_output, status='end')
        self.end()

    def end(self):
        # link file to workflow output_dir
        os.mkdir(self.output_dir + '/' + 'wgcna_prepare')
        os.mkdir(self.output_dir + '/' + 'wgcna_module')
        os.mkdir(self.output_dir + '/' + 'wgcna_relate')
        os.mkdir(self.output_dir + '/' + 'wgcna_network')

        for each in glob.glob(self.wgcna_prepare.output_dir+'/*'):
            base_name = os.path.basename(each)
            target = os.path.join(self.output_dir + '/' + 'wgcna_prepare', base_name)
            os.link(each, target)
        for each in glob.glob(self.wgcna_module.output_dir+'/*'):
            base_name = os.path.basename(each)
            if 'Network-heatmap' in base_name:
                target = os.path.join(self.output_dir, base_name)
            else:
                target = os.path.join(self.output_dir + '/' + 'wgcna_module', base_name)
            os.link(each, target)
        for each in glob.glob(self.wgcna_relate.output_dir+'/*'):
            base_name = os.path.basename(each)
            target = os.path.join(self.output_dir + '/' + 'wgcna_relate', base_name)
            os.link(each, target)
        # link network result
        for each_tool in self.wgcna_network:
            for each in glob.glob(each_tool.output_dir+'/*'):
                base_name = each_tool.option("module") + '.' + os.path.basename(each)
                # base_name = os.path.basename(each)
                target = os.path.join(self.output_dir, 'wgcna_network', base_name)
                os.link(each, target)
        # upload result
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["06 Advanced_Analysis", "", "高级分析结果目录",0],
            ["06 Advanced_Analysis/01 WGCNA", "", "WGCNA分析", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "WGCNA一键化分析结果目录",0,"211595"],
            ["./*Network-heatmap*", "", "TOM矩阵热图", 0,],
            ["./wgcna_prepare", "", "WGCNA数据预处理文件",0,"211596"],
            ["./wgcna_prepare/powerEstimate_*", "", "soft power 值(power为文件名中的数字)",0,"211597"],
            ["./wgcna_prepare/ignored_gene.list", "", "没有用于聚类的基因列表",0,"211598"],
            ["./wgcna_prepare/sampleClustering.pdf", "", "样本聚类", 0],
            ["./wgcna_prepare/pick_power.pdf", "", "无尺度容适曲线及无尺度平均连通度曲线", 0],
            # ["./wgcna_prepare/sample.cluster.*.txt", "", "聚类树形状和分支长度等信息",0,"211599"],
            ["./wgcna_prepare/scale_free_analysis.xls", "", "无尺度分析结果",0,"211600"],
            ["./wgcna_prepare/exp_matrix_after_filtering.txt", "", "过滤后的表达量表",0,"211601"],
            ["./wgcna_prepare", "", "WGCNA预处理分析结果目录",0,"211602"],
            ["./wgcna_relate", "", "WGCNA模块与表型相关性文件",0,"211603"],
            ["./wgcna_relate/block_*_gene_*pdf", "", "基因与表型相关性热图pdf",0,"211604"],
            ["./wgcna_relate/block_*_gene_*png", "", "基因与表型相关性热图png",0,"211605"],
            ["./wgcna_relate/module_trait.correlation.xls", "", "模块与表型相关性系数表",0,"211606"],
            ["./wgcna_relate/module_trait.correlation_pvalues.xls", "", "模块与表型相关显著性统计表",0,"211607"],
            ["./wgcna_relate/gene_trait.correlation.xls", "", "基因与表型相关性系数表",0,"211608"],
            ["./wgcna_module", "", "WGCNA模块识别文件",0,"211609"],
            ["./wgcna_module/block_1_dendrogram.pdf", "", "模块分类树pdf图",0,"211610"],
            ["./wgcna_module/block_1_dendrogram.png", "", "模块分类树png图",0,"211611"],
            ["./wgcna_module/seq_id2gene_name.txt", "", "模块成员基因列表",0,"211612"],
            ["./wgcna_module/module_size.stat.xls", "", "模块成员统计表",0,"211613"],
            ["./wgcna_module/module_corr.matrix.xls", "", "模块相关性关系表",0,"211614"],
            ["./wgcna_module/membership.xls", "", "模块成员聚类详情表",0,"211615"],
            ["./wgcna_module/eigengenes.txt", "", "各模块特征基因列表",0,"211616"],
            ["./wgcna_module/*RData", "", "模块分析R数据, 可用R打开查看",0,"211617"],
            ["./wgcna_network", "", "WGCNA可视化分析文件",0,"211618"],
            ["./wgcna_network/*.json", "", "该模块网络图数据, 文本格式",0,"211619"],
            ["./wgcna_network/*network.nodes.txt", "", "该模块网络图节点列表",0,"211620"],
            ["./wgcna_network/*network.edges.txt", "", "该模块网络图边列表",0,"211621"],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],

        ])
        super(WgcnaPipelineWorkflow, self).end()

#   run tools
    def run_exp(self):
        self.exp_genelist = os.path.join(self.work_dir, 'exp_genelist.txt')
        gene_list = list()
        with open(self.option('gene_list').path, 'r') as g:
            for line in g.readlines():
                gene_list.append(line.strip())

        exp = pd.read_table(self.exp_file, header=0, index_col=0, sep='\t')
        exp_new = exp.loc[gene_list]
        exp_new.to_csv(self.exp_genelist, header=True, index=True, sep='\t')

    def run_wgcna_prepare(self):
        if self.option('gene_list').is_set:
            options = dict(
                exp=self.exp_genelist,
                me=self.option('me'),
                cv=self.option('cv'),
            )
        else:
            options = dict(
                exp=self.exp_file,
                me=self.option('me'),
                cv=self.option('cv'),
            )
        self.wgcna_prepare.set_options(options)
        self.wgcna_prepare.run()

    def run_wgcna_module(self):
        exp_matrix = os.path.join(self.wgcna_prepare.output_dir, "exp_matrix_after_filtering.txt")
        if self.option('power') == '自确定' or float(self.option("power")) <= 0:
            self.power = glob.glob(self.wgcna_prepare.output_dir + '/powerEstimate_*')[0].split("_")[-1]
            if self.power == "NA":
                self.power = 6
        else:
            self.power = self.option("power")
        options = dict(
            datExpr=exp_matrix,
            mergeCutHeight=float(self.option('mergeCutHeight')),
            power=int(self.power),
            networkType=self.option('networkType'),
            minModuleSize=int(self.option('minModuleSize')),
            minKMEtoStay=float(self.option('minKMEtoStay')),
        )
        self.wgcna_module.set_options(options)
        self.wgcna_module.run()

    def run_wgcna_relate(self):
        exp_matrix = os.path.join(self.wgcna_prepare.output_dir, "exp_matrix_after_filtering.txt")
        eigengenes = os.path.join(self.wgcna_module.output_dir, "eigengenes.txt")
        module_Rdata = os.path.join(self.wgcna_module.output_dir, "blockwiseModules_result.RData")
        options = dict(
            datExpr=exp_matrix,
            MEs=eigengenes,
            traits=self.trait_file,
            corType=self.option('corr_method'),
            block_Rdata=module_Rdata
        )
        self.wgcna_relate.set_options(options)
        self.wgcna_relate.run()

    def run_wgcna_network(self):
        trait_corr_f = glob.glob(self.wgcna_relate.output_dir+'/module_trait.correlation.xls')[0]
        corr_pvalue_f = glob.glob(self.wgcna_relate.output_dir+'/module_trait.correlation_pvalues.xls')[0]
        trait_corr = pd.read_table(trait_corr_f, index_col=0, header=0)
        corr_pvalue = pd.read_table(corr_pvalue_f, index_col=0, header=0)
        if "MEgrey" in trait_corr.index:
            trait_corr = trait_corr.loc[trait_corr.index.drop("MEgrey"), :]
            corr_pvalue = corr_pvalue.loc[corr_pvalue.index.drop("MEgrey"), :]
        exp_matrix = os.path.join(self.wgcna_prepare.output_dir,"exp_matrix_after_filtering.txt")
        self.gid2gname = self.export_wgcna_g2n_matrix_new(exp_matrix)
        # self.export_wgcna_g2n_matrix(exp_matrix)
        target_modules = set()
        traits = trait_corr.columns
        # for trait in traits:
        #     p_small_modules = corr_pvalue[corr_pvalue[trait] <= 0.1].index
        #     if not p_small_modules.shape[0]:
        #         p_small_modules = corr_pvalue[corr_pvalue[trait] <= 1].index
        #     select_module = trait_corr.loc[p_small_modules, [trait]].abs().idxmax()[0]
        for trait in traits:
            p_small_modules = corr_pvalue[corr_pvalue[trait] <= 1].index
            # select_module = trait_corr.loc[p_small_modules, [trait]].abs().idxmax()[0]
            for i in range(len(p_small_modules)):
                select_module=p_small_modules[i].split('ME')[1]
                target_modules.add(select_module)
        for module in target_modules:
            tool = self.add_tool("tool_lab.wgcna.wgcna_network")
            options = dict(
                module=module,
                threshold=self.option('threshold'),
                top=self.option("top"),
                step3output=self.wgcna_relate.output_dir,
                step2output=self.wgcna_module.output_dir,
            )
            tool.set_options(options)
            self.wgcna_network.append(tool)
        self.on_rely(self.wgcna_network, self.set_db, 'wgcnanetwork')
        for tool in self.wgcna_network:
            tool.run()

#   dump data to db for each step
    def add_prepare_result(self):
        # add result info
        self.dump_tool.add_prepare_detail(self.wgcna_prepare.output_dir, main_id=self.option('main_id'),)

    def add_module_result(self):
        gene2name = self.gid2gname
        exp_matrix = self.wgcna_prepare.output_dir + '/exp_matrix_after_filtering.txt'
        self.dump_tool.add_module_detail(self.wgcna_module.work_dir, gene2name, exp_matrix, self.option('main_id'))

    def add_relate_result(self):
        # add result info
        seq_annot = self.wgcna_module.work_dir + '/gene_module_detail.xls'
        self.dump_tool.add_relate_detail(self.wgcna_relate.work_dir, seq_annot, self.option('main_id'))

    def add_network_result(self):
        module_list = list()
        for each_tool in self.wgcna_network:
            # add network detail
            self.dump_tool.add_network_detail(each_tool.output_dir, self.option('main_id'), each_tool.option('module'))
            module_list.append(each_tool.option('module'))
        self.dump_tool.update_db_record('sg_wgcna_pipeline', self.option('main_id'), network_module=module_list)

    def export_wgcna_g2n_matrix_new(self, exp_matrix):
        target_seqs = pd.read_table(exp_matrix, header=0, index_col=0).index.tolist()
        gene2name = pd.DataFrame({'gene_id': target_seqs, 'gene_name': None})
        output = os.path.join(self.wgcna_module.output_dir, "seq_id2gene_name.txt")
        gene2name.to_csv(output, sep='\t', header=True, index=False)
        return output

    def check_file_path(self):
        collection_name = 'tool_thurl'
        project_type = 'tool_lab'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        conn_upset = db[collection_name]
        status = 'start'
        count_time = 0
        while status == 'start':
            if count_time > 600:
                self.set_error('超过十分钟还没有结果文件生成，请检查是否生成文件时报错')
                break
            time.sleep(10)
            print 'sleep 10s'
            try:
                upset = conn_upset.find_one(
                    {'task_id': self.option('project_task_id'), 'relate_id': ObjectId(self.option('relate_id'))})
                status = upset['status']
            except:
                pass
            count_time += 10
        upset = conn_upset.find_one(
            {'task_id': self.option('project_task_id'), 'relate_id': ObjectId(self.option('relate_id'))})
        exp_path = upset['exp_path']
        group_path = upset['group_path']
        return exp_path, group_path

    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        if os.path.exists(path):
            to_path = path
        elif exists(path):
            download(path, to_path)
        else:
            self.set_error('file can not find %s', variables=(path,), code='13700502')
        return to_path


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.wgcna_pipeline import WgcnaPipelineWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'wgcna_pipeline_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'tool_lab.wgcna_pipeline',
            'options': {
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/WGCNA/exp_matrix',
                'trait_path': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/WGCNA/all_66_3.txt.checked.txt',
                "networkType": "signed",
                "threshold": "0.02",
                "trait_type": "continuous",
                "cv": "0.1",
                "minModuleSize": "30",
                "top": "30",
                "power": "14",
                "minKMEtoStay": "0.3",
                "me": "0.5",
                "corr_method": "pearson",
                "mergeCutHeight": "0.25",
            }
        }
        wsheet = Sheet(data=data)
        wf =WgcnaPipelineWorkflow(wsheet)
        wf.sheet.id = 'batch_effect'
        wf.sheet.project_sn = 'batch_effect'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
