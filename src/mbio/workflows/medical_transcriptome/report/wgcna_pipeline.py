# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
# from mainapp.models.mongo.medical_transcriptome import MedicalTranscriptome
from collections import OrderedDict
from bson.objectid import ObjectId
import pandas as pd
import glob
import datetime
import os
import json
import re
from mbio.packages.lnc_rna.copy_file import CopyFile
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.medical_transcriptome.wgcna.wgcna import Wgcna
from mbio.packages.ref_rna_v2.chart_advance import ChartAdvance


class WgcnaPipelineWorkflow(Workflow):
    """
    666
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(WgcnaPipelineWorkflow, self).__init__(wsheet_object)
        options = list()
        for each in self._sheet.options():
            options.append(dict(name=each, type="string"))
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.wgcna_prepare = self.add_tool("ref_rna_v2.wgcna.wgcna_prepare")
        self.wgcna_module = self.add_tool("ref_rna_v2.wgcna.wgcna_module")
        self.wgcna_relate = self.add_tool("ref_rna_v2.wgcna.wgcna_relate")
        # self.wgcna_network = self.add_tool("rna.wgcna.wgcna_network")
        self.wgcna_network = list()
        self.dump_tool = self.api.api("medical_transcriptome.wgcna")
        self.medical_transcriptome = Wgcna()  # 用于创建主表  #add by fwy 20210108
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/05 Advanced_Analysis/01 WGCNA')
        self.inter_dirs = []

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(WgcnaPipelineWorkflow, self).send_log(data)

    def run(self):
        self.wgcna_prepare.on("end", self.run_wgcna_module)
        self.wgcna_module.on("end", self.run_wgcna_relate)
        self.wgcna_relate.on("end", self.run_wgcna_network)
        self.get_run_log()
        self.run_wgcna_prepare()
        super(WgcnaPipelineWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("medical_transcriptome", table="sg_wgcna_pipeline", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        # save workflow output path--建议标配
        workflow_output = self.get_workflow_output_dir()
        pipeline_params = self.medical_transcriptome.get_main_info(self.option("main_id"),'sg_wgcna_pipeline', self.option("task_id"))['params']
        params_dict = json.loads(pipeline_params)
        params_dict['power'] = self.power
        new_params = json.dumps(params_dict,sort_keys=True,separators=(',',':'))
        self.params = params_dict
        self.dump_tool.update_db_record('sg_wgcna_pipeline', self.option('main_id'), main_id=ObjectId(self.option('main_id')),
                                        output_dir=workflow_output, status='end', params=new_params)
        # add result info
        wgcna_prepare_id = self.add_prepare_result()
        wgcna_module_id = self.add_module_result(wgcna_prepare_id)
        wgcna_relate_id = self.add_relate_result(wgcna_module_id)
        self.add_network_result(wgcna_module_id, wgcna_relate_id)
        self.end()

    def end(self):
        # link file to workflow output_dir
        os.mkdir(self.output_dir + '/' + 'prepare')
        os.mkdir(self.output_dir + '/' + 'module')
        os.mkdir(self.output_dir + '/' + 'relate')
        os.mkdir(self.output_dir + '/' + 'network')
        self.chart()
        for each in glob.glob(self.wgcna_prepare.output_dir+'/*'):
            base_name = os.path.basename(each)
            target = os.path.join(self.output_dir + '/' + 'prepare', base_name)
            os.link(each, target)
        for each in glob.glob(self.wgcna_module.output_dir+'/*'):
            base_name = os.path.basename(each)
            target = os.path.join(self.output_dir + '/' + 'module', base_name)
            os.link(each, target)
        for each in glob.glob(self.wgcna_relate.output_dir+'/*'):
            base_name = os.path.basename(each)
            target = os.path.join(self.output_dir + '/' + 'relate', base_name)
            os.link(each, target)
        # link network result
        for each_tool in self.wgcna_network:
            for each in glob.glob(each_tool.output_dir+'/*'):
                base_name = each_tool.option("module") + '.' + os.path.basename(each)
                target = os.path.join(self.output_dir, 'network', base_name)
                os.link(each, target)
        # upload result

        self.inter_dirs = [
            ["05 Advanced_Analysis", "", "高级分析结果目录",0],
            ["05 Advanced_Analysis/01 WGCNA", "", "WGCNA分析", 0]
        ]

        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "WGCNA一键化分析结果目录",0,"211595"],
            ["./prepare", "", "WGCNA预处理分析结果目录",0,"211596"],
            ["./prepare/powerEstimate_*", "", "soft power 值",0,"211597"],
            ["./prepare/ignored_gene.list", "", "没有用于聚类的基因列表",0,"211598"],
            ["./prepare/sample.cluster.*.txt", "", "聚类树形状和分支长度等信息",0,"211599"],
            ["./prepare/scale_free_analysis.xls", "", "无尺度分析结果",0,"211600"],
            ["./prepare/exp_matrix_after_filtering.txt", "", "过滤后的表达量表",0,"211601"],
            ["./prepare", "", "WGCNA预处理分析结果目录",0,"211602"],
            ["./relate", "", "WGCNA模块与表型分析结果目录",0,"211603"],
            ["./relate/block_*_gene_*pdf", "", "基因与表型相关性热图pdf",0,"211604"],
            ["./relate/block_*_gene_*png", "", "基因与表型相关性热图png",0,"211605"],
            ["./relate/module_trait.correlation.xls", "", "模块与表型相关性系数表",0,"211606"],
            ["./relate/module_trait.correlation_pvalues.xls", "", "模块与表型相关显著性统计表",0,"211607"],
            ["./relate/gene_trait.correlation.xls", "", "基因与表型相关性系数表",0,"211608"],
            ["./module", "", "WGCNA模块识别分析结果目录",0,"211609"],
            ["./module/block_*_dendrogram.pdf", "", "模块分类树pdf图",0,"211610"],
            ["./module/block_*_dendrogram.png", "", "模块成员详情表", "模块分类树png图",0,"211611"],
            ["./module/seq_id2gene_name.txt", "", "模块成员基因列表",0,"211612"],
            ["./module/module_size.stat.xls", "", "模块成员统计表",0,"211613"],
            ["./module/module_corr.matrix.xls", "", "模块相关性关系表",0,"211614"],
            ["./module/membership.xls", "", "模块成员聚类详情表",0,"211615"],
            ["./module/eigengenes.txt", "", "各模块特征基因列表",0,"211616"],
            ["./module/*RData", "", "模块分析R数据, 可用R打开查看",0,"211617"],
            ["./network", "", "WGCNA可视化分析结果目录",0,"211618"],
            ["./network/*.json", "", "该模块网络图数据, 文本格式",0,"211619"],
            ["./network/*network.nodes.txt", "", "该模块网络图节点列表",0,"211620"],
            ["./network/*network.edges.txt", "", "该模块网络图边列表",0,"211621"],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(WgcnaPipelineWorkflow, self).end()


#   run tools
    def run_wgcna_prepare(self):
        options = dict(
            exp=self.option('exp_matrix'),
            me=self.option('me'),
            cv=self.option('cv'),
        )
        self.wgcna_prepare.set_options(options)
        self.wgcna_prepare.run()

    def run_wgcna_module(self):
        exp_matrix = os.path.join(self.wgcna_prepare.output_dir, "exp_matrix_after_filtering.txt")
        if float(self.option("power")) <= 0:
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
            traits=self.option("trait_path"),
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
        # self.gid2gname = self.export_wgcna_g2n_matrix(exp_matrix)
        if os.path.exists(os.path.join(self.wgcna_module.output_dir, "seq_id2gene_name.txt")):
            os.remove(os.path.join(self.wgcna_module.output_dir, "seq_id2gene_name.txt"))

        os.link(self.option("gene_name"), os.path.join(self.wgcna_module.output_dir, "seq_id2gene_name.txt"))
        self.gid2gname = os.path.join(self.wgcna_module.output_dir, "seq_id2gene_name.txt")
        # self.export_wgcna_g2n_matrix(exp_matrix)
        target_modules = set()
        traits = trait_corr.columns
        for trait in traits:
            p_small_modules = corr_pvalue[corr_pvalue[trait] <= 0.1].index
            if not p_small_modules.shape[0]:
                p_small_modules = corr_pvalue[corr_pvalue[trait] <= 1].index
            select_module = trait_corr.loc[p_small_modules, [trait]].abs().idxmax()[0]
            target_modules.add(select_module[2:])
        for module in target_modules:
            tool = self.add_tool("ref_rna_v2.wgcna.wgcna_network")
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
        # create main table
        exp_info = self.medical_transcriptome.get_main_info(self.option("exp_id"), 'sg_exp', self.option("task_id"))
        project_sn = exp_info["project_sn"]
        group_dict = json.loads(self.option('raw_group_dict'), object_pairs_hook=OrderedDict)
        # create main table record
        try:
            exp_info = json.loads(exp_info['params'])
            level,exp_type,quant_method = self.option("level")[0].upper(),exp_info.get('exp_type', "TPM"),exp_info.get('method', "RSEM")
        except:
            level = self.option("level")[0].upper()
            exp_type = 'TPM'
            quant_method = 'RSEM'
        params = dict(
            task_id=self.option("task_id"),
            submit_location="wgcnaprepare",
            task_type=2,
            group_id=self.option("group_id"),
            level=self.option("level"),
            # category=self.option("category"),
            group_dict=group_dict,
            me=self.option("me"),
            cv=self.option("cv"),
            geneset_id=self.option("geneset_id"),
            exp_id=self.option("exp_id")
        )
        params = json.dumps(params,sort_keys=True,separators=(',',':'))
        name = "PipelinePrepare" + '_' + level + '_' + quant_method + '_' + exp_type.upper() + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=self.option("task_id"),
            name=name,
            version="v1",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            exp_id=self.option("exp_id"),
            desc='wgcna pre-processing analysis main table',
            level=self.option("level"),
            # category=self.option("category"),
            params=params,
            pipeline_id=ObjectId(self.option("main_id")),
            status="start")
        main_id = self.medical_transcriptome.insert_main_table('sg_wgcna_prepare', main_info)
        # add result info
        self.dump_tool.add_prepare_detail(self.wgcna_prepare.output_dir, main_id=str(main_id),)
        # update main table
        with open(self.work_dir + '/group_info.txt') as f:
            _ = f.readline()
            group_info = dict()
            for line in f:
                s, g = line.strip().split('\t')
                group_info[s] = g
        workflow_output = self.get_workflow_output_dir()
        exp_matrix = workflow_output + '/prepare/exp_matrix_after_filtering.txt'
        self.dump_tool.update_db_record('sg_wgcna_prepare', main_id, group_info=group_info, exp_matrix=exp_matrix)
        return main_id

    def chart(self):
        chart = ChartAdvance()
        chart.work_dir = self.work_dir + "/"
        sample_tree = self.wgcna_prepare.work_dir + "/sample.cluster.dendrogram.txt"
        group_dict = json.loads(self.option("group_dict"))
        chart.chart_wgcna_sample_tree(sample_tree, group_dict)
        module_tree = self.wgcna_module.work_dir + "/unmerged.module_corr.dendrogram.txt"
        chart.chart_wgcna_module_tree(module_tree)

        module_corr = self.wgcna_module.work_dir + "/module_corr.matrix.xls"
        module_corr_tree = self.wgcna_module.work_dir + "/module_corr.tree.txt"
        chart.chart_wgcna_module_corr(module_corr, module_corr_tree)
        module_stat = self.wgcna_module.work_dir + "/module_size.stat.xls"
        chart.chart_wgcna_module_column(module_stat)
        chart.to_pdf()

        relation_corr = self.wgcna_relate.work_dir + "/module_trait.correlation.xls"
        relation_corr_pvalue = self.wgcna_relate.work_dir + "/module_trait.correlation_pvalues.xls"

        chart.chart_wgcna_relation_corr(relation_corr, relation_corr_pvalue, module_stat)

        gene_trait_corr = self.wgcna_relate.work_dir + "/gene_trait.correlation.xls"
        seq_annot = self.wgcna_module.work_dir + "/gene_module_detail.xls"

        chart.chart_wgcna_relation_ms(gene_trait_corr, seq_annot, relation_corr, relation_corr_pvalue)

        chart.to_pdf()
        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*sample_tree*.pdf")
        for p in pdf_file:
            os.link(p, self.output_dir + "/prepare/" + os.path.basename(p))
        pdf_file = glob.glob(self.work_dir + "/wgcna.module*.pdf")
        for p in pdf_file:
            os.link(p, self.output_dir + "/module/" + os.path.basename(p))
        pdf_file = glob.glob(self.work_dir + "/wgcna*relation*.pdf")
        for p in pdf_file:
            os.link(p, self.output_dir + "/relate/" + os.path.basename(p))


    def add_module_result(self, wgcna_prepare_id):
        # create main table
        params = dict(
            task_id=self.option("task_id"),
            submit_location="wgcnamodule",
            task_type=2,
            wgcna_prepare_id=str(wgcna_prepare_id),
            level=self.option("level"),
            mergeCutHeight=self.option("mergeCutHeight"),
            power=self.power,
            minModuleSize=self.option("minModuleSize"),
            networkType=self.option("networkType"),
            minKMEtoStay=self.option("minKMEtoStay"),
        )
        params = json.dumps(params, sort_keys=True, separators=(',',':'))
        name = "PipelineModule" + '_' + self.option("level")[0].upper() + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=self._sheet.project_sn,
            task_id=self.option("task_id"),
            name=name,
            version="v1",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            wgcna_prepare_id=wgcna_prepare_id,
            desc='wgcna module identification analysis main table',
            level=self.option("level"),
            # category=self.option("category"),
            params=params,
            status="start",
            pipeline_id=ObjectId(self.option("main_id")),
            output_dir=self.get_workflow_output_dir() + '/module',
        )
        main_id = self.medical_transcriptome.insert_main_table('sg_wgcna_module', main_info)
        # add result info
        gene2name = self.gid2gname
        exp_matrix = self.wgcna_prepare.output_dir + '/exp_matrix_after_filtering.txt'
        self.dump_tool.add_module_detail(self.wgcna_module.work_dir, gene2name, exp_matrix, str(main_id))
        return main_id

    def add_relate_result(self, wgcna_module_id):
        # create main table record
        params = dict(
            task_id=self.option("task_id"),
            submit_location="wgcnarelate",
            task_type=2,
            wgcna_module_id=str(wgcna_module_id),
            level=self.option("level"),
            corr_method=self.option("corr_method"),
            trait_type=self.option("trait_type"),
            trait_path=self.params["trait_path"],
            file_id=self.option("file_id"),
        )
        params = json.dumps(params, sort_keys=True, separators=(',',':'))
        name = "PipelineRelate" + '_' + self.option("level")[0].upper() + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=self._sheet.project_sn,
            task_id=self.option("task_id"),
            name=name,
            version="v1",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            wgcna_module_id=wgcna_module_id,
            desc='wgcna relation analysis main table',
            level=self.option("level"),
            # category=self.option("category"),
            params=params,
            status="start",
            pipeline_id=ObjectId(self.option("main_id")),
            output_dir=self.get_workflow_output_dir() + '/relate',
        )
        main_id = self.medical_transcriptome.insert_main_table('sg_wgcna_relate', main_info)
        # add result info
        seq_annot = self.wgcna_module.work_dir + '/gene_module_detail.xls'
        self.dump_tool.add_relate_detail(self.wgcna_relate.work_dir, seq_annot, str(main_id))
        return main_id

    def add_network_result(self, wgcna_module_id, wgcna_relate_id):
        for each_tool in self.wgcna_network:
            params = dict(
                task_id=self.option("task_id"),
                submit_location="wgcnanetwork",
                task_type=2,
                wgcna_relate_id=str(wgcna_relate_id),
                level=self.option("level"),
                threshold=self.option("threshold"),
                module=each_tool.option("module"),
                top=each_tool.option("top"),
            )
            packed_params = json.dumps(params, sort_keys=True, separators=(',',':'))
            name = "PipelineNetwork" + '_' + each_tool.option("module").replace(",","_") + '_'
            # name = "PipelineNetwork" + '_' + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=self._sheet.project_sn,
                task_id=self.option("task_id"),
                name=name,
                version="v1",
                level=self.option("level"),
                # category=self.option("category"),
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                wgcna_module_id=wgcna_module_id,
                wgcna_relate_id=wgcna_relate_id,
                desc='wgcna network analysis',
                params=packed_params,
                pipeline_id=ObjectId(self.option("main_id")),
                output_dir=self.get_workflow_output_dir() + '/network/' + each_tool.option("module") +'.network.json',
                # output_dir=self.get_workflow_output_dir() + '/network.json',
                status="start")
            main_id = self.medical_transcriptome.insert_main_table("sg_wgcna_network", main_info)
            # add network detail
            self.dump_tool.add_network_detail(each_tool.output_dir, main_id)

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if re.match(r'tsanger:',workflow_output):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        return workflow_output

#   export file for tool
    def export_wgcna_g2n_matrix(self, exp_matrix):
        target_seqs = pd.read_table(exp_matrix, header=0, index_col=0).index
        task_id = self.option('task_id')
        annot_table = self.dump_tool.db['sg_annotation_query']
        try:
            annot_main = annot_table.find_one({"task_id": task_id, "type": "latest", "status": "end"})
        except:
            self.set_error("cannot find sg_annotation_query main table", code = "13703401")
        else:
            if annot_main is None:
                annot_main = annot_table.find_one({"task_id": task_id, "type": "origin", "status": "end"})
        if "main_id" not in annot_main:
            annot_main_id = annot_main['_id']
        else:
            annot_main_id = annot_main['main_id']
        annot_detail = self.dump_tool.db['sg_annotation_query_detail']
        level = self.option('level')
        if level[0].upper() == 'G':
            query_dict = dict(query_id=annot_main_id, is_gene=True)
            result_dict = dict(_id=0, gene_name=1, gene_id=1)
            result = annot_detail.find(query_dict, result_dict)
            gene2name = pd.DataFrame(list(result))
            gene2name.set_index('gene_id', inplace=True)
        if level[0].upper() == 'T':
            query_dict = dict(query_id=annot_main_id, )
            result_dict = dict(_id=0, gene_name=1, gene_id=1, transcript_id=1)
            result = annot_detail.find(query_dict, result_dict)
            gene2name = pd.DataFrame(list(result))
            gene2name.set_index('transcript_id', inplace=True)
            gene2name = gene2name.loc[:, ["gene_id","gene_name"]]
        gene2name = gene2name.loc[list(target_seqs), :]
        gene2name.reset_index(inplace=True)
        output = os.path.join(self.wgcna_module.output_dir, "seq_id2gene_name.txt")
        gene2name.to_csv(output, sep='\t', header=True, index=False)
        gene2name = pd.read_table(output, header=0)
        gene2name["gene_name"] = gene2name["gene_name"].fillna(gene2name["gene_id"])
        # gene2name = gene2name.fillna(method="pad", axis=1)
        gene2name.to_csv(output, sep='\t', header=True, index=False)
        return output
