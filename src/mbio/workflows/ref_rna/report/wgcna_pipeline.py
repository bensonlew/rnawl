# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
from mainapp.models.mongo.ref_rna import RefRna
from collections import OrderedDict
from bson.objectid import ObjectId
import pandas as pd
import glob
import datetime
import os
import json
import re


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
        self.wgcna_prepare = self.add_tool("rna.wgcna.wgcna_prepare")
        self.wgcna_module = self.add_tool("rna.wgcna.wgcna_module")
        self.wgcna_relate = self.add_tool("rna.wgcna.wgcna_relate")
        # self.wgcna_network = self.add_tool("rna.wgcna.wgcna_network")
        self.wgcna_network = list()
        self.dump_tool = self.api.api("ref_rna.wgcna")
        self.ref_rna = RefRna()  # 用于创建主表

    def run(self):
        self.wgcna_prepare.on("end", self.run_wgcna_module)
        self.wgcna_module.on("end", self.run_wgcna_relate)
        self.wgcna_relate.on("end", self.run_wgcna_network)
        self.run_wgcna_prepare()
        super(WgcnaPipelineWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        # save workflow output path--建议标配
        self.logger.info("step1")
        workflow_output = self.get_workflow_output_dir()
        pipeline_params = self.ref_rna.get_main_info(self.option("main_id"),'sg_wgcna_pipeline')['params']
        params_dict = json.loads(pipeline_params)
        params_dict['power'] = self.power
        new_params = json.dumps(params_dict,sort_keys=True,separators=(',',':'))
        self.dump_tool.update_db_record('sg_wgcna_pipeline', self.option('main_id'), output_dir=workflow_output, params=new_params)
        # add result info
        wgcna_prepare_id = self.add_prepare_result()
        wgcna_module_id = self.add_module_result(wgcna_prepare_id)
        wgcna_relate_id = self.add_relate_result(wgcna_module_id)
        self.add_network_result(wgcna_module_id, wgcna_relate_id)
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
            target = os.path.join(self.output_dir + '/' + 'wgcna_module', base_name)
            os.link(each, target)
        for each in glob.glob(self.wgcna_relate.output_dir+'/*'):
            base_name = os.path.basename(each)
            target = os.path.join(self.output_dir + '/' + 'wgcna_relate', base_name)
            os.link(each, target)
        # link network result
        for each_tool in self.wgcna_network:
            for each in glob.glob(each_tool.output_dir+'/*'):
                base_name = os.path.basename(each)
                target = os.path.join(self.output_dir, 'wgcna_network', base_name)
                os.link(each, target)
        # upload result
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "wgcna_pipeline"],
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
        options = dict(
            datExpr=exp_matrix,
            MEs=eigengenes,
            traits=self.option("trait_path"),
            corType=self.option('corr_method'),
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
        self.gid2gname = self.export_wgcna_g2n_matrix(exp_matrix)
        target_modules = set()
        traits = trait_corr.columns
        for trait in traits:
            p_small_modules = corr_pvalue[corr_pvalue[trait] <= 0.1].index
            if not p_small_modules.shape[0]:
                p_small_modules = corr_pvalue[corr_pvalue[trait] <= 1].index
            select_module = trait_corr.loc[p_small_modules, [trait]].abs().idxmax()[0]
            target_modules.add(select_module[2:])
        for module in target_modules:
            tool = self.add_tool("rna.wgcna.wgcna_network")
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
        self.logger.info("step2")
        exp_info = self.ref_rna.get_main_info(self.option("exp_id"), 'sg_express')
        project_sn = exp_info["project_sn"]
        group_dict = json.loads(self.option('raw_group_dict'), object_pairs_hook=OrderedDict)
        # create main table record
        exp_info = json.loads(exp_info['params'])
        exp_level,exp_type,quant_method = self.option("exp_level")[0].upper(),exp_info['type'],exp_info['express_method']
        params = dict(
            task_id=self.option("task_id"),
            submit_location="wgcnaprepare",
            task_type=2,
            exp_id=self.option("exp_id"),
            group_id=self.option("group_id"),
            exp_level=self.option("exp_level"),
            group_dict=group_dict,
            me=self.option("me"),
            cv=self.option("cv"),
            geneset_id=self.option("geneset_id"),
        )
        params = json.dumps(params,sort_keys=True,separators=(',',':'))
        name = "PipelinePrepare" + '_' + exp_level + '_' + quant_method + '_' + exp_type.upper() + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=self.option("task_id"),
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            exp_id=self.option("exp_id"),
            desc='wgcna pre-processing analysis main table',
            type=self.option("exp_level"),
            params=params,
            pipeline_id=ObjectId(self.option("main_id")),
            status="start")
        main_id = self.ref_rna.insert_main_table('sg_wgcna_prepare', main_info)
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
        exp_matrix = workflow_output + '/wgcna_prepare/exp_matrix_after_filtering.txt'
        self.dump_tool.update_db_record('sg_wgcna_prepare', main_id, group_info=group_info, exp_matrix=exp_matrix)
        self.logger.info(main_id)
        return main_id

    def add_module_result(self, wgcna_prepare_id):
        # create main table
        self.logger.info("step3")
        params = dict(
            task_id=self.option("task_id"),
            submit_location="wgcnamodule",
            task_type=2,
            wgcna_prepare_id=str(wgcna_prepare_id),
            exp_level=self.option("exp_level"),
            mergeCutHeight=self.option("mergeCutHeight"),
            power=self.power,
            minModuleSize=self.option("minModuleSize"),
            networkType=self.option("networkType"),
            minKMEtoStay=self.option("minKMEtoStay"),
        )
        params = json.dumps(params, sort_keys=True, separators=(',',':'))
        name = "PipelineModule" + '_' + self.option("exp_level")[0].upper() + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=self._sheet.project_sn,
            task_id=self.option("task_id"),
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            wgcna_prepare_id=wgcna_prepare_id,
            desc='wgcna module identification analysis main table',
            type=self.option("exp_level"),
            params=params,
            status="start",
            pipeline_id=ObjectId(self.option("main_id")),
            output_dir=self.get_workflow_output_dir() + '/wgcna_module',
        )
        main_id = self.ref_rna.insert_main_table('sg_wgcna_module', main_info)
        # add result info
        gene2name = self.gid2gname
        exp_matrix = self.wgcna_prepare.output_dir + '/exp_matrix_after_filtering.txt'
        self.dump_tool.add_module_detail(self.wgcna_module.work_dir, gene2name, exp_matrix, str(main_id))
        self.logger.info(main_id)
        return main_id

    def add_relate_result(self, wgcna_module_id):
        # create main table record
        self.logger.info("step4")
        params = dict(
            task_id=self.option("task_id"),
            submit_location="wgcnarelate",
            task_type=2,
            wgcna_module_id=str(wgcna_module_id),
            exp_level=self.option("exp_level"),
            corr_method=self.option("corr_method"),
            trait_type=self.option("trait_type"),
            trait_path=self.option("trait_path"),
            file_id=self.option("file_id"),
        )
        params = json.dumps(params, sort_keys=True, separators=(',',':'))
        name = "PipelineRelate" + '_' + self.option("exp_level")[0].upper() + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=self._sheet.project_sn,
            task_id=self.option("task_id"),
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            wgcna_module_id=wgcna_module_id,
            desc='wgcna relation analysis main table',
            type=self.option("exp_level"),
            params=params,
            status="start",
            pipeline_id=ObjectId(self.option("main_id")),
            output_dir=self.get_workflow_output_dir() + '/wgcna_relate',
        )
        main_id = self.ref_rna.insert_main_table('sg_wgcna_relate', main_info)
        # add result info
        seq_annot = self.wgcna_module.work_dir + '/gene_module_detail.xls'
        self.dump_tool.add_relate_detail(self.wgcna_relate.work_dir, seq_annot, str(main_id))
        self.logger.info(main_id)
        return main_id

    def add_network_result(self, wgcna_module_id, wgcna_relate_id):
        self.logger.info("step5")
        for each_tool in self.wgcna_network:
            params = dict(
                task_id=self.option("task_id"),
                submit_location="wgcnanetwork",
                task_type=2,
                wgcna_relate_id=str(wgcna_relate_id),
                exp_level=self.option("exp_level"),
                threshold=self.option("threshold"),
                module=each_tool.option("module"),
                top=each_tool.option("top"),
            )
            packed_params = json.dumps(params, sort_keys=True, separators=(',',':'))
            name = "PipelineNetwork" + '_' + each_tool.option("module").replace(",","_") + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=self._sheet.project_sn,
                task_id=self.option("task_id"),
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                wgcna_module_id=wgcna_module_id,
                wgcna_relate_id=wgcna_relate_id,
                desc='wgcna network analysis',
                params=packed_params,
                pipeline_id=ObjectId(self.option("main_id")),
                output_dir=self.get_workflow_output_dir() + '/wgcna_network/' + each_tool.option("module") +'.network.json',
                status="start")
            main_id = self.ref_rna.insert_main_table("sg_wgcna_network", main_info)
            self.logger.info(main_id)
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
        sg_task = self.dump_tool.db['sg_task']
        sg_task_info = sg_task.find_one({"task_id": task_id})
        if 'is_demo' in sg_task_info and int(sg_task_info['is_demo']) != 0:
            task_id = sg_task_info['demo_id']
        annot_main = annot_table.find_one({"task_id": task_id})
        if "main_id" not in annot_main:
            annot_main_id = annot_main['_id']
        else:
            annot_main_id = annot_main['main_id']
        annot_detail = self.dump_tool.db['sg_annotation_query_detail']
        anno_type = self.option('exp_level')
        query_dict = dict(query_id=annot_main_id,anno_type=anno_type,)
        result_dict = dict(_id=0,gene_name=1,gene_id=1)
        if anno_type == 'transcript':
            result_dict.update({"transcript_id": 1})
        result = annot_detail.find(query_dict,result_dict)
        print(query_dict)
        print(result_dict)
        gene2name = pd.DataFrame(list(result))
        print(gene2name.head())
        if anno_type == 'transcript':
            gene2name.set_index('transcript_id',inplace=True)
            gene2name = gene2name.loc[:,["gene_id","gene_name"]]
        else:
            gene2name.set_index('gene_id',inplace=True)
        gene2name = gene2name.loc[list(target_seqs), :]
        gene2name.reset_index(inplace=True)
        output = os.path.join(self.wgcna_module.output_dir, "seq_id2gene_name.txt")
        gene2name.to_csv(output, sep='\t', header=True, index=False)
        gene2name = pd.read_table(output, header=0)
        gene2name.fillna(method="pad", axis=1, inplace=True)
        gene2name.to_csv(output, sep='\t', header=True, index=False)
        return output
