# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
from mainapp.models.mongo.itraq_and_tmt import ItraqTmt
from collections import OrderedDict
from bson.objectid import ObjectId
import pandas as pd
import glob
import datetime
import os
import json
import re
from mbio.packages.itraq_and_tmt.wgcna.wgcna import Wgcna
from mbio.packages.itraq_and_tmt.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder

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
        self.wgcna_prepare = self.add_tool("itraq_and_tmt.wgcna.wgcna_prepare")
        self.wgcna_module = self.add_tool("itraq_and_tmt.wgcna.wgcna_module")
        self.wgcna_relate = self.add_tool("itraq_and_tmt.wgcna.wgcna_relate")
        # self.wgcna_network = self.add_tool("rna.wgcna.wgcna_network")
        self.wgcna_network = list()
        self.dump_tool = self.api.api("itraq_and_tmt.wgcna")
        self.itraq_and_tmt = Wgcna() # add by fwy 20210108 mongo库升级
        # self.itraq_and_tmt = ItraqTmt()  # 用于创建主表
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/6_Wgcna')
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
        self.run_wgcna_prepare()
        super(WgcnaPipelineWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        # save workflow output path--建议标配
        workflow_output = self.get_workflow_output_dir()
        pipeline_params = self.itraq_and_tmt.get_main_info(self.option("main_id"),'sg_wgcna_pipeline', self.option("task_id"))['params']
        params_dict = json.loads(pipeline_params)
        params_dict['power'] = self.power
        new_params = json.dumps(params_dict,sort_keys=True,separators=(',',':'))
        self.dump_tool.update_db_record('sg_wgcna_pipeline', self.option('main_id'), main_id=ObjectId(self.option('main_id')),
                                        output_dir=workflow_output, status='end', params=new_params)
        # add result info
        wgcna_prepare_id = self.add_prepare_result()
        wgcna_module_id = self.add_module_result(wgcna_prepare_id)
        wgcna_relate_id = self.add_relate_result(wgcna_module_id)
        self.add_network_result(wgcna_module_id, wgcna_relate_id)
        self.end()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir+'/'
        chart.output_dir = self.work_dir+'/'
        
        # group_dict = {"An10": ["An10_1", "An10_2", "An10_3"], "An5": ["An5_1", "An5_2", "An5_3"], "An_cPP": ["An_cPP_1", "An_cPP_2", "An_cPP_3"], "Con": ["Con_1", "Con_2", "Con_3"]}
        group_dict = json.loads(self.option('raw_group_dict'), object_pairs_hook=OrderedDict)
        sample_tree = os.path.join(self.work_dir, "WgcnaPrepare", "sample.cluster.dendrogram.txt")
        #
        module_tree = os.path.join(self.work_dir, "WgcnaModule", "unmerged.module_corr.dendrogram.txt")
        #
        module_corr = os.path.join(self.work_dir, "WgcnaModule", "module_corr.matrix.xls")
        module_corr_tree = os.path.join(self.work_dir, "WgcnaModule", "module_corr.tree.txt")
        #
        module_stat = os.path.join(self.work_dir, "WgcnaModule", "module_size.stat.xls")
        #
        relation_corr = os.path.join(self.work_dir, "WgcnaRelate", "module_trait.correlation.xls")
        relation_corr_pvalue = os.path.join(self.work_dir, "WgcnaRelate", "module_trait.correlation_pvalues.xls")
        #
        gene_trait_corr = os.path.join(self.work_dir, "WgcnaRelate", "protein_trait.correlation.xls")
        seq_annot = os.path.join(self.work_dir, "WgcnaModule", "protein_module_detail.xls")
        #
        wgcna_prepare_curve = os.path.join(self.work_dir, "WgcnaPrepare", "scale_free_analysis.xls")

        #
        if os.path.exists(sample_tree):
            chart.chart_wgcna_sample_tree(sample_tree,group_dict)
        if os.path.exists(wgcna_prepare_curve):
            chart.chart_wgcna_prepare_curve(wgcna_prepare_curve)
        #
        if os.path.exists(module_stat):
            chart.chart_wgcna_module_column(module_stat)
        if os.path.exists(module_corr) and os.path.exists(module_corr_tree):
            chart.chart_wgcna_module_corr(module_corr,module_corr_tree)
        if os.path.exists(module_tree):
            chart.chart_wgcna_module_tree(module_tree)
        #
        if os.path.exists(relation_corr) and os.path.exists(relation_corr_pvalue) and os.path.exists(module_stat):
            chart.chart_wgcna_relation_corr(relation_corr,relation_corr_pvalue,module_stat)
        if os.path.exists(gene_trait_corr) and os.path.exists(seq_annot) and os.path.exists(relation_corr) and os.path.exists(relation_corr_pvalue):
            chart.chart_wgcna_relation_ms(gene_trait_corr,seq_annot,relation_corr,relation_corr_pvalue)


        chart.to_pdf()
        # # move pdf to result dir
        def os_link(ori_filename, dirname, filename):
            if os.path.exists(os.path.join(self.work_dir, ori_filename)):
                os.link(os.path.join(self.work_dir, ori_filename), os.path.join(self.output_dir, dirname, filename))
        for i,ii,iii in [\
        ["wgcna.sample_tree.heat_corr.pdf","wgcna_prepare","cluster.pdf"],#chart_wgcna_sample_tree\
        ["wgcna_prepare_average.wgcnascatter.pdf","wgcna_prepare","mean.pdf"],#chart_wgcna_prepare_curve\
        ["wgcna_prepare_adapt.wgcnascatter.pdf","wgcna_prepare","scale.pdf"],#chart_wgcna_prepare_curve\
        #
        # ["?????","wgcna_module","dendrogram.pdf"],\
        ["wgcna.module_stat.column.pdf","wgcna_module","num.pdf"],#chart_wgcna_module_column\
        ["wgcna.mofule_corr.heat_corr.pdf","wgcna_module","relate.pdf"],#chart_wgcna_module_corr\
        ["wgcna.module_tree.heat_corr.pdf","wgcna_module","modulecluster.pdf"],#chart_wgcna_module_tree\
        #
        ["wgcna.relation_heat.heat_corr.pdf","wgcna_relate","module_trait_cor.pdf"],#chart_wgcna_relation_corr\
        # ["?????","wgcna_relate","protein_trait_cor.pdf"]\
        ]:
            os_link(i, ii, iii)
        for pdf_file in glob.glob(self.work_dir + "/wgcna.*.relation_ms.column_conf.pdf"):#chart_wgcna_relation_ms
            file_name = os.path.basename(pdf_file)
            os_link(file_name, "wgcna_relate", file_name[6:-28]+"_MS_bar.pdf")
        for pdf_file in glob.glob(self.work_dir + "/wgcnarelation_ms.*.scatter.pdf"):#chart_wgcna_relation_ms
            file_name = os.path.basename(pdf_file)
            os_link(file_name, "wgcna_relate", file_name[17:-12]+"_MSGS_plot.pdf")
        
        pdf_file = glob.glob(self.work_dir + "/WgcnaModule/*dendrogram.pdf")[0]
        os.link(pdf_file, os.path.join(self.output_dir, "wgcna_module", "dendrogram.pdf"))
        pdf_file = glob.glob(self.work_dir + "/WgcnaRelate/*trait.correlation.pdf")[0]
        os.link(pdf_file, os.path.join(self.output_dir, "wgcna_relate", "protein_trait_cor.pdf"))

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
                base_name = each_tool.option("module") + '.' + os.path.basename(each)
                target = os.path.join(self.output_dir, 'wgcna_network', base_name)
                os.link(each, target)
        #
        self.chart()
        # upload result
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["6_Wgcna", "", "WGCNA",0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "wgcna_pipeline", 0,  "211248"],
            ["./wgcna_prepare", "", "数据预处理", 0],
            ["./wgcna_prepare/cluster.pdf", "", "样本聚类", 0],
            ["./wgcna_prepare/scale.pdf", "", "无尺度容适曲线", 0],
            ["./wgcna_prepare/mean.pdf", "", "无尺度平均连通度曲线", 0],
            ["./wgcna_module", "", "模块识别", 0],
            ["./wgcna_module/dendrogram.pdf", "", "模块分类树", 0],
            ["./wgcna_module/num.pdf", "", "模块成员统计图", 0],
            ["./wgcna_module/relate.pdf", "", "模块相关性图", 0],
            ["./wgcna_module/modulecluster.pdf", "", "模块聚类图", 0],
            ["./wgcna_relate", "", "模块分析", 0],
            ["./wgcna_relate/module_trait_cor.pdf", "", "模块与表型相关性热图", 0],
            ["./wgcna_relate/protein_trait_cor.pdf", "", "蛋白与表型相关性热图", 0],
            ["./wgcna_relate/*MS_bar.pdf", "", "MS分析柱图", 0],
            ["./wgcna_relate/*MSGS_plot.pdf", "", "MM-GS散点图", 0],
            ["./wgcna_network", "", "可视化分析", 0],
            ["./wgcna_network/Network.pdf", "", "可视化网络图", 0],
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
        eigenproteins = os.path.join(self.wgcna_module.output_dir, "eigenproteins.txt")
        module_Rdata = os.path.join(self.wgcna_module.output_dir, "blockwiseModules_result.RData")
        options = dict(
            datExpr=exp_matrix,
            MEs=eigenproteins,
            traits=self.option("trait_path"),
            corType=self.option('corr_method'),
            block_Rdata=module_Rdata
        )
        self.wgcna_relate.set_options(options)
        self.wgcna_relate.run()

    def run_wgcna_network(self):
        import shutil
        trait_corr_f = glob.glob(self.wgcna_relate.output_dir+'/module_trait.correlation.xls')[0]
        corr_pvalue_f = glob.glob(self.wgcna_relate.output_dir+'/module_trait.correlation_pvalues.xls')[0]
        trait_corr = pd.read_table(trait_corr_f, index_col=0, header=0)
        corr_pvalue = pd.read_table(corr_pvalue_f, index_col=0, header=0)
        if "MEgrey" in trait_corr.index:
            trait_corr = trait_corr.loc[trait_corr.index.drop("MEgrey"), :]
            corr_pvalue = corr_pvalue.loc[corr_pvalue.index.drop("MEgrey"), :]
        exp_matrix = os.path.join(self.wgcna_prepare.output_dir,"exp_matrix_after_filtering.txt")
        self.gid2gname = self.export_wgcna_g2n_matrix(exp_matrix)
        # self.export_wgcna_g2n_matrix(exp_matrix)
        target_modules = set()
        traits = trait_corr.columns
        for trait in traits:
            p_small_modules = corr_pvalue[corr_pvalue[trait] <= 0.1].index
            if not p_small_modules.shape[0]:
                p_small_modules = corr_pvalue[corr_pvalue[trait] <= 1].index
            select_module = trait_corr.loc[p_small_modules, [trait]].abs().idxmax()[0]
            target_modules.add(select_module[2:])
        step3output = os.path.join(self.work_dir, 'step3output')
        if os.path.exists(step3output):
            shutil.rmtree(step3output)
        os.makedirs(step3output)
        for file in glob.glob(self.wgcna_relate.output_dir + '/*'):
            source = file
            base = os.path.basename(file).replace('protein', 'gene')
            shutil.copy(source, os.path.join(step3output, base))
        step2output = os.path.join(self.work_dir, 'step2output')
        if os.path.exists(step2output):
            shutil.rmtree(step2output)
        os.makedirs(step2output)
        for file in glob.glob(self.wgcna_module.output_dir + '/*'):
            source = file
            base = os.path.basename(file).replace('protein', 'gene')
            shutil.copy(source, os.path.join(step2output, base))
        for module in target_modules:
            tool = self.add_tool("itraq_and_tmt.wgcna.wgcna_network")
            options = dict(
                module=module,
                threshold=self.option('threshold'),
                top=self.option("top"),
                step3output=step3output,
                step2output=step2output,
            )
            tool.set_options(options)
            self.wgcna_network.append(tool)
        self.on_rely(self.wgcna_network, self.set_db, 'wgcnanetwork')
        for tool in self.wgcna_network:
            tool.run()

#   dump data to db for each step
    def add_prepare_result(self):
        # create main table
        exp_info = self.itraq_and_tmt.get_main_info(self.option("exp_id"), 'sg_express', self.option("task_id"))
        project_sn = exp_info["project_sn"]
        group_dict = json.loads(self.option('raw_group_dict'), object_pairs_hook=OrderedDict)
        # create main table record
        params = dict(
            task_id=self.option("task_id"),
            submit_location="wgcnaprepare",
            task_type=2,
            exp_id=self.option("exp_id"),
            group_id=self.option("group_id"),
            group_dict=group_dict,
            me=self.option("me"),
            cv=self.option("cv"),
            proteinset_id=self.option("proteinset_id"),
        )
        params = json.dumps(params,sort_keys=True,separators=(',',':'))
        name = "PipelinePrepare" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=self.option("task_id"),
            name=name,
            version="v2",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            exp_id=self.option("exp_id"),
            desc='wgcna pre-processing analysis main table',
            params=params,
            pipeline_id=ObjectId(self.option("main_id")),
            status="start")
        main_id = self.itraq_and_tmt.insert_main_table('sg_wgcna_prepare', main_info)
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
        return main_id

    def add_module_result(self, wgcna_prepare_id):
        # create main table
        params = dict(
            task_id=self.option("task_id"),
            submit_location="wgcnamodule",
            task_type=2,
            wgcna_prepare_id=str(wgcna_prepare_id),
            mergeCutHeight=self.option("mergeCutHeight"),
            power=self.power,
            minModuleSize=self.option("minModuleSize"),
            networkType=self.option("networkType"),
            minKMEtoStay=self.option("minKMEtoStay"),
        )
        params = json.dumps(params, sort_keys=True, separators=(',',':'))
        name = "PipelineModule" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=self._sheet.project_sn,
            task_id=self.option("task_id"),
            name=name,
            version="v2",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            wgcna_prepare_id=wgcna_prepare_id,
            desc='wgcna module identification analysis main table',
            params=params,
            status="start",
            pipeline_id=ObjectId(self.option("main_id")),
            output_dir=self.get_workflow_output_dir() + '/wgcna_module',
        )
        main_id = self.itraq_and_tmt.insert_main_table('sg_wgcna_module', main_info)
        # add result info
        protein2name = self.gid2gname
        exp_matrix = self.wgcna_prepare.output_dir + '/exp_matrix_after_filtering.txt'
        self.dump_tool.add_module_detail(self.wgcna_module.work_dir, protein2name, exp_matrix, str(main_id))
        return main_id

    def add_relate_result(self, wgcna_module_id):
        # create main table record
        params = dict(
            task_id=self.option("task_id"),
            submit_location="wgcnarelate",
            task_type=2,
            wgcna_module_id=str(wgcna_module_id),
            corr_method=self.option("corr_method"),
            trait_type=self.option("trait_type"),
            trait_path=self.option("trait_path"),
            file_id=self.option("file_id"),
        )
        params = json.dumps(params, sort_keys=True, separators=(',',':'))
        name = "PipelineRelate" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=self._sheet.project_sn,
            task_id=self.option("task_id"),
            name=name,
            version="v2",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            wgcna_module_id=wgcna_module_id,
            desc='wgcna relation analysis main table',
            params=params,
            status="start",
            pipeline_id=ObjectId(self.option("main_id")),
            output_dir=self.get_workflow_output_dir() + '/wgcna_relate',
        )
        main_id = self.itraq_and_tmt.insert_main_table('sg_wgcna_relate', main_info)
        # add result info
        seq_annot = self.wgcna_module.work_dir + '/protein_module_detail.xls'
        self.dump_tool.add_relate_detail(self.wgcna_relate.work_dir, seq_annot, str(main_id))
        return main_id

    def add_network_result(self, wgcna_module_id, wgcna_relate_id):
        for each_tool in self.wgcna_network:
            params = dict(
                task_id=self.option("task_id"),
                submit_location="wgcnanetwork",
                task_type=2,
                wgcna_relate_id=str(wgcna_relate_id),
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
                version="v2",
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                wgcna_module_id=wgcna_module_id,
                wgcna_relate_id=wgcna_relate_id,
                desc='wgcna network analysis',
                params=packed_params,
                pipeline_id=ObjectId(self.option("main_id")),
                output_dir=self.get_workflow_output_dir() + '/wgcna_network/' + each_tool.option("module") +'.network.json',
                status="start")
            main_id = self.itraq_and_tmt.insert_main_table("sg_wgcna_network", main_info)
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
        query_dict = dict(query_id=annot_main_id)
        result_dict = dict(_id=0, description=1, accession_id=1)
        result = annot_detail.find(query_dict, result_dict)
        protein2name = pd.DataFrame(list(result))
        protein2name.set_index('accession_id', inplace=True)
        protein2name = protein2name.loc[list(target_seqs), :]
        protein2name.reset_index(inplace=True)
        output = os.path.join(self.wgcna_module.output_dir, "seq_id2protein_name.txt")
        protein2name.to_csv(output, sep='\t', header=True, index=False)
        protein2name = pd.read_table(output, header=0)
        protein2name = protein2name.fillna(method="pad", axis=1)
        protein2name.to_csv(output, sep='\t', header=True, index=False)
        return output
