# -*- coding: utf-8 -*-
# __author__ = 'fwy'
from biocluster.config import Config
import datetime
import unittest
import json
import web
import os
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.medical_transcriptome import *
from mbio.api.to_file.medical_transcriptome import *
from bson.objectid import ObjectId
from collections import OrderedDict
from biocluster.config import Config
import pandas as pd
from biocluster.file import download


class SnpSearchAction(MedicalTranscriptomeController):

    def __init__(self):
        super(SnpSearchAction, self).__init__(instant=False)

    def check_target_geneset(self):
        # project_type = 'medical_transcriptome'
        # db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        target_genes = self.input_data.target_genes
        geneset_snp = self.input_data.geneset_snp
        id_type = self.input_data.id_type
        if geneset_snp == "" and target_genes == "":
            final_target_genes = "all"
        else:
            if geneset_snp != "":
                collection = self.db['sg_geneset_detail']
                results = collection.find_one({"geneset_id": ObjectId(geneset_snp)})
                try:
                    seq_list = results["seq_list"]
                except:
                    info = {'success': False, 'info': "No geneset: %s" % str(str(geneset_snp))}
                    return json.dumps(info)
                if target_genes != "":
                    target_genes = target_genes.split(",")
                    final_target_genes = set(target_genes) & set(seq_list)
                    if not final_target_genes:
                        info = {'success': False, 'info': "输入的基因不在选定的基因集中"}
                        return [False,json.dumps(info)]
                    else:
                        final_target_genes = list(final_target_genes)
                else:
                    final_target_genes = list(seq_list)
            else:
                target_genes = target_genes.split(",")
                final_target_genes = list(target_genes)
        return [True,final_target_genes]
    
    def check_filter_snp(self,final_target_genes):
        project_type = 'medical_transcriptome'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        task_id = self.input_data.task_id
        connect = db['sg_snp']
        result = connect.find_one({'task_id': task_id, 'status': 'end'})
        snp_method = result["params"]["method_type"]
        type = self.input_data.type
        if type == "snp":
            s3_snp_path = result["result_dir"] + "/snp_anno.xls"
            inter_dir = self.create_tmp_dir(self.input_data.task_id, "snp/")
            snp_path = os.path.join(inter_dir, "snp", snp_method + "snp_anno.xls")
            download(s3_snp_path, snp_path)
            # snp_path = self.download_from_s3(s3_snp_path, inter_dir=inter_dir)
        else:
            s3_snp_path = result["result_dir"] + "/indel_anno.xls"
            inter_dir = self.create_tmp_dir(self.input_data.task_id, "indel/")
            snp_path = os.path.join(inter_dir, snp_method + "indel_anno.xls")
            download(s3_snp_path, snp_path)
        id_type = self.input_data.type
        sample = self.input_data.sample
        region = self.input_data.region
        depth_comare = self.input_data.depth_comare
        depth = int(self.input_data.depth)
        snp_detail = pd.read_table(snp_path)
        # 过滤条件较多,分步过滤
        # step1:按照区域和种类分类
        search_df = snp_detail[
            (snp_detail["type"] == type) & (snp_detail["Anno"] == region)]
        # step2 按照基因过滤,如果客户没选择基因过滤,则会传入一个all，否则将是一个文件。
        if final_target_genes ==  "all":
            pass
        else:
            search_df = search_df[search_df["GENE(in or nearby)"].isin(final_target_genes)]
        # step3 按照样本过滤
        if sample == "all":
            search_df = search_df
            sample_depths =self.get_depth_columns(search_df)
            if depth_comare == "greater":
                search_df["min_depth"] = search_df.apply(lambda x: x[sample_depths].min(), axis=1)
                search_df = search_df[search_df["min_depth"] > depth]
                search_df = search_df.drop("min_depth", axis=1)
            elif depth_comare == "greateroreq":
                search_df["min_depth"] = search_df.apply(lambda x: x[sample_depths].min(), axis=1)
                search_df = search_df[search_df["min_depth"] >= depth]
                search_df = search_df.drop("min_depth", axis=1)
            elif depth_comare == "less":
                search_df["max_depth"] = search_df.apply(lambda x: x[sample_depths].max(), axis=1)
                search_df = search_df[search_df["max_depth"] < depth]
                search_df = search_df.drop("max_depth", axis=1)
            elif depth_comare == "lessoreq":
                search_df["max_depth"] = search_df.apply(lambda x: x[sample_depths].max(), axis=1)
                search_df = search_df[search_df["max_depth"] <= depth]
                search_df = search_df.drop("max_depth", axis=1)
            elif depth_comare == "equal":
                search_df["reamin"] = search_df.apply(self.get_all,args=(depth,sample_depths,), axis=1)
                search_df = search_df[search_df["reamin"] == 1]
                search_df = search_df.drop("reamin", axis=1)
        else:
            search_df = search_df[
                ["GENE(in or nearby)", "Gene name", "Gene description", "Chrom", "Start", "End", "Ref", "Alt",
                 "Total depth", "QUAL", "Anno", "MUT type", "MUT info", "type", sample + "_sampledepth",
                 sample + "genotype"]]
            selected_column = sample + "_sampledepth"
            if depth_comare == "greater":
                search_df = search_df[search_df[selected_column] > depth]
            elif depth_comare == "greateroreq":
                search_df = search_df[search_df[selected_column] >= depth]
            elif depth_comare == "less":
                search_df = search_df[search_df[selected_column] < depth]
            elif depth_comare == "lessoreq":
                search_df = search_df[search_df[selected_column] <= depth]
            elif depth_comare == "equal":
                search_df = search_df[search_df[selected_column] == depth]
        if search_df.shape[0] == 0:
            info = {'success': False, 'info': "没有满足符合要求的snp事件"}
            return [False, json.dumps(info)]
        elif search_df.shape[0] > 50000:
            info = {'success': False, 'info': "snp事件大于5w条不予展示"}
            return [False, json.dumps(info)]
        else:
            snp_path = os.path.join(inter_dir, snp_method + "filter_snp_anno.xls")
            return [True, snp_path]

    @check_sig
    def POST(self):
        # project_type = 'medical_transcriptome'
        # db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        data = web.input()
        print data
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ["snp_id","target_genes", 'geneset_snp', 'type',"id_type","sample","region",'depth_comare',"depth"]

        self.input_data = data
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2901101', 'variables': variables}
                print info
                return json.dumps(info)

        print("bengkuikaishi")
        print(str(datetime.datetime.now()))
        record_counts = self.db["sg_snp_search"].find(
            {'task_id': data.task_id, "status": {"$in": ["end","start"]}}).count()
        # if record_counts >= 4 :
        #     info = {"success": False, "info": u"已有四条运行记录，不可再投递"}
        #     return json.dumps(info)
        check_geneset = self.check_target_geneset()
        if check_geneset[0] is not True:
            return check_geneset[1]
        # check_filter_snp = self.check_filter_snp(check_geneset[1])
        # if check_filter_snp[0] is not True:
        #     return check_filter_snp[1]
        print("bengkuijieshu")
        print(str(datetime.datetime.now()))
        check_geneset = self.check_target_geneset()
        task_id = data.task_id
        task_info = self.medical_transcriptome.get_task_info(task_id=task_id)
        project_sn = task_info['project_sn']
        params_json = {
            'task_id': data.task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'snp_id':str(data.snp_id),
            'geneset_snp': data.geneset_snp,
            'target_genes': data.target_genes,
            'type':  data.type,
            'id_type': data.id_type,
            'sample': data.sample,
            'region': data.region,
            'depth_comare': data.depth_comare,
            'depth': data.depth,
        }
        # if hasattr(data,'sample'):
        #     params_json.update({
        #         "sample": data.sample
        #     })
        # if hasattr(data, 'group_dict'):
        #     group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        #     params_json.update({
        #         'group_dict': group_dict
        #     })
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))

        # self.medical_transcriptome.delete_search_records("sg_snp_search",data.task_id,10)#这个10的意思是最多保留10条
        name = "Snp_search" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=data.task_id,
            name=name,
            version="v1",
            desc='Snp Search main table',
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            params=params,
            status="start",
        )

        main_id = self.medical_transcriptome.insert_main_table('sg_snp_search', main_info)
        new_task_id = self.medical_transcriptome.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}

        connect = self.db["sg_snp"]
        result = connect.find_one({'task_id': task_id,'main_id':ObjectId(data.snp_id), 'status': 'end'})
        result_dir = result["result_dir"]
        options = {
            'result_path': result_dir,
            'target_genes':data.target_genes,
            "geneset_snp":data.geneset_snp,
            "type": data.type,
            "id_type": data.id_type,
            "sample": data.sample,
            "anno": task_id,
            'region': data.region,
            'depth_comare': data.depth_comare,
            'depth': data.depth,
            'update_info': json.dumps({str(main_id): "sg_snp_search"}),
            'main_id': str(main_id),
            'main_table_data': main_table_data,
        }
        to_files = ["medical_transcriptome.get_gene_detail_whole(anno)",
                    "medical_transcriptome.export_multi_gene_list_search(geneset_snp)"]


        # prepare to file
        task_name = 'medical_transcriptome.report.snp_search'

        self.set_sheet_data(name=task_name,
                            options=options,
                            to_file=to_files,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            project_sn=project_sn,
                            new_task_id=new_task_id,
                            task_id=data.task_id)
        task_info = super(SnpSearchAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        # task_info['group_dict'] = group_dict
        return json.dumps(task_info)

        # 更新基因集的使用信息
        # self.whole_transcriptome.insert_geneset_info(data.geneset_id, 'geneset_cluster', str(main_id))
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.medical_transcriptome.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.medical_transcriptome.update_group_compare_is_use(data.task_id, data.control_id)
        return json.dumps(task_info)



class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/medical_transcriptome/snp_search "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="8ju59of6lkdojsbreb482jveeb",
            task_type="2",
            submit_location="snp_search",
            sample = "all",
            target_genes ="",
            geneset_snp = "5facb20d17b2bf65566b8881",
            type ="snp",
            id_type="gene_id",
            depth_comare ="greater",
            depth = "0",
            region ="all",
            snp_id ="5facb47017b2bf655673e27c",
            # sample='BY4741_1,Ni_BY_1,H4K5R_1',
            # group_dict=json.dumps({"H1": ["H1581_1", "H1581_2", "H1581_3"], "H2": ["H1581_4","H1581_5", "H1581_6"]}).replace('"', '\\"')
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
