# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import web
import json
import datetime
import unittest
import os
from mainapp.libs.signature import check_sig
from mainapp.libs.param_pack import *
from mainapp.models.mongo.denovo_rna_v2 import DenovoRnaV2
from mainapp.controllers.project.denovo_rna_v2_controller import DenovoRnaV2Controller


class DenovoAnnotationAction(DenovoRnaV2Controller):
    """
    denovo 注释重运行接口
    """
    def __init__(self):
        super(DenovoAnnotationAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get("HTTP_CLIENT")
        print data
        return_result = self.check_options(data)

        if return_result:
            var = []
            var.append(return_result)
            info = {"success": False, "info": return_result, "code": 'C1600301', "variables": var}
            return json.dumps(info)
        stat_info = self.denovo_rna_v2.get_main_info_by_record("sg_annotation_stat", task_id=str(data.task_id), type="origin")
        if not stat_info:
            info = {"success": False, "info": "stat_id不存在,请确认参数是否正确", "code": 'C1600302', "variables": ''}
            return json.dumps(info)
        print stat_info

        run_info = self.denovo_rna_v2.get_main_info_by_record("sg_annotation_stat", task_id=str(data.task_id), status="start", type="latest")
        if run_info:
            var = []
            var.append(run_info['name'])
            info = {"success": False, "info": "任务%s正在计算，请等待其完成在下次运行"%(run_info['name']), "code": 'C1600303', "variables": var}
            return json.dumps(info)

        latest_info = self.denovo_rna_v2.get_main_info_by_record("sg_annotation_stat", task_id=str(data.task_id), status="end",  type="latest")
        latest_main_id = ""
        try:
            group_info = self.denovo_rna_v2.get_main_info_by_record("sg_specimen_group", task_id=str(data.task_id))
            groups = map(str, group_info["category_names"])
            species = list()
            for specie in group_info["specimen_names"]:
                species.append(map(str, specie))
            group_dict = dict(zip(groups, species))
        except:
            group_dict = dict()
        gene_exp = self.denovo_rna_v2.get_main_info_by_record("sg_exp", task_id=str(data.task_id), status="end",  exp_level="G", exp_type="TPM")["main_id"]

        trans_exp = self.denovo_rna_v2.get_main_info_by_record("sg_exp", task_id=str(data.task_id), status="end",  exp_level="T", exp_type="TPM")["main_id"]
        task_info = self.denovo_rna_v2.get_task_info(data.task_id)
        trans2gene = task_info['assemble_t2g']
        if latest_info:
            latest_main_id = latest_info["_id"]

        taxon = "Animals"
        if stat_info.has_key("taxonomy"):
            taxon = stat_info["taxonomy"]
        else:
            pass

        params_json = {
            "task_id": str(data.task_id),
            "nr_evalue": str(data.nr_evalue),
            "nr_similarity": str(data.nr_similarity),
            "nr_identity": str(data.nr_identity),
            "swissprot_evalue": str(data.swissprot_evalue),
            "swissprot_similarity": str(data.swissprot_similarity),
            "swissprot_identity": str(data.swissprot_identity),
            "cog_evalue": str(data.cog_evalue),
            "cog_similarity": str(data.cog_similarity),
            "cog_identity": str(data.cog_identity),
            "kegg_evalue": str(data.kegg_evalue),
            "kegg_similarity": str(data.kegg_similarity),
            "kegg_identity": str(data.kegg_identity),
            "pfam_evalue": str(data.pfam_evalue),
            "submit_location": data.submit_location,
            "task_type": int(data.task_type)
        }

        if hasattr(data, "exclude_taxon"):
            params_json.update({
                "exclude_taxon": str(data.exclude_taxon)
            })

        main_table_name = "AnnotationStat_latest_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        mongo_data = [
            ("project_sn", stat_info["project_sn"]),
            ("task_id", stat_info["task_id"]),
            ("seq_type", "new"),
            ("status", "start"),
            ("name", main_table_name),
            ("database", "nr,swissprot,cog,kegg,pfam"),
            ("desc", "注释统计主表"),
            ("result_dir", ""),
            ("taxonomy", taxon),
            ("version", "v2"),
            ("type", "latest"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.denovo_rna_v2.insert_main_table("sg_annotation_stat", mongo_data)
        update_info = {str(main_table_id): "sg_annotation_stat"}
        origin_result = str(stat_info["result_dir"])
        if not origin_result.endswith('/'):
            origin_result += '/'
        origin_param = stat_info["params"]
        print "group_dict is {}".format(str(group_dict))
        options = {
            "nr_evalue": float(data.nr_evalue),
            "nr_similarity": float(data.nr_similarity),
            "nr_identity": float(data.nr_identity),
            "swissprot_evalue": float(data.swissprot_evalue),
            "swissprot_similarity": float(data.swissprot_similarity),
            "swissprot_identity": float(data.swissprot_identity),
            "cog_evalue": float(data.cog_evalue),
            "cog_similarity": float(data.cog_similarity),
            "cog_identity": float(data.cog_identity),
            "kegg_evalue": float(data.kegg_evalue),
            "kegg_similarity": float(data.kegg_similarity),
            "kegg_identity": float(data.kegg_identity),
            "pfam_evalue": float(data.pfam_evalue),
            "task_id": str(data.task_id),
            "stat_id": str(main_table_id),
            "last_id": str(latest_main_id),
            "trans2gene": str(trans2gene),
            "origin_result": self.use_s3(origin_result),
            "origin_param": origin_param,
            "taxonomy": taxon,
            "gene_exp": gene_exp,
            "trans_exp": trans_exp,
            "group_dict": json.dumps(group_dict),
            "update_info": json.dumps(update_info)
        }
        if hasattr(data, "exclude_taxon"):
            options.update({
                "exclude_taxon": str(data.exclude_taxon)
            })
        to_files = ["denovo_rna_v2.export_exp_matrix(gene_exp)",
                    "denovo_rna_v2.export_exp_matrix(trans_exp)"]

        self.set_sheet_data(name="denovo_rna_v2.report.denovo_annotation",
                            options=options,
                            #main_table_name="AnnotationStat/" + main_table_name,
                            main_table_name=main_table_name,
                            task_id=stat_info["task_id"],
                            project_sn=stat_info["project_sn"],
                            to_file=to_files, )
        task_info = super(DenovoAnnotationAction, self).POST()
        task_info["content"] = {"ids": {"id": str(main_table_id), "name": main_table_name}}
        return json.dumps(task_info)

    def check_options(self, data):
        """
        检查网页端传来的参数是否正确
         ["task_id", "nr_evalue", "nr_score", "nr_similarity", "nr_identity",
          "swissprot_evalue", "swissprot_score", "swissprot_similarity",
           "swissprot_identity", "submit_location"]
        """
        warn_info = list()
        if not (hasattr(data, "task_id")):
            warn_info.append("缺少参数task_id")

        if not (hasattr(data, "task_type")):
            warn_info.append("缺少参数task_type")

        if not (hasattr(data, "nr_evalue")):
            data.nr_evalue = 1e-3
        else:
            try:
                if float(data.nr_evalue) > 1e-3:
                    warn_info.append("NR E-value值必须小于1e-3")
            except:
                warn_info.append('输入的NR E-value不是数字 ')

        if not (hasattr(data, "nr_similarity")):
            data.nr_similarity = 0
        else:
            try:
                if float(data.nr_similarity) < 0 or float(data.nr_similarity) > 100:
                    warn_info.append("NR Similarity值需在0-100范围内")
            except:
                warn_info.append('NR similarity 值必须是数字')

        if not (hasattr(data, "nr_identity")):
            data.nr_identity = 0
        else:
            try:
                if float(data.nr_identity) < 0 or float(data.nr_identity) > 100:
                    warn_info.append("NR Identity值需在0-100范围内")
            except:
                warn_info.append('NR Identity值必须是数字')

        if not (hasattr(data, "swissprot_similarity")):
            data.swissprot_similarity = 0
        else:
            try:
                if float(data.swissprot_similarity) < 0 or float(data.swissprot_similarity) > 100:
                    warn_info.append("Swiss-Prot Similarity值需在0-100范围内")
            except:
                warn_info.append('Swiss-port Similarity 必须是数字')

        if not (hasattr(data, "swissprot_evalue")):
            data.swissprot_evalue = 1e-3
        else:
            try:
                if float(data.swissprot_evalue) > 1e-3:
                    warn_info.append("Swiss-Prot E-value值需小于1e-3")
            except:
                warn_info.append('Swiss-Prot E-value必须是数字')

        if not (hasattr(data, "swissprot_identity")):
            data.swissprot_identity = 0
        else:
            try:
                if float(data.swissprot_identity) < 0 or float(data.swissprot_identity) > 100:
                    warn_info.append("Swiss-Prot Identity值需在0-100范围内")
            except:
                warn_info.append('Swiss-Prot Identity 必须是数字')


        if not (hasattr(data, "cog_evalue")):
            data.cog_evalue = 1e-3
        else:
            try:
                if float(data.cog_evalue) > 1e-3:
                    warn_info.append("COG E-value值必须小于1e-3")
            except:
                warn_info.append('输入的COG E-value不是数字 ')

        if not (hasattr(data, "cog_similarity")):
            data.cog_similarity = 0
        else:
            try:
                if float(data.cog_similarity) < 0 or float(data.cog_similarity) > 100:
                    warn_info.append("COG Similarity值需在0-100范围内")
            except:
                warn_info.append('COG similarity 值必须是数字')

        if not (hasattr(data, "cog_identity")):
            data.cog_identity = 0
        else:
            try:
                if float(data.cog_identity) < 0 or float(data.cog_identity) > 100:
                    warn_info.append("COG Identity值需在0-100范围内")
            except:
                warn_info.append('COG Identity值必须是数字')

        if not (hasattr(data, "kegg_evalue")):
            data.kegg_evalue = 1e-3
        else:
            try:
                if float(data.kegg_evalue) > 1e-3:
                    warn_info.append("KEGG E-value值必须小于1e-3")
            except:
                warn_info.append('输入的KEGG E-value不是数字 ')

        if not (hasattr(data, "kegg_similarity")):
            data.kegg_similarity = 0
        else:
            try:
                if float(data.kegg_similarity) < 0 or float(data.kegg_similarity) > 100:
                    warn_info.append("KEGG Similarity值需在0-100范围内")
            except:
                warn_info.append('KEGG similarity 值必须是数字')

        if not (hasattr(data, "kegg_identity")):
            data.kegg_identity = 0
        else:
            try:
                if float(data.kegg_identity) < 0 or float(data.kegg_identity) > 100:
                    warn_info.append("KEGG Identity值需在0-100范围内")
            except:
                warn_info.append('KEGG Identity值必须是数字')

        if not (hasattr(data, "pfam_evalue")):
            data.pfam_evalue = 1e-3
        else:
            try:
                if float(data.pfam_evalue) > 1e-3:
                    warn_info.append("PFAM E-value值必须小于1e-3")
            except:
                warn_info.append('输入的PFAM E-value不是数字 ')

        if hasattr(data, "exclude_taxon"):
            taxons = str(data.exclude_taxon).split(",")
            for taxon in taxons:
                if taxon == "":
                    continue
                try:
                    int(taxon)
                except:
                    warn_info.append('taxon 必需为数字')
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.denovo_rna_v2.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.denovo_rna_v2.update_group_compare_is_use(data.task_id, data.control_id)
        return ";".join(warn_info)

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/denovo_rna_v2/denovo_annotation "
        cmd += "-b http://192.168.12.101:9090 "
        args = dict(
            task_id='denovo_rna_v2_upgrade',
            nr_evalue=0.00001,
            nr_similarity=40,
            nr_identity=80,
            swissprot_evalue=0.000001,
            swissprot_similarity=40,
            swissprot_identity=80,
            cog_evalue=0.00001,
            cog_similarity=40,
            cog_identity=80,
            kegg_evalue=1e-6,
            kegg_similarity=90,
            kegg_identity=90,
            exclude_taxon='946362,81529',
            pfam_evalue=1e-4,
            submit_location="annotationstat",
            task_type=2,
        )
        arg_names, arg_values = args.keys(), args.values()
        arg_values = [str(x) for x in arg_values]
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
