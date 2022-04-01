# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20190320
import os
import web
import json
import datetime
from biocluster.config import Config
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from bson.objectid import ObjectId
from mainapp.controllers.project.dna_controller import DnaController


class ChromosomeDistributionAction(DnaController):
    """
    染色体分布图接口
    """
    def __init__(self):
        super(ChromosomeDistributionAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        # subname = []
        data = web.input()
        print data
        params = ["marker_type", "data_type", "win_step", "analysis_object", "graphic_style", "task_id", "task_type",
                  "submit_location", "chongmingming_result"]
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数!" % param}
                return json.dumps(info)
        if data.marker_type not in ['Gene', 'SNP+InDel', "SNP", "InDel", "CNV", "SV", "SSR"]:
            info = {"success": False, "info": "marker_type：%s类型不合法" % data.marker_type}
            return json.dumps(info)
        if data.graphic_style not in ['heatmap', 'area']:
            info = {"success": False, "info": "图形样式参数%s不正确" % data.graphic_style}
            return json.dumps(info)
        if data.data_type not in ['before', 'after', 'ref']:
            info = {"success": False, "info": "分析数据来源参数%s不正确" % data.data_type}
            return json.dumps(info)
        if data.marker_type != 'SNP+InDel' and data.analysis_object == 'all':
            info = {"success": False, "info": "标记类型不为'SNP+InDel'时选择对象参数不能为all"}
            return json.dumps(info)
        if data.data_type == 'ref' and data.marker_type not in ["Gene", "SSR"]:
            info = {"success": False, "info": "当数据类型是ref的时候，标记类型必须为Gene 或者SSR！"}
            return json.dumps(info)
        if data.marker_type == 'CNV' and data.data_type == 'before' and data.analysis_object == '':
            info = {"success": False, "info": "CNV原始vcf进行分析的时候，必须要选择样本！"}
            return json.dumps(info)
        db = self.set_project_type(data)
        result = Dna(db).find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        # noinspection PyBroadException
        try:
            pop_final_vcf = result["pop_final_vcf"]
            project_sn = result["project_sn"]
            member_id = result["member_id"]
            pop_sv_vcf = result['pop_sv_vcf']
            # pop_ssr_vcf = result['pop_ssr_vcf']
            cnv_anno_path = result['cnv_anno_path']
            genome_version_id = result['genome_version_id']
        except:
            info = {"success": False, "info": "sg_task表里没有找到对应信息，请检查!"}
            return json.dumps(info)
        if data.data_type == 'before':
            if data.marker_type in ['SNP+InDel', "SNP", "InDel"]:  # pop_final_vcf
                pos_file = Dna(db).set_file_path(data.task_id, pop_final_vcf, data.client)
            elif data.marker_type == 'CNV':
                cnv_file = os.path.join(cnv_anno_path, "{}.cnv.anno.xls".format(data.analysis_object))
                pos_file = Dna(db).set_file_path(data.task_id, cnv_file, data.client)
            elif data.marker_type == 'SSR':
                result = Dna(db).find_one(collection="sg_ssr_marker",
                                          query_dic={"task_id": data.task_id, "name": {'$ne': "Reference"}})
                # noinspection PyBroadException
                try:
                    pop_ssr_vcf = result['ssr_vcf']
                except:
                    info = {"success": False, "info": "请确认是否做了ssr标记分析!--正常结束"}
                    return json.dumps(info)
                pos_file = Dna(db).set_file_path(data.task_id, pop_ssr_vcf, data.client)
            elif data.marker_type == 'SV':
                pos_file = Dna(db).set_file_path(data.task_id, pop_sv_vcf, data.client)
            else:
                info = {"success": False, "info": "%s不合法！" % data.marker_type}
                return json.dumps(info)
        elif data.data_type == 'after':
            collection = 'sg_variant_site_table'
            # if data.marker_type in ['SNP+InDel', "SNP", "InDel"]:
            #     origin = "variant"
            # elif data.marker_type == 'CNV':
            #     origin = "cnv"
            # elif data.marker_type == 'SV':
            #     origin = 'sv'
            # elif data.marker_type == 'SSR':
            #     origin = 'ssr'
            # else:
            #     info = {"success": False, "info": "%s不合法！" % data.marker_type}
            #     return json.dumps(info)
            result = Dna(db).find_one(collection=collection, query_dic={"_id": ObjectId(data.analysis_object)})
            # noinspection PyBroadException
            try:
                vcf_path = result["vcf_path"]
                # if data.subname:
                #     subname = result['subname']
            except:
                info = {"success": False, "info": "%s表里没有找到对应信息，请检查!" % collection}
                return json.dumps(info)
            # if data.subname:
            #     if data.subname in subname:
            #         index_ = subname.index(data.subname)
            #         pos_file = Dna(db).set_file_path(data.task_id, vcf_path[index_], data.client)
            #     else:
            #         info = {"success": False, "info": "%s不在%s中！" % (data.subname, subname)}
            #         return json.dumps(info)
            # else:
            #     if isinstance(vcf_path, list):
            #         pos_file = Dna(db).set_file_path(data.task_id, vcf_path[0], data.client)
            #     else:
            #         pos_file = Dna(db).set_file_path(data.task_id, vcf_path, data.client)
            pos_file = Dna(db).set_file_path(data.task_id, vcf_path, data.client)
        else:
            result = Dna(db).find_one(collection="sg_species_version", query_dic={"_id": genome_version_id})
            # noinspection PyBroadException
            try:
                genome_path = result["ssr_path"]
            except:
                info = {"success": False, "info": "sg_species_version表里没有找到对应信息，请检查!"}
                return json.dumps(info)
            if data.marker_type == 'Gene':
                gff_file = Config().SOFTWARE_DIR + "/database/dna_geneome/{}ref/genes.gff".format(genome_path)
                if not os.path.exists(gff_file):
                    gff_file = Config().SOFTWARE_DIR + "/database/dna_geneome/{}ref/genes.gtf".format(genome_path)
                pos_file = Dna(db).set_file_path(data.task_id, gff_file, data.client)
            else:
                ssr_path = Config().SOFTWARE_DIR + "/database/dna_geneome/{}ssr.ref.result.xls".format(genome_path)
                pos_file = Dna(db).set_file_path(data.task_id, ssr_path, data.client)
        win_step = int(data.win_step)
        params_json = {
            "marker_type": data.marker_type,
            "data_type": data.data_type,
            "win_step": str(win_step),
            "analysis_object": data.analysis_object,
            "graphic_style": data.graphic_style,
            "task_type": int(data.task_type),
            "submit_location": data.submit_location,
            "task_id": data.task_id,
            "chongmingming_result": data.chongmingming_result
            # "subname": data.subname
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        main_table_name = Dna(db).set_main_table_name("Chromosome_distribution", data.chongmingming_result)
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "chromosome_distribution主表！"),
            ('member_id', member_id),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna(db).insert_main_table(collection="sg_marker_distribution", data=mongo_data)
        Dna(db).update_db_record(collection="sg_marker_distribution", query_dict={"_id": main_id},
                                 update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_marker_distribution"}
        options = {
            "pos_file": pos_file,
            "marker_type": data.marker_type,
            "data_type": data.data_type,
            "win_step": win_step,
            "update_info": json.dumps(update_info),
            "main_id": str(main_id),
            "analysis_object": data.analysis_object,
            "graphic_style": data.graphic_style,
            "project_sn": project_sn,
            "task_id": data.task_id
        }
        if hasattr(data, "project_type"):
            options["project_type"] = data.project_type
        self.set_sheet_data(name="wgs_v2.chromosome_distribution", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id, main_table_name="chromosome_distribution/" + main_table_name,
                            options=options, params=params, db_type=db, target_output=True,
                            analysis_name='chromosomedistribution', module_type='tool')
        task_info = super(ChromosomeDistributionAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

    def set_project_type(self, data):
        if not hasattr(data, "project_type"):
            db = "dna_wgs_v2"
        elif data.project_type == "dna_gmap":
            db = "dna_gmap"
        elif data.project_type == "dna_wgs":
            db = 'dna_wgs'
        else:
            db = data.project_type
        return db


if __name__ == '__main__':
    """
    python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post s/wgs_v2/chromosome_distribution
    -c client03 -b http://192.168.12.101:9090
    -n "marker_type;data_type;win_step;analysis_object;graphic_style;task_id;task_type;submit_location;
    chongmingming_result;subname" -d "Gene;ref;100000000;;heatmap;wgs_v2;1;chr_distributon;;"
    """
    import argparse
    parser = argparse.ArgumentParser(description="用于对染色体分布图接口测试")
    parser.add_argument("-m", "--marker_type", type=str, help="Gene, SNP+InDel',SNP, InDel, CNV, SV, SSR",
                        default='Gene')
    parser.add_argument("-d", "--data_type", type=str, help="before, after, ref", default='ref')
    parser.add_argument("-w", "--win_step", type=int, help="window size", default=50000000)
    parser.add_argument("-a", "--analysis_object", type=str, help="all, sample_id")
    parser.add_argument("-g", "--graphic_style", type=str, help="heatmap, area", default='heatmap')
    parser.add_argument("-t", "--task_id", type=str, help="task_id", default='wgs_v2')
    parser.add_argument("-s", "--subname", type=str, help="比较分析批量的子分析名字")
    args = parser.parse_args()
    cmd = "python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py post s/wgs_v2/chromosome_distribution " \
          "-c client03 -b http://192.168.12.101:9090 "
    cmd += '-n "{}" -d "{}" '\
        .format("marker_type;data_type;win_step;analysis_object;graphic_style;task_id;task_type;submit_location;"
                "chongmingming_result;subname", ';'.join([args.marker_type, args.data_type, str(args.win_step),
                                                          args.analysis_object, args.graphic_style, args.task_id, "1",
                                                          "chr_distributon", args.subname]))
    print cmd
    os.system(cmd)

