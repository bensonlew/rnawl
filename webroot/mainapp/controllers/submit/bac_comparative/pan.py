# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang' @20191122
import web,re
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bac_comparative_controller import BacComparativeController
from mainapp.libs.signature import check_sig
from bson import SON
from mainapp.libs.param_pack import group_detail_sort


class PanAction(BacComparativeController):
    def __init__(self):
        super(PanAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = ['task_id', 'gene_dir', 'submit_location', 'task_type','group_detail', 'group_id', 'cluster']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        task_name = 'bac_comp_genome.report.pan'
        module_type = 'workflow'  # 可以不配置
        project_sn = self.bac_comparative.get_projectsn(data.task_id)
        if data.cluster not in ['usearch', 'cd_hit', 'mmseq', 'blast+mcl', 'cdhit+blast+mcl']:
            info = {'success': False, 'info': '%s参数缺少!' % data.cluster}
            return json.dumps(info)
        ###放入params的参数，需要根据传入的参数进行判断
        params = {
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'task_id': data.task_id
        }
        all_sample_list = []
        if data.group_detail:
            group_detail = json.loads(data.group_detail)
            for key in group_detail:
                group_string = group_detail[key]
                # group_list = group_string.split(",")
                for x in group_string:
                    if x not in all_sample_list:
                        all_sample_list.append(x)
        if len(all_sample_list) == 1:
            info = {'success': False, 'info': '样本数不能为1，请重新选择样本进行分析!'}
            return json.dumps(info)

        if data.cluster in ['usearch']:
            params["cluster"] = data.cluster
            if not hasattr(data, "identity"):
                info = {'success': False, 'info': 'parameters missing:%s' % data.identity}
                return json.dumps(info)
            else:
                params["identity"] = data.identity

        if data.cluster in ['cd_hit','mmseq']:
            params["cluster"] = data.cluster
            if not hasattr(data, "identity"):
                info = {'success': False, 'info': 'parameters missing:%s' % data.identity}
                return json.dumps(info)
            else:
                params["identity"] = data.identity
            if not hasattr(data, "coverage"):
                info = {'success': False, 'info': 'parameters missing:%s' % data.coverage}
                return json.dumps(info)
            else:
                params["coverage"] = data.coverage

        if data.cluster in ['blast+mcl']:
            params["cluster"] = data.cluster
            if not hasattr(data, "software"):
                info = {'success': False, 'info': 'parameters missing:%s' % data.identity}
                return json.dumps(info)
            else:
                if data.software not in ['PGAP', 'Orthofinder', 'Orthomcl', 'GET_HOMOLOGUS']:
                    info = {'success': False, 'info': '%s参数缺少!' % data.software}
                    return json.dumps(info)
                else:
                    params["software"] = data.software
                    if not hasattr(data, "coverage"):
                        info = {'success': False, 'info': 'parameters missing:%s' % data.coverage}
                        return json.dumps(info)
                    else:
                        params["coverage"] = data.coverage
                    if not hasattr(data, "identity"):
                        info = {'success': False, 'info': 'parameters missing:%s' % data.identity}
                        return json.dumps(info)
                    else:
                        params["identity"] = data.identity
                    if not hasattr(data, "inflation"):
                        info = {'success': False, 'info': 'parameters missing:%s' % data.inflation}
                        return json.dumps(info)
                    else:
                        params["inflation"] = data.inflation
                    if data.software in ['PGAP']:
                        if not hasattr(data, "method"):
                            info = {'success': False, 'info': 'parameters missing:%s' % data.method}
                            return json.dumps(info)
                        else:
                            params["method"] = data.method
                            if not data.method in ['GF', 'MP']:
                                info = {'success': False, 'info': 'mehtod not in ["GF", "MP"]:%s' % data.method}
                                return json.dumps(info)
                            else:
                                if data.method in ['GF']:
                                    if not hasattr(data, "evalue"):
                                        info = {'success': False, 'info': 'parameters missing:%s' % data.evalue}
                                        return json.dumps(info)
                                    else:
                                        params["evalue"] = data.evalue
                                    if not hasattr(data, "score"):
                                        info = {'success': False, 'info': 'parameters missing:%s' % data.score}
                                        return json.dumps(info)
                                    else:
                                        params["score"] = data.score
                                elif data.method in ['MP']:
                                    if not hasattr(data, "local"):
                                        info = {'success': False, 'info': 'parameters missing:%s' % data.local}
                                        return json.dumps(info)
                                    else:
                                        params["local"] = data.local
                                    if not hasattr(data, "globa"):
                                        info = {'success': False, 'info': 'parameters missing:%s' % data.globa}
                                        return json.dumps(info)
                                    else:
                                        params["globa"] = data.globa

                    elif data.software in ['Orthofinder']:
                        if not hasattr(data, "coverage"):
                            info = {'success': False, 'info': 'parameters missing:%s' % data.coverage}
                            return json.dumps(info)
                        if not hasattr(data, "identity"):
                            info = {'success': False, 'info': 'parameters missing:%s' % data.identity}
                            return json.dumps(info)
                        if not hasattr(data, "inflation"):
                            info = {'success': False, 'info': 'parameters missing:%s' % data.inflation}
                            return json.dumps(info)
                    elif data.software in ['Orthomcl']:
                        if not hasattr(data, "evalue"):
                            info = {'success': False, 'info': 'parameters missing:%s' % data.evalue}
                            return json.dumps(info)
                        else:
                            params["evalue"] = data.evalue

                    elif data.software in ['GET_HOMOLOGUS']:
                        if not hasattr(data, "method"):
                            info = {'success': False, 'info': 'parameters missing:%s' % data.method}
                            return json.dumps(info)
                        else:
                            if not data.method in ['BDBH', 'OMCL', 'OCOG']:
                                info = {'success': False, 'info': 'mehtod not in ["BDBH", "OMCL", "OCOG"]:%s' % data.method}
                                return json.dumps(info)
                            else:
                                params["method"] = data.method
                            if not hasattr(data, "identity"):
                                info = {'success': False, 'info': 'parameters missing:%s' % data.identity}
                                return json.dumps(info)
                            else:
                                params["identity"] = data.identity
                            if not hasattr(data, "coverage"):
                                info = {'success': False, 'info': 'parameters missing:%s' % data.coverage}
                                return json.dumps(info)
                            else:
                                params["coverage"] = data.coverage
                            if not hasattr(data, "inflation"):
                                info = {'success': False, 'info': 'parameters missing:%s' % data.inflation}
                                return json.dumps(info)
                            else:
                                params["inflation"] = data.inflation
                    else:
                        info = {'success': False, 'info': 'parameters missing:%s' % data.software}
                        return json.dumps(info)

        if data.cluster in ['cdhit+blast+mcl']:
            params["cluster"] = data.cluster
            if not hasattr(data, "coverage"):
                info = {'success': False, 'info': 'parameters missing:%s' % data.coverage}
                return json.dumps(info)
            else:
                params["coverage"] = data.coverage
            if not hasattr(data, "identity"):
                info = {'success': False, 'info': 'parameters missing:%s' % data.identity}
                return json.dumps(info)
            else:
                params["identity"] = data.identity
            if not hasattr(data, "inflation"):
                info = {'success': False, 'info': 'parameters missing:%s' % data.inflation}
                return json.dumps(info)
            else:
                params["inflation"] = data.inflation

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_table_name = 'Pan_'+ datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        if hasattr(data, "summary_path"):
            anno_file = data.summary_path
        else:
            anno_file = self.bac_comparative.get_annofile(data.task_id)
        update_info = {}
        mongo_data = [
            ('project_sn', project_sn),
            ('status', 'start'),
            ('task_id', data.task_id),
            ('desc', '正在计算'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", params)
        ]
        main_id = self.bac_comparative.insert_main_table('pan', mongo_data)
        update_info[str(main_id)] = 'pan'
        options = {
                   'update_info': json.dumps(update_info),
                   'params': params,
                   'main_id': str(main_id),
                   'group': data.group_id,
                   "group_detail": data.group_detail,
                   'fasta_dir': data.gene_dir,
                   'anno_file': anno_file,
                   }
        to_file = "bac_comparative.export_group_table_by_detail(group)"

        if hasattr(data, "cluster"): ###选择的聚类方式
            if data.cluster in ["cd_hit"]:
                options["cluster_method"] = 'cdhit'
            else:
                options["cluster_method"] = data.cluster

        if hasattr(data, "software"): ###选择的聚类软件
            options["method"] = (data.software).lower()
            if data.software in ["PGAP"]:
                if hasattr(data, "method"): ###选择的聚类方法
                    options["pgap_mehod"] = (data.method)
            elif data.software in ["GET_HOMOLOGUS"]:
                if hasattr(data, "method"): ###选择的聚类方法
                    options["homologus_mehod"] = (data.method).lower()

        if hasattr(data, "coverage"): ###选择的比对的coverage
            options["coverage"] = float(data.coverage) / 100
            if (int(float(data.coverage)) < 0) and (int(float(data.coverage)) > 100):
                info = {'success': False, 'info': 'coverage 必须在[0,100]中'}
                return json.dumps(info)

        if hasattr(data, "identity"): ###选择的比对的identity
            options["identity"] = float(data.identity) / 100
            if (int(float(data.identity)) < 20) and (int(float(data.identity)) > 100):
                info = {'success': False, 'info': 'identity 必须在[20,100]中'}
                return json.dumps(info)

        if hasattr(data, "inflation"):###选择的聚类膨胀系数
            options["inflation"] = float(data.inflation)
            if (int(float(data.inflation)) < 1) and (int(float(data.inflation)) > 5):
                info = {'success': False, 'info': 'inflation 必须在[1,5]中'}
                return json.dumps(info)

        if hasattr(data, "evalue"):###选择的evalue值
            options["evalue"] = data.evalue

        if hasattr(data, "score"):###选择的得分系数
            options["score"] = int(float(data.score))
            if (int(float(data.score)) < 10):
                info = {'success': False, 'info': 'score 必须在大于等于10中'}
                return json.dumps(info)

        if hasattr(data, "local"): ###选择的local系数
            options["local"] = float(data.local)

        if hasattr(data, "globa"):###选择的global系数
            options["globa1"] = float(data.globa)

        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="PAN/" + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params,
                            to_file=to_file)
        task_info = super(PanAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)