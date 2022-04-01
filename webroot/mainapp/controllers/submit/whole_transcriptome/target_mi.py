# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import os
import unittest

import web

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.param_pack import *
from mainapp.libs.signature import check_sig


class TargetMiAction(WholeTranscriptomeController):
    """
    whole_transcriptome 靶基因预测与注释重运行接口
    """

    def __init__(self):
        super(TargetMiAction, self).__init__(instant=False)

    def check_params(self, task_info, data):
        spe_tax = {
            "Mus_musculus": 10090,
            "Rattus_norvegicus": 10116,
            "Monodelphis_domestica": 13616,
            "Xenopus_tropicalis": 8364,
            "Gallus_gallus": 9031,
            "Macaca_mulatta": 9544,
            "Pan_troglodytes": 9598,
            "Homo_sapiens": 9606,
            "Canis_lupus_familiaris": 9615,
            "Bos_taurus": 9913
        }

        if not self.genome_info['organism_name'] in spe_tax and hasattr(data, 'targetscan') and data[
            'targetscan'] == "yes":
            info = {"success": False, "info": "该项目物种 %s 不支持targetscan" % self.genome_info['organism_name']}
            return json.dumps(info)

    @check_sig
    def POST(self):
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get("HTTP_CLIENT")
        print  "data is {}".format(data)
        return_result = self.check_options(data)

        if return_result:
            variables = []
            variables.append(return_result)
            info = {"success": False, "info": '%s' % return_result}
            return json.dumps(info)

        task_info = self.whole_transcriptome.get_main_info_by_record("sg_task", task_id=str(data.task_id))

        run_info = self.whole_transcriptome.get_main_info_by_record("target_mi", task_id=str(data.task_id),
                                                                    status="start", type="latest")
        '''
        if run_info:
            variables = []
            variables.append(run_info['name'])
            info = {"success": False, "info": "任务%s正在计算，请等待其完成在下次运行" % run_info['name']}
            return json.dumps(info)
        '''

        task_info = self.whole_transcriptome.get_task_info(data.task_id)
        self.genome_info = self.whole_transcriptome.get_genome_info_by_genome_id(task_info['genome_id'])
        self.check_params(task_info, data)

        latest_info_target = self.whole_transcriptome.get_main_info_by_record("target_mi", task_id=str(data.task_id),
                                                                              status="end", type="latest")
        latest_main_id_target = ""

        try:
            target_info = self.whole_transcriptome.get_main_info_by_record("target_mi", task_id=str(data.task_id),
                                                                           status="end", type="origin")

            result = target_info["result_dir"]
        except:
            result = ""

        '''
        if 'target' in task_info and os.path.exists(task_info["target_mi"]):
            target = task_info["target_mi"]
        else:
            target = result + "/" + "target.fa"
        '''
        if latest_info_target:
            latest_main_id_target = latest_info_target["_id"]

        params_json = {
            "task_id": str(data.task_id),
            "submit_location": data.submit_location,
            "task_type": int(data.task_type)
        }
        for method in [
            'm_miranda', 'm_targetscan', 'm_psrobot', 'm_targetfinder', 'm_rnahybrid',
            'l_miranda', 'l_targetscan', 'l_psrobot', 'l_targetfinder', 'l_rnahybrid',
            'c_miranda', 'c_targetscan', 'c_psrobot', 'c_targetfinder', 'c_rnahybrid']:
            if hasattr(data, method):
                params_json.update({
                    method: data[method]
                })
        # 参数检查
        for category in ['mRNA', 'lncRNA', 'circRNA']:
            m_list = list()
            abr = category[0]
            for method in [
                abr + '_miranda', abr + '_targetscan', abr + '_psrobot', abr + '_targetfinder', abr + '_rnahybrid']:

                if hasattr(data, method):
                    if data[method] == "yes":
                        m_list.append(method)
            if int(data[abr + "_min_support"]) > len(m_list):
                info = {"success": False, "info": "靶{}最少支持的数量多于选择的预测方法".format(category)}
                return json.dumps(info)





        for par in ['m_miranda_score', 'm_miranda_energy', 'm_miranda_strict',
                    'm_rnahybird_num', 'm_rnahybird_energy', 'm_rnahybird_pvalue',
                    'm_ps_robot_score', 'm_targetfinder_score',
                    'l_miranda_score', 'l_miranda_energy', 'l_miranda_strict',
                    'l_rnahybird_num', 'l_rnahybird_energy', 'l_rnahybird_pvalue',
                    'l_ps_robot_score', 'l_targetfinder_score',
                    'c_miranda_score', 'c_miranda_energy', 'c_miranda_strict',
                    'c_rnahybird_num', 'c_rnahybird_energy', 'c_rnahybird_pvalue',
                    'c_ps_robot_score', 'c_targetfinder_score'
                    ]:
            if hasattr(data, par):
                params_json.update({
                    par: data[par]
                })


        for par in ['m_min_support', 'l_min_support', 'c_min_support']:
            if hasattr(data, par):
                params_json.update({
                    par: data[par]
                })

        for par in ['miset_id', 'targetset_id']:
            if hasattr(data, par):
                params_json.update({
                    par: data[par]
                })
            else:
                setattr(data, par, "All")

        main_table_name2 = "SmallTarget_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))

        mongo_data2 = [
            ("project_sn", task_info["project_sn"]),
            ("task_id", task_info["task_id"]),
            ("status", "start"),
            ("name", main_table_name2),
            ("desc", "靶基因预测主表"),
            ("result_dir", ""),
            ("taxonomy", task_info['options']['mirbase_category']),
            ("type", "latest"),
            ("version", "v1"),
            ("species_name", self.genome_info['organism_name']),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]

        main_table_id2 = self.whole_transcriptome.insert_main_table("target_mi", mongo_data2)
        update_info = {str(main_table_id2): "target_mi"}

        options = {
            "novol": data.task_id + ";" + "miRNA;new" + ";" + data.miset_id,
            "known": data.task_id + ";" + "miRNA;ref" + ";" + data.miset_id,
            "m_novol": data.task_id + ";" + "mRNA;new" + ";" + data.targetset_id,
            "m_known": data.task_id + ";" + "mRNA;ref" + ";" + data.targetset_id,
            "task_id": str(data.task_id),
            "last_id_target": str(latest_main_id_target),
            "target_id": str(main_table_id2),
            "species_name": self.genome_info['organism_name'],
            "taxonomy": task_info['options']['mirbase_category'],
            "update_info": json.dumps(update_info),
            "m_min_support": data.m_min_support,
            "l_min_support": data.l_min_support,
            "c_min_support": data.c_min_support,
            "circ_detail": data.task_id,
            "anno_detail": data.task_id

        }

        if 'circRNA' in task_info['rna']:
            options.update({
                "c_novol": data.task_id + ";" + "circRNA;new" + ";" + data.targetset_id,
                "c_known": data.task_id + ";" + "circRNA;ref" + ";" + data.targetset_id
            })

        if 'lncRNA' in task_info['rna']:
            options.update({
                "l_novol": data.task_id + ";" + "lncRNA;new" + ";" + data.targetset_id,
                "l_known": data.task_id + ";" + "lncRNA;ref" + ";" + data.targetset_id,
                "lnc_ref_detail": task_info["output"] + '/lncrna/01_lncRNA_Analysis/01_Known_lncRNA/known_lncrna_detail.xls',
                "lnc_new_detail": task_info["output"] + '/lncrna/01_lncRNA_Analysis/02_Novel_lncRNA/novel_lncrna_predict_detail.xls'
            })

        for method in [
            'm_miranda', 'm_targetscan', 'm_psrobot', 'm_targetfinder', 'm_rnahybrid',
            'l_miranda', 'l_targetscan', 'l_psrobot', 'l_targetfinder', 'l_rnahybrid',
            'c_miranda', 'c_targetscan', 'c_psrobot', 'c_targetfinder', 'c_rnahybrid']:
            if hasattr(data, method):
                options.update({
                    method: data[method]
                })

        for par in ['m_miranda_score', 'm_miranda_energy', 'm_miranda_strict',
                    'm_rnahybird_num', 'm_rnahybird_energy', 'm_rnahybird_pvalue',
                    'm_ps_robot_score', 'm_targetfinder_score',
                    'l_miranda_score', 'l_miranda_energy', 'l_miranda_strict',
                    'l_rnahybird_num', 'l_rnahybird_energy', 'l_rnahybird_pvalue',
                    'l_ps_robot_score', 'l_targetfinder_score',
                    'c_miranda_score', 'c_miranda_energy', 'c_miranda_strict',
                    'c_rnahybird_num', 'c_rnahybird_energy', 'c_rnahybird_pvalue',
                    'c_ps_robot_score', 'c_targetfinder_score'
                    ]:
            if hasattr(data, par):
                options.update({
                    par: data[par].lstrip("\\")
                })

        to_file = [
            "whole_transcriptome.advance.export_mrna_detail(anno_detail)",
            "whole_transcriptome.advance.get_seq(novol)",
            "whole_transcriptome.advance.get_seq(known)",
            "whole_transcriptome.advance.get_seq(m_novol)",
            "whole_transcriptome.advance.get_seq(m_known)",

        ]

        if 'circRNA' in task_info['rna']:
            to_file += [
                "whole_transcriptome.advance.export_circ_detail(circ_detail)",
                "whole_transcriptome.advance.get_seq(c_novol)",
                "whole_transcriptome.advance.get_seq(c_known)"
            ]
        if 'lncRNA' in task_info['rna']:
            to_file += [
                "whole_transcriptome.advance.get_seq(l_novol)",
                "whole_transcriptome.advance.get_seq(l_known)"
            ]

        self.set_sheet_data(name="whole_transcriptome.report.target_mi",
                            options=options,
                            main_table_name=main_table_name2,
                            task_id=task_info["task_id"],
                            project_sn=task_info["project_sn"],
                            to_file=to_file, )
        task_info = super(TargetMiAction, self).POST()
        task_info["content"] = {"ids": {"id": str(main_table_id2), "name": main_table_name2}}
        return json.dumps(task_info)

    def check_options(self, data):
        """
        检查网页端传来的参数是否正确
        """

        if len(data.miset_id) == 24 or data.miset_id.lower() == "all":
            pass
        else:
            return "请选择mirna基因集"

        if len(data.targetset_id) == 24 or data.targetset_id.lower() == "all":
            pass
        else:
            return "请选择target基因集"
        warn_info = list()
        method_num = 0
        return False


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/whole_transcriptome/target_mi "
        cmd += "-b http://bcl.tsg.com "
        args = dict({
            "task_id": "tsg_36088",
            "m_miranda": "yes",
            "c_miranda": "yes",
            "l_miranda": "yes",
            "m_miranda_score": "170",
            "m_miranda_energy": '\-20',
            "m_miranda_strict": "on",
            "l_miranda_score": "160",
            "l_miranda_energy": '\-20',
            "l_miranda_strict": "on",
            "c_miranda_score": "165",
            "c_miranda_energy": '\-20',
            "c_miranda_strict": "on",
            "m_min_support": 1,
            "l_min_support": 1,
            "c_min_support": 1,
            "submit_location": "target_mi",
            "task_type": 2,
        })
        arg_names, arg_values = args.keys(), args.values()
        arg_values = [str(x) for x in arg_values]
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
