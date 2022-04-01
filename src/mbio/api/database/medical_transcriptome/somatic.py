# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# 20201102

from __future__ import division

import datetime
import json
import os
import unittest

import pandas as pd
from bson.objectid import ObjectId

from mbio.api.database.medical_transcriptome.api_base import ApiBase


class Somatic(ApiBase):
    def __init__(self, bind_object):
        super(Somatic, self).__init__(bind_object)
        self._project_type = 'medical_transcriptome'
        # self._db_name = Config().MONGODB + '_ref_rna_v2'

    def add_somatic_main(self, snp_anno=None, project_sn=None, task_id=None, main_id=None, params=None, method_type=None,
                         predict_stat= None,new_output=None):
        """
        导入SNP主表信息
        :param snp_anno: snp模块的结果文件夹，如果传就直接导入detail信息~/output_dir/
        :return:
        """
        if main_id is None:
            # prepare main table info
            name = "Somatic" + '_' + method_type + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                # project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='somatic analysis main table',
                params=params,
                status="start",
                version='v1'
            )
            snp_id = self.create_db_table('sg_somatic', [main_info])
        if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
            snp_id = ObjectId(main_id)

        if snp_anno:
            self.add_snp_detail(snp_anno, snp_id, new_output=new_output)
            predict_pdf = os.path.join(self.bind_object.sheet.output, '03Gene_structure_analysis','02Somatic_Analysis', 'predict')
            self.add_pathogenicity_stat(predict_stat,snp_id)
            self.update_db_record('sg_somatic', snp_id, status="end",graph_path = predict_pdf ,
                                  main_id=snp_id)
        return snp_id

    def add_snp_detail(self, snp_anno, snp_id, new_output=None):
        """
        导入SNP详情表的函数
        :param snp_anno: snp模块的结果文件夹，即~/output_dir/
        :param snp_id: SNP主表的ID
        :return:
        """
        data_anno = os.path.join(snp_anno, "data_anno_pre.xls")
        with open(data_anno, "r") as f:
            headers = f.readline().strip().split("\t")
            data_list = list()
            for row, line in enumerate(f):
                data = dict()
                line = line.strip().split('\t')
                if len(headers) == len(line):
                    for n in range(len(headers)):
                        try:
                            data.update({headers[n]: float(line[n])})
                        except:
                            data.update({headers[n]: line[n]})
                    data.update({"somatic_id": snp_id})
                    data_list.append(data)
                    if row != 0 and row % 100000 == 0:
                        try:
                            self.create_db_table('sg_somatic_detail', data_list)
                        except Exception as e:
                            self.bind_object.set_error("导入SNP详情表出错:%s", variables=(e), code="55600504")
                        else:
                            self.bind_object.logger.info("导入SNP详情表成功，已经导入了%s条记录" % str(row))
                            data_list = list()
            if data_list:
                try:
                    self.create_db_table('sg_somatic_detail', data_list)
                except Exception as e:
                    self.bind_object.set_error("导入SNP详情表出错:%s", variables=(e), code="55600505")
                else:
                    self.bind_object.logger.info("导入SNP详情表成功，已经导入了%s条记录" % str(row))
                    data_list = list()

        df_depth_pre = pd.read_table(os.path.join(snp_anno, "snp_depth_statistics.xls"))
        df_depth_pre = df_depth_pre.fillna('')
        df_depth_pre['somatic_id'] = snp_id
        depth_data = df_depth_pre.to_dict('records')
        for depth in depth_data:
            for key in depth:
                try:
                    depth[key] = float(depth[key])
                except:
                    pass
        df_depth_pre.rename(columns={"range_key": "Depth"}, inplace=True)
        df_depth_pre.drop(['somatic_id', 'stat_type'], axis=1, inplace=True)
        df_depth_pre_col = df_depth_pre.columns[0:-1].tolist()
        df_depth_pre_col.insert(0, df_depth_pre.columns[-1])
        df_depth = df_depth_pre[df_depth_pre_col]
        df_depth.to_csv(new_output + "/snp_depth_statistics.xls", sep="\t", header=True, index=False)


        type_data_pre = pd.read_table(os.path.join(snp_anno, "snp_transition_tranversion_statistics.xls"))
        type_data_pre['somatic_id'] = snp_id
        type_data_pre = type_data_pre.fillna('')
        type_data = type_data_pre.to_dict('records')
        for type in type_data:
            for key in type:
                try:
                    type[key] = float(type[key])
                except:
                    pass
        type_data_pre.rename(columns={"range_key": "type"}, inplace=True)
        type_data_pre.drop(['somatic_id', 'stat_type'], axis=1, inplace=True)
        type_data_pre_col = type_data_pre.columns[0:-1].tolist()
        type_data_pre_col.insert(0, type_data_pre.columns[-1])
        df_type_data = type_data_pre[type_data_pre_col]
        df_type_data.loc['col_sum'] = df_type_data.apply(lambda x: x.sum())
        df_type_data.loc['col_sum', 'type'] = "total"
        # df_type_data.loc['col_sum', 'stat_type'] = "type_stat"
        df_type_data.to_csv(new_output + "/snp_transition_tranversion_statistics.xls", sep="\t", header=True,
                            index=False)

        snp_pos_data_pre = pd.read_table(os.path.join(snp_anno, "snp_position_distribution.xls"))
        snp_pos_data_pre = snp_pos_data_pre.fillna('')
        snp_pos_data_pre['somatic_id'] = snp_id
        snp_pos_data = snp_pos_data_pre.to_dict('records')
        for pos in snp_pos_data:
            for key in pos:
                try:
                    pos[key] = float(pos[key])
                except:
                    pass
        snp_pos_data_pre.rename(columns={"range_key": "Region"}, inplace=True)
        snp_pos_data_pre.drop(['somatic_id', 'stat_type'], axis=1, inplace=True)
        snp_pos_data_pre_col = snp_pos_data_pre.columns[0:-1].tolist()
        snp_pos_data_pre_col.insert(0, snp_pos_data_pre.columns[-1])
        df_snp_pos_data = snp_pos_data_pre[snp_pos_data_pre_col]
        df_snp_pos_data.loc['col_sum'] = df_snp_pos_data.apply(lambda x: x.sum())
        df_snp_pos_data.loc['col_sum', 'Region'] = "total"
        # df_snp_pos_data.loc['col_sum', 'stat_type'] = "snp_distribution"
        df_snp_pos_data.to_csv(new_output + "/snp_position_distribution.xls", sep="\t", header=True, index=False)

        indel_pos_data_pre = pd.read_table(os.path.join(snp_anno, "indel_position_distribution.xls"))
        indel_pos_data_pre = indel_pos_data_pre.fillna('')
        indel_pos_data_pre['somatic_id'] = snp_id
        indel_pos_data = indel_pos_data_pre.to_dict('records')
        for indel in indel_pos_data:
            for key in indel:
                try:
                    indel[key] = float(indel[key])
                except:
                    pass
        indel_pos_data_pre.rename(columns={"range_key": "Region"}, inplace=True)
        indel_pos_data_pre.drop(['somatic_id', 'stat_type'], axis=1, inplace=True)
        indel_pos_data_pre_col = indel_pos_data_pre.columns[0:-1].tolist()
        indel_pos_data_pre_col.insert(0, indel_pos_data_pre.columns[-1])
        df_indel_pos_data = indel_pos_data_pre[indel_pos_data_pre_col]
        df_indel_pos_data.loc['col_sum'] = df_indel_pos_data.apply(lambda x: x.sum())
        df_indel_pos_data.loc['col_sum', 'Region'] = "total"
        # df_indel_pos_data.loc['col_sum', 'stat_type'] = "indel_distribution"
        df_indel_pos_data.to_csv(new_output + "/indel_position_distribution.xls", sep="\t", header=True, index=False)
        with open(os.path.join(snp_anno, "main_info")) as main_r:
            main_dict = dict()
            for line in main_r:
                line = line.strip().split('\t')
                if line:
                    main_dict.update({line[0]: line[1].split(';')})

        try:
            stat_collection = self.db["sg_somatic_stat"]
            stat_collection.insert_many(depth_data)
            # stat_collection.insert_many(freq_data)
            stat_collection.insert_many(type_data)
            stat_collection.insert_many(snp_pos_data)
            stat_collection.insert_many(indel_pos_data)
            main_collection = self.db["sg_somatic"]
            main_collection.update({"_id": ObjectId(snp_id)}, {"$set": main_dict})
        except Exception as e:
            self.bind_object.set_error("导入SNP统计信息出错:%s", variables=(e), code="55600506")
        else:
            self.bind_object.logger.info("导入SNP统计信息成功")
        with open(os.path.join(snp_anno, "snp_annotation_statistics.xls"), "r")as all_info, open(
                os.path.join(new_output, "indel_anno.xls"), "w") as ide, open(os.path.join(new_output, "snp_anno.xls"),
                                                                              "w") as sp:
            header = all_info.readline()
            ide.write(header)
            sp.write(header)
            for line in all_info.readlines():
                lines = line.split("\t")
                if lines[-1].startswith("snp"):
                    sp.write(line)
                else:
                    ide.write(line)

    # os.link(os.path.join(snp_anno, "snp_annotation_statistics.xls"), os.path.join(new_output, "snp_anno.xls"))
    def add_pathogenicity_stat(self, predict_stat, somatic_id):
        predict_stat_df = pd.read_table(predict_stat)
        predict_stat_df["somatic_id"] = somatic_id
        predict_dict = predict_stat_df.to_dict("r")
        self.create_db_table('sg_somatic_predict', predict_dict)


    def run(self):
        params = {"method_type":"sentieon","submit_location":"somatic","task_id":"medical_transcriptome","task_type":2}
        new_snp_rewrite = "/mnt/ilustre/users/sanger-dev/workspace/20201013/MedicalTranscriptome_gbnijnfcuslii6rieegpbdgpt5/CallSnpIndelSentieon/output/predeal"
        self.add_somatic_main(snp_anno=new_snp_rewrite, method_type="sentieon",
                          task_id="medical_transcriptome", project_sn="medical_transcriptome",params=params,predict_stat ="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/test1102/test1",
                          new_output="/mnt/ilustre/users/sanger-dev/workspace/20201013/MedicalTranscriptome_gbnijnfcuslii6rieegpbdgpt5/CallSnpIndelSentieon/output")


class TestFunction(unittest.TestCase):
    """
    测试导表函数
    """

    def test_mongo(test):
        from mbio.workflows.medical_transcriptome.medical_transcriptome_test_api import MedicalTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet

        data = {
            # "id": "denovo_rna_v2" + str(random.randint(1,10000)),
            "id": "medical_transcriptome",
            "project_sn": "medical_transcriptome",
            "type": "workflow",
            "name": "medical_transcriptome.medical_transcriptome_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = MedicalTranscriptomeTestApiWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("medical_transcriptome.somatic")
        wf.test_api.run()


if __name__ == '__main__':
    unittest.main()
