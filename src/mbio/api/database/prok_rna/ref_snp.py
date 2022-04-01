# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
from __future__ import division
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
from bson.objectid import ObjectId
# from cStringIO import StringIO
from collections import Counter
import datetime
import glob
from mbio.api.database.prok_rna.api_base import ApiBase
import json
import pandas as pd


class RefSnp(ApiBase):
    def __init__(self, bind_object):
        super(RefSnp, self).__init__(bind_object)
        self._project_type = 'prok_rna'
        #self._db_name = Config().MONGODB + '_prok_rna'

    def add_snp_main(self, snp_anno=None, project_sn=None, task_id=None, main_id=None,params=None, method_type=None, new_output=None, group=None):
        """
        导入SNP主表信息
        :param snp_anno: snp模块的结果文件夹，如果传就直接导入detail信息~/output_dir/
        :return:
        """
        if main_id is None:
            # prepare main table info
            name = "SNP" + '_' + method_type + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='snp_indel analysis main table',
                params=params,
                status="start",
            )
            snp_id = self.create_db_table('sg_snp', [main_info])
        if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
            snp_id = ObjectId(main_id)

        if snp_anno:
            self.add_snp_detail(snp_anno, snp_id, new_output=new_output, group=group)
            self.update_db_record('sg_snp', snp_id, status="end",
                                  main_id=snp_id)
        return snp_id

    def add_snp_detail(self, snp_anno, snp_id, new_output=None, group=None):
        """
        导入SNP详情表的函数
        :param snp_anno: snp模块的结果文件夹，即~/output_dir/
        :param snp_id: SNP主表的ID
        :return:
        """
        if group:
            sample_list_ = list()
            with open(group, 'r') as g:
                for line in g.readlines():
                    if line.startswith('#'):
                        continue
                    sample_list_.append(line.strip().split('\t')[0])
        snp_anno = glob.glob('{}/*'.format(snp_anno))[0]
        snp_type_stat = {}
        snp_pos_stat = {}
        indel_pos_stat = {}
        all_depth_stat = {}
        all_freq_stat = {}
        depth_list = ["<=30", "31-100", "101-200", "201-300", "301-400", "401-500", ">500"]
        # graph_data_list = []
        chroms = set()
        distributions = set()
        data_list = []
        sample_names = []
        sample_old_index = []
        sample_old_gene = []

        with open(snp_anno, "r") as f:
            sample_names = f.readline().strip().split("\t")[11:]
            for s in sample_names:
                snp_type_stat[s] = {}
                snp_pos_stat[s] = {}
                indel_pos_stat[s] = {}
                # all_freq_stat[s] = [0]
                all_freq_stat[s] = []   # modified by zhangyitong on 20211021
                all_depth_stat[s] = [0, 0, 0, 0, 0, 0, 0]
                sample_old_index.append(-1)
                sample_old_gene.append('')
            # print f.next()
            for line in f:
                line = line.strip().split("\t")
                sample_infos = line[11:]
                new_gene_name = line[7]
                chroms.add(line[0])
                distributions.add(line[6])
                snp_type = line[3] + "/" + line[4]
                data = {
                    "snp_id": snp_id,
                    "type": "snp" if len(line[3]) + len(line[4]) == 2 and "-" not in snp_type else "indel",
                    "chrom": line[0],
                    "start": line[1],
                    "end": line[2],
                    "ref": line[3],
                    "alt": line[4],
                    "qual": float(line[5]),
                    "reads_num": int(line[8]),
                    # "mut_rate": 0.33,
                    "anno": line[6],
                    "gene": line[7],
                    "mut_type": line[9],
                    "mut_info": line[10],
                    "snp_type": snp_type
                }
                for n, s in enumerate(sample_names):
                    rate = sample_infos[n]
                    # print sample_infos[n]
                    mut_rate = 0
                    depth_num = -1
                    single_and_all = '.'
                    if rate != './.' and rate != '0/0':
                        single_and_all = rate.split("/")[1]
                        mut_rate = round(int(rate.split("/")[0])/int(rate.split("/")[1]), 4)
                        # 统计各样本的突变频率数目
                        if line[6] == "exonic":
                            sample_old_index[n], all_freq_stat[s] = self.gene_num_stat(new_gene_name, sample_old_gene[n], sample_old_index[n], all_freq_stat[s])
                            sample_old_gene[n] = new_gene_name
                        #  统计各样本的SNP类型数目
                        if not '-' in snp_type and len(snp_type) == 3:
                            snp_type_stat[s] = self.type_stat(snp_type, snp_type_stat[s])
                        depth_num = int(rate.split("/")[1])

                        # 统计SNP/Indel位置信息
                        if data["type"] == "snp":
                            snp_pos_stat[s] = self.type_stat(data["anno"], snp_pos_stat[s])
                        else:
                            indel_pos_stat[s] = self.type_stat(data["anno"], indel_pos_stat[s])
                    data[s + "_mut_rate"] = mut_rate
                    data[s + "_reads_rate"] = single_and_all
                    all_depth_stat[s] = self.get_depth_stat(depth_num, all_depth_stat[s])
                data_list.append(data)
                # #统计reads深度
        df_data_list_pre = pd.DataFrame(data_list)
        df_data_list_pre.to_csv(new_output + "/data_anno_pre.xls", sep="\t", header=True, index=False)
        # snp_types = ["A/G", "A/C", "C/T", "G/A", "G/C", "C/A", "A/T", "C/G", "G/T", "T/C", "T/A", "T/G"]

        # new_freq_stat = self.freq_stat(all_freq_stat)
        # depth_data = self.get_stat_data(sample_names, depth_list, all_depth_stat, snp_id, "depth_stat")
        # freq_data = self.get_stat_data(sample_names, [1, 2, 3, 4, ">=5"], new_freq_stat, snp_id, "freq_stat")
        # type_data = self.get_stat_data(sample_names, snp_types, snp_type_stat, snp_id, "type_stat")
        # snp_pos_data = self.get_stat_data(sample_names, list(distributions), snp_pos_stat, snp_id, "snp_distribution")
        # indel_pos_data = self.get_stat_data(sample_names, list(distributions), indel_pos_stat, snp_id, "indel_distribution")
        depth_data = self.get_stat_data(sample_names, depth_list, all_depth_stat, snp_id, "depth_stat")
        df_depth_pre = pd.DataFrame(depth_data)
        df_depth_pre.rename(columns={"range_key": "Depth"}, inplace=True)
        df_depth_pre.drop(['snp_id', 'stat_type'], axis=1, inplace=True)
        df_depth_pre_col = df_depth_pre.columns[0:-1].tolist()
        df_depth_pre_col.insert(0, df_depth_pre.columns[-1])
        df_depth = df_depth_pre[df_depth_pre_col]
        df_depth.to_csv(new_output + "/snp_depth_statistics.xls", sep="\t", header=True, index=False)

        new_freq_stat = self.freq_stat(all_freq_stat)
        freq_data = self.get_stat_data(sample_names, [1, 2, 3, 4, ">=5"], new_freq_stat, snp_id, "freq_stat")
        freq_data_pre = pd.DataFrame(freq_data)
        freq_data_pre.rename(columns={"range_key": "Freq"}, inplace=True)
        freq_data_pre.drop(['snp_id', 'stat_type'], axis=1, inplace=True)
        freq_data_pre_col = freq_data_pre.columns[0:-1].tolist()
        freq_data_pre_col.insert(0, freq_data_pre.columns[-1])
        df_freq_data = freq_data_pre[freq_data_pre_col]
        df_freq_data.to_csv(new_output + "/snp_freq_statistics.xls", sep="\t", header=True, index=False)

        snp_types = ["A/G", "A/C", "C/T", "G/A", "G/C", "C/A", "A/T", "C/G", "G/T", "T/C", "T/A", "T/G"]
        type_data = self.get_stat_data(sample_names, snp_types, snp_type_stat, snp_id, "type_stat")
        type_data_pre = pd.DataFrame(type_data)
        type_data_pre.rename(columns={"range_key": "type"}, inplace=True)
        type_data_pre.drop(['snp_id', 'stat_type'], axis=1, inplace=True)
        type_data_pre_col = type_data_pre.columns[0:-1].tolist()
        type_data_pre_col.insert(0, type_data_pre.columns[-1])
        df_type_data = type_data_pre[type_data_pre_col]
        df_type_data.to_csv(new_output + "/snp_transition_tranversion_statistics.xls", sep="\t", header=True, index=False)

        snp_pos_data = self.get_stat_data(sample_names, list(distributions), snp_pos_stat, snp_id, "snp_distribution")
        snp_pos_data_pre = pd.DataFrame(snp_pos_data)
        print(snp_pos_data_pre.head(3))
        print(99999999999)
        snp_pos_data_pre.rename(columns={"range_key": "Region"}, inplace=True)
        snp_pos_data_pre.drop(['snp_id', 'stat_type'], axis=1, inplace=True)
        snp_pos_data_pre_col = snp_pos_data_pre.columns[0:-1].tolist()
        snp_pos_data_pre_col.insert(0, snp_pos_data_pre.columns[-1])
        df_snp_pos_data = snp_pos_data_pre[snp_pos_data_pre_col]
        df_snp_pos_data.to_csv(new_output + "/snp_position_distribution.xls", sep="\t", header=True, index=False)

        indel_pos_data = self.get_stat_data(sample_names, list(distributions), indel_pos_stat, snp_id, "indel_distribution")
        indel_pos_data_pre = pd.DataFrame(indel_pos_data)
        indel_pos_data_pre.rename(columns={"range_key": "Region"}, inplace=True)
        indel_pos_data_pre.drop(['snp_id', 'stat_type'], axis=1, inplace=True)
        indel_pos_data_pre_col = indel_pos_data_pre.columns[0:-1].tolist()
        indel_pos_data_pre_col.insert(0, indel_pos_data_pre.columns[-1])
        df_indel_pos_data = indel_pos_data_pre[indel_pos_data_pre_col]
        df_indel_pos_data.to_csv(new_output + "/indel_position_distribution.xls", sep="\t", header=True, index=False)

        try:
            collection = self.db["sg_snp_detail"]
            collection.insert_many(data_list)
            # graph_collection = self.db["sg_snp_graphic"]
            # graph_collection.insert_many(graph_data_list)
            stat_collection = self.db["sg_snp_stat"]
            stat_collection.insert_many(depth_data)
            stat_collection.insert_many(freq_data)
            stat_collection.insert_many(type_data)
            stat_collection.insert_many(snp_pos_data)
            stat_collection.insert_many(indel_pos_data)
            main_collection = self.db["sg_snp"]
            main_collection.update({"_id": ObjectId(snp_id)}, {"$set": {"snp_types": snp_types, "specimen": sample_list_, "distributions": list(distributions), "chroms": list(chroms)}})
        except Exception as e:
            print("导入SNP统计信息出错:%s" % e)
        else:
            print("导入SNP统计信息成功")

    def get_depth_stat(self, depth_num, target_list):
        if depth_num == -1:
            pass
        else:
            if depth_num < 31:
                target_list[0] += 1
            elif 30 < depth_num < 101:
                target_list[1] += 1
            elif 100 < depth_num < 201:
                target_list[2] += 1
            elif 200 < depth_num < 301:
                target_list[3] += 1
            elif 300 < depth_num < 401:
                target_list[4] += 1
            elif 400 < depth_num < 501:
                target_list[5] += 1
            else:
                target_list[6] += 1
        # print target_dict
        return target_list

    def type_stat(self, dict_key, target_dict):
        if dict_key in target_dict:
            target_dict[dict_key] += 1
        else:
            target_dict[dict_key] = 1
        return target_dict

    def get_stat_data(self, sample_names, range_list, value_dict, snp_id, stat_type):
        stat_data = []
        for n, ds in enumerate(range_list):
            data = {
                "snp_id": snp_id,
                "range_key": ds,
                "stat_type": stat_type
            }
            for s in sample_names:
                if type(value_dict[s]) is list:
                    data["{}".format(s)] = value_dict[s][n]
                else:
                    if ds not in value_dict[s].keys():
                        value_dict[s][ds] = 0
                    data["{}".format(s)] = value_dict[s][ds]
            stat_data.append(data)
        # print stat_data
        return stat_data

    def gene_num_stat(self, new_gene_name, old_gene_name, old_index, new_list):
        if new_gene_name == old_gene_name:
            new_list[old_index] += 1
        else:
            new_list.append(1)
            old_index += 1
        return old_index, new_list

    def freq_stat(self, all_freq_stat):
        for a_s in all_freq_stat:
            c = Counter(all_freq_stat[a_s])
            values = c.values()
            keys = c.keys()
            new_d = {">=5": 0}
            for n, s in enumerate(keys):
                # print s
                if s < 5:
                    new_d[s] = values[n]
                else:
                    new_d[">=5"] += values[n]
            all_freq_stat[a_s] = new_d
        return all_freq_stat


if __name__ == '__main__':
     anno = RefSnp(None)
     new_snp_rewrite = "/mnt/ilustre/users/sanger-dev/sg-users/litangjian/merge"
     new_output = "/mnt/ilustre/users/sanger-dev/sg-users/litangjian/merge"
     anno.add_snp_main(snp_anno=new_snp_rewrite, method_type="gatk",
                       task_id="test1", project_sn="test1", new_output=new_output)
