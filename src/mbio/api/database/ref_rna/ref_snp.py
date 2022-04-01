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


class RefSnp(Base):
    def __init__(self, bind_object):
        super(RefSnp, self).__init__(bind_object)
        self._project_type = 'ref_rna'
        #self._db_name = Config().MONGODB + '_ref_rna'

    def add_snp_main(self, snp_anno=None):
        """
        导入SNP主表信息
        :param snp_anno: snp模块的结果文件夹，如果传就直接导入detail信息~/output_dir/
        :return:
        """
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "name": "snp_origin",
            "status": "end",
            "desc": "SNP主表",
            "snp_type": ["A/T", "A/C", "A/G", "C/G", "C/T", "G/C"],
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_snp"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        if snp_anno:
            self.add_snp_detail(snp_anno, inserted_id)
        # add_coverage_detail(coverage, inserted_id)
        return inserted_id

    def add_snp_detail(self, snp_anno, snp_id):
        """
        导入SNP详情表的函数
        :param snp_anno: snp模块的结果文件夹，即~/output_dir/
        :param snp_id: SNP主表的ID
        :return:
        """
        snp_anno = glob.glob('{}/*'.format(snp_anno))[0]
        snp_type_stat = {}
        snp_pos_stat = {}
        indel_pos_stat = {}
        all_depth_stat = {}
        all_freq_stat = {}
        depth_list = ["<=30", "31-100", "101-200", "201-300", "301-400", "401-500", ">501"]
        # graph_data_list = []
        chroms = set()
        distributions = set()
        data_list = []
        sample_names = []
        sample_old_index = []
        sample_old_gene = []

        with open(snp_anno, "r") as f:
            sample_names = f.readline().strip().split("\t")[10:]
            for s in sample_names:
                snp_type_stat[s] = {}
                snp_pos_stat[s] = {}
                indel_pos_stat[s] = {}
                all_freq_stat[s] = [0]
                all_depth_stat[s] = [0, 0, 0, 0, 0, 0, 0]
                sample_old_index.append(-1)
                sample_old_gene.append('')
            # print f.next()
            for line in f:
                line = line.strip().split("\t")
                sample_infos = line[10:]
                new_gene_name = line[6]
                chroms.add(line[0])
                distributions.add(line[5])
                snp_type = line[3] + "/" + line[4]
                data = {
                    "snp_id": snp_id,
                    "type": "snp" if len(line[3]) + len(line[4]) == 2 and "-" not in snp_type else "indel",
                    "chrom": line[0],
                    "start": line[1],
                    "end": line[2],
                    "ref": line[3],
                    "alt": line[4],
                    "reads_num": int(line[7]),
                    # "mut_rate": 0.33,
                    "anno": line[5],
                    "gene": line[6],
                    "mut_type": line[8],
                    "mut_info": line[9],
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
                        if line[5] == "exonic":
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
        snp_types = ["A/G", "A/C", "C/T", "G/A", "G/C", "C/A", "A/T", "C/G", "G/T", "T/C", "T/A", "T/G"]

        new_freq_stat = self.freq_stat(all_freq_stat)
        depth_data = self.get_stat_data(sample_names, depth_list, all_depth_stat, snp_id, "depth_stat")
        freq_data = self.get_stat_data(sample_names, [1, 2, 3, 4, ">=5"], new_freq_stat, snp_id, "freq_stat")
        type_data = self.get_stat_data(sample_names, snp_types, snp_type_stat, snp_id, "type_stat")
        snp_pos_data = self.get_stat_data(sample_names, list(distributions), snp_pos_stat, snp_id, "snp_distribution")
        indel_pos_data = self.get_stat_data(sample_names, list(distributions), indel_pos_stat, snp_id, "indel_distribution")

        try:
            collection = self.db["sg_snp_detail"]
            if len(data_list) > 5000:
                for i in range(0, len(data_list), 3000):
                    tmp_list = data_list[i: i+3000]
                    collection.insert_many(tmp_list)
            else:
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
            main_collection.update({"_id": ObjectId(snp_id)}, {"$set": {"snp_types": snp_types, "specimen": sample_names, "distributions": list(distributions), "chroms": list(chroms)}})
        except Exception, e:
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
