# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# modified 2018.0830

from api_base import ApiBase
from collections import defaultdict
import os
import datetime
import json


class Hapmap(ApiBase):
    """
    群体进化，单倍体图谱
    """
    def __init__(self, bind_object):
        super(Hapmap, self).__init__(bind_object)
        self._project_type = "dna_evolution"

    def add_sg_hapmap(self, project_sn, task_id, params=None, name=None):
        """
        单倍体图谱主表
        需在modu:hapmap内更新主表的header和region_path字段，见updata_sg_hapmap
        """
        if params:
            new_params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        else:
            new_params = "null"
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "params": new_params,
            "name": name if name else "origin_sg_hapmap",
            "desc": "单倍体图谱接口主表",
        }
        main_id = self.db['sg_hapmap'].insert_one(insert_data).inserted_id
        self.update_db_record("sg_hapmap", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_haploid(self, main_id, file_path, analysis, trait, task_id, area):
        """
        单倍体图普导表
        analysis:分析方法
        :param main_id:
        :param file_path: 每个区域的pop.ld文件
        e.g.["/mnt/ilustre/users/sanger-dev/workspace/20181025/Single_hapmap_7_20181025/Hapmap/RegionHapmap/output/chr3_1_217536/pop.ld,
        /mnt/ilustre/users/sanger-dev/workspace/20181025/Single_hapmap_7_20181025/Hapmap/RegionHapmap/output/chr4_1_218215/pop.ld']
        """
        origin_id = self.check_objectid(main_id)  # 检查id是否是OBjectID
        insert_data = {
            "origin_id": origin_id,
            "name": '单倍体图谱',
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "task_id": task_id,
            "type": 1,
            "trait": trait,
            "analysis": analysis,
            "area": area
        }
        heatmap_id = self.db["sg_heatmap"].insert_one(insert_data).inserted_id
        for i in area:
            new_file_path = os.path.join(os.path.join(file_path, "region_dir"), (i + '/' + i + '.LDheatmap'))
            self.add_sg_haploid_detail(file_path=new_file_path, heatmap_id=heatmap_id, area=i)

    def add_sg_haploid_detail(self, file_path, heatmap_id, area):
        insert_data = {
                    "heatmap_id": heatmap_id,
                    "area": area,
                    "file_path": file_path
                }
        self.db["sg_heatmap_detail"].insert_one(insert_data)

    # def add_sg_haploid_detail(self, file_path, heatmap_id):
    #     """
    #     单倍体图细节表，一个区间一条记录。
    #     """
    #     x = []
    #     y = []
    #     data_dic = {}
    #     datas = []
    #     mongo_list = []
    #     with open(file_path, 'r')as rd:
    #         lines = rd.readlines()
    #         for line in lines[1:]:
    #             temp = line.strip().split()
    #             x.append(int(temp[1]))
    #             y.append(int(temp[4]))
    #             data_dic[str(int(temp[1]))+'||'+str(int(temp[4]))] = float(temp[6])
    #     all = x + y
    #     c = list(set(all))
    #     c.sort()
    #     n = 0
    #     print '2222---------------'
    #     print len(c)
    #     for i in c:
    #         data_list_append = []
    #         for l in range(n, len(c)):
    #             list_in = str(c[l])+'||'+str(i)
    #             list_out = str(i)+'||'+str(c[l])
    #             append_value = 0
    #             for key, value in data_dic.items():
    #                 if key == list_in or key == list_out:
    #                     append_value = value
    #             data_list_append.append(append_value)
    #         n += 1
    #         datas.append(data_list_append)
    #     print datas
    #     temp = file_path.strip().split("/")
    #     area = temp[-2]
    #     if len(datas) > 100:
    #         data_splist = self.splist(datas, 10)
    #         for i in range(10):
    #             insert_data = {
    #                 "heatmap_id": heatmap_id,
    #                 "area": area,
    #                 "value": data_splist[i],
    #                 "sort": i
    #             }
    #             mongo_list.append(insert_data)
    #         print '4444---------------'
    #         if len(mongo_list) == 0:
    #             self.bind_object.logger.info("{}的结果为空！".format(file_path))
    #         else:
    #             self.col_insert_data("sg_heatmap_detail", mongo_list)
    #     else:
    #         insert_data = {
    #             "heatmap_id": heatmap_id,
    #             "area": area,
    #             "value": datas,
    #             "sort": 0
    #         }
    #         self.db["sg_heatmap_detail"].insert_one(insert_data)

    def add_sg_hapmap_stat(self, hapmap_id, header, region_dir):
        """
        单倍体图谱--关联区域统计
        /mnt/ilustre/users/sanger-dev/workspace/20180829/Single_region_hapmap_220180829/RegionHapmap/output/TRAIT2.GLM.select.region
        补pl to TRAIT2.GLM.select.region.xls
        """
        hapmap_id = self.check_objectid(hapmap_id)   # 检查id是否是OBjectID
        self.check_exists(region_dir)   # 检查文件是否存在
        region_list = []
        region_list_2 = os.listdir(region_dir)
        for i in region_list_2:
            if i.endswith('.xls'):
                region_list.append(i)
        for file in region_list:
            trait_name = file.strip().split('.')[0]
            analysis = file.strip().split('.')[1]
            data_list = []
            with open(region_dir + "/" + file, 'r') as r:
                lines = r.readlines()
                for line in lines:
                    if line.startswith("#") or line.startswith("@"):
                        pass
                    else:
                        tmp = line.strip().split('\t')
                        insert_data = {
                            "hapmap_id": hapmap_id,
                            "trait": trait_name,
                            "chr": tmp[0],
                            "start": int(tmp[1]),
                            "end": int(tmp[2]),
                            "min_pvalue": float(tmp[3]),
                            "gene_num": int(tmp[4]),
                            "analysis": analysis
                        }
                        if header == 1:
                            insert_data['snp_num'] = int(tmp[6])
                        else:
                            insert_data['snp_num'] = int(tmp[6])
                            insert_data['indel_num'] = int(tmp[5])
                        data_list.append(insert_data)
                if len(data_list) == 0:
                    self.bind_object.logger.info("为空！{}".format(file))
                    print "为空！"
                else:
                    self.col_insert_data("sg_hapmap_stat", data_list)

    def splist(self, l, x):
        """
        用于自动等分list
        :param l:list
        :param x:按多少等分
        :return:
        """
        length = (len(l) / x if len(l) % x == 0 else len(l) / x + 1)
        return [l[m:m + length] for m in range(0, len(l), length)]


if __name__ == "__main__":
    a = Hapmap(None)
    # project_sn = 'dna_evolution_cui'
    # task_id = 'dna_evolution_cui'
    # params = {"sg_gwas_analysis_id": "5adecf95a4e1af4fa5c55351"}
    # hapmap_id = a.add_sg_hapmap(project_sn, task_id, params)
    # header = 'snp,indel'
    # region_dir = '/mnt/ilustre/users/sanger-dev/home/cuiqingmei/1.project/4.evolution/api.hapmap/region_dir'
    # a.add_sg_hapmap_stat(hapmap_id, header, region_dir)
    main_id = '5adecf95a4e1af4fa5c55351'
    file_path = ['/mnt/ilustre/users/sanger-dev/workspace/20181113/Single_tsg_32120_Hapmap_1113152720629524/Hapmap/RegionHapmap/output/chr3_1_117536/pop.ld',
                 '/mnt/ilustre/users/sanger-dev/workspace/20181113/Single_tsg_32120_Hapmap_1113152720629524/Hapmap/RegionHapmap/output/chr4_1_118215/pop.ld',
                 '/mnt/ilustre/users/sanger-dev/workspace/20181113/Single_tsg_32120_Hapmap_1113152720629524/Hapmap/RegionHapmap/output/chr7_1_135456/pop.ld',
                 '/mnt/ilustre/users/sanger-dev/workspace/20181113/Single_tsg_32120_Hapmap_1113152720629524/Hapmap/RegionHapmap/output/chr9_1_103584/pop.ld']
    analysis = "FarmCPU"
    trait = "FA"
    task_id = 'dna_evolution'
    a.add_sg_haploid(main_id, file_path, analysis, trait, task_id)
