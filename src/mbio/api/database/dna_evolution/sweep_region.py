# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify:20180920

from api_base import ApiBase
from collections import defaultdict
import os
import datetime


class SweepRegion(ApiBase):
    """
    受选择区域筛选导表
    """
    def __init__(self, bind_object):
        super(SweepRegion, self).__init__(bind_object)
        self._project_type = "dna_evolution"
        self.project_sn = self.bind_object.sheet.project_sn
        self.task_id = self.bind_object.sheet.id.split("_SweepRegion")[0]

    def add_sg_sweep_region(self, params=None, name=None):
        """
        主表sg_sweep_region
        """
        insert_data = {
            "project_sn": self.project_sn,
            "task_id": self.task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "params": params if params else "",
            "name": name if name else "origin_sweep_region",
            "desc": "受选择区域筛选",
        }
        main_id = self.db['sg_sweep_region'].insert_one(insert_data).inserted_id
        self.update_db_record("sg_sweep_region", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_sweep_region_detail(self, sweep_region_id, diff_group_name, select_region_stat):
        """
        sg_sweep_region_detail
        select_region_stat: 1-2.1.pi_tajimaD_fst.select.region.stat
        diff_group_name: 差异分组名称，A_vs_B
        """
        sweep_region_id = self.check_objectid(sweep_region_id)
        self.check_exists(select_region_stat)
        data_list = []
        diff_group = os.path.basename(select_region_stat).split(".")[0]
        pop1_pi = diff_group.split("-")[0] + "_pi"
        pop2_pi = diff_group.split("-")[1] + "_pi"
        pop1_tajima = diff_group.split("-")[0] + "_tajima"
        pop2_tajima = diff_group.split("-")[1] + "_tajima"
        with open(select_region_stat, "r") as f:
            lines = f.readlines()
            header = lines[0].strip().split("\t")
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "sweep_region_id": sweep_region_id,
                    "diff_group_name": diff_group_name,
                    "chr": item[0],
                    "start": item[1],
                    "end": item[2],
                    "gene_num": int(item[3]),
                    "snp_num": int(item[4]),
                    "indel_num": int(item[5]),
                    "average_fst": float(item[6]),
                    pop1_pi: float(item[7]),
                    pop2_pi: float(item[8]),
                    pop1_tajima: float(item[9]),
                    pop2_tajima: float(item[10])
                }
                data_list.append(insert_data)
        if len(data_list) == 0:
            self.bind_object.logger.info("{}的受选择区域结果为空！".format(diff_group_name))
        else:
            self.col_insert_data("sg_sweep_region_detail", data_list)

    def update_diff_group(self, sweep_region_id, diff_group, sweep_region_dir):
        """
        更新sg_sweep_regin的diff_group
        diff_group: 差异分组list,["A_vs_B", "B_vs_C"]
        sweep_region_dir: 上传的结果文件夹
        """
        sweep_region_id = self.check_objectid(sweep_region_id)
        self.update_db_record("sg_sweep_region", {"main_id": sweep_region_id}, {"diff_group": diff_group, "sweep_region_dir": sweep_region_dir})
        self.bind_object.logger.info("更新差异分组:{}成功！".format(diff_group))


if __name__ == "__main__":
    a = SweepRegion(None)
    project_sn = "test_zj"
    task_id = "tsg_32120"
    sweep_region_id = a.add_sg_sweep_region(project_sn, task_id)
    diff_group_name = "1_vs_2"
    select_region_stat = "/mnt/ilustre/users/sanger-dev/workspace/20180920/Single_sweep_region_single3/SweepRegionSingle/SweepRegionStat/output/1-2.1.pi_tajimaD_fst.select.region.anno"
    a.add_sg_sweep_region_detail(sweep_region_id, diff_group_name, select_region_stat)
    diff_group_name = "1_vs_3"
    select_region_stat = "/mnt/ilustre/users/sanger-dev/workspace/20180920/Single_sweep_region_single3/SweepRegionSingle/SweepRegionStat/output/1-2.2.pi_tajimaD_fst.select.region.anno"
    a.add_sg_sweep_region_detail(sweep_region_id, diff_group_name, select_region_stat)
    diff_group = ["1_vs_2", "1_vs_3"]
    a.update_diff_group(sweep_region_id, diff_group)
