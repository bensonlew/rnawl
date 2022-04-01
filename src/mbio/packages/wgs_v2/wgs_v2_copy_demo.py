# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 20190505

from biocluster.api.database.base import Base


class WgsV2CopyDemo(Base):
    """
    复制DNA wgs_V2 demo数据
    """
    def __init__(self, old_task_id, new_task_id, new_project_sn, new_member_id):
        """
        """
        super(WgsV2CopyDemo, self).__init__()
        self._project_type = "dna_wgs_v2"
        self._old_task_id = old_task_id
        self._new_member_id = new_member_id
        self._new_project_sn = new_project_sn
        self._new_task_id = new_task_id

    def run(self):
        self.copy_sg_task()
        self.copy_sg_specimen()
        # self.copy_sg_version()
        self.copy_main_collection(collection="sg_specimen_group")    # 软件列表
        self.copy_main_collection(collection="sg_software")    # 软件列表
        self.copy_main_collection(collection="sg_specimen_qc")
        self.copy_main_collection(collection="sg_mapping")			# 基因组比对主表
        self.copy_main_collection(collection="sg_coverage_window")			# 基因组比对主表
        self.copy_main_collection(collection="sg_variant_call")			# snp变异检测主表
        self.copy_main_collection(collection="sg_variant_anno")			# snp功能注释
        self.copy_main_collection(collection="sg_variant_compare")   # indel变异检测主表
        self.copy_main_collection(collection="sg_variant_compare_filter")
        # self.copy_main_collection(collection="sg_chromosome_window")		# indel变异检测主表
        self.copy_main_collection(collection="sg_assembly")		# indel功能注释
        self.copy_main_collection(collection="sg_circos")           # circos圈图主表
        self.copy_main_collection(collection="sg_cnv_call")             # 基因注释主表
        self.copy_main_collection(collection="sg_cnv_compare")
        self.copy_main_collection(collection="sg_crispr_analysis")
        self.copy_main_collection(collection="sg_gene_seq")
        self.copy_main_collection(collection="sg_marker_distribution")
        self.copy_main_collection(collection="sg_marker_position_table")
        self.copy_main_collection(collection="sg_primer")
        self.copy_main_collection(collection="sg_sequence")
        self.copy_main_collection(collection="sg_region_anno")
        self.copy_main_collection(collection="sg_ssr_compare")
        self.copy_main_collection(collection="sg_sv_call")
        self.copy_main_collection(collection="sg_sv_compare")
        self.copy_main_collection(collection="sg_transgenosis")
        self.copy_main_collection(collection="sg_variant_site_table")
        self.copy_main_collection(collection="sg_results")
        self.copy_main_collection(collection="sg_ssr_marker")



    def copy_sg_task(self):
        """
        复制sg_task的数据
        """
        coll = self.db["sg_task"]
        result = coll.find_one({"task_id": self._old_task_id})
        if not result:
            raise Exception("sg_task表里没有找到demo的task_id：%s" % self._old_task_id)
        result.pop("_id")
        result["task_id"] = self._new_task_id
        result["project_sn"] = self._new_project_sn
        result["member_id"] = self._new_member_id
        result["is_demo"] = 1  # 该任务是demo
        try:  # 原始demo任务
            result["demo_id"] = self._old_task_id
        except:
            result.pop("demo_id")
            result["demo_id"] = self._old_task_id
        self.db["sg_task"].insert_one(result)
        print "复制sg_task表成功"

    def copy_sg_specimen(self):
        """
        复制sg_specimen的数据
        """
        results = self.db["sg_specimen"].find({"task_id": self._old_task_id})
        news = []
        for s in results:
            s.pop("_id")
            s["task_id"] = self._new_task_id
            s["project_sn"] = self._new_project_sn
            news.append(s)
        if news:
            new_result = self.db["sg_specimen"].insert_many(news)
        else:
            raise Exception("sg_specimen里没有该任务{}的样本信息,请核对该任务结果是否完整")
        print "复制sg_specimen表成功"

    """
    需要确保这个sg_version是否还需要"""
    # def copy_sg_version(self):
    #     """
    #     复制sg_version的数据
    #     """
    #     result = self.db["sg_version"].find_one({"task_id": self._old_task_id})
    #     result.pop("_id")
    #     result["task_id"] = self._new_task_id
    #     self.db["sg_version"].insert_one(result)
    #     print "复制sg_version表成功"

    def copy_main_collection(self, collection):
        """
        复制需要返回旧表与新表关系的数据表
        :params collection: 复制的表的名称
        """
        results = self.db[collection].find({"task_id": self._old_task_id})
        if not results:
            raise Exception("没有在表:{}中找到task_id:{}，请检查".format(collection, self._old_task_id))
        news = []
        for result in results:
            result["task_id"] = self._new_task_id
            if "project_sn" in result.keys():
                result["project_sn"] = self._new_project_sn
            if "member_id" in result.keys():
                result["member_id"] = self._new_member_id
            if collection == "sg_ssr":  # 将拉取的demo里的原始的下载路径修改成新的项目里的路径，modifed by zengjing 20180627
                try:
                    new = "/".join(result["download_path"].split("/")[5:])
                    result["download_path"] = "rerewrweset/files/" + self._new_member_id + "/" + self._new_project_sn + "/" + self._new_task_id + "/" + new
                except:
                    result["download_path"] = "rerewrweset/files/" + self._new_member_id + "/" + self._new_project_sn + "/" + self._new_task_id + "/workflow_results/08.ssr/"
            if collection == "sg_task":  # sg_task表增加原始demo的相关信息
                result["is_demo"] = 1
                try:
                    result["demo_id"] = self._old_task_id
                except:
                    result.pop("demo_id")
                    result["demo_id"] = self._old_task_id
            result.pop("_id")
            news.append(result)
        if news:
            self.db[collection].insert_many(news)
        print "复制{}的数据成功".format(collection)


if __name__ == "__main__":
    a = WgsV2CopyDemo("sanger_85433", "wgs_v2_new1", "wgs_v2_new1-1", "wgs_v2_new1-1new")
    a.run()
