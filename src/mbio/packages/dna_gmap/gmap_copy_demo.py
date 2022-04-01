# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 20180717

from biocluster.api.database.base import Base


class GmapCopyDemo(Base):
    """
    复制DNA GMAP demo数据
    """
    def __init__(self, old_task_id, new_task_id, new_project_sn, new_member_id):
        """
        """
        super(GmapCopyDemo, self).__init__()
        self._project_type = "dna_gmap"
        self._old_task_id = old_task_id
        self._new_member_id = new_member_id
        self._new_project_sn = new_project_sn
        self._new_task_id = new_task_id

    def run(self):
        self.copy_main_collection(collection="sg_task")				# sg_task 任务表！
        self.copy_main_collection(collection="sg_specimen")	        # sg_specimen 样本表
        self.copy_main_collection(collection="sg_mapping")			# 基因组比对主表
        self.copy_main_collection(collection="sg_snp_call")			# snp变异检测主表
        self.copy_main_collection(collection="sg_snp_anno")			# snp功能注释
        self.copy_main_collection(collection="sg_indel_call")		# indel变异检测主表
        self.copy_main_collection(collection="sg_indel_anno")		# indel功能注释
        self.copy_main_collection(collection="sg_circos")           # circos圈图主表
        self.copy_main_collection(collection="sg_anno")             # 基因注释主表
        # self.copy_main_collection(collection="sg_cnv_call")			# cnv变异统计主表 task_id project_sn
        # self.copy_main_collection(collection="sg_sv_call")			# sv变异统计主表
        # self.copy_main_collection(collection="sg_ssr")              # ssr标准基因组主表
        self.copy_main_collection(collection="sg_igv")
        # self.copy_main_collection(collection="sg_specimen_group")
        self.copy_sg_specimen_child()
        # self.copy_main_collection(collection="sg_specimen_child")
        self.copy_main_collection(collection="sg_specimen_qc")
        self.copy_main_collection(collection="sg_specimen_other")
        try:
            other_task_id = self._new_task_id.split('_')[1]
        except:
            other_task_id = self._new_task_id
        self.db['sg_specimen_other'].update({"task_id": self._new_task_id},
                                            {'$set': {"other_task_id": other_task_id}}, upsert=True, multi=True)
        self.copy_main_collection(collection="sg_marker")
        self.copy_main_collection(collection="sg_binmarker")
        self.copy_sgsg_feature_file()
        self.copy_main_collection(collection="sg_lg")
        self.copy_main_collection(collection="sg_tree")
        self.copy_main_collection(collection="sg_evalutaion")
        self.copy_main_collection(collection="sg_subtype_matrix")

    def copy_sg_specimen_child(self):
        """
        复制sg_specimen_child
        """
        results = self.db["sg_specimen_child"].find({"task_id": self._old_task_id})
        if not results:
            raise Exception("没有在表:{}中找到task_id:{}，请检查".format("sg_specimen_child", self._old_task_id))
        news = []
        for result in results:
            result["task_id"] = self._new_task_id
            result["child_task_id"] = self._new_task_id.split("_")[-1]
            result.pop("_id")
            news.append(result)
        if news:
            self.db["sg_specimen_child"].insert_many(news)
        print "复制{}的数据成功".format("sg_specimen_child")

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
                    result["download_path"] = "rerewrweset/files/" + self._new_member_id + "/" +\
                                              self._new_project_sn + "/" + self._new_task_id + "/" + new
                except:
                    result["download_path"] = "rerewrweset/files/" + self._new_member_id + "/" + self._new_project_sn \
                                              + "/" + self._new_task_id + "/workflow_results/08.ssr/"
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

    def copy_sgsg_feature_file(self):
        """
        复制表sg_feature_file、sg_feature_file_detail
        """
        results = self.db["sg_feature_file"].find({"task_id": self._old_task_id})
        if not results:
            raise Exception("没有在表:sg_feature_file中找到task_id:{}，请检查".format(self._old_task_id))
        for result in results:
            result["task_id"] = self._new_task_id
            result["project_sn"] = self._new_project_sn
            old_feature_file_id = result["_id"]
            result.pop("_id")
            feature_file_id = self.db["sg_feature_file"].insert_one(result).inserted_id
            results_ = self.db["sg_feature_file_detail"].find({"feature_file_id": old_feature_file_id})
            print "复制sg_feature_file的数据成功"
            news = []
            for result_ in results_:
                result_["feature_file_id"] = feature_file_id
                result_.pop("_id")
                news.append(result_)
            if news:
                self.db["sg_feature_file_detail"].insert_many(news)
                print "复制sg_feature_file_detail的数据成功"


if __name__ == "__main__":
    a = GmapCopyDemo("tsg_29745", "wgs_copy_demo1", "wgs_copy_demo2", "wgs_copy_demo3")
    a.run()
