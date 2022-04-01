# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.02.28

from biocluster.api.database.base import Base


class BsaCopyDemo(Base):
    """
    复制BSA demo数据
    """
    def __init__(self, old_task_id, new_task_id, new_project_sn, new_member_id):
        super(BsaCopyDemo, self).__init__()
        self._project_type = "bsa"
        self._old_task_id = old_task_id
        self._new_task_id = new_task_id
        self._new_project_sn = new_project_sn
        self._new_member_id = new_member_id

    def run(self):
        self.copy_sg_task()
        self.copy_sg_specimen()
        self.copy_sg_version()
        self.copy_main_collection(collection="sg_specimen_qc")
        self.copy_main_collection(collection="sg_mapping")
        self.copy_main_collection(collection="sg_snp_call")
        self.copy_main_collection(collection="sg_snp_anno")
        self.copy_main_collection(collection="sg_indel_call")
        self.copy_main_collection(collection="sg_indel_anno")
        self.copy_main_collection(collection="sg_anno")
        self.copy_main_collection(collection="sg_slidingwin")
        self.copy_main_collection(collection="sg_region")
        self.copy_main_collection(collection="sg_region_anno")  # region_id
        self.copy_main_collection(collection="sg_genome")

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

    def copy_sg_version(self):
        """
        复制sg_version的数据
        """
        result = self.db["sg_version"].find_one({"task_id": self._old_task_id})
        result.pop("_id")
        result["task_id"] = self._new_task_id
        self.db["sg_version"].insert_one(result)
        print "复制sg_version表成功"

    def copy_sg_manhattan(self, slidingwin_ids):
        """
        复制sg_manhattan的数据
        """
        results = self.db["sg_region"].find({"task_id": self._old_task_id})
        data_list = []
        for result in results:
            s = self.db["sg_manhattan"].find_one({"origin_id": result["slidingwin_id"]})
            s["origin_id"] = slidingwin_ids[s["origin_id"]]
            s.pop["_id"]
            data_list.append(s)
        self.db["sg_manhattan"].insert_many(data_list)
        print "复制sg_manhattan表成功"

    def copy_main_collection(self, collection):
        """
        复制需要返回旧表与新表关系的数据表
        :params collection: 复制的表的名称
        """
        results = self.db[collection].find({"task_id": self._old_task_id})
        for result in results:
            result["task_id"] = self._new_task_id
            result["project_sn"] = self._new_project_sn
            if "member_id" in result.keys():
                result["member_id"] = self._new_member_id
            result.pop("_id")
            main_id = self.db[collection].insert_one(result).inserted_id
        print "复制{}的数据成功".format(collection)


if __name__ == "__main__":
    a = BsaCopyDemo("Oryza_sativa_mbmpwpwb", "demo_test3", "demo_test", "demo_test")
    a.run()
