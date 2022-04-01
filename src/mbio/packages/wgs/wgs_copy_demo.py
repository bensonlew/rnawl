# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# modified 2018.0508


from biocluster.api.database.base import Base


class WgsCopyDemo(Base):
    """
    复制WGS demo数据
    """
    def __init__(self, old_task_id, new_task_id, new_project_sn, new_member_id):
        """
        一个demo需要连接三个字段：用户信息——项目id——任务id
        lasted modifed by hd 20180609
        1）完成pep8报错
        2）解决拉取demo时fastq列表会出现样本重复的现象，以及fastq.list打不来
        """
        super(WgsCopyDemo, self).__init__()
        self._project_type = "dna_wgs"
        self._old_task_id = old_task_id
        self._new_member_id = new_member_id     # 用户信息
        self._new_project_sn = new_project_sn   # 项目id
        self._new_task_id = new_task_id         # 任务id

    def run(self):
        self.copy_sg_task()					# sg_task 任务表！
        self.copy_sg_specimen()				# sg_specimen 样本表
#        self.copy_sg_version()				# sg_version 版本表
        self.copy_main_collection(collection="sg_mapping")			# 基因组比对主表
        self.copy_main_collection(collection="sg_snp_call")			# snp变异检测主表
        self.copy_main_collection(collection="sg_snp_anno")			# snp功能注释
        self.copy_main_collection(collection="sg_indel_call")		# indel变异检测主表
        self.copy_main_collection(collection="sg_indel_anno")		# indel功能注释
        self.copy_main_collection(collection="sg_circos")           # circos圈图主表
        self.copy_main_collection(collection="sg_anno")             # 基因注释主表
        self.copy_main_collection(collection="sg_cnv_call")			# cnv变异统计主表 task_id project_sn
        self.copy_main_collection(collection="sg_sv_call")			# sv变异统计主表
        self.copy_main_collection(collection="sg_ssr")              # ssr标准基因组主表
        # sg_task里面有genome_version_id 已经关联到基因组版本信息！！
        self.copy_main_collection(collection="sg_igv")
        self.copy_main_collection(collection="sg_specimen_group")
        self.copy_main_collection(collection="sg_specimen_qc")
        self.copy_main_collection(collection="sg_specimen_other")
        try:
            other_task_id = self._new_task_id.split('_')[1]
        except:
            other_task_id = self._new_task_id
        self.db['sg_specimen_other'].update({"task_id": self._new_task_id},
                                            {'$set': {"other_task_id": other_task_id}}, upsert=True, multi=True)
        # modified by hd 20180609

    def copy_main_collection(self, collection):					    # 赋予新的任务id 项目id
        """
        复制需要返回旧表与新表关系的数据表
        :params collection: 复制的表的名称
        """
        results = self.db[collection].find({"task_id": self._old_task_id})  # 字典
        for result in results:
            result["task_id"] = self._new_task_id
            # 每一个旧的任务id都赋值为传进来的新的任务new_task_id，结果输出为一个字典
            if "project_sn" in result.keys():
                result["project_sn"] = self._new_project_sn		# 每一个旧的项目id都赋值为传进来的新的项目new_task_id
            if "member_id" in result.keys():
                result["member_id"] = self._new_member_id
            if collection == "sg_ssr":  # 将拉取的demo里的原始的下载路径修改成新的项目里的路径，modifed by zengjing 20180627
                try:
                    new = "/".join(result["download_path"].split("/")[5:])
                    result["download_path"] = "rerewrweset/files/" + self._new_member_id + "/" + self._new_project_sn + "/" + self._new_task_id + "/" + new
                except:
                    result["download_path"] = "rerewrweset/files/" + self._new_member_id + "/" + self._new_project_sn + "/" + self._new_task_id + "/workflow_results/08.ssr/"
            result.pop("_id")								# 删掉id信息
            main_id = self.db[collection].insert_one(result).inserted_id  # import pymongo;insert_one
        print "复制{}的数据成功".format(collection)

    def copy_sg_task(self):
        """
        复制sg_task的数据
        """
        coll = self.db["sg_task"]
        result = coll.find_one({"task_id": self._old_task_id})  # demo的old task_id  &&
        if not result:
            raise Exception("sg_task表里没有找到demo的task_id：%s" % self._old_task_id)
            # raise Exception("抛出一个异常") ，raise捕捉异常，
        result.pop("_id")
        result["task_id"] = self._new_task_id
        result["project_sn"] = self._new_project_sn
        result["member_id"] = self._new_member_id
        result["is_demo"] = 1  # 该任务是demo
        #  赋值一对键值对&证明这个是demo；前端可能会用到判定；不是demo的时候，不为1或者不存在，看流程判定！！
        try:  # 原始demo任务
            result["demo_id"] = self._old_task_id
        except:											# 处理异常；
            result.pop("demo_id")
            result["demo_id"] = self._old_task_id
        self.db["sg_task"].insert_one(result)
        print "复制sg_task表成功"

    def copy_sg_specimen(self):
        """
        复制sg_specimen的数据
        """
        results = self.db["sg_specimen"].find({"task_id": self._old_task_id})     # .find || .find_one 插入整个new的字典
        news = []
        for s in results:
            s.pop("_id")
            s["task_id"] = self._new_task_id
            s["project_sn"] = self._new_project_sn
            news.append(s)          				# 遍历结束新增进去
        if news:
            self.db["sg_specimen"].insert_many(news)   # .insert_many || .insert_one
        else:
            raise Exception("sg_specimen里没有该任务{}的样本信息,请核对该任务结果是否完整")
        print "复制sg_specimen表成功"

if __name__ == "__main__":
    a = WgsCopyDemo("tsg_29745", "wgs_copy_demo1", "wgs_copy_demo2", "wgs_copy_demo3")
    a.run()
