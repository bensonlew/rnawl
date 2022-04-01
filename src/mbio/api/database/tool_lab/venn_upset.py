# !usr/bin/python
# -*- coding: utf-8 -*-
from api_base import ApiBase
import json


class VennUpset(ApiBase):
    """..."""
    def __init__(self, bind_object):
        super(VennUpset, self).__init__(bind_object)

    def sg_manhattan(self, task_id, project_sn, params):  # 用于测试建立主表
        """
        曼哈顿图主表
        :return:
        """
        main_id = self.add_main_table("sg_manhattan", task_id, project_sn, params, "origin_sg_manhattan",
                                      "曼哈顿图导表", "开始进行曼哈顿图导表")
        self.update_db_record("sg_manhattan", {"_id": main_id}, {"main_id": main_id, "status": "end"})
        return main_id

    def add_sg_venn_upset(self, main_id, path):
        id = self.check_objectid(main_id)
        venn_upset_data = {"name": "name",
                           # "data": "value",
                           "category": "category",
                           "condition": {"type": "venn_upset"}}
        self.update_db_record("sg_upset", {"_id": id}, {"venn_upset_data": json.dumps(venn_upset_data)})
        data_list = []
        with open(path) as f:
            lines = f.readlines()
            for line in lines:
                insert_data = {
                    "name": line.strip().split("\t")[0],
                    # "value": line.strip().split("\t")[1:],
                    "data": line.strip().split("\t")[1:],
                    "type": "venn_upset",
                    "category": line.strip().split("\t")[0],
                    "upset_id": id
                }
                data_list.append(insert_data)
        if len(data_list) == 0:
            self.bind_object.logger.info("{}文件为空！".format(path))
        else:
            self.col_insert_data("sg_upset_detail", data_list)









if __name__ == "__main__":
    a = VennUpset(None)
    task_id = "tool_lab"
    project_sn = "tool_lab_20200416"
    # a.sg_manhattan(task_id, project_sn, "manhattan")
    a.add_sg_venn_upset("5e97eeb517b2bf6d8f8f35f8", "/mnt/ilustre/users/sanger-dev/workspace/20200507/VennUpset_tsg_3421_0507165024849807_5569/output/venn_upset/venn_upset.txt")

