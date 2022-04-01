# !usr/bin/python
# -*- coding: utf-8 -*-
from api_base import ApiBase
import os


class Circos(ApiBase):
    """..."""
    def __init__(self, bind_object):
        super(Circos, self).__init__(bind_object)

    def sg_manhattan(self, task_id, project_sn, params):  # 用于测试建立主表
        """
        曼哈顿图主表
        :return:
        """
        main_id = self.add_main_table("sg_manhattan", task_id, project_sn, params, "origin_sg_manhattan",
                                      "曼哈顿图导表", "开始进行曼哈顿图导表")
        self.update_db_record("sg_manhattan", {"_id": main_id}, {"main_id": main_id, "status": "end"})
        return main_id

    def add_sg_circos_path(self, main_id, **params_dict):
        id = self.check_objectid(main_id)
        self.update_db_record("sg_circos", {"_id": id}, {"circos_path": params_dict["path"],
                                                         "outerRadius": params_dict["outerRadius"],
                                                         "innerRadius": params_dict["innerRadius"],
                                                         "circle_num": params_dict["circle_num"],
                                                         "circos_type": params_dict["circos_type"]})


if __name__ == "__main__":
    a = Circos(None)
    task_id = "tool_lab"
    project_sn = "tool_lab_20200416"
    # a.sg_manhattan(task_id, project_sn, "manhattan")
    # a.add_sg_circos_chr("5e97eeb517b2bf6d8f8f35f8", "/mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/tool_lab/circos/Circos1.txt")
    # a.add_sg_circos_chart("5e97eeb517b2bf6d8f8f35f8", "/mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/tool_lab/circos/circos2.txt")
    a.add_sg_circos_path("5e97eeb517b2bf6d8f8f35f8",
                         **{"outerRadius": 252, 'innerRadius': 246, 'path': "AAA", 'circle_num': 4}
                         )



