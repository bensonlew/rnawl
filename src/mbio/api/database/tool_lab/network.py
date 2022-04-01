# -*- coding: utf-8 -*-
# __author__ = 'binbinzhao'

import json
from api_base import ApiBase


class Network(ApiBase):
    def __init__(self, bind_object):
        super(Network, self).__init__(bind_object)

    def add_sg_network(self, main_id, node_path, link_path):
        data_list_node = []
        data_list_link = []
        with open(node_path) as f, open(link_path) as f1:
            lines = f.readlines()
            rows = f1.readlines()
            self.bind_object.logger.info(lines)
            self.bind_object.logger.info(rows)

            index = 0
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                  "network_id": self.check_objectid(main_id),
                  "name_id": item[0],
                  "group": "",
                  "type": "node",
                }
                if item[1] != "NA":
                    insert_data["size"] = item[1]
                else:
                    insert_data["size"] = ""
                insert_data["node_id"] = index
                index += 1
                data_list_node.append(insert_data)
            for row in rows[1:]:
                ele = row.strip().split("\t")
                insert_data2 = {
                    "network_id": self.check_objectid(main_id),
                    "name": str(ele[0]) + "_" + str(ele[1]),
                    "group": "",
                    "type": "link",
                    "source": int(ele[0]),
                    "target": int(ele[1]),
                    "value": "",
                }
                if ele[2] != "NA":
                    insert_data2["size"] = ele[2]
                else:
                    insert_data2["size"] = ""
                data_list_link.append(insert_data2)
            if len(data_list_node) == 0:
                self.bind_object.logger.info("{}文件为空！".format(node_path))
            else:
                self.col_insert_data("sg_network_node", data_list_node)
            if len(data_list_link) == 0:
                self.bind_object.logger.info("{}文件为空！".format(link_path))
            else:
                self.col_insert_data("sg_network_link", data_list_link)

        update_dict = {
            "node_data": json.dumps({
                "name": "node_name",
                "condition": {'type': "node"}
            }),
            "link_data": json.dumps({
                "name": "name",
                "value": "value",
                "condition": {'type': "link"}
            }),
        }
        self.bind_object.logger.info("&&&&&&&&&&&&&&&&&&导表结束&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
        self.update_db_record("sg_network", {"_id": self.check_objectid(main_id)}, update_dict)


if __name__ == '__main__':
    a = Network(None)
    mainid = "5eead08a17b2bf7d0fd2d207"
    out_dir = '/mnt/ilustre/users/sanger-dev/workspace/20200714/Network_tsg_3421_0714130030659182_2285/output/network'
    a.add_sg_network(mainid, out_dir, out_dir)
