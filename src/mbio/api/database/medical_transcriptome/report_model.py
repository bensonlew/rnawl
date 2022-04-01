# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from mbio.api.database.medical_transcriptome.api_base import ApiBase
from biocluster.api.database.base import report_check
import json
from bson.objectid import ObjectId
import os
import sys
import datetime

class ReportModel(ApiBase):
    '''
    last_modify: 2021.08.06
    '''
    def __init__(self, bind_object):
        super(ReportModel, self).__init__(bind_object)


    def get_report_model(self, task_id):
        one_touch_dict = dict()
        collection = self.db["sg_one_touch"]
        results = collection.find({"task_id": task_id}, {"_id": 0})
        deposits_dict = dict()
        for result in results:
            deposits_dict[result["table_deposit_id"]] = result["id"]
            result["table_deposit_id"] = str(result["table_deposit_id"])
            one_touch_dict[result['id']] = result

        with open("one_touch.json", 'w') as f:
            f.write(json.dumps(one_touch_dict, indent=4, ensure_ascii=False).encode('utf8').decode())

        deposit_dict = dict()
        collection = self.db["sg_result_table_deposit"]
        results = collection.find({"_id": {"$in": deposits_dict.keys()}})
        for result in results:
            a = result['_id']
            a_id = deposits_dict[a]
            result_clean = result.pop("_id")
            deposit_dict[a_id] = result

        with open("deposit.json", 'w') as f:
            f.write(json.dumps(deposit_dict, indent=4, ensure_ascii=False).encode('utf8').decode())

    def get_geneset_name2id(self, task_id):
        coll = self.db["sg_geneset"]
        genesets_records = coll.find({"task_id": task_id})
        self.geneset_name2id = {record["name"]: str(record["main_id"]) for record in genesets_records}

    # def convert_main_dict2query_dict(self, main_dict):
    #     query_update_dict = dict()
    #     for k, v in main_dict:
    #         if k == "sample":
    #             pass
    #         if k == "gene_set":

    def choose_one_from_record(self, main_dict, records):
        if "gene_set" in main_dict:
            for record in records:
                geneset_name = main_dict["gene_set"]
                if geneset_name in self.geneset_name2id:
                    pass
                else:
                    if "_up" in geneset_name:
                        geneset_name = geneset_name.split("_up")[0]
                    if "_down" in geneset_name:
                        geneset_name = geneset_name.split("_down")[0]
                        
                if self.geneset_name2id[geneset_name] in record["params"]:
                    record_choose = record
                    break

            if "record_choose" in locals():
                pass
            else:
                # print(main_dict)
                # print("not find geneset {}".format(self.geneset_name2id[main_dict["gene_set"]]                ))
                records.rewind()
                record_choose = records.next()

        else:
            record_choose = records.next()
        return record_choose

    def get_result_table_id(self, task_id, p_dict):
        '''
        获取结果表ID
        '''
        query_dict = {"task_id": task_id}
        coll = self.db[p_dict["main_coll"]]
        if "main_dict" in p_dict:
            query_dict.update(query_dict)
            main_records = coll.find(query_dict)
            if main_records.count() == 0:
                print("not find {} {}".format(query_dict, p_dict))
                result_table_id = "None"
                result_table_name = "undefined"
            else:
                record_choose = self.choose_one_from_record(p_dict["main_dict"], main_records)
                result_table_id = record_choose["main_id"]
                result_table_name = record_choose.get("name", "undefined")
        else:
            main_record = coll.find_one({"task_id": task_id})
            if not main_record:
                print("not find {} {}".format(query_dict, p_dict))
                result_table_id = None
                result_table_name = "undefined"
            else:
                result_table_id = main_record["main_id"]
                result_table_name = main_record.get("name", "undefined")
        return result_table_id, result_table_name

    def add_report_image(self, task_id, report_config, img_path, member_id=""):
        '''
        倒入结题报告图片表格
        '''

        self.get_geneset_name2id(task_id)
        one_touch_collection = self.db["sg_one_touch"]
        deposit_collection = self.db["sg_result_table_deposit"]

        script_path = os.path.split(os.path.realpath(__file__))[0]

        one_touch_dicts = self.load_one_touch(os.path.join(script_path, "one_touch.json"))
        deposit_dicts = self.load_deposit(os.path.join(script_path, "deposit.json"))

        with open(report_config, 'r') as f:
            report_dict = json.load(f)
        serial_no = 0
        for p_id, p_dict in report_dict.items():
            if p_id in one_touch_dicts:
                serial_no += 1
                one_touch_dict = one_touch_dicts[p_id]
                one_touch_dict_new = one_touch_dict.copy()

                # 获取deposit模板信息
                # print one_touch_dict["table_deposit_id"]
                if p_id not in deposit_dicts:
                    continue
                deposit_dict = deposit_dicts[p_id]
                deposit_dict_new = deposit_dict.copy()

                # 获取分析主表
                # print("p_dict {}".format(p_dict))
                result_table_id, result_table_name = self.get_result_table_id(task_id, p_dict)

                # print p_dict

                # 生成该项目deposit数据 
                if "with_table" in deposit_dict_new:
                    deposit_dict_table = deposit_dict_new.copy()
                    deposit_dict_table.pop("img_url")
                    deposit_dict_table.pop("img_type")
                    deposit_dict_table.update({
                        "task_id": task_id,
                        "result_table_id": result_table_id,
                        "created_time": datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                        "result_table_name": result_table_name,
                        "type": "table",
                        "serial_no": serial_no
                    })
                    serial_no += 1
                    deposit_collection.insert_one(deposit_dict_table).inserted_id

                deposit_dict_new.update({
                    "task_id": task_id,
                    "img_url": os.path.join(img_path, p_dict["img"]),
                    "result_table_id": result_table_id,
                    "result_table_name": result_table_name,
                    "created_time": datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                    "serial_no": serial_no
                })
                deposit_id = deposit_collection.insert_one(deposit_dict_new).inserted_id
                # print("insert {}".format(deposit_dict))

                # 生成项目one_touch_dict数据
                
                
                one_touch_dict_new.update({
                    "table_deposit_id": deposit_id,
                    "task_id": task_id,
                    "img_url": os.path.join(img_path, p_dict["img"]),
                    "result_table_id": result_table_id,
                    "result_table_name": result_table_name,
                    "created_time": datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                    "create_id": member_id

                })
                # print one_touch_dict_new
                print("insert {}".format(one_touch_dict_new))
                one_touch_collection.insert_one(one_touch_dict_new)


    def load_one_touch(self, json_f):
        with open(json_f, 'r') as f:
            one_touch_dict = json.load(f)
            return one_touch_dict

    def load_deposit(self, json_f):
        with open(json_f, 'r') as f:
            deposit_dict = json.load(f)
            return deposit_dict

if __name__ == '__main__':
    os.environ["current_mode"]="workflow"
    os.environ["NTM_PORT"]="7322"
    os.environ["WFM_PORT"]="7321"

    report = ReportModel(None)
    report._config.DBVersion = 0
    # # report._config.ProjectID = "v0afleo6p1211n885l0oo0av7v"
    # report.get_report_model(task_id=sys.argv[1])
    if len(sys.argv) == 2:
        report._config.ProjectID = "v0afleo6p1211n885l0oo0av7v"
        report._config.DBVersion = 1
        report.get_report_model(task_id=sys.argv[1])
    else:
        report.add_report_image(
            task_id=sys.argv[1],
            report_config=sys.argv[2],
            img_path= "test"
        )
