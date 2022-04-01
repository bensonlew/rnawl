# -*- coding: utf-8 -*-
## __author__: guhaidong
from biocluster.api.database.base import Base
from bson.objectid import ObjectId
import types
import ConfigParser

class CleanMongo(Base):
    def __init__(self, bind_object, project_type):
        super(CleanMongo, self).__init__(bind_object)
        self._project_type = project_type
        self.rcf = ConfigParser.RawConfigParser()
        self.status_db = self.db["sg_status"]  # sg_record_status? sg_result_table_deposit?
        self.record_db = self.db["sg_record_status"]
        self.report_db = self.db["sg_result_table_deposit"]
        self.main_table_dic = dict()
        self.main_detail_map = list()
        self.task_id = ""

    def clean_all_var(self):
        self.main_table_dic = dict()
        self.main_detail_map = list()
        self.task_id = ""

    def get_main_table_info(self, task_id):
        self.task_id = task_id
        infos = self.status_db.find({"task_id": task_id})
        if infos and infos.count() != 0:
            for info in infos:
                main_table_name = info["type_name"]
                main_table_id = info["table_id"]
                if self.main_table_dic.has_key(main_table_name):
                    self.main_table_dic[main_table_name].append(main_table_id)
                else:
                    self.main_table_dic[main_table_name] = [main_table_id]
        else:
            raise Exception("找不到任务: %s" % task_id)

    def get_detail_table_info(self, table, from_status=True):
        """
        [detail_list]
        main_table_name = detail_name1,detail_name2...
        [main_id_map]
        default_main_name = _id
        ***_main_name = main_id_map
        [detail_info]
        detail_name1 = main_table_head
        detail_name2 = main_table_head
        """
        self.rcf.read(table)
        if from_status:
            main_collections = self.main_table_dic.keys()
        else:
            main_collections = self.db.collection_names()
        for name in main_collections:
            try:
                detail_str = self.rcf.get("detail_list", name)
            except:
                if from_status:
                    print "主表%s在交互分析运行记录中找到，但未配置相关明细表信息，请检查" % name
                else:
                    print "主表%s在数据库中找到，但未配置相关明细表信息，请检查" % name
                continue
            try:
                main_id_map = self.rcf.get("main_id_map", name + "_main_name")
            except:
                main_id_map = self.rcf.get("main_id_map", "default_main_name")
            detail_table_info = {
                "name": name,
                "main_id_name": main_id_map,
                "info": {}
            }
            detail_table_list = [] if "" in detail_str.split(",") else detail_str.split(",")
            for detail_table in detail_table_list:
                main_table_head = self.rcf.get("detail_info", detail_table)
                detail_table_info['info'][detail_table] = main_table_head
            self.main_detail_map.append(detail_table_info)

    def get_detail_table_info_from_mongo(self, coll_name, main_id_name="main_id", main_id_loc="asso_id"):
        """
        coll_name: 用于记录mongo表间对应关系的collection name
        main_id_name: 默认在主表中对应main_id的字段，默认为main_id
        main_id_loc: 在对应关系表中特别规定主表中main_id的字段，默认为asso_id
        """
        self.rcf = None
        result = self.db[coll_name].find({"is_detail": "n"})
        for one in result:
            analysis = one["analysis"]
            name = one["table"]
            if one[main_id_loc] == "n":
                main_id_name = main_id_name
            else:
                main_id_name = one[main_id_loc]
            detail_result = self.db[coll_name].find({"is_detail": "y", "analysis": analysis})
            info = {}
            for one_detail in detail_result:
                info[one_detail['table']] = one_detail['asso_id']
            self.main_detail_map.append({
                "name": name,
                "main_id_name":main_id_name,
                "info": info
            })

    def insert_detail_table_info_to_mongo(self, coll_name, default_main_name="_id"):
        is_first_list = self.rcf.get("other conditions", "is_first").split(",")
        is_first_list = [] if is_first_list == [""] else is_first_list
        for one_ana in self.main_detail_map:
            analysis = table = one_ana['name']
            main_id_name = "n" if one_ana["main_id_name"] == default_main_name else one_ana["main_id_name"]
            is_first = "y" if table in is_first_list else "n"
            is_workflow = "y"
            data = {
                "is_first": is_first,
                "is_workflow": is_workflow,
                "is_detail": "n",
                "analysis": analysis,
                "asso_id": main_id_name,
                "table": table
            }
            self.db[coll_name].insert(data)
            for detail_table in one_ana["info"].keys():
                data = {
                    "is_first": "y" if detail_table in is_first_list else "n",
                    "is_workflow": is_workflow,
                    "is_detail": "y",
                    "analysis": analysis,
                    "asso_id": one_ana["info"][detail_table],
                    "table": detail_table
                }
                self.db[coll_name].insert(data)

    def add_record(self, record_file):
        pass

    def check_object_id(self, input_str):
        if not isinstance(input_str, ObjectId):
            if isinstance(input_str, types.StringTypes):
                input_str = ObjectId(input_str)
            else:
                raise Exception("main_id is no object id:%s" % input_str)
        return input_str


    def remove_detail(self, record_file=None):
        if record_file:
            file = open(record_file, "a")
            file.write("INPUT detail table remove record:\n")
        for one_main_coll in self.main_detail_map:
            main_coll_name = one_main_coll['name']
            main_id_name = one_main_coll['main_id_name']
            id_list = self.main_table_dic.setdefault(main_coll_name, [])
            for one_detail_coll in one_main_coll["info"].keys():
                this_detail_main_head = one_main_coll["info"][one_detail_coll]
                for main_coll_id in id_list:
                    result = self.db[main_coll_name].find_one({"_id": main_coll_id})
                    if result:
                        main_id = result[main_id_name]
                        if main_id_name != "_id" and main_id != main_coll_id:
                            continue
                        if record_file:
                            file.write("db.getCollection('%s').remove({%s: ObjectId('%s')})\n" % (one_detail_coll, this_detail_main_head, main_id))
                        else:
                            self.remove_one_detail(one_detail_coll, main_id, this_detail_main_head)
                    else:
                        print "主表%s找不到_id: %s" % (main_coll_name, main_coll_id)
        if record_file:
            file.close()

    def remove_one_detail(self, detail_table_name, main_table_id, main_head):
        main_table_id = self.check_object_id(main_table_id)
        result = self.db[detail_table_name].remove({main_head: main_table_id})
        rm_num = result['n']
        if rm_num == 0:
            print "TAKE CARE: detail coll: %s, find 0 records: {%s: ObjectId(\"%s\")}" % (detail_table_name, main_head, main_table_id)
        else:
            print "collection %s remove %s records" % (detail_table_name, rm_num)

    def remove_main(self, record_file=None):
        if record_file:
            file = open(record_file, "a")
            file.write("INPUT main table remove record:\n")
        for one_main_table in self.main_table_dic.keys():
            main_table_id_list = self.main_table_dic[one_main_table]
            if record_file:
                file.write("%s\t%s\n" % (one_main_table, main_table_id_list))
            else:
                self.db[one_main_table].remove({"_id": {"$in": main_table_id_list}})
                self.remove_record_db(main_table_id_list)
                self.remove_report_db(main_table_id_list)
                self.remove_status_db(main_table_id_list)
        if record_file:
            file.close()

    def rm_rubbish_main(self):
        # 与任务无关，只是清除掉所有没有参数的主表
        pass

    def rm_task(self, task_id, record_file=None):
        if record_file:
            file = open(record_file, "a")
            file.write("INPUT detail table remove record:\n")
        for one_main_coll in self.main_detail_map:
            main_coll_name = one_main_coll['name']
            main_id_name = one_main_coll['main_id_name']
            coll_result = self.db[main_coll_name].find({"task_id": task_id})
            print "%s find %s records" % (main_coll_name, coll_result.count())
            for one_detail_coll in one_main_coll["info"].keys():
                this_detail_main_head = one_main_coll["info"][one_detail_coll]
                for one_result in coll_result.clone():
                    try:
                        main_id = one_result[main_id_name]
                    except:
                        print "DEBUG %s please check main_id is %s or not,\n more_detail: %s" % (main_coll_name, main_id_name, one_result)
                        continue
                    if main_id_name != "_id" and one_result["_id"] != main_id:
                        print "collection:%s" % main_coll_name
                        print "main_id_name:%s != _id, _id %s != main_id %s" % (main_id_name, one_result["_id"], main_id)
                        continue
                    if record_file:
                        file.write("db.getCollection('%s').remove({%s: ObjectId('%s')})\n" % (one_detail_coll, this_detail_main_head, main_id))
                    else:
                        self.remove_one_detail(one_detail_coll, main_id, this_detail_main_head)
            if record_file:
                file.write("rm main_table %s, task_id %s\n" % (main_coll_name, task_id))
            else:
                result = self.db[main_coll_name].remove({"task_id": task_id})
                rm_num = result['n']
                if rm_num == 0:
                    print "TAKE CARE: main coll: %s, find 0 records: {task_id: '%s'}" % (main_coll_name, task_id)
                else:
                    print "collection %s remove %s records" % (main_coll_name, rm_num)
        table_list = self.rcf.get("other conditions", "have_no_detail").split(",") if self.rcf else []
        if record_file:
            for table in table_list + ["sg_task", "sg_result_table_deposit", "sg_status"]:
                if table:
                    file.write("db.getCollection('%s').remove({task_id: '%s'})\n" % (table, task_id))
            file.close()
        else:
            for table in table_list + ["sg_task", "sg_result_table_deposit", "sg_status"]:
                if table:
                    result = self.db[table].remove({"task_id": task_id})
                    rm_num = result['n']
                    if rm_num == 0:
                        print "TAKE CARE: main coll: %s, find 0 records: {task_id: '%s'}" % (table, task_id)
                    else:
                        print "collection %s remove %s records" % (table, rm_num)

    def rm_rubbish_detail(self, record_file=None):
        # 与任务无关，只是清除掉所有没有主表的明细表
        for one_main_coll in self.main_detail_map:
            main_coll_name = one_main_coll['name']
            main_id_name = one_main_coll['main_id_name']
            if record_file:
                file = open(record_file, "a")
                file.write("rum_rubbish_detail record:\n")
            else:
                file = None
            for one_detail_coll in one_main_coll["info"].keys():
                this_detail_main_head = one_main_coll["info"][one_detail_coll]
                self.check_and_remove(one_detail_coll, this_detail_main_head, main_coll_name, main_id_name, file)
            if record_file:
                file.close()

    def check_and_remove(self, detail_coll_name, detail_main, main_coll_name, main_id_name, record_file=None):
        result = self.db[detail_coll_name].find({})
        print_record = True
        if result.count() == 0:
            return
        tmp_rm_main_id = ""
        for one in result:
            try:
                main_id = one[detail_main]
            except:
                if print_record:
                    print_record = False
                    print "明细表%s没有字段%s" % (detail_coll_name, detail_main)
                continue
            if main_id == tmp_rm_main_id:
                continue
            else:
                tmp_rm_main_id = main_id
                main_result = self.db[main_coll_name].find_one({main_id_name: main_id})
                if not main_result:
                    if record_file:
                        record_file.write("db.getCollection('%s').remove({%s: ObjectId('%s')})\n" % (detail_coll_name, detail_main, main_id))
                    else:
                        result = self.db[detail_coll_name].remove({detail_main: main_id})
                        rm_num = result["n"]
                        print "collection %s remove %s records" % (detail_coll_name, rm_num)

    def remove_record_db(self, table_id_list):
        result = self.record_db.remove({"table_id": {"$in": table_id_list}, "task_id": self.task_id})
        print "remove record_db %s records" % result['n']

    def remove_report_db(self, result_table_id_list):
        new_list = [str(i) for i in result_table_id_list]
        result = self.report_db.remove({"result_table_id": {"$in": new_list}, "task_id": self.task_id})
        print "remove report_db %s records" % result['n']

    def remove_status_db(self, table_id_list):
        result = self.status_db.remove({"table_id": {"$in": table_id_list}, "task_id": self.task_id})
        print "remove status_db %s records" % result['n']

    def run(self, task_id, conf_table=None, conf_db=None, detail_log=None, main_log=None):
        if detail_log and not main_log:
            raise Exception("删主表前必须先删除明细表")
        if conf_table and conf_db:
            raise Exception("不能同时引用两种配置文件")
        if not conf_table and not conf_db:
            raise Exception("必须输入配置文件或对应mongo表")
        self.clean_all_var()
        self.get_main_table_info(task_id)
        if conf_table:
            self.get_detail_table_info(conf_table) # 根据配置文件获取主表与明细表对应关系，格式见方法说明
        elif conf_db:
            self.get_detail_table_info_from_mongo(coll_name=conf_db)
        # self.get_detail_table_info_from_mongo(collection_infos)  # 从某个mongo库获取主表与明细表对应关系，方法需自行定义
        self.remove_detail(record_file=detail_log)  # 先删detail再删main, 使用record参数不进行真正的删除，仅记录要删除的内容
        self.remove_main(record_file=main_log)
        # self.remove_detail()  # 真正的清除
        # self.remove_main()  # 真正的清除
        #################################
        #    HOW TO REMOVE ALL TASK
        #################################
        #  self.get_detail_table_info(conf_table, from_status=False)
        #  self.rm_task(task_id)
        ################################################################
        #   HOW TO REMOVE DETAIL TABLE WHICH MAIN TABLE DOESN'T EXISTS
        ################################################################
        #  self.get_detail_table_info(conf_table, from_status=False)
        #  self.rm_rubbish_detail()
        ################################################################

if __name__ == '__main__':
    pass
    # examples:
    # cl = clean_mongo(None, "metagenomic")
    # cl.run("tsg_32489", "metag_conf.txt", detail_log="/mnt/ilustre/users/sanger-dev/sg-users/guhaidong/debug/rm_record/tsg_32489_detail_table.log", main_log="/mnt/ilustre/users/sanger-dev/sg-users/guhaidong/debug/rm_record/tsg_32489_main_tabe.log")
    # cl = clean_mongo(None, "bacgenome")
    # cl.run("tsg_32645", "bacgenome_mongo_conf.txt", detail_log="tsg_32645_detail_table.log", main_log="tsg_32645_main_table.log")
    # cl.get_detail_table_info("bacgenome_mongo_conf.txt",from_status=False)
    # cl.rm_rubbish_detail(record_file="bacgenome_need_to_rm.txt")
    # cl.rm_task("tsg_32667", "tsg_32667_need_to_rm.txt")
    # cl = clean_mongo(None, "fungigenome")
    # cl.get_detail_table_info("fungigenome_mogo_conf.txt", from_status=False)
    # cl.rm_rubbish_detail(record_file="fungi_need_to_rm.txt")
    # cl.rm_task("tsg_30845", "tsg_30845_need_to_rm.txt")
    # cl.run("tsg_30845", "fungigenome_mogo_conf.txt", detail_log="tsg_30845_detail_table.log", main_log="tsg_30845_main_table.log")
    # m = clean_mongo(None, "meta")
    # m.get_detail_table_info("config/meta_conf.txt", from_status=False)
    # # m.insert_detail_table_info_to_mongo("table_relation", default_main_name="_id")
    # task_list = ["sanger_136735"]
    # for task_id in task_list:
    #     print "===============TASK ID %s===============" % task_id
    #     m.rm_task(task_id, "task/sanger_136735.txt")
    # m.insert_detail_table_info_to_mongo("table_relation", default_main_name="_id")
    # m.get_detail_table_info("metag_conf.txt", from_status=False)
    # m.rm_rubbish_detail(record_file="metag_need_to_rm.txt")
    # m.rm_task("sanger_136863", "task/sanger_136863_need_to_rm.txt")
    # m.run("sanger_136863", "config/metag_conf.txt", detail_log="jiaohu_log/sanger_136863_detail.log", main_log="jiaohu_log/sanger_136863_main.log")
    # b = clean_mongo(None, "metabolome")
    # b.get_detail_table_info_from_mongo(coll_name="table_relation")
    # b.rm_rubbish_detail("metab_need_to_rm.txt")
    # b.rm_task("tsg_31456", "tsg_31456_need_to_rm.txt")
    # b.run("tsg_32644", conf_db="table_relation", detail_log="tsg_32644_detail.log", main_log="tsg_32644_main.log")
