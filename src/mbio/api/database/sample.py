# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
# last modified guhaidong 20171116
from biocluster.api.database.base import Base, report_check
import re
# from biocluster.config import Config


class Sample(Base):
    def __init__(self, bind_object):
        super(Sample, self).__init__(bind_object)
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB
        self.sample_table_ids = list()
        self.main_task_id = "_".join(self.bind_object.sheet.id.split('_')[0:2])  # add main_task_id by guhaidong 20171116

    @report_check
    def add_samples_info(self, file_path):
        self.bind_object.logger.info("add_samples_info start")
        data_list = []
        with open(file_path, 'r') as f:
            l = f.readline()
            if not re.match(r"sample\s+reads\s+bases\s+avg\s+min\s+max", l):
                self.bind_object.logger.error("文件%s格式不正确，请选择正确的样品信息文件" % file_path)
                self.bind_object.set_error("结果文件格式不正确", code="51007001")
            while True:
                line = f.readline().strip('\n')
                if not line:
                    break
                line_data = line.split("\t")
                data = {
                    "project_sn": self.bind_object.sheet.project_sn,
                    "task_id": self.bind_object.sheet.id,
                    "specimen_name": line_data[0],
                    "read_number": int(line_data[1]),
                    "base_number": int(line_data[2]),
                    "average_length": float(line_data[3]),
                    "min_length": int(line_data[4]),
                    "max_length": int(line_data[5]),
                    "is_initial": 1,
                    # by houshuang 20190815>>>
                    "new_name": line_data[0],
                    "desc": '-'
                    # <<<
                }
                data_list.append(data)
        try:
            collection = self.db["sg_specimen"]
            result = collection.insert_many(data_list)
            self.sample_table_ids = result.inserted_ids[:]
        except Exception, e:
            self.bind_object.logger.error("导入样品信息数据出错:%s" % e)
            self.bind_object.set_error("导入样品信息数据出错", code="51007002")
        else:
            self.bind_object.logger.info("导入样品信息数据成功:%s" % result.inserted_ids)
        return self.sample_table_ids

    @report_check
    def add_base_info(self, sample_name, file_path):
        collection = self.db["sg_specimen"]
        results = collection.find({"specimen_name": sample_name, "task_id": self.main_task_id})  # add task_id by guhaidong 20171116
        if results.count():
            specimen_id = self._find_specimen_id(results)
            if not specimen_id:
                self.bind_object.set_error("没有找到对应的样品%s信息，请先导入样品信息表!", variables=(sample_name), code="51007003")
        else:
            self.bind_object.set_error("没有找到对应的样品%s信息，请先导入样品信息表!", variables=(sample_name), code="51007003")
        data_list = []
        with open(file_path, 'r') as f:
            l = f.readline()
            if not re.match(r"column\tcount\tmin\tmax\tsum\tmean\tQ1\tmed\tQ3\tIQR\tlW\trW\t"
                            r"A_Count\tC_Count\tG_Count\tT_Count\tN_Count\tMax_count", l):
                self.bind_object.logger.error("文件%s格式不正确，请选择正确的碱基统计文件" % file_path)
                self.bind_object.set_error("碱基统计文件格式不正确", code="51007004")
            while True:
                line = f.readline().strip('\n')
                if not line:
                    break
                line_data = line.split("\t")
                data = {
                    "project_sn": self.bind_object.sheet.project_sn,
                    "task_id": self.bind_object.sheet.id,
                    "specimen_id": specimen_id,
                    "column": int(line_data[0]),
                    "min": int(line_data[2]),
                    "max": int(line_data[3]),
                    "q1": int(line_data[6]),
                    "q3": int(line_data[8]),
                    "median": int(line_data[7]),
                    "average": round(float(line_data[5]), 2),
                    "lw": int(line_data[10]),
                    "rw": int(line_data[11])
                }
                data_list.append(data)
        try:
            collection = self.db["sg_specimen_sequence"]
            collection.insert_many(data_list)
            main_collection = self.db["sg_specimen"]
            #main_collection.update({"_id": specimen_id},{"$set": { "main_id": specimen_id}})

        except Exception, e:
            self.bind_object.logger.error("导入样品%s的碱基统计信息出错:%s" % (sample_name, e))
            self.bind_object.set_error("导入样品的碱基统计信息出错", code="51007005")
        else:
            self.bind_object.logger.info("导入样品%s的碱基统计信息成功" % sample_name)

    @report_check
    def add_reads_len_info(self, step_length, file_path):
        data_list = []
        all_specimen_id = []
        with open(file_path, 'r') as f:
            l = f.readline()
            if not re.match(r"^sample\t", l):
                self.bind_object.logger.error("文件%s格式不正确，请选择正确的碱基统计文件" % file_path)
                self.bind_object.set_error("碱基统计文件格式不正确", code="51007004")
            else:
                length_list = l.strip("\n").split("\t")
                length_list.pop(0)
            while True:
                line = f.readline().strip('\n')
                if not line:
                    break
                line_data = line.split("\t")
                collection = self.db["sg_specimen"]
                results = collection.find({"specimen_name": line_data[0], "task_id": self.main_task_id})  # add task_id by guhaidong 20171116
                if results.count():
                    specimen_id = self._find_specimen_id(results)
                    if not specimen_id:
                        self.bind_object.set_error("没有找到对应的样品%s信息，请先导入样品信息表!", variables=(line_data[0]), code="51007003")
                else:
                    self.bind_object.set_error("没有找到对应的样品%s信息，请先导入样品信息表!", variables=(line_data[0]), code="51007003")
                step_data = {}
                i = 0
                for step in length_list:
                    i += 1
                    step_data[step] = int(line_data[i])
                data = {
                    "project_sn": self.bind_object.sheet.project_sn,
                    "task_id": self.bind_object.sheet.id,
                    "specimen_id": specimen_id,
                    "step": step_length,
                    "value": step_data
                }
                data_list.append(data)
                all_specimen_id.append(str(specimen_id))
        params = all_specimen_id
        main_info = dict(
            project_sn=self.bind_object.sheet.project_sn,
            task_id=self.bind_object.sheet.id,
            name='sg_specimen_step_' + self.bind_object.sheet.id,
            desc='sg_specimen_step',
            params=params,
            status="end",
            step=20
        )
        try:
            collection = self.db["sg_specimen_step"]
            collection.insert_many(data_list)
            if step_length == 20:
                main = self.db["sg_specimen_step_main"]
                main.insert_many([main_info])
        except Exception, e:
            self.bind_object.logger.error("导入步%s的步长序列长度统计出错:%s" % (step_length, e))
            self.bind_object.set_error("步长序列长度统计出错", code="51007006")
        else:
            self.bind_object.logger.info("导入步%s的步长序列长度统计成功" % step_length)

    @report_check
    def get_spname_spid(self):
        if not self.sample_table_ids:
            self.bind_object.logger.error("样本id列表为空，请先调用add_samples_info产生sg_speciem的id")
            self.bind_object.set_error("没有找到样本信息", code="51007007")
        collection = self.db["sg_specimen"]
        spname_spid = dict()
        for id_ in self.sample_table_ids:
            result = collection.find_one({"_id": id_, "task_id": self.main_task_id})  # add task_id by guhaidong 20171116
            if not result:
                self.bind_objecct.set_error("意外错误，无法找到样本id", code="51007008")
            spname_spid[result["specimen_name"]] = id_
        return spname_spid

    def _find_specimen_id(self, results):
        specimen_id = ""
        for result in results:
            if result["_id"] in self.sample_table_ids:
                specimen_id = result["_id"]
                break
        return specimen_id
        
    @report_check
    def add_valid_sequence_info(self, valid_sequence_path):
        self.bind_object.logger.info("add_valid_sequence_info start")
        data_list = []
        with open(valid_sequence_path, 'r') as f:
            #f.readline()
            line = f.readline()
            line_data = line.strip().split("\t")
            data = {
                "project_sn": self.bind_object.sheet.project_sn,
                "task_id": self.bind_object.sheet.id,
                "amplified_region": line_data[0],
                "samples": line_data[1],
                "sequences": line_data[2],
                "bases": line_data[3],
                "average_length": line_data[4]
            }
            data_list.append(data)
        try:
            collection = self.db["sg_valid_sequence_info"]
            result = collection.insert_many(data_list)
            self.sample_table_ids = result.inserted_ids[:]

            # for table_id in self.sample_table_ids:
            #     collection.update({"_id": table_id},{"$set" : {"main_id": table_id}})

        except Exception, e:
            self.bind_object.logger.error("导入样品有效信息数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入样品有效信息数据成功:%s" % result.inserted_ids)
    
    @report_check
    def add_raw_sequence_info(self, raw_sequence_path):
        print("add_raw_sequence_info start")
        data_list = []
        with open(raw_sequence_path, 'r') as f:
            f.readline()
            line = f.readline()
            line_data = line.strip().split("\t")
            data = {
                "project_sn": self.bind_object.sheet.project_sn,
                "task_id": self.bind_object.sheet.id,
                "amplified_region": line_data[0],
                "insert_size": line_data[1],
                "sequencing_length": line_data[2],
                "raw_reads": line_data[3],
                "total_base": line_data[4]
            }
            data_list.append(data)
        try:
            collection = self.db["sg_raw_sequence_info"]
            result = collection.insert_many(data_list)
            self.sample_table_ids = result.inserted_ids[:]

            # for table_id in self.sample_table_ids:
            #     collection.update({"_id": table_id},{"$set" : {"main_id": table_id}})

        except Exception, e:
            print("导入样品原始信息数据出错:%s" % e)
        else:
            print("导入样品原始信息数据成功:%s" % result.inserted_ids)
