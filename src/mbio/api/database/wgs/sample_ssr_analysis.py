# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.26

from api_base import ApiBase
import datetime
import os


class SampleSsrAnalysis(ApiBase):
    def __init__(self, bind_object):
        """
        WGS SSR分析导表
        """
        super(SampleSsrAnalysis, self).__init__(bind_object)
        self._project_type = "dna_wgs"

    def add_sg_ssr_specimen(self, project_sn, task_id, params=None, name=None):
        """
        sg_ssr_specimen
        """
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_ssr_specimen",
            "params": params if params else "null",
            "desc": "SSR样本基因组主表"
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_ssr_specimen", data_list)
        self.update_db_record("sg_ssr_specimen", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_ssr_specimen_stat(self, ssr_specimen_id, specimen_id, ssr_stat):
        """
        sg_ssr_specimen_stat
        ssr_stat: ssr.stat
        """
        ssr_specimen_id = self.check_objectid(ssr_specimen_id)
        self.check_exists(ssr_stat)
        data_list = []
        with open(ssr_stat, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "ssr_specimen_id": ssr_specimen_id,
                    "specimen_id": specimen_id,
                    "chr": item[0],
                    "ssr_num": int(item[1]),
                    "c": int(item[2]),
                    "c_star": int(item[3]),
                    "p1": int(item[4]),
                    "p2": int(item[5]),
                    "p3": int(item[6]),
                    "p4": int(item[7]),
                    "p5": int(item[8]),
                    "p6": int(item[9]),
                }
                data_list.append(insert_data)
        self.col_insert_data("sg_ssr_specimen_stat", data_list)

    def add_sg_ssr_specimen_detail(self, ssr_specimen_id, specimen_id, ssr_detail):
        """
        sg_ssr_specimen_detail
        ssr_detail: specimen.result
        """
        ssr_specimen_id = self.check_objectid(ssr_specimen_id)
        self.check_exists(ssr_detail)
        data_list = []
        # num = 0
        with open(ssr_detail, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                try:
                    insert_data = {
                        "ssr_specimen_id": ssr_specimen_id,
                        "specimen_id": specimen_id,
                        "chr": item[0],
                        "ssr_type": item[2],
                        "ssr": item[3],
                        "start": int(item[5]),
                        "end": int(item[6]),
                        "forward_primer": item[7],
                        "forward_tm": item[8],
                        "forward_gc": item[9],
                        "forward_len": item[10],
                        "reverse_primer": item[11],
                        "reverse_tm": item[12],
                        "reverse_gc": item[13],
                        "reverse_len": item[14],
                        "product_size": item[15]
                    }
                except:
                    insert_data = {
                        "ssr_specimen_id": ssr_specimen_id,
                        "specimen_id": specimen_id,
                        "chr": item[0],
                        "ssr_type": item[2],
                        "ssr": item[3],
                        "start": int(item[5]),
                        "end": int(item[6]),
                        "forward_primer": "-",
                        "forward_tm": "-",
                        "forward_gc": "-",
                        "forward_len": "-",
                        "reverse_primer": "-",
                        "reverse_tm": "-",
                        "reverse_gc": "-",
                        "reverse_len": "-",
                        "product_size": "-"
                    }
                data_list.append(insert_data)
                # insert_data = {
                #     "ssr_specimen_id": ssr_specimen_id,
                #     "specimen_id": specimen_id,
                #     "chr": item[0],
                #     "ssr_nr": item[1],
                #     "ssr_type": item[2],
                #     "ssr": item[3],
                #     "size": item[4],
                #     "start": int(item[5]),
                #     "end": int(item[6])
                # }
                # origin_id = self.db["sg_ssr_specimen_detail"].insert_one(insert_data).inserted_id
                # try:
                #     num = len(item[7].split(";"))
                #     for i in range(num):
                #         primer_data = {
                #             "ssr_id": ssr_specimen_id,
                #             "origin_id": origin_id,
                #             "chr": item[0],
                #             "ssr_type": item[2],
                #             "size": item[4],
                #             "start": int(item[5]),
                #             "end": int(item[6]),
                #             "forward_primer": item[7].split(";")[i],
                #             "forward_tm": float(item[8].split(";")[i]),
                #             "forward_gc": float(item[9].split(";")[i]),
                #             "forward_len": int(item[10].split(";")[i]),
                #             "reverse_primer": item[11].split(";")[i],
                #             "reverse_tm": float(item[12].split(";")[i]),
                #             "reverse_gc": float(item[13].split(";")[i]),
                #             "reverse_len": int(item[14].split(";")[i]),
                #             "product_size": int(item[15].split(";")[i])
                #         }
                #         data_list.append(primer_data)
                #         # self.db["sg_ssr_primer"].insert_one(primer_data)
                # except:
                #     for i in range(num):
                #         primer_data = {
                #             "ssr_id": ssr_specimen_id,
                #             "origin_id": origin_id,
                #             "chr": item[0],
                #             "ssr_type": item[2],
                #             "size": item[4],
                #             "start": int(item[5]),
                #             "end": int(item[6]),
                #             "forward_primer": "-",
                #             "forward_tm": "-",
                #             "forward_gc": "-",
                #             "forward_len": "-",
                #             "reverse_primer": "-",
                #             "reverse_tm": "-",
                #             "reverse_gc": "-",
                #             "reverse_len": "-",
                #             "product_size": "-"
                #         }
                #         data_list.append(primer_data)
                #         # self.db["sg_ssr_primer"].insert_one(primer_data)
        # self.col_insert_data("sg_ssr_primer", data_list)
        self.col_insert_data("sg_ssr_specimen_detail", data_list)

    def add_sg_ssr_specimen_detail_new(self, ssr_specimen_id, specimen_id, ssr_detail, download_file):
        """
        sg_ssr_specimen_detail
        ssr_detail: specimen.result
        """
        ssr_specimen_id = self.check_objectid(ssr_specimen_id)
        self.check_exists(ssr_detail)
        data_list = []
        self.update_db_record("sg_ssr_specimen", {"main_id": ssr_specimen_id}, {"download_path": download_file})
        with open(ssr_detail, "r") as f:
            head = f.readline().strip().split("\t")
            num = len(head[6:]) / 9
            for line in f:
                item = line.strip().split("\t")
                forward_primer, forward_tm, forward_gc, forward_len = [], [], [], []
                reverse_primer, reverse_tm, reverse_gc, reverse_len = [], [], [], []
                product_size = []
                for i in range(num):
                    j = 9 * i
                    try:
                        forward_primer.append(item[7+j])
                        forward_tm.append(item[8+j])
                        forward_gc.append(item[9+j])
                        forward_len.append(item[10+j])
                        reverse_primer.append(item[11+j])
                        reverse_tm.append(item[12+j])
                        reverse_gc.append(item[13+j])
                        reverse_len.append(item[14+j])
                        product_size.append(item[15+j])
                    except:
                        break
                insert_data = {
                    "ssr_specimen_id": ssr_specimen_id,
                    "specimen_id": specimen_id,
                    "chr": item[0],
                    "ssr_nr": item[1],
                    "ssr_type": item[2],
                    "ssr": item[3],
                    "size": item[4],
                    "start": int(item[5]),
                    "end": int(item[6]),
                    "forward_primer": ";".join(forward_primer),
                    "forward_tm": ";".join(forward_tm),
                    "forward_gc": ";".join(forward_gc),
                    "forward_len": ";".join(forward_len),
                    "reverse_primer": ";".join(reverse_primer),
                    "reverse_tm": ";".join(reverse_tm),
                    "reverse_gc": ";".join(reverse_gc),
                    "reverse_len": ";".join(reverse_len),
                    "product_size": ";".join(product_size)
                }
                data_list.append(insert_data)
        self.col_insert_data("sg_ssr_specimen_detail", data_list)

    def get_specimen_fastq_list(self, task_id, specimen_names, fastq_list, target_dir=None):
        """
        SSR交互分析，根据task_id、specimen_names得到fastq_list
        """
        print task_id
        # result = self.col_find_one("sg_task", query_dic={"task_id": task_id})
        # fastq_dir = result["clean_fastq_path"]
        # self.bind_object.logger.info(fastq_dir)
        # if fastq_dir.startswith("rerewrweset"):
        #     fastq_dir = os.path.join(target_dir, fastq_dir)
        with open(fastq_list, "w") as w:
            for s in specimen_names.split(","):
                s_results = self.col_find("sg_specimen", query_dic={"task_id": task_id, "old_name": s})
                for s_result in s_results:
                    clean_dir = s_result["clean_path"].split("|")
                    for clean_path in clean_dir:
                        try:
                            fastq_l = os.path.join(target_dir, clean_path.split(",")[0])
                            fastq_r = os.path.join(target_dir, clean_path.split(",")[1])
                            # if not os.path.exists(fastq_l):   # 这里要根据具体的s3文件检查
                            #     raise Exception("文件{}不存在，请检查".format(fastq_l))
                            # if not os.path.exists(fastq_r):
                            #     raise Exception("文件{}不存在，请检查".format(fastq_r))
                            w.write(s + "\t" + fastq_l + "\t" + fastq_r + "\n")
                        except:
                            raise Exception("任务{}的sg_specimen表没有clean_path信息，请检查".format(task_id))

    def update_ssr_output_dir(self, ssr_specimen_id, output_dir, download_file):
        """
        更新样本基因组SSR分析的output_dir，用于样本基因组SSR比较分析
        """
        ssr_specimen_id = self.check_objectid(ssr_specimen_id)
        self.update_db_record("sg_ssr_specimen", {"_id": ssr_specimen_id}, {"output_dir": output_dir, "download_path": download_file})


if __name__ == "__main__":
    a = SampleSsrAnalysis(None)
    project_sn = 'wgs_test'
    task_id = 'wgs_test'
    # ssr_specimen_id = a.add_sg_ssr_specimen(project_sn, task_id)
    # # ssr_specimen_id = "5aea6131a4e1af344973ad89"
    # specimen_id = "GC_bulk"
    # ssr_stat = "/mnt/ilustre/users/sanger-dev/workspace/20180502/SampleSsr_wgs_test_0502181855_3886_6295/SsrPrimer/output/GC_bulk.ssr.stat"
    # ssr_detail = "/mnt/ilustre/users/sanger-dev/workspace/20180502/SampleSsr_wgs_test_0502181855_3886_6295/SsrPrimer/output/GC_bulk.result"
    # a.add_sg_ssr_specimen_stat(ssr_specimen_id, specimen_id, ssr_stat)
    # a.add_sg_ssr_specimen_detail(ssr_specimen_id, specimen_id, ssr_detail)
    # specimen_id = "JY102"
    # ssr_stat = "/mnt/ilustre/users/sanger-dev/workspace/20180502/SampleSsr_wgs_test_0502181855_3886_6295/SsrPrimer/output/JY102.ssr.stat"
    # ssr_detail = "/mnt/ilustre/users/sanger-dev/workspace/20180502/SampleSsr_wgs_test_0502181855_3886_6295/SsrPrimer/output/JY102.result"
    # a.add_sg_ssr_specimen_stat(ssr_specimen_id, specimen_id, ssr_stat)
    # a.add_sg_ssr_specimen_detail(ssr_specimen_id, specimen_id, ssr_detail)
    ssr_specimen_id = "5b1a4f61ffec6054b5b648c0"
    # specimen_id = "yellow"
    # ssr_stat = "/mnt/ilustre/users/sanger-test/workspace/20180604/SampleSsr_tsanger_30180_0604134221570747_3942/output/yellow.ssr.stat"
    # ssr_detail = "/mnt/ilustre/users/sanger-test/workspace/20180604/SampleSsr_tsanger_30180_0604134221570747_3942/output/yellow.result"
    # a.add_sg_ssr_specimen_detail(ssr_specimen_id, specimen_id, ssr_detail)
    # a.add_sg_ssr_specimen_stat(ssr_specimen_id, specimen_id, ssr_stat)
    specimen_id = "dark1"
    # ssr_stat = "/mnt/ilustre/users/sanger-test/workspace/20180604/SampleSsr_tsanger_30180_0604134221570747_3942/output/dark1.ssr.stat"
    # ssr_detail = "/mnt/ilustre/users/sanger-test/workspace/20180604/SampleSsr_tsanger_30180_0604134221570747_3942/output/dark1.result"
    # a.add_sg_ssr_specimen_detail(ssr_specimen_id, specimen_id, ssr_detail)
    # a.add_sg_ssr_specimen_stat(ssr_specimen_id, specimen_id, ssr_stat)
    # ssr_detail = "/mnt/ilustre/users/sanger-test/workspace/20180608/SampleSsr_tsanger_30180_0608174153706508_7679/SsrPrimer/output/yellow_result.xls"
    download_file = "rerewrweset/files/m_188/188_5b03d16580da8/tsanger_30180/interaction_results/sg_ssr_specimen"
    ssr_detail = "/mnt/ilustre/users/sanger-test/workspace/20180608/SampleSsr_tsanger_30180_0608174153706508_7679/SsrPrimer/output/dark1_result.xls"
    a.add_sg_ssr_specimen_detail_new(ssr_specimen_id, specimen_id, ssr_detail, download_file)
