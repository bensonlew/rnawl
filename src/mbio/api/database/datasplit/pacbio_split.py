# -*- coding: utf-8 -*-
# __author__ = 'Xue Qinwen'
# last_modify: 20190124

import os
import re
import json
import types
import datetime
from bson.objectid import ObjectId
from biocluster.api.database.base import Base, report_check

class PacbioSplit(Base):
    def __init__(self, bind_object):
        super(PacbioSplit, self).__init__(bind_object)
        self._project_type = 'datasplit'

    def update_sg_pacbio_specimen(self,split_id,s3_upload,stat_dir):
        print(1111)
        md5_dict = {}
        with open(os.path.join(os.path.dirname(stat_dir),"data/md5sum.txt"),"r") as md5file:
            while 1:
                line = md5file.readline()
                if not line:
                    break
                fd = line.rstrip().split("  ")
                md5_dict[fd[1]] = fd[0]
        for i in os.listdir(stat_dir):
            sample_name = i.split('_value.txt')[0]
            statfile = os.path.join(stat_dir,i)
            stat = {}
            with open(statfile,'r') as sf:
                line = sf.readline()
                title = line.rstrip().split('\t')
                rline = sf.readline()
                result = rline.rstrip().split('\t')
                for j in range(len(title)):
                    stat[title[j]] = result[j]
            query_dict = {
                "import_id": ObjectId(split_id),
                "majorbio_name": sample_name
                }
            if stat['atype'] == "diversity":
                update_dict = {
                    "reads": stat["reads"],
                    "qc_reads":stat["qc_reads"],
                    "reads_base": stat["reads_base"],
                    "total_len": stat["total_len"],
                    "ave_len": stat["ave_len"],
                    "raw_path": os.path.join(s3_upload,'PacbioQcStat/data/{}.ccs.fastq.gz'.format(sample_name)),
                    "clean_path": os.path.join(s3_upload,'PacbioQcStat/data/{}_value.fastq.gz'.format(sample_name)),
                    "raw_md5sum":md5_dict['{}.ccs.fastq.gz'.format(sample_name)],
                    "clean_md5sum":md5_dict['{}_value.fastq.gz'.format(sample_name)]
                }
            else:
                update_dict = {
                    "reads": stat["reads"],
                    "qc_reads":stat["qc_reads"],
                    "reads_base": stat["reads_base"],
                    "total_len": stat["total_len"],
                    "ave_len": stat["ave_len"],
                    "raw_path": os.path.join(s3_upload,'PacbioQcStat/data/{}.ccs.bam'.format(sample_name)),
                    "raw_md5sum":md5_dict['{}.ccs.bam'.format(sample_name)],
                    # "clean_path": os.path.join(s3_upload,'PacbioQcStat/{}_value.fastq'.format(sample_name)),
                }
            self.db["sg_pacbio_specimen"].update(query_dict, {"$set": update_dict})
    
    def update_sg_pacbio(self, split_id, stat_file):
        print(1111)
        stat = {}
        with open(stat_file,'r') as sf:
            line = sf.readline()
            title = line.rstrip().split('\t')
            rline = sf.readline()
            result = rline.rstrip().split('\t')
            for i in range(len(title)):
                stat[title[i]] = result[i]
        query_dict = {
                "_id": ObjectId(split_id),
                }
        update_dict = {
                    "subreads": stat["subreads_num"],
                    "ccs_reads":stat["ccsreads_num"],
                    "get_ratio": stat["get_ratio"],
                    "reconize_reads": stat["reconize_reads"],
                    "reconize_retio": stat["reconize_retio"],
                    "avg_len_ccs": stat["avg_len_ccs"]
                }
        self.db["sg_pacbio"].update(query_dict, {"$set": update_dict})

    def add_pacbio_bar(self, split_id, stat_file):
        bar_id = self.add_sg_bar(split_id,"ccs长度分布图","pacbio_split",[])
        with open(stat_file,'r') as sf:
            while 1:
                line = sf.readline()
                if not line:
                    break
                fd = line.rstrip().split('\t')
                self.add_sg_bar_detail(bar_id,fd[0],int(fd[1]))

    def add_sg_bar(self, split_id, name, location, categories, types=1):
        """
        导入柱状图数据
        origin_id: sg_split_clean_qc的_id
        name: 样本名称
        """
        print(1111)
        # if not isinstance(origin_id, ObjectId):
        #     if isinstance(origin_id, StringTypes):
        #         origin_id = ObjectId(origin_id)
        #     else:
        #         raise Exception("origin_id必须为ObjectId对象或其对应的字符串!")
        insert_data = {
            "task_id": ObjectId(split_id),
            "origin_id": ObjectId(split_id),
            "name": name,
            "categories": categories,
            "type": types,
            "location": location,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        return self.db["sg_bar"].insert_one(insert_data).inserted_id

    def add_sg_bar_detail(self, bar_id, name, value):
        """
        导入柱状图细节表
        """
        print(1111)
        data_list = []
        insert_data = {
            "bar_id": bar_id,
            "name": name,
            "value": value
        }
        data_list.append(insert_data)
        self.db["sg_bar_detail"].insert_many(data_list)
        
if __name__ == "__main__":
    a = PacbioSplit(None)
    split_id = "6136f57017b2bf11119ac9c3"
    s3_upload = "s3://rerewrweset/files/datasplit/2019/20190516sXten/5d09bed39dc6d6565858e0d4_20190621_131348"
    statread_new = "/mnt/ilustre/users/sanger-dev/workspace/20210917/PacbioSplit_SP20210915-1631687156_20210917_124744/PacbioStat/output/statReads.new.txt"
    stat_origin = "/mnt/ilustre/users/sanger-dev/workspace/20210917/PacbioSplit_SP20210915-1631687156_20210917_124744/PacbioStat/output/m64236_210806_111052"
    statistic = "/mnt/ilustre/users/sanger-dev/workspace/20210915/Single_qc_stat_233/PacbioQcStat/output/statistic"
    a.update_sg_pacbio_specimen(split_id,s3_upload,statistic)
    a.update_sg_pacbio(split_id,stat_origin)
    a.add_pacbio_bar(split_id,statread_new)
    