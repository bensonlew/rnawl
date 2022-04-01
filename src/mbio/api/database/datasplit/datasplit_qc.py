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

class DatasplitQc(Base):
    def __init__(self, bind_object):
        super(DatasplitQc, self).__init__(bind_object)
        self._project_type = 'datasplit'

    def update_sg_qc_specimen(self,qc_id,fastq_dir,s3_upload_dir):
        qc_id = ObjectId(qc_id)
        fq_list = {}
        md5sum_info = {}
        sample_list = []
        sample = ""
        md5sum_file = os.path.join(fastq_dir,"md5sum.txt")
        md5sum_dict = {}
        if os.path.exists(md5sum_file):
            with open(md5sum_file,"r") as mf:
                while 1:
                    line = mf.readline()
                    if not line:
                        break
                    fd = line.rstrip().split("  ")
                    md5sum_dict[fd[1]] = fd[0]
        for f in os.listdir(fastq_dir):
            if f.endswith(".fq.gz") or f.endswith(".fastq.gz") or f.endswith(".fastq"):
                    # lib = f.split("--")[0]
                sample = f.split("--")[1].split(".")[0]
                if sample not in sample_list:
                    sample_list.append(sample)
                    md5sum_info[sample] = []
                    fq_list[sample] = []
                fa_path = s3_upload_dir+f if s3_upload_dir.endswith("/") else s3_upload_dir+"/"+f
                    # fa_work_path = os.path.join(fastq_dir, f)
                    # clean_bytes = str(os.path.getsize(os.path.join(fastq_dir, f)))
                md5sum_info[sample].append(md5sum_dict[f])
                fq_list[sample].append(fa_path)
        for sample in sample_list:
            if sample != "":
                query_dict = {"qc_id": qc_id, "specimen_name": sample}
                update_dict = {"clean_path": ";".join(fq_list[sample]),"clean_md5sum":";".join(md5sum_info[sample])}
                self.db["sg_qc_specimen"].update(query_dict, {"$set": update_dict})
    
    # def update_sg_qc(self, qc_id):
    #     print(1111)
    #     stat = {}
    #     query_dict = {
    #             "_id": ObjectId(split_id),
    #             }
    #     update_dict = {
                    
    #             }
    #     self.db["sg_pacbio"].update(query_dict, {"$set": update_dict})

    
if __name__ == "__main__":
    a = DatasplitQc(None)
    split_id = "6136f57017b2bf11119ac9c3"
    s3_upload = "s3://rerewrweset/files/datasplit/2019/20190516sXten/5d09bed39dc6d6565858e0d4_20190621_131348"
    statread_new = "/mnt/ilustre/users/sanger-dev/workspace/20210917/DatasplitQc_SP20210915-1631687156_20210917_124744/PacbioStat/output/statReads.new.txt"
    stat_origin = "/mnt/ilustre/users/sanger-dev/workspace/20210917/DatasplitQc_SP20210915-1631687156_20210917_124744/PacbioStat/output/m64236_210806_111052"
    statistic = "/mnt/ilustre/users/sanger-dev/workspace/20210915/Single_qc_stat_233/PacbioQcStat/output/statistic"
    a.update_sg_pacbio_specimen(split_id,s3_upload,statistic)
    a.update_sg_pacbio(split_id,stat_origin)
    a.add_pacbio_bar(split_id,statread_new)
    