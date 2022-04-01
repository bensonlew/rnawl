# -*- coding: utf-8 -*-
# __author__ = 'liulinmeng'
import datetime
from biocluster.api.database.base import Base, report_check
from bson import SON


class BacgTaskInfo(Base):
    def __init__(self, bind_object):
        super(BacgTaskInfo, self).__init__(bind_object)
        self._project_type = 'bacgenome'

    def add_task_info(self, db_name=None):
        if db_name:
            self._db_name = db_name
        json_data = [
            ('task_id', self.bind_object.sheet.id),
            ('member_id', self.bind_object.sheet.member_id),
            ('project_sn', self.bind_object.sheet.project_sn),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('member_type', int(self.bind_object.sheet.member_type)),
            ('cmd_id', int(self.bind_object.sheet.cmd_id)),
            ('is_demo', 0),
            ('demo_id', self.bind_object.sheet.id)
        ]
        self.db['sg_task'].insert_one(SON(json_data))
        self.bind_object.logger.info('任务信息导入sg_task成功。')

    #@report_check
    def add_assem_task_info(self, db_name=None, task_type='mix'):
        """
        细菌基因组拼接工作流的任务导表
        :param db_name:
        :param task_type: 纯二代拼接draft, 纯三代拼接chr,均含mix
        :return:
        """
        self._project_type = "bac_assem"
        json_data = [
            ('task_id', self.bind_object.sheet.id),
            ('member_id', self.bind_object.sheet.member_id),
            ('project_sn', self.bind_object.sheet.project_sn),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('member_type', int(self.bind_object.sheet.member_type)),
            ('cmd_id', int(self.bind_object.sheet.cmd_id)),
            ('is_demo', 0),
            ('demo_id', self.bind_object.sheet.id),
            ('task_type', task_type)
        ]
        if self.bind_object.sheet.option("qc"):
            if self.bind_object.sheet.option("qc_tool") == "fastp":
                json_data.append(("qc_tool", "fastp"))
            else:
                json_data.append(("qc_tool", "Trimmomatic, SeqPrep, Sickle, FastqTotalHighQualityBase.jar;bamtools"))
        else:
            json_data.append(("qc_tool", "-"))
        if task_type == "draft":
            assem_tool = self.get_draft_assem_tool() + ", GapCloser"
            assem_para = "-" # 以后再加
        else:
            assem_tool = self.get_draft_assem_tool() + ", GapCloser"
            assem_para = "-"
        json_data.append(("assem_tool", assem_tool))
        json_data.append(("assem_para", assem_para))
        self.db['sg_task'].insert_one(SON(json_data))
        self.bind_object.logger.info('任务信息导入sg_task成功')

    def get_draft_assem_tool(self):
        if self.bind_object.option("pe_assem_tool") == "soapdenovo":
            return "SOAPdenovo"
        elif self.bind_object.option("pe_assem_tool") == "velvet":
            return "Velvet"
"""
    def get_chr_assem_para(self):
        if self.bind_object.option("css_assem_tool") == "canu":
            para = "error_rate %s" % self.bind_object.option("error_rate")
            para += ", cor_min_coverage %s" % self.bind_object.option("cor_min_coverage")
            para += ", cor_mhap_sensitivity %s" % self.bind_object.option("cor_mhap_sensitivity")
        else:
            para = "software default"
        return para
"""