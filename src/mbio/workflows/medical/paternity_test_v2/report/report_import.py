# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class ReportImportWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        用途及运行逻辑：
        用于点击存入报告的时候，将用于出报告的数据存入到数据库中，并整合了该家系其他的一些信息，如果该家系的报告信息已经存入，
        就直接改变is_report为no，页面默认显示types为1的结果信息。同时要去更新下sg_family中的已经出报告与未出报告的字段。
        __auther__: hongdong
        last modified by hongdong@20171206
        :param wsheet_object:
        """
        self._sheet = wsheet_object
        super(ReportImportWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "dad_id", "type": "string"},
            {"name": "mom_id", "type": "string"},
            {"name": "son_id", "type": "string"},
            {"name": "update_info", "type": "string"},  # 用于后面更新主表
            {"name": "member_id", "type": "string"},  # 会员ID
            {"name": "father_err_id", "type": "string"},  # 旧的错配位点主表id
            {"name": "err_min", 'type': "int"},
            {"name": "report_time", "type": "string"},  # 报告中存入的日期
            {"name": "report_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.api_pt = self.api.api('medical.paternity_test_v2')

    def check_options(self):
        """检查参数设置"""
        if not self.option('dad_id'):
            raise OptionError('缺少参数dad_id')
        if not self.option('mom_id'):
            raise OptionError('缺少参数mom_id')
        if not self.option('son_id'):
            raise OptionError('缺少参数son_id')
        if not self.option('father_err_id'):
            raise OptionError('缺少参数father_err_id')
        if not self.option('err_min'):
            raise OptionError('缺少参数err_min')
        return True

    def run(self):
        """
        运行逻辑：
        1)检查F_M_S的信息在sg_father_summary中是否存在，存在的话，说明该家系的报告信息已经存入到库中了，
        那么该次我们就只要更新is_report为no，然后重新插入一个新的相关信息记录到库中，并定义is_report为yes
        2)如果检查F_M_S的信息在sg_father_summary中不存在，那么就,直接查表整合相关的信息到库中
        3)上述都完成后，要去更新下sg_family中的已经出报告与未出报告的字段
        :return:
        """
        self.start_listener()
        self.fire("start")
        name = '_'.join([self.option("dad_id"), self.option("mom_id"), self.option("son_id")])
        self.api_pt.check_report_result(name)
        self.api_pt.import_report_summary_result(self.option("dad_id"), self.option("mom_id"),
                                                 self.option("son_id"), self.option("father_err_id"), self.output_dir,
                                                 self.option('err_min'), self.option("report_time"), name,
                                                 self.option("report_id"))
        self.api_pt.update_report_done(self.option("dad_id"), self.option("father_err_id"))
        self.end()

    def end(self):
        super(ReportImportWorkflow, self).end()
