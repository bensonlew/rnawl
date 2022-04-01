# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.26

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
from biocluster.api.file.lib.transfer import MultiFileTransfer
import json


class VennAnalysisAgent(Agent):
    """
    WGS项目对比较分析的结果做venn分析
    """
    def __init__(self, parent):
        self._sheet = parent
        super(VennAnalysisAgent, self).__init__(parent)
        options = [
            {"name": "set_ids", "type": "string"},
            {"name": "compare_type", "type": "string"},  # 比较分析结果的类型,snp/indel
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "project_type", "type": "string"},
            {'name': "is_bucket", "type": "string", 'default': "false"},
            {"name": "rerewrwese_path", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("set_ids"):
            raise OptionError("请设置set_ids", code="34507601")
        if not self.option("compare_type"):
            raise OptionError("请设置比较分析的类型", code="34507602")
        if self.option("compare_type") not in["snp", "indel"]:
            raise OptionError("compare_type类型必须为snp/indel", code="34507603")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(VennAnalysisAgent, self).end()


class VennAnalysisTool(Tool):
    def __init__(self, config):
        super(VennAnalysisTool, self).__init__(config)
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.venn_path = self.config.PACKAGE_DIR + "/wgs/venn.pl"

    def get_variant_result_list(self):
        """
        得到每个set_id对应的variant结果文件
        """
        venn_api = self.api.api("wgs.venn_analysis")
        set_ids = self.option("set_ids").split(",")
        type = self.option("compare_type")
        self.diff_list = os.path.join(self.work_dir, "diff_result.list")
        self.logger.info(set_ids)
        # data = os.path.join(os.path.abspath(os.path.join(self.work_dir, "..")), "data.json")
        # f = open(data, "r")
        # data_json = json.loads(f.read())
        # if data_json["client"] == 'client01':
        #     rerewrwese_path = "/mnt/ilustre/data"
        # else:
        #     rerewrwese_path = "/mnt/ilustre/tsanger-data"
        rerewrwese_path = self.option("rerewrwese_path")
        if self.option("project_type"):
            venn_api._project_type = self.option("project_type")
        venn_api.get_variant_result_list(set_ids, type, self.diff_list, rerewrwese_path)
        if self.option('is_bucket') == 'true':
            self.download_from_s3_()

    def download_from_s3_(self):
        """
        从对象存储中下载文件到指定路径
        :return:
        """
        if not os.path.exists(self.work_dir + "/temp"):
            os.mkdir(self.work_dir + "/temp")
        if os.path.exists(self.work_dir + "/new_diff.list"):
            os.remove(self.work_dir + "/new_diff.list")
        self.logger.info("开始下载对象存储中的文件！")
        transfer = MultiFileTransfer()
        with open(self.diff_list, 'r') as r, open(self.work_dir + "/new_diff.list", 'w') as w:
            data = r.readlines()
            for line in data:
                tmp = line.strip().split('\t')
                transfer.add_download(tmp[1], "{}/temp/{}/".format(self.work_dir, tmp[0]))
                w.write("{}\t{}\n".format(tmp[0], self.work_dir + "/temp/{}/{}"
                                          .format(tmp[0], os.path.basename(tmp[1]))))
        transfer.perform()
        self.logger.info("下载对象存储中的文件成功！")

    def run_venn(self):
        """
        venn.pl
        """
        venn_result = os.path.join(self.work_dir, "venn_output")
        if not os.path.exists(venn_result):
            os.mkdir(venn_result)
        infile = self.work_dir + "/new_diff.list" if self.option('is_bucket') == 'true' else self.diff_list
        cmd = "{} {} -pop {} -out {}".format(self.perl_path, self.venn_path, infile, venn_result)
        command = self.add_command("venn_analysis", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("venn_analysis完成")
        else:
            self.set_error("venn_analysis失败", code="34507601")

    def set_output(self):
        venn_result = os.path.join(self.work_dir, "venn_output")
        for f in os.listdir(venn_result):
            if f == "data.xls":
                os.link(os.path.join(venn_result, "data.xls"), os.path.join(self.output_dir, "venn_stat.xls"))
            else:
                os.link(os.path.join(venn_result, f), os.path.join(self.output_dir, f))

    def set_db(self):
        """
        保存结果到mongo
        """
        venn_api = self.api.api("wgs.venn_analysis")
        if self.option("project_type"):
            venn_api._project_type = self.option("project_type")
        main_id = self.option("main_id")
        set_ids = self.option("set_ids").split(",")
        table_name, table_title = [], []
        i = 0
        file_dir = os.listdir(self.output_dir)
        for f in file_dir:
            if f == "venn_stat.xls":
                venn_stat = os.path.join(self.output_dir, "venn_stat.xls")
                venn_api.add_sg_venn_analysis_stat(main_id, venn_stat)
            else:
                name = f.split(".result")[0]
                table_name.append({"name" + str(i): name})
                venn_detail = os.path.join(self.output_dir, f)
                title = venn_api.add_sg_venn_analysis_detail(main_id, "name" + str(i), venn_detail)
                table_title.append(title)
                i += 1
        venn_api.update_table_name(main_id, {"table_name": table_name})
        venn_api.update_table_name(main_id, {"table_title": table_title})
        # venn_id = venn_api.add_sg_venn(main_id)  # 取消venn图
        # venn_api.add_sg_venn_detail(venn_id, set_ids, self.diff_list)

    def run(self):
        super(VennAnalysisTool, self).run()
        self.get_variant_result_list()
        self.run_venn()
        self.set_output()
        self.set_db()
        self.end()
