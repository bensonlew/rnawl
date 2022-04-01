# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180409

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class SampleCompareWorkflow(Workflow):
    """
    交互分析：样本比较分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SampleCompareWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "sample1", "type": "string"},
            {"name": "sample2", "type": "string"},
            {"name": "is_same", "type": "string", "default": "true"},  # 页面样本间拷贝数变异是否相同
            {"name": "vcf_file", "type": "infile", 'format': 'bsa.vcf'},
            {"name": "funtype", "type": "string"},  # HIGH,LOW,MODIFIER,MODERATE
            {"name": "efftype", "type": "string"},
            {"name": "len1", "type": "string", "default": "1"},  # snp的时候默认len1与len2是1
            {"name": "len2", "type": "string", "default": "1"},
            {"name": "dep1", "type": "string"},  # 2,5
            {"name": "dep2", "type": "string"},  # 2,6
            {"name": "location", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "project_type", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.cnv_diff = self.add_tool("wgs.sample_compare")

    def check_options(self):
        if not self.option("sample1"):
            raise OptionError("缺少样本1参数", code="14500801")
        if not self.option("sample2"):
            raise OptionError("缺少样本2参数", code="14500802")
        if not self.option("is_same"):
            raise OptionError("缺少is_same参数", code="14500803")
        if not self.option("vcf_file"):
            raise OptionError("缺少vcf_file参数", code="14500804")
        return True

    def run_sample_compare(self):
        options = {
            "sample1": self.option("sample1"),
            "sample2": self.option("sample2"),
            "is_same": self.option("is_same"),
            "vcf_file": self.option("vcf_file").prop['path'],
            "funtype": self.option("funtype"),
            "efftype": self.option("efftype"),
            "len1": self.option("len1"),
            "len2": self.option("len2"),
            "dep1": self.option("dep1"),
            "dep2": self.option("dep2"),
            "location": self.option("location"),
        }
        self.cnv_diff.set_options(options)
        self.cnv_diff.on("end", self.set_output, "sample_compare")
        self.cnv_diff.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'sample_compare':
            self.linkdir(obj.output_dir, 'sample_compare')
        self.set_db()

    def linkdir(self, dirpath, dirname):
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
                    # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_db(self):
        self.logger.info("保存结果到mongo")
        api_path = self.api.api("wgs.snp_indel_compare")
        if self.option("project_type"):
            api_path._project_type = self.option("project_type")
        types = "snp" if self.option("submit_location") == "snpcompare_specimen" else "indel"
        api_path.add_sg_distribution_new(self.option("main_id"), self.cnv_diff.work_dir + '/win.stat.xls', types)  # 染色体分布图导表
        api_path.add_sg_snp_indel_compare_detail(self.option("main_id"), self.cnv_diff.output_dir + '/pos.variant',
                                                 types)
        api_path.add_sg_snp_indel_compare_eff_stat(self.option("task_id"), self.option("main_id"),
                                                   self.cnv_diff.output_dir + '/Eff.stat', types)
        api_path.add_sg_snp_indel_compare_stat(self.option("task_id"), self.option("main_id"),
                                               self.cnv_diff.output_dir + '/Ann.stat', types)
        download_path = self._sheet.output.rstrip('/') + "/sample_compare/pos.variant"  # 主表增加详情表上传磁盘的路径，用于下载表格，modified by zengjing 20180612
        api_path.update_download_path(self.option("main_id"), types, download_path)
        self.end()

    def run(self):
        self.run_sample_compare()
        super(SampleCompareWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(SampleCompareWorkflow, self).end()
