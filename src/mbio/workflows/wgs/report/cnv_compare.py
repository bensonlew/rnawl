# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180409

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class CnvCompareWorkflow(Workflow):
    """
    交互分析：cnv的差异比较分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CnvCompareWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "sample1", "type": "infile", "format": "bsa.vcf"},
            {"name": "sample2", "type": "infile", "format": "bsa.vcf"},
            {"name": "is_same", "type": "string"},
            {"name": "variation_type", "type": "string", "default": ""},
            {"name": "variation_len", "type": "string", "default": ""},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "project_type", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.cnv_diff = self.add_tool("wgs.cnv_compare")

    def check_options(self):
        if not self.option('sample1'):
            raise OptionError('必须输入样本1', code="14500301")
        if not self.option('sample2'):
            raise OptionError('必须输入样本2', code="14500302")
        if not self.option('is_same'):
            raise OptionError('必须输入相同或者不同的参数', code="14500303")
        if not self.option('variation_type'):
            raise OptionError('必须输入variation_type', code="14500304")
        if not self.option('variation_len'):
            raise OptionError('必须输入variation_len', code="14500305")
        return True

    def run_cnv_diff(self):
        options = {
            "sample1": self.option("sample1").prop['path'],
            "sample2": self.option("sample2").prop['path'],
            "is_same": self.option("is_same"),
            "variation_type": self.option("variation_type"),
            "variation_len": self.option("variation_len")
        }
        self.cnv_diff.set_options(options)
        self.cnv_diff.on("end", self.set_output, "cnv_diff")
        self.cnv_diff.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'cnv_diff':
            self.linkdir(obj.output_dir, 'cnv_diff')
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
        cnv_diff = self.api.api("wgs.cnv_compare")
        if self.option("project_type"):
            cnv_diff._project_type = self.option("project_type")
        n = 0
        for m in os.listdir(self.cnv_diff.output_dir):
            cnv_diff.add_cnv_diff_result(self.option("main_id"), os.path.join(self.cnv_diff.output_dir, m))
            n += 1
        if n > 1:
            self.set_error("output_dir路径中xls文件只能有一个，请检查！", code="14500303")
        self.end()

    def run(self):
        self.run_cnv_diff()
        super(CnvCompareWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(CnvCompareWorkflow, self).end()
