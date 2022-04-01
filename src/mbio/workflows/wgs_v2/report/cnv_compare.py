# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20190307

import re
import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.file import getsize, exists, list_dir


class CnvCompareWorkflow(Workflow):
    """
    交互分析：cnv的差异比较分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CnvCompareWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "infile_path", "type": "string"},
            {"name": "region", "type": "string"},  # chr1:1-500,chr2:2-4
            {"name": "marktype", "type": "string"},  # same or diff or all same,diff
            {"name": "samples", "type": "string"},  # 样本对a|b,b|c
            {"name": "region_type", "type": "string", 'default': 'real_region'},  # 用于后面区域过滤，是过滤
            # 范围还是具体的点
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "project_type", "type": "string"},
            {"name": "analysis_model", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.cnv_diff = self.add_tool("wgs_v2.cnv_compare")

    def check_options(self):
        if not self.option("samples"):
            raise OptionError("缺少samples参数")
        if not self.option("infile_path"):
            raise OptionError("缺少infile_path参数")
        if not self.option("region"):
            raise OptionError("缺少region参数")
        if not self.option("marktype"):
            raise OptionError("缺少marktype参数")
        return True

    def run_cnv_diff(self):
        options = {
            "infile_path": self.work_dir + '/temp/',
            "region": self.option("region"),
            "marktype": self.option("marktype"),
            "samples": self.option("samples"),
            "region_type": self.option("region_type"),
            "analysis_model": self.option("analysis_model")
        }
        self.cnv_diff.set_options(options)
        self.cnv_diff.on("end", self.set_output, "cnv_diff")
        self.cnv_diff.run()

    def dowmload_from_s3(self):
        """
        从对象存储中下载文件到指定路径
        :return:
        """
        if not os.path.exists(self.work_dir + "/temp"):
            os.mkdir(self.work_dir + "/temp")
        self.logger.info("开始下载对象存储中的文件！")
        samples = []
        for m in self.option("samples").split(','):
            for n in m.split('|'):
                samples.append(n)
        transfer = MultiFileTransfer()
        for sample in list(set(samples)):
            source = os.path.join(self.option("infile_path"), "{}.cnv.anno.xls".format(sample))
            if not exists(source):
                self.set_error("文件%s不存在！" % source)
            transfer.add_download(source, '{}/temp/'.format(self.work_dir))
        transfer.perform()
        self.logger.info("下载对象存储中的文件成功！")

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'cnv_diff':
            self.linkdir(obj.output_dir, '')
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
        cnv_diff = self.api.api("wgs_v2.cnv_compare")
        if self.option("project_type"):
            cnv_diff._project_type = self.option("project_type")
        subname = []
        vcf_path = []
        for m in os.listdir(self.output_dir):
            if re.match('.*\.stat\.xls$', m):
                nname = m.split('.')[0]
                subname.append("{}".format(nname))
                vcf_path.append(self._sheet.output.rstrip('/') + '/{}.xls'.format(nname))
                cnv_diff.add_sg_cnv_compare_detail(self.option('main_id'),
                                                   os.path.join(self.output_dir, '{}.xls'.format(nname)),
                                                   "{}".format(nname), self.option("analysis_model"))
                cnv_diff.add_sg_cnv_compare_stat(self.option('main_id'), os.path.join(self.output_dir, m),
                                                 "{}".format(nname))
        cnv_diff.update_subname(self.option('main_id'), subname, vcf_path)
        self.end()

    def run(self):
        self.dowmload_from_s3()
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
