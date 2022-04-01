# -*- coding: utf-8 -*-
# __author__ = 'wentian.liu'
# modified 20190308

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class CnvCallWorkflow(Workflow):
    """
    cnv_call
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CnvCallWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "ref_gff", "type": "string"},
            {"name": "bam_list", "type": "infile", "format": "bsa.dir"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "target_path", 'type': "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.cnv_call = self.add_module("wgs.cnv_call")
        self.ref_gff = os.path.join(self.config.SOFTWARE_DIR, ("database/dna_geneome/" + self.option("ref_gff")))

    def check_options(self):
        if not self.option('ref_gff'):
            raise OptionError('必须输入ref_gff', code="14500301")
        if not self.option('bam_list'):
            raise OptionError('必须输入bam_list', code="14500302")
        if not self.option('main_id'):
            raise OptionError('必须输入main_id', code="14500303")
        if not self.option('task_id'):
            raise OptionError('必须输入task_id', code="14500304")
        return True

    def run_cnv_call(self):
        options = {
            "ref_gff": self.ref_gff,
            "bam_list": os.path.join(self.work_dir, "bam.list")
        }
        self.cnv_call.set_options(options)
        self.cnv_call.on("end", self.set_output, "cnv_call")
        self.cnv_call.run()

    def make_bam_list(self):
        file_list = os.listdir(self.option('bam_list').prop['path'])
        write_lines = ""
        for i in file_list:
            if i.endswith(".sort.bam"):
                tmp = i.strip().split(".sort.bam")
                sample = tmp[0]
                bam_path = os.path.join(self.option('bam_list').prop['path'], i)
                write_lines += sample + "\t" + bam_path + "\n"
        with open(os.path.join(self.work_dir, "bam.list"), "w")as fw:
            fw.write(write_lines)

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'cnv_call':
            self.linkdir(obj.output_dir, 'cnv_call')
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
        cnv_api = self.api.api('wgs.cnv')
        cnv_api.project_type = "dna_wgs_v2"
        call_id = self.option("main_id")
        cnv_api.add_sg_cnv_call_stat(call_id, os.path.join(self.output_dir, "cnv_call/cnv.stat.xls"))
        self.logger.info("cnv统计表导表成功！")
        self.logger.info("开始进行cnv长度统计导表")
        cnv_api.add_cnv_length_bar(self.option("task_id"), call_id, os.path.join(self.output_dir, "cnv_call/length"))  # 导入length文件夹
        cnv_api.update_db_record("sg_task", {"task_id": self.option("task_id")},
                                 {"cnv_anno_path": (self.option('target_path') + "/cnv_call/anno/")})
        self.logger.info("cnv长度统计导表成功！")
        self.end()

    def run(self):
        self.make_bam_list()
        self.run_cnv_call()
        super(CnvCallWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(CnvCallWorkflow, self).end()
