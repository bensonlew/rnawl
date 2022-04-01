# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180409

import re
import os
from bson.objectid import ObjectId
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class DrawCircosWorkflow(Workflow):
    """
    交互分析：circos图
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DrawCircosWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "windows", "type": "int", "default": 100000},
            {"name": "snp", "type": "infile", 'format': 'bsa.vcf'},  # snp.anno.primary.vcf
            {"name": "indel", "type": "infile", 'format': 'bsa.vcf'},  # indel.anno.primary.vcf
            # {"name": "gff", "type": "string"},   # 这个是gene.gff文件，不是ref.gff文件  都不进行传了直接通过基因组获取
            {"name": "sv_path", "type": "infile", 'format': 'bsa.dir'},
            {"name": "cnv_path", "type": "infile", 'format': 'bsa.dir'},
            # {"name": "chrlist", "type": "string"},
            {"name": "chrs", "type": "string", "default": "all"},  # 传进来是染色体列，如：chr1,chr2,chr3 逗号分隔
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "genome_version_id", "type": "string"},
            {"name": "project_type", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.circos = self.add_module("wgs.draw_circos")

    def check_options(self):
        if not self.option('snp'):
            raise OptionError('必须提供snp结果表', code="14500501")
        if not self.option('indel'):
            raise OptionError('必须提供indel结果表', code="14500502")
        if not self.option("sv_path"):
            raise OptionError("必须提供sv_path文件！", code="14500503")
        if not self.option("cnv_path"):
            raise OptionError("必须提供cnv_path文件！", code="14500504")
        return True

    def run_circos(self):
        options = {
            "windows": self.option("windows"),
            "snp": self.option("snp").prop['path'],
            "indel": self.option("indel").prop['path'],
            "genome_version_id": self.option("genome_version_id"),
            # "gff": self.option("gff"),
            "sv_path": self.option("sv_path").prop['path'],
            "cnv_path": self.option("cnv_path").prop['path'],
            # "chrlist": self.option("chrlist"),
            "chrs": self.option("chrs")
        }
        self.circos.set_options(options)
        self.circos.on("end", self.set_output, "circos")
        self.circos.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'circos':
            self.linkdir(obj.output_dir, self.output_dir)
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
        self.logger.info("设置circos的结果路径")
        sample_id = []
        for file_ in os.listdir(self.output_dir):
            m = re.match(r"(.*)\.png$", file_)
            if m:
                sample_id.append(m.group(1))
        api = self.api.api("wgs.api_base")
        if self.option("project_type"):
            api._project_type = self.option("project_type")
        target_dir = self._sheet.output.rstrip('/') + '/'
        api.update_db_record("sg_circos", {"_id": ObjectId(self.option("main_id"))},
                             {"circos_result_path": target_dir, "sample_list": sample_id})
        self.end()

    def run(self):
        self.run_circos()
        super(DrawCircosWorkflow, self).run()

    def end(self):
        """
        这里后面要重新定义下文件名字
        :return:
        """
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(DrawCircosWorkflow, self).end()
