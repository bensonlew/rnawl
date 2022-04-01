# -*- coding: utf-8 -*-
# __author__ = 'wentian.liu'
# modified 20190319

import os
import json
import gevent
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class PrimerDesignWorkflow(Workflow):
    """
    引物设计
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PrimerDesignWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "ref_fa", "type": "string"},
            {"name": "marker_type", "type": "string"},  # 标记类型snp/indel,ssr
            {"name": "data_type", "type": "string"},  # 选择标记位置信息 vcf,location(位置文件),custom(自定义)
            {"name": "marker_detail", "type": "string"},
            # 选择标记位置信息 location:[{location:chr:start-end},{}] custom:[chr:start-end,....]
            {"name": "condition", "type": "string"},  # 设置引物设计条件 {snp:{tm1:,tm2:....},indel"{}}
            {"name": "vcf", "type": "infile", "format": "wgs_v2.vcf"},  # vcf文件:原始vcf或者运行结果vcf。
            {"name": "ssr", "type": "infile", "format": "wgs_v2.bam"},  # 参考基因组中的ssr.ref.result.xls。
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.cut_vcf = self.add_tool("wgs_v2.cut_vcf")
        # self.primer_design = self.add_tool("wgs.primer_design")
        self.primer_design_tools = []
        self.primer_num = 0

    def check_options(self):
        if not self.option('ref_fa'):
            raise OptionError('必须输入ref_fa', code="14500301")
        if not self.option('marker_type'):
            raise OptionError('必须输入marker_type', code="14500301")
        if not self.option('data_type'):
            raise OptionError('必须输入data_type', code="14500302")
        # if not self.option('marker_detail'):
        #     raise OptionError('必须输入marker_detail', code="14500302")
        if not self.option('condition'):
            raise OptionError('必须输入condition', code="14500302")
        if not self.option('vcf'):
            raise OptionError('必须输入vcf', code="14500302")
        # if not self.option('ssr'):
        #     raise OptionError('必须输入ssr', code="14500303")
        if not self.option('main_id'):
            raise OptionError('必须输入main_id', code="14500302")
        return True

    def run_cut_vcf(self):
        """
        生成primer_design所需的输入文件。
        :return:
        """
        options = {
            "marker_type": self.option('marker_type'),
            "data_type": self.option('data_type'),
            "vcf": self.option('vcf').prop['path']
        }
        if not self.option('data_type') == "vcf":
            options["marker_detail"] = self.option('marker_detail')
        if self.option('marker_type') == "ssr" and self.option('data_type') == "custom":
            options["ssr"] = self.option('ssr').prop['path']
        self.cut_vcf.set_options(options)
        self.cut_vcf.on("end", self.run_primer_design)
        self.cut_vcf.on("end", self.set_output, "cut_vcf")
        self.cut_vcf.run()

    def run_primer_design(self):
        if self.option('marker_type') == "snp/indel":
            type_list = ["snp", "indel"]
            for i in type_list:
                condition = json.loads(self.option('condition'))[i]
                if int(condition["primer_num"]) > self.primer_num:
                    self.primer_num = int(condition["primer_num"])
                self.primer_design = self.add_tool("wgs_v2.primer_design")
                options = {
                    "ref_fa": self.option('ref_fa'),
                    "diff_variant": os.path.join(self.output_dir, "cut_vcf/" + i + "_infile.txt"),
                    "tm1": condition["tm1"],
                    "tm2": condition["tm2"],
                    "product_size": condition["product_size"],
                    "primer_num": condition["primer_num"],
                    "project": "wgs_v2"
                }
                self.primer_design.set_options(options)
                self.primer_design_tools.append(self.primer_design)
            self.on_rely(self.primer_design_tools, self.conbine)
            for tool in self.primer_design_tools:
                gevent.sleep(1)
                tool.run()
        else:
            self.primer_design = self.add_tool("wgs_v2.primer_design")
            condition = json.loads(self.option('condition'))["ssr"]
            if int(condition["primer_num"]) > self.primer_num:
                self.primer_num = int(condition["primer_num"])
            options = {
                "ref_fa": self.option('ref_fa'),
                "diff_variant": os.path.join(self.output_dir, "cut_vcf/infile.txt"),
                "tm1": condition["tm1"],
                "tm2": condition["tm2"],
                "product_size": condition["product_size"],
                "primer_num": condition["primer_num"],
                "project": "wgs_v2"
            }
            self.primer_design.set_options(options)
            self.primer_design.on("end", self.set_output, "primer_design")
            self.primer_design.run()

    def conbine(self):
        """
        当marker_type为snp/indel时，将两个结果文件合并。
        """
        write_lines = ""
        line_list = []
        for i in range(0, 2):
            with open(os.path.join(self.primer_design_tools[i].output_dir, "variation.result"), "r")as fr:
                lines = fr.readlines()
                line_list += lines[i:]
        for line in line_list:
            write_lines += line
        newdir = os.path.join(self.output_dir, "primer_design")
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        with open(os.path.join(self.output_dir, "primer_design/variation.result"), "w")as fw:
            fw.write(write_lines)
        self.set_db()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'cut_vcf':
            self.linkdir(obj.output_dir, 'cut_vcf')
        if event['data'] == 'primer_design':
            self.linkdir(obj.output_dir, 'primer_design')
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
        primer_design_api = self.api.api('wgs_v2.primer_design')
        self.logger.info("开始进行primer_design导表")
        primer_id = self.option('main_id')
        primer_result = os.path.join(self.output_dir, "primer_design/variation.result")
        primer_design_api.add_sg_primer_detail(primer_id, primer_result, self.primer_num)
        self.logger.info("primer_design导表成功！")
        self.end()

    def run(self):
        self.run_cut_vcf()
        super(PrimerDesignWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(PrimerDesignWorkflow, self).end()
