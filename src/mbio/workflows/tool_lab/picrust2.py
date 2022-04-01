# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir
import os
import subprocess

class Picrust2Workflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(Picrust2Workflow, self).__init__(wsheet_object)
        options = [
            {"name": "data_table", "type": "infile", 'format': "meta.otu.otu_table"},
            {"name": "in_fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "style", "type": "string", "default": "16S"},  ##选择方法16S、18S、ITS
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.predict = self.add_tool("tool_lab.picrust2_predict")

    def check_options(self):
        if not self.option("data_table").is_set:
            raise OptionError("请传入数据表！", code="")
        if not self.option("in_fasta").is_set:
            raise OptionError("请传入数据表！", code="")

    def run(self):
        self.get_fasta()
        self.run_predict()
        super(Picrust2Workflow, self).run()

    def get_fasta(self):
        raw_fasta = self.option("in_fasta").prop['path']
        self.new_fasta = self.work_dir + '/new_fasta.fasta'
        with open(raw_fasta) as f, open(self.new_fasta,"w") as t:
            data = f.readlines()
            for i in data:
                if i.startswith(">"):
                    t.write(i.strip().split(" ")[0].split("\t")[0] + "\n")
                else:
                    t.write(i)

    def run_predict(self):
        """
        用picrust2进行预测计算
        :return:
        """
        opts = {
            "otu_fasta": self.new_fasta,
            "otu_table": self.option("data_table"),
            "analysis_type": self.option("style")
        }
        self.predict.set_options(opts)
        self.predict.on('end', self.set_db)
        self.predict.run()

    def zip_file(self):
        os.chdir(self.predict.output_dir)
        cmd = "zip -q -r {} ./*".format(self.output_dir + "/picrust2_results.zip")
        try:
            self.logger.info("开始对文件%s进行压缩！" % self.output_dir)
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError:
            self.set_error("对文件%s压缩失败！", variables=(self.output_dir), code="")

    def set_db(self):
        #link_dir(self.predict.output_dir, self.output_dir)
        self.zip_file()
        api_main = self.api.api("tool_lab.common")
        api_main.add_main_table("picrust2", main_id = self.option('main_id'))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "dir", "PICRUSt2功能预测结果目录", 0, "110273"],
            ["./picrust2_results.zip", "zip", "所有结果的压缩文件", 0],
        ])
        super(Picrust2Workflow, self).end()

    def check_data(self):
        data_file = self.option("data_table").prop['path']
        fasta_file = self.option("in_fasta").prop['path']
        data = []
        fasta = []
        num = 0
        self.data_info = {}
        with open(data_file) as f,open(fasta_file) as t:
            data_tmp = f.readlines()
            fasta_tmp = t.read()
            for i in data_tmp[1:]:
                data.append(i.split("\t")[0])
            for i in fasta_tmp.split(">")[1:]:
                fasta.append(i.split("\n"))
            for i in data:
                if i in fasta:
                    num += 1
            self.data_info["data"] = len(data)
            self.data_info["fasta"] = len(fasta)
            self.data_info["fasta"] = len(fasta)