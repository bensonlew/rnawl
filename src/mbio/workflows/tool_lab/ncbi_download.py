# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_file,link_dir
import os

class NcbiDownloadWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(NcbiDownloadWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "accsession_file", "type": "infile", 'format': "bacgenome.simple_file"},
            {"name": "accsession", "type": "string"},
            {"name": "method", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.download = self.add_tool('toolapps.ncbi_download')

    def check_options(self):
        if not self.option("accsession_file").is_set and not self.option("accsession"):
            raise OptionError("必须传入登录号！", code="")

    def run(self):
        self.run_download()
        super(NcbiDownloadWorkflow, self).run()

    def run_download(self):
        if self.option("accsession"):
            all_sample = []
            for i in self.option("accsession").strip().split("\n"):
                if i:
                    all_sample.append(i)
            self.list = ";".join(all_sample)
            self.download.set_options({
                "sample_list": self.list,
            })
        else:
            self.list = ''
            with open(self.option("accsession_file").prop['path']) as v:
                data = v.readlines()
                for i in data:
                    if i.strip():
                        self.list += (i.strip()+";")
            self.list = self.list.strip(";")
            self.download.set_options({
                "sample_list": self.list,
            })
        self.download.on("end", self.set_output)
        self.download.run()

    def set_output(self):
        self.logger.info("set_output")
        gca_gcf_dict = {}
        with open(self.download.work_dir + "/ncbi/gca_gcf.txt") as f:
            data = f.readlines()
            if data:
                for lines in data:
                    line = lines.strip().split("\t")
                    gca_gcf_dict[line[0]] = line[1]
        sample_list1 = []
        for i in self.list.split(";"):
            sample_list1.append(i)
        with open(self.output_dir + "/README.xls","w") as t:
            t.write("genomeID" + "\t" + "stat" + "\n")
            for x in sample_list1:
                if os.path.exists(self.download.work_dir + "/ncbi"):
                    for v in os.listdir(self.download.work_dir + "/ncbi"):
                        if x in v:
                            has = "true1"
                            link_dir(self.download.work_dir + "/ncbi/" + v, self.output_dir + "/" + v)
                            break
                        else:
                            if x in gca_gcf_dict.keys() and gca_gcf_dict[x] in v:
                                has = "true2"
                                link_dir(self.download.work_dir + "/ncbi/" + v, self.output_dir + "/" + v)
                                break
                            else:
                                has = "false"
                    if has == "true1":
                        t.write(x + "\t" + "下载完成" + "\n")
                    elif has == "true2":
                        t.write(x + "\t" + "下载完成(对应GCF号为{})".format(gca_gcf_dict[x]) + "\n")
                    else:
                        t.write(x + "\t" + "下载失败" + "\n")
                else:
                    t.write(x + "\t" + "下载失败" + "\n")
        self.set_db()

    def set_db(self):
        self.logger.info("开始导表")
        api_main = self.api.api("tool_lab.common")
        api_main.add_main_table("ncbi_download", main_id = self.option('main_id'))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "NCBI参考基因组批量下载结果目录", 0],
            ["./*", "", "下载结果", 0],
            ["./README.xls", "xls", "下载情况统计表", 0],
        ])
        super(NcbiDownloadWorkflow, self).end()