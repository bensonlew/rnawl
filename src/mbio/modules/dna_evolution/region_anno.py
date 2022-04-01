# -*- coding: utf-8 -*-
# __author__ = 'liuwentian'
# modified 2018.1123


from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from biocluster.config import Config
import os, re
import json
import time, datetime
from bson.objectid import ObjectId


class RegionAnnoModule(Module):
    """
    群体进化，GWAS关联分析接口的module
    """
    def __init__(self, work_id):
        super(RegionAnnoModule, self).__init__(work_id)
        options = [
            {"name": "main_id", "type": "string"},  #
            {"name": "task_id", "type": "string"},  #
            {"name": "pop_summary_path", "type": "infile", "format": "dna_evolution.pop_summary"},
            {"name": "region_path", "type": "string"},  # sg_anno_params的region_path
            {"name": "update_info", "type": "string"},
            {"name": "region_select", "type": "string"},  # chr1/1/10000
            {"name": "genome_version_id", "type": "string"}
        ]
        self.add_option(options)
        self.region_anno = self.add_tool("dna_evolution.region_anno")
        self.go_summary = self.add_tool("dna_evolution.go_summary")
        self.pop_summary_path = 0

    def check_options(self):
        if not self.option("main_id"):
            raise OptionError("请设置main_id")
        if not self.option("pop_summary_path").prop['path']:
            raise OptionError("请上传pop_summary_path")
        if not self.option("task_id"):
            raise OptionError("请设置task_id")
        if not self.option("region_path"):
            raise OptionError("请设置region_path")
        if not self.option("genome_version_id"):
            raise OptionError("请设置genome_version_id")
        if not self.option("region_select"):
            raise OptionError("请设置region_select")

    def run_region_anno(self):
        """
        """
        go_summary_path = os.path.join(self.go_summary.output_dir, "pop.2.enrich")
        self.region_anno.set_options({
            "main_id": self.option("main_id"),
            "task_id": self.option("task_id"),
            "pop_summary_path": self.pop_summary_path,
            "region_path": self.option("region_path"),
            "update_info": self.option("update_info"),
            "genome_version_id": self.option("genome_version_id"),
            "go_summary_path": go_summary_path,
            "pathway_path": self.parent._sheet.output
        })
        self.region_anno.on('end', self.set_output, "region_anno")
        self.region_anno.on('end', self.end)
        self.region_anno.run()
        print("\n##### region_anno运行结束\n")

    def run_go_summary(self):
        """
        """
        options = {
            "pop_summary": self.pop_summary_path
        }
        self.go_summary.set_options(options)
        self.go_summary.on('end', self.set_output, "go_summary")
        self.go_summary.on("end", self.run_region_anno)
        self.go_summary.run()

    def new_pop_summary(self):
        """
        根据区域过滤pop.summary
        """
        if json.loads(self.option("region_select"))[0] == "All Region":
            self.pop_summary_path = self.option("pop_summary_path").prop["path"]
        else:
            new_pop_summary_path = os.path.join(self.work_dir, "pop.summary")
            n = 0
            with open(self.option("pop_summary_path").prop['path'], 'r')as fr:
                lines = fr.readlines()
                with open(new_pop_summary_path, 'w')as fw:
                    fw.write(lines[0])
                with open(new_pop_summary_path, 'a')as fa:
                    for line in lines[1:]:
                        tmp = line.strip().split('\t')
                        for i in json.loads(self.option("region_select")):
                            temp = i.strip("\"").split("-")
                            chr = temp[0]
                            start = temp[1]
                            end = temp[2]
                            if tmp[4] == chr and int(tmp[5]) >= int(start) and int(tmp[6]) <= int(end):
                                fa.write(line)
                                n = n + 1
                            continue
            self.pop_summary_path = new_pop_summary_path
            if n == 0:
                self.set_error("过滤后的pop.summary文件为空！")

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'region_anno':
            self.linkdir(obj.output_dir, 'region_anno')
        elif event['data'] == 'go_summary':
            self.linkdir(obj.output_dir, 'go_summary')
        else:
            pass

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
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
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def run(self):
        super(RegionAnnoModule, self).run()
        self.new_pop_summary()
        self.run_go_summary()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        super(RegionAnnoModule, self).end()
