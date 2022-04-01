# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# last_modify:20180612

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import json


class MlodCalcModule(Module):
    """
    用于拆分pop.filter.marker等文件，随后mlod 计算和汇总
    lasted modified by hongdong@20180628
    1)解决pep8错误
    2)完善infile与outfile
    3)解决当只有一个split文件的时候，合并的时候会出错
    """
    def __init__(self, work_id):
        super(MlodCalcModule, self).__init__(work_id)
        options = [
            {"name": "marker", "type": 'infile', "format": "dna_gmap.marker"},
            {"name": "only", "type": 'int', "default": 200},
            {"name": "total_mlod", "type": "outfile", "format": "dna_gmap.mlod"}
        ]
        self.add_option(options)
        self.split_num = 0
        self.mlod_calc_tools = []

    def check_options(self):
        if not self.option("marker"):
            raise OptionError("请输入pop.filtered.marker或者Total.bin.marker文件", code="24800301")
        return True

    def split_mark_run(self):
        """
        用于将marker进行分群
        修复了当marker行数小于200的时候，文件为空的问题，modified by hongdong 20180714
        """
        if not os.path.exists(self.work_dir + "/sub"):
            os.mkdir(self.work_dir + "/sub")
        else:
            os.system("rm -r {}".format(self.work_dir + "/sub"))
            os.mkdir(self.work_dir + "/sub")
        with open(self.option("marker").prop['path'], "r") as f:
            line = f.readlines()
            if len(line) < self.option("only"):
                with open(self.work_dir + "/sub/sub.0.genotype", "w") as w:
                    for m in line:
                        w.write(m)
            else:
                self.split_num = int(len(line) / self.option("only"))
                i = 0
                while i < self.split_num:
                    name = "sub." + str(i) + ".genotype"
                    f = os.path.join(self.work_dir + "/sub/", name)
                    with open(f, "w") as w:
                        w.write(line[0])
                        j = i * self.option("only")
                        while j < len(line):
                            w.write(line[j])
                            j += 1
                    i += 1

    def mlod_calc_run(self):
        for i in os.listdir(self.work_dir + "/sub"):
            mlod_calc = self.add_tool("dna_gmap.mlod_calc")
            m = i.split(".")[1]
            path = os.path.join(self.work_dir + "/sub/" + i)
            if m == self.split_num - 1:
                mlod_calc.set_options({
                        "sub_genotype": path
                    })
                self.mlod_calc_tools.append(mlod_calc)
            else:
                mlod_calc.set_options({
                    "sub_genotype": path,
                    "only": self.option("only")
                })
                self.mlod_calc_tools.append(mlod_calc)
        for j in range(len(self.mlod_calc_tools)):
            self.mlod_calc_tools[j].on("end", self.set_output, 'mlod_calc_dir')
        if self.mlod_calc_tools:
            if len(self.mlod_calc_tools) > 1:
                self.on_rely(self.mlod_calc_tools, self.mlod_merge_run)
            elif len(self.mlod_calc_tools) == 1:
                self.mlod_calc_tools[0].on('end',  self.mlod_merge_run)
            for tool in self.mlod_calc_tools:
                gevent.sleep(1)
                tool.run()
        else:
            self.set_error("mlod_calc_tools列表为空！", code="24800304")

    def mlod_merge_run(self):
        if len(self.mlod_calc_tools) > 1:
            result_path = self.output_dir + "/mlod_calc_dir"
        else:
            result_path = self.mlod_calc_tools[0].output_dir
        code = os.system("cat {}/sub.*.mlod > {}/Total.mlod".format(result_path, self.output_dir))
        if code == 1:
            self.set_error("merge失败", code="24800301")
        elif code == 0:
            self.logger.info("merge 成功")
        self.option("total_mlod").set_path(self.output_dir + "/Total.mlod")
        self.end()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'mlod_calc_dir':
            self.linkdir(obj.output_dir, 'mlod_calc_dir')
        else:
            pass

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

        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def run(self):
        super(MlodCalcModule, self).run()
        self.split_mark_run()
        self.mlod_calc_run()

    def end(self):
        super(MlodCalcModule, self).end()
