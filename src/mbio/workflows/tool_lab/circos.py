# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 20200508

import os
import gevent
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class CircosWorkflow(Workflow):
    """
    相关性小工具工作流
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CircosWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "chr_file", "type": "string"},
            # {"name": "circos_file", "type": "string"},
            # {"name": "circos_type", "type": "string"},
            {"name": "circos_file_1", "type": "string"},
            {"name": "circos_file_2", "type": "string"},
            {"name": "circos_file_3", "type": "string"},
            {"name": "circos_file_4", "type": "string"},
            {"name": "circos_file_5", "type": "string"},
            {"name": "circos_file_6", "type": "string"},
            {"name": "circos_file_7", "type": "string"},
            {"name": "circos_file_8", "type": "string"},
            {"name": "circos_file_9", "type": "string"},
            {"name": "circos_file_10", "type": "string"},
            {"name": "circos_file_11", "type": "string"},
            {"name": "circos_file_12", "type": "string"},
            {"name": "circos_file_13", "type": "string"},
            {"name": "circos_file_14", "type": "string"},
            {"name": "circos_file_15", "type": "string"},
            {"name": "circos_type_1", "type": "string"},
            {"name": "circos_type_2", "type": "string"},
            {"name": "circos_type_3", "type": "string"},
            {"name": "circos_type_4", "type": "string"},
            {"name": "circos_type_5", "type": "string"},
            {"name": "circos_type_6", "type": "string"},
            {"name": "circos_type_7", "type": "string"},
            {"name": "circos_type_8", "type": "string"},
            {"name": "circos_type_9", "type": "string"},
            {"name": "circos_type_10", "type": "string"},
            {"name": "circos_type_11", "type": "string"},
            {"name": "circos_type_12", "type": "string"},
            {"name": "circos_type_13", "type": "string"},
            {"name": "circos_type_14", "type": "string"},
            {"name": "circos_type_15", "type": "string"},
            {'name': 'outerRadius', 'type': 'int'},
            {'name': 'thickness', 'type': 'int'},
            {'name': 'circle_num', 'type': 'int'},
            {'name': 'main_id', 'type': 'string'},
            {'name': "update_info", 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.file = {}
        self.type = {}
        self.circos_tools = list()

    def check_options(self):
        if not self.option("chr_file"):
            raise OptionError("请设置坐标文件")
        # if not self.option("circos_file"):
        #     raise OptionError("请设置每一圈画图数据")
        if not self.option("circos_file_1"):
            raise OptionError("请设置画图数据")
        # if not self.option("circos_type"):
        #     raise OptionError("请设置每一圈画图类型")
        if not self.option("circos_type_1"):
            raise OptionError("请设置画图类型")
        # if not len(self.option("circos_type").strip().split(",")) == len(self.option("circos_file").strip().split(",")):
        #     raise OptionError("circos圈文件与画图类型数量不一致，请指定每一圈画图类型")

    def run_circos(self):
        self.circos_file = []
        self.circos_type = []
        with open(self.work_dir + "/circos.list","w") as bl:
            for il in range(1,self.option("circle_num")+1):
                bl.write("{}\t{}\n".format(self.option("circos_file_{}".format(il)),self.option("circos_type_{}".format(il))))
        file_list = os.path.join(self.work_dir,"circos.list")
        # files = self.option("circos_file").strip().split(",")
        # types = self.option("circos_type").strip().split(",")
        # self.circos_tools = []
        with open (file_list,"r") as lf:
            while 1:
                line = lf.readline()
                if not line:
                    break
                fd = line.rstrip().split("\t")
                self.file = fd[0]
                self.type = fd[1]
                self.circos_file.append(fd[0])
                self.circos_type.append(fd[1])
                # self.logger.info(self.circos_file)
                # self.logger.info(self.circos_type)
        files = self.circos_file
        print files
        types = self.circos_type
        print types
        files.append(self.option("chr_file"))
        types.append("c_chr")
        for i in range(len(files)):
            circos = self.add_tool("tool_lab.circos")
            circos.set_options({
                "circos_file": files[i].strip(),
                "circos_type": types[i].strip(),
                "index": str(i),
            })
            self.circos_tools.append(circos)
        for j in range(len(self.circos_tools)):
            self.circos_tools[j].on("end", self.set_output, "circos")
        if self.circos_tools:
            if len(self.circos_tools) > 1:
                self.on_rely(self.circos_tools, self.wait_end)
            elif len(self.circos_tools) == 1:
                self.circos_tools[0].on('end', self.wait_end)
        else:
            self.set_error("variant_compare_tools为空！")
        for tool in self.circos_tools:
            gevent.sleep(1)
            tool.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'circos':
            self.linkdir(obj.output_dir, 'circos')

    def wait_end(self):
        gevent.sleep(1)
        self.end()

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
        files = os.listdir(os.path.join(self.output_dir, "circos"))
        file_list = []
        for file in files:
            file_list.append(os.path.join(self._sheet.output, file))
        self.logger.info(file_list)
        new_path = ",".join(file_list)
        self.logger.info(new_path)
        self.logger.info("保存结果到mongo")
        params_dict = {}
        params_dict["outerRadius"] = self.option("outerRadius")
        params_dict["innerRadius"] = abs(self.option("outerRadius") - self.option("thickness"))
        params_dict["circle_num"] = self.option("circle_num")
        params_dict["path"] = new_path
        params_dict["circos_type"] = self.option("circos_type")
        api_manhattan = self.api.api("tool_lab.circos")
        self.logger.info(params_dict)
        api_manhattan.add_sg_circos_path(self.option("main_id"), **params_dict)

    def run(self):
        self.run_circos()
        super(CircosWorkflow, self).run()

    def end(self):
        self.set_db()
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(CircosWorkflow, self).end()