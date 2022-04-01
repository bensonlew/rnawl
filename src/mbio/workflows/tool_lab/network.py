# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 20200623

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class NetworkWorkflow(Workflow):
    """
    相关性小工具工作流
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(NetworkWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "link_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "node_table", "type": "infile", "format": "tool_lab.simple"},
            {'name': 'main_id', 'type': 'string'},
            {'name': "update_info", 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.network = self.add_tool("tool_lab.network")

    def check_options(self):
        if not self.option("link_table").is_set:
            raise OptionError("必须输入link_table文件")
        if not self.option("node_table").is_set:
            raise OptionError("必须输入node_table文件")
        return True

    def run_network(self):
        options = {
            "link_table": self.option("link_table"),
            "node_table": self.option("node_table"),
        }
        self.network.set_options(options)
        self.network.on("end", self.set_output, "network")
        self.network.on("end", self.end)
        self.network.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'network':
            self.linkdir(obj.output_dir, 'network')

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
        self.logger.info("保存network结果到mongo")
        api_network = self.api.api("tool_lab.network")
        api_network.add_sg_network(self.option("main_id"), self.work_dir + "/Network/output/newnodes",
                                   self.work_dir + "/Network/output/newlinks")
        self.logger.info("-------------------------------准备结束------------------------------------------------------")

    def run(self):
        self.run_network()
        # self.set_db()
        super(NetworkWorkflow, self).run()

    def end(self):
        self.set_db()
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(NetworkWorkflow, self).end()
