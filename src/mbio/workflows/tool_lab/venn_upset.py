# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 20200416
from bson.objectid import ObjectId
import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import time
from biocluster.file import getsize, exists
from biocluster.file import download


class VennUpsetWorkflow(Workflow):
    """
    相关性小工具工作流
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(VennUpsetWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "upset_table", "type": "infile", "format": "tool_lab.upset_table"},
            {'name': 'task_id', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
            {'name': "update_info", 'type': 'string'},
            {'name': 'project_task_id', 'type': 'string'},
            {"name": 'source', 'type': 'string', 'default': 'tool_lab'},
            {'name': 'relate_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.venn_upset = self.add_tool("tool_lab.venn_upset")

    def check_options(self):
        if not self.option("upset_table"):
            raise OptionError("必须输入upset_table文件")
        return True

    def run_venn_upset(self):
        if self.option('source') == 'project':
            upset = self.upset_venn_path
        if self.option('source') == 'tool_lab':
            upset = self.option('upset_table').path
        options = {
            "upset_table": upset,
        }
        self.venn_upset.set_options(options)
        self.venn_upset.on("end", self.set_output, "venn_upset")
        self.venn_upset.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'venn_upset':
            self.linkdir(obj.output_dir, 'venn_upset')
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
        api_venn_upset = self.api.api("tool_lab.venn_upset")
        api_venn_upset.add_sg_venn_upset(self.option("main_id"),
                                            self.output_dir + "/venn_upset/venn_upset.txt")
        self.end()

    def run(self):
        if self.option('source') == 'project':
            self.file_path = self.check_file_path()
            self.upset_venn_path = self.download_s3_file(self.file_path, 'upset_venn.txt')
        self.run_venn_upset()
        super(VennUpsetWorkflow, self).run()

    def check_file_path(self):
        collection_name = 'tool_thurl'
        project_type = 'tool_lab'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        conn_upset = db[collection_name]
        status = 'start'
        count_time = 0
        while status == 'start':
            if count_time > 600:
                self.set_error('超过十分钟还没有结果文件生成，请检查是否生成文件时报错')
                break
            time.sleep(10)
            print 'sleep 10s'
            try:
                upset = conn_upset.find_one(
                    {'task_id': self.option('project_task_id'), 'relate_id': ObjectId(self.option('relate_id'))})
                status = upset['status']
            except:
                pass
            count_time += 10
        upset = conn_upset.find_one(
            {'task_id': self.option('project_task_id'), 'relate_id': ObjectId(self.option('relate_id'))})
        file_path = upset['file_path']
        return file_path

    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        elif os.path.exists(path):
            to_path = path
        elif exists(path):
            download(path, to_path)
        else:
            self.set_error('file can not find %s', variables=(path), code = '13700502')
        return to_path

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(VennUpsetWorkflow, self).end()
