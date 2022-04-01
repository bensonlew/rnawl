# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 20200416

import os
import re
import math
import time
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from biocluster.file import getsize, exists
from biocluster.file import download
from biocluster.config import Config
from bson.objectid import ObjectId



class PrimerDesignWorkflow(Workflow):
    """
    相关性小工具工作流
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PrimerDesignWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "seq", "type": "infile", "format": "tool_lab.seq_table"},
            {"name": "primer_length", "type": "int", "default": 20},
            {"name": "tm1", "type": "float", "default": 55},  # float, Tm1 (℃)
            {"name": "tm2", "type": "float", "default": 60},  # float, Tm2 (℃),要大于tm1
            {"name": "gc1", "type": "float", "default": 40},  # float, gc1 (℃)
            {"name": "gc2", "type": "float", "default": 60},  # float, gc2 (℃),要大于gc1
            {"name": "product_size", "type": "string", "default": "100-300"},  # Product Size(bp),数值要为int,范围间用-分隔，多个用,分隔,如：100-300
            {"name": "primer_num", "type": "int", "default": 3},  # Max Pairs Primer Number,范围:[1,5]
            {'name': 'main_id', 'type': 'string'},
            {'name': "update_info", 'type': 'string'},
            {"name": 'source', 'type': 'string', 'default': 'tool_lab'}, #['tool_lab', 'project'],
            {"name": 'relate_id', 'type': 'string'},
            {'name': 'task_id', 'type': 'string'},
            {'name': 'project_task_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.primer = self.add_tool("tool_lab.primer_design")
        self.seq = ''

    def check_options(self):
        if not self.option("primer_length"):
            raise OptionError("请输入引物长度")
        if self.option("tm1") >= self.option("tm2"):
            raise OptionError("tm1要小于tm2")
        if self.option("gc1") >= self.option("gc2"):
            raise OptionError("gc1要小于gc2")
        if not self.option("product_size"):
            raise OptionError("请设置product_size")
        if not re.search(r"(\d+)-(\d+).*", self.option("product_size")):
            raise OptionError("%s product_size格式不正确,用-分隔", variables=(self.option("product_size")))
        if not self.option("primer_num"):
            raise OptionError("请设置primer_num", code="34504406")
        if self.option("primer_num") > 5 or self.option("primer_num") < 1:
            raise OptionError("primer_num:%s范围不为[1,5]", variables=(self.option("primer_num")))

    def run_primer(self):
        if self.option('source') == 'tool_lab':
            seq_file = self.option('seq').path
        elif self.option('source') == 'project':
            seq_file = self.seq_path
        with open(seq_file, 'r') as s:
            lines = s.readlines()
            if len(lines) == 1:
                for line in lines:
                    self.seq = str(line.rstrip().split('\t')[0])
            else:
                self.seq = str(lines[1].rstrip().split('\t')[0])
        self.logger.info(self.seq)
        options = {
            "seq": self.seq,
            "primer_length": self.option('primer_length'),
            "tm1": self.option('tm1'),
            'tm2': self.option('tm2'),
            "gc1": self.option('gc1'),
            "gc2": self.option('gc2'),
            "product_size": self.option('product_size'),
            "primer_num": self.option('primer_num'),
        }
        self.primer.set_options(options)
        self.primer.on("end", self.set_output, "primer_design")
        self.primer.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'primer_design':
            self.linkdir(obj.output_dir, 'primer_design')
        primer_input = os.path.join(self.work_dir, 'seq.txt')
        if self.option('source') == 'project' and os.path.exists(primer_input):
            os.link(primer_input, os.path.join(self.output_dir, 'primer_design', 'primer_input.txt'))
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

    # def set_db(self):
    #     self.logger.info("保存结果到mongo")
    #     api_primer = self.api.api("tool_lab.primer_design")
    #     primer_result = os.path.join(self.output_dir, "primer_design/variation.result")
    #     api_primer.add_sg_primer_detail(self.option('main_id'), primer_result)
    #     self.end()

    def run(self):
        if self.option('source') == 'project':
            self.file_path = self.check_file_path()
            self.seq_path = self.download_s3_file(self.file_path, 'seq.txt')
        self.run_primer()
        super(PrimerDesignWorkflow, self).run()

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
        elif exists(path):
            download(path, to_path)
        else:
            self.set_error('file can not find %s', variables=(path), code = '13700502')
        return to_path

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ['primer_results.txt', 'txt', '引物设计结果文件', 0],
            ['primer_input.txt', 'txt', '引物设计输入文件', 0],
        ])
        super(PrimerDesignWorkflow, self).end()
