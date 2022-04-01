# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
# last_modify:2017.05.24

from biocluster.module import Module
import os
import shutil
import time
from biocluster.core.exceptions import OptionError


class HmmscanModule(Module):
    def __init__(self, work_id):
        super(HmmscanModule, self).__init__(work_id)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "database", "type": "string", "default": "cazy"},  # 数据库：cazy, pfam  add by zhujuan10180124
            {"name": "lines", "type": "int", "default": 200000},  # 将fasta序列拆分此行数的多个文件
            {"name": "align_result", "type": "outfile", "format": "meta_genomic.hmmscan_table"} , # 设置结果文件后面要用
        ]
        self.add_option(options)  #####检查option是否list格式，其中每个opt是否字典格式
        self.split_fasta = self.add_tool("sequence.split_fasta")
        self.step.add_steps('cazy_align', 'split_fasta')
        self.cazy_align_tools = []
        self.cat_out = self.add_tool("align.cat_hmmscanout")
        self.align_result_path = ''

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="21100301")
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数", code="21100302")
        return True

    def run_split_fasta(self):
        self.split_fasta.set_options({  #######批量设置参数值，并传给tool的option
                                        "fasta": self.option("query"),
                                        "lines": self.option("lines"),
                                        })
        self.split_fasta.on('end', self.run_align)  #########split发生end事件触发run.align
        self.split_fasta.run()

    def run_align(self):
        self.align_result_path = os.path.join(self.output_dir, "align_result")
        if os.path.exists(self.align_result_path):
            pass
        else:
            os.mkdir(self.align_result_path)
        n = 0
        for f in os.listdir(self.split_fasta.output_dir):
            file_path = os.path.join(self.split_fasta.output_dir, f)
            align_hmmscan = self.add_tool('align.hmmscan')
            align_hmmscan.set_options({
                "faa_file": file_path,
                "database": self.option("database")
            })
            if len(os.listdir(
                    self.split_fasta.output_dir)) > 1:  # last_modify by zhujuan 20180124 当 self.cazy_align_tools == 1时， self.set_output 和 self.run_cat_result同时执行，造成结果报错
                align_hmmscan.on('end', self.set_output, 'cazy_align')
            self.cazy_align_tools.append(align_hmmscan)
            n += 1
        if len(self.cazy_align_tools) > 1:
            self.on_rely(self.cazy_align_tools, self.run_cat_result)
        else:
            self.cazy_align_tools[0].on('end', self.run_cat_result)
        for tool in self.cazy_align_tools:
            tool.run()

    def run_cat_result(self):
        if len(self.cazy_align_tools) > 1:
            path =self.align_result_path
        else:
            path = self.cazy_align_tools[0].output_dir
        self.cat_out.set_options({
            "hmmscan_out": path,
        })
        self.cat_out.on('end', self.set_output, 'cazy_out')
        self.cat_out.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'cazy_align':
            file_path = obj.option('hmmscan_out_dm').prop['path']
            if os.path.exists(os.path.join(self.align_result_path, os.path.basename(file_path))):
                os.remove(os.path.join(self.align_result_path, os.path.basename(file_path)))
            os.link(file_path, os.path.join(self.align_result_path, os.path.basename(file_path)))  #####创建硬链
        if event['data'] == 'cazy_out':
            if os.path.exists(self.output_dir + '/align_result.txt'):
                os.remove(self.output_dir + '/align_result.txt')
            os.link(obj.option('hmmscan_result').prop['path'], self.output_dir + '/align_result.txt')
            self.option("align_result", self.output_dir + "/align_result.txt")
            self.end()

    def run(self):
        super(HmmscanModule, self).run()
        self.run_split_fasta()

    def end(self):
        super(HmmscanModule, self).end()

    def linkdir(self, dirpath, dirname):  # 暂时无用
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
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                file_name = os.listdir(oldfiles[i])
                os.mkdir(newfiles[i])
                for file_name_ in file_name:
                    os.link(os.path.join(oldfiles[i], file_name_), os.path.join(newfiles[i], file_name_))
