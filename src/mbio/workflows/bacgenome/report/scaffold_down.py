# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modifies 20190320


import os,re
from biocluster.workflow import Workflow
import datetime
from biocluster.api.file.remote import RemoteFileManager
import shutil
from biocluster.file import download
from biocluster.api.file.lib.transfer import MultiFileTransfer


class ScaffoldDownWorkflow(Workflow):
    """
    供细菌基因组组装部分的序列下载接口
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(ScaffoldDownWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "main_table_name", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "client", "type": "string"},
            {"name": "type", "type": "string"},
            {"name": "seq_path", "type": "string"},
            {"name": "specimen_name", "type": "string"},
            # 'GBB01:aaab;GBB02:bbbc;GBB03:CCC'
            {"name": "seq_list", "type": "string"},
            # GBB01:Chromosome;GBB02:Chromosome;GBB03:Chromos
            {"name": "params", "type": "string"},
            {"name": "seq_type", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool_list = []
        self.file_path = self._sheet.output

    def run_seq_down(self):
        transfer = MultiFileTransfer()
        self.down_dir = {}
        specimen = self.get_dic(self.option("specimen_name"))
        specimen_list = self.get_dic(self.option("seq_list"))
        sample_seq_path = eval(self.option('seq_path'))
        for k in sorted(specimen_list.keys()):
            self.logger.info(sample_seq_path[k])
            if self.option('type') in ['uncomplete']:
                self.down = self.add_tool("bacgenome.extract_seq_ass")
                seq_path = sample_seq_path[k]
                options = {
                    'seq_path': seq_path,
                    'seq_list': specimen_list[k],
                    'sample_new': specimen[k],
                    'seq_type': self.option('seq_type')
                }
                self.down.set_options(options)
                self.down_dir[k] = self.down.output_dir
                self.tool_list.append(self.down)

            elif self.option('type') in ['complete']:
                self.down = self.add_tool("bacgenome.extract_seq_assemble")
                transfer.add_download(sample_seq_path[k], self.work_dir + '/' + k + '/')
                transfer.perform()
                seq_path = self.work_dir + '/' + k + '/'
                options = {
                    'seq_path': seq_path,
                    'type': self.option('type'),
                    'seq_list': specimen_list[k],
                    'sample': k,
                    'sample_new': specimen[k],
                }
                self.down.set_options(options)
                self.down_dir[k] = self.down.output_dir
                self.tool_list.append(self.down)
        self.on_rely(self.tool_list, self.set_db)
        for tool in self.tool_list:
            tool.run()

    def run(self):
        self.run_seq_down()
        super(ScaffoldDownWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        main_id = self.option('main_id')
        specimen_list = self.get_dic(self.option("seq_list"))
        specimen = self.get_dic(self.option("specimen_name"))
        if self.option('type') in ['uncomplete']:
            for k in sorted(specimen_list.keys()):
                self.linkdir(self.down_dir[k], self.output_dir)
            if len(os.listdir(self.output_dir)) > 1:
                file_name = " ".join(os.listdir(self.output_dir))
                os.system("cd " + self.output_dir + ";tar -czvf sequence_" + datetime.datetime.now().strftime(
                        "%Y%m%d_%H%M%S%f")[:-3] + ".tar.gz " + file_name)
                os.system("cd " + self.output_dir + ";rm -rf " + file_name)
        elif self.option('type') in ['complete']:
            self.logger.info(sorted(specimen_list.keys()))
            for k in sorted(specimen_list.keys()):
                self.logger.info(self.down_dir[k] + '/' + specimen[k])
                if os.path.exists(self.output_dir + '/' + specimen[k]):
                    shutil.rmtree(self.output_dir + '/' + specimen[k])
                shutil.copytree(self.down_dir[k] + '/' + specimen[k], self.output_dir + '/' + specimen[k])
                self.logger.info(self.output_dir)
            self.logger.info("aaaaaaaaaa")
            file_name = " ".join(os.listdir(self.output_dir))
            self.logger.info(file_name)
            os.system("cd " + self.output_dir + ";tar -zcvf sequence_" + datetime.datetime.now().strftime(
                "%Y%m%d_%H%M%S%f")[:-3] + ".tar.gz " + file_name)
            os.system("cd " + self.output_dir + ";rm -rf " + file_name)
        api_seqdown = self.api.api('bacgenome.seq_down')
        self.logger.info('dao biao ')
        self.logger.info(self.option('main_id'))
        file = os.listdir(self.output_dir)[0]
        self.logger.info(self.file_path + file)
        api_seqdown.add_seq_down(main_id=main_id, link_path=self.file_path + file)
        self.end()

    def end(self):
        repaths = [
            [".", "", ""],
        ]
        regexps = [
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(ScaffoldDownWorkflow, self).end()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件的output目录
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        self.logger.info(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])

    def get_dic(self, file):
        dict = {}
        list = file.split(';')
        for i in list:
            arry = i.split(':')
            dict[arry[0]] = arry[1]
        return dict

    def get_path_dic(self, file):
        dict = {}
        list = file.split(';')
        for i in list:
            arry = i.split(':')
            dict[arry[0]] = arry[1] + ':' + arry[2]
        return dict