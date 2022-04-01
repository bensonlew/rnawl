# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last_modifies 20180326


import os,re
from biocluster.workflow import Workflow
import datetime
from biocluster.api.file.remote import RemoteFileManager
import shutil
from biocluster.config import Config
from biocluster.file import download
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.file import getsize, exists

class SequenceDownWorkflow(Workflow):
    """
    供细菌基因组预测部分的序列下载接口
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(SequenceDownWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "client", "type": "string"},
            {"name": "seq_path", "type": "string"},  # 预测结果的序列和二级机构文件
            {"name": "type", "type": "string"},
            {"name": "specimen_gene", "type": "string"},
            # '{"B": "gene001,gene002,gene003", "A": "gene001,gene002,gene003", "C": "gene001,gene002,gene003"}'
            {"name": "specimen_gene_new", "type": "string"},
            # '{"B": {"B":"B_new","gene": "gene_new"}, "A": {"A":"A_new","gene": "gene_new"}, "C": {"C":"C_new","gene": "gene_new"}}'
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool_list = []
        self.file_path = self._sheet.output

    def run_seq_down(self):
        transfer = MultiFileTransfer()
        self.down_dir = {}
        specimen_gene = eval(self.option("specimen_gene"))
        specimen_gene_new = eval(self.option("specimen_gene_new"))
        for k in sorted(specimen_gene.keys()):
            self.logger.info(k)
            self.logger.info(self.option('seq_path'))
            sample_seq_path = eval(self.option('seq_path'))
            self.down = self.add_tool("bacgenome.extract_seq_bygenenu")
            if re.search(r'://', sample_seq_path[0]):
                seq_file = sample_seq_path[0] + k + sample_seq_path[1]
                self.logger.info(seq_file)
                transfer.add_download(seq_file, self.work_dir + '/' + self.option('type') + '/' + k + '/')
                transfer.perform()
                seq_path = self.work_dir + '/' + self.option('type') + '/' + k + '/' + k + sample_seq_path[2]
            else:
                seq_file = sample_seq_path[0] + k + sample_seq_path[1]
                self.logger.info(seq_file)
                transfer.add_download(seq_file, self.work_dir + '/' + self.option('type') + '/' + k + '/')
                transfer.perform()
                seq_path = self.work_dir + '/' + self.option('type') + '/' + k + '/' + k + sample_seq_path[2]
            options = {
                'seq_prefix_path': seq_path,
                'type': self.option('type'),
                'gene_list': specimen_gene[k],
                'sample': k,
                'sample_gene_new': str(specimen_gene_new[k])
            }
            self.down.set_options(options)
            self.down_dir[k] = self.down.output_dir
            self.tool_list.append(self.down)
        self.on_rely(self.tool_list, self.set_db)
        for tool in self.tool_list:
            tool.run()

    def run(self):
        self.run_seq_down()
        super(SequenceDownWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        main_id = self.option('main_id')
        specimen_gene = eval(self.option("specimen_gene"))
        for k in sorted(specimen_gene.keys()):
            self.linkdir(self.down_dir[k], self.output_dir)
        if len(os.listdir(self.output_dir)) > 1:
            file_name = " ".join(os.listdir(self.output_dir))
            os.system(
                "cd " + self.output_dir + ";tar czvf " + self.option('type') + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[
                                                                  :-3] + ".tar.gz " + file_name)
            os.system("cd " + self.output_dir + ";rm -rf " + file_name)
        api_seqdown = self.api.api('bacgenome.seq_down')
        self.logger.info('dao biao ')
        self.logger.info(self.option('main_id'))
        file =os.listdir(self.output_dir)[0]
        self.logger.info(self.file_path + file)
        api_seqdown.add_seq_down(main_id=main_id,link_path=self.file_path + file)
        self.end()

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

    def end(self):
        repaths = [
            [".", "", ""],
        ]
        regexps = [
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(SequenceDownWorkflow, self).end()
