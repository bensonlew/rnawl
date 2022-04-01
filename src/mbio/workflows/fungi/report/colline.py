# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
# last_modifies 20180621


import os
from biocluster.workflow import Workflow
import datetime
from biocluster.api.file.remote import RemoteFileManager
from biocluster.config import Config
import re
import shutil
from biocluster.file import exists
from biocluster.file import download

class CollineWorkflow(Workflow):
    """
    供真菌基因组共性分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(CollineWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "specimen_id", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "main_name", "type": "string"},
            {"name": "ref_path", "type": "string", "default" :""},
            {"name": "ref_id", "type": "string", "default" :""},
            {"name": "genbank", "type":"string", "default" :""},
            {"name": "update_info", "type": "string"},
            {"name": "middle_path", "type":"string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool_list = []
        self.config = Config()
        self.colline = self.add_tool("fungi_genome.colline")
        self.down_tool = self.add_tool('fungi_genome.download')
        self.down_list = []


    def run_colline(self):
        ##快速测试用
        #self.query = '/mnt/ilustre/users/sanger-dev/sg-users/zouguanqing/colline/R01.close.300.fna.all'
        #self.ref = '/mnt/ilustre/users/sanger-dev/sg-users/zouguanqing/colline/GCF_000002495.2_MG8_genomic.fna.all'
        if not os.path.getsize(self.query):
            self.logger.set_error('query is null')
        self.colline.set_options({
            "qfna": self.query,
            "rfna": self.ref,
            "sample": self.option('specimen_id')
            })
        self.colline.on('end', self.set_db)
        self.colline.run()


    def set_db(self):
        api_colline = self.api.api('fungi_genome.colline')
        ref = self.ref.split('/')[-1].split('.')[0]
        #data_path = '/fungi/' + self.option('task_id') + '/colline/' + self.option('main_id')
        data_path = self._sheet.output + '/circos.png'
        api_colline.add_detail(self.colline.output_dir, main_id=self.option('main_id'),specimen_id=self.option('specimen_id'),ref=ref,data_path=data_path)
        self.end()

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.run_get_seq()
        #self.run_colline()
        super(CollineWorkflow, self).run()

    def end(self):
        repaths = [
            [".", "", "共性性结果目录",0,'140506'],
        ]
        regexps = [
            [r'.*_blocks_coords.txt', 'txt', '结果文件',0,'140501'],
            [r'.*_coverage_report.txt', 'txt','结果文件',0,'140502'],
            [r'.*circos.png', 'png','结果文件',0,'140503'],
            [r'.*circos.svg', 'svg','结果文件',0,'140504'],
            [r'.*block_new.xls', 'xls','结果文件',0,'140505']
        ]
        sdir = self.add_upload_dir(self.colline.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(CollineWorkflow, self).end()

    def get_seq_path(self, file_path):
        client = Config().get_mongo_client('fungigenome')
        db = client[Config().get_mongo_dbname('fungigenome')]
        try:
            ass_info = db['assemble'].find_one({'task_id': self.option('task_id')})
            if 'origin_id' in ass_info:
                ass_id = ass_info['origin_id']
            else:
                ass_id = ass_info['_id']
            seq_info = db['assemble_seq'].find_one({'assemble_id': ass_id})
            file_path = seq_info['seq_path'] +\
                "/{0}/assembly_predict/assembly/{0}_scaf.fna".format(self.option("specimen_id"))
            return file_path
        except Exception as e:
            self.set_error("找不到qeury文件:{}\n{}".format(file_path, e))

    def run_get_seq(self):
        #transfer = MultiFileTransfer()
        pre_assemble_file =self.option('middle_path')
        self.sanger_prefix = Config().get_project_region_bucket(project_type="fungi")
        self.database = '/mnt/ilustre/tsanger-data/'
        if self._sheet.client == 'client01':
            self.database = '/mnt/ilustre/data/'
        self.query_1 = self.database +'rerewrweset/' + self.option('middle_path') + "/{0}/assembly_predict/assembly/{0}_scaf.fna".format(self.option("specimen_id"))
        self.logger.info("self.query_1: " + self.query_1 )
        if exists(self.query_1):
            self.query = self.work_dir + "/{0}_scaf.fna".format(self.option("specimen_id"))
            download(self.query_1, self.query)
            #transfer.add_download(self.query_1, self.query)
            #transfer.perform()
        else:
            query_tmp = pre_assemble_file + "/{0}/assembly_predict/assembly/{0}_scaf.fna".format(self.option("specimen_id"))
            query_base_name = query_tmp.split('/')[-1]
            self.query = os.path.join(self.down_tool.work_dir, query_base_name)
            remote_query_path = os.path.join(self.sanger_prefix, query_tmp)
            if not exists(remote_query_path):
                remote_query_path = self.get_seq_path(remote_query_path)
            self.down_list.append([remote_query_path, './' + query_base_name])

        if self.option("genbank")!= "":
            self.ref = self.config.SOFTWARE_DIR + '/database/NCBI_fungi/fna/'+self.option('genbank')+'.fna'
        elif self.option("ref_path")!= "" :
            ref_base_name = self.option('ref_path').split('/')[-1]
            if '://' in self.option('ref_path'):
                self.ref = os.path.join(self.down_tool.work_dir, ref_base_name)
                self.down_list.append([self.option('ref_path'), './' + ref_base_name])
            else:
                self.ref_1 = self.database +'/'+ self.option("ref_path").split(':')[-1]
                self.ref = self.work_dir + '/' + ref_base_name
                download(self.ref_1, self.ref)
        else:
            self.ref_1 = self.database + 'rerewrweset/' +  pre_assemble_file + "/{0}/assembly_predict/assembly/{0}_scaf.fna".format(self.option("ref_id"))
            if exists(self.ref_1):
                self.ref = self.work_dir + '/{0}_scaf.fna'.format(self.option("ref_id"))
                download(self.ref_1, self.ref)
            else:
                ref_tmp = pre_assemble_file + "/{0}/assembly_predict/assembly/{0}_scaf.fna".format(self.option("ref_id"))
                ref_base_name = ref_tmp.split('/')[-1]
                self.ref = os.path.join(self.down_tool.work_dir, ref_base_name)
                self.down_list.append([os.path.join(self.sanger_prefix, ref_tmp), './' + ref_base_name])
        self.logger.info("self.down_lsit: {}".format(self.down_list))

        self.down_tool.set_options({'item':str(self.down_list)})
        self.down_tool.on('end',self.run_colline)
        self.down_tool.run()
