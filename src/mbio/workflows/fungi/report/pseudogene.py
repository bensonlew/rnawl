# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
# last_modifies 20180621


import os
from biocluster.workflow import Workflow
import datetime
from biocluster.api.file.remote import RemoteFileManager
from biocluster.config import Config
import re
from biocluster.file import exists
from biocluster.file import download

class PseudogeneWorkflow(Workflow):
    """
    供真菌基因组假基因预测
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(PseudogeneWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "specimen_list", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "main_name", "type": "string"},
            {"name": "ref_path", "type": "string", "default":""},
            {"name": "ref", "type": "string", "default":""},
            {"name": "update_info", "type": "string"},
            {"name": "middle_path", "type":"string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool_list = []
        self.config = Config()
        self.down_tool = self.add_tool('fungi_genome.download')
        self.down_list = []

    def run_pseudogene(self):

        self.pse_dir = {}
        for i in self.spe_list:
            self.pseudogene = self.add_tool("gene_structure.pseudogene")
            local_q = self.query_pre_path + '/' + i + "_scaf.fna"

            self.pseudogene.set_options({
                "query": local_q,
                "ref": self.ref,
                "sample": i
            })

            self.pse_dir[i] = self.pseudogene.output_dir

            self.tool_list.append(self.pseudogene)
        if len(self.tool_list) > 1:
            self.on_rely(self.tool_list, self.set_db)
            self.logger.info(self.tool_list)
        else:
            self.tool_list[0].on('end', self.set_db)

        for tool in self.tool_list:
            tool.run()


    def set_db(self):

        api_pseudogene = self.api.api('fungi_genome.pseudogene')

        id = self.option('main_id').split(",")
        dir_name = self.option('main_name').split(",")
        n = 0
        for i in self.spe_list:
            api_pseudogene.add_detail(self.pse_dir[i], main_id=id[n],specimen_id=i)

            self.linkdir(self.pse_dir[i],dir_name[n])
            n += 1
        self.end()

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.run_get_seq()
        super(PseudogeneWorkflow, self).run()


    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
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
            [".", "", "假基因预测结果目录",0,'140510'],
        ]
        regexps = [
            [r'.*_pseudogene.xls', 'xls', '假基因预测结果文件',0,'140511']
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(PseudogeneWorkflow, self).end()


    def get_seq_v1(self):

        self.spe_list = self.option('specimen_list').split(",")
        database = '/mnt/ilustre/tsanger-data/'
        if self._sheet.client == 'client01':
            database = '/mnt/ilustre/data/'
        assemble =os.path.join(database, 'rerewrweset', self.option('middle_path'), "{0}/assembly_predict/assembly/{0}_scaf.fna")
        self.query =[assemble.format(i) for i in self.spe_list]

        ##测试用：
        #self.ref = "/mnt/ilustre/users/sanger-dev/sg-users/zouguanqing/pseudogene/t300.faa"
        #self.query = ["/mnt/ilustre/users/sanger-dev/sg-users/zouguanqing/pseudogene/t1.fna"]

    def get_seq_v2(self):
        self.sanger_prefix = Config().get_project_region_bucket(project_type="fungi")
        self.spe_list = self.option('specimen_list').split(",")
        assemble =os.path.join(self.sanger_prefix,self.option('middle_path'),"{0}/assembly_predict/assembly/{0}_scaf.fna")
        self.query =[assemble.format(i) for i in self.spe_list]
        self.logger.info(self.query)

        for n , i in enumerate(self.spe_list):
            self.down_list.append([self.query[n], './'+ i + "_scaf.fna"])

    def run_get_seq(self):
        self.get_seq_v1()
        if exists(self.query[0]):
            for n , i in enumerate(self.spe_list):
                local_q = self.work_dir + '/seq/' + i + "_scaf.fna"
                if not os.path.exists(local_q):
                    download(self.query[n] ,local_q)
            self.query_pre_path = self.work_dir
        else:
            self.get_seq_v2()
            self.query_pre_path = self.down_tool.work_dir

        if self.option('ref') != "":
            self.ref = self.config.SOFTWARE_DIR + '/database/NCBI_fungi/faa/'+self.option('ref')+'.faa'
        else:
            base_name = self.option('ref_path').split('/')[-1]
            self.ref = os.path.join(self.down_tool.work_dir ,base_name)
            if '://' in self.option('ref_path'):
                self.down_list.append([self.option('ref_path'), './' + base_name])
            else:
                database = '/mnt/ilustre/tsanger-data/'
                if self._sheet.client == 'client01':
                    database = '/mnt/ilustre/data/'
                ref_tmp = database + '/rerewrweset/'+ self.option('ref_path').split('rerewrweset')[-1]
                download(ref_tmp, self.ref)
        self.down_tool.set_options({'item':str(self.down_list)})
        self.down_tool.on('end',self.run_pseudogene)
        self.down_tool.run()