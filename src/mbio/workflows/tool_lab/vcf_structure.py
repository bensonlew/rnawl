# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = "XueQinwen"

import os
import re 
import datetime
from biocluster.workflow import Workflow
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.core.exceptions import OptionError

class VcfStructureWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(VcfStructureWorkflow, self).__init__(wsheet_object)
        options = [
            {"name":"vcf_file","type":"infile","format":"dna_gmap.vcf"},
            {"name":"k_min","type":"int","default":2},
            {"name":"k_max","type":"int","default":20},
            {"name":"max_missing","type":"float","default":0.3},
            {"name":"maf","type":"float","default":0.05},
            {"name":"main_id","type":"string"},
            {"name":"update_info","type":"string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.filter = self.add_tool("tool_lab.vcf_filter")
        self.plink = self.add_tool("tool_lab.vcf_plink")
        self.structure = self.add_module('tool_lab.pop_structure')
    
    def check_option(self):
        """
        参数检查
        """
        if not self.option("vcf_file"):
            raise OptionError("必须输入vcf文件")

    def run_filter(self):
        self.filter.set_options({
            "vcf_file":self.option("vcf_file"),
            "max_missing":self.option("max_missing"),
            "maf":self.option("maf")
        })
        self.filter.on('end',self.run_plink)
        self.filter.run()

    def run_plink(self):
        self.filter_vcf = os.path.join(self.filter.output_dir,"pop.recode.vcf")
        self.plink.set_options({
            "vcf_file":self.filter_vcf,
        })
        self.plink.on("end",self.run_structure)
        self.plink.run()
    
    def run_structure(self):
        fam = os.path.join(self.plink.output_dir,"pop.fam")
        bed = os.path.join(self.plink.output_dir,"pop.bed")
        self.structure.set_options({
            "pop_fam":fam,
            "pop_bed":bed,
            "k_min":self.option("k_min"),
            "k_max":self.option("k_max")
        })
        self.structure.on("end",self.set_output)
        self.structure.run()
    
    def set_output(self):
        self.linkdir(self.structure.output_dir, 'structure')
        self.set_db()
    
    def linkdir(self, dirpath, dirname):
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
                    # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_db(self):
        self.logger.info("开始导表")
        api_structure = self.api.api("tool_lab.vcf_structure")
        for root,_,files in os.walk(os.path.join(self.structure.output_dir,"structure")):
            for name in files:
                if len(name.split(".")) == 3:
                    k  = name.split(".")[1]
                    if name.split(".")[2] == "log":
                        continue
                else:
                    continue
                api_structure.add_structrue(self.option("main_id"),k,
                                os.path.join(root,name))
        api_structure.add_cverror_detail(self.option("main_id"),
                                os.path.join(self.structure.output_dir,"cverror/cv.error"))
        self.logger.info("导表结束")
        self.end()

    def run(self):
        self.run_filter()
        super(VcfStructureWorkflow, self).run()
    
    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(VcfStructureWorkflow, self).end()

if __name__ == "__main__":
    from biocluster.wsheet import Sheet
    import random
    add_time = datetime.datetime.now().strftime("%Y%m%d%H%M%S%f") +\
        str(random.randint(1000,10000))
    data = {
        "name":"test_vcf_structure",
        "id":"vcf_structure_" + str(random.randint(1,10000)),
        "type":"workflow",
        "options":{
            "vcf_file": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/test_WGS/4.structure/pop.recode.vcf",
            "k_max": 5,
            "main_id":add_time
        }

    }
    wsheet = Sheet(data=data)
    wf = VcfStructureWorkflow(wsheet)
    wf.run()
