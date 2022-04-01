# -*- coding: utf-8 -*-
# __author__ = "XueQinwen"
# last_modify: 20210908

import re
import os
import json
import time
import shutil
# from src.mbio.workflows.datasplit_v2.submit import Submit
# from submit import Submit
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
# from upload_s5cmd import UploadS5cmd

class PacbioSplitWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PacbioSplitWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "pacbio_params", "type": "infile", "format": "datasplit.library_params", "required": True},  # 进行三代拆分参数的参数文件
            {"name": "split_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.end_times = 0

    def check_options(self):
        if not self.option("pacbio_params").is_set:
            raise OptionError("缺少三代数据拆分参数,请检查")
        f = open(self.option("pacbio_params").prop["path"],"rb")
        try:
            json_dict = json.loads(f.read())
            self.params_json = json_dict["pacbio_split"]
        except:
            raise OptionError("JSON格式不正确")
    
    def run_ccs(self):
        # self.bam = ""
        self.ccs_workdir = os.path.join(self.work_dir,"ccs_bam")
        if not os.path.exists(self.ccs_workdir):
            os.mkdir(self.ccs_workdir)
        self.ccs_tools = []
        i = 1
        while i <= 10:
            ccs = self.add_tool("datasplit_v2.ccs")
            options = {
                "input_bam":self.new_bam_path,
                "chunk": str(i),
                "ccs_bam": os.path.join(self.ccs_workdir,"ccs{}.bam".format(i))
            }
            ccs.set_options(options)
            self.ccs_tools.append(ccs)
            i+=1
        if len(self.ccs_tools) == 1:
            self.ccs_tools[0].on('end',self.run_pbmerge)
        else:
            self.on_rely(self.ccs_tools,self.run_pbmerge)
        for tool in self.ccs_tools:
            tool.run()
        

    def run_pbmerge(self):
        self.pbmerge = self.add_tool("datasplit_v2.pbmerge")
        options = {
            "ccs_bams":self.ccs_workdir
        }
        self.pbmerge.set_options(options)
        self.pbmerge.on("end",self.run_lima)
        self.pbmerge.run()

    def run_lima(self):
        self.lima = self.add_tool("datasplit_v2.lima")
        if self.data_type == "subreads":
            options = {
                "ccs_bam":os.path.join(self.pbmerge.work_dir,"ccs.new.bam"),
                "barcode_path":self.barcode_path,
                "index_path":self.index_path
            }
        else:
            options = {
                "ccs_bam":self.new_bam_path,
                "barcode_path":self.barcode_path,
                "index_path":self.index_path
            }
        self.lima.set_options(options)
        self.lima.on("end",self.run_ccs_statistic)
        self.lima.run()

    def run_ccs_statistic(self):
        self.ccs_sta = self.add_tool("datasplit_v2.pacbio_stat")
        if self.data_type == "subreads":
            options = {
                "ccsbam":os.path.join(self.pbmerge.work_dir,"ccs.new.bam"),
                "split_dir":os.path.join(self.lima.output_dir,"split"),
                "origin_bam":self.new_bam_path,
                # "origin_bam":"/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/checkSmrtLink/ccs10.bam",
                "index_path":self.index_path
            }
        else:
            options = {
                "ccsbam":self.new_bam_path,
                "split_dir":os.path.join(self.lima.output_dir,"split"),
                "origin_bam":self.new_bam_path,
                "index_path":self.index_path
            }
        self.ccs_sta.set_options(options)
        self.ccs_sta.on("end",self.run_qc_stat)
        self.ccs_sta.run()

    def run_qc_stat(self):
        self.qc_stat = self.add_module("datasplit_v2.pacbio_qc_stat")
        options = {
            "split_dir":os.path.join(self.lima.output_dir,"split"),
            "index_path":self.index_path
        }
        self.qc_stat.set_options(options)
        self.qc_stat.on("end",self.set_output)
        self.qc_stat.run()

    def set_output(self):
        """
        pacbio_stat
        pacbio_qc_stat
        """
        self.linkdir(os.path.join(self.qc_stat.output_dir,""), 'PacbioQcStat')
        # self.linkdir(self.ccs_sta.output_dir,"PacbioStat")
        self.run_md5sum()
        # self.set_db()


    def linkdir(self, dirpath, dirname):
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath,i) for i in allfiles]
        newfiles = [os.path.join(newdir,i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i],newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.system('cp -r %s %s ' % (oldfiles[i], newdir))

    def set_db(self):
        s3_upload_dir = self._sheet.output
        app_pacbio = self.api.api("datasplit.pacbio_split")
        app_pacbio.update_sg_pacbio_specimen(self.option("split_id"),s3_upload_dir,os.path.join(self.output_dir,"PacbioQcStat/statistic"))
        app_pacbio.update_sg_pacbio(self.option("split_id"),os.path.join(self.ccs_sta.output_dir,os.path.basename(self.params_json["bam_path"]).split('.')[0]))
        app_pacbio.add_pacbio_bar(self.option("split_id"),os.path.join(self.ccs_sta.output_dir,"statReads.new.txt"))
        self.end()
        

    def get_index_table(self):
        csv_sheet = self.params_json["sample_sheet"]
        self.data_type = ""
        self.index_path = os.path.join(self.work_dir,"sample_list.txt")
        barcode = {}
        barcode_list = []
        with open(csv_sheet,'r') as cs, open(self.index_path,'w') as ip:
            line = cs.readline()
            title = line.rstrip().split(',')
            sn_info = {}
            temp_barcode = []
            while 1:
                line = cs.readline()
                if not line:
                    break
                fd = line.rstrip().split(',')
                sample_name = fd[1]
                for i in range(len(title)):
                    sn_info[title[i]] = fd[i]
                barcode[sn_info["f_name"]] = sn_info["f_barcode"]
                barcode[sn_info["r_name"]] = sn_info["r_barcode"]
                temp_barcode.append(sn_info["f_name"])
                temp_barcode.append(sn_info["r_name"])
                self.data_type = sn_info["data_type"]
                ip.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(sample_name,sn_info["f_name"],sn_info["r_name"],sn_info["type"],sn_info["primer_name"],sn_info["sample_name"]))
            barcode_list = sorted(set(temp_barcode))
        self.barcode_path = os.path.join(self.work_dir,"barcode.fasta")
        with open(self.barcode_path,"w") as bp:
            for i in barcode_list:
                bp.write(">{}\n".format(i))
                bp.write("{}\n".format(barcode[i]))

    def get_bam_pbi(self):
        if not os.path.exists(self.params_json["bam_path"]+".pbi"):
            origin_bam_name = os.path.basename(self.params_json["bam_path"])
            self.new_bam_path = os.path.join(self.work_dir,origin_bam_name)
            if not os.path.exists(self.new_bam_path):
                shutil.copyfile(self.params_json["bam_path"],self.new_bam_path)
            options = {
                "input_bam": self.new_bam_path
            }
            self.pbindex_tool = self.add_tool("datasplit_v2.pbindex")
            self.pbindex_tool.set_options(options)
            self.pbindex_tool.on("end",self.run_pacbio)
            self.pbindex_tool.run()
        else:
            self.new_bam_path = self.params_json["bam_path"]
            self.run_pacbio()

    def run_pacbio(self):
        self.get_index_table()
        if self.data_type == "subreads":
            self.run_ccs()
        elif self.data_type == "ccs":
            self.run_lima()
        else:
            self.set_error("三代拆分给的数据类型不对：{}".format(self.data_type))

    def run_md5sum(self):
        self.logger.info("开始进行md5校验")
        # self.md_tool = []
        # for f in os.listdir(self.output_dir):
        #     # if f in ["meta_qc", "dna_raw"]:
        #     if f in ["meta_raw", "meta_clean", "dna_raw"]:
        options = {"fastq_dir": os.path.join(self.output_dir, "PacbioQcStat/data")}
        self.md5sum = self.add_tool("datasplit_v2.md5sum")
        self.md5sum.set_options(options)
        self.md5sum.on("end", self.set_db)
        self.md5sum.run()
        # self.md_tool.append(self.md5sum)
        # self.logger.info(len(self.md_tool))
        # if len(self.md_tool) > 1:
        #     self.on_rely(self.md_tool, self.set_db)
        #     for md_tool in self.md_tool:
        #         md_tool.run()
        # elif len(self.md_tool) == 1:
        #     self.md_tool[0].on("end", self.set_db)
        #     self.md_tool[0].run()
        # else:
        # self.end()

    def run(self):
        self.get_bam_pbi()
        super(PacbioSplitWorkflow,self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".","","结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            ["","",""]
        ])
        super(PacbioSplitWorkflow,self).end()

                
            