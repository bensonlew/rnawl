# -*- coding: utf-8 -*-
# __author__ = '...'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import types
import subprocess
import shutil
from biocluster.core.exceptions import OptionError


class RnaeditingAgent(Agent):
    
    def __init__(self, parent):
        super(RnaeditingAgent, self).__init__(parent)
        options = [
            {"name": "rna_bam_dir", "type": "infile","format":"align.bwa.bam_dir"},
            {"name": "ref_hg19.fa", "type": "infile", "format": "sequence.fasta"}
            ]
        self.add_option(options)
        self.on('start', self.step_start)
        self.on('end', self.step_end)
      
    def check_options(self):
        if not self.option("rna_bam_dir").is_set:
            raise OptionError("必须提供bam文件！")
        if not self.option("ref_hg19.fa").is_set:
            raise OptionError("请提供参考基因组文件")
        return True 
        
    def set_resource(self):
        
        self._cpu = 10
        self._memory = '500G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
                [".", "", "RNA编辑结果输出目录"],
                ["./RDD_output.Step2.BamList.txt", "txt", "RNA编辑输入bam文件列表"],
                ["./RDD_output.Step2.BamHeader.txt", "txt", "RNA编辑输入bam头文件"],
                ["./RDD_output.RDD.RawList.txt", "txt", "RNA编辑位点原始列表"],
                ["./RDD_output.Step4.Predictor.Training.Log", "log", "训练日志"],
                ["./RDD_output.Step4.Predictor.Training.Log", "log", "evaluation日志"],
                ["./RDD_output.Step5.Results.Summary.txt", "txt", "RNA编辑位点概览表"],
                ["./RDD_output.RDDpred.results_report.txt", "txt", "RNA编辑位点详表"],
                ["./RDD_output.VafList.txt","txt","Vaf表"],
                ])
        print self.get_upload_files()
        super(RnaeditingAgent, self).end()

class RnaeditingTool(Tool):
    def __init__(self, config):
        super(RnaeditingTool, self).__init__(config)
        self._version = "1.1"
        self.cmd_path = self.config.SOFTWARE_DIR + '/bioinfo/rna/RDDpred_v1.1/RDDpred.py'
         
    def run(self):
        
        super(RnaeditingTool, self).run()
        self.run_rna_editing_py()

    def run_rna_editing_py(self):
        txt = os.path.join(self.work_dir,)
        if not os.path.exists("bam_list.txt"):
            os.mknod("bam_list.txt")
    
        with open("self.work_dir/bam_list.txt","w+") as bam_list:
            for bam_files in os.listdir(self.options("rna_bam_dir")):
                bam_list.write(str(bam_files)+"\n")
        
        
        cmd = self.config.SOFTWARE_DIR + '/program/Python/bin/python '
        cmd += self.cmd_path
        cmd += ' -rbl %s -rsf %s -tdp %s -ops %s -psl %s -nsl %s' % (self.work_dir/bam_list.txt,self.options("ref_hg19.fa").prop["path"],self.config.SOFTWARE_DIR + '/bioinfo/rna/RDDpred_v1.1/ToolBox.Dir',self.work_dir/RNA_Editing/RDD_prefix,self.config.SOFTWARE_DIR + '/bioinfo/rna/RDDpred_v1.1/PriorData/hg19.PublicSites.txt',self.config.SOFTWARE_DIR + '/bioinfo/rna/RDDpred_v1.1/PriorData/hg19.MES_Sites.txt')
        self.logger.info('开始运行RDDpred并检测编辑位点')
        
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('RDDpred检测完成')
        except subprocess.CalledProcessError:
            self.logger.info('RDDpred检测失败')
            self.set_error('运行RDDpred.py失败')
        self.end()
        
        
