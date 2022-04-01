# -*- coding: utf-8 -*-
# __author__ = '...'
#__modified__ = 'moli.zhou'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError


class RnaeditingAgent(Agent):
    def __init__(self, parent):
        super(RnaeditingAgent, self).__init__(parent)
        options = [
            {"name": "rna_bam_dir", "type": "infile", "format": "align.bwa.bam_dir"},
            {"name": "ref_hg19.fa", "type": "infile", "format": "sequence.fasta"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("rna_bam_dir").is_set:
            raise OptionError("必须提供bam文件！")
        if not self.option("ref_hg19.fa").is_set:
            raise OptionError("请提供参考基因组文件")
        return True

    def set_resource(self):
        self._cpu = 10
        self._memory = '200G'

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
            ["./RDD_output.VafList.txt", "txt", "Vaf表"],
        ])
        print self.get_upload_files()
        super(RnaeditingAgent, self).end()


class RnaeditingTool(Tool):
    def __init__(self, config):
        super(RnaeditingTool, self).__init__(config)
        self._version = "1.1"
        self.cmd_path = self.config.SOFTWARE_DIR + '/bioinfo/rna/RDDpred_v1.1.Dir/RDDpred.py'
        self.script_path = self.config.SOFTWARE_DIR + '/bioinfo/rna/scripts'
        self.out_files = ["RDD_info", "RDD_output.Argument.Log", "RDD_output.Step2.Rescue.GenomicList.txt", "RDD_output.Step2.BamList.txt", "RDD_output.Step2.BamHeader.txt", "RDD_output.RDD.RawList.txt", "RDD_output.Step4.Predictor.Training.Log", "RDD_output.Step4.Attribute.Evaluation.Log", "RDD_output.Step5.Results.Summary.txt", "RDD_output.Prediction.ResultList.txt", "RDD_output.ModelDir", "RDD_output.RDDpred.results_report.txt", "RDD_output.VafList.txt"]
        self.remove_files = ["RDD_output.Step2.BamList.txt", "RDD_output.Step2.BamHeader.txt", "RDD_output.RDD.RawList.txt", "RDD_output.Step4.Predictor.Training.Log", "RDD_output.Step4.Predictor.Training.Log", "RDD_output.ModelDir", "RDD_info", "RDD_output.Step2.Rescue.GenomicList.txt", "RDD_output.Argument.Log", "RDD_output.Step4.Attribute.Evaluation.Log"]

    def run(self):
        super(RnaeditingTool, self).run()
        self.run_rna_editing_py()
        self.reassign_run()
        self.end()

    def run_rna_editing_py(self):
        if not os.path.exists("RNA_Editing"):
            os.mkdir("RNA_Editing")
        with open("bam_list.txt", "wb") as bam_list:
            for bam_files in os.listdir(self.option("rna_bam_dir").prop["path"]):
                bam_list.write(os.path.join(self.option("rna_bam_dir").prop["path"], bam_files) + "\n")

        cmd = 'program/Python/bin/python {} -rbl {} -rsf {} -tdp {} -ops {} -psl {} ' \
               '-nsl {}'.format(self.cmd_path, os.path.join(self.work_dir, "bam_list.txt"),
                                self.option("ref_hg19.fa").prop["path"],
                                self.config.SOFTWARE_DIR + '/bioinfo/rna/RDDpred_v1.1.Dir/ToolBox.Dir',
                                os.path.join(self.work_dir, "output/RDD_prefix"),
                                self.config.SOFTWARE_DIR + '/bioinfo/rna/RDDpred_v1.1.Dir/PriorData/hg19.PublicSites.txt',
                                self.config.SOFTWARE_DIR + '/bioinfo/rna/RDDpred_v1.1.Dir/PriorData/hg19.MES_Sites.txt')

        self.logger.info('开始运行RDDpred并检测编辑位点')
        self.logger.debug(cmd)

        cmd_obj = self.add_command("rddpred_cmd", cmd)
        cmd_obj.run()
        self.wait()
        if cmd_obj.return_code == 0:
            self.logger.info('RDDpred检测完成')
            self.result = os.path.join(self.work_dir,'output/RDD_prefix.RDDpred.results_report.txt')
        else:
            self.logger.error('RDDpred检测失败')
            self.set_error("运行RDDpred.py失败")
            raise Exception("运行RDDpred.py失败")

    def reassign_run(self):
        cmd = 'program/Python/bin/python {} {} {}'.format(os.path.join(self.script_path,'rna_edit_reassign.py'), self.result, self.output_dir)

        self.logger.info('开始运行RDDpred并检测编辑位点')
        self.logger.debug(cmd)

        cmd_obj = self.add_command("reassign_cmd", cmd)
        cmd_obj.run()
        self.wait()
        if cmd_obj.return_code == 0:
            self.logger.info('重排完成')
        else:
            self.logger.error('重排失败')
            self.set_error("运行rna_edit_reassign.py失败")
            raise Exception("运行rna_edit_reassign.py失败")


    def linkfile(self, oldfile, newname):
        newpath = os.path.join(self.output_dir, newname)
        if os.path.exists(newpath):
            os.remove(newpath)
        os.link(oldfile, newpath)
