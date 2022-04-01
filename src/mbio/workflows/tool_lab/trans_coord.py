# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import os
from biocluster.workflow import Workflow
from biocluster.config import Config
import datetime
import unittest
import types
from bson.objectid import ObjectId


class TransCoordWorkflow(Workflow):
    """
   基因组坐标转换：
   提供两种可选方式：
   方式1：
   当前物种:others
   此时允许客户上传当前物种fasta文件和目前物种fasta文件以及基因组坐标转换文件
   通过序列比对和多个ucsc小工具生成fasta文件之间相关的chain文件，后针对基因组坐标转换文件进行转换
   方式2：
   针对已有物种和已知物种版本的基因组之间进行基因组坐标转换
   chain文件从ucsc上下载
   基因组坐标文件格式要求如下
    NW_101466.1     4157    4281    1
    NW_101466.1     4000    4103    2
    NW_101466.1     3366    3937    3
    NW_101466.1     2632    3320    4

    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TransCoordWorkflow, self).__init__(wsheet_object)
        options = [
            # 选择当前物种["Human","Mouse","Zebrafish","C.elegans","D.melanogaster","Pig","Others"]
            {"name": "current_species", "type": "string"},
            #当前物种选择others时，由客户上传fasta文件(当前物种)
            {"name": "current_fasta", "type": "infile", "format": "ref_rna_v2.common"},
            # 当前物种选择others时，由客户上传fasta文件(目标物种)
            {"name": "target_fasta", "type": "infile", "format": "ref_rna_v2.common"},
            # 无论哪种模式均需要上传待转化bed文件
            {"name": "coordinate_file", "type": "infile", "format": "ref_rna_v2.common"},  # fasta文件
            #当前版本
            #[{"Human":"GRCh38 / hg38","GRCh37 / hg19"},{"Mouse":["GRCm38 / mm10",""NCBI37 / mm9]},"Zebrafish":["GRCz11 / danRer11","GRCz10 / danRer10"],
            #{"C.elegans":[" WS220 / ce11","WS220 / ce10"]},{"D.melanogaster":["BDGP Release 6 + ISO1 MT / dm6"]},{"Pig",["SGSC Sscrofa10.2 / susScr3","SGSC Sscrofa11.1 / susScr11"]}
            {"name": "current_version", "type": "string"},
            #目标物种
            {"name": "target_species", "type": "string", "default": "1"},
            #目标物种版本
            # [{"Human":"GRCh38 / hg38","GRCh37 / hg19"},{"Mouse":["GRCm38 / mm10",""NCBI37 / mm9]},"Zebrafish":["GRCz11 / danRer11","GRCz10 / danRer10"],
            # {"C.elegans":[" WS220 / ce11","WS220 / ce10"]},{"D.melanogaster":["BDGP Release 6 + ISO1 MT / dm6"]},{"Pig",["SGSC Sscrofa10.2 / susScr3","SGSC Sscrofa11.1 / susScr11"]}
            {"name": "target_version", "type": "string", "default": "5"},
            {"name": "main_id", "type": "string",},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.trans_coord.lift_over")
        self.set_options(self._sheet.options())
        if self.option("current_species").lower() != "others":
            self.map_chain = os.path.join(self.config.SOFTWARE_DIR, 'database/lnc_rna/liftOver/{}To{}{}.over.chain.gz'.format(self.option('current_version'),self.option('target_version')[0].upper(),self.option('target_version')[1:]))
            self.logger.info("{}".format(self.map_chain))


    def run(self):
        if self.option("current_species").lower() == "others":
            self.run_custom()
        else:
            self.run_liftover()

        super(TransCoordWorkflow, self).run()

    def run_custom(self):
        self.chain_create=self.add_module("tool_lab.trans_coord")
        opts = {
            'raw_fasta': self.option("current_fasta"),
            'target_fasta': self.option('target_fasta'),
        }
        self.chain_create.set_options(opts)
        self.chain_create.on("end", self.run_liftover)
        self.chain_create.run()


    def run_liftover(self):
        if self.option("current_species").lower() == "others":
            self.map_chain=os.path.join(self.chain_create.output_dir,"lift.chain")
        opts = {
            'coordinate_file': self.option("coordinate_file").prop["path"],
            'chain_file': self.map_chain,
        }
        self.tool.set_options(opts)
        self.tool.on('end', self.set_output)
        self.tool.run()

    def set_output(self):
        for file in os.listdir(self.tool.output_dir):
            os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "基因组坐标转换转换结果文件",0],
            # [r'.*\.cpm2tpm\.xls', 'xls', '定量指标cpm转tpm结果文件', 0],
            # [r'.*\.count2tpm\.xls', 'xls', '定量指标count转tpm结果文件', 0],
            # [r'.*\.count2cpm\.xls', 'xls', '定量指标count转cpm结果文件', 0],
            # [r'.*\.count2fpkm\.xls', 'xls', '定量指标count转fpkm结果文件', 0],
            # [r'.*\.fpkm2tpm\.xls', 'xls', '定量指标fpkm转tpm结果文件', 0],
            # [r'.*\.cpm2fpkm\.xls', 'xls', '定量指标cpm转fpkm结果文件', 0],
            # [r'.*\.fpkm2cpm\.xls', 'xls', '定量指标fpkm转cpm结果文件', 0],
        ])
        super(TransCoordWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.trans_coord import TransCoordWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "trans_coord" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.trans_coord",
            "options": dict(
                # genome=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                # fasta_seq=">scaffold1\nATTATAAAATGAGGAGATGTAAATTTTAAGGAAAAAAATAAAGCTGTCGAGATTTTCTCGACAGCCTGAAGAGCACCCATCATAAAAGGGTGCTCTTTTTCTTTTCTCTTCCCGATTTGTTCACAGGCTGAATCCTCTCCTCATATGCTGAAAGGGAGATTCAGAAAATATTACAGACTTTTTTAAATAAAGAAGGTGAACGGTCAGTATGTTGAGCGGTTTAACGGTTGCGGTGATCGGGGGAGATGCAAGGCAGCTTGAAATCATTCGCAAGCTGTCACAGCAGCATGCCAAAGTGTTTTTGGTCGGATTTGATCAGCTGGATCATGGGTTTATCGGTGCTGAAAAGCTTAAAATGTCAGAACTTCCATTTGAACAGGTAGACAGTATGATTCTGCCGGTATCAGGTGCAACAGATGAAGGCGTCGTCGCCACAGTTTTCTCAAATGAGCAGGTCGTGCTGGAAGCAGAATATTTAGAAAGAACTCCAGCACATTGTACCTTGTACTCAGGTATTTCTAATACGTACTTAGACAATCTGGCAAAGCAGGTGAACCGGAAGCTTGTGAAGCTGTTTGAGCGCGATGATATTGCCATATATAACTCTATTCCAACAGTTGAAGGGATTATCATGATGGCCATTCAGCAAACGGACTATACGATTCATGGATCACATGTCGCTGTCCTCGGGCTTGGGAGAACAGGGCTCACAATTGCCCGCACAT",
                # current_species="Human",
                # current_version="hg38",
                # target_species="Human",
                # target_version='hg19',
                # coordinate_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/TransCoord/input/test.bed"
                current_fasta="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/TransCoord/input/ASM14920v2.fasta",
                target_fasta="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/TransCoord/input/ASM18445v3.fasta",
                current_species="Others",
                coordinate_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/TransCoord/input/lncRNA.bed"
            )
        }
        wsheet = Sheet(data=data)
        wf =TransCoordWorkflow(wsheet)
        wf.sheet.id = 'trans_coord'
        wf.sheet.project_sn = 'trans_coord'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
