# -*- coding: utf-8 -*-
# __author__ = 'ysh'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os, shutil
import re
import unittest
import pandas as pd

class MethylationAgent(Agent):
    """
    三代数据甲基化位点预测
    """

    def __init__(self, parent):
        super(MethylationAgent, self).__init__(parent)
        options = [
            {"name": "input", "type": "string"}, # 多文件逗号分隔
            #{"name": "ref_input", "type": "infile", "format": "bacgenome.methy_file","check":"check_ref_file"},
            {"name": "ref_input", "type": "string", "default": ""},
            {"name": "result", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "input_type", "type": "string"},
            {"name": "sample", "type": "string"} , # 样品名称
            {"name": "ref_type", "type":"string","default":"fasta"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("input"):
            raise OptionError("必须设置输入文件，bam或者h5文件")
        if not self.option("input_type"):
            raise OptionError("必须设置输入文件类型")
        if not self.option("sample"):
            raise OptionError("必须设置输入样本名称")
        if not self.option("ref_input"):
            raise OptionError("必须设置输入参考序列文件")
        return True

    def set_resource(self):
        self._cpu = 16
        self._memory = '90G'

    def end(self):
        super(MethylationAgent, self).end()


class MethylationTool(Tool):
    def __init__(self, config):
        super(MethylationTool, self).__init__(config)
        self._version = "1.0"
        #self.pbsmrtpipe = self.config.SOFTWARE_DIR + '/bioinfo//Genomic/Sofware/smrttools/smrtcmds/bin/pbsmrtpipe'
        self.pbsmrtpipe = '/bioinfo//Genomic/Sofware/smrttools/smrtcmds/bin/pbsmrtpipe'
        self.creat_xml = '/bioinfo//Genomic/Sofware/smrttools/smrtcmds/bin/dataset'
        self.bax2bam = "/bioinfo/Genomic/Sofware/smrttools/smrtcmds/bin/bax2bam"
        self.fasta2reference = "/bioinfo/Genomic/Sofware/smrttools/smrtcmds/bin/fasta-to-reference"
        self.perl_script = self.config.PACKAGE_DIR + '/bacgenome/get_methy.pl '
        self.perl_path = "/program/perl-5.24.0/bin/perl"

    def run(self):
        super(MethylationTool, self).run()
        self.data_xml = os.path.join(self.work_dir, "data.xml")
        sample_dir = os.path.join(self.output_dir, self.option("sample"))
        if not os.path.exists(sample_dir):
            os.mkdir(sample_dir)
        ref_type = self.option("ref_type")
        #ref_type = self.option("ref_input").prop["ref_type"]
        if ref_type == "xml":
            self.ref_xml = self.option("ref_input")    ####.prop["path"]
        elif ref_type == "fasta":
            self.creat_fasta_to_ref(self.option("ref_input"))   ##.prop["path"])
        input_type = self.option("input_type")
        if input_type == "h5":
            self.h5_to_bam(self.option("input"))
            self.creat_xml_file(self.bam_dir + "sample.subreads.bam", self.data_xml)
        elif input_type == "bam":
            self.creat_xml_file("subread", self.option("input"), self.data_xml)
        self.methy_predict()
        self.get_motif_detail()
        self.set_output()
        self.end()

    def methy_predict(self):
        cmd = '{} pipeline-id --local-only pbsmrtpipe.pipelines.ds_modification_motif_analysis -e eid_subread:{} -e eid_ref_dataset:{} -o {}'.\
            format(self.pbsmrtpipe, self.data_xml, self.ref_xml, self.work_dir)
        command = self.add_command('predict', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("甲基化预测成功")
        else:
            self.set_error("甲基化预测失败")

    def get_motif_detail(self):
        motif_detail = self.output_dir + "/" + self.option("sample") + "/{}.motif_detail.xls".format(self.option("sample"))
        self.motif_gff = self.work_dir + "/tasks/motif_maker.tasks.reprocess-0/motifs.gff"
        if not os.path.exists(self.motif_gff):
            self.set_error("motif_gff文件生成失败".format(self.motif_gff))
        cmd = '{} {} {} {} {} '.format(self.perl_path, self.perl_script, self.motif_gff, self.option('sample'),motif_detail)
        command = self.add_command('get_motif_detail', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("get_motif_detail文件生成成功")
        else:
            self.set_error("get_motif_detail文件生成失败")


    def creat_xml_file(self, inputfile, xml_name):
        '''
        # 用于为输入文件构建smrttool使用的xml文件
        filetype: subread, reference
        '''
        cmd = "{} create --type SubreadSet --force {} {} ".format(self.creat_xml,xml_name,inputfile)
        command = self.add_command('creat_xml', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("creat xml file sucess!")
        else:
            self.set_error("creat xml file failed!")

    def h5_to_bam(self, inputfiles):
        self.bam_dir = self.work_dir + "/" + self.option("sample") + "/data/"
        if not os.path.exists(self.bam_dir):
            os.makedirs(self.bam_dir)
        #cmd = self.bax2bam + "  --output-xml -o {} ".format(bam_dir)
        cmd = self.bax2bam + " -o {} ".format(self.bam_dir + "sample")
        inputfiles = inputfiles.split(",")
        for each in inputfiles:
            cmd += " {} ".format(each)
        command = self.add_command('h5_to_bam', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("h5_to_bam sucess!")
        else:
            self.set_error("h5_to_bam failed!")

    def creat_fasta_to_ref(self, ref_fasta):
        self.logger.info("start creat fasta fai")
        out_name = "reference"
        if os.path.exists(self.work_dir + "/reference"):
            shutil.rmtree(self.work_dir + "/reference")
        cmd = self.fasta2reference + " --skip-ngmlr {} {} {}".format(ref_fasta, self.work_dir, out_name)
        command = self.add_command('creat_fasta_fai', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("creat_fasta_fai sucess!")
        else:
            self.set_error("creat_fasta_fai failed!")
        self.ref_xml = self.work_dir + "/" +out_name +"/referenceset.xml"

    def set_output(self):
        self.motif_csv = self.work_dir + "/tasks/motif_maker.tasks.find_motifs-0/motifs.csv"
        motif_csv = pd.read_csv(self.motif_csv)
        link_motif_csv = self.output_dir + "/" +self.option("sample") + "/{}.motifs.csv".format(self.option("sample"))
        if len(motif_csv) > 1:
            if os.path.exists(link_motif_csv):
                os.remove(link_motif_csv)
            if os.path.exists(self.motif_csv):
                self.logger.info(self.motif_csv)
                os.link(self.motif_csv, link_motif_csv)


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            #"id": "CreatTable" + str(random.randint(1, 10000)),
            "id": "methylation_ref",
            "type": "tool",
            "name": "bacgenome.methylation",
            "instant": True,
            "options": dict(
                input="/mnt/ilustre/users/isanger/sg-users/guhaidong/old/bacgenome/pbdata/m181029_145741_42228_c101465842550000001823290111181857_s1_p0.1.bax.h5,/mnt/ilustre/users/isanger/sg-users/guhaidong/old/bacgenome/pbdata/m181029_145741_42228_c101465842550000001823290111181857_s1_p0.2.bax.h5,/mnt/ilustre/users/isanger/sg-users/guhaidong/old/bacgenome/pbdata/m181029_145741_42228_c101465842550000001823290111181857_s1_p0.3.bax.h5",
                ref_input="/mnt/ilustre/users/sanger-dev/sg-users/yuanshaohua/measurement/smrt/fasta_file/reference.fasta",
                sample="BH",
                input_type="h5",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()