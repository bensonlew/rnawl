# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# last_modify:20201021

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import unittest
import gevent.subprocess as subprocess
from biocluster.config import Config


class SomaticPredictModule(Module):
    """
    bam文件分析前处理
    """
    def __init__(self, work_id):
        super(SomaticPredictModule, self).__init__(work_id)
        options = [
            {"name": "ref_fasta", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "combine_vcf", "type": "infile", "format": "ref_rna_v2.common"},  # 输入用于做前处理的bam文件
        ]
        self.add_option(options)
        self.dbnsfp_tools = self.add_tool("medical_transcriptome.somatic.dbnsfp_tools")
        self.predict_stat = self.add_tool("medical_transcriptome.somatic.predict_stat")

    def check_options(self):
        # if not self.option("ref_fasta").is_set:
        #     raise OptionError("缺少ref_fasta参数")
        # if not self.option("in_bam"):
        #     raise OptionError("请输入bam_文件")
        # if self.option("call_type"):
        #     if self.option("call_type") not in ["samtools", "gatk", "freebayes", "sentieon"]:
        #         raise OptionError("call_type方式不合法，必须是samtools, gatk, freebayes, sentieon", code="25600307")
        # else:
        #     raise OptionError("缺少call_type参数", code="25600308")
        # ref_file = os.path.dirname(self.option("ref_fasta").prop['path'])
        # ref_file_name = os.path.basename(self.option("ref_fasta").prop['path']).split(".fa")[0]



        #if not os.path.isfile(os.path.join(ref_file, "{}.dict".format(ref_file_name))) \
                #or not os.path.isfile(os.path.join(ref_file, "{}.fa.fai".format(ref_file_name))):
            #raise OptionError("参考组配置文件缺少dcit与fa.fai文件", code="24500406")
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run_vcf_split(self):
        self.vcf_split=self.add_tool("medical_transcriptome.somatic.vcf_split")
        self.vcf_split.set_options({
            "combine_vcf": self.option("combine_vcf").prop["path"],
            "line_num": 10000,
        })
        if self.option("analysis_format").lower() == "bam":
            self.bam_addrg.on('end', self.run_bam_sort_index)
        elif self.option("analysis_format").lower() == "cram":
            self.bam_addrg.on('end', self.run_bam2cram)
        self.bam_addrg.run()

    def run_somatic_predict(self):

        self.dbnsfp_tools.set_options({
            "combine_vcf": self.option("combine_vcf").prop["path"],
        })
        self.dbnsfp_tools.on("end", self.run_predict_stat)
        self.dbnsfp_tools.run()

    def run_predict_stat(self):
        self.predict_stat.set_options({
            "predict_result":os.path.join(self.dbnsfp_tools.output_dir,"result"),
        })
        self.predict_stat.on("end", self.set_output)
        self.predict_stat.run()


    def set_output(self):
        # allfiles = os.listdir(self.sentieon_haplptyper.output_dir)
        # resultfiles = [os.path.join(self.sentieon_haplptyper.output_dir, i) for i in allfiles]
        # outfiles = [os.path.join(self.output_dir, i) for i in allfiles]
        # for outfile in outfiles:
        #     if os.path.exists(outfile):
        #         if os.path.isfile(outfile):
        #             os.remove(outfile)
        #         else:
        #             os.system('rm -r %s' % outfile)
        #         # self.logger.info('rm -r %s' % newfile)
        # for i in range(len(allfiles)):
        #     os.link(resultfiles[i], outfiles[i])
        self.end()


    def run(self):
        super(SomaticPredictModule, self).run()
        self.logger.info("开始准备参考基因组索引文件及相关切割工作")
        # self.genome_split()
        self.logger.info("参考基因组索引文件准备及切割工作完成")
        self.run_somatic_predict()


    def set_bam_list(self):
        """
        samtools和freebayes与gatk的bamlist格式不一致，这里进行处理下
        """
        self.logger.info("开始设置bamlist！")
        with open(self.option("bam_list"), "r") as r, open("{}/bam.list".format(self.work_dir), "w") as w:
            data = r.readlines()
            for line in data:
                temp = line.strip().split("\t")
                w.write(temp[1] + "\n")
        self.logger.info("设置bamlist成功！")

    def end(self):
        super(SomaticPredictModule, self).end()



class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna'
        data = {
            "id": "snp" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "ref_rna_v3.pre4sentieon",
            "instant": False,
            "options": dict(
                ref_fasta="/mnt/ilustre/users/sanger-dev/workspace/20200114/Single_sentieonnew_majorbio_2319378889/CallSnpIndelSentieon/Pre4sentieon2/Lcu_new_hic.fasta",
                in_bam="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/piplinetest/test_data/test6/Refrna_majorbio_231937/output/bam/LC_PE_2.bam",
                call_type="sentieon",
                analysis_format="bam"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()