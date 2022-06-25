# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# modified 2020.10.26

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
import os
import subprocess
import pandas as pd 
from biocluster.config import Config

class PredictStatAgent(Agent):
    """
    基因比对接口
    """

    def __init__(self, parent):
        super(PredictStatAgent, self).__init__(parent)
        options = [
            {"name": "predict_result", "type": "infile", "format": "ref_rna_v2.common"},  # 输入的bam
        ]
        self.add_option(options)
        self._memory_increase_step = 200

    def check_options(self):
        # if not self.option("bam_file"):
        #     raise OptionError("请设置bam路径")
        # if not self.option("fa_file"):
        #     raise OptionError("请设置ref.fa路径")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        super(PredictStatAgent, self).end()


class PredictStatTool(Tool):
    def __init__(self, config):
        super(PredictStatTool, self).__init__(config)
        # self.R = self.config.SOFTWARE_DIR + "/bioinfo/miniconda2/bin/Rscript"
        self.R = self.config.SOFTWARE_DIR + "/bioinfo/rna/miniconda2/bin/Rscript"
        self.Rscrip = self.config.PACKAGE_DIR + '/medical_transcriptome/predict_scatter.r'
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)

    def run_predict_stat(self):
        """
        对somatic预测结果进行统计

        """
        predict_detail_df = pd.read_table(self.option("predict_result").prop["path"],usecols=['#CHROM', 'POS', 'REF','ALT','ref','alt','SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred','Ensembl_geneid','genename'])
        predict_detail_df["SIFT_score"] = predict_detail_df["SIFT_score"].apply(lambda x: x.split(";")[0])
        predict_detail_df["SIFT_pred"] = predict_detail_df["SIFT_pred"].apply(lambda x: x.split(";")[0])
        predict_detail_df["Polyphen2_HDIV_score"] = predict_detail_df["Polyphen2_HDIV_score"].apply(lambda x: x.split(";")[0])
        predict_detail_df["Polyphen2_HDIV_pred"] = predict_detail_df["Polyphen2_HDIV_pred"].apply(lambda x: x.split(";")[0])
        stat_df = pd.DataFrame()
        stat_df["TYPE"] = ["SIFT_pred", "Polyphen2_HDIV_pred", "common"]
        SIFT_num = predict_detail_df.shape[0] - predict_detail_df.loc[:, "SIFT_score"].value_counts()["."]
        Polyphen2_HDIV_num = predict_detail_df.shape[0] - predict_detail_df.loc[:, "Polyphen2_HDIV_score"].value_counts()["."]
        Common_num = predict_detail_df[(predict_detail_df["Polyphen2_HDIV_score"] != ".") & (predict_detail_df["SIFT_score"] != ".") ].shape[0]
        stat_df["Number"] = [SIFT_num,Polyphen2_HDIV_num,Common_num]
        SIFT_score = predict_detail_df.loc[:, "SIFT_pred"].value_counts()["D"]
        Polyphen2_HDIV_score = predict_detail_df.loc[:, "Polyphen2_HDIV_pred"].value_counts()["D"]+ predict_detail_df.loc[:, "Polyphen2_HDIV_pred"].value_counts()["P"]
        Common_score = len ( set(predict_detail_df[predict_detail_df["SIFT_pred"] == "D"].index) |set(predict_detail_df[predict_detail_df["Polyphen2_HDIV_pred"] == "D"].index)|set(predict_detail_df[predict_detail_df["Polyphen2_HDIV_pred"] == "P"].index))
        stat_df["Score"] = [SIFT_score, Polyphen2_HDIV_score, Common_score]
        stat_df.to_csv(os.path.join(self.output_dir,"predict_stat"),sep="\t",index=False)
        predict_detail_df.to_csv(os.path.join(self.output_dir,"predict_detail"),sep="\t",index=False)

    def run_plot_predeal(self):
        """
        对s统计结果进行过滤和我预处理准备画图

        """
        detail_df  = pd.read_table(os.path.join(self.output_dir,"predict_detail"))
        draw_info_df = detail_df[[ "Polyphen2_HDIV_score","SIFT_score"]]
        draw_info_df_filter = draw_info_df[
            (draw_info_df["SIFT_score"] != ".") & (draw_info_df["Polyphen2_HDIV_score"] != ".")]
        def type(series):
            if series["SIFT_score"] < 0.05:
                if series["Polyphen2_HDIV_score"] > 0.453:
                    type = "SIFT and Polyphen"
                else:
                    type = "SIFT"
            else:
                if series["Polyphen2_HDIV_score"] > 0.453:
                    type = "Polyphen2_HDIV"
                else:
                    type = "other"
            return type
        draw_info_df_filter["SIFT_score"] = draw_info_df_filter["SIFT_score"].astype("float")
        draw_info_df_filter["Polyphen2_HDIV_score"] = draw_info_df_filter["Polyphen2_HDIV_score"].astype("float")
        draw_info_df_filter["type"] = draw_info_df_filter.apply(type, axis=1)
        draw_info_df_filter.to_csv(os.path.join(self.work_dir,"predict_matrix"), sep="\t", index=False)

    def run_scatter_plot(self):
        '''绘制散点图'''
        self.logger.info("开始绘制散点图")
        cmd = "{} ".format(self.R)
        cmd += "{} ".format(self.Rscrip)
        cmd += "-i {} ".format(os.path.join(self.work_dir,"predict_matrix"))
        cmd += "-o {} ".format(os.path.join(self.output_dir, "predict"))
        self.logger.info("使用R进行散点图绘制,代码如下")
        self.logger.info(cmd)
        command = self.add_command("scatter_plot", cmd, ignore_error=True, shell=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("散点图绘制完成!")
        elif command.return_code in [1, -9]:  # 当返回码为1或-9，加内存重试
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("散点图绘制失败!")

    def run(self):
        super(PredictStatTool, self).run()
        self.run_predict_stat()
        self.run_plot_predeal()
        self.run_scatter_plot()
        self.end()



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
            "id": "dbnsfp" + str(random.randint(1, 10000))+"yyyy",
            "type": "tool",
            "name": "medical_transcriptome.somatic.predict_stat",
            "instant": False,
            "options": dict(
                predict_result="/mnt/ilustre/users/sanger-dev/workspace/20201123/MedicalTranscriptome_e4dcldqfiur7iecm4s2ok0apv9/CallSomaticSentieon/SomaticPredictNew/DbnsfpToolsNew/output/final_result",
                # fa_file="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/dna/Homo_sapiens.GRCh38.dna.toplevel.fa",
                # name="add_sort",
                # file_format="bam",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()