# -*- coding: utf-8 -*-

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil
import re
import subprocess
import itertools
import pandas as pd


class BugbasePredictAgent(Agent):
    """
    多样性进行BugBase预测
    现有开发思路是不允许老师对上传样本的标签进行分类分析
    """
    def __init__(self, parent):
        super(BugbasePredictAgent, self).__init__(parent)
        options = [
            {"name": "otu_fasta", "type": "infile", "format": "sequence.fasta"}, # otu代表序列
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table"}, # otu丰度表
            {"name": "mapping_file", "type": "infile", "format": "bacgenome.simple_file"},  # mapping_file
        ]

        self.add_option(options)
        self.step.add_steps("bugbase_predict")
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.bugbase_predict.start()
        self.step.update()

    def step_end(self):
        self.step.bugbase_predict.finish()
        self.step.update()

    def check_options(self):
        """
        进行参数二次检查
        :return:
        """
        if not self.option("otu_fasta").is_set:
            raise OptionError("请传入otu_reps序列！")
        if not self.option("otu_table").is_set:
            raise OptionError("请传入otu丰度文件！")

    def set_resource(self):
        """
        设置资源
        :return:
        """
        self._cpu = 4
        self._memory = '60G'

    def end(self):
        """
        结束啦
        :return:
        """
        super(BugbasePredictAgent, self).end()


class BugbasePredictTool(Tool):
    def __init__(self, config):
        super(BugbasePredictTool, self).__init__(config)
        """
        设置软件、脚本、数据库的路径
        """
        self.set_environ(PYTHONPATH = self.config.SOFTWARE_DIR + "/bioinfo/meta/picrust-1.1.0/")
        self.python_scripts = "/program/Python/bin/"
        self.scripts = self.config.SOFTWARE_DIR + "/bioinfo/meta/16s_scripts/"
        self.picrust_path = self.config.SOFTWARE_DIR + "/bioinfo/meta/picrust-1.1.0/scripts/"
        #self.Fundb = self.config.SOFTWARE_DIR + "/database/taxon_db/train_RDP_taxon/RDP_trained_greengenes135_16s_bacteria/RdpTaxonAssigner_training_seqs.fasta"  # 设置数据库文件路径，97_otus.fasta
        self.Fundb = self.config.SOFTWARE_DIR + "/database/align/ncbi/db/bugbase_greengene_database.fasta"  # 设置数据库文件路径，97_otus.fasta  /bioinfo/meta/16sFundb/rep_set/97_otus.fasta
        self.Fundb_tax = self.config.SOFTWARE_DIR + "/database/align/ncbi/db/bugbase_greengene_database.tax"
        self.bugbase = "{0}/bioinfo/meta/BugBase-master:{0}/program/R-3.3.1/bin:{0}/bioinfo/meta/BugBase-master/bin:{0}/program/R-3.3.1/bin:{0}/bioinfo/meta/BugBase-master/R_lib".format(self.config.SOFTWARE_DIR,self.config.SOFTWARE_DIR,self.config.SOFTWARE_DIR,self.config.SOFTWARE_DIR)
        self.r_home = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/"
        self.library = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/library: {0}/bioinfo/meta/BugBase-master/R_lib" .format(self.config.SOFTWARE_DIR)

        self.set_environ(PATH=self.bugbase, LD_LIBRARY_PATH=self.library, R_HOME=self.r_home)
        self.set_environ(BUGBASE_PATH="{0}/bioinfo/meta/BugBase-master".format(self.config.SOFTWARE_DIR))
        self.bugbase_scripts = "{0}/bioinfo/meta/BugBase-master/bin/".format(self.config.SOFTWARE_DIR)
        self.r_path = "/program/R-3.3.1/bin/Rscript"
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.db_path = self.config.SOFTWARE_DIR + "/database/align/ncbi/db/bugbase_greengene_database"
        self.cmd_path = "bioinfo/align/ncbi-blast-2.3.0+/bin/blastn"
        self.set_environ(BLASTDB=self.db_path)

    def run_blast(self):
        """
        直接与greengene 16s 细菌库比对得到greengene id
        :return:
        """
        cmd = "{} -query {} -db {} -out {} -evalue 1e-5 -num_threads 8 -max_target_seqs 1 -outfmt 6".format(self.cmd_path,
            self.option("otu_fasta").prop["path"], self.db_path, self.work_dir+"/blast.txt")
        self.logger.info("开始运行blast")
        blast_command = self.add_command("blast", cmd)
        blast_command.run()
        self.wait()
        if blast_command.return_code == 0:
            self.logger.info("运行blast完成")
        else:
            self.set_error("blast运行出错!", code="31100202")

    def blast_replace_otu_name(self):
        self.logger.info("对原始otu表替换otu的id")
        otu_table_path = self.option("otu_table").prop['path']
        ref_table_path = self.work_dir+"/blast.txt"
        # 获取对应otu和greengene的对应关系
        with open(ref_table_path, 'r') as m:
            lins = m.readlines()
            self.otu_greengene = {}
            self.greengene_otu = {}
            for lin in lins:
                lin = lin.strip().split("\t")
                if lin[0] not in self.otu_greengene.keys():
                    self.otu_greengene[lin[0]] = lin[1]
                if lin[1] not in self.greengene_otu.keys():
                    self.greengene_otu[lin[1]] = lin[0]

        new_otu_table_path = os.path.join(self.work_dir, 'new_otu_table.xls')
        if os.path.exists(new_otu_table_path):
            os.remove(new_otu_table_path)
        ww = open(new_otu_table_path, 'w')
        with open(otu_table_path, 'r') as f:
            lines = f.readlines()
            head_list = lines[0].strip().split("\t")
            head_list[0] = "GreenGene_ID"
            new_headline = '\t'.join(head_list)
            ww.write(new_headline + "\n")
            for line in lines[1:]:
                if not line.startswith("#"):
                    data_list = []
                    line = line.strip().split("\t")
                    otu_id = line[0]
                    if otu_id in self.otu_greengene.keys():
                        data_list.append(self.otu_greengene[otu_id])
                        if len(data_list) != 0:
                            data_list = data_list + line[1:]
                            ww.write('\t'.join(data_list))
                            ww.write('\n')
                        else:
                            self.set_error("未能成功获取正确的otu与greengene的对应关系")
                else:
                    ww.write(line)
        ww.close()
        # 对相同greengenes数据库中相同id的数据进行合并
        data = pd.read_table(new_otu_table_path, sep='\t', header=0)
        data = data.groupby("GreenGene_ID").sum()
        bugbase_table_path = os.path.join(self.work_dir, 'otu_input_table.txt')
        data.to_csv(bugbase_table_path, sep='\t')
        self.logger.info("替换otu文件完成")



    def picrust_predict(self):
        """
        根据代表序列获取otu与greengene数据库中的对应关系
        通过比对GreenGene--16s数据库，得到97%相似性的物种或者otu对应关系
        """
        self.logger.info("开始用picrust获取greengene数据库中的对应关系")
        if os.path.exists(os.path.join(self.work_dir, "Results")):
            shutil.rmtree(os.path.join(self.work_dir, "Results"))
        else:
            os.mkdir(os.path.join(self.work_dir, "Results"))
        result_path = os.path.join(self.work_dir, "Results")
        self.logger.info("设置picrust的结果文件完成")

        if not os.path.exists(os.path.join(self.work_dir,"otu_picking_params_97.txt")):
            self.logger.info("判断txt文件是否存在")
            f = open("otu_picking_params_97.txt", "wb")
            f.write("pick_otus:enable_rev_strand_match True\n")
            f.write("pick_otus:similarity 0.97")
            f.close()
        otu_params_path = os.path.join(self.work_dir, "otu_picking_params_97.txt")
        cmd = "{}pick_closed_reference_otus.py -i {} -o {} -f -p {} -r {} -t {}".format(self.python_scripts, self.option("otu_fasta").prop["path"], result_path, otu_params_path, self.Fundb, self.Fundb_tax)
        self.logger.info(cmd)
        command = self.add_command("picrust_predict", cmd).run()
        self.wait(command)
        if command.return_code in [0]:
            self.logger.info("picrust运行完成！")
        else:
            self.set_error("picrust预测失败！")

    def replace_otu_name(self):
        """
        根据picrust预测得到的otu和序列的对应关系替换otu表的otu的id，并对原始otu表进行整理合并
        为什么要这样做？？是因为bugbase中只识别greengene数据库中的id，对otu为序列号的表并不识别
        """
        self.logger.info("对原始otu表替换otu的id")
        otu_table_path = self.option("otu_table").prop['path']
        ref_table_path = os.path.join(self.work_dir, "Results", "uclust_ref_picked_otus/otu_fasta_rep_otus.txt")
        #获取对应otu和greengene的对应关系
        with open(ref_table_path, 'r') as m:
            lins = m.readlines()
            self.otu_greengene = {}
            self.greengene_otu = {}
            for lin in  lins:
                lin = lin.strip().split("\t")
                if lin[1] not in self.otu_greengene.keys():
                    self.otu_greengene[lin[1]] = lin[0]
                if lin[0] not in self.greengene_otu.keys():
                    self.greengene_otu[lin[0]] = lin[1]

        new_otu_table_path = os.path.join(self.work_dir, 'new_otu_table.xls')
        if os.path.exists(new_otu_table_path):
            os.remove(new_otu_table_path)
        ww = open(new_otu_table_path, 'w')
        with open(otu_table_path, 'r') as f:
            lines = f.readlines()
            head_list = lines[0].strip().split("\t")
            head_list[0] = "GreenGene_ID"
            new_headline = '\t'.join(head_list)
            ww.write(new_headline + "\n")
            for line in lines[1:]:
                if not line.startswith("#"):
                    data_list = []
                    line = line.strip().split("\t")
                    otu_id = line[0]
                    if otu_id in self.otu_greengene.keys():
                        data_list.append(self.otu_greengene[otu_id])
                        if len(data_list) != 0:
                            data_list = data_list + line[1:]
                            ww.write('\t'.join(data_list))
                            ww.write('\n')
                        else:
                            self.set_error("未能成功获取正确的otu与greengene的对应关系")
                else:
                    ww.write(line)
        ww.close()
        # 对相同greengenes数据库中相同id的数据进行合并
        data = pd.read_table(new_otu_table_path, sep='\t', header=0)
        data = data.groupby("GreenGene_ID").sum()
        bugbase_table_path = os.path.join(self.work_dir, 'otu_input_table.txt')
        data.to_csv(bugbase_table_path, sep='\t')
        self.logger.info("替换otu文件完成")

    def run_bugbase_predict(self):
        """
        运行bugbase软件进行预测结果，没有kegg结果
        """
        self.logger.info("开始运行bugbase软件进行预测")
        bugbase_out = os.path.join(self.work_dir, "bugbase")
        if os.path.exists(bugbase_out):
            shutil.rmtree(bugbase_out)
        input_bugbase = os.path.join(self.work_dir, 'otu_input_table.txt')
        cmd = "{} {}run.bugbase.r -i {}  -m {} -c BODY_SITE -o {}".format(self.r_path, self.bugbase_scripts, input_bugbase, self.option("mapping_file").prop["path"], bugbase_out)

        self.logger.info(cmd)
        command = self.add_command("bugbase_predict", cmd).run()
        self.wait(command)
        if command.return_code in [0]:
            self.logger.info("bugbase运行完成！")
        else:
            self.set_error("与bugbase数据库比对成功的greengene id过少！")

    def run_kegg_predict(self):
        """
        运行Bugbase软件对kegg结果进行预测，
        :return:
        """
        self.logger.info("开始运行bugbase软件对kegg通路进行预测")
        bugbase_out = os.path.join(self.work_dir, "bugbase_kegg")
        if os.path.exists(bugbase_out):
            shutil.rmtree(bugbase_out)
        input_bugbase = os.path.join(self.work_dir, 'otu_input_table.txt')
        cmd = "{} {}/run.bugbase.r -i {} -a -x -k -o {}".format(self.r_path, self.bugbase_scripts, input_bugbase, bugbase_out)

        self.logger.info(cmd)
        command = self.add_command("bugbase_predict_kegg", cmd).run()
        self.wait(command)
        if command.return_code in [0]:
            self.logger.info("bugbase预测kegg通路完成！")
        else:
            self.set_error("bugbase预测kegg通路失败！")

    def translate_table(self, table, out_table):
        """
        转置一个表
        :param table:输入表path
        :param out_table: 输出表path
        :return:
        """
        self.logger.info("开始对table表进行转置和修改结果表！")
        data = pd.read_table(table, sep='\t', header=0)
        name_list = data.columns
        name_list2 = ['Phenotypes']
        for i in name_list[1:]:
            name_list2.append(i)
        data.columns = name_list2
        data2 = data.astype('str').T
        data2.to_csv(out_table, sep='\t', header=0)

    def translate_kegg_table(self, table, out_table):
        """
        对kegg得到的结果进行处理，转置和切分为两列
        :param table: 输入表
        :param out_table: 输出表
        :return:
        """
        self.logger.info("对kegg预测结果表进行转置和修改结果表！")
        data = pd.read_table(table, sep='\t', header=0)
        name_list = data.columns
        name_list[0] = 'Phenotypes'
        data.columns = name_list
        data2 = (data.set_index('Phenotypes')).T
        all_name_list = list(data2.columns)
        data2['Module'], data2['Description'] = data2['index'].str.split("_", 1).str
        data2['Description'] = data2['Description'].str.rstrip("_")
        all_name_list.remove("index")
        all_name_list.remove("Module")
        all_name_list.remove("Description")
        select_list = ['Module', 'Description'] + all_name_list
        data3 = data2[select_list]
        data3.to_csv(out_table, sep='\t', index=0)
        self.logger.info("对kegg的统计结果完成")

    def add_otu_id(self, table, out_table):
        """
        对bugbase预测得到的结果进行增加替换的名称
        主要是16s的otu注释的结果表
        :return:
        """
        self.logger.info("开始对预测的结果增加替换结果的名称")
        with open(table, 'r') as f, open(out_table, 'w') as w:
            lines = f.readlines()
            head_list = lines[0].strip().split("\t")
            new_head_list = ['OTU_ID', 'GreenGene_ID'] + head_list
            w.write("{}\n".format("\t".join(new_head_list)))
            for line in lines[1:]:
                line = line.strip().split("\t")
                greengene_id = line[0]
                if greengene_id in self.greengene_otu.keys():
                    otu_id = self.greengene_otu[greengene_id]
                    new_line_list = [otu_id] + line
                    new_line_string = '\t'.join(new_line_list)
                    w.write(new_line_string + "\n")
        self.logger.info("替换16s预测结果完成！")

    def run(self):
        """
        运行
        """
        super(BugbasePredictTool, self).run()
        self.logger.info("开始运行预测啦")
        self.picrust_predict()
        self.replace_otu_name()
        self.run_bugbase_predict()
        #self.run_kegg_predict()
        self.set_output()
        self.logger.info("运行tool结束")
        self.end()

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始设置结果文件目录")
        self.logger.info("开始预测bugbase的文件连接")
        bugbase_path = os.path.join(self.work_dir, "bugbase")
        normalize_table_path = os.path.join(bugbase_path, "normalized_otus/16s_normalized_otus.txt")
        out_normalize_table_path = os.path.join(self.output_dir, "16s_normalized_asv.txt")
        if os.path.exists(out_normalize_table_path):
            os.remove(out_normalize_table_path)
            os.link(normalize_table_path, out_normalize_table_path)
        else:
            os.link(normalize_table_path, out_normalize_table_path)
        new_out_normalize_table_path = os.path.join(self.output_dir, "16s_normalized_asv.txt_new")
        self.add_otu_id(out_normalize_table_path, new_out_normalize_table_path)
        os.remove(out_normalize_table_path)
        os.rename(new_out_normalize_table_path, out_normalize_table_path)

        bugbase_prediction = os.path.join(bugbase_path, "predicted_phenotypes/predictions.txt")
        out_bugbase_prediction = os.path.join(self.output_dir, "bugbase_predictions.txt")
        if os.path.exists(out_bugbase_prediction):
            os.remove(out_bugbase_prediction)
            os.link(bugbase_prediction, out_bugbase_prediction)
        else:
            os.link(bugbase_prediction, out_bugbase_prediction)
        new_out_bugbase_prediction = os.path.join(self.output_dir, "bugbase_predictions.txt_new")
        self.translate_table(out_bugbase_prediction, new_out_bugbase_prediction)
        os.remove(out_bugbase_prediction)
        os.rename(new_out_bugbase_prediction, out_bugbase_prediction)
        os.system("sed -i 's/NA/0/g' {}".format(out_bugbase_prediction))

        bugbase_thresholds = os.path.join(bugbase_path, "thresholds/thresholds_used.txt")
        out_bugbase_thresholds = os.path.join(self.output_dir, "bugbase_thresholds.txt")
        if os.path.exists(out_bugbase_thresholds):
            os.remove(out_bugbase_thresholds)
            os.link(bugbase_thresholds, out_bugbase_thresholds)
        else:
            os.link(bugbase_thresholds, out_bugbase_thresholds)
        bugbase_variances = os.path.join(bugbase_path, "thresholds/variances.txt")
        out_bugbase_variances = os.path.join(self.output_dir, "bugbase_variances.txt")
        if os.path.exists(out_bugbase_variances):
            os.remove(out_bugbase_variances)
            os.link(bugbase_variances, out_bugbase_variances)
        else:
            os.link(bugbase_variances, out_bugbase_variances)

        contributing_file = os.path.join(bugbase_path, "otu_contributions/contributing_otus.txt")
        contributing_file_new = os.path.join(self.output_dir, "contributing_asv.txt")
        if os.path.exists(contributing_file_new):
            os.remove(contributing_file_new)
            os.link(contributing_file, contributing_file_new)
        else:
            os.link(contributing_file, contributing_file_new)
        contributing_file_new1 = os.path.join(self.output_dir, "contributing_asv.txt_new")
        self.add_otu_id(contributing_file_new, contributing_file_new1)
        os.remove(contributing_file_new)
        os.rename(contributing_file_new1, contributing_file_new)

        self.logger.info("完成预测bugbase的文件连接")

        """
        self.logger.info("开始预测bugbase的kegg文件连接")
        bugbase_kegg_path = os.path.join(self.work_dir, "bugbase_kegg")
        normalize_table_kegg_path = os.path.join(bugbase_kegg_path, "normalized_otus/16s_normalized_otus.txt")
        out_normalize_table_kegg_path = os.path.join(self.output_dir, "16s_normalized_otus_kegg.txt")
        if os.path.exists(out_normalize_table_kegg_path):
            os.remove(out_normalize_table_kegg_path)
            os.link(normalize_table_kegg_path, out_normalize_table_kegg_path)
        else:
            os.link(normalize_table_kegg_path, out_normalize_table_kegg_path)
        new_out_normalize_table_kegg_path = os.path.join(self.output_dir, "16s_normalized_otus_kegg.txt_new")
        self.add_otu_id(out_normalize_table_kegg_path, new_out_normalize_table_kegg_path)
        os.remove(out_normalize_table_kegg_path)
        os.rename(new_out_normalize_table_kegg_path, out_normalize_table_kegg_path)

        bugbase_kegg_prediction = os.path.join(bugbase_kegg_path, "predicted_phenotypes/predictions.txt")
        out_bugbase_kegg_prediction = os.path.join(self.output_dir, "bugbase_kegg_predictions.txt")
        if os.path.exists(out_bugbase_kegg_prediction):
            os.remove(out_bugbase_kegg_prediction)
            os.link(bugbase_kegg_prediction, out_bugbase_kegg_prediction)
        else:
            os.link(bugbase_kegg_prediction, out_bugbase_kegg_prediction)
        new_out_bugbase_kegg_prediction = os.path.join(self.output_dir, "bugbase_kegg_predictions.txt_new")
        self.translate_kegg_table(out_bugbase_kegg_prediction, new_out_bugbase_kegg_prediction)
        os.remove(out_bugbase_kegg_prediction)
        os.rename(new_out_bugbase_kegg_prediction, out_bugbase_kegg_prediction)


        bugbase_kegg_thresholds = os.path.join(bugbase_kegg_path, "thresholds/thresholds_used.txt")
        out_bugbase_kegg_thresholds = os.path.join(self.output_dir, "bugbase_kegg_thresholds.txt")
        if os.path.exists(out_bugbase_kegg_thresholds):
            os.remove(out_bugbase_kegg_thresholds)
            os.link(bugbase_kegg_thresholds, out_bugbase_kegg_thresholds)
        else:
            os.link(bugbase_kegg_thresholds, out_bugbase_kegg_thresholds)

        bugbase_kegg_variances = os.path.join(bugbase_kegg_path, "thresholds/variances.txt")
        out_bugbase_kegg_variances = os.path.join(self.output_dir, "bugbase_kegg_variances.txt")
        if os.path.exists(out_bugbase_kegg_variances):
            os.remove(out_bugbase_kegg_variances)
            os.link(bugbase_kegg_variances, out_bugbase_kegg_variances)
        else:
            os.link(bugbase_kegg_variances, out_bugbase_kegg_variances)

        self.logger.info("完成预测bugbase的kegg文件连接")
        """
