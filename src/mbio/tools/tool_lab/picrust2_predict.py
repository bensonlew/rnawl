# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'@20191224


from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil
import re
import subprocess
import itertools
import pandas as pd


class Picrust2PredictAgent(Agent):
    """
    picrust2预测（主要功能是预测16s、真菌和ITS序列）
    这里因为picrust2的计算结果真菌和ITS与16s的结果有差异，只有酶的结果和pathway的结果
    """
    def __init__(self, parent):
        super(Picrust2PredictAgent, self).__init__(parent)
        options = [
            {"name": "otu_fasta", "type": "infile", "format": "sequence.fasta"}, ##代表序列OTU或者ASV序列
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table"},##代表序列的otu表
            {"name": "thread", "type": "string", "default": "4"}, ##线程数
            {"name": "database", "type": "string", "default": "COG,EC,KO"},##注释参考库
            {"name": "hsp", "type": "string", "default": "mp"}, ##HSP方法选择，主要有5中方法{mp, emp_prob, subtree_average,pic,scp}
            {"name": "analysis_type", "type": "string", "default": "16S"},##选择方法16S、18S、ITS
        ]
        self.add_option(options)
        self.step.add_steps("picrust_predict")
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self._memory_increase_step = 20

    def step_start(self):
        self.step.picrust_predict.start()
        self.step.update()

    def step_end(self):
        self.step.picrust_predict.finish()
        self.step.update()

    def check_options(self):
        """
        进行参数二次检查
        :return:
        """
        if not self.option("otu_fasta").is_set:
            raise OptionError("请传入otu_reps序列！", code="32707401")
        if not self.option("otu_table").is_set:
            raise OptionError("请传入otu_table.xls文件！", code="32707402")
        if (not self.option("database")) and (self.option("analysis_type") not in ["16S"]):
            raise OptionError("必须提供参考的数据库database！", code="32707403")
        if self.option("hsp") not in ["mp", "emp_prob", "pic", "scp", "subtree_average"]:
            raise OptionError("请传入HSP方法", code="32707404")

    def set_resource(self):
        """
        设置初始化资源的利用
        :return:
        """
        otu_path = self.option("otu_fasta").prop["path"]
        number = 0
        with open(otu_path, 'r') as f:
            for line in f:
                if line[0] == ">":
                    number += 1
        if number >= 1000:
            self._cpu = 8
            self._memory = '50G'
        else:
            self._cpu = 8
            self._memory = '40G'

    def end(self):
        self.logger.info("运行结束啦")
        super(Picrust2PredictAgent, self).end()


class Picrust2PredictTool(Tool):
    def __init__(self, config):
        super(Picrust2PredictTool, self).__init__(config)
        self.activate = self.config.SOFTWARE_DIR + "/program/miniconda3/bin/activate"
        self.picrust = self.config.SOFTWARE_DIR + "/program/miniconda3/envs/picrust2" ##设置picrust2的conda环境
        self.shell_path = "/program/sh"
        self.picrust2_script = self.config.SOFTWARE_DIR + "/bioinfo/meta/picrust2/picrust2-2.2.0-b/scripts/picrust2_pipeline.py"
        self.deactivate = self.config.SOFTWARE_DIR + "/program/miniconda3/bin/deactivate"
        self.pick_annotation = self.config.PACKAGE_DIR + "/meta/scripts/pick_annotation.py"
        self.python = "/miniconda2/bin/python"
        self.ref_dir_18s = self.config.SOFTWARE_DIR + "/bioinfo/meta/picrust2/picrust2-2.2.0-b/picrust2/default_files/fungi/fungi_18S"
        self.custom_trait_tables_18s = self.config.SOFTWARE_DIR + "/bioinfo/meta/picrust2/picrust2-2.2.0-b/picrust2/default_files/fungi/ec_18S_counts.txt.gz"
        self.marker_gene_table_18s = self.config.SOFTWARE_DIR + "/bioinfo/meta/picrust2/picrust2-2.2.0-b/picrust2/default_files/fungi/18S_counts.txt.gz"
        self.pathway_map_18s =  self.config.SOFTWARE_DIR + "/bioinfo/meta/picrust2/picrust2-2.2.0-b/picrust2/default_files/pathway_mapfiles/metacyc_path2rxn_struc_filt_euk.txt"
        self.ref_dir_its = self.config.SOFTWARE_DIR + "/bioinfo/meta/picrust2/picrust2-2.2.0-b/picrust2/default_files/fungi/fungi_ITS"
        self.custom_trait_tables_its = self.config.SOFTWARE_DIR + "/bioinfo/meta/picrust2/picrust2-2.2.0-b/picrust2/default_files/fungi/ec_ITS_counts.txt.gz"
        self.marker_gene_table_its = self.config.SOFTWARE_DIR + "/bioinfo/meta/picrust2/picrust2-2.2.0-b/picrust2/default_files/fungi/ITS_counts.txt.gz"
        self.pathway_map_its =  self.config.SOFTWARE_DIR + "/bioinfo/meta/picrust2/picrust2-2.2.0-b/picrust2/default_files/pathway_mapfiles/metacyc_path2rxn_struc_filt_fungi.txt"
        # self.set_environ(PYTHONPATH = self.config.SOFTWARE_DIR + "/bioinfo/meta/picrust-1.1.0/")
        self.unzip_file_dict = {
            "pred_metagenome_unstrat.tsv.gz": "pred_metagenome_unstrat.tsv",
            "weighted_nsti.tsv.gz": "weighted_nsti.tsv",
            "COG_predicted.tsv.gz": "COG_predicted.tsv",
            "EC_predicted.tsv.gz": "EC_predicted.tsv",
            "KO_predicted.tsv.gz": "KO_predicted.tsv",
            "marker_predicted_and_nsti.tsv.gz": "marker_predicted_and_nsti.tsv",
            "path_abun_unstrat.tsv.gz": "path_abun_unstrat.tsv",
            "ec_18S_counts.txt_predicted.tsv.gz": "ec_18S_counts.txt_predicted.tsv",
            "ec_ITS_counts.txt_predicted.tsv.gz": "ec_ITS_counts.txt_predicted.tsv",
            "seqtab_norm.tsv.gz": "seqtab_norm.tsv"
        }

    def func_predict_analysis(self):
        """
        进行picrust2预测
        :return:
        """
        self.logger.info("开始运行picrust2软件进行预测啦！")
        create_cmd = ''
        if self.option("analysis_type") in ["16s", "16S"]:
            shell_cmd = self.create_cmd("16S")
            create_cmd = "{} {}".format(self.shell_path, shell_cmd)
        elif self.option("analysis_type") in ["18S", "18s"]:
            shell_cmd = self.create_cmd("18S")
            create_cmd = "{} {}".format(self.shell_path, shell_cmd)
        elif self.option("analysis_type") in ["ITS", "its", "Its"]:
            shell_cmd = self.create_cmd("ITS")
            create_cmd = "{} {}".format(self.shell_path, shell_cmd)
        self.logger.info(create_cmd)
        command = self.add_command("predict", create_cmd).run()
        self.wait(command)
        if command.return_code in [0]:
            self.logger.info("picrust2运行完成！")
        else:
            self.set_error("picrust2预测失败！", code="32707401")

    def create_cmd(self, type):
        """
        用于生成picrust2的运行shell命令
        :param type:
        :return:
        """
        shell_cmd = os.path.join(self.work_dir, "picrust2.sh")
        fasta_path = self.option("otu_fasta").prop['path']
        fasta_path = self.trim_N(fasta_path)  # by xieshichang 20200426
        otu_table_path = self.option("otu_table").prop['path']
        self.picrust2_out = os.path.join(self.work_dir, "picrust2")
        if os.path.exists(self.picrust2_out):
            shutil.rmtree(self.picrust2_out)
        f = open(shell_cmd, "w")
        if type in ["16s", "16S"]:
            f.write("#! /bin/bash\n")
            f.write("source {} {}\n".format(self.activate, self.picrust))
            f.write("python {} -s {} -i {} -o {} -p {} --in_traits {} --verbose\n".format(self.picrust2_script, fasta_path, otu_table_path,self.picrust2_out, 8 , self.option("database")))
            f.write("source {} {}\n".format(self.deactivate, self.picrust))
        elif type in ["18s", "18S"]:
            f.write("#! /bin/bash\n")
            f.write("source {} {}\n".format(self.activate, self.picrust))
            f.write("python {} -s {} -i {} -o {} --ref_dir {} --custom_trait_tables {} --marker_gene_table {} --reaction_func {} --pathway_map {} -p {} --verbose\n".format(self.picrust2_script, fasta_path, otu_table_path, self.picrust2_out, self.ref_dir_18s, self.custom_trait_tables_18s, self.marker_gene_table_18s, self.custom_trait_tables_18s, self.pathway_map_18s, 8))
            f.write("source {} {}\n".format(self.deactivate, self.picrust))
        elif type in ["ITS", "its", "Its"]:
            f.write("#! /bin/bash\n")
            f.write("source {} {}\n".format(self.activate, self.picrust))
            f.write("python {} -s {} -i {} -o {} --ref_dir {} --custom_trait_tables {} --marker_gene_table {} --reaction_func {} --pathway_map {} -p {} --verbose\n".format(self.picrust2_script, fasta_path, otu_table_path, self.picrust2_out, self.ref_dir_its, self.custom_trait_tables_its, self.marker_gene_table_its, self.custom_trait_tables_its, self.pathway_map_its, 8))
            f.write("source {} {}\n".format(self.deactivate, self.picrust))
        else:
            self.set_error("不存在的type类型%s", variables=(type), code="32707402")
        f.close()
        return shell_cmd

    def trim_N(self, path):
        """
        去除otu fasta文件中序列末尾的N
        by xieshichang 20200426
        """
        rt_path = self.work_dir + '/trim_N.' + os.path.basename(path)
        cmd = self.config.SOFTWARE_DIR + "/program/perl-5.24.0/bin/perl -lne 's/N+$// unless /^\>/;print' " + "{} > {}".format(path, rt_path)
        command = self.add_command('fasta_trim_n', cmd, shell=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("成功去除otu序列结尾的N")
        return rt_path

    def unzip_file(self):
        """
        对picrust2的预测结果进行解压
        :return:
        """
        self.logger.info("开始对所有文件结果进行解压！")
        all_files = os.listdir(self.picrust2_out)
        for file in all_files:
            file_path = os.path.join(self.picrust2_out, file)
            if os.path.isfile(file_path):
                if file.endswith(".gz"):
                    file_name = os.path.basename(file_path)
                    out_file_name = self.unzip_file_dict[file_name]
                    out_file_path = os.path.join(self.picrust2_out, out_file_name)
                    try:
                        self.logger.info("开始对文件%s进行解压！"%file_name)
                        cmd = "zcat {} > {}".format(file_path, out_file_path)
                        self.logger.info(cmd)
                        subprocess.check_output(cmd, shell=True)
                        self.logger.info("运行%s解压成功！"%file_name)
                    except subprocess.CalledProcessError:
                        self.set_error("对文件%s解压失败！", variables=(file_name), code="32707403")
            elif os.path.isdir(file_path):
                second_files = os.listdir(file_path)
                for f in second_files:
                    f_path = os.path.join(file_path, f)
                    self.logger.info("要对文件夹下的文件进行解压了")
                    if os.path.isfile(f_path):
                        if f.endswith(".gz"):
                            f_name = os.path.basename(f_path)
                            out_file_name = self.unzip_file_dict[f_name]
                            out_file_path = os.path.join(file_path, out_file_name)
                            try:
                                self.logger.info("开始对文件%s进行解压！"%f_name)
                                cmd = "zcat {} > {}".format(f_path, out_file_path)
                                self.logger.info(cmd)
                                subprocess.check_output(cmd, shell=True)
                                self.logger.info("运行%s解压成功！"%f_name)
                            except subprocess.CalledProcessError:
                                self.set_error("对文件%s解压失败！", variables=(f_name), code="32707404")
                    elif os.path.isdir(f_path):
                        self.logger.info("不再对文件夹下的文件夹进行解压了")
            else:
                self.set_error("不存在的文件类型：%s" , variables=(file_path), code="32707405")
        self.logger.info("运行文件解压成功！")

    def run_cog_anno(self):
        """
        进行cog功能的注释
        :return:
        """
        self.logger.info("开始对cog功能进行注释啦")
        out_cog = os.path.join(self.work_dir, "COG")
        if os.path.exists(out_cog):
            shutil.rmtree(out_cog)
        os.mkdir(out_cog)
        cog_abundance_path = os.path.join(self.picrust2_out, "COG_metagenome_out/pred_metagenome_unstrat.tsv")
        nsti_path = os.path.join(self.picrust2_out, "COG_metagenome_out/weighted_nsti.tsv")
        out_nsti_path = os.path.join(out_cog, "weighted_nsti.xls")
        os.link(nsti_path, out_nsti_path)
        if os.path.exists(cog_abundance_path):
            cog_cmd = "{} {} -i {} -database {} -o {}".format(self.python, self.pick_annotation, cog_abundance_path, "cog", out_cog)
            self.logger.info(cog_cmd)
            command = self.add_command("anno_cog", cog_cmd).run()
            self.wait(command)
            if command.return_code in [0]:
                self.logger.info("anno cog 运行完成！")
            else:
                self.set_error("anno cog 运行失败！", code="32707406")
        else:
            self.logger.info("cog 预测失败")

    def run_kegg_anno(self):
        """
        进行kegg功能的注释
        :return:
        """
        self.logger.info("开始对kegg功能进行注释啦")
        out_kegg = os.path.join(self.work_dir, "KEGG")
        if os.path.exists(out_kegg):
            shutil.rmtree(out_kegg)
        os.mkdir(out_kegg)
        ko_abundance_path = os.path.join(self.picrust2_out, "KO_metagenome_out/pred_metagenome_unstrat.tsv")
        EC_metagenome_out = os.path.join(self.picrust2_out, "EC_metagenome_out/pred_metagenome_unstrat.tsv")
        nsti_path = os.path.join(self.picrust2_out, "KO_metagenome_out/weighted_nsti.tsv")
        out_ko_nsti_path = os.path.join(out_kegg, "weighted_nsti.xls")
        os.link(nsti_path, out_ko_nsti_path)
        if os.path.exists(ko_abundance_path) and os.path.exists(EC_metagenome_out):
            kegg_cmd = "{} {} -i {} -database {} -e {} -o {}".format(self.python, self.pick_annotation, ko_abundance_path,
                                                                   "kegg", EC_metagenome_out, out_kegg)
            self.logger.info(kegg_cmd)
            command = self.add_command("anno_kegg", kegg_cmd).run()
            self.wait(command)
            if command.return_code in [0]:
                self.logger.info("anno kegg 运行完成！")
            else:
                self.set_error("anno kegg 运行失败！", code="32707407")
        else:
            self.logger.info("kegg 预测失败")

    def run_metacyc_anno(self):
        """
        进行metacyc功能的注释
        :return:
        """
        self.logger.info("开始对metacyc功能进行注释啦")
        out_metacyc = os.path.join(self.work_dir, "MetaCyc")
        if os.path.exists(out_metacyc):
            shutil.rmtree(out_metacyc)
        os.mkdir(out_metacyc)
        metacyc_abundance_path = os.path.join(self.picrust2_out, "pathways_out/path_abun_unstrat.tsv")
        # os.link(metacyc_abundance_path, metacyc_file)
        if os.path.exists(metacyc_abundance_path):
            metacyc_cmd = "{} {} -i {} -database {} -o {}".format(self.python, self.pick_annotation,
                                                                  metacyc_abundance_path, "metacyc", out_metacyc)
            self.logger.info(metacyc_cmd)
            command = self.add_command("anno_metacyc", metacyc_cmd).run()
            self.wait(command)
            if command.return_code in [0]:
                self.logger.info("anno metacyc 运行完成！")
            else:
                self.set_error("anno metacyc 运行失败！", code="32707408")
        else:
            self.logger.info("metacyc 预测失败")

    def run_ec_anno(self):
        """
        对18S注释的酶进行注释
        :return:
        """
        self.logger.info("开始对kegg的酶进行注释啦")
        out_ec = os.path.join(self.work_dir, "KEGG")
        if os.path.exists(out_ec):
            shutil.rmtree(out_ec)
        os.mkdir(out_ec)
        if self.option("analysis_type") in ["18s", "18S"]:
            ec_abundance_path = os.path.join(self.picrust2_out,"ec_18S_counts.txt_metagenome_out/pred_metagenome_unstrat.tsv")
            nsti_path = os.path.join(self.picrust2_out, "ec_18S_counts.txt_metagenome_out/weighted_nsti.tsv")
        else:
            ec_abundance_path = os.path.join(self.picrust2_out,"ec_ITS_counts.txt_metagenome_out/pred_metagenome_unstrat.tsv")
            nsti_path = os.path.join(self.picrust2_out, "ec_ITS_counts.txt_metagenome_out/weighted_nsti.tsv")
        out_nsti_path = os.path.join(out_ec, "weighted_nsti.xls")
        os.link(nsti_path, out_nsti_path)
        if os.path.exists(ec_abundance_path):
            metacyc_cmd = "{} {} -i {} -database {} -o {}".format(self.python, self.pick_annotation,
                                                                  ec_abundance_path, "enzyme", out_ec)
            self.logger.info(metacyc_cmd)
            command = self.add_command("anno_ec", metacyc_cmd).run()
            self.wait(command)
            if command.return_code in [0]:
                self.logger.info("anno ec 运行完成！")
            else:
                self.set_error("anno ec 运行失败！", code="32707409")
        else:
            self.logger.info("ec 预测失败")

    def check_file(self):
        """
        对预测的结果进行检查
        :return:
        """
        self.logger.info("开始对picrust2的结果进行检查！")
        if self.option("analysis_type") in ["16s", "16S"]:
            out_path = self.picrust2_out
            files_dir = ["COG_metagenome_out", "EC_metagenome_out", "KO_metagenome_out", "pathways_out"]
            files = ["COG_predicted.tsv.gz", "EC_predicted.tsv.gz", "KO_predicted.tsv.gz","marker_predicted_and_nsti.tsv.gz", "out.tre"]
            for f in files:
                f_path = os.path.join(out_path, f)
                if not os.path.exists(f_path):
                    self.logger.info("%s文件不存在，picrust2运行的结果失败！请检查选择的数据类型是否正确！"%f_path)
            for j in files_dir:
                j_path = os.path.join(out_path, j)
                if not os.path.exists(j_path):
                    self.logger.info("%s文件不存在，picrust2运行的结果失败！请检查选择的数据类型是否正确！"%j_path)
        elif self.option("analysis_type") in ["18s", "18S"]:
            out_path = self.picrust2_out
            files_dir = ["ec_18S_counts.txt_metagenome_out", "pathways_out"]
            files = ["marker_predicted_and_nsti.tsv.gz", "ec_18S_counts.txt_predicted.tsv.gz", "out.tre"]
            for f in files:
                f_path = os.path.join(out_path, f)
                if not os.path.exists(f_path):
                    self.set_error("%s文件不存在，picrust2运行的结果失败！请检查选择的数据类型是否正确！", variables=(f_path), code="32707410")
            for j in files_dir:
                files_dir_path = os.path.join(out_path, j)
                if not os.path.exists(files_dir_path):
                    self.set_error("%s文件不存在，picrust2运行的结果失败！请检查选择的数据类型是否正确！", variables=(files_dir_path), code="32707411")
                if not os.path.isdir(files_dir_path):
                    self.set_error("%s文件夹不存在，picrust2运行的结果失败！请检查选择的数据类型是否正确！", variables=(files_dir_path), code="32707412")
        elif self.option("analysis_type") in ["ITS", "its", "Its"]:
            out_path = self.picrust2_out
            files_dir = ["ec_ITS_counts.txt_metagenome_out", "pathways_out"]
            files = ["ec_ITS_counts.txt_predicted.tsv.gz", "marker_predicted_and_nsti.tsv.gz", "out.tre"]
            for f in files:
                f_path = os.path.join(out_path, f)
                if not os.path.exists(f_path):
                    self.set_error("%s文件不存在，picrust2运行的结果失败！请检查选择的数据类型是否正确！", variables=(f_path), code="32707413")
            for j in files_dir:
                files_dir_path = os.path.join(out_path, j)
                if not os.path.exists(files_dir_path):
                    self.set_error("%s文件不存在，picrust2运行的结果失败！请检查选择的数据类型是否正确！", variables=(files_dir_path), code="32707414")
                if not os.path.isdir(files_dir_path):
                    self.set_error("%s文件夹不存在，picrust2运行的结果失败！请检查选择的数据类型是否正确！", variables=(files_dir_path), code="32707415")
        self.logger.info("对picrust2的结果检查完毕")

    def rename(self, infile, outfile):
        """
        对中间过程文件（EC_predicted.xls）进行重命令，去除酶的EC
        :param infile:
        :param outfile:
        :return:
        """
        with open(infile, 'r') as f, open(outfile, 'w') as w:
            lines = f.readlines()
            header_list = lines[0].strip().split("\t")
            head = []
            head.append(header_list[0])
            for line in header_list[1:]:
                ec_name = line.split(":")[1]
                head.append(ec_name)
            head_list = "\t".join(head)
            w.write(head_list + "\n")
            for n_line in lines[1:]:
                w.write(n_line)

    def run(self):
        """
        运行
        """
        super(Picrust2PredictTool, self).run()
        self.logger.info("开始运行啦")
        if self.option("analysis_type") in ["16s", "16S"]:
            self.func_predict_analysis()
            self.check_file()
            self.unzip_file()
            self.run_cog_anno()
            self.run_kegg_anno()
            self.run_metacyc_anno()
        elif self.option("analysis_type") in ["18s", "18S", "ITS", "its", "Its"]:
            self.func_predict_analysis()
            self.check_file()
            self.unzip_file()
            self.run_ec_anno()
            self.run_metacyc_anno()
        self.set_output()
        self.end()

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        if self.option("analysis_type") in ["16s", "16S"]:
            cog_dir = os.path.join(self.output_dir, "COG")
            kegg_dir = os.path.join(self.output_dir, "KEGG")
            metacyc_dir = os.path.join(self.output_dir, "MetaCyc")
            if os.path.exists(cog_dir):
                shutil.rmtree(cog_dir)
            if os.path.exists(kegg_dir):
                shutil.rmtree(kegg_dir)
            if os.path.exists(metacyc_dir):
                shutil.rmtree(metacyc_dir)
            cog_list = os.listdir(os.path.join(self.work_dir, "COG"))
            if len(cog_list) != 0:
                shutil.copytree(os.path.join(self.work_dir, "COG"), cog_dir)
            kegg_list = os.listdir(os.path.join(self.work_dir, "KEGG"))
            if len(kegg_list) != 0:
                shutil.copytree(os.path.join(self.work_dir, "KEGG"), kegg_dir)
            meta_list = os.listdir(os.path.join(self.work_dir, "MetaCyc"))
            if len(meta_list) != 0:
                shutil.copytree(os.path.join(self.work_dir, "MetaCyc"), metacyc_dir)
            cog_file = os.path.join(self.output_dir, "COG_predicted.xls")
            ec_file = os.path.join(self.output_dir, "EC_predicted.xls")
            ko_file = os.path.join(self.output_dir, "KO_predicted.xls")
            marker_file = os.path.join(self.output_dir, "marker_predicted_and_nsti.xls")
            tree = os.path.join(self.output_dir, "out.tre")
            if os.path.exists(cog_file):
                os.remove(cog_file)
                os.link(self.picrust2_out + "/COG_predicted.tsv", cog_file)
            else:
                os.link(self.picrust2_out + "/COG_predicted.tsv", cog_file)
            if os.path.exists(ec_file):
                os.remove(ec_file)
                self.rename(self.picrust2_out + "/EC_predicted.tsv",self.picrust2_out + "/EC_predicted2.tsv")
                os.link(self.picrust2_out + "/EC_predicted2.tsv", ec_file)
            else:
                self.rename(self.picrust2_out + "/EC_predicted.tsv",self.picrust2_out + "/EC_predicted2.tsv")
                os.link(self.picrust2_out + "/EC_predicted2.tsv", ec_file)
            if os.path.exists(ko_file):
                os.remove(ko_file)
                os.link(self.picrust2_out + "/KO_predicted.tsv", ko_file)
            else:
                os.link(self.picrust2_out + "/KO_predicted.tsv", ko_file)
            if os.path.exists(marker_file):
                os.remove(marker_file)
                os.link(self.picrust2_out + "/marker_predicted_and_nsti.tsv", marker_file)
            else:
                os.link(self.picrust2_out + "/marker_predicted_and_nsti.tsv", marker_file)
            if os.path.exists(tree):
                os.remove(tree)
                os.link(self.picrust2_out + "/out.tre", tree)
            else:
                os.link(self.picrust2_out + "/out.tre", tree)
        else:
            kegg_dir = os.path.join(self.output_dir, "KEGG")
            metacyc_dir = os.path.join(self.output_dir, "MetaCyc")
            if os.path.exists(kegg_dir):
                shutil.rmtree(kegg_dir)
            if os.path.exists(metacyc_dir):
                shutil.rmtree(metacyc_dir)
            shutil.copytree(os.path.join(self.work_dir, "KEGG"), kegg_dir)
            shutil.copytree(os.path.join(self.work_dir, "MetaCyc"), metacyc_dir)
            if self.option("analysis_type") in ["18s", "18S"]:
                ec_file = os.path.join(self.picrust2_out, "ec_18S_counts.txt_predicted.tsv")
            else:
                ec_file = os.path.join(self.picrust2_out, "ec_ITS_counts.txt_predicted.tsv")
            ec_out = os.path.join(self.output_dir, "EC_predicted.xls")
            if os.path.exists(ec_out):
                os.remove(ec_out)
                self.rename(ec_file, self.picrust2_out + "/EC_predicted2.tsv")
                os.link(self.picrust2_out + "/EC_predicted2.tsv", ec_out)
            else:
                self.rename(ec_file, self.picrust2_out + "/EC_predicted2.tsv")
                os.link(self.picrust2_out + "/EC_predicted2.tsv", ec_out)
            marker_file = os.path.join(self.output_dir, "marker_predicted_and_nsti.xls")
            tree = os.path.join(self.output_dir, "out.tre")
            if os.path.exists(marker_file):
                os.remove(marker_file)
                os.link(self.picrust2_out + "/marker_predicted_and_nsti.tsv", marker_file)
            else:
                os.link(self.picrust2_out + "/marker_predicted_and_nsti.tsv", marker_file)
            if os.path.exists(tree):
                os.remove(tree)
                os.link(self.picrust2_out + "/out.tre", tree)
            else:
                os.link(self.picrust2_out + "/out.tre", tree)
