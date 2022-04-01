# -*- coding: utf-8 -*-
# __author__ = 'chenyanyan'
# last_modifiy:2017.12.27 by zhengyuan, add NSTI
# modified 2016.12.06

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil
import re
import subprocess
import itertools
from mbio.packages.annotation.kegg_annotation import KeggAnnotation
import pandas as pd


class FunctionPredictAgent(Agent):
    """
    16s功能预测
    """
    def __init__(self, parent):
        super(FunctionPredictAgent, self).__init__(parent)

        options = [
            {"name": "otu_reps.fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "otu_table.xls", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "db", "type": "string", "default": "both"}
        ]

        self.add_option(options)
        self.step.add_steps("func_predict")
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.func_predict.start()
        self.step.update()

    def step_end(self):
        self.step.func_predict.finish()
        self.step.update()

    def check_options(self):
        if not self.option("otu_reps.fasta").is_set:
            raise OptionError("请传入otu_reps序列！", code="32704401")
        if not self.option("otu_table.xls").is_set:
            raise OptionError("请传入otu_table.xls文件！", code="32704402")
        if not self.option("db") in ["kegg", 'cog', 'both']:
            raise OptionError("数据库必须为kegg或者cog,both！", code="32704403")

    def set_resource(self):
        self._cpu = 10
        self._memory = '100G'

    def end(self):

        super(FunctionPredictAgent, self).end()


class FunctionPredictTool(Tool):
    def __init__(self, config):
        super(FunctionPredictTool, self).__init__(config)
        """
        设置软件、脚本、数据库的路径
        """
        self.set_environ(PYTHONPATH = self.config.SOFTWARE_DIR + "/bioinfo/meta/picrust-1.1.0/")
        self.python_scripts = self.config.SOFTWARE_DIR + "/program/Python/bin/"
        self.perl = self.config.SOFTWARE_DIR + "/program/perl/perls/perl-5.24.0/bin/perl"
        self.scripts = self.config.SOFTWARE_DIR + "/bioinfo/meta/16s_scripts/"
        self.picrust_path = self.config.SOFTWARE_DIR + "/bioinfo/meta/picrust-1.1.0/scripts/"
        self.Fundb = self.config.SOFTWARE_DIR + "/bioinfo/meta/16sFundb/rep_set/97_otus.fasta"  # 设置数据库文件路径，97_otus.fasta
        self.category_db = self.config.SOFTWARE_DIR + "/bioinfo/meta/16sFundb/database/"
        self.image_magick = self.config.SOFTWARE_DIR + "/program/ImageMagick/bin/convert"
        self.cmds = []
        self.kegg_desc_add = self.config.PACKAGE_DIR + "/meta/scripts/kegg_desc_add.py"
        self.kegg_json = self.config.SOFTWARE_DIR + '/bioinfo/meta/16sFundb/database/kegg_v94.2_ko.detail.json'
        self.kegg_category = self.config.SOFTWARE_DIR + '/bioinfo/meta/16sFundb/database/kegg_v94.2_ko.category.txt'
        self.file_convert = self.config.PACKAGE_DIR + "/meta/scripts/file_convert.py"#add by qingchen.zhang@20191121
        self.ec_database = self.config.SOFTWARE_DIR + '/bioinfo/meta/16sFundb/database/kegg_enzyme_v94.2.xls' ## add by qingchen.zhang @ 20201126
        self.ko_database = self.config.SOFTWARE_DIR + '/bioinfo/meta/16sFundb/database/kegg_v94.2_KO.xls' ## add by qingchen.zhang @ 20201126

    def func_predict_analysis(self):

        keggs = []
        cogs = []
        if not os.path.exists(os.path.join(os.getcwd(), "Results")):
            os.mkdir("Results")
        result_path = os.path.join(os.getcwd(), "Results")

        if not os.path.exists("otu_picking_params_97.txt"):
            self.logger.info("判断txt文件是否存在")
            f = open("otu_picking_params_97.txt", "wb")
            f.write("pick_otus:enable_rev_strand_match True\n")
            f.write("pick_otus:similarity 0.97")
            f.close()
        otu_params_path = os.path.join(os.getcwd(), "otu_picking_params_97.txt")


        cmd1 = "{}pick_closed_reference_otus.py -i {} -o Results -f -p otu_picking_params_97.txt -r {}".format(self.python_scripts, self.option("otu_reps.fasta").prop["path"], self.Fundb)
        self.logger.info(cmd1)

        cmd2 = "{} {}parse_OTU.pl -a {} -b Results/uclust_ref_picked_otus/otu_reps_otus.txt -o new.otu_table.xls".format(self.perl, self.scripts, self.option("otu_table.xls").prop["path"])
        print cmd2

        self.cmds.append(cmd1)
        self.cmds.append(cmd2)
        cmd3 = '{}biom convert -i new.otu_table.xls -o otu_table.biom --table-type "OTU table" --to-hdf5'.format(self.python_scripts)
        self.cmds.append(cmd3)
        cmd4 = '{}python {}normalize_by_copy_number.py -i otu_table.biom -o normalized_otus.biom'.format(self.python_scripts, self.picrust_path)
        self.cmds.append(cmd4)
        cmd5 = '{}biom convert -i normalized_otus.biom -o normalized_otus.txt  --table-type "OTU table" --to-tsv'.format(self.python_scripts)
        self.cmds.append(cmd5)
        cmd6 = "cat normalized_otus.txt|sed -n '2p'|sed 's/#//' >normalized_otus.xls"
        self.cmds.append(cmd6)
        cmd7 = "cat normalized_otus.txt|sed -n '3,$p'|sed 's/\.0	/	/g'|sed 's/\.0$//g'>>normalized_otus.xls"
        self.cmds.append(cmd7)

        cmd8_ko = 'python {}predict_metagenomes.py -i normalized_otus.biom -o predictions_ko.biom -t ko -a NSTI_KEGG.xls'.format(self.picrust_path)
        cmd8_convert = '{}biom convert -i predictions_ko.biom -o predictions_ko.txt --table-type "OTU table" --to-tsv'.format(self.python_scripts)
        cmd8_trans1 = "cat predictions_ko.txt|sed -n '2p'|sed 's/#OTU ID/KO ID/' >predictions_ko.xls"
        cmd8_trans2 = "cat predictions_ko.txt|sed -n '3,$p'|sort -V |sed 's/\.0//g'>>predictions_ko.xls"


        cmd8_trans3 = "{}python {} predictions_ko.xls".format(self.python_scripts, self.file_convert)


        cmd9_cog = '{}python {}predict_metagenomes.py -i normalized_otus.biom -o predictions_cog.biom -t cog -a NSTI_COG.xls'.format(self.python_scripts, self.picrust_path)
        #self.cmds.append(cmd9_cog)
        cmd9_convert = '{}biom convert -i predictions_cog.biom -o predictions_cog.txt --table-type "OTU table" --to-tsv'.format(self.python_scripts)
        #self.cmds.append(cmd9_convert)
        cmd9_trans1 = "cat predictions_cog.txt|sed -n '2p'|sed 's/#OTU ID/COG ID/' >predictions_cog.xls"
        cmd9_trans2 = "cat predictions_cog.txt|sed -n '3,$p'|sort -V |sed 's/\.0//g'>>predictions_cog.xls"
        cmd9_trans3 = '{}python {} predictions_cog.xls'.format(self.python_scripts, self.file_convert)#add by qingchen.zhang@20191121

        """
        计算丰度表
        """

        cmd10 = '{}python {}categorize_by_function.py -i predictions_ko.biom -c KEGG_Pathways -l 1 -o predictions_ko.L1.biom --ignore none'.format(self.python_scripts, self.picrust_path)
        cmd11 = '{}python {}categorize_by_function.py -i predictions_ko.biom -c KEGG_Pathways -l 2 -o predictions_ko.L2.biom'.format(self.python_scripts, self.picrust_path)
        cmd12 = '{}python {}categorize_by_function.py -i predictions_ko.biom -c KEGG_Pathways -l 3 -o predictions_ko.L3.biom'.format(self.python_scripts, self.picrust_path)

        cmd10_convert = '{}biom convert -i predictions_ko.L1.biom -o predictions_ko.L1.txt --table-type "OTU table" --to-tsv'.format(self.python_scripts)
        cmd10_trans1 = "cat predictions_ko.L1.txt|sed -n '2p'|sed 's/#OTU ID/CategoryL1/' >predictions_ko.L1.xls"
        cmd10_trans2 = "cat predictions_ko.L1.txt|sed -n '3,$p'|sed 's/\.0//g'>>predictions_ko.L1.xls"
        cmd10_trans3 = '{}python {} predictions_ko.L1.xls'.format(self.python_scripts, self.file_convert)#add by qingchen.zhang@20191121

        cmd11_convert = '{}biom convert -i predictions_ko.L2.biom -o predictions_ko.L2.txt --table-type "OTU table" --to-tsv'.format(self.python_scripts)
        cmd11_trans1 = "cat predictions_ko.L2.txt|sed -n '2p'|sed 's/#OTU ID/CategoryL2/' >predictions_ko.L2.xls"
        cmd11_trans2 = "cat predictions_ko.L2.txt|sed -n '3,$p'|sed 's/\.0//g'>>predictions_ko.L2.xls"
        cmd11_trans3 = '{}python {} predictions_ko.L2.xls'.format(self.python_scripts, self.file_convert)#add by qingchen.zhang@20191121

        cmd12_convert = '{}biom convert -i predictions_ko.L3.biom -o predictions_ko.L3.txt --table-type "OTU table" --to-tsv'.format(self.python_scripts)
        cmd12_trans1 = "cat predictions_ko.L3.txt|sed -n '2p'|sed 's/#OTU ID/CategoryL3/' >predictions_ko.L3.xls"
        cmd12_trans2 = "cat predictions_ko.L3.txt|sed -n '3,$p'|sed 's/\.0//g'>>predictions_ko.L3.xls"
        cmd12_trans3 = '{}python {} predictions_ko.L3.xls'.format(self.python_scripts, self.file_convert)#add by qingchen.zhang@20191121

        """
        生成kegg.pathway.profile.xls和kegg.enzyme.profile.xls文件
        """

        cmd13 = '{} {}kegg_anno.pl -i predictions_ko.xls -m {} -o kegg -e {} -k {} > kegg_anno.log'.format(self.perl, self.scripts,self.category_db, self.ec_database, self.ko_database)
        cmd14 = '{} {}cog_anno.pl -i predictions_cog.xls -o cog'.format(self.perl, self.scripts)

        cmd14_trans1 = '{}python {} cog.category.function.xls'.format(self.python_scripts, self.file_convert) #add by qingchen.zhang@20191121

        cmd15 = '{}python {}cog-boxplot.py -f cog.category.function.xls -c cog.descrip.table.xls -o cog.box.pdf'.format(self.python_scripts, self.scripts)
        cmd16 = '{}python {}cog_category.py {} {} predictions_cog.xls'.format(self.python_scripts, self.scripts, self.category_db + 'NOG.funccat.txt', self.category_db + 'eggnog.fun_cate.txt')



        #add v4 201910
        cmd_add_kegg_desc = '{}python {} -k predictions_ko.xls -L2 predictions_ko.L2.xls -L3 predictions_ko.L3.xls -k_json {} -k_database {}'.format(self.python_scripts, self.kegg_desc_add,self.kegg_json, self.kegg_category)
        if self.option("db") == "cog" or self.option("db") == "both":
            self.cmds.append(cmd9_cog)
            self.cmds.append(cmd9_convert)
            self.cmds.append(cmd9_trans1)
            self.cmds.append(cmd9_trans2)
            self.cmds.append(cmd9_trans3)
            self.cmds.append(cmd14)
            self.cmds.append(cmd14_trans1)
            self.cmds.append(cmd15)
            self.cmds.append(cmd16)
            self.logger.info("cog分析命令行添加完毕！")
        if self.option("db") == "kegg" or self.option("db") == "both":
            self.cmds.append(cmd8_ko)
            self.cmds.append(cmd8_convert)
            self.cmds.append(cmd8_trans1)
            self.cmds.append(cmd8_trans2)
            self.cmds.append(cmd8_trans3)
            self.cmds.append(cmd10)
            self.cmds.append(cmd10_convert)
            self.cmds.append(cmd10_trans1)
            self.cmds.append(cmd10_trans2)
            self.cmds.append(cmd10_trans3)
            self.cmds.append(cmd11)
            self.cmds.append(cmd11_convert)
            self.cmds.append(cmd11_trans1)
            self.cmds.append(cmd11_trans2)
            self.cmds.append(cmd11_trans3)
            self.cmds.append(cmd12)
            self.cmds.append(cmd12_convert)
            self.cmds.append(cmd12_trans1)
            self.cmds.append(cmd12_trans2)
            self.cmds.append(cmd12_trans3)
            self.cmds.append(cmd13)
            self.cmds.append(cmd_add_kegg_desc)  # add 20191016 v4
            self.logger.info("kegg分析命令行添加完毕！")

    def run_cmds(self, cmds):
            """
            循环运行命令行，传入命令行列表
            """

            self.logger.info("开始运行命令！")
            self.func_predict_analysis()
            count = 0
            for cmd in cmds:
                self.logger.info(cmd)
                try:
                    subprocess.check_output(cmd, shell=True)
                    self.logger.info(cmd + " 运行完成！")
                    count += 1
                except subprocess.CalledProcessError:
                    import traceback
                    print traceback.format_exc()
                    self.logger.info('CMD:{}'.format(cmd))
                    break
            if count == len(cmds):
                self.logger.info("全部命令运行完成！！！")
                self.setoutput()
            else:
                if count == 0:
                    log_path = os.path.join(os.getcwd(), "Results/uclust_ref_picked_otus/otu_reps_otus.log")
                    log_file = open(log_path)
                    for i in log_file:
                            if i.strip() == "Num OTUs:0":
                                self.set_error('16S功能预测需测序区域为16S rRNA，请检测测序引物序列。', code="32704401")
                else:
                    self.set_error('命令cmd{} ERROR，运行失败！', variables=(count - 1), code="32704402")


    def run(self):
        """
        """
        super(FunctionPredictTool, self).run()
        self.logger.info("运行函数run_cmds")
        self.run_cmds(self.cmds)
        self.logger.info("函数run_cmds运行结束！")



    def setoutput(self):
        if self.option("db") == 'kegg' or self.option("db") == "both" :
            if os.path.exists("KEGG"):
                shutil.rmtree("KEGG")
            os.mkdir("KEGG")

            new_file = open("pid.txt", "w+")
            f = open("kegg.pathway.profile.xls")
            for line in itertools.islice(open("kegg.pathway.profile.xls"), 1, None):
                pid = line.split()[0]
                url = line.strip().split()[-1]
                kos = url.split('?')[-1]
                kos = ";".join(kos.split("+")[1:])
                line = pid + "\t" + kos
                new_file.write(line)
                new_file.write('\n')
            new_file.close()
            for f in os.listdir(os.getcwd()):
                if re.match("predictions_ko.*", f) or re.match("predictions_module.xls", f):
                    shutil.move(f, "KEGG")
                if re.match("kegg.*.xls", f):
                    shutil.move(f, "KEGG")
            shutil.move("NSTI_KEGG.xls", "KEGG")
            shutil.move("KEGG", self.output_dir)
            ##fix by qingchen.zhang @20191127
            # os.makedirs("output/KEGG/Pathway_pdf")
            # pdf_dir = os.path.join(os.getcwd(), "output/KEGG/Pathway_pdf")
            # KeggAnnotation().getPic("pid.txt", pdf_dir)
            # os.makedirs("output/KEGG/Pathway_png")
            # png_dir = os.path.join(os.getcwd(), "output/KEGG/Pathway_png")
            # if os.path.exists(pdf_dir):
            #     if not os.path.exists(png_dir):
            #         os.mkdir(png_dir)
            #     for f in os.listdir(pdf_dir):
            #         m = re.match(r"(.*)pdf$", f)
            #         if m:
            #             q = m.group(1)
            #         fp = os.path.join(pdf_dir, f)
            #         fq = png_dir + '/' + q + 'png'
            #         cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' + fp + ' ' + fq
            #         try:
            #             subprocess.check_output(cmd, shell=True)
            #         except subprocess.CalledProcessError:
            #             self.set_error('图片格式pdf转png出错', code="32704403")
            for file in os.listdir(os.path.join(self.output_dir, "KEGG")):
                file_path = os.path.join(self.output_dir + "/KEGG", file)
                if file.endswith("_old"):
                    os.remove(file_path)
                if file.endswith(".txt"):
                    os.remove(file_path)
        if self.option("db") == "cog" or self.option("db") == "both":
            if os.path.exists("COG"):
                shutil.rmtree("COG")
            os.mkdir("COG")
            for f in os.listdir(os.getcwd()):
                if re.match("predictions_cog.*", f):
                    shutil.move(f, "COG")
                if re.match("cog.*.xls", f):
                    shutil.move(f, "COG")
            shutil.move("NSTI_COG.xls", "COG")
            shutil.move("COG", self.output_dir)
        self.logger.info("成功移动文件夹！")
        self.end()