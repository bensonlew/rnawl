# -*- coding: utf-8 -*-
# __author__ = 'zzg'
# version 1.0
# last_modify: 2020.03.01

import os, time
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir

class PlasmidPredictAgent(Agent):
    """
    质粒预测小工具
    采用conda进行
    """
    def __init__(self, parent):
        super(PlasmidPredictAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},# 序列文件
            {"name": "sample_name", "type": "string"},# 样本名
            {"name": "threshold", "type": "float", "default": 0.7},
            {"name": "outfmt", "type": "int", "default": 6},
            {"name": "evalue", "type": "float", "default": 1e-5},
            {"name": "max_target_seqs", "type": "int", "default": 10},
        ]
        self.add_option(options)
        #self._memory_increase_step = 20

    def check_options(self):
        if not self.option("fasta").is_set:
            raise OptionError("必须输入序列文件")
        if not self.option("fasta").is_set:
            raise OptionError("必须输入样品名称")

    def set_resource(self):
        self._cpu = 4
        self._memory = '20G'

    def end(self):
        super(PlasmidPredictAgent, self).end()

class PlasmidPredictTool(Tool):
    def __init__(self, config):
        super(PlasmidPredictTool, self).__init__(config)
        #self.seqstat = self.config.SOFTWARE_DIR + "/bioinfo/seq/scripts/seqstat"
        self.conda = self.config.SOFTWARE_DIR + "/program/miniconda3/"
        self.set_environ(PATH = self.conda +"bin")
        self.path = self.config.SOFTWARE_DIR + "/program/Python35/bin:"
        self.lib = self.config.SOFTWARE_DIR + "/program/Python35/lib"
        self.set_environ(PATH=self.path, LD_LIBRARY_PATH=self.lib)
        self.shell_path = "/program/sh"
        #self.blast = self.config.SOFTWARE_DIR + "bioinfo/align/ncbi-blast-2.3.0+/bin/"
        self.database = self.config.SOFTWARE_DIR + "/database"


    def creat_shell(self, input_file, threshold):
        """
        重写一个shell，将运行的命令直接写入文件里面去
        :return:
        """
        self.fix_sh = os.path.join(self.work_dir, "plas_flow.sh")
        with open(self.fix_sh, "w") as w:
            w.write("#! /bin/bash\n")
            w.write("source activate plasflow-1.1\n")
            #w.write("source {} {}\n".format(self.activate, self.roary))
            w.write("PlasFlow.py --input {} --output plasflow_predictions.tsv --threshold {} > plasflow.log\n".format(input_file,threshold))
            w.write("conda deactivate\n")

    def run_predict(self):
        """
        质粒鉴定
        :return:
        """
        self.creat_shell(self.option("fasta").prop['path'], self.option("threshold"))
        #cmd = "{} --input {} --output plasflow_predictions.tsv --threshold {}" .format("bioinfo/tool_lab/PlasFlow-1.1/PlasFlow-master/PlasFlow.py", self.option("fasta").prop['path'], self.option("threshold"))
        cmd = "{} {}".format(self.shell_path, self.fix_sh)
        self.logger.info(cmd)
        command = self.add_command("run_predict", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_predict运行完成！")
        else:
            self.set_error("run_predict运行运行出错!")
        self.plas = os.path.getsize(self.work_dir + "/plasflow_predictions.tsv_plasmids.fasta")
        return self.plas

    def run_anno(self):
        """
        质粒注释
        :return:
        """
        cmd = "{} -query {}/plasflow_predictions.tsv_plasmids.fasta -db {} -out plasmid_blast.txt -outfmt {} -evalue {} -max_target_seqs {}".format(
            "bioinfo/align/ncbi-blast-2.3.0+/bin/blastn", self.work_dir, self.database + "/PLSDB/plsdb.fna",
            self.option("outfmt"), self.option("evalue"), self.option("max_target_seqs"))
        self.logger.info(cmd)
        command = self.add_command("run_anno", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_anno运行完成！")
        else:
            self.set_error("run_anno运行运行出错!")

    def run_stat(self):
        self.plsdb_file = self.output_dir + '/' + self.option("sample_name") + '_' + 'plasmid_plsdb.tsv'
        all_acc = []
        with open(self.work_dir + '/plasmid_blast.txt') as f, open(self.database + '/PLSDB/plsdb.tsv')as v, open(self.plsdb_file,'w') as t:
            t.write("UID_NUCCORE" + "\t" + "ACC_NUCCORE" + "\t" + "Description_NUCCORE" + "\t" + "CreateDate_NUCCORE" + "\t" + "Topology_NUCCORE" + "\t" + "Completeness_NUCCORE" + "\t" + "TaxonID_NUCCORE" + "\t" + "Genome_NUCCORE" + "\t" + "Length_NUCCORE" + "\t" + "Source_NUCCORE" + "\t" + "UID_ASSEMBLY" + "\t" + "Status_ASSEMBLY" + "\t" + "SeqReleaseDate_ASSEMBLY" + "\t" + "SubmissionDate_ASSEMBLY" + "\t" + "Latest_ASSEMBLY" + "\t" + "UID_BIOSAMPLE" + "\t" + "ACC_BIOSAMPLE" + "\t" + "Location_BIOSAMPLE" + "\t" + "Coordinates_BIOSAMPLE" + "\t" + "IsolationSource_BIOSAMPLE" + "\t" + "Host_BIOSAMPLE" + "\t" + "SamplType_BIOSAMPLE" + "\t" + "taxon_name" + "\t" + "taxon_rank" + "\t" + "lineage" + "\t" + "taxon_species_id" + "\t" + "taxon_species_name" + "\t" + "taxon_genus_id" + "\t" + "taxon_genus_name" + "\t" + "taxon_family_id" + "\t" + "taxon_family_name" + "\t" + "taxon_order_id" + "\t" + "taxon_order_name" + "\t" + "taxon_class_id" + "\t" + "taxon_class_name" + "\t" + "taxon_phylum_id" + "\t" + "taxon_phylum_name" + "\t" + "taxon_superkingdom_id" + "\t" + "taxon_superkingdom_name" + "\t" + "loc_lat" + "\t" + "loc_lng" + "\t" + "loc_parsed" + "\t" + "GC_NUCCORE" + "\t" + "Identical" + "\t" + "OldVersion" + "\t" + "hits_rMLST" + "\t" + "hitscount_rMLST" + "\t" + "D1" + "\t" + "D2" + "\t" + "plasmidfinder" + "\t" + "pmlst" + "\n")
            data1 = f.readlines()
            for i in data1:
                all_acc.append(i.split("\t")[1])
            data2 = v.readlines()
            for k in all_acc:
                for x in data2:
                    if k == x.split("\t")[1]:
                        t.write(x)

    def set_output(self):
        """
        设置结果文件目录
        """
        file_list = ["plasflow_predictions.tsv","plasflow_predictions.tsv_chromosomes.fasta","plasflow_predictions.tsv_plasmids.fasta","plasflow_predictions.tsv_unclassified.fasta","plasmid_blast.txt","plasflow.log"]
        for file in file_list:
            if os.path.exists(self.output_dir + '/' +  self.option("sample_name") + '_' + file):
                os.remove(self.output_dir + '/' + self.option("sample_name") + '_' + file)
            if os.path.exists(self.work_dir + '/' + file):
                os.link(self.work_dir + '/' + file,self.output_dir + '/' + self.option("sample_name") + '_' +file)
        self.logger.info("生成结果文件完成")

    def run(self):
        """
        运行
        """
        super(PlasmidPredictTool, self).run()
        plas = self.run_predict()
        if plas == 0:
            self.end()
        else:
            self.run_anno()
            self.run_stat()
        self.set_output()
        self.end()