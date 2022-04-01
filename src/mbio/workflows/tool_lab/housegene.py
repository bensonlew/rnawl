# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir
from mbio.packages.tool_lab.common_function import fasta_dir_sort
import os
import subprocess
import gevent

class HousegeneWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(HousegeneWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "file_style", "type": "string"},
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "fasta_dir", "type": "infile", "format": "sequence.fasta_dir"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.modules = []
        #self.gene_predict = self.add_module("bac_comp_genome.gene_predict")
        #self.core_gene = self.add_module("bac_comp_genome.get_coregene")

    def check_options(self):
        if (not self.option("fasta").is_set) and (not self.option("fasta_dir").is_set):
            raise OptionError("请传入序列文件或序列文件夹！", code="")

    def run(self):
        if self.option("fasta_dir").is_set:
            fasta_dir_sort(self.option("fasta_dir").prop['path'],self.work_dir + "/fasta_dir_sort")
        self.run_hgene()
        super(HousegeneWorkflow, self).run()

    def run_hgene(self):
        """
        持家基因预测
        :return:
        """
        if self.option("fasta").is_set:
            self.hgene = self.add_module("tool_lab.housegene")
            opts = {
                "fasta": self.option("fasta"),
                "sample_name": self.option("fasta").prop['path'].split("/")[-1].split(".")[0]
            }
            self.hgene.set_options(opts)
            self.hgene.on('end', self.set_output)
            self.hgene.run()
        else:
            files = os.listdir(os.path.join(self.work_dir, "fasta_dir_sort"))
            for file in files:
                self.hgene = self.add_module("tool_lab.housegene")
                opts= {
                    "fasta": os.path.join(os.path.join(self.work_dir,"fasta_dir_sort"),file),
                    "sample_name": file.strip(".fasta")
                }
                self.hgene.set_options(opts)
                self.modules.append(self.hgene)
            self.logger.info(self.modules)
            if self.modules:
                if len(self.modules) > 1:
                    self.on_rely(self.modules, self.set_output)
                elif len(self.modules) == 1:
                    self.modules[0].on('end', self.set_output)
                for module in self.modules:
                    self.logger.info(module)
                    gevent.sleep(1)
                    module.run()

    def set_output(self):
        stat_file = self.output_dir + "/all.hgene.xls"
        if self.option("fasta").is_set:
            link_dir(self.hgene.output_dir ,self.output_dir)
            print self.output_dir + "/stat_file.txt"
            with open(stat_file,"w") as w,open(self.output_dir + "/stat_file.txt","r") as f:
                data = f.readlines()[0]
                w.write("Genome" + "\t" + "Gene num" + "\n")
                w.write(data.split("\t")[0] + "\t" +data.split("\t")[1] + "\n")
                w.close()
                f.close()
            os.remove(self.output_dir + "/stat_file.txt")
        else:
            with open(stat_file,"w") as w:
                w.write("Genome" + "\t" + "Gene num" + "\n")
                for module in self.modules:
                    self.logger.info(module)
                    with open(module.output_dir + "/stat_file.txt", "r") as f:
                        self.logger.info(module.output_dir + "/stat_file.txt")
                        sample = f.readlines()[0].split("\t")
                        w.write(sample[0] + "\t" + sample[1] + "\n")
                        f.close()
                    if os.path.exists(module.output_dir + '/' + sample[0]):
                        link_dir(module.output_dir + '/' + sample[0], self.output_dir + '/' + sample[0])

                w.close()
        self.set_db()

    def set_db(self):
        api_main = self.api.api("tool_lab.housegene")
        api_main.add_hgene_detail(main_id=self.option('main_id'), file_path=self.output_dir)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "dir", "持家基因预测结果目录", 0],
            ["./all.hgene.xls", "xls", "持家基因预测结果统计表", 0],
            ["./*", "dir", "持家预测结果目录", 0],
            ["./*/*blast.xls", "xls", "持家基因预测结果详情表", 0],
            ["./*/*_hgene.fa", "fa", "持家基因氨基酸序列", 0],
            ["./*/*all.hgene.fa", "fa", "持家基因按照规定顺序连接的蛋白序列文件", 0]
        ])
        super(HousegeneWorkflow, self).end()