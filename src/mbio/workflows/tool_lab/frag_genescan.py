# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir
from mbio.packages.tool_lab.common_function import fasta_dir_sort
import os
import re, shutil
import subprocess
import gevent

class FragGenescanWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(FragGenescanWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "seq_style", "type": "string", "default":"full"},
            {"name": "file_style", "type": "string"},
            {'name': 'fasta', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'fasta_dir', 'type': 'infile', 'format': 'sequence.fasta_dir'},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.tools = []

    def check_options(self):
        if not self.option("fasta").is_set and not self.option("fasta_dir").is_set:
            raise OptionError("请传入序列文件或序列文件夹！")

    def run(self):
        if self.option("fasta_dir").is_set:
            fasta_dir_sort(self.option("fasta_dir").prop["path"], self.work_dir + "/fasta_dir_sort")
        self.run_fraggenescan()
        super(FragGenescanWorkflow, self).run()

    def run_fraggenescan(self):
        """
        FragGeneScan基因预测
        :return:
        """
        if self.option("fasta").is_set:
            self.all_sample = [os.path.basename(self.option("fasta").prop["path"]).split(".")[0]]
            self.fraggene = self.add_tool("tool_lab.frag_genescan")
            opts = {
                "fasta": self.option("fasta"),
                "sample_name": os.path.basename(self.option("fasta").prop["path"]).split(".")[0]
            }
            self.fraggene.set_options(opts)
            self.fraggene.on('end', self.set_output)
            self.fraggene.run()
        else:
            self.all_sample = []
            files = os.listdir(os.path.join(self.work_dir, "fasta_dir_sort"))
            for file in files:
                self.all_sample.append(file.replace(".fasta",""))
                self.fraggene = self.add_tool("tool_lab.frag_genescan")
                opts = {
                    "fasta": os.path.join(os.path.join(self.work_dir, "fasta_dir_sort"), file),
                    "sample_name": file.replace(".fasta","")
                }
                self.fraggene.set_options(opts)
                self.tools.append(self.fraggene)
            self.logger.info(self.tools)
            if self.tools:
                if len(self.tools) > 1:
                    self.on_rely(self.tools, self.set_output)
                elif len(self.tools) == 1:
                    self.tools[0].on('end', self.set_output)
                for tool in self.tools:
                    self.logger.info(tool)
                    gevent.sleep(1)
                    tool.run()

    def set_output(self):
        if self.option("fasta").is_set:
            for file in os.listdir(self.fraggene.output_dir):
                if os.path.exists(self.output_dir + "/" + file):
                    os.remove(self.output_dir + "/" + file)
                os.link(self.fraggene.output_dir + "/" + file, self.output_dir + "/" + file)
        else:
            for tool in self.tools:
                for file in os.listdir(tool.output_dir):
                    if os.path.exists(self.output_dir + "/" + file):
                        os.remove(self.output_dir + "/" + file)
                    os.link(tool.output_dir + "/" + file, self.output_dir + "/" + file)
        sample_result = []
        with open(self.output_dir + "/GenePredictionDetail.xls","w") as t, open(self.output_dir + "/GenePredictionStat.xls","w") as v:
            t.write("Gene ID\tSample Name\tLocation\tStrand\tStart\tEnd\tGene Len（bp）\tProtein Len\n")
            v.write("Sample Name\tGene No.\tGene Total Len（bp）\tGene Average Len（bp）\tGC Content  (%)\tGene/Genome (%)\n")
            for file in os.listdir(self.output_dir):
                num = 0
                total_len = 0
                gc = 0
                all = 0
                genome_all_len = 0
                if file.endswith(".gff"):
                    sample_result.append(file.replace(".gff",""))
                    if self.option("fasta").is_set:
                        fasta_file = self.option("fasta").prop["path"]
                    else:
                        fasta_file = os.path.join(os.path.join(self.work_dir, "fasta_dir_sort"), file.replace(".gff","") + ".fasta")
                    with open(self.output_dir + "/" + file) as f,open(self.output_dir + "/" + file.replace(".gff","") + ".ffn") as g, open(fasta_file) as y:
                        data = f.readlines()
                        data2 = g.readlines()
                        data3 = y.readlines()
                        for i in data[1:]:
                            num +=1
                            total_len += (abs(int(i.strip().split("\t")[4]) - int(i.strip().split("\t")[3])) + 1)
                            t.write("ORF" + str(num) + "\t" + file.replace(".gff","") + "\t" + i.strip().split("\t")[0]
                                    + "\t" + i.strip().split("\t")[6] + "\t" + i.strip().split("\t")[3] + "\t" + i.strip().split("\t")[4]
                                    + "\t" + str(abs(int(i.strip().split("\t")[4]) - int(i.strip().split("\t")[3])) + 1)
                                    + "\t" + str(int((abs(int(i.strip().split("\t")[4]) - int(i.strip().split("\t")[3])) + 1)/3)) + "\n")

                        for x in data2:
                            if x.startswith(">"):
                                pass
                            else:
                                gc += (x.count("G") + x.count("g") + x.count("c") + x.count("C"))
                                all += len(x.strip())

                        for b in data3:
                            if b.startswith(">"):
                                pass
                            else:
                                genome_all_len += len(b.strip())

                    v.write(file.replace(".gff","") + "\t" + str(num) + "\t" + str(total_len) + "\t" +str(round(total_len/num,2))
                            + "\t" + str(round((float(gc)/all)*100,2)) + "\t" + str(round((float(all)/genome_all_len)*100,2)) + "\n")
            for sample in self.all_sample:
                if sample not in sample_result:
                    v.write(sample + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\n")
        self.set_db()

    def set_db(self):
        api_main = self.api.api("tool_lab.frag_genescan")
        if self.option('main_id'):
            main_id =  self.option('main_id')
        else:
            main_id = api_main.add_frag_main()
        api_main.add_frag_detail(main_id=main_id, detail_path=self.output_dir + "/GenePredictionDetail.xls", stat_file=self.output_dir + "/GenePredictionStat.xls")
        api_main.add_frag_main(main_id)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "dir", "FragGeneScan基因预测结果目录", 0],
            ["./GenePredictionStat.xls", "xls", "FragGeneScan基因预测统计表", 0],
            ["./GenePredictionDetail.xls", "xls", "FragGeneScan基因预测详情表", 0],
            ["./*.ffn", "ffn", "FragGeneScan基因预测ffn文件", 0],
            ["./*.faa", "faa", "FragGeneScan基因预测faa文件", 0],
            ["./*.gff", "gff", "FragGeneScan基因预测gff文件", 0],
        ])
        super(FragGenescanWorkflow, self).end()