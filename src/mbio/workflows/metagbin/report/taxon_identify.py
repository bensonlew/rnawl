# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
# last_modifies 20180621

import os
from biocluster.workflow import Workflow
from biocluster.file import download

class TaxonIdentifyWorkflow(Workflow):
    """
    bins 共线性分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(TaxonIdentifyWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "bin_name", "type": "string"},
            {"name": "method", "type": "string"},
            {"name": "query", "type": "infile", "format": "sequence.fasta"},
            {"name": "ref", "type": "infile", "format": "sequence.fasta,sequence.fasta_dir"},
            {"name": "ref_str", "type": "string"},  # 种水平鉴定可能传来一组基因组
            {"name": "blasr_db", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "task_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.colline = self.add_tool("fungi_genome.colline")  # 这个tool换一个，改下结果文件的格式
        self.other = self.add_module("annotation.bin_tax")  # 做blasr,pocp,ani分析
        self.api_obj = self.api.api("metagbin.common_api")

    def run_colline(self):
        self.colline.set_options({
            "qfna": self.option("query").path,
            "rfna": self.option("ref").path,
            "sample": self.option("bin_name")
        })
        self.colline.on("end", self.set_colline_db)
        self.colline.run()

    def run_blasr(self):
        self.other.set_options({
            "query": self.option("query"),
            "blasr_db": self.option("blasr_db"),
            "method": "blasr",
            "sample": self.option("bin_name"),
            "task_id": self.option("task_id")
        })
        self.other.on("end", self.set_blasr_db)
        self.other.run()

    def run_pocp(self):
        self.other.set_options({
            "query": self.option("query"),
            "ref": self.option("ref"),
            "method": "pocp",
            "task_id": self.option("task_id")
        })
        self.other.on("end", self.set_pocp_db)
        self.other.run()


    def run_ani(self):
        opts = {
            "query": self.option("query"),
            "method": "ani",
            "task_id": self.option("task_id")
        }
        if self.option("ref") and self.option("ref").is_set:
            opts["ref"] = self.option("ref")
        else:
            ref_path = os.path.join(self.work_dir, "download_file")
            for seq in self.option("ref_str").split(","):
                file_name = os.path.basename(seq)
                if os.path.isfile(os.path.join(ref_path, file_name)):
                    continue
                download(seq, os.path.join(ref_path, file_name))
            opts["ref"] = ref_path
        self.other.set_options(opts)
        self.other.on("end", self.set_ani_db)
        self.other.run()

    def set_colline_db(self):
        file = os.path.join(self.colline.output_dir, self.option("bin_name") + "_block_new.xls")
        if not os.path.isfile(file) or os.path.getsize(file) == 0:
            # self.set_error("circos结果为空")
            self.end()
        new_file = os.path.join(self.work_dir, "export_colline.xls")
        ref_name = self.option("ref").path.split("/")[-1].split('.')[0].split("._CDS")[0]
        with open(file, "r") as file1, open(new_file, "w") as file2:
            lines = file1.readlines()
            if len(lines) == 0:
                self.end()
            for line in lines:
                line = line.strip().split("\t")
                self.logger.info(line)
                line[0] = line[0].replace(" ", "")
                line[1] = line[1][2:]
                if line[6] == 'seq1':
                    line[6] = self.option("bin_name")
                else:
                    line[6] = ref_name
                file2.write("\t".join(line) + "\n")
        self.api_obj.add_main_detail(new_file, "colline_detail", self.option("main_id"),
                                     "block,location,strand,start,end,len,genome_id,loc_start,loc_end", has_head=False,
                                     main_name="circos_id", main_table="colline",
                                     update_dic={"cir_path": self.sheet.output + "circos.png"})
        # 注意结果表第一列需要去空格，第二列需要去掉前两个字符
        self.end()

    def set_blasr_db(self):
        file = os.path.join(self.other.output_dir, "blasr_result.xls")
        if os.path.isfile(file):
            self.api_obj.add_main_detail(file, "identif_16s_detail", self.option("main_id"), "genome_id,taxon,identify", has_head=True, main_name="s16_id", convert_dic={"identify": float})
            self.api_obj.update_genome_taxon(self.option("bin_name"), file, task_id=self.option("task_id"))
        self.end()

    def set_pocp_db(self):
        file = os.path.join(self.other.output_dir, "pocp_result.xls")
        if os.path.isfile(file):
            self.api_obj.add_main_detail(file, "identif_genus_detail", self.option("main_id"), "genome_id,ref,ref_genus,prop", has_head=True, main_name="genus_id", convert_dic={"prop": float})
        self.end()

    def set_ani_db(self):
        file = os.path.join(self.other.output_dir, "ani_result.xls")
        if os.path.isfile(file):
            self.api_obj.add_main_detail(file, "identif_species_detail", self.option("main_id"), "genome_id,ref,ref_species,ani", has_head=True, main_name="species_id", convert_dic={"ani": float})
        # else:
        #     self.set_error("ani 计算结果为空")
        self.end()

    def end(self):
        bin = self.option("bin_name")
        repaths = []
        if self.option("method") == "colline":
            repaths.append([".", "", "共线性结果目录"])
            regexps = [
                [r'.*_blocks_coords.txt', 'txt', '%s共线性区块的结果信息表' % bin],
                [r'.*_coverage_report.txt', 'txt','%s共线性coverage的结果信息表' % bin],
                [r'.*circos.png', 'png', "%s的基因组共线性SVG统计图" % bin],
                [r'.*circos.svg', 'svg','%s的基因组共线性PNG统计图' % bin],
                [r'.*block_new.xls', 'xls','%s的基因组共线性分析结果表' % bin]
            ]
            sdir = self.add_upload_dir(self.colline.output_dir)
            sdir.add_relpath_rules(repaths)
            sdir.add_regexp_rules(regexps)
        else:
            if self.option("method") == "blasr":
                repaths.append([".", "", "16S物种分类结果目录"])
            elif self.option("method") == "pocp":
                repaths.append([".", "", "POCP分析结果目录"])
            elif self.option("method") == "ani":
                repaths.append([".", "", "ANI分析结果目录"])
            regexps = [
                [r'.*result.xls', 'xls', '结果文件']
            ]
            sdir = self.add_upload_dir(self.other.output_dir)
            sdir.add_relpath_rules(repaths)
            sdir.add_regexp_rules(regexps)
        super(TaxonIdentifyWorkflow, self).end()

    def run(self):
        if self.option("method") == "colline":
            self.run_colline()
        elif self.option("method") == "blasr":
            self.run_blasr()
        elif self.option("method") == "pocp":
            self.run_pocp()
        elif self.option("method") == "ani":
            self.run_ani()
        self.api_obj.update_status_genome_id(self.option("main_id"), self.option("bin_name"))
        super(TaxonIdentifyWorkflow, self).run()