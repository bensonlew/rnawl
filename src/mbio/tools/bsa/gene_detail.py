# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.02.27

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class GeneDetailAgent(Agent):
    """
    BSA基因详情页，根据chrom、start、end得到基因序列信息或变异位点信息
    """
    def __init__(self, parent):
        super(GeneDetailAgent, self).__init__(parent)
        options = [
            # {"name": "index_file", "type": "infile", "format": "bsa.vcf"},  # 0X-0X的index-calc.result.index文件
            {"name": "variant_file", "type": "infile", "format": "bsa.vcf"},  # 0X-0X的index-calc.result.variant文件
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},  # ref.fa
            {"name": "gene_id", "type": "string"},  # 基因 ID
            {"name": "chrom", "type": "string"},  # 基因所在的染色体
            {"name": "start", "type": "int"},  # 基因起始位置
            {"name": "end", "type": "int"},  # 基因终止位置
            {"name": "type", "type": "string"},  # 查找的基因类型, seq/detail
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("type"):
            raise OptionError("请设置查找的基因类型", code="31500301")
        if self.option("type") == "detail":
            if not self.option("variant_file").is_set:
                raise OptionError("请设置0X-0X的index-calc.result.variant文件", code="31500302")
            if not self.option("gene_id"):
                raise OptionError("请设置基因ID", code="31500303")
        elif self.option("type") == "seq":
            if not self.option("ref_fa").is_set:
                raise OptionError("请设置ref_fa文件", code="31500304")
            if not self.option("chrom"):
                raise OptionError("请设置查找的染色体", code="31500305")
            if not self.option("start"):
                raise OptionError("请设置基因起始位置", code="31500306")
            if not self.option("end"):
                raise OptionError("请设置基因终止位置", code="31500307")
        else:
            raise OptionError("类型只能为detail/seq", code="31500308")

    def set_resource(self):
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        super(GeneDetailAgent, self).end()


class GeneDetailTool(Tool):
    def __init__(self, config):
        super(GeneDetailTool, self).__init__(config)
        self.gene_detail = self.config.PACKAGE_DIR + "/bsa/gene_detail.py"
        self.python = "miniconda2/bin/python"

    def run_find_gene_detail(self):
        """
        查找变异位点信息
        """
        cmd = "{} {} -i {} ".format(self.python, self.gene_detail, self.option("variant_file").prop["path"])
        # cmd = "{} {} -i {} ".format(self.python, self.gene_detail, self.option("index_file").prop["path"])
        # cmd += "-t {} -chr {} -start {} ".format("detail", self.option("chrom"), self.option("start"))
        # cmd += "-end {} -o {}".format(self.option("end"), self.output_dir + "/filter.result.index")
        cmd += "-t {} -g {} -o {} ".format("detail", self.option("gene_id"), self.output_dir + "/filter.result.index")
        command = self.add_command("gene_detail", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("查找基因变异位点信息运行完成")
        else:
            api_bgd = self.api.api("bsa.bsa_gene_detail")
            api_bgd.update_sg_gene_status(main_id=self.option("main_id"), collection="sg_gene_index")
            self.set_error("查找基因变异位点信息运行出错!", code="31500301")

    def find_gene_seq(self):
        """
        查找基因序列信息
        """
        cmd = "{} {} -i {} ".format(self.python, self.gene_detail, self.option("ref_fa").prop["path"])
        cmd += "-t {} -chr {} -start {} ".format("seq", self.option("chrom").lower(), self.option("start"))
        cmd += "-end {} -o {}".format(self.option("end"), self.output_dir + "/seq.fa")
        command = self.add_command("gene_seq", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("查找基因序列信息运行完成")
        else:
            api_bgd = self.api.api("bsa.bsa_gene_detail")
            api_bgd.update_sg_gene_status(main_id=self.option("main_id"), collection="sg_gene_seq")
            self.set_error("查找基因序列信息运行出错!", code="31500302")

    def set_db(self):
        self.logger.info("保存结果到mongo")
        api_bgd = self.api.api("bsa.bsa_gene_detail")
        if self.option("type") == "detail":
            api_bgd.add_sg_gene_index_detail(index_id=self.option("main_id"), index_path=self.output_dir + "/filter.result.index")
        else:
            api_bgd.add_sg_gene_seq(main_id=self.option("main_id"), seq_path=self.output_dir + "/seq.fa")

    def run(self):
        super(GeneDetailTool, self).run()
        if self.option("type") == "detail":
            self.run_find_gene_detail()
        else:
            self.find_gene_seq()
        self.set_db()
        self.end()
