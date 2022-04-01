# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
from mbio.packages.ref_rna.filter_gtf import FilterGtf
import re
import pandas as pd
from BCBio import GFF
from BCBio.GFF import GFFExaminer
from gtfparse import read_gtf
import unittest


class ExtractGffFastaAgent(Agent):
    """
    用于提取gff/gtf文件中的mRNA、cds以及pep序列
    """
    def __init__(self, parent):
        super(ExtractGffFastaAgent, self).__init__(parent)
        options = [
            {"name": "genome", "type": "infile", "format": "ref_rna_v2.fasta"},  # 参考基因文件
            {"name": "gff", "type": "infile", "format": "ref_rna_v2.gtf"},  # 参考基因的注释文件
            {"name": "strand", "type": "string", "default": ''},  # + or -
            {"name": "chr", "type": "string", "default": ''},  # 指定染色体名称
            {"name": "start", "type": "int", "default": ''},  # 起始位置
            {"name": "end", "type": "int", "default": ''},  # 终止为止
            {"name": "seq_type", "type": "string", "default": ''},  # 序列类型, mRNA/CDS/protein/all
            {"name": "reverse", "type": "bool", "default": ''},
            {"name": "no_pseudo", "type": "bool", "default": ''},  # 过滤掉含有 'pseudo' 的注释信息
            {"name": "no_cds", "type": "bool", "default": ''},  # 丢弃掉无CDS的转录本
        ]
        self.add_option(options)
        self.step.add_steps("extract_seq")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.extract_seq.start()
        self.step.update()

    def stepfinish(self):
        self.step.extract_seq.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('genome').is_set:
            raise OptionError('必须输入参考基因组序列文件')
        if not self.option('gff').is_set:
            raise OptionError('必须输入参考注释文件gff文件')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = "15G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(ExtractGffFastaAgent, self).end()


class ExtractGffFastaTool(Tool):
    def __init__(self, config):
        super(ExtractGffFastaTool, self).__init__(config)
        self._version = "v1.0.1"
        self.cufflinks_path = '/bioinfo/rna/cufflinks-2.2.1/'
        self.bioawk_path = Config().SOFTWARE_DIR + '/bioinfo/seq/bioawk/'
        self.bedtools_path = Config().SOFTWARE_DIR + '/bioinfo/seq/bedtools-2.25.0/bin/'
        self.script_path = '/bioinfo/rna/scripts/'
        self.gffread_path = "bioinfo/rna/cufflinks-2.2.1/gffread"

    def run(self):
        """
        运行
        :return:
        """
        super(ExtractGffFastaTool, self).run()
        self.run_gff_to_gtf()
        if self.option("chr"):
            df = read_gtf(self.gtf)
            self.chrs = list(set(df['seqname']))
            if self.option("chr") not in self.chrs:
                raise OptionError('染色体名称不存在')
            if not (self.option("start") and self.option("end")):
                filter_gtf = os.path.join(self.work_dir, os.path.basename(self.option("gff").prop['path']) + ".filter.gtf")
                with open(self.option("gff").prop['path'], "r") as f, open(filter_gtf, "w") as w:
                    for line in f:
                        if line.startswith("#"):
                            w.write(line)
                        else:
                            if line.strip().split("\t")[0] == self.option("chr"):
                                w.write(line)
                self.gtf = filter_gtf
        self.run_gtf_to_fa()
        self.set_output()
        self.end()

    def run_gff_to_gtf(self):
        self.logger.info("转换gff文件为gtf文件")
        gff = self.option("gff").prop['path']
        tmp_gtf = os.path.join(self.work_dir, os.path.basename(gff) + ".tmp.gtf")
        self.gtf = os.path.join(self.work_dir, os.path.basename(gff) + ".gtf")
        genome = self.option("genome").prop["path"]
        to_gtf_cmd = '%s %s -g %s -T -O -o %s' % (self.gffread_path, gff, genome, tmp_gtf)
        command = self.add_command("gff_to_gtf", to_gtf_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行gff_to_gtf完成")
        else:
            self.set_error("运行gff_to_gtf运行出错!")
            return False
        ## 过滤gtf文件，去掉不含gene_id、transcript_id、链信息的注释行
        with open(tmp_gtf, "r") as f1, open(self.gtf, "w") as w1:
            for line in f1:
                items = line.strip().split("\t")
                if len(items) != 9:
                    continue
                if not re.search(r'transcript_id', items[8]):
                    continue
                if not re.search(r'gene_id', items[8]):
                    continue
                if not re.search(r'[+|-]', items[6]):
                    continue
                if not line.strip().endswith(";"):
                    w1.write(line.strip() + ";\n")
                else:
                    w1.write(line)

    def run_gtf_to_fa(self):
        """
        运行gtf_to_fasta，转录本gtf文件转fa文件
        """
        cmd = self.cufflinks_path + "gffread"
        if self.option("no_cds"):
            cmd += " -C"
        if self.option("no_pseudo"):
            cmd += " --no-pseudo"
        limit = ""
        if self.option("chr") and self.option("start") and self.option("end"):
            if self.option("strand"):
                limit = self.option("strand") + self.option("chr") + ":" + str(self.option("start")) + ".." + str(self.option("end"))
            else:
                limit = self.option("chr") + ":" + str(self.option("start")) + ".." + str(self.option("end"))
        if limit:
            if self.option("reverse"):
                cmd += " -R"
            cmd += " -r %s %s -g %s -w %s" % (limit, self.gtf, self.option('genome').prop['path'], "mRNA.fa")
        else:
            cmd += " %s -g %s -w %s" % (self.option("gff").prop['path'], self.option('genome').prop['path'], "mRNA.fa")

        self.logger.info('运行gtf_to_transcripts，形成fasta文件')
        command = self.add_command("gtf_to_fa_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gtf_to_fasta运行完成")
        else:
            self.set_error("gtf_to_fasta运行出错!")

        cmd = self.cufflinks_path + "gffread"
        limit = ""
        if self.option("chr") and self.option("start") and self.option("end"):
            if self.option("strand"):
                limit = self.option("strand") + self.option("chr") + ":" + str(self.option("start")) + ".." + str(self.option(
                    "end"))
            else:
                limit = self.option("chr") + ":" + str(self.option("start")) + ".." + str(self.option("end"))
        if limit:
            if self.option("reverse"):
                cmd += " -R"
            cmd += " -r %s %s -g %s -x %s -y %s" % (limit, self.gtf, self.option('genome').prop['path'], "cds.fa", "protein.fa")
        else:
            cmd += " %s -g %s -x %s -y %s" % (self.gtf, self.option('genome').prop['path'], "cds.fa", "protein.fa")
        self.logger.info('运行gtf_to_cds，形成cds pep文件')
        command = self.add_command("gtf_to_fa2cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gtf_to_cds运行完成")
        else:
            self.set_error("gtf_to_cds运行出错!")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            if self.option("seq_type").lower() == "all":
                all_files = ['mRNA.fa', 'cds.fa', 'protein.fa']
            elif self.option("seq_type").lower() == "mRNA":
                all_files = ['mRNA.fa']
            elif self.option("seq_type").lower() == "cds":
                all_files = ['cds.fa']
            elif self.option("seq_type").lower() == "protein":
                all_files = ['protein.fa']
            for each in all_files:
                fname = os.path.basename(each)
                link = os.path.join(self.output_dir, fname)
                if os.path.exists(link):
                    os.remove(link)
                os.link(each, link)
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "extract_gff_fasta" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.extract_gff_fasta",
            "instant": True,
            "options": dict(
                genome=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                gff=test_dir + "/" + "gtf/Saccharomyces_cerevisiae.R64-1-1.39.gtf",
                seq_type='all',
                reverse=False,
                no_cds=False,
                chr="avs",
                strand="+",
                start=10,
                end=10000
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
