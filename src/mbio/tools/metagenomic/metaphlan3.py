# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
# __last_modify__ = '20210304'


from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import json
from biocluster.core.exceptions import OptionError
import shutil


class Metaphlan3Agent(Agent):
    """
    metaphlan2基于reads进行物种分类以及丰度统计
    """

    def __init__(self, parent):
        super(Metaphlan3Agent, self).__init__(parent)
        options = [
            {"name": "read1", "type": "infile", "format": "sequence.fastq"},
            {"name": "read2", "type": "infile", "format": "sequence.fastq"},
            {"name": "single", "type": "infile", "format": "sequence.fastq"},
            {"name": "sample", "type": "string"},
            {"name": "read_min_len", "type": "int", "default": 70},
            {"name": "bt2_ps", "type": "string", "default": "very-sensitive"},#sensitive,very-sensitive,sensitive-local,very-sensitive-local
            # {"name": "tax_lev", "type": "string", "default": "a"},#a,k,p,c,o,f,g,s
            {"name": "min_cu_len", "type": "int", "default": 2000},
            {"name": "stat", "type": "string", "default": "avg_g"},#avg_g,avg_l,tavg_g,tavg_l,wavg_g,wavg_l,med
            {"name": "stat_q", "type": "float", "default": 0.2},  #
            {"name": "threads", "type": "int", "default": 10},
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("read1").is_set:
            raise OptionError("请输入read1文件！")
        if not self.option("read2").is_set:
            raise OptionError("请输入read2文件！")
        if not self.option("sample"):
            raise OptionError("请输入sample的样品名称！")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = self.option('threads')
        self._memory = '50G'


class Metaphlan3Tool(Tool):
    def __init__(self, config):
        super(Metaphlan3Tool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/metaphlan3/bin/"
        self.db_path = self.path +\
            '../lib/python3.7/site-packages/metaphlan/metaphlan_databases/mpa_v30_CHOCOPhlAn_201901'
        self.set_environ(PATH=self.path)
        self.out_sam = os.path.join(self.work_dir, self.option('sample') + ".sam")
        self.out = os.path.join(self.output_dir, self.option('sample') + ".taxon.xls")

    def run(self):
        super(Metaphlan3Tool, self).run()
        #self.run_bowtie2()
        self.run_metaphlan3()
        self.set_output()
        self.end()

    def run_bowtie2(self):
        cmd = "{}/bowtie2 --sam-no-hd --sam-no-sq --no-unal --sensitive -S {}" +\
            " -x {} -p {}"
        cmd = cmd.format(self.path, self.out_sam, self.db_path, self.option('threads'))
        if self.option('single').is_set:
            cmd += ' -U ' + self.option("single").path
        else:
            cmd += " -1 {} -2 {}".format(self.option("read1").path, self.option("read2").path)
        command = self.add_command('run_bowtiew2', cmd, shell=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_bowtie2运行完成")
            #os.remove('umapped.fastq')
        else:
            self.set_error("run_bowtie2运行出错")

    def run_metaphlan3(self):
        fq_path = ''
        if self.option('single').is_set:
            fq_path = self.option("single").path
        else:
            fq_path = self.option("read1").path + ',' + self.option("read2").path
        if self.option('single').is_set: fq_path = self.option("single").path
        if os.path.exists("metagenome.bowtie2.bz2"):
            os.remove("metagenome.bowtie2.bz2")
        cmd = "{}/metaphlan {} --bowtie2out metagenome.bowtie2.bz2 --nproc {} -o {} --input_type fastq -x mpa_v30_CHOCOPhlAn_201901" +\
            " --sample_id {} --min_cu_len {} --stat {} --stat_q {}"
        cmd = cmd.format(self.path, fq_path, self.option('threads'), self.out,
                         self.option("sample"), self.option("min_cu_len"),
                         self.option("stat"), self.option('stat_q'))
        command = self.add_command("run_metaphlan3", cmd, shell=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_metaphlan3运行完成！")
        else:
            self.set_error("run_metaphlan3运行完成运行出错!")

    def set_output(self):
        level_labs = {'K': 2, 'P': 3, 'C': 4, 'O': 5, 'F': 6, 'G': 7, 'S': 8}
        out_hander = {}
        for h in level_labs:
            out_f = os.path.join(self.output_dir, "{}_{}.xls".format(self.option("sample"), h))
            if os.path.exists(out_f):
                os.remove(out_f)
            out_h = open(out_f, 'w')
            out_h.write("taxonomy_id\trelative_abundance\tlinkage\n")
            out_hander[h] = out_h
        with open(self.out, 'r') as r:
            for line in r:
                if line.startswith("#"):
                    continue
                line = line.strip().split('\t')
                clades = line[0].split('|')
                level_lab = clades[-1][0].upper()
                linkage = {s[0].upper(): s[3:] for s in clades}
                taxon = line[0].split('|')[-1]
                relative_abundance = line[2]
                if level_lab in level_labs:
                    out_hander[level_lab].write("{}\t{}\t{}\n".format(taxon, relative_abundance, json.dumps(linkage)))
