# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20190222

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class PfamAgent(Agent):
    """
    参考基因组的pfam数据库注释
    """
    def __init__(self, parent):
        super(PfamAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "sample_name", "type": "string"},
            {"name": "e_seq ", "type": "float", "default": 1e-5},
            {"name": "e_dom ", "type": "float", "default": 1e-5},
            {"name": "b_seq ", "type": "float", "default": 1e-5},
            {"name": "b_dom ", "type": "float", "default": 1e-5},
        ]
        self.add_option(options)
        self.step.add_steps('analysis')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.analysis.start()
        self.step.update()

    def step_end(self):
        self.step.analysis.finish()
        self.step.update()

    def check_options(self):
        if not self.option("fasta"):
            raise OptionError("请输入fasta文件！")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 8
        self._memory = '30G'

    def end(self):
        super(PfamAgent, self).end()


class PfamTool(Tool):
    def __init__(self, config):
        super(PfamTool, self).__init__(config)
        self.set_environ(PERL5LIB=self.config.SOFTWARE_DIR + "/bioinfo/wgs_v2/PfamScan")
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/program/perl/perls/perl-5.24.0/bin")
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/wgs_v2")
        self.pfam_scan_path = "bioinfo/wgs_v2/PfamScan/pfam_scan.pl"
        self.pfam_db = self.config.SOFTWARE_DIR + "/database/pfam_33.1/"

    def script_run(self):
        """
        运行 pfam_scan
        :return:
        """
        self.outfile = os.path.join(self.output_dir, self.option('sample_name') + '.pfam')
        cmd = "{} -fasta {} -dir {} -outfile {} -cpu 8"\
            .format(self.pfam_scan_path, self.option('fasta').prop['path'], self.pfam_db, self.outfile)
        self.logger.info(cmd)
        self.logger.info("start pfam_scan")
        command = self.add_command("pfam_scan", cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("pfam_scan script finish！")
        else:
            if not self.check_error:
                self.set_error("pfam_scan出错！")
            else:
                self.logger.info("pfam_scan，正常退出！")

    def run_stat(self):
        """
        整理输出
        :return:
        """
        desc_dict = {}
        num = 0
        seq_num = 0
        desc_table_path = self.pfam_db + 'Pfam-A.clans.tsv'
        with open(self.outfile) as f, open(desc_table_path) as g, open(self.output_dir + "/" + self.option('sample_name') + ".result.xls","w") as t:
            t.write("seq id" + "\t" + "alignment start" + "\t" + "alignment end" + "\t" + "envelope start" + "\t" + "envelope end" + "\t" + "hmm acc" + "\t" + "hmm name" + "\t" + "type" + "\t" + "hmm start" + "\t" + "hmm end" + "\t" + "hmm length" + "\t" + "bitscore" + "\t" + "E-value" + "\t" + "significance" + "\t" + "clan" + "\t" + "Description"+ "\n")
            data1 = f.readlines()
            data2 = g.readlines()

            for x in data2:
                desc_dict[x.split("\t")[0]] = x.strip().split("\t")[4]
            for i in data1:
                if i.startswith("#"):
                    pass
                else:
                    if i.strip():
                        num += 1
                        if i.strip().split()[5].split(".")[0] in desc_dict.keys():
                            t.write("\t".join(i.strip().split()) + "\t" + desc_dict[i.strip().split()[5].split(".")[0]] + "\n")
                        else:
                            t.write("\t".join(i.strip().split()) + "\t" + "-" + "\n")
        with open(self.option('fasta').prop['path']) as v, open(self.output_dir + "/stat","w") as j:
            data3 = v.readlines()
            for i in data3:
                if i.startswith(">"):
                    seq_num += 1
            j.write(self.option('sample_name') + '\t' + str(seq_num) + '\t' + str(num))
        os.remove(self.outfile)

    def run(self):
        super(PfamTool, self).run()
        self.script_run()
        self.run_stat()
        self.end()

    def check_error(self, file_path):
        with open(file_path, 'r') as r:
            for line in r:
                if re.match('Unrecognised character \[\*\] in.*', line):
                    return True
        return False
