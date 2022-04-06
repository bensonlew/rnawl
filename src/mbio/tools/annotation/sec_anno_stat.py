# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2018/11/15'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os, re
from biocluster.core.exceptions import OptionError
import pandas as pd
from mbio.packages.metabolome.common import link_file
from mbio.packages.taxon.mask_taxon import mask_taxon


class SecAnnoStatAgent(Agent):

    def __init__(self, parent):
        super(SecAnnoStatAgent, self).__init__(parent)
        options = [
            {"name": "predict_dir", "type": "string"},
            {"name": "gene_nr_anno", "type": "infile", "format": "sequence.profile_table"},
            {"name": "nr_type", "type": "string", "default": "besthit"},  # besthit, lca, deunclassify
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("predict_dir"):
            raise OptionError("must input predict_dir", code="31204501")
        elif not os.path.isdir(self.option("predict_dir")):
            raise OptionError("predict_dir is not correct directory path", code="31204502")
        if not self.option("gene_nr_anno").is_set:
            raise OptionError("must input gene_nr_anno", code="31204503")
        if not self.option("reads_profile_table").is_set:
            raise OptionError("must input reads_profile_table", code="31204504")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = "10G"


class SecAnnoStatTool(Tool):
    def __init__(self, config):
        super(SecAnnoStatTool, self).__init__(config)
        self._version = "1.0"
        self.sec_stat_path = self.config.PACKAGE_DIR + "/metagenomic/sec_stat.py "
        self.python_path = 'miniconda2/bin/python'
        self.r_path = 'program/R-3.3.1/bin/Rscript'
        self.type_list = []
        self.name_to_name = {}

    def run_all(self):
        """
        description
        :return:
        """
        gram_pos = gram_neg = euk = ttss = None
        for file in os.listdir(self.option('predict_dir')):
            if "Gram+" in file:
                gram_pos = os.path.join(self.option('predict_dir'), file)
            elif "Gram-" in file:
                gram_neg = os.path.join(self.option('predict_dir'), file)
            elif "Euk" in file:
                euk = os.path.join(self.option('predict_dir'), file)
            elif "ttss" in file:
                ttss = os.path.join(self.option('predict_dir'), file)  # 是否提供一個文件名稱前綴？
        if gram_pos:
            self.stat("grampos", gram_pos)
            self.type_list.append("grampos")
        if gram_neg:
            self.stat("gramneg", gram_neg)
            self.type_list.append("gramneg")
        if gram_pos and gram_neg:
            self.stat("bac", gram_pos + ',' + gram_neg)
            self.type_list.append("bac")
        if euk:
            self.stat("euk", euk)
            self.type_list.append("euk")
        if gram_pos and gram_neg and euk:
            self.stat("all", gram_pos + ',' + gram_neg + ',' + euk)
            self.type_list.append("all")
        if ttss:
            self.stat("ttss_" + self.option("nr_type"), ttss)
            self.type_list.append("ttss_" + self.option("nr_type"))

    def dashrepl(self, matchobj):
        """
        add func by guhaidong 20171031
        """
        return self.name_to_name[matchobj.groups()[0]]

    def add_taxon(self, old_result, taxon_result):
        """
        add func by guhaidong 20171031
        description: 将旧注释的名称，根据词典替换成新注释名称
        """
        with open(old_result, "r") as f, open(taxon_result, "w") as w:
            # w.write(old_result)
            for i in f.readlines():
                #line = i.strip()
                new_line = re.sub(r"(name\d+)", self.dashrepl, i)
                w.write(new_line)

    def stat(self, mytype, table):
        mytype = str(mytype)
        out = os.path.join(self.work_dir, mytype)
        cmd = '{} {} -p {} -t {} -l Species -s {} -o {}'.format(self.python_path, self.sec_stat_path, self.option('reads_profile_table').path,
                                         self.option('gene_nr_anno').path, table, out)
        cmd += ' -e'
        command = self.add_command(mytype, cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("%s統計成功" % mytype)
        elif command.return_code == -11:
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("%s stat error" , variables=( mytype), code="31204501")
        self.name_to_name = mask_taxon(out + "/fisher.input", out + "/fisher_tmp_input")
        os.rename(out + "/fisher.input", out + "/fisher_input_old")
        os.rename(out + "/fisher_tmp_input", out + "/fisher.input")
        new_r_script = os.path.join(self.work_dir, "%s_run_fisher_test.r" % mytype)
        if os.path.exists(new_r_script):
            os.remove(new_r_script)
        os.link(os.path.join(self.work_dir, "run_fisher_test.r"), new_r_script)
        cmd = "%s %s_run_fisher_test.r" % (self.r_path, mytype)
        command_name = "fisher_%s_cmd" % mytype
        command2 = self.add_command(command_name, cmd).run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("%s success" % command_name)
        else:
            self.set_error("%s failed" , variables=( command_name), code="31204502")

    def mv_child(self, type):
        old_summary = os.path.join(self.work_dir, type, "summary.txt")
        new_summary = os.path.join(self.output_dir, type + "_summary.txt")
        link_file(old_summary, new_summary)
        fisher_input = os.path.join(self.work_dir, type, "fisher_input_old")
        fisher_output = os.path.join(self.work_dir, type, "fisher.output")
        fisher_output_with_taxon = os.path.join(self.work_dir, type, "fisher.output_taxon")
        self.add_taxon(fisher_output, fisher_output_with_taxon)
        new_fisher_output = os.path.join(self.output_dir, type + "_fisher.txt")
        data1 = pd.read_table(fisher_input, index_col=0)
        data2 = pd.read_table(fisher_output_with_taxon, index_col=0)
        data3 = pd.merge(data1,data2[["odds_ratio", "pvalue", "corrected_pvalue"]], left_index=True, right_index=True).sort_values("odds_ratio",ascending=False)
        data3.to_csv(new_fisher_output, sep="\t", index=True)

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        for i in self.type_list:
            self.mv_child(i)

    def run(self):
        super(SecAnnoStatTool, self).run()
        self.run_all()
        self.set_output()
        self.end()
