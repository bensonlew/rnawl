# -*- coding: utf-8 -*-
# __author__ = 'hongyu.chen'
# last modify 20180814

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import subprocess
import os
import time
import random


class DpCorrectionAgent(Agent):
    """
    MS降噪
    """
    def __init__(self, parent):
        super(DpCorrectionAgent, self).__init__(parent)
        options = [
            {"name": "mom_tab", "type": "string"},
            {"name": "preg_tab", "type": "string"},
            {"name": "new_id", "type": "string"}        # example: WQ182579-S-T1-M
        ]
        self.add_option(options)
        self.step.add_steps('dp_correction')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self.queue = "gada"

    def step_start(self):
        self.step.dp_correction.start()
        self.step.update()

    def step_end(self):
        self.step.dp_correction.finish()
        self.step.update()

    def check_options(self):
        if not self.option("mom_tab"):
            raise OptionError("缺少mom_tab参数")
        if not self.option("preg_tab"):
            raise OptionError("缺少preg_tab参数")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        super(DpCorrectionAgent, self).end()


class DpCorrectionTool(Tool):
    def __init__(self, config):
        super(DpCorrectionTool, self).__init__(config)
        self.mom_tab = self.option("mom_tab")
        self.preg_tab = self.option("preg_tab")
        self.new_id = self.option("new_id")
        self.pos = self.config.SOFTWARE_DIR + "/bioinfo/medical/pos.txt"
        self.ref_all = self.config.SOFTWARE_DIR + "/database/human/pt_ref/tab_all"

    def add_gene(self, tab_file):
        dic = {}
        with open(tab_file) as FA:
            for line in FA:
                if len(line) < 8:  # 过滤掉列数小于8的位点
                    continue
                else:
                    line = line.strip()
                    tmp = line.split("\t")
                    sample_name, chrom, pos, ref, alt, dp, ref_dp, alt_dp = tmp
                    if int(ref_dp) == 0 and int(alt_dp) == 0:
                        continue
                    rf = float(ref_dp) / (float(ref_dp) + float(alt_dp))
                    if rf >= 0.8:
                        tmp.append("0/0")
                    elif rf <= 0.2:
                        tmp.append("1/1")
                    else:
                        tmp.append("0/1")
                    sample_name, chrom, pos, ref, alt, dp, ref_dp, alt_dp, genetype = tmp
                    dic[pos] = sample_name + '\t' + ref + '\t' + alt + '\t' + dp + '\t' + ref_dp + '\t' + alt_dp + '\t' + genetype
        return dic

    def merge_family(self):
        dic_M = self.add_gene(self.mom_tab)
        dic_S = {}
        output_file = self.new_id + ".merge_family.xls"
        with open(self.preg_tab, "r") as FB:
            for line in FB:
                line = line.strip()
                line = line.split('\t')
                sample_name, chrom, pos, ref, alt, dp, ref_dp, alt_dp = line
                if int(ref_dp) == 0 and int(alt_dp) == 0:
                    continue
                dic_S[pos] = sample_name + '\t' + chrom + '\t' + pos + '\t' + ref + '\t' + alt + '\t' + dp + '\t' + ref_dp + '\t' + alt_dp

        with open(self.pos, 'r') as FC, open(output_file, "w") as f:
            for key in FC:
                key = key.strip()
                M_value = dic_M.get(str(key), 'NA')
                S_value = dic_S.get(str(key), 'NA')
                content = S_value + '\t' + M_value + '\n'
                f.write(content)
            f.close()
        subprocess.call("sed -i '/NA/d' {}".format(output_file), shell=True)

    def dp_correction_run(self):
        self.merge_family()
        data = pd.read_table(self.new_id + ".merge_family.xls", header=None, sep='\t')
        data.columns = ['S_id', 'chro', 'pos', 'S_ref', 'S_alt', 'S_dp', 'S_ref_dp', 'S_alt_dp', 'M_id', 'M_ref',
                        'M_alt', 'M_dp', 'M_ref_dp', 'M_alt_dp', 'M_gene']
        dat = data.copy()
        for index in range(0, len(dat) - 1):
            M_rf = float(dat.loc[index, 'M_alt_dp']) / float(dat.loc[index, 'M_alt_dp'] + dat.loc[index, 'M_ref_dp'])
            S_rf = float(dat.loc[index, 'S_alt_dp']) / float(dat.loc[index, 'S_alt_dp'] + dat.loc[index, 'S_ref_dp'])
            if abs(M_rf - S_rf) > 0.5 or S_rf > 0.995 or S_rf < 0.005:
                continue
            else:
                if dat.loc[index, 'M_gene'] == "0/0":  # 母本0/0  父本1/1 or 0/1是,背噪值以母本1-rf为准
                    BZ = float(dat.loc[index, 'M_alt_dp']) / float(
                        dat.loc[index, 'M_alt_dp'] + dat.loc[index, 'M_ref_dp'])
                    preg_alt_dp = dat.loc[index, 'S_alt_dp'] - BZ * (
                            dat.loc[index, 'S_alt_dp'] + dat.loc[index, 'S_ref_dp'])
                    if preg_alt_dp < 0:
                        dat.loc[index, 'S_alt_dp'] = 0
                    else:
                        dat.loc[index, 'S_alt_dp'] = int(round(preg_alt_dp))

                if dat.loc[index, 'M_gene'] == "1/1":
                    BZ = float(dat.loc[index, 'M_ref_dp']) / float(
                        dat.loc[index, 'M_ref_dp'] + dat.loc[index, 'M_alt_dp'])
                    preg_ref_dp = dat.loc[index, 'S_ref_dp'] - BZ * (
                            dat.loc[index, 'S_ref_dp'] + dat.loc[index, 'S_alt_dp'])
                    if preg_ref_dp < 0:
                        dat.loc[index, 'S_ref_dp'] = 0
                    else:
                        dat.loc[index, 'S_ref_dp'] = int(round(preg_ref_dp))
        new_S_dat = dat.iloc[:, 0:8]
        new_S_dat["S_id"] = self.new_id
        new_S_dat.to_csv(self.new_id + '.mem.sort.hit.vcf.tab', sep='\t', index=False, header=None)

    def set_output(self):
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        file0 = '/{}.mem.sort.hit.vcf.tab'.format(self.new_id)
        os.link(self.work_dir + file0, self.output_dir + "/{}.tab".format(self.new_id))
        self.logger.info("设置文件夹路径成功")
        self.logger.info("开始link样本到tab_all文件夹！")
        # if self.config.RGW_ENABLE:
        #     self.upload_to_s3("output/*", "s3://app/database/human/pt_ref/tab_all/")
        # else:
        time.sleep(random.randint(3, 25))   # 随机停避免多个tool同时复制胎儿矫正后tab造成的文件冲突20191119 by hd
        if not os.path.exists(self.ref_all + "/{}.tab".format(self.new_id)):
            self.logger.info("移动时间：{}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
            os.link(self.output_dir + "/{}.tab".format(self.new_id), self.ref_all + "/{}.tab".format(self.new_id))
            self.logger.info("在tab_all中备份该样本：{}".format(self.new_id))
        else:
            self.logger.info("tab_all中已有该样本：{}".format(self.new_id))

    def run(self):
        super(DpCorrectionTool, self).run()
        self.dp_correction_run()
        self.set_output()
        self.api.api('medical.paternity_test_v3.paternity_test_v3').add_tab(
            self.output_dir + "/{}.tab".format(self.new_id), "sg_sample_tab_ms")
        self.api.api('medical.paternity_test_v3.paternity_test_v3').update_sample_tab(
            self.new_id, "", "555555555555555555555555", "wqcf")
        self.end()
