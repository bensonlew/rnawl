# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import os
import re
import subprocess
import shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
from mbio.files.meta.otu.otu_table import OtuTableFile
import numpy as np


class AsvTaxonStatAgent(Agent):
    """
    metaasv
    """
    def __init__(self, parent):
        super(AsvTaxonStatAgent, self).__init__(parent)
        options = [
            {'name': 'in_otu_table', 'type': 'infile', 'format': 'meta.otu.otu_table'},  # 输入的otu表
            {'name': 'taxon_file', 'type': 'infile', 'format': 'taxon.seq_taxon'},  # 输入的taxon文件
            {'name': 'sub_otu_table', 'type': 'infile', 'format': 'meta.otu.otu_table'},  # 输入抽平后的otu_taxon表
            {'name': 'otu_taxon_biom', 'type': 'outfile', 'format': 'meta.otu.biom'},  # 输出的biom文件
            {'name': 'otu_taxon_table', 'type': 'outfile', 'format': 'meta.otu.otu_table'},  # 输出的otu表文件
            {'name': 'otu_taxon_dir', 'type': 'outfile', 'format': 'meta.otu.tax_summary_dir'}, # 输出的otu_taxon_dir(absolute)文件夹
            {'name': 'otu_taxon_dir', 'type': 'outfile', 'format': 'meta.otu.tax_summary_dir'}]  # 输出的otu_taxon_dir文件夹
        self.add_option(options)
        self.step.add_steps('OtuTaxonStat')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self.otu_taxon_otu_r = ''

    def step_start(self):
        self.step.OtuTaxonStat.start()
        self.step.update()

    def step_end(self):
        self.step.OtuTaxonStat.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        """
        if self.option("sub_otu_table").is_set:
            pass
        else:
            if not self.option("in_otu_table").is_set:
                raise OptionError("输入的OTU文件不能为空")
            if not self.option("taxon_file").is_set:
                raise OptionError("输入的taxon文件不能为空")
            self.option("in_otu_table").get_info()
            if self.option("in_otu_table").prop['metadata'] == "taxonomy":
                raise OptionError("otu表不应该有taxonomy信息")
        return True

    def end(self):
        super(AsvTaxonStatAgent, self).end()

    def set_resource(self):
        """
        设置所需要的资源
        """
        self._cpu = 2
        self._memory_increase_step = 20
        if self.option("sub_otu_table").is_set:
            self._memory = '20G'
        else:
            table_size = self.option("in_otu_table").get_size() / 1024  # 单位K/get_size是b单位
            if table_size < 1:
                self._memory = '20G'
            else:
                multiple = table_size / 1024.00 / 5.00    # 文件每增加5M，mem增加1G
                self._memory = str(20 + int(multiple)) + 'G'


class AsvTaxonStatTool(Tool):
    """
    otu taxon stat tool
    需要软件biom
    需要脚本make_otu_table.py,summarize_taxa.py,sum_tax.pl,
    """
    def __init__(self, config):
        super(AsvTaxonStatTool, self).__init__(config)
        self._version = 1.0
        self._biom_path = "miniconda2/bin/biom"
        self._make_otu_table_path = "miniconda2/bin/make_otu_table.py"
        self._summarize_taxa_path = "miniconda2/bin/summarize_taxa.py"
        self._sum_tax_path = os.path.join(Config().SOFTWARE_DIR, "bioinfo/taxon/scripts/sum_tax.fix.pl")
        self.otu_taxon_dir = os.path.join(self.work_dir, "output", "tax_summary_a")

    def split_subsample_table(self):
        """
        将抽平后的otu表拆分成in_otu_table和seq_taxon文件

        :return: 生成in_otu_table和seq_taxon文件
        """
        with open(self.option("sub_otu_table").prop['path'], 'r') as sr:
            head = sr.next().strip()
            w1 = open("in_otu_table", 'w')
            w2 = open("seq_taxon", 'w')
            w1.write(head + "\n")
            for line in sr:
                lst = re.split("\t", line.strip())
                tmp = lst[0].split(";")
                otu = tmp[-1]
                tax = ";".join(tmp[0:-1])
                w1.write(otu + "\t" + "\t".join(lst[1:]) + "\n")
                w2.write(otu + "\t" + tax + "\n")
            w1.close()
            w2.close()

    def get_biom_otu(self):
        """
        根据in_otu_table和seq_taxon文件生成biom表和otu表

        :return: 生成的biom和otu文件的路径
        """
        if self.option("sub_otu_table").is_set:
            self.split_subsample_table()
            tax_file = "seq_taxon"
            otu_file = "in_otu_table"
        else:
            tax_file = self.option("taxon_file").prop['path']
            otu_file = self.option("in_otu_table").prop['path']
        otu_tax = dict()
        with open(tax_file, 'r') as r1:
            for line in r1:
                line = re.sub("; ",";",line)
                line = re.split('\t', line.rstrip('\n'))
                tmp = re.split(';', re.sub(r'\s+', '_', line[1]))
                new_line = '; '.join(tmp)
                otu_tax[line[0]] = new_line
        taxon_otu = os.path.join(self.work_dir, "output", "asv_taxon.xls")
        otu_stat = os.path.join(self.work_dir, "output", "asv_summary.xls")
        self.logger.info("正在生成otu_taxon.xls和otu_summary.xls")
        with open(otu_file, 'r') as r2:
            with open(taxon_otu, 'w') as w:
                with open(otu_stat, 'w') as w2:
                    line1 = r2.next().rstrip('\n')
                    if re.search(r'Constructed from biom', line1):
                        line1 = r2.next().rstrip('\n')
                    line1 = line1.strip().split("\t")
                    line1[0] = "ASV ID"
                    new_line1 = "\t".join(line1)
                    w.write(new_line1 + "\t" + "taxonomy" + "\n")
                    w2.write(new_line1 + "\t" + "taxonomy" + "\n")
                    for line in r2:
                        line = line.rstrip('\n')
                        name = re.split('\t', line)[0]
                        line = re.sub(r'\.0', '', line)
                        tmp_lst = re.split('\t', line)[1:]
                        line_stat = ["1" if int(i) else "0" for i in tmp_lst]
                        line_stat = name + "\t" + "\t".join(line_stat)
                        if name in otu_tax.keys():
                            w.write(line + '\t' + otu_tax[name] + "\n")
                            w2.write(line_stat + '\t' + otu_tax[name] + "\n")
                        else:
                            add_line = "d__unclassified;k__unclassified;p__unclassified;c__unclassified;o__unclassified;f__unclassified;g__unclassified;s__unclassified"
                            w.write(line + '\t' + add_line + "\n")
                            w2.write(line_stat + '\t' + add_line + "\n")
        taxon_otu_obj = OtuTableFile()
        taxon_otu_obj.set_path(taxon_otu)
        taxon_otu_obj.get_info()
        new_taxon_otu = taxon_otu + ".new"
        taxon_otu_obj.complete_taxonomy(taxon_otu, new_taxon_otu)
        os.remove(taxon_otu)
        shutil.copy2(new_taxon_otu, taxon_otu)
        os.remove(new_taxon_otu)
        taxon_stat_obj = OtuTableFile()
        taxon_stat_obj.set_path(otu_stat)
        taxon_stat_obj.get_info()
        new_otu_stat = otu_stat + ".new"
        taxon_stat_obj.complete_taxonomy(otu_stat, new_otu_stat)
        os.remove(otu_stat)
        shutil.copy2(new_otu_stat, otu_stat)
        os.remove(new_otu_stat)

        biom = os.path.join(self.work_dir, "output", "asv_taxon.biom")
        cmd = self._biom_path + " convert -i " + taxon_otu + " -o " + biom\
            + " --process-obs-metadata taxonomy --table-type \"OTU table\" --to-hdf5"
        create_taxon_biom = self.add_command("create_taxon_biom", cmd)
        self.logger.info("由otu开始转化biom")
        create_taxon_biom.run()
        self.wait(create_taxon_biom)
        if create_taxon_biom.return_code == 0:
            self.logger.info("taxon_biom生成成功")
        else:
            self.set_error("taxon_biom生成失败")
        return(biom, taxon_otu)

    def get_diff_level(self, biom):
        """
        :param biom: biom文件路径
        """
        tax_summary_a_dir = os.path.join(self.work_dir, "output", "tax_summary_a")
        tax_summary_dir = os.path.join(self.work_dir, "output", "tax_summary_r")
        if os.path.exists(tax_summary_dir):
            shutil.rmtree(tax_summary_dir)
        if os.path.exists(tax_summary_a_dir):
            shutil.rmtree(tax_summary_a_dir)
        cmd = self._summarize_taxa_path + " -i " + biom + ' -o ' + tax_summary_a_dir\
            + " -L 1,2,3,4,5,6,7,8 -a "
        cmd2 = self._summarize_taxa_path + " -i " + biom + ' -o ' + tax_summary_dir\
            + " -L 1,2,3,4,5,6,7,8 "  # modify by zhouxuan 2016.11.29 (add 10 line)
        create_tax_summary_ = self.add_command("create_tax_summary_r", cmd2)
        self.logger.info("开始生成tax_summary文件夹")
        create_tax_summary_ .run()
        self.wait(create_tax_summary_)
        if create_tax_summary_.return_code == 0:
            self.logger.info("tax_summary文件夹生成成功")
        else:
            self.set_error("tax_summary文件夹生成失败")
            raise Exception("tax_summary文件夹生成失败")

        create_tax_summary = self.add_command("create_tax_summary", cmd)
        self.logger.info("开始生成tax_summary_a文件夹")
        create_tax_summary.run()
        self.wait(create_tax_summary)
        if create_tax_summary.return_code == 0:
            self.logger.info("文件夹生成成功")
        else:
            self.set_error("文件夹生成失败")
            raise Exception("文件夹生成失败")

        list_ = os.listdir(tax_summary_a_dir)
        for my_otu_table in list_:
            if re.search(r"txt", my_otu_table):
                my_otu_table = os.path.join(tax_summary_a_dir, my_otu_table)
                otu_basename = os.path.basename(my_otu_table)
                otu_basename = re.sub(r'\.txt$', r'.xls', otu_basename)
                otu_name = os.path.join(tax_summary_a_dir, otu_basename)
                cmd = self._sum_tax_path + " -i " + my_otu_table + " -o " + otu_name
                self.logger.info(cmd)
                try:
                    subprocess.check_call(cmd, shell=True)
                except subprocess.CalledProcessError:
                    self.set_error("运行sum_tax.pl出错")
                    raise Exception("运行sum_tax.pl出错")

        list__ = os.listdir(tax_summary_dir)
        for my_otu_table in list__:
            if re.search(r"txt", my_otu_table):
                my_otu_table = os.path.join(tax_summary_dir, my_otu_table)
                otu_basename = os.path.basename(my_otu_table)
                otu_basename = re.sub(r'\.txt$', r'.xls', otu_basename)
                otu_name = os.path.join(tax_summary_dir, otu_basename)
                cmd = self._sum_tax_path + " -i " + my_otu_table + " -o " + otu_name
                self.logger.info(cmd)
                try:
                    subprocess.check_call(cmd, shell=True)
                except subprocess.CalledProcessError:
                    self.set_error("运行sum_tax.pl出错")
                    raise Exception("运行sum_tax.pl出错")

        list_ = os.listdir(tax_summary_a_dir)
        for table in list_:
            if re.search(r"txt$", table):
                file_ = os.path.join(tax_summary_a_dir, table)
                os.remove(file_)
            if re.search(r"new$", table):
                name = re.sub(r"txt\.new$", r"full.xls", table)
                file_ = os.path.join(tax_summary_a_dir, table)
                new_file = os.path.join(tax_summary_a_dir, name)
                os.rename(file_, new_file)
        self.logger.info("tax_summary_a_dir开始整理输出文件夹")

        list__ = os.listdir(tax_summary_dir)  # modify by zhouxuan 2016.11.29 (add 11 line) 222
        for table in list__:
            if re.search(r"txt$", table):
                file_ = os.path.join(tax_summary_dir, table)
                os.remove(file_)
            if re.search(r"new$", table):
                name = re.sub(r"txt\.new$", r"full.xls", table)
                file_ = os.path.join(tax_summary_dir, table)
                new_file = os.path.join(tax_summary_dir, name)
                os.rename(file_, new_file)
        self.logger.info("tax_summary_dir开始整理输出文件夹")

        self.rename()

    def rename(self):
        """
        将各级文件重命名，将原始的OTU和biom放入到tax_summary文件夹中
        """
        level_level = {
            "L1": "Domain",
            "L2": "Kingdom",
            "L3": "Phylum",
            "L4": "Class",
            "L5": "Order",
            "L6": "Family",
            "L7": "Genus",
            "L8": "Species"
        }
        tax_summary_dir = os.path.join(self.work_dir, "output", "tax_summary_r")
        list__ = os.listdir(tax_summary_dir)
        for table in list__:
            match = re.search(r"(.+)(L\d)(.+)", table)
            prefix = match.group(1)
            suffix = match.group(3)
            level = match.group(2)
            newname = prefix + level_level[level] + ".percent" + suffix
            table = os.path.join(tax_summary_dir, table)
            newname = os.path.join(tax_summary_dir, newname)
            os.rename(table, newname)

        tax_summary_a_dir = os.path.join(self.work_dir, "output", "tax_summary_a")
        list_ = os.listdir(tax_summary_a_dir)
        for table in list_:
            match = re.search(r"(.+)(L\d)(.+)", table)
            prefix = match.group(1)
            suffix = match.group(3)
            level = match.group(2)
            newname = prefix + level_level[level] + suffix
            table = os.path.join(tax_summary_a_dir, table)
            newname = os.path.join(tax_summary_a_dir, newname)
            os.rename(table, newname)

        otu_taxon_otu = os.path.join(tax_summary_a_dir, "asv_taxon_asv.xls")  # 获得otu序列信息表(绝对丰度表)
        if self.option('sub_otu_table').is_set:
            in_table = "in_otu_table"
        else:
            in_table = self.option('in_otu_table').prop['path']
        with open(in_table, 'r') as r:
            with open(otu_taxon_otu, 'w') as w:
                line1 = r.next()
                if re.search(r'Constructed from biom', line1):
                    line1 = r.next()
                w.write(line1)
                for line in r:
                    line = line
                    line = re.sub(r'\.0', '', line)
                    w.write(line)

        otu_taxon_otu_ = os.path.join(tax_summary_dir, "asv_taxon_asv.percent.xls")  # 相对丰度表
        self.percent(otu_taxon_otu, otu_taxon_otu_)
        self.otu_taxon_otu_r = otu_taxon_otu_

        biom2 = os.path.join(tax_summary_dir, "asv_taxon_asv.percent.biom")
        cmd2 = self._biom_path + " convert -i " + otu_taxon_otu_ + " -o " + biom2\
            + " --table-type \"OTU table\" --to-hdf5"
        create_taxon_biom_otu_ = self.add_command("create_taxon_otu_biom_", cmd2)
        self.logger.info("由otu开始转化biom")
        create_taxon_biom_otu_.run()
        self.wait(create_taxon_biom_otu_)
        if create_taxon_biom_otu_.return_code == 0:
            self.logger.info("taxon_biom_otu生成成功")
        else:
            self.set_error("taxon_biom_otu生成失败")


        biom = os.path.join(tax_summary_a_dir, "asv_taxon_asv.biom")  # otu绝对丰度表的biom格式文件
        cmd = self._biom_path + " convert -i " + otu_taxon_otu + " -o " + biom\
            + " --table-type \"OTU table\" --to-hdf5"
        create_taxon_biom_otu = self.add_command("create_taxon_otu_biom", cmd)
        self.logger.info("由otu开始转化biom")
        create_taxon_biom_otu.run()
        self.wait(create_taxon_biom_otu)
        if create_taxon_biom_otu.return_code == 0:
            self.logger.info("taxon_biom_otu生成成功")
        else:
            self.set_error("taxon_biom_otu生成失败")

    def percent(self, origin_file, percent_file):
        tmp = np.loadtxt(origin_file, dtype=np.str, delimiter="\t")
        data = tmp[1:, 1:].astype(np.float)
        array_data = data / data.sum(0)
        array = np.c_[tmp[1:, 0], array_data]
        array_all = np.row_stack((tmp[0, :], array))
        s = np.char.encode(array_all, 'utf-8')
        np.savetxt(percent_file, s, fmt='%s', delimiter='\t', newline='\n')

    def get_full_otu(self, otu_table):
        """
        生成名称为全场的OTU表，放到tax_summary_dir里面
        """
        otu_full_path = os.path.join(self.work_dir, "output", "tax_summary_a", "asv_taxon_asv.full.xls")
        with open(otu_table, 'rb') as r, open(otu_full_path, 'wb') as w:
            lines = r.readlines()
            head = lines[0].rstrip("\r\n").split("\t")
            head.pop(-1)
            w.write("\t".join(head) + "\tTotal\tPercent")
            w.write("\n")
            total_num = 0
            self.total_dict = {}
            for lin in lines[1:]:
                asv_number = 0
                lin = lin.strip().split("\t")
                asv_name = lin[0]
                lin[0] = "{};{}".format(lin[-1], lin[0])
                for number in lin[1:-1]:
                    asv_number += int(number)
                total_num += asv_number
                self.total_dict[asv_name] = asv_number

            self.total_percent_dict = {}
            for line in lines[1:]:
                line = line.rstrip("\r\n").split("\t")
                asv_names = line[0]
                line[0] = "{};{}".format(line[-1], line[0])
                line.pop(-1)
                if asv_names in self.total_dict:
                    asv_num = self.total_dict[asv_names]
                else:
                    asv_num = 0
                percent = float(asv_num) / total_num
                self.total_percent_dict[asv_names] = percent
                w.write("\t".join(line) + "\t{}\t{}".format(asv_num, percent))
                w.write("\n")

        otu_full_path_r = os.path.join(self.work_dir, "output", "tax_summary_r", "asv_taxon_asv.percent.full.xls")
        self.logger.info("开始生成相对全局otu表")
        with open(otu_full_path, 'rb') as r1, open(otu_full_path_r, 'wb') as ws, open(self.otu_taxon_otu_r, 'rb') as r2:
            line = r1.readlines()
            line2 = r2.readlines()
            col = len(line)
            self.logger.info("行数为{}".format(col))
            for x in range(0, col):
                if line[x].startswith("ASV ID"):
                    name = line[x].strip().split("\t")[0:-2]
                    ws.write("\t".join(name) + "\tPercent\n")
                else:
                    f = line[x].rstrip("\r\n").split(";")
                    m = line2[x].rstrip("\r\n").split("\t")
                    m_asv_name = m[0]
                    ws.write(";".join(f[0:-1]))
                    ws.write(";" + "\t".join(m) + "\t{}\n".format(self.total_percent_dict[m_asv_name]))


    def run(self):
        """
        运行
        """
        super(AsvTaxonStatTool, self).run()
        (biom, otu_table) = self.get_biom_otu()
        self.option("otu_taxon_biom").set_path(biom)
        self.option("otu_taxon_biom").check()
        self.option("otu_taxon_table").set_path(otu_table)
        self.option("otu_taxon_table").check()
        self.get_diff_level(biom)
        self.option("otu_taxon_dir").set_path(self.otu_taxon_dir)
        self.option("otu_taxon_table").check()
        self.get_full_otu(otu_table)
        self.logger.info("otu_taxon完成，即将退出程序")
        self.end()
