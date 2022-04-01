# -*- coding: utf-8 -*-
# __author__ = 'wentian.liu'
# modified 2019.04.9

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import json
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.file import getsize, exists, list_dir


class CircosAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(CircosAgent, self).__init__(parent)
        options = [
            # {"name": "snp", "type": "string"},  # snp输入文件
            # {"name": "indel", "type": "string"},  # indel输入文件
            # {"name": "cnv", "type": "string"},  # cnv输入文件
            # {"name": "sv", "type": "string"},  # sv输入文件
            # {"name": "snpplusindel", "type": "string"},  # snpplusindel输入文件
            # {"name": "ssr", "type": "string"},  # ssr输入文件
            {"name": "gff", "type": "string"},  # ref.gff
            {"name": "chrlist", "type": "string"},  # 所有样本ref.chrlist
            {"name": "color", "type": "int"},  # 颜色方案1,2,3,4,5
            {"name": "chromosome", "type": "string"},  # 选择基因组区域，按逗号分隔。
            {"name": "variant", "type": "string"},  # 数据
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "target_path", 'type': "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("gff"):
            raise OptionError("请设置gff")
        if not self.option("chrlist"):
            raise OptionError("请设置chrlist")
        if not self.option("color"):
            raise OptionError("请设置color")
        if not self.option("chromosome"):
            raise OptionError("请设置chromosome")
        if not self.option("variant"):
            raise OptionError("请设置variant")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(CircosAgent, self).end()


class CircosTool(Tool):
    def __init__(self, config):
        super(CircosTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/WGS/circos-0.69-6/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/perl/perls/perl-5.24.0/bin')
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.circos = self.config.PACKAGE_DIR + '/wgs_v2/circos.pipeline.pl'
        self.variant = json.loads(self.option("variant"))

    def make_list(self):
        """
        生成circos需要的输入文件circos.list
        :return:
        """
        write_lines = "#chromosomes\tcolors\tcircosnumber\ttype\ttypeplot\tsampleid\twindowssize\tpwd\n"
        circosnumber_list5 = ["0.95r,0.85r", "0.75r,0.65r", "0.55r,0.45r", "0.35r,0.25r", "0.20r,0.10r"]
        circosnumber_list4 = ["0.95r,0.85r", "0.75r,0.65r", "0.55r,0.45r", "0.35r,0.25r"]
        circosnumber_list3 = ["0.75r,0.65r", "0.55r,0.45r", "0.35r,0.25r"]
        circosnumber_list2 = ["0.75r,0.50r", "0.35r,0.20r"]
        circosnumber_list1 = ["0.75r,0.35r"]
        num = 0
        variant_len = len(self.variant)
        print self.variant
        self.logger.info(self.variant)
        for i in self.variant:
            type = i["variant"]
            typeplot = i["style"]
            pwd = i["pwd"]
            file_path = self.dowmload_from_s3(pwd)
            self.logger.info(file_path)
            if i["type"] == "ref" or i["type"] == "after":
                sampleid = "NA"
            elif i["analysis_object"] == "all":
                sampleid = "NA"
            else:
                sampleid = i["analysis_object"]
            windowssize = i["win_step"]
            if variant_len == 1:
                circosnumber_list = circosnumber_list1
            elif variant_len == 2:
                circosnumber_list = circosnumber_list2
            elif variant_len == 3:
                circosnumber_list = circosnumber_list3
            elif variant_len == 4:
                circosnumber_list = circosnumber_list4
            else:
                circosnumber_list = circosnumber_list5
            circosnumber = circosnumber_list[num]
            num += 1
            write_lines += "all" + "\t" + str(self.option("color")) + "\t" + str(circosnumber) + "\t" + type + "\t"\
                           + str(typeplot) + "\t" + str(sampleid) + "\t" + str(windowssize) + "\t" + str(file_path) + "\n"
        with open(os.path.join(self.work_dir, "circos.list"), "w")as fw:
            fw.write(write_lines)

    def cut_vcf(self, vcf):
        """
        将vcf切成snp的vcf和indel的vcf
        :return:
        """
        snp_write_lines = ""
        indel_write_lines = ""
        with open(vcf, "r")as fr:
            # lines = fr.readlines()
            for line in fr:
                if line.startswith("#"):
                    snp_write_lines += line
                    indel_write_lines += line
                else:
                    temp = line.strip().split("\t")
                    conbine = ",".join([temp[3], temp[4]])
                    tmp = conbine.strip().split(",")
                    line_type = "snp"
                    for i in tmp:
                        if len(i) > 1:
                            line_type = "indel"
                    if line_type == "indel":
                        indel_write_lines += line
                    if line_type == "snp":
                        snp_write_lines += line
        with open(os.path.join(self.work_dir, "snp.vcf"), "w")as fw1:
            fw1.write(snp_write_lines)
        with open(os.path.join(self.work_dir, "indel.vcf"), "w")as fw2:
            fw2.write(indel_write_lines)

    def run_circos(self):
        """
        chrlist:ref.chrlist
        gff:ref.gff
        outfile:输出文件名称
        outdir：输出路径
        paramlist：circos.list
        chrom：选中的染色体
        """
        cmd = "{} {} -chrlist {}  -gff {} -outfile {} -outdir {} -paramlist {} -chrom {}"\
            .format(self.perl_path, self.circos, self.option("chrlist"), self.option("gff"), "circos", self.output_dir,
                    os.path.join(self.work_dir, "circos.list"), self.option("chromosome"))
        # for i in self.variant:
        #     if i["variant"] == "snp" or i["variant"] == "indel":
        #         self.cut_vcf(self.option(i["variant"]))
        #         cmd += " -{} {}".format(i["variant"], os.path.join(self.work_dir, (i["variant"] + ".vcf")))
        #     elif i["variant"] == "gene":
        #         pass
        #     elif i["variant"] == "ssr" and i["type"] == "ref":
        #         cmd += " -{} {}".format(i["variant"], self.option("ssr"))
        #     else:
        #         cmd += " -{} {}".format(i["variant"], self.option(i["variant"]))
        command = self.add_command("circos", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("circos运行成功")
        else:
            self.set_error("circos运行失败")

    def dowmload_from_s3(self, pwd):
        """
        从对象存储中下载文件到指定路径
        :return:
        """
        if not pwd.startswith("/mnt"):
            if not os.path.exists(self.work_dir + "/temp"):
                os.mkdir(self.work_dir + "/temp")
            self.logger.info("开始下载对象存储中的文件！")
            transfer = MultiFileTransfer()
            if not exists(pwd):
                self.set_error("文件%s不存在！" % pwd)
            transfer.add_download(pwd, '{}/temp/'.format(self.work_dir))
            transfer.perform()
            self.logger.info("下载对象存储中的文件成功！")
            file_name_list = pwd.strip().split("/")
            file_name = file_name_list[-1]
            file_path = os.path.join(self.work_dir, ("temp/" + file_name))
            return file_path
        else:
            return pwd

    def set_db(self):
        if os.path.exists(os.path.join(self.output_dir, "circos.png")):
            self.logger.info("circos运行成功,png已产生！")
        else:
            self.set_error("png生成失败")
        print self.option('target_path')
        circos = self.api.api("wgs_v2.circos")
        self.logger.info("开始进行circos导表")
        if self.option("main_id"):
            circos.add_sg_circos(self.option("main_id"), (self.option('target_path') + "/circos.png"),
                                 (self.option('target_path') + "/circos.svg"))
        self.logger.info("circos导表成功！")

    def run(self):
        super(CircosTool, self).run()
        self.make_list()
        self.run_circos()
        self.set_db()
        self.end()
