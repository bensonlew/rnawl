# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os
import re
import shutil
import gevent
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir,link_file


class FramebotModule(Module):
    """
    metaasv heatmap
    """
    def __init__(self, work_id):
        super(FramebotModule, self).__init__(work_id)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "database", "type": "string"},  # 数据库选择
            {"name": "acid_length", "type": "int", "default": 80},  # 氨基酸长度阈值
            {"name": "seq_identity", "type": "float", "default": 0.4},
            {"name": "ref_acid", "type": "infile", "format": "sequence.fasta"},  # 参考氨基酸fasta文件
            {"name": "nucl_fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "prot_fasta", "type": "infile", "format": "sequence.fasta"},
            {'name': 'framebot_nucl_fasta', 'type': 'outfile', 'format': 'sequence.fasta'},  # 输出核酸序列
            {'name': 'framebot_prot_fasta', 'type': 'outfile', 'format': 'sequence.fasta'},  # 输出核酸序列
            {"name": "path", "type": "string"},  # 降噪结果路径
            {"name": "taxon_file", "type": "outfile", "format": "taxon.seq_taxon"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option('fasta').is_set:
            raise OptionError('必须提供输入fasta序列')
        if self.option("database") == "custom_mode":
            if not self.option('ref_acid').is_set:
                raise OptionError('必须提供输入核酸序列')
        else:
            if self.option("database") not in ['fgr/amoA_archaea_202012', 'fgr/amoA_bacteria_202012',
                                               'fgr/amoA_AOB_like_202012', 'fgr/amoA_comammox_202012',
                                               'fgr/nosZ_202012',
                                               'fgr/nosZ_atypical_1_202012', 'fgr/nosZ_atypical_2_202012',
                                               'fgr/nirK_202012',
                                               'fgr/nirS_202012', 'fgr/mcrA_202012', 'fgr/nifH_202012',
                                               'fgr/pmoA_202012',
                                               'fgr/mmoX_202012']:
                raise OptionError("数据库%s不被支持", variables=(self.option("database")), code="12700104")

    def run_framebot(self):
        self.framebot = self.add_tool("meta.framebot")
        opts = {
            "fasta": self.option("fasta"),
            "database": self.option("database"),
            "acid_length": self.option("acid_length"),
            "seq_identity": self.option("seq_identity"),
            "ref_acid": self.option("ref_acid"),
        }
        self.framebot.set_options(opts)
        self.framebot.on("end", self.get_result)
        self.framebot.run()

    def get_result(self):
        if os.path.getsize(self.framebot.output_dir + '/framebot_nucl.fasta') == 0:
            self.set_error("FrameBot矫正后序列数为0")
        nucl_fasta = os.path.join(self.output_dir, "framebot_nucl.fasta")
        prot_fasta = os.path.join(self.output_dir, "framebot_prot.fasta")
        for i in nucl_fasta, prot_fasta:
            if os.path.exists(i):
                os.remove(i)
        os.link(self.framebot.output_dir + '/framebot_nucl.fasta', nucl_fasta)
        os.link(self.framebot.output_dir + '/framebot_prot.fasta', prot_fasta)
        md5file = self.option('path') + "/ASV_md5.xls"
        tablefile = self.option('path') + "/ASV_table.xls"
        if os.path.getsize(nucl_fasta) > 0:
            asv_chang_file = self.output_dir + "/ASV_change.xls"
            asv_fasta_file = self.output_dir + "/ASV_reps.fasta"
            new_md5file = self.output_dir + "/ASV_md5.xls"
            new_tablefile = self.output_dir + "/ASV_table.xls"
            with open(asv_chang_file, "w") as t, open(nucl_fasta, "r") as f, open(tablefile, "r") as v, open(md5file, "r") as m, \
                    open(asv_fasta_file, "w") as fasta, open(new_md5file, "w") as md5, open(new_tablefile, "w") as table:
                data1 = f.readlines()
                data2 = v.readlines()
                data4 = m.readlines()
                asv_list = []
                asv_dict = {}  # raw_asvid : new_asvid
                num = 0
                for i in data1:
                    if i.startswith(">"):
                        asv_list.append(i.strip().lstrip(">"))
                for x in asv_list:
                    num += 1
                    asv_dict[x] = "ASV" + str(num)
                t.write("raw_asvid\tnew_asvid\t")
                for n in asv_dict:
                    t.write(n + "\t" + asv_dict[n] + "\n")
                for num in range(len(data1)):
                    if data1[num].startswith(">"):
                        if data1[num].strip().lstrip(">") in asv_dict:
                            fasta.write(data1[num] + data1[num + 1])
                md5.write("ASV ID\tmd5\n")
                for line2 in data4:
                    if line2.split("\t")[0] in asv_dict:
                        md5.write(line2)
                table.write(data2[0])
                for line3 in data2:
                    if line3.split("\t")[0] in asv_dict:
                        table.write(line3)
            # 转换的表是用md5 值作为 ID 的
            if os.path.exists(self.work_dir + "/qza"):
                shutil.rmtree(self.work_dir + "/qza")
            os.mkdir(self.work_dir + "/qza")
            md5_table = self.work_dir + "/qza/ASV_table.xls"
            md5_fasta = self.work_dir + "/qza/ASV_reps.fasta"
            md5_dict = {}
            with open(md5_fasta, "w") as mf, open(md5_table, "w") as mt, open(asv_fasta_file) as f1, open(new_md5file) as f2, open(new_tablefile) as f3:
                data5 = f1.readlines()
                data6 = f2.readlines()
                data7 = f3.readlines()
                for xx in data6[1:]:
                    md5_dict[xx.strip().split("\t")[0]] = xx.strip().split("\t")[1]
                for xxx in data5:
                    if xxx.startswith(">"):
                        mf.write(">" + md5_dict[xxx.strip().split(">")[1]] + "\n")
                    else:
                        mf.write(xxx)
                mt.write(data7[0])
                for xxxx in data7[1:]:
                    mt.write(md5_dict[xxxx.strip().split("\t")[0]] + "\t".join(xxxx.strip().split("\t")[1:]) + "\n")
        self.run_convert_qza()

    def run_convert_qza(self):
        """
        转table转为qza文件
        :return:
        """
        self.file_list = []
        allfiles = ["ASV_table.xls", "ASV_reps.fasta"]
        for file in allfiles:
            file_path = os.path.join(self.work_dir+ "/qza", file)
            self.convert_qza = self.add_tool("metaasv.file_to_qzaqzv")
            if re.search(r'ASV_table', file):
                options = {
                    "input_file": file_path,
                    "type": "table",
                    "prefix": "ASV_table"
                }
            else:
                options = {
                    "input_file": file_path,
                    "type": "fasta",
                    "prefix": "ASV_reps"
                }
            self.convert_qza.set_options(options)
            self.file_list.append(self.convert_qza)
        if len(self.file_list) > 1:
            self.on_rely(self.file_list, self.set_output)
        else:
            self.file_list[0].on("end", self.set_output)
        for tool in self.file_list:
            tool.run()
            gevent.sleep(0)

    def run(self):
        super(FramebotModule, self).run()
        self.run_framebot()

    def set_output(self):
        """
        生成结果文件目录
        :return:
        """
        for tool in self.file_list:
            for file in os.listdir(tool.output_dir):
                file_path = os.path.join(tool.output_dir, file)
                out_path = os.path.join(self.output_dir, file)
                link_file(file_path, out_path)
        self.option('framebot_nucl_fasta', self.output_dir + "/ASV_reps.fasta")
        self.option('framebot_prot_fasta', self.output_dir + '/framebot_prot.fasta')
        self.end()

    def end(self):
        super(FramebotModule, self).end()
