# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

"""axe-demux barcode拆分工具 """
import re
import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class AxeDemuxAgent(Agent):
    """
    axe-demux
    """
    def __init__(self, parent=None):
        super(AxeDemuxAgent, self).__init__(parent)
        options = [
            {'name': 'fq1', 'type': "infile", "format": "sequence.fastq"},
            {'name': 'fq2', 'type': "infile", "format": "sequence.fastq"},
            {'name': 'ziplevel', 'type': "string", "default": "6"},  # 压缩级别
            {'name': 'combinatorial', "type": "string", "default": "2"},  # Use combinatorial barcode matching
            {'name': 'mismatch', "type": "string", "default": "0"},  # mismatch
            {'name': "library_info", "type": "infile", "format": "sequence.barcode_info"},  # 文库信息，包含文库中样本酶的信息，文库类型等
        ]
        self.add_option(options)
        self.step.add_steps("axe_demux")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.axe_demux.start()
        self.step.update()

    def stepfinish(self):
        self.step.axe_demux.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        """
        if not self.option('fq1'):
            raise OptionError('必须输入1端序列')
        if not self.option('fq2'):
            raise OptionError('必须输入2端序列')
        if not self.option('library_info'):
            raise OptionError('必须输入文库信息表')
        return True

    def set_resource(self):
        """
        设置所需要的资源
        """
        self._cpu = 1
        self._memory = '2G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(AxeDemuxAgent, self).end()


class AxeDemuxTool(Tool):
    """
    """
    def __init__(self, config):
        super(AxeDemuxTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/library/gsl23/lib')
        self.axe_demux_path = 'bioinfo/seq/axe-demux'
        self.lib_type = ''
        self.fq_name = ''

    def run(self):
        super(AxeDemuxTool, self).run()
        self.get_barcode()
        self.run_axe_demux()
        self.set_output()
        self.end()

    def get_barcode(self):
        """
        根据输入的酶的信息，获取对应的barcode并建立文件
        :return:
        """
        # enzyme_dict = {"10": "GGCTAC", "16": "ACCTCT", "18": "ATAATC", "46": "TGCTTA", "52": "TTCCAC", "53": "CTCC",
        #                "61": "AGGC", "62": "GATC", "63": "TCAC", "66": "TCACC", "78": "CATCT", "85": "TCGTT",
        #                "96": "ACGTGTT", "100": "CGCTGAT", "101": "CGGTAGA", "104": "TAGCGGA", "108": "ACGACTAC",
        #                "111": "TGCAAGGA", "114": "CCGGATAT", "119": "CCATGGGT", "121": "CGTGTGGT", "122": "GCTGTGGA",
        #                "123": "GGATTGGT", "124": "GTGAGGGT", "M1": "AACG", "M3": "TTACA", "M6": "CCTCAG",
        #                "M9": "GAAGATCC", "T1": "AGAACATA", "T5": "TTCAATG", "T7": "GAACTA", "T12": "CCTA"}
        enzyme_dict = {"10": "GGCTAC", "16": "ACCTCT", "18": "ATAATC", "46": "TGCTTA", "52": "TTCCAC", "53": "CTCC",
                       "61": "AGGC", "62": "GATC", "63": "TCAC", "66": "TCACC", "78": "CATCT", "85": "TCGTT",
                       "96": "ACGTGTT", "100": "CGCTGAT", "101": "CGGTAGA", "104": "TAGCGGA", "108": "ACGACTAC",
                       "111": "TGCAAGGA", "114": "CCGGATAT", "119": "CCATGGGT", "121": "CGTGTGGT", "122": "GCTGTGGA",
                       "123": "GGATTGGT", "124": "GTGAGGGT", "M1": "AACG", "M3": "TTACA", "M6": "CCTCAG",
                       "M9": "GAAGATCC", "T1": "AGAACATA", "T5": "TTCAATG", "T7": "GAACTA", "T12": "CCTA",
                       "PstI124": "GTGAGGGT", "PstI123": "GGATTGGT", "PstI122": "GCTGTGGA", "PstI121": "CGTGTGGT",
                       "PstI119": "CCATGGGT", "PstI114": "CCGGATAT", "PstI111": "TGCAAGGA", "PstI108": "ACGACTAC",
                       "PstI104": "TAGCGGA", "PstI101": "CGGTAGA", "PstI100": "CGCTGAT", "PstI96": "ACGTGTT",
                       "PstI85": "TCGTT", "PstI78": "CATCT", "PstI62": "GATC", "PstI61": "AGGC", "PstI53": "CTCC",
                       "PstI52": "TTCCAC", "PstI46": "TGCTTA", "PstI18": "ATAATC", "PstI16": "ACCTCT", "PstI10": "GGCTAC",
                       "EcoRI124": "GTGAGGGT", "EcoRI123": "GGATTGGT", "EcoRI122": "GCTGTGGA", "EcoRI121": "CGTGTGGT",
                       "EcoRI119": "CCATGGGT", "EcoRI114": "CCGGATAT", "EcoRI111": "TGCAAGGA", "EcoRI108": "ACGACTAC",
                       "EcoRI104": "TAGCGGA", "EcoRI101": "CGGTAGA", "EcoRI100": "CGCTGAT", "EcoRI96": "ACGTGTT",
                       "EcoRI85": "TCGTT", "EcoRI78": "CATCT", "EcoRI62": "GATC", "EcoRI61": "AGGC", "EcoRI53": "CTCC",
                       "EcoRI52": "TTCCAC", "EcoRI46": "TGCTTA", "EcoRI18": "ATAATC", "EcoRI16": "ACCTCT", "EcoRI10": "GGCTAC",
                       "TaqI124": "GTGAGGGT", "TaqI123": "GGATTGGT", "TaqI122": "GCTGTGGA", "TaqI121": "CGTGTGGT",
                       "TaqI119": "CCATGGGT", "TaqI114": "CCGGATAT", "TaqI111": "TGCAAGGA", "TaqI108": "ACGACTAC",
                       "TaqI104": "TAGCGGA", "TaqI101": "CGGTAGA", "TaqI100": "CGCTGAT", "TaqI96": "ACGTGTT", "TaqI85": "TCGTT",
                       "TaqI78": "CATCT", "TaqI62": "GATC", "TaqI61": "AGGC", "TaqI53": "CTCC", "TaqI52": "TTCCAC", "TaqI46": "TGCTTA",
                       "TaqI18": "ATAATC", "TaqI16": "ACCTCT", "TaqI10": "GGCTAC", "MseI124": "GTGAGGGT", "MseI123": "GGATTGGT",
                       "MseI122": "GCTGTGGA", "MseI121": "CGTGTGGT", "MseI119": "CCATGGGT", "MseI114": "CCGGATAT", "MseI111": "TGCAAGGA",
                       "MseI108": "ACGACTAC", "MseI104": "TAGCGGA", "MseI101": "CGGTAGA", "MseI100": "CGCTGAT", "MseI96": "ACGTGTT",
                       "MseI85": "TCGTT", "MseI78": "CATCT", "MseI62": "GATC", "MseI61": "AGGC", "MseI53": "CTCC", "MseI52": "TTCCAC",
                       "MseI46": "TGCTTA", "MseI18": "ATAATC", "MseI16": "ACCTCT", "MseI10": "GGCTAC"}
        lib_type = set()
        with open(self.option("library_info").prop["path"])as fr, open(self.work_dir + '/barcode.txt', "w+")as fw:
            lines = fr.readlines()
            for line in lines[1:]:
                tmp = line.strip().split("\t")
                if re.search("RAD", tmp[4]):
                    self.lib_type = "RAD"
                elif re.search("GBS", tmp[4]):
                    self.lib_type = "GBS"
                else:
                    self.set_error("文库：{}类型：{}有问题，请检查".format(tmp[3], tmp[4]))
                barcode_dict = {}
                self.fq_name = tmp[1] + ':' + tmp[3]
                if tmp[7] in enzyme_dict.keys():
                    barcode_dict[tmp[7]] = enzyme_dict[tmp[7]]
                if tmp[8] in enzyme_dict.keys():
                    barcode_dict[tmp[8]] = enzyme_dict[tmp[8]]
                enzyme = tmp[7]
                # for key in enzyme_dict.keys():
                #     m = re.match(r'[A-Z]*[a-z]*[A-Z]*([0-9]+)[a-z]*[A-Z]*', enzyme)
                #     if str(m.group(1)) == key:
                #         barcode_dict[key] = enzyme_dict[key]
                if len(barcode_dict) == 0:
                    self.set_error("样本{}酶的信息有误，请核实".format(tmp[5]))
                elif len(barcode_dict) == 1:
                    if self.lib_type == 'GBS':
                        self.set_error("样本{}酶的信息与文库类型{}不对应，请核实".format(tmp[5], self.lib_type ))
                    fw.write(str(barcode_dict.values()[0]) + "\t" + str(tmp[5]) + "\n")  # 如果是RAD，单barcode
                elif len(barcode_dict) == 2:
                    if self.lib_type == 'RAD':
                        self.set_error("样本{}酶的信息与文库类型{}不对应，请核实".format(tmp[5], self.lib_type ))
                    else:
                        barcode1 = str(barcode_dict.keys()[0])  # 如果是GBS，双barcode，需要根据酶对barcode排序
                        barcode2 = str(barcode_dict.keys()[1])
                        info1 = enzyme.find(barcode1)
                        info2 = enzyme.find(barcode2)
                        if info1 > info2:
                            fw.write(str(barcode_dict.values()[0]) + "\t" + str(
                                barcode_dict.values()[1]) + "\t" + str(tmp[5]) + "\n")
                        else:
                            fw.write(str(barcode_dict.values()[1]) + "\t" + str(
                                barcode_dict.values()[0]) + "\t" + str(tmp[5]) + "\n")
                else:
                    self.set_error("酶的信息{}有误，请核实".format(enzyme))

    def run_axe_demux(self):
        """
        运行axe-demux
        """
        if self.lib_type == 'RAD':
            cmd = self.axe_demux_path + ' -m %s -z %s -b %s -t %s -f %s -r %s -F %s -R %s' % (
                self.option("mismatch"), self.option('ziplevel'), self.work_dir + '/barcode.txt', self.work_dir + '/' + self.fq_name + '.split.stat',
                self.option('fq1').prop['path'], self.option('fq2').prop['path'], self.work_dir + '/' + self.fq_name + ":", self.work_dir + '/' + self.fq_name + ":")
        else:
            cmd = self.axe_demux_path + ' -m %s -z %s -c %s -b %s -t %s -f %s -r %s -F %s -R %s' % (
                self.option("mismatch"), self.option('ziplevel'), self.option('combinatorial'), self.work_dir + '/barcode.txt', self.work_dir + '/' + self.fq_name + '.split.stat',
                self.option('fq1').prop['path'], self.option('fq2').prop['path'], self.work_dir + '/' + self.fq_name + ":",
                self.work_dir + '/' + self.fq_name + ":")
        self.logger.info('运行axe-demux，进行二次拆分')
        command = self.add_command("axe-demux_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("axe-demux运行完成")
        else:
            self.set_error("axe-demux运行出错!")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            all_files = os .listdir(self.work_dir)
            for files in all_files:
                if files.endswith(".fastq.gz"):
                    if re.search(r'unknown_R', files):
                        pass
                    else:
                        os.link(self.work_dir + "/" + files, self.output_dir + '/' + files)
            self.logger.info("设置axe-demux分析结果目录成功")

        except Exception as e:
            self.logger.info("设置axe-demux分析结果目录失败{}".format(e))
            self.set_error("设置axe-demux分析结果目录失败{}".format(e))
