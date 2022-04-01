# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

"""多样性去接头 """
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError

import os


class GetvaildBybarcodeAgent(Agent):
    """
    GetvaildBybarcode
    """
    def __init__(self, parent=None):
        super(GetvaildBybarcodeAgent, self).__init__(parent)
        options = [
            {'name': 'fq1', 'type': "infile", "format": "sequence.fastq"},  # 一次拆分后的文库序列
            {'name': 'fq2', 'type': "infile", "format": "sequence.fastq"},
            {'name': "barcode_info", "type": "infile", "format": "sequence.barcode_info"},  # 样本总信息
            {'name': 'out_fq1', 'type': "outfile", "format": "sequence.fastq"},
            {'name': 'out_fq2', 'type': "outfile", "format": "sequence.fastq"},  # 去嵌合体之后的fq序列
            {'name': 'lib_name', 'type': "string"},  # 文库名称
            # {'name': "barcode_file", "type": "outfile", "format": "sequence.barcode_info"},  # barcode信息
            # {'name': "primer_file", "type": "outfile", "format": "sequence.barcode_info"},  # primer信息
        ]
        self.add_option(options)
        self.step.add_steps("getvaild_bybarcode")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.getvaild_bybarcode.start()
        self.step.update()

    def stepfinish(self):
        self.step.getvaild_bybarcode.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        """
        if not self.option('fq1'):
            raise OptionError('必须输入1端序列')
        if not self.option('fq2'):
            raise OptionError('必须输入2端序列')
        if not self.option('barcode_info'):
            raise OptionError('必须输入barcode信息表')
        return True

    def set_resource(self):
        """
        设置所需要的资源
        """
        self._cpu = 2
        self._memory = '2G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(GetvaildBybarcodeAgent, self).end()


class GetvaildBybarcodeTool(Tool):
    """
    """
    def __init__(self, config):
        super(GetvaildBybarcodeTool, self).__init__(config)
        # self.barcode_list = self.config.SOFTWARE_DIR + '/database/datasplit/barcode.list'
        # self.primer_list = self.config.SOFTWARE_DIR + '/database/datasplit/primer.list'
        self.getvaild_bybarcode_path = '/bioinfo/seq/scripts/getValid_rawSeq_byBarcode.2.pl '
        # self.lib_name = ''  # 文库名
        # self.barcode_path = ''  # barcode信息表
        # self.primer_path = ''  # 引物信息表
        # self.sample_info = defaultdict(list)  # 存放样本的引物、酶、签订测序量、插入片段等信息

    def run(self):
        super(GetvaildBybarcodeTool, self).run()
        # self.get_barcode()
        self.run_getvaild_bybarcode()
        self.set_output()
        self.end()

    def run_getvaild_bybarcode(self):
        """
        运行GetvaildBybarcode,进行质控
        """
        cmd = self.getvaild_bybarcode_path + ' %s %s %s %s' % (
            self.option('fq1').prop['path'], self.option('fq2').prop['path'], self.option('barcode_info').prop['path'],
            self.work_dir + "/" + self.option('lib_name') + ".all.raw")
        self.logger.info('运行getvaild_bybarcode，进行去低值')
        command = self.add_command("getvaild_bybarcode_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("getvaild_bybarcode运行完成")
        else:
            self.set_error("getvaild_bybarcode运行出错!")

    # def get_barcode(self):
    #     """
    #     获得barcode文件及引物文件
    #     :return:
    #     """
    #     lib = set()
    #     with open(self.option('barcode_info').prop['path'])as fr:
    #             lines = fr.readlines()
    #             for line in lines[1:]:
    #                 tmp = line.strip().split("\t")
    #                 lib.add(tmp[1])
    #                 # self.sample_info[tmp[0]].append(tmp[2])  # 引物
    #                 # self.sample_info[tmp[0]].append(tmp[3])  # 酶
    #                 # self.sample_info[tmp[0]].append(tmp[4])  # 签订测序量
    #                 # self.sample_info[tmp[0]].append(tmp[5])  # 插入片段长度
    #     if len(lib) != 1:
    #         self.set_error('样本信息表中的文库名称不一致，{}，请核实！'.format(lib))
    #     self.lib_name = list(lib)[0]
        # self.barcode_path = self.output_dir + '/' + self.lib_name + '.barcode.config'
        # self.primer_path = self.output_dir + '/' + self.lib_name + '.barcode.primer.config'
        # with open(self.barcode_path, 'w+')as fw1, open(self.primer_path, 'w+')as fw2, open(self.barcode_list)as fr1, open(
        #         self.primer_list)as fr2:
        #     lines1 = fr1.readlines()
        #     lines2 = fr2.readlines()
        #     fw1.write('#Sample\tBarcode-tag\tFbarcode\tRbarcode\n')
        #     fw2.write('#Sample\tF-barcode\tLinkPrimer\tR-barcode\tReversePrimer\n')
        #     for sample in self.sample_info.keys():
        #         new_line = sample
        #         for line in lines1[1:]:
        #             tmp1 = line.strip().split("\t")
        #             if self.sample_info[sample][1] == tmp1[0]:
        #                 fw1.write(sample + '\t' + tmp1[1] + '\t' + tmp1[2] + '\t' + tmp1[3] + '\n')
        #                 for line2 in lines2[1:]:
        #                     tmp2 = line2.strip().split("\t")
        #                     two_primer = self.sample_info[sample][0].strip().split("_")
        #                     if two_primer[0] == tmp2[0]:
        #                         new_line = new_line + '\t' + tmp1[2] + '\t' + tmp2[1]
        #                     elif two_primer[1] == tmp2[0]:
        #                         new_line = new_line + '\t' + tmp1[3] + '\t' + tmp2[1]
        #         fw2.write(new_line + '\n')

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        try:
            os.link(self.work_dir + "/" + self.option('lib_name') + ".all.raw.valid.1.fq",
                    self.output_dir + '/' + self.option('lib_name') + ".all.raw.valid.1.fq")
            os.link(self.work_dir + "/" + self.option('lib_name') + ".all.raw.valid.2.fq",
                    self.output_dir + '/' + self.option('lib_name') + ".all.raw.valid.2.fq")
            self.option('out_fq1').set_path(self.output_dir + '/' + self.option('lib_name') + ".all.raw.valid.1.fq")
            self.option('out_fq2').set_path(self.output_dir + '/' + self.option('lib_name') + ".all.raw.valid.2.fq")
            # self.option('barcode_file').set_path(self.output_dir + '/' + self.lib_name + '.barcode.config')
            # self.option('primer_file').set_path(self.output_dir + '/' + self.lib_name + '.barcode.primer.config')
            self.logger.info("设置getvaild_bybarcode分析结果目录成功")

        except Exception as e:
            self.logger.info("设置getvaild_bybarcode分析结果目录失败{}".format(e))
            self.set_error("设置getvaild_bybarcode分析结果目录失败{}".format(e))
