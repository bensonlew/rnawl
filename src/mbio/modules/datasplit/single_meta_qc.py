# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

import os
import re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict


class SingleMetaQcModule(Module):
    """
    多样性对单个文库进行二次拆分及质控
    author: wangzhaoyue
    last_modify: 2017.12.05
    """

    def __init__(self, work_id):
        super(SingleMetaQcModule, self).__init__(work_id)
        options = [
            {"name": "fq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # 输入文件, 输入的单个文库文件夹，存放序列两个gz文件
            {'name': "barcode_info", "type": "infile", "format": "sequence.barcode_info"},  # 文库中样本barcode及引物信息等
            {'name': 'lib_insert_size', "type": "int"},  # 文库插入片段长度
            {'name': 'fq_type', 'type': "string", "default": "PE"},  # PE or SE
            {'name': 'leading', 'type': "string", "default": "0"},  # 切除首端碱基质量小于0的碱基或者N
            {'name': 'tailing', 'type': "string", "default": "20"},  # 切除末端碱基质量小于20的碱基或者N
            {'name': 'sliding_window', 'type': "string", "default": "50:20"},
            # 例50:20  Windows的size是50个碱基，其平均碱基质量小于20，则切除
            {'name': 'minlen', 'type': "string", "default": "50"},  # 最低reads长度
            {'name': 'valid_len', "type": "int"},  # -l,长度过滤阈值
            {'name': 'min_lenth', 'type': "string", "default": "10"},  # -m,两个reads之间所需的最小重叠长度，以提供可靠的重叠
            {'name': 'max_lenth', "type": "string", "default": "100"},  # -M,两个reads之间的最大重叠长度
            {'name': "mismatch_rate", "type": "string", "default": "0.2"},  # -x,错配和重叠长度允许的最大比率
            {'name': "pred", "type": "string", "default": "33"},  # -p,FASTQ文件中的碱基的质量值，Pred33/Pred64.
            {'name': 'thread', "type": "string", "default": "6"},  # -t,线程数
            {'name': "min_len", "type": "int"},  # -m,最小长度
            {'name': 'split_type', 'type': "string", "default": "Auto"},  # 拆分样本序列类型 Pair or Single or Auto

        ]
        self.add_option(options)
        self.fq1_file = ''  # 文库中的序列R1端
        self.fq2_file = ''  # 文库中的序列R2端
        self.lib_name = ''  # 文库名
        self.trim_length = 0  # 长度过滤阈值
        self.min_length = 0  # 最小长度
        self.tools = []  # 拼接前长度过滤的tool
        self.contract_sample_num = defaultdict(list)  # 存放 合同样本的签订测序量 合同：[样本，数据量]
        self.sample_primer = defaultdict(list)  # 存放样本的引物信息
        self.reach_num = True  # 是否达到数据量
        self.barcode_path = ''  # barcode信息表
        self.primer_path = ''  # 引物信息表
        self.R1_primer_path = ''  # 需要单端拆分的样本引物信息表
        self.R1_split_samples = []
        self.getvaild_bybarcode = self.add_tool('datasplit.getvaild_bybarcode')
        self.trimmomatic = self.add_tool('datasplit.trimmomatic')
        self.single_primer_list = ["Arch344F_Arch915R", "bamoA1F_bamoA2R", "amoAF_amoAR", "MLfF_MLrR", "A189F_mb661R"]

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('fq_dir'):
            raise OptionError('必须输入文库文件夹文件')
        if not self.option('lib_insert_size'):
            raise OptionError('必须输入文库插入片段长度')
        if not self.option('barcode_info'):
            raise OptionError('必须输入文库的样本信息表')
        if self.option('split_type') not in ["Pair", "Single", "Auto"]:
            raise OptionError('拆分类型只能是Pair或Single或自动拆分')
        return True

    def getvaild_bybarcode_run(self):
        opts = ({
            "fq1": self.fq1_file,
            "fq2": self.fq2_file,
            "barcode_info": self.barcode_path,
            "lib_name": self.lib_name,
        })
        self.getvaild_bybarcode.set_options(opts)
        self.getvaild_bybarcode.on('end', self.set_output, 'getvaild_bybarcode')
        self.getvaild_bybarcode.run()

    def trimmomatic_run(self):
        self.trimmomatic.set_options({
            "fq1": self.getvaild_bybarcode.option('out_fq1'),
            "fq2": self.getvaild_bybarcode.option('out_fq2'),
            "fq_type": self.option('fq_type'),
            "leading": self.option('leading'),
            "tailing": self.option('tailing'),
            "sliding_window": self.option('sliding_window'),
            "minlen": self.option('minlen'),
            "lib_name": self.lib_name,
        })
        self.trimmomatic.on('end', self.set_output, 'trimmomatic')
        self.trimmomatic.run()

    def trim_fqseq_front_run(self):  # 拼接前长度过滤
        n = 0
        fq_list = [self.trimmomatic.option('out_fq1'), self.trimmomatic.option('out_fq2')]
        for fq in fq_list:
            self.trim_fqseq_front = self.add_tool('datasplit.trim_fqseq')  # 拼接前的长度过滤
            self.trim_fqseq_front.set_options({
                "fq": fq,
                "valid_len": (str(self.option('valid_len')) if self.option('valid_len') else str(self.trim_length)),
                "lib_name": self.lib_name,
            })
            self.trim_fqseq_front.on('end', self.set_output, 'trim_fqseq_front{}')
            self.tools.append(self.trim_fqseq_front)
            n += 1
        self.on_rely(self.tools, self.flash_run)
        for tool in self.tools:
            tool.run()

    def flash_run(self):
        self.flash = self.add_tool('datasplit.flash')
        self.flash.set_options({
            "fq1": (self.tools[0].option('out_fq') if self.option('lib_insert_size') <= 380 else self.trimmomatic.option('out_fq1')),
            "fq2": (self.tools[1].option('out_fq') if self.option('lib_insert_size') <= 380 else self.trimmomatic.option('out_fq2')),
            "min_lenth": self.option('min_lenth'),
            "max_lenth": self.option('max_lenth'),
            "mismatch_rate": self.option('mismatch_rate'),
            "pred": self.option('pred'),
            "thread": self.option('thread'),
            "lib_name": self.lib_name,
        })
        self.flash.on('end', self.split_by_barcode_run)
        self.flash.on('end', self.set_output, 'flash')
        self.flash.run()

    def split_by_barcode_run(self):  # 拼接，拆样本
        self.split_by_barcode = self.add_tool('datasplit.split_by_barcode')
        self.split_by_barcode.set_options({
            "fq": self.flash.option('out_fq'),
            "barcode_info": self.primer_path,
            "lib_name": self.lib_name,
        })
        self.split_by_barcode.on('end', self.trim_fqseq_run)
        self.split_by_barcode.on('end', self.set_output, 'flash_split_by_barcode')
        self.split_by_barcode.run()

    def trim_fqseq_run(self):  # 拼接，长度过滤
        self.trim_fqseq = self.add_tool('datasplit.trim_fqseq')
        self.trim_fqseq.set_options({
            "fq":  self.split_by_barcode.option('out_fq'),
            "min_len": (str(self.option('min_len')) if self.option('min_len') else str(self.min_length)),
            "lib_name": self.lib_name,
        })
        self.trim_fqseq.on('end', self.fastq_extract_run)
        self.trim_fqseq.on('end', self.set_output, 'flash_trim_fqseq')
        self.trim_fqseq.run()

    def fastq_extract_run(self):  # 拼接 拆样本，统计样本的数据量
        self.fastq_extract = self.add_tool('datasplit.meta_fastq_extract')
        self.fastq_extract.set_options({
            "in_fastq": self.trim_fqseq.option('out_fq'),
        })
        self.fastq_extract.on('end', self.set_output, 'flash_fastq_extract')
        if self.option('split_type') == 'Auto':  # 自动拆分时，需根据情况判断是否需要用单端序列重新拆分
            self.fastq_extract.on('end', self.judge_data)
        else:
            self.fastq_extract.on('end', self.end)
        self.fastq_extract.run()

    def R1_split_by_barcode_run(self):  # R1,拆样本
        self.R1_split_by_barcode = self.add_tool('datasplit.split_by_barcode')
        self.R1_split_by_barcode.set_options({
            "fq": self.trimmomatic.option('out_fq1'),
            "barcode_info": (self.R1_primer_path if self.option('split_type') == 'Auto' else self.primer_path),
            "lib_name": self.lib_name,
            "split_type": "Single",
        })
        self.R1_split_by_barcode.on('end', self.R1_trim_fqseq_run)
        self.R1_split_by_barcode.on('end', self.set_output, 'R1_split_by_barcode')
        self.R1_split_by_barcode.run()

    def R1_trim_fqseq_run(self):  # 长度过滤
        self.R1_trim_fqseq = self.add_tool('datasplit.trim_fqseq')
        self.R1_trim_fqseq.set_options({
            "fq":  self.R1_split_by_barcode.option('out_fq'),
            "min_len": (str(self.option('min_len')) if self.option('min_len') else str(self.min_length)),
            "lib_name": self.lib_name,
        })
        self.R1_trim_fqseq.on('end', self.R1_fastq_extract_run)
        self.R1_trim_fqseq.on('end', self.set_output, 'R1_trim_fqseq')
        self.R1_trim_fqseq.run()

    def R1_fastq_extract_run(self):  # 拆样本，统计样本的数据量
        self.R1_fastq_extract = self.add_tool('datasplit.meta_fastq_extract')
        self.R1_fastq_extract.set_options({
            "in_fastq": self.R1_trim_fqseq.option('out_fq'),
        })
        self.R1_fastq_extract.on('end', self.set_output, 'R1_fastq_extract')
        self.R1_fastq_extract.run()

    def run(self):
        """
        运行
        :return:
        """
        super(SingleMetaQcModule, self).run()
        self.get_info()
        # time.sleep(2)
        self.getvaild_bybarcode.on('end', self.trimmomatic_run)
        if self.option('split_type') != 'Single':  # 不是单端拆分，则需要拼接
            if self.option('lib_insert_size') < 200:
                self.min_length = self.option('lib_insert_size') - 20
            else:
                self.min_length = 200
            if self.option('lib_insert_size') <= 380:
                if 0 < self.option('lib_insert_size') <= 220:
                    self.trim_length = 160
                elif 220 < self.option('lib_insert_size') <= 320:
                    self.trim_length = 200
                elif 320 < self.option('lib_insert_size') <= 380:
                    self.trim_length = 250
                self.trimmomatic.on('end', self.trim_fqseq_front_run)
            else:
                self.trimmomatic.on('end', self.flash_run)

        else:
            self.trimmomatic.on('end', self.R1_split_by_barcode_run)
        self.getvaild_bybarcode_run()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.work_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.link(oldfiles[i], newdir)

    def set_output(self, event):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if event["data"] == "getvaild_bybarcode":
            self.linkdir(self.getvaild_bybarcode.output_dir, self.output_dir + '/getvaild_bybarcode')
        if event["data"] == "trimmomatic":
            self.linkdir(self.trimmomatic.output_dir, self.output_dir + '/trimmomatic')
        if event["data"] == "trim_fqseq_front0":
            self.linkdir(self.tools[0].output_dir, self.output_dir + '/trim_fqseq_front')
        if event["data"] == "trim_fqseq_front1":
            self.linkdir(self.tools[1].output_dir, self.output_dir + '/trim_fqseq_front')
        if event["data"] == "flash":
            self.linkdir(self.flash.output_dir, self.output_dir + '/flash')
        if event["data"] == "flash_split_by_barcode":
            self.linkdir(self.split_by_barcode.output_dir, self.output_dir + '/flash_split_by_barcode')
        if event["data"] == "flash_trim_fqseq":
            self.linkdir(self.trim_fqseq.output_dir, self.output_dir + '/flash_trim_fqseq')
        if event["data"] == "flash_fastq_extract":
            if not os.path.exists(self.output_dir + '/fastq_extract'):
                os.mkdir(self.output_dir + '/fastq_extract')
            for f in os.listdir(self.fastq_extract.output_dir + '/fastq/'):
                if f.endswith("gz"):
                    f1 = self.fastq_extract.output_dir + '/fastq/' + f
                    s = f.split(".fastq.gz")[0]
                    f2 = self.output_dir + '/fastq_extract/' + self.lib_name + ":" + self.project_specimen[s] + ":" + s + ".fq.gz"
                    if os.path.exists(f2):
                        os.remove(f2)
                    os.link(f1, f2)
        if event["data"] == "R1_split_by_barcode":
            self.linkdir(self.R1_split_by_barcode.output_dir, self.output_dir + '/R1_split_by_barcode')
        if event["data"] == "R1_trim_fqseq":
            self.linkdir(self.R1_trim_fqseq.output_dir, self.output_dir + '/R1_trim_fqseq')
        if event["data"] == "R1_fastq_extract":
            all_fq = os.listdir(self.R1_fastq_extract.output_dir + '/fastq')
            if not os.path.exists(self.output_dir + '/R1_fastq_extract/'):
                os.mkdir(self.output_dir + '/R1_fastq_extract/')
            for fq in all_fq:
                if fq.endswith("gz"):
                    fq_path = self.R1_fastq_extract.output_dir + '/fastq/' + fq
                    s = fq.split(".fastq.gz")[0]
                    f2 = self.output_dir + '/R1_fastq_extract/' + self.lib_name + ":" + self.project_specimen[s] + ":" + s + ".R1.fq.gz"
                    if os.path.exists(f2):
                        os.remove(f2)
                    os.link(fq_path, f2)
            self.logger.info("设置结果目录成功")
            self.end()

    def end(self):
        super(SingleMetaQcModule, self).end()

    def get_info(self):
        """
        获取文库中的两端序列，获取barcode文件及primer文件
        """
        self.project_specimen = {}
        with open(self.option('barcode_info').prop['path'])as fr:
            lines = fr.readlines()
            tmp = lines[1].strip().split("\t")
            self.lib_name = tmp[1]
            for line in lines[1:]:
                item = line.strip().split("\t")
                self.project_specimen[item[0]] = item[2]
        self.barcode_path = self.work_dir + '/' + self.lib_name + '.barcode.config'
        self.primer_path = self.work_dir + '/' + self.lib_name + '.barcode.primer.config'
        all_fq_files = os.listdir(self.option('fq_dir').prop['path'])
        for fq in all_fq_files:
            if re.search(r'_R1_', fq):
                self.fq1_file = os.path.join(self.option('fq_dir').prop['path'], fq)
            if re.search(r'_R2_', fq):
                self.fq2_file = os.path.join(self.option('fq_dir').prop['path'], fq)
        if self.fq1_file == '' or self.fq2_file == '':
            raise Exception('文库文件夹中没有对应的R1,R2序列，请核实！')
        with open(self.option('barcode_info').prop['path'])as fr, open(self.barcode_path, 'w+')as fw1, open(
            self.primer_path, 'w+')as fw2:
            lines = fr.readlines()
            fw1.write('#Sample\tBarcode-tag\tFbarcode\tRbarcode\n')
            fw2.write('#Sample\tF-barcode\tLinkPrimer\tR-barcode\tReversePrimer\n')
            for line in lines[1:]:
                tmp = line.strip().split("\t")
                self.contract_sample_num[tmp[2]].append([tmp[0], tmp[6]])
                self.sample_primer[tmp[0]] = [tmp[3], tmp[4], tmp[7]]  # 样本的引物，引物类型，插入片段长度
                fw1.write(tmp[0] + '\t' + tmp[8] + '\t' + tmp[9] + '\t' + tmp[11] + '\n')
                fw2.write(tmp[0] + '\t' + tmp[9] + '\t' + tmp[10] + '\t' + tmp[11] + '\t' + tmp[12] + '\n')

    def judge_data(self):
        """
        拼接之后的样本序列判断数据量是否达标，达标返回True,否则返回Flase,
        当合同中的样本数据量都不达标时，我们需要判断引物，
        通用引物：单端序列，取R1端拆出来的样本
                  不是单端序列，取拼接后拆出来的样本
        外来引物：insert_size >=550,取R1端拆出来的样本
                  insert_size <550，取拼接后拆出来的样本
        :return:
        """
        real_num = {}
        error_sample = []
        info_file = self.fastq_extract.work_dir + '/info.txt'
        with open(info_file)as fr:
            lines = fr.readlines()
            for line in lines[1:]:
                tmp = line.strip().split("\t")
                real_num[tmp[1]] = tmp[3]      # 统计拼接后拆分样本的数据量
        for key in self.contract_sample_num.keys():
            reach_num = 0
            for sample in self.contract_sample_num[key]:
                try:
                    if float(real_num[sample[0]]) >= float(sample[1]):
                        reach_num += 1
                    else:
                        reach_num += 0
                except:
                    self.logger.info("没有拆出样本:{}".format(sample[0]))
            if reach_num == 0:    # 为零证明整个合同中没有一个达标的，则需进一步判断引物
                for i in self.contract_sample_num[key]:
                    error_sample.append(i)
        if len(error_sample) != 0:
            for sample in error_sample:
                if self.sample_primer[sample[0]][1] == 'mj':
                    if self.sample_primer[sample[0]][0] in self.single_primer_list:
                        self.R1_split_samples.append(sample[0])
                else:
                    if float(self.sample_primer[sample[0]][2]) >= 550:
                        self.R1_split_samples.append(sample[0])
            if len(self.R1_split_samples) > 0:
                self.R1_primer_path = self.work_dir + '/R1_' + self.lib_name + '.barcode.primer.config'
                with open(self.primer_path)as fr, open(
                        self.R1_primer_path, 'w+')as fw:
                    lines = fr.readlines()
                    fw.write('#Sample\tF-barcode\tLinkPrimer\tR-barcode\tReversePrimer\n')
                    for line in lines[1:]:
                        tmp = line.strip().split("\t")
                        if tmp[0] in self.R1_split_samples:
                            fw.write(line)
                self.R1_split_by_barcode_run()
            else:
                self.end()
        else:
            self.end()
