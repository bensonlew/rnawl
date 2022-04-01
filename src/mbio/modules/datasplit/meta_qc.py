# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

import os
import re
import time
import json
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict


class MetaQcModule(Module):
    """
    多样性对多个文库进行二次拆分及质控
    author: wangzhaoyue
    last_modify: 2017.12.11
    """

    def __init__(self, work_id):
        super(MetaQcModule, self).__init__(work_id)
        options = [
            {"name": "lib_path", "type": "infile", "format": "datasplit.path"},  # list,存放文库文件夹对应的路径信息
            {'name': "barcode_info", "type": "infile", "format": "sequence.barcode_info"},  # 文库中样本barcode及引物信息等
            {'name': "lib_specimen_id", "type": "string"},  # 文库-样本-样本的ObjectID,便于出现相同样本名称的处理
            {'name': 'lib_insert_size', "type": "string"},  # 文库插入片段长度
            {'name': 'fq_type', 'type': "string", "default": "PE"},  # PE or SE
            {'name': 'leading', 'type': "string", "default": "0"},  # 切除首端碱基质量小于0的碱基或者N
            {'name': 'tailing', 'type': "string", "default": "20"},  # 切除末端碱基质量小于20的碱基或者N
            {'name': 'sliding_window', 'type': "string", "default": "50:20"},
            # 例50:20  Windows的size是50个碱基，其平均碱基质量小于20，则切除
            {'name': 'minlen', 'type': "string", "default": "50"},  # 最低reads长度
            {'name': 'valid_len', "type": "string"},  # -l,长度过滤阈值
            {'name': 'min_lenth', 'type': "string", "default": "10"},  # -m,两个reads之间所需的最小重叠长度，以提供可靠的重叠
            {'name': 'max_lenth', "type": "string", "default": "100"},  # -M,两个reads之间的最大重叠长度
            {'name': "mismatch_rate", "type": "string", "default": "0.2"},  # -x,错配和重叠长度允许的最大比率
            {'name': "pred", "type": "string", "default": "33"},  # -p,FASTQ文件中的碱基的质量值，Pred33/Pred64.
            {'name': 'thread', "type": "string", "default": "6"},  # -t,线程数
            {'name': "min_len", "type": "string"},  # -m,最小长度
            {'name': 'split_type', 'type': "string", "default": "Auto"},  # 拆分样本序列类型 Pair or Single or Auto
        ]
        self.lib_name = {}
        self.modules = []
        self.lib_info = defaultdict(list)
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('lib_path'):
            raise OptionError('必须输入文库文件夹对应的路径信息')
        if not self.option('lib_insert_size'):
            raise OptionError('必须输入文库插入片段长度')
        if not self.option('barcode_info'):
            raise OptionError('必须输入文库的样本信息表')
        if self.option('split_type') not in ["Pair", "Single", "Auto"]:
            raise OptionError('拆分类型只能是Pair或Single或自动拆分')
        return True

    def single_meta_qc_run(self):
        n = 0
        for i in self.lib_specimen.keys():
            self.single_meta_qc = self.add_module('datasplit.single_meta_qc')
            config = self.lib_specimen[i]["barcode_info"]
            lib = self.lib_specimen[i]["lib"]
            opts = {
                "fq_dir": self.lib_name[lib],
                "barcode_info": config,
                "lib_insert_size": self.option("lib_insert_size"),
                "fq_type": self.option('fq_type'),
                "leading": self.option('leading'),
                "tailing": self.option('tailing'),
                "sliding_window": self.option('sliding_window'),
                "minlen": self.option('minlen'),
                "min_lenth": self.option('min_lenth'),
                "max_lenth": self.option('max_lenth'),
                "mismatch_rate": self.option('mismatch_rate'),
                "pred": self.option('pred'),
                "thread": self.option('thread'),
                "split_type": self.option('split_type'),
            }
            if self.option('valid_len'):
                opts.update({'valid_len': self.option('valid_len')})
            if self.option('min_len'):
                opts.update({'min_len': self.option('min_len')})
            self.single_meta_qc.set_options(opts)
            self.modules.append(self.single_meta_qc)
            n += 1
        self.on_rely(self.modules, self.set_output)
        for module in self.modules:
            module.run()

    def run(self):
        """
        运行
        :return:
        """
        super(MetaQcModule, self).run()
        self.get_info()
        f = open(self.option("lib_specimen_id"), "r")
        self.lib_specimen_id = json.loads(f.read())
        time.sleep(2)
        self.single_meta_qc_run()

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
        for file_dir in allfiles:
            new_path = os.path.join(newdir, file_dir)
            if not os.path.exists(new_path):
                os.mkdir(new_path)
            dir_path = os.path.join(dirpath, file_dir)
            allfqs = os.listdir(dir_path)
            if file_dir == "fastq_extract":
                for fq in allfqs:
                    oldfile = os.path.join(dir_path, fq)
                    lib_name = fq.split(":")[0]
                    project_sn = fq.split(":")[1]
                    specimen_name = fq.split(":")[2].split(".fq")[0]
                    id = self.lib_specimen_id[lib_name][project_sn][specimen_name]
                    newfile = os.path.join(new_path, lib_name + ":" + specimen_name + ":" + id + ".fq.gz")
                    if os.path.exists(newfile):
                        os.remove(newfile)
                    os.link(oldfile, newfile)
            elif file_dir == "R1_fastq_extract":
                for fq in allfqs:
                    oldfile = os.path.join(dir_path, fq)
                    lib_name = fq.split(":")[0]
                    project_sn = fq.split(":")[1]
                    specimen_name = fq.split(":")[2].split(".R1.fq")[0]
                    id = self.lib_specimen_id[lib_name][project_sn][specimen_name]
                    newfile = os.path.join(new_path, lib_name + ":" + specimen_name + ".R1:" + id + ".fq.gz")
                    if os.path.exists(newfile):
                        os.remove(newfile)
                    os.link(oldfile, newfile)
            else:
                for fq in allfqs:
                    oldfile = os.path.join(dir_path, fq)
                    newfile = os.path.join(new_path, fq)
                    if os.path.exists(newfile):
                        os.remove(newfile)
                    os.link(oldfile, newfile)

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for module in self.modules:
            module_name = str(module.work_dir).strip().split("/")[-1]
            self.linkdir(module.output_dir, self.output_dir + '/' + module_name)
        self.logger.info("设置结果目录成功")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(MetaQcModule, self).end()

    def get_info(self):
        """
        按照文库，将对应的文库序列及文库信息分开
        :return:
        """
        self.lib_specimen = {}
        self.all_lib = []
        i = 0
        with open(self.option('barcode_info').prop['path'])as fr:
            lines = fr.readlines()
            for line in lines[1:]:
                tmp = line.strip().split('\t')
                if tmp[1] not in self.all_lib:
                    i += 1
                    self.all_lib.append(tmp[1])
                    self.lib_specimen[str(i)] = {"lib": tmp[1], "specimen": [], "info": []}
                if tmp[0] in self.lib_specimen[str(i)]["specimen"]:
                    i += 1
                    self.lib_specimen[str(i)] = {"lib": tmp[1], "specimen": [], "info": []}
                self.lib_specimen[str(i)]["specimen"].append(tmp[0])
                self.lib_specimen[str(i)]["info"].append(line)
        for key in self.lib_specimen.keys():
            barcode_info = self.work_dir + '/' + self.lib_specimen[key]["lib"] + "." + key + '.all.barcode.primer.config'
            self.lib_specimen[key]["barcode_info"] = barcode_info
            with open(barcode_info, 'w+')as fw:
                fw.write('#Sample\tLibrary\tContract\tPrimer\tPrimer_type\tBarcode\tMinSeqNum\tInsertSize\tBarcode-tag\tF-barcode\tLinkPrimer\tR-barcode\tReversePrimer\n')
                for i in self.lib_specimen[key]["info"]:
                    fw.write(i)    # 将文库信息按照文库拆开
        with open(self.option('lib_path').prop['path'])as f:
            for line in f:
                # tmp = line.strip().split('/')
                # lib = tmp[-1].strip().split("Sample_")
                # if lib[-1] in self.all_lib:
                #     self.lib_name[lib[-1]] = line.strip()  # 换行符要去掉
                # else:
                #     raise Exception('文库路径中的文库名{}在样本信息表中不存在，请核实！'.format(lib[-1]))
                tmp = line.strip().split('\t')
                lib = tmp[0]
                if lib in self.all_lib:
                    self.lib_name[lib] = tmp[1]
                else:
                    raise Exception('文库路径中的文库名{}在样本信息表中不存在，请核实！'.format(lib))
