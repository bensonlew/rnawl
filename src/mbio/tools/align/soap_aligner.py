# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess
import shutil
from mbio.files.sequence.file_sample import FileSampleFile


class SoapAlignerAgent(Agent):
    """
    mapping by SOAP2
    author: zouxuan
    modified at date 20170911
    """

    def __init__(self, parent):
        super(SoapAlignerAgent, self).__init__(parent)
        options = [
            {"name": "fafile", "type": "infile", "format": "sequence.fasta"},  # 非冗余基因集
            {"name": "sample", "type": "string"},  # sample的名称
            {"name": "insertSize", "type": "int"},  # 插入片段长度
            {"name": "index", "type": "infile", "format": "align.bwt_index_dir"},  # build_gene生成的索引文件
            {"name": "fq_r", "type": "infile", "format": "sequence.fastq"},  # 右端fastq文件
            {"name": "fq_l", "type": "infile", "format": "sequence.fastq"},  # 左端fastq文件
            {"name": "fq_s", "type": "infile", "format": "sequence.fastq"},  # 单端fastq文件
            {"name": "map_dir", "type": "outfile", "format": "align.map_dir"},  # map结果
            {"name": "repeat", "type": "int", "default": 1},  # how to report repeat hits, 0=none, 1=random one, 2=all
            {"name": "seed", "type": "int", "default": 35},
            # align the initial n bps as a seed means whole lengths of read
            {"name": "mode", "type": "int", "default": 4},
            # match mode for each read or the seed part of read, which shouldn't contain more than 2 mismatches: 0 for exact mathc only; 1 for 1 mismatch; 2 for 2 mismatch; 4 for find the best hits
            {"name": "processors", "type": "int", "default": 6},
            {"name": "mismatch", "type": "int", "default": 20},  # maximum number of mismatches allowed on a read
            {"name": "identity", "type": "float", "default": 0.95}  # identity
        ]
        self.add_option(options)
        self.step.add_steps('soapaligner')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.soapaligner.start()
        self.step.update()

    def step_end(self):
        self.step.soapaligner.finish()
        self.step.update()

    def check_options(self):
        if not self.option("repeat") in [0, 1, 2]:
            raise OptionError("repeat必须为0,1,或2", code="31102101")
        if not self.option("mode") in [0, 1, 2, 4]:
            raise OptionError("repeat必须为0,1,2,或4", code="31102102")
        if not 0 < self.option("seed") <= 256:
            raise OptionError('seed参数必须设置在1-256之间: %s', variables=(self.option('seed')), code="31102103")
        if not 0 < self.option("identity") <= 1:
            raise OptionError("identity必须在0，1之间", code="31102104")
        #if (not self.option("fq_s").is_set) or (not (self.option("fq_l").is_set and self.option("fq_r").is_set)):
        #   raise OptionError("必须输入两条双端序列或一条单端序列或三条序列都输入", code="31102105")
        if (not (self.option("fq_l").is_set and self.option("fq_r").is_set)):
           raise OptionError("必须输入两条双端序列或一条单端序列或三条序列都输入", code="31102105")
        if not self.option("sample"):
            raise OptionError("必须输入样品名", code="31102106")
        if not self.option("insertSize"):
            raise OptionError("必须输入插入片段长度", code="31102107")

    def set_resource(self):
        self._cpu = 15
        if os.path.getsize(self.option("fafile").prop['path'])/100000000< 15:
            self._memory = '20G'
        else:
            self._memory = str(os.path.getsize(self.option("fafile").prop['path'])/150000000 + 10) + 'G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(SoapAlignerAgent, self).end()


class SoapAlignerTool(Tool):
    def __init__(self, config):
        super(SoapAlignerTool, self).__init__(config)
        self._version = "1.1"
        self.soap_path = 'bioinfo/uniGene/soap2.21release/soap'
        self.index = self.option("index").prop['path'] + '/' +  os.path.basename(self.option("fafile").prop['path'])+'.index'

    def run(self):
        super(SoapAlignerTool, self).run()
        self.make_file()
        if self.option("fq_s").is_set:
            self.single_map()
        if self.option("fq_l").is_set and self.option("fq_r").is_set:
            self.pair_map()
        self.set_output()

    def make_file(self):
        """
        生成以样品名为文件名的文件夹
        """
        if os.path.exists(os.path.join(self.work_dir, self.option("sample"))):
            pass
        else:
            os.mkdir(os.path.join(self.work_dir, self.option("sample")))

    def pair_map(self):
        thesample = self.option("sample")
        out_dir = self.work_dir + '/' + thesample
        pe_cmd = '%s -a %s -b %s -D %s -o %s -2 %s' % (
            self.soap_path, self.option("fq_r").prop['path'], self.option("fq_l").prop['path'], self.index,
            out_dir + '/' + thesample + '.soap.pair.pe',
            out_dir + '/' + thesample + '.soap.pair.se')  ##setting input and output file
        pe_cmd += ' -r %s -l %s -M %s -p %s -v %s -c %s -m %s -x %s 2' % (
            self.option("repeat"), self.option("seed"), self.option("mode"), self.option("processors"),
            self.option("mismatch"), self.option("identity"), self.option("insertSize") - 100,
            self.option("insertSize") + 100)  ##setting parameters
        self.logger.info(pe_cmd)
        self.logger.info(self.soap_path)
        command1 = self.add_command('pair_map', pe_cmd).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("pair_map succeed")
        else:
            self.set_error("pair_map failed", code="31102101")
            raise Exception("pair_map failed")

    def single_map(self):
        thesample = self.option("sample")
        out_dir = self.work_dir + '/' + thesample
        se_cmd = '%s -a %s -D %s -o %s' % (self.soap_path, self.option("fq_s").prop['path'], self.index,
                                           out_dir + '/' + thesample + '.soap.single.se')  ##setting input and output file
        se_cmd += ' -r %s -l %s -M %s -p %s -v %s -c %s -m %s -x %s 2' % (
            self.option("repeat"), self.option("seed"), self.option("mode"), self.option("processors"),
            self.option("mismatch"), self.option("identity"), self.option("insertSize") - 100,
            self.option("insertSize") + 100)  ##setting parameters
        self.logger.info(se_cmd)
        command2 = self.add_command('single_map', se_cmd).run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("single_map succeed")
        else:
            self.set_error("single_map failed", code="31103102")
            raise Exception("single_map failed")

    def set_output(self):
        self.logger.info("set output")
        self.linkdir(os.path.join(self.work_dir, self.option("sample")), self.option("sample"))
        self.option('map_dir', self.output_dir)
        self.end()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])
