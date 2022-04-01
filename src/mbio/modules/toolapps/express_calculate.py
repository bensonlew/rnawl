# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import os, re
import shutil
import gevent
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
from mbio.packages.toolapps.common import link_dir


class ExpressCalculateModule(Module):
    """
    小工具：宏转录组 表达量计算
    运用软件：RSEM软件
    kalisto 和 salmon软件参数与转录组不同，转录组不可直接用
    """
    def __init__(self, work_id):
        super(ExpressCalculateModule, self).__init__(work_id)
        options = [
            {"name": "fq_type", "type": "string"},  # PE OR SE
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},  # trinity.fasta文件
            {"name": "ref_dir", "type": "infile", "format": "toolapps.fasta_dir"},  # trinity.fasta文件夹
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"}, # PE数据或者SE数据，包含所有的左端和右端的fq文件或者single端文件
            {"name": "exp_way", "type": "string", "default": "fpkm"}, #计算表达量的指标或者方法
            {"name": "software", "type": "string", 'default': "rsem"}, ## 计算使用的软件
        ]
        self.add_option(options)
        self.sequence = self.add_module('metagenome.reads_unzip')
        self.bowtie_build = self.add_tool("toolapps.express_rsem")
        self.merge_rsem = self.add_tool("toolapps.merge_rsem")
        self.statistics = self.add_tool("toolapps.express_statistics")
        self.run_tools = []
        self.bam_path = self.work_dir + '/bowtie2_bam_dir/'
        self.rsem_path = os.path.join(self.work_dir, "result_out")

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        if not self.option('fq_type'):
            raise OptionError('必须设置测序类型：PE OR SE！', code="24400501")
        if self.option('fq_type') not in ['PE', 'SE']:
            raise OptionError('测序类型不在所给范围内！', code="24400502")
        if self.option("exp_way") not in ['fpkm', 'tpm']:
            raise OptionError("所设表达量的代表指标不在范围内，请检查！", code="24400503")
        if not self.option('ref_fa').is_set and not self.option("ref_dir").is_set:
            raise OptionError('必须设置参考基因组！', code="24400504")
        if not self.option("fastq_dir").is_set:
            raise OptionError('必须要设置测序序列文件夹！', code="24400505")
        if self.option("software") not in ["rsem", "kallisto", "salmon"]:
            raise OptionError('所设置的软件名称不在["rsem", "kallisto", "salmon"]中', code="24400506")
        if self.option("software") in ['rsem']:
            if self.option("exp_way") not in ["fpkm", "tpm"]:
                raise OptionError('所设置的指标名称不在["fpkm", "tpm"]中', code="24400507")
        elif self.option("software") in ['kallisto', 'salmon']:
            if self.option("exp_way") not in ["tpm"]:
                raise OptionError('所设置的指标名称不在["tpm"]中', code="24400508")
        else:
            raise OptionError('错误的软件名称，请检查软件名称', code="24400509")
        return True

    def run_sequence(self):
        """
        对数据进行解压，上传序列为压缩格式
        :return:
        """
        self.logger.info("开始进行解压文件")
        opts = {
            'fastq_dir': self.option('fastq_dir'),
        }
        self.sequence.set_options(opts)
        self.sequence.run()

    def run_bowtie_build(self):
        """
        bowtie2建索引
        :return:
        """
        self.logger.info("开始用bowtie2建立索引！")
        opts = {
            'fq_type': self.option('fq_type'),
            'rsem_fa': self.option('ref_fa'),
            'only_bowtie_build': True
            }
        self.bowtie_build.set_options(opts)
        self.bowtie_build.run()

    def run_rsem(self):
        """
        运行RSEM软件计算表达量
        :return:
        """
        self.logger.info("开始用RSEM软件计算表达量")
        samples = self.get_list()
        if not os.path.exists(self.bam_path):
            os.mkdir(self.bam_path)
        else:
            shutil.rmtree(self.bam_path)
            os.mkdir(self.bam_path)

        for sample in samples:
            type_dict = samples[sample]
            self.rsem = self.add_tool('toolapps.express_rsem')
            if not self.option("ref_dir").is_set:
                opts = {
                    'fq_type': self.option('fq_type'),
                    'only_bowtie_build': False,
                    'rsem_fa': self.bowtie_build.option('fa_build'),
                    'sample': sample,
                    }
            else:
                fasta_samples = self.get_fasta_list()
                self.logger.info("samples:{}".format(fasta_samples))
                opts = {
                    'fq_type': self.option('fq_type'),
                    'only_bowtie_build': False,
                    'bowtie_build_rsem': True,
                    'rsem_fa': fasta_samples[sample],
                    'sample': sample,
                    }
            for direct in type_dict.keys():
                if direct in ['l']:
                    opts["fq_l"] = type_dict[direct]
                elif direct in ['r']:
                    opts["fq_r"] = type_dict[direct]
                elif direct in ['s']:
                    opts["fq_s"] = type_dict[direct]
                else:
                    raise OptionError('list.txt文件含有第三列中有错误的类型%s', variables=(direct), code="24400510")
            self.rsem.set_options(opts)
            self.run_tools.append(self.rsem)
        if len(self.run_tools) > 1:
            self.on_rely(self.run_tools, self.run_merge_rsem)
        else:
            self.run_tools[0].on("end", self.run_merge_rsem)
        for tool in self.run_tools:
            tool.run()
            gevent.sleep(0)

    def run_merge_rsem(self):
        """
        将RSEM软件计算的结果进行合并，并按照不同的指标方法进行计算
        :return:
        """
        self.logger.info("合并RSEM软件的结果")
        if os.path.exists(self.rsem_path):
            shutil.rmtree(self.rsem_path)
            os.mkdir(self.rsem_path)
        else:
            os.mkdir(self.rsem_path)
        samples = self.get_list()
        if len(samples.keys()) > 1:
            for tool in self.run_tools:
                link_dir(tool.output_dir, self.rsem_path)
            self.merge_rsem.set_options({
                'rsem_files': self.rsem_path,
                'exp_way': self.option('exp_way')
            })
            self.merge_rsem.on('end', self.set_output)
            self.merge_rsem.run()
        else:
            self.set_output()

    def run_kallisto(self):
        """
        运行软件kallisto软件计算表达量
        :return:
        """
        self.logger.info("开始用Kallisto软件计算表达量")
        samples = self.get_list()
        for sample in samples:
            type_dict = samples[sample]
            self.kallisto = self.add_tool('toolapps.kallisto')
            if not self.option("ref_dir").is_set:
                opts = {
                    'fq_type': self.option('fq_type'),
                    'ref_fa': self.option("ref_fa"),
                    'sample': sample,
                    }
            else:
                fasta_samples = self.get_fasta_list()
                opts = {
                    'fq_type': self.option('fq_type'),
                    'ref_fa': fasta_samples[sample],
                    'sample': sample,
                    }
            for direct in type_dict.keys():
                if direct in ['l']:
                    opts["fq_l"] = type_dict[direct]
                elif direct in ['r']:
                    opts["fq_r"] = type_dict[direct]
                elif direct in ['s']:
                    opts["fq_s"] = type_dict[direct]
                else:
                    raise OptionError('list.txt文件含有第三列中有错误的类型%s', variables=(direct), code="24400511")
            self.kallisto.set_options(opts)
            self.run_tools.append(self.kallisto)
        if len(self.run_tools) > 1:
            self.on_rely(self.run_tools, self.run_statistic)
        else:
            self.run_tools[0].on("end", self.run_statistic)
        for tool in self.run_tools:
            tool.run()
            gevent.sleep(0)

    def run_salmon(self):
        """
        运行salmon软件计算表达量
        :return:
        """
        self.logger.info("开始用Kallisto软件计算表达量")
        samples = self.get_list()
        for sample in samples:
            type_dict = samples[sample]
            self.salmon = self.add_tool('toolapps.salmon')
            if not self.option("ref_dir").is_set:
                opts = {
                    'fq_type': self.option('fq_type'),
                    'ref_fa': self.option("ref_fa"),
                    'sample': sample,
                    }
            else:
                fasta_samples = self.get_fasta_list()
                opts = {
                    'fq_type': self.option('fq_type'),
                    'ref_fa': fasta_samples[sample],
                    'sample': sample,
                    }
            for direct in type_dict.keys():
                if direct in ['l']:
                    opts["fq_l"] = type_dict[direct]
                elif direct in ['r']:
                    opts["fq_r"] = type_dict[direct]
                elif direct in ['s']:
                    opts["fq_s"] = type_dict[direct]
                else:
                    raise OptionError('list.txt文件含有第三列中有错误的类型%s', variables=(direct), code="24400512")
            self.salmon.set_options(opts)
            self.run_tools.append(self.salmon)
        if len(self.run_tools) > 1:
            self.on_rely(self.run_tools, self.run_statistic)
        else:
            self.run_tools[0].on("end", self.run_statistic)
        for tool in self.run_tools:
            tool.run()
            gevent.sleep(0)

    def run_statistic(self):
        """
        对salmon软件和kallisto软件的结果进行合并
        :return:
        """
        self.logger.info("合并软件计算的结果")
        if os.path.exists(self.rsem_path):
            shutil.rmtree(self.rsem_path)
            os.mkdir(self.rsem_path)
        else:
            os.mkdir(self.rsem_path)

        for tool in self.run_tools:
            link_dir(tool.output_dir, self.rsem_path)
        self.statistics.set_options({
            'merge_files': self.rsem_path,
            'exp_way': self.option('exp_way'),
            'software': self.option("software")
        })
        self.statistics.on('end', self.set_output)
        self.statistics.run()

    def get_list(self):
        """
        根据上传的fq文件夹list获取对应关系
        list：文件名称 样本名称 类型（l、r、s）
        :return:
        """
        self.logger.info("根据fq文件夹的list文件获取对应关系")
        list_path = os.path.join(self.sequence.output_dir,"data","list.txt")
        if os.path.exists(list_path):
            self.logger.info(list_path)
        sample = {}
        with open(list_path, "rb") as l:
            for line in l:
                line = line.strip().split("\t")
                if len(line) == 3:
                    sample_path = os.path.join(self.sequence.output_dir,"data", line[0])
                    if line[1] not in sample:
                        sample[line[1]] = {line[2]: sample_path}
                    elif line[2] in sample[line[1]].keys(): ##相同样本名称，相同类型的文件以空格相隔
                        sample[line[1]][line[2]] = sample[line[1]][line[2]] + " " + sample_path
                    else:
                        sample[line[1]][line[2]] = sample_path
                else:
                    raise OptionError('list.txt文件格式有误', code="24400513")
        return sample

    def get_fasta_list(self):
        """
        根据上传的fasta文件夹获取样本名称与文件名称的对应关系
        :return:
        list 格式 文件名称 样本名称
        """
        self.logger.info("根据fq文件夹的list文件获取对应关系")
        sample = {}
        list_path = os.path.join(self.option("ref_dir").prop["fasta_dir"],"list.txt")
        if os.path.exists(list_path):
            self.logger.info(list_path)
            with open(list_path, "rb") as l:
                for line in l:
                    line = line.strip().split("\t")
                    sample_path = line[0]##line[0]为path
                    if len(line) == 2:
                        if line[1] not in sample:
                            sample[line[1]] = sample_path
                        else:
                            raise OptionError('样本名称存在重复: %s', variables=(line[1]), code="24400514")
                    else:
                        raise OptionError('list.txt文件格式有误', code="24400515")
        else:
            all_files = os.listdir(self.option("ref_dir").prop["fasta_dir"])
            for file in all_files:
                file_path = os.path.join(self.option("ref_dir").prop["fasta_dir"], file)
                if re.search(r'\.fa', file):
                    sample_name = file.strip('.fa')
                elif re.search(r'\.fasta', file):
                    sample_name = file.strip('.fasta')
                if sample_name not in sample.keys():
                    sample[sample_name] = file_path
                else:
                    raise OptionError('样本名称存在重复: %s', variables=(sample_name), code="24400516")
        return sample

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始设置结果文件目录！")
        if self.option("software") in ["rsem"]:
            for tool in self.run_tools:
                link_dir(tool.output_dir, os.path.join(self.output_dir, 'express_out'))
            if os.path.exists(os.path.join(self.output_dir, 'merge_statistics')):
                shutil.rmtree(os.path.join(self.output_dir, 'merge_statistics'))
            samples = self.get_list()
            if len(samples.keys()) > 1:
                link_dir(self.merge_rsem.output_dir, os.path.join(self.output_dir, 'merge_statistics'))
        else:
            for tool in self.run_tools:
                link_dir(tool.output_dir, os.path.join(self.output_dir, 'express_out'))
            if os.path.exists(os.path.join(self.output_dir, 'merge_statistics')):
                shutil.rmtree(os.path.join(self.output_dir, 'merge_statistics'))
            link_dir(self.statistics.output_dir, os.path.join(self.output_dir, 'merge_statistics'))
        self.logger.info("设置结果文件目录完成！")
        self.end()

    def run(self):
        """
        运行
        :return:
        """
        self.logger.info("start running！")
        super(ExpressCalculateModule, self).run()
        if self.option("software") in ["rsem"]:
            if self.option("ref_dir").is_set:
                self.sequence.on("end", self.run_rsem)
                self.run_sequence()
            else:
                self.sequence.on("end", self.run_bowtie_build)
                self.bowtie_build.on('end', self.run_rsem)
                self.run_sequence()
        elif self.option("software") in ["kallisto"]:
            self.sequence.on("end", self.run_kallisto)
            self.run_sequence()
        elif self.option("software") in ["salmon"]:
            self.sequence.on("end", self.run_salmon)
            self.run_sequence()
        else:
            raise OptionError('错误的软件名称，请检查软件名称', code="24400517")

    def end(self):
        """
        结束
        :return:
        """
        self.logger.info("开始上传结果文件目录")
        repaths = [
            [".", "", "表达量分析模块结果输出目录"]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        super(ExpressCalculateModule, self).end()
