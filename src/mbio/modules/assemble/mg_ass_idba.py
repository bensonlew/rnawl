# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import os
import re
import math
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from biocluster.config import Config
from bson.objectid import ObjectId
import gevent


class MgAssIdbaModule(Module):
    """
    宏基因运用idba/megait组装,单拼&混拼
    author: guhaidong
    last_modify: 2018.10.23
    """

    def __init__(self, work_id):
        super(MgAssIdbaModule, self).__init__(work_id)
        options = [
            # {"name": "data_id", "type": "string"},  # 主表任务ID，导出测序量与样品类型
            {"name": "raw_stat", "type": "infile", "format": "sequence.profile_table"},  # 输入文件，原始序列统计文件
            {"name": "QC_dir", "type": "infile", "format": "sequence.fastq_dir"},  # 输入文件，质控后的文件夹
            {"name": "qc_stat", "type": "infile", "format": "sequence.profile_table"},  # 输入文件，质控或去宿主后的统计结果
            {"name": "assemble_tool", "type": "string", "default": "idba"},  # 拼接工具选择{idba|megahit|automatic}
            {"name": "method", "type": "string"},  # 拼接方法选择{simple|multiple}
            {"name": "use_newbler", "type": "string", "default": "True"},
            {"name": "min_contig", "type": "int", "default": 300},  # 输入最短contig长度，默认300
            {"name": "contig", "type": "outfile", "format": "sequence.fasta_dir"},  # 输出contig路径, sample.contig.fa
            {"name": "contig_stat", "type": "outfile", "format": "sequence.profile_table"},  # 输出contig质量统计结果表
        ]
        self.add_option(options)
        self.sample = []
        self.qc_file = {}  # 质控数据信息
        self.sum_tools = []
        # self.tools = []
        self.single_module = []  # idba拼接，混拼后加入第一次混拼过程
        self.bowtie_module = []  # bowtie2 map reads
        self.extract_fq_module = []  # 根据sam文件提取reads
        self.zip_length = []
        self.cat_file = []
        self.mix_module = []
        self.number = 1
        '''
        if self.option('assemble_tool') == 'idba':
            self.idba = self.add_tool('assemble.idba')
            self.step.add_steps("idba")
        elif self.option('assemble_tool') == 'megahit':
            self.megahit = self.add_tool('assemble.megahit')
            self.step.add_steps('megahit')
        if self.option('method') == 'simple':
            self.contig_stat = self.add_tool("assemble.contig_stat")
            self.len_distribute = self.add_tool("sequence.length_distribute")
            self.step.add_steps("contig_stat", "length_distribute")
        elif self.option('method') == 'multiple':
            self.bowtie2 = self.add_tool("align.bowtie2")
            self.extract_fq = self.add_tool('sequence.extract_fastq_by_sam')
            self.cat_reads = self.add_tool("sequence.cat_reads")
            self.mix_assem = self.add_tool("assemble.megahit")
            self.cut_length = self.add_tool('sequence.cut_length')
            self.newbler = self.add_tool("assemble.newbler")
            self.sort_result = self.add_tool("assemble.sort_idba_result")
            self.contig_stat = self.add_tool("assemble.contig_stat")
            self.len_distribute = self.add_tool("sequence.length_distribute")
            self.step.add_steps("bowtie2", "extract_fq", "cat_reads", "mix_assem", "cut_length", "newbler",
                                "sort_result", "contig_stat", "length_distribute")
        '''
        self.idba = self.add_tool('assemble.idba')
        #self.megahit = self.add_tool('assemble.megahit')
        self.bowtie2 = self.add_tool("align.bowtie2")
        self.extract_fq = self.add_tool('sequence.extract_fastq_by_sam')
        self.cut_length = self.add_tool('sequence.cut_length')
        self.newbler = self.add_tool("assemble.newbler")
        #self.mix_assem = self.add_tool("assemble.megahit")
        self.sort_result = self.add_tool("assemble.sort_idba_result")
        self.contig_stat = self.add_tool("assemble.contig_stat")
        self.len_distribute = self.add_tool("sequence.length_distribute")
        self.yasuo = self.add_tool("sequence.zip")
        self.cat_max_reads = self.add_tool("metagenomic.merge_file")
        self.step.add_steps("idba", "megahit", "bowtie2", "extract_fq", "cat_reads", "mix_assem", "cut_length",
                            "newbler", "sort_result", "contig_stat", "length_distribute")

    def check_options(self):
        """
        检查参数
        :return:
        """
        '''
        if not self.option('data_id'):
            raise OptionError("必须输入任务ID")
        elif len(self.option('data_id')) != 24:
            raise OptionError("任务ID长度错误")
        '''
        if not self.option('raw_stat'):
            raise OptionError('必须输入原始序列统计文件', code="21300101")
        if not self.option('QC_dir'):
            raise OptionError('必须输入质控后的fq文件夹', code="21300102")
        if not self.option('qc_stat'):
            raise OptionError('必须输入质控后的统计文件', code="21300103")
        if not self.option('method'):
            raise OptionError('必须输入拼接方法', code="21300104")
        if self.option('method') not in ['simple', "multiple"]:
            raise OptionError('拼接方法必须为simple/multiple之一', code="21300105")
        if self.option("use_newbler") not in ['True', 'False']:
            raise OptionError("option use_newbler must in range ['True', 'False']")
        elif self.option("use_newbler") == "True" and self.option("method") == "simple":
            raise OptionError("Simple assemble can't support newbler")
        return True

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def idba_run(self):
        """
        idba组装
        :return:
        """
        n = 0
        # db = Config().mongo_client.tsanger_metagenomic
        # collection = db['mg_data_stat']
        # object_id = ObjectId(self.option('data_id'))
        self.qc_file = self.get_list()
        '''
        results = collection.find({'data_stat_id': object_id})
        if not results.count():
            raise Exception('没有找到样品集数据')
        if results is None:
            raise Exception('没有找到样品集数据2')
        raw_rd_len, base_num, insert_dic, sample_type = self.get_dic(results)
        '''
        raw_rd_len, base_num, insert_dic, sample_type = self.get_dic()
        for key in insert_dic.keys():
            assem_mem, split_num = self.get_mem(sample_type[key], base_num[key])  # 计算运行内存及是否需要拆分
            # self.logger.info('type is ' + sample_type[key] + '; base_num is ' + base_num[key] + " assem_mem is " +
            # str(assem_mem) + '\n')
            self.sample.append(key)
            self.idba = self.add_tool('assemble.idba')
            self.step.add_steps('IDBA_{}'.format(n))
            opts = ({
                "fastq1": self.option('QC_dir').prop['path'] + '/' + self.qc_file[key]['l'],
                "fastq2": self.option('QC_dir').prop['path'] + '/' + self.qc_file[key]['r'],
                "sample_name": key,
                "split_num": split_num,
                "mem": 100 if assem_mem < 100 else assem_mem,  # assem_mem，测试提供
            })
            if self.option('method') == 'simple':
                opts['min_contig'] = self.option('min_contig')
            else:
                opts['min_contig'] = self.option('min_contig')
            #if 's' in self.qc_file[key].keys():
                #opts['fastqs'] = self.option('QC_dir').prop['path'] + '/' + self.qc_file[key]['s']
            self.idba.set_options(opts)
            step = getattr(self.step, 'IDBA_{}'.format(n))
            step.start()
            self.idba.on('end', self.finish_update, 'IDBA_{}'.format(n))
            self.single_module.append(self.idba)
            self.sum_tools.append(self.idba)
            n += 1
        if len(self.single_module) == 1:
            if self.option('method') == 'multiple':
                self.single_module[0].on('end', self.bowtie2_run)
            else:
                self.single_module[0].on('end', self.contig_stat_run)
        else:
            if self.option('method') == 'multiple':
                self.on_rely(self.single_module, self.bowtie2_run)
                self.step.bowtie2.start()
            else:
                self.on_rely(self.single_module, self.contig_stat_run)
                self.step.contig_stat.start()
            self.step.update()
        for module in self.single_module:
            module.run()
            gevent.sleep(0)

    def megahit_run(self):
        """
        进行megahit拼接
        :return:
        """
        self.number = 1
        self.qc_file = self.get_list()
        '''
        db = Config().mongo_client.tsanger_metagenomic
        collection = db['mg_data_stat']
        object_id = ObjectId(self.option('data_id'))
        self.qc_file = self.get_list()
        results = collection.find({'data_stat_id': object_id})
        if not results.count():
            raise Exception('没有找到样品集数据')
        if results is None:
            raise Exception('没有找到样品集数据2')
        raw_rd_len, base_num, insert_dic, sample_type = self.get_dic(results)
        '''
        raw_rd_len, base_num, insert_dic, sample_type = self.get_dic()
        for key in insert_dic.keys():
            assem_mem, split_num = self.get_mem(sample_type[key], base_num[key], tool="megahit")  # 计算运行内存及是否需要拆分
            self.sample.append(key)
            self.megahit = self.add_tool('assemble.megahit')
            self.step.add_steps('MEGAHIT_{}'.format(self.number))
            opts = ({
                "fastq1": self.option('QC_dir').prop['path'] + '/' + self.qc_file[key]['l'],
                "fastq2": self.option('QC_dir').prop['path'] + '/' + self.qc_file[key]['r'],
                "sample_name": key,
                "mem": 50 if assem_mem < 50 else assem_mem,  # assem_mem，测试提供
                "mem_mode": 'moderate'  # 使用普通内存参数运行megahit @ 20180328
            })
            if self.option('method') == 'simple':
                opts['min_contig'] = self.option('min_contig')
            else:
                opts['min_contig'] = self.option('min_contig')
            #if 's' in self.qc_file[key].keys():
                #opts['fastqs'] = self.option('QC_dir').prop['path'] + '/' + self.qc_file[key]['s']
            self.megahit.set_options(opts)
            step = getattr(self.step, 'MEGAHIT_{}'.format(self.number))
            step.start()
            self.megahit.on('end', self.finish_update, 'MEGAHIT_{}'.format(self.number))
            self.single_module.append(self.megahit)
            self.sum_tools.append(self.megahit)
            self.number += 1
        if len(self.single_module) == 1:
            if self.option('method') == 'multiple':
                self.single_module[0].on('end', self.bowtie2_run)
            else:
                self.single_module[0].on('end', self.contig_stat_run)
        else:
            if self.option('method') == 'multiple':
                self.on_rely(self.single_module, self.bowtie2_run)
                self.step.bowtie2.start()
            else:
                self.on_rely(self.single_module, self.contig_stat_run)
                self.step.contig_stat.start()
            self.step.update()
        for module in self.single_module:
            module.run()
            gevent.sleep(0)

    def bowtie2_run(self):
        """
        进行bowtie2比对
        :return:
        """
        n = 0
        self.get_contig_file()
        for samples in self.sample:
            self.bowtie2 = self.add_tool("align.bowtie2")
            self.step.add_steps('bowtie2_{}'.format(n))
            opts = ({
                'ref_fasta': self.work_dir + '/contig_dir/' + samples + '.contig.fa',
                'fastq1': self.option('QC_dir').prop['path'] + '/' + self.qc_file[samples]['l'],
                'fastq2': self.option('QC_dir').prop['path'] + '/' + self.qc_file[samples]['r'],
            })
            #if 's' in self.qc_file[samples].keys():
                #opts['fastqs'] = self.option('QC_dir').prop['path'] + '/' + self.qc_file[samples]['s']
            self.bowtie2.set_options(opts)
            step = getattr(self.step, 'bowtie2_{}'.format(n))
            step.start()
            self.bowtie2.on('end', self.finish_update, 'bowtie2_{}'.format(n))
            self.bowtie_module.append(self.bowtie2)
            self.sum_tools.append(self.bowtie2)
            n += 1
        if len(self.bowtie_module) == 1:
            self.bowtie_module[0].on('end', self.extract_fq_run)
        else:
            self.on_rely(self.bowtie_module, self.extract_fq_run)
            self.step.extract_fq.start()
            self.step.update()
        for module in self.bowtie_module:
            module.run()
            gevent.sleep(0)

    def extract_fq_run(self):
        """
        根据bowtie2比对结果，挑选出fastq序列
        :return:
        """
        n = 0
        for module in self.bowtie_module:
            self.extract_fq = self.add_tool('sequence.extract_fastq_by_sam')
            self.step.add_steps('extract_fq_{}'.format(n))
            opts = ({
                'sam': module.option('sam_file'),  # 测试一下这样传参
            })
            file_dir = os.listdir(module.option('sam_file').prop['path'])
            # name = os.path.basename(module.option('sam_file').prop['path']).split('.')[0]
            for files in file_dir:
                if 'pair.sam' in files:
                    opts['fq_type'] = 'PE'
                    break
            self.extract_fq.set_options(opts)
            step = getattr(self.step, 'extract_fq_{}'.format(n))
            step.start()
            self.extract_fq.on('end', self.finish_update, 'extract_fq_{}'.format(n))
            self.extract_fq_module.append(self.extract_fq)
            self.sum_tools.append(self.extract_fq)
            n += 1
        if len(self.extract_fq_module) == 1:
            self.extract_fq_module[0].on('end', self.cat_reads_run)
        else:
            self.on_rely(self.extract_fq_module, self.cat_reads_run)
            self.step.cat_reads.start()
            self.step.update()
        for module in self.extract_fq_module:
            module.run()
            gevent.sleep(0)

    def cat_reads_run(self):
        """
        将bowtie2比对得到的fastq按照pes的关系归为三类，按照统一的样品次序分别cat在一起，作为混合拼接的输入文件
        :return:
        """
        self.get_fq_file()
        n = self.get_fq_sort_file()
        self.logger.info("切分文件个数：%s" % n)
        for i in range(1, n+1):
            self.cat_reads = self.add_tool("sequence.cat_reads")
            self.cat_reads.set_options({
                'map_dir': self.work_dir + '/sort_size_%s' % i,
            })
            self.cat_file.append(self.cat_reads)
        if len(self.cat_file) == 1:
            self.cat_file[0].on('end', self.mix_assem_run)
        else:
            self.on_rely(self.cat_file, self.mix_assem_run)
        for tool in self.cat_file:
            tool.run()
            gevent.sleep(0)
        #self.step.cat_reads.finish()
        #self.step.mix_assem.start()
        #self.step.update()

    def mix_assem_run(self):
        """
        对没有mapping到contig的reads进行混合拼接
        """
        self.logger.info("megahit第几个tool：{}".format(self.number))
        #self.number += 2
        for module in self.cat_file:
            self.mix_assem = self.add_tool("assemble.megahit")
            self.step.add_steps('MEGAHIT_{}'.format(self.number))
            opts = ({
                'mem': 200,  # assem_mem 混拼使用大内存
                'mem_mode': 'moderate',  # 内存使用模式不再指定内存大小，便于Agent修改内存 @ 20180328
                'sample_name': 'Megahit_Mix',
            })
            file_list = os.listdir(module.output_dir)
            for file in file_list:
                if file.endswith('l.fq'):
                    opts['fastq1'] = os.path.join(module.output_dir, file)
                elif file.endswith('r.fq'):
                    opts['fastq2'] = os.path.join(module.output_dir, file)
            #elif file.endswith('s.fq'):
                #opts['fastqs'] = os.path.join(self.cat_reads.output_dir, file)
            self.mix_assem.set_options(opts)
            self.mix_module.append(self.mix_assem)
            step = getattr(self.step, 'MEGAHIT_{}'.format(self.number))
            step.start()
            self.mix_assem.on('end', self.finish_update, 'MEGAHIT_{}'.format(self.number))
            self.number += 1
        if len(self.mix_module) == 1:
            if self.option("use_newbler") == "False":  # 增加newbler参数
                self.mix_module[0].on("end", self.cat_file_reads)
            else:
                self.mix_module[0].on('end', self.cut_length_run)
        else:
            if self.option("use_newbler") == "False":  # 增加newbler参数
                self.on_rely(self.mix_module, self.cat_file_reads)
            else:
                self.on_rely(self.mix_module, self.cut_length_run)
        # self.mix_assem.on('end', self.cut_length_run)
        #self.single_module.append(self.mix_assem)
        for tool in self.mix_module:
            tool.run()
            gevent.sleep(0)
        #self.mix_assem.run()
        #self.step.mix_assem.finish()
        #self.step.cut_length.start()  # step无效，所以不进行修改
        #self.step.update()

    def cat_file_reads(self):
        """
        将多个混拼结果cat到一起，方便后续进行统计
        :return:
        """
        number = 0
        file_dir = os.path.join(self.work_dir, 'tmp_mix')
        if os.path.exists(file_dir):
            shutil.rmtree(file_dir)
            os.mkdir(file_dir)
        else:
            os.mkdir(file_dir)
        for tool in self.mix_module:
            number += 1
            for file in os.listdir(tool.output_dir):
                if file.endswith('contig.fa'):
                    old_file = os.path.join(tool.output_dir, file)
                    new_file = file_dir + '/mix_file_{}.fa'.format(str(number))
                    os.link(old_file, new_file)
        self.cat_max_reads.set_options({
            "fa_dir": file_dir
        })
        self.cat_max_reads.on("end", self.contig_stat_run)
        self.cat_max_reads.run()
        self.single_module.append(self.cat_max_reads)

    def cut_length_run(self):
        """
        按照统一长度对contig进行拆分
        :return:
        """
        self.get_contig_file()
        self.cut_length.set_options({
            'contig': self.work_dir + '/contig_dir'
        })
        self.cut_length.on('end', self.newbler_run)
        self.cut_length.run()
        self.step.cut_length.finish()
        self.step.cut_length.start()
        self.step.update()

    def newbler_run(self):
        """
        进行newbler拼接
        :return:
        """
        self.newbler.set_options({
            'contig': self.cut_length.option('short_contig'),
            'all_length': self.option('min_contig'),
            'mem': 150,
            'cpu': 1
        })  # 增加cpu参数为1 by GHD @ 20180502
        # 增加内存至150G by GHD @ 20180720
        self.newbler.on('end', self.sort_idba_result_run)
        self.newbler.run()
        self.step.newbler.finish()
        self.step.update()

    def sort_idba_result_run(self):
        """
        整合混拼结果
        :return:
        """
        self.sort_result.set_options({
            'idba_contig': self.cut_length.option('cut_contig'),
            'newbler': self.newbler.output_dir,
            'min_contig': self.option('min_contig'),
        })
        self.sort_result.on('end', self.contig_stat_run)
        self.sort_result.run()
        self.step.sort_result.finish()
        self.step.update()
        self.single_module.append(self.sort_result)

    def contig_stat_run(self):
        """
        汇总信息并统计
        :return:
        """
        self.get_contig_file()
        self.contig_stat.set_options({
            "contig_dir": self.work_dir + '/contig_dir',
            "assembly_stat": "assembly.stat",
            "min_contig": str(self.option("min_contig"))
        })
        self.step.contig_stat.start()
        self.step.update()
        # self.contig_stat.on('end', self.len_distribute_run)
        self.contig_stat.on('end', self.run_zip_length)
        self.contig_stat.run()
        self.step.contig_stat.finish()
        self.step.length_distribute.start()
        self.step.update()

    def run_zip_length(self):
        """
        并行运行zip和length_distribute
        :return:
        """
        self.zip_run()
        self.len_distribute_run()
        for tool in self.zip_length:
            tool.run()
        self.on_rely(self.zip_length, self.set_output)
        #self.step.length_distribute.finish()
        #self.step.update()

    def zip_run(self):
        """
        对单拼各样品或者混拼总样品进行压缩
        :return:
        """
        if self.option('method') == 'multiple' and self.option("use_newbler") == "True":
            opts = {
                "file_path": self.sort_result.output_dir + '/Newbler_Mix.contig.fa'
            }
        else:
            opts = {
                "file_dir": self.work_dir + '/contig_dir'
            }
        self.yasuo.set_options(opts)
        self.zip_length.append(self.yasuo)

    def len_distribute_run(self):
        """
        长度分布
        :return:
        """
        self.len_distribute.set_options({
            "fasta_dir": self.contig_stat.output_dir,
            "len_range": "200,400,500,600,800,1000,1500,2000,3000",
        })
        self.zip_length.append(self.len_distribute)

    def get_list(self):
        """
        根据QC路径下list.txt，将文件信息转换成字典
        :return: dic:file_dic[sample_name][file_type]
        """
        file_dic = dict()
        self.logger.info("QC_dir路径：{}".format(self.option('QC_dir').path))
        ab_rout = self.option('QC_dir').prop['path'] + '/list.txt'
        with open(ab_rout, 'r') as list_file:
            for line in list_file:
                info = line.strip('\r\n').split('\t')
                name = info[1]
                type = info[2]
                if type not in ['l', 's', 'r']:
                    raise OptionError('质控样品的类型错误，必须为l/r/s之一', code="21300106")
                if name in file_dic.keys():
                    if type in file_dic[name].keys():
                        raise OptionError('质控list表中包含重复的样品及其pse类型，请检查质控list.txt ', code="21300107")
                    else:
                        file_dic[name][type] = info[0]
                else:
                    file_dic[name] = {type: info[0]}
        return file_dic

    def get_mem(self, type, base_number, tool="idba"):
        """
        from biocluster.config import Config
        sanger_type = self.parent._sheet.output.split(':')[0]
        sanger_prefix = Config().get_netdata_config(sanger_type)
        config_path = os.path.join(sanger_prefix[sanger_type + '_path'], "rerewrweset/metag/pipeline.config")
        # config_path = "/mnt/ilustre/tsanger-data/rerewrweset/metag/pipeline.config"  # module 测试用 add by GHD @20180313
        f = open(config_path, "r")
        content = f.readlines()
        for line in content[1:]:
            if self.parent._sheet.id in line and "assem" in line:
                tmp = line.strip().split("\t")
                mem = int(float(tmp[2]))
                split_n = 1
                return mem, split_n
        """
        # do follow lines if can't match task_id
        mem, split_n = self.calculate_mem(type, base_number)
        if tool == "idba":
            return mem, split_n
        elif tool == "megahit":
            return mem / 2, split_n

    def calculate_mem(self, type, base_number):
        """
        根据样品类型和此样品测序量，计算拼接所用内存, idba计算方法，megahit要除2
        :param type: 样品类型
        :param base_number: 样品测序量
        :return: int:mem int:split_n
        """
        type_coefficient = {
            "human": "0.6",
            "human gut": "0.6",
            "gut": "0.6",
            "soil": "1.5",
            "water": "1",
            "sludge": "1.2",
        }
        mem_base_ratio = 10
        if type in type_coefficient.keys():
            mem = float(base_number) / 1000000000 * float(type_coefficient[type]) * mem_base_ratio
            split_n = int(mem / 250) + 1
        else:
            mem = 250
            split_n = 1
        if mem > 250:
            mem = 250
        else:
            mem = int(mem)
        if float(base_number) / 1000000000 > 15:  # modified @ 20180309
            mem = 200
        else:
            mem = 100
        return mem, split_n
        # return 100, split_n

    # def get_dic(self, stat_file, insert_file):
    def get_dic(self):
        """
        根据reads_stat和insert_size 获得每个样本的最大读长和平均插入片段长度的字典
        :return:
        """
        insert_size = dict()
        base_number = dict()
        raw_read_len = dict()
        samp_type = dict()
        '''
        for one in sample_data:
            insert_size[one['sample_name']] = one['insert_size']
            raw_read_len[one['sample_name']] = one['raw_read_len']
            samp_type[one['sample_name']] = one['sample_source']
            if one['type'] == 'clean':
                base_number[one['sample_name']] = one['clean_base']
            elif one['type'] == 'optimised':
                base_number[one['sample_name']] = one['opt_base']
            else:
                raise OptionError("拼接前样品集type必须为clean或optimised")
        return raw_read_len, base_number, insert_size, samp_type

    '''
        with open(self.option('raw_stat').prop['path']) as fr1:
            lines = fr1.readlines()
            for line in lines[1:]:
                line_list = line.strip().split('\t')
                sample_name = line_list[0]
                samp_type[sample_name], insert_size[sample_name], raw_read_len[sample_name] = line_list[1:4]
        with open(self.option('qc_stat').prop['path']) as fr2:
            lines = fr2.readlines()
            for line in lines[1:]:
                line_list = line.strip().split("\t")
                sample_name = line_list[0]
                base = line_list[2]
                base_number[sample_name] = base
        return raw_read_len, base_number, insert_size, samp_type

    def get_contig_file(self):
        """
        将所有拼接结果文件整理成一个文件夹，为统计做准备
        :return:
        """
        if os.path.exists(self.work_dir + '/contig_dir'):
            shutil.rmtree(self.work_dir + '/contig_dir')
        os.mkdir(self.work_dir + '/contig_dir')
        for module in self.single_module:
            all_files = os.listdir(module.output_dir)
            for files in all_files:
                if files.endswith('.contig.fa'):
                    fa_file = os.path.join(module.output_dir, files)
                    link_fa_file = os.path.join(self.work_dir + '/contig_dir', files)
                    if os.path.exists(link_fa_file):
                        os.remove(link_fa_file)
                    os.link(fa_file, link_fa_file)
                elif files.endswith('.scaf.fa'):
                    fa_file = os.path.join(module.output_dir, files)
                    link_fa_file = os.path.join(self.work_dir + '/contig_dir',"Megahit_Mix.contig.fa")
                    if os.path.exists(link_fa_file):
                        os.remove(link_fa_file)
                    os.link(fa_file, link_fa_file)
                else:
                    pass

    def get_fq_file(self):
        """
        将所有未map的reads整理成一个文件夹
        :return:
        """
        if not os.path.exists(self.work_dir + '/unmap_dir'):
            os.mkdir(self.work_dir + '/unmap_dir')
        new_list = open(self.work_dir + '/unmap_dir/list.txt', 'w')
        for module in self.extract_fq_module:
            all_files = os.listdir(module.output_dir)
            list_file = os.path.join(module.output_dir, 'list.txt')
            with open(list_file, 'r') as lf:
                for line in lf:
                    new_list.write(line)
            for files in all_files:
                if files.endswith('fq'):
                    fa_file = os.path.join(module.output_dir, files)
                    link_fa_file = os.path.join(self.work_dir + '/unmap_dir', files)
                    if os.path.exists(link_fa_file):
                        os.remove(link_fa_file)
                    os.link(fa_file, link_fa_file)
                else:
                    pass

    def get_fq_sort_file(self):
        """
        将所有未map的reads整理成多个文件夹（根据样本大小100G进行判断）
        :return:
        """
        self.get_fq_file()
        new_number = 0
        new_num = 0
        file_size = {}
        file_numbler = 1
        with open(self.work_dir + '/unmap_dir/list.txt', "r") as f, open(self.work_dir + '/size_list.txt', 'w') as w:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split("\t")
                if re.search(r'1.fq', line[0]):
                    file = os.path.join(self.work_dir + '/unmap_dir/', line[0])
                    number = os.path.getsize(file)
                    file_size[line[0]] = number
                    new_number += number
                    w.write('{}\t{}\t{}\n'.format(file, line[1], number)) #统计左端序列的大小，并写入新文件size_list.txt中

        n = (float(new_number) / (100*1024*1024*1024)) #可以分几份
        n = math.ceil(n)
        n = int(n)
        if n <= file_numbler:
            n = file_numbler
        self.logger.info('%s'% n)
        os.system('sort -n -r -k 3 {} > {}'.format(self.work_dir + '/size_list.txt', self.work_dir + '/sort_size.txt')) # 按照左端大小进行排序，排除掉样本差异比较大造成切分文件夹个数错误
        with open(self.work_dir + '/sort_size.txt', 'r') as f1:
            for line in f1:
                line = line.strip().split("\t")
                num = int(line[2])
                new_num += num
                if new_num > (100*1024*1024*1024):# 如果大于100G则重新连接到新的文件夹
                    new_num = num
                    file_numbler += 1
                if os.path.exists(self.work_dir + '/sort_size_%s' % str(file_numbler)):
                    pass
                else:
                    os.mkdir(self.work_dir + '/sort_size_%s' % str(file_numbler))
                if not os.path.exists(self.work_dir + '/sort_size_%s/%s' % (str(file_numbler), str(line[1])+".l.fq")):
                    os.link(line[0], self.work_dir + '/sort_size_%s/%s' % (str(file_numbler), str(line[1])+".l.fq"))
                    os.link(self.work_dir + '/unmap_dir/%s' % (str(line[1])+".2.fq"), self.work_dir + '/sort_size_%s/%s' % (str(file_numbler), str(line[1])+".r.fq"))
        for i in range(1, n+1):
            if not os.path.exists(os.path.join(self.work_dir + '/sort_size_%s' % i, 'list.txt')):
                self.get_file_list(self.work_dir + '/sort_size_%s' % i)
        return n

    def get_file_list(self, path):
        """
        用于获得指定路径下的list文件
        :param path: 指定路径
        :return:
        """
        self.logger.info("指定文件夹的路径:{}".format(path))
        all_files = os.listdir(path)
        list_file = os.path.join(path, 'list.txt')
        with open(list_file, 'w') as w:
            for files in all_files:
                files_path = os.path.join(path, files)
                files_name = os.path.basename(files_path)
                if files.endswith('l.fq'):
                    name = files_name.split('.')[0]
                    type = 'l'
                elif files.endswith('r.fq'):
                    name = files_name.split('.')[0]
                    type = 'r'
                else:
                    name =""
                    type = ""
                w.write('{}\t{}\t{}\n'.format(files_name, name, type))

    def run(self):
        """
        运行
        :return:
        """
        super(MgAssIdbaModule, self).run()
        if self.option('assemble_tool') == 'idba':
            self.idba_run()
        elif self.option('assemble_tool') == 'megahit':
            self.megahit_run()
        elif self.option('assemble_tool') == 'automatic' and self.option('method') == 'multiple':
            raw_rd_len, base_num, insert_dic, sample_type = self.get_dic()
            choose_assemble_tool = 'idba'
            for key in insert_dic.keys():
                if float(base_num[key]) / 1000000000 > 20:
                    choose_assemble_tool = 'megahit'
                    break
            if choose_assemble_tool == 'idba':
                self.idba_run()
            elif choose_assemble_tool == 'megahit':
                self.megahit_run()
                # if self.option('method') == 'multiple':
                # self.bowtie2_run()
                # self.cat_reads_run()
                # self.mix_assem_run()
                # self.cut_length_run()
                # self.newbler_run()
                # self.sort_idba_result_run()
        # self.contig_stat_run()

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
                # os.link(oldfiles[i], newdir)
                oldfile_basename = os.path.basename(oldfiles[i])
                self.linkdir(oldfiles[i], os.path.join(newdir, oldfile_basename))

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        os.mkdir(self.output_dir + '/predict')
        stat_file = '/assembly.stat'
        if os.path.exists(self.output_dir + stat_file):
            os.remove(self.output_dir + stat_file)
        os.link(self.contig_stat.work_dir + stat_file, self.output_dir + stat_file)  # 如果重运行，contig_stat.output_dir下无此文件
        if os.path.exists(self.contig_stat.output_dir + stat_file):
            os.remove(self.contig_stat.output_dir + stat_file)  # 不用newbler条件下，contig_stat路径拷贝到predict下不需要此文件
        self.linkdir(self.yasuo.output_dir, self.output_dir)
        # if self.option('method') == 'simple':
        if self.option("use_newbler") == 'False':  # 不用newbler进行此种复制
            # self.linkdir(self.contig_stat.output_dir, self.output_dir)
            self.linkdir(self.contig_stat.output_dir, self.output_dir + '/predict')
            # self.option('contig').set_path(self.output_dir)
            # self.option('contig', self.output_dir)
        elif self.option('method') == 'multiple':
            # stat_file = '/assembly.stat'
            # # new_mix_file = '/Newbler_Mix.contig.fa'
            # if os.path.exists(self.output_dir + stat_file):
            #     os.remove(self.output_dir + stat_file)
            # self.linkdir(self.sort_result.output_dir, self.output_dir)
            self.linkdir(self.sort_result.output_dir + '/predict', self.output_dir + '/predict')
            # os.link(self.contig_stat.output_dir + stat_file, self.output_dir + stat_file)
            # self.option('contig').set_path(self.output_dir + '/predict')
            # self.option('conitg', self.output_dir + '/predict')
        self.option('contig').set_path(self.output_dir + '/predict')
        self.linkdir(self.len_distribute.output_dir, self.output_dir + '/len_distribute')
        self.logger.info("设置结果目录")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(MgAssIdbaModule, self).end()
