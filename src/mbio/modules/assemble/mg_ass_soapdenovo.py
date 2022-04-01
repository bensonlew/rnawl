# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import os
import re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from biocluster.config import Config
from bson.objectid import ObjectId


class MgAssSoapdenovoModule(Module):
    """
    宏基因运用soapdenovo2组装
    author: guhaidong
    last_modify: 2017.09.15
    """
    def __init__(self, work_id):
        super(MgAssSoapdenovoModule, self).__init__(work_id)
        options = [
            # {"name": "data_id", "type": "string"},  # 主表任务ID
            {"name": "raw_stat", "type": "infile", "format": "sequence.profile_table"},  # 输入文件，原始序列统计文件
            {"name": "QC_dir", "type": "infile", "format": "sequence.fastq_dir"},  # 输入文件，质控或去宿主后的文件夹
            {"name": "qc_stat", "type": "infile", "format": "sequence.profile_table"},  # 输入文件，质控或去宿主后的统计结果
            #{"name": "reads_stat", "type": "string"},  # read最大读长,质控后的统计文件
            #{"name": "insert_size", "type": "string"},  # 平均插入片段长度
            {"name": "reverse_seq", "type": "string", "default": "0"},   # 配置文件的其他参数
            {"name": "asm_flags", "type": "string", "default": "3"},  # 配置文件的其他参数
            {"name": "rank", "type": "string", "default": "1"},  # 配置文件的其他参数
            {"name": "min_contig", "type": "int", "default": 500},  # 输入最短contig长度，默认500
            # {"name": "scafSeq", "type": "outfile", "format": "sequence.fasta"},  # 输出文件,sample.scafSeq
            # {"name": "scaftig", "type": "outfile", "format": "sequence.fasta"},  # 输出文件，scaffold去掉N后的序列
            {"name": "contig", "type": "outfile", "format": "sequence.fasta"},
            # 输出文件，去掉小于最短contig长度的序列
        ]
        self.add_option(options)
        self.qc_file = {}  # 质控数据信息
        self.sum_tools = []
        self.tools = []
        self.single_module = []
        self.step.add_steps("SOAPdenovo2", "contig_stat","length_distribute")
        self.kmer_list = ['39', '43', '47']
        self.zip_length = []

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
            raise OptionError('必须输入原始序列统计文件', code="21300201")
        if not self.option('QC_dir'):
            raise OptionError('必须输入质控后的fq文件夹', code="21300202")
        if not self.option('qc_stat'):
            raise OptionError('必须输入质控后的统计文件', code="21300203")
        #if not self.option('reads_stat'):
        #    raise OptionError('必须输入质控后的统计文件，需要read最大读长')
        #if not self.option('insert_size'):
        #    raise OptionError('必须输入平均插入片段的长度')
        return True

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def soapdenovo2_run(self):
        """
        SOAPdenovo组装
        :return:
        """
        n = 0
        # db = Config().mongo_client.tsanger_metagenomic
        # collection = db['mg_data_stat']
        # collection = db['data_stat']
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
        # self.logger.info(self.option('min_contig'))
        for key in insert_dic.keys():
            assem_mem = self.get_mem(sample_type[key], base_num[key]) # 计算运行内存
            # self.logger.info('type is ' + sample_type[key] + '; base_num is ' + base_num[key] + " assem_mem is " + str(assem_mem) + '\n')
            for kmer in self.kmer_list:
                self.SOAPdenovo2 = self.add_module('assemble.single_soap_denovo')
                self.step.add_steps('SOAPdenovo2_{}'.format(n))
                opts = ({
                    "fastq1": self.option('QC_dir').prop['path'] + '/' + self.qc_file[key]['l'],
                    "fastq2": self.option('QC_dir').prop['path'] + '/' + self.qc_file[key]['r'],
                    #"fastqs": self.option('QC_dir').prop['path']+'/' + key + '.sickle.s.fastq',
                    "sample_name": key,
                    "mem": 100 if assem_mem < 100 else assem_mem,
                    "max_rd_len": raw_rd_len[key],
                    "insert_size": insert_dic[key],
                    "reverse_seq": self.option('reverse_seq'),
                    "asm_flags": self.option('asm_flags'),
                    "rank": self.option('rank'),
                    "kmer": kmer,
                    "min_contig": str(self.option('min_contig'))
                })
                if 's' in self.qc_file[key].keys():
                    opts['fastqs'] = self.option('QC_dir').prop['path'] + '/' + self.qc_file[key]['s']
                self.SOAPdenovo2.set_options(opts)
                step = getattr(self.step, 'SOAPdenovo2_{}'.format(n))
                step.start()
                self.SOAPdenovo2.on('end', self.finish_update, 'SOAPdenovo2_{}'.format(n))
                self.single_module.append(self.SOAPdenovo2)
                self.sum_tools.append(self.SOAPdenovo2)
                n += 1
        if len(self.single_module) == 1:
            self.single_module[0].on('end', self.contig_stat_run)
        else:
            self.on_rely(self.single_module, self.contig_stat_run)
            self.step.contig_stat.start()
            self.step.update()
        for module in self.single_module:
            module.run()

    def contig_stat_run(self):
        """
        汇总信息并统计
        :return:
        """
        kmer = ",".join(self.kmer_list)
        self.get_file()
        self.contig_stat = self.add_tool("assemble.contig_stat")
        self.contig_stat.set_options({
            "contig_dir": self.work_dir + '/scaftig_dir',
            "choose_kmer": kmer,
            "assembly_stat": "assembly.stat",
        })
        self.contig_stat.on('end', self.run_zip_length)
        self.contig_stat.run()
        self.step.contig_stat.finish()
        self.step.length_distribute.start()
        self.step.update()

    def zip_run(self):
        """
        对单拼各样品或者混拼总样品进行压缩
        :return:
        """
        self.yasuo = self.add_tool('sequence.zip')
        self.yasuo.set_options({
                "file_dir": self.contig_stat.output_dir
            })
        #self.yasuo.on('end', self.len_distribute_run)
        self.zip_length.append(self.yasuo)

    def len_distribute_run(self):
        """
        长度分布
        :return:
        """
        self.len_distribute = self.add_tool("sequence.length_distribute")
        self.len_distribute.set_options({
            #"fasta_dir": self.work_dir + '/scaftig_dir',
            "fasta_dir": self.contig_stat.output_dir,
            "len_range": "200,400,500,600,800,1000,1500,2000,3000",
        })
        #self.len_distribute.on('end', self.run_zip_length)
        self.zip_length.append(self.len_distribute)

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

    def get_list(self):
        """
        根据QC路径下list.txt，将文件信息转换成字典
        :return: dic:file_dic[sample_name][file_type]
        """
        file_dic = dict()
        ab_rout = self.option('QC_dir').prop['path'] + '/list.txt'
        with open(ab_rout, 'r') as list_file:
            for line in list_file:
                info = line.split('\t')
                name = info[1]
                type = info[2].split('\n')[0]
                if type not in ['l', 's', 'r']:
                    raise OptionError('质控样品的类型错误，必须为l/r/s之一', code="21300204")
                if name in file_dic.keys():
                    if type in file_dic[name].keys():
                        raise OptionError('质控list表中包含重复的样品及其pse类型，请检查质控list.txt ', code="21300205")
                    else:
                        file_dic[name][type] = info[0]
                else:
                    file_dic[name] = {type: info[0]}
        return file_dic

    def get_mem(self, type, base_number):
        """
        from biocluster.config import Config
        sanger_type = self.parent._sheet.output.split(':')[0]
        sanger_prefix = Config().get_netdata_config(sanger_type)
        config_path = os.path.join(sanger_prefix[sanger_type + '_path'], "rerewrweset/metag/pipeline.config")
        f = open(config_path, "r")
        content = f.readlines()
        for line in content[1:]:
            if self.parent._sheet.id in line and "assem" in line:
                tmp = line.strip().split("\t")
                mem = int(float(tmp[2]))
                return mem
        """
        # do follow lines if can't match task_id
        mem = self.calculate_mem(type, base_number)
        return mem

    def calculate_mem(self,type,base_number):
        """
        根据样品类型和此样品测序量，计算拼接所用内存
        :param type: 样品类型
        :param base_number: 样品测序量
        :return: int:mem
        """
        type_coefficient = {
            "human": "0.6",
            "human gut": "0.6",
            "soil": "1.5",
            "water": "1",
            "sludge": "1.2"
        }
        mem_base_ratio = 10
        if type in type_coefficient.keys():
            mem = float(base_number) / 1000000000 * float(type_coefficient[type]) * mem_base_ratio
        else:
            mem = 250
        if mem > 250:
            mem = 250
        else:
            mem = int(mem)
        return mem

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

    def get_file(self):
        """
        将所有拼接结果文件整理成一个文件夹，为统计做准备
        :return:
        """
        if not os.path.exists(self.work_dir + '/scaftig_dir'):
            os.mkdir(self.work_dir + '/scaftig_dir')
        for module in self.single_module:
            all_files = os.listdir(module.output_dir)
            for files in all_files:
                if files.endswith('.contig.fa'):
                    fa_file = os.path.join(module.output_dir, files)
                    link_fa_file = os.path.join(self.work_dir + '/scaftig_dir', files)
                    if os.path.exists(link_fa_file):
                        os.remove(link_fa_file)
                    os.link(fa_file, link_fa_file)
                else:
                    pass

    def run(self):
        """
        运行
        :return:
        """
        super(MgAssSoapdenovoModule, self).run()
        self.soapdenovo2_run()

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
        stat_path = self.contig_stat.output_dir
        stat_file = '/assembly.stat'
        self.linkdir(self.yasuo.output_dir, self.output_dir)
        if os.path.exists(self.output_dir + stat_file):
            os.remove(self.output_dir + stat_file)
        if os.path.exists(self.output_dir + stat_file + '.tar.gz'):
            os.remove(self.output_dir + stat_file + '.tar.gz')
        os.link(self.contig_stat.work_dir + stat_file, self.output_dir + stat_file)
        if os.path.exists(stat_path + stat_file):
            os.remove(stat_path + stat_file)
        self.linkdir(self.contig_stat.output_dir, self.output_dir + '/predict')
        self.linkdir(self.len_distribute.output_dir, self.output_dir + '/len_distribute')
        self.logger.info("设置结果目录")
        self.option('contig').set_path(self.output_dir + '/predict')
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(MgAssSoapdenovoModule, self).end()