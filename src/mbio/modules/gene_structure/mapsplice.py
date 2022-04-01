# -*- coding: utf-8 -*-

# __author__ = 'linfang.jin'
# time: 2017/1/25 9:04

import importlib
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.core.exceptions import FileError
from biocluster.module import Module
from mbio.packages.ref_rna.trans_step import *
import re
from mbio.files.sequence.file_sample import FileSampleFile
from mbio.files.ref_rna.gene_structure.mapsplice_fa_dir import MapspliceFaDirFile


class MapspliceModule(Module):
    def __init__(self, work_id):
        super(MapspliceModule, self).__init__(work_id)
        options = [
            {"name": "have_single_seq", "type": "int", "default": 0},  ##只在0或1  0代表还没有单序列fasta文件夹 1则代表有
            {"name": "have_index", "type": "int", "default": 0},  ##只在0或1  0代表还没有bowtie1 index 1则代表有
            {"name": "single_fa_dir", "type": "infile", "format": "ref_rna.gene_structure.mapsplice_fa_dir"},
            {"name": "single_fa_index_dir", "type": "string", "default": ""},
            # split fasta 参数
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},
            # bowtie index 参数
            {"name": "bowtie_version", "type": "string", "default": "bowtie1"},
            # mapsplice map option
            {"name": "seq_type", "type": "string", "default": "paired"},  # ‘paired’ 或者 ‘single’
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # fastq文件夹
            {"name": "threads", "type": "int", "default": 10},
            # {"name": "ref_gtf", "type": "infile", "format": "ref_rna.assembly.gtf"},  # 一定要设置
            {"name": "segments_len", "type": "int", "default": 25},  # 限定在[18,25]之间
            {"name": "min_map_len", "type": "int", "default": 50},
            {"name": "min_intron_len", "type": "int", "default": 50},
            {"name": "max_intron_len", "type": "int", "default": 300000},
            {"name": "non_cano_double_anchor", "type": "int", "default": 1},  # 0:不设定，1：设定
            {"name": "non_cano_single_anchor", "type": "int", "default": 1},  # 0:不设定，1：设定
            {"name": "splice_mismatch_max_num", "type": "int", "default": 1},  # 必须在[0,2]之间
            {"name": "max_append_mismatch", "type": "int", "default": 3},
            {"name": "max_insertion_len", "type": "int", "default": 3},  # 范围为[0,10]
            {"name": "max_deletion_len", "type": "int", "default": 6},
            {"name": "double_anchor_option", "type": "int", "default": 0},
            {"name": "single_anchor_option", "type": "int", "default": 0},
            {"name": "fusion_option", "type": "int", "default": 0},
            {"name": "fusion_non_canonical", "type": "int", "default": 0},
            {"name": "min_fusion_distance", "type": "int", "default": 10000}
        ]
        self.add_option(options)
        self.ref_link = ""
        self.index_tools = []
        self.samples = {}
        self.mapsplice_map_tools = []
        self.ref_single_seq_fa = os.path.join(self.output_dir, "ref_single_seq_fa")
        self.ref_seq_fa_index = os.path.join(self.output_dir, "ref_seq_fa_index")

    def check_options(self):
        if self.option('have_single_seq') not in (0, 1):
            raise OptionError('have_single_seq 必须设置为0或1')
        if self.option('have_index') not in (0, 1):
            raise OptionError('have_index 必须设置为0或1')
        if self.option('have_single_seq') == 0:
            if not self.option('ref_fa'):
                raise OptionError("必须设置输入的fasta文件")
            if self.option('bowtie_version') != 'bowtie1':
                raise OptionError("bowtie_version只使用bowtie1做索引")
        if self.option('have_single_seq') == 1 and self.option('have_index') == 0:
            if not self.option('single_fa_dir'):
                raise OptionError("必须设置输入的single_fa_dir")
            if self.option('bowtie_version') != 'bowtie1':
                raise OptionError("bowtie_version只使用bowtie1做索引")
        if not self.option('fastq_dir'):
            raise OptionError('必须设置fastq文件夹')
        if not self.option('seq_type'):
            raise OptionError("必须设置测序类型： paired 或者 single")
        if self.option("fastq_dir").is_set:
            # self.logger.info(self.samples)
            list_path = os.path.join(self.option("fastq_dir").path, "list.txt")
            row_num = len(open(list_path, "r").readline().split())
            self.logger.info(row_num)
            if self.option('seq_type') == "paired" and row_num != 3:
                raise OptionError("paired序列list文件应该包括文件名、样本名和左右端说明三列")
            elif self.option('seq_type') == "single" and row_num != 2:
                raise OptionError("single序列list文件应该包括文件名、样本名两列")
        # if not self.option("fastq_dir").is_set and self.option('seq_type') in ["paired"]:
        #     if self.option("single_end_reads").is_set:
        #         raise OptionError("您上传的是单端测序的序列，请上传双端序列")
        #     elif not (self.option("left_reads").is_set and self.option("right_reads").is_set):
        #         raise OptionError("您漏了某端序列")
        # if not self.option("fastq_dir").is_set and self.option('seq_type') == "single":
        #     if not self.option("single_end_reads").is_set:
        #         raise OptionError("请上传单端序列")
        #     elif self.option("left_reads").is_set or self.option("right_reads").is_set:
        #         raise OptionError("有单端的序列就够啦")

        if self.option('segments_len') < 18 or self.option('segments_len') > 25:
            raise OptionError('reads 被分为的小片段（segments）长度应在18bp 和25bp 之间')
        if self.option('splice_mismatch_max_num') < 0 or self.option('splice_mismatch_max_num') > 2:
            raise OptionError('splice_mismatch_max_num 长度应在0bp 和2bp 之间')
        if self.option('min_intron_len') > self.option('max_intron_len'):
            raise OptionError('splice junctions 的最小长度不应该大于splice junction的最大长度')
        if self.option('max_insertion_len') > 10 or self.option('max_insertion_len') < 0:
            raise OptionError('max_insertion_len 应介于[0,10]之间')
        if self.option('double_anchor_option') not in (0, 1):
            raise OptionError('double_anchor_option 必须设置为0或1')
        if self.option('single_anchor_option') not in (0, 1):
            raise OptionError('single_anchor_option 必须设置为0或1')
        if self.option('fusion_option') not in (0, 1):
            raise OptionError(' fusion_option必须设置为0或1')
        if self.option('fusion_non_canonical') not in (0, 1):
            raise OptionError(' fusion_non_canonical 必须设置为0或1')
        return True

    def run(self):
        self.split_ref_fasta_run()
        super(MapspliceModule, self).run()

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def split_ref_fasta_run(self):
        if self.option('have_single_seq') == 0:
            self.link_ref()
            self.step.add_steps('mapsplice_split_fa')
            self.step.mapsplice_split_fa.start()
            self.step.update()
            self.logger.info("模块开始启动mapsplice_split_fa程序")
            self.mapsplice_split_fa = self.add_tool("ref_rna.gene_structure.mapsplice_split_fa")
            self.mapsplice_split_fa.set_options({
                "ref_fa": self.ref_link
            })
            """绑定下一个将要运行的步骤"""
            self.mapsplice_split_fa.on("end", self.finish_update, 'mapsplice_split_fa')
            self.mapsplice_split_fa.on('end', self.bowtie_index_run)
            self.mapsplice_split_fa.run()
        else:
            self.bowtie_index_run()

    def bowtie_index_run(self):
        if self.option('have_index') == 0:
            mapsplice_fa_dir_obj = MapspliceFaDirFile()
            if self.option('have_single_seq') == 0:
                mapsplice_fa_dir = self.mapsplice_split_fa.output_dir
                mapsplice_fa_dir_obj.set_path(mapsplice_fa_dir)
            if self.option('have_single_seq') == 1:
                mapsplice_fa_dir_obj = self.option('single_fa_dir')
            if mapsplice_fa_dir_obj.check():
                fa_lst = [os.path.join(mapsplice_fa_dir_obj.path, f) for f in os.listdir(mapsplice_fa_dir_obj.path) if
                          re.match(r'^(\S+)\.fa$', f.strip())]
                n = 0
                self.logger.info("模块开始启动bowtie_index程序")
                for fa in fa_lst:
                    n = n + 1
                    self.step.add_steps('bowtie_index_{}'.format(n))
                    step = getattr(self.step, 'bowtie_index_{}'.format(n))
                    step.start()
                    self.step.update()
                    bowtie_index = self.add_tool('ref_rna.gene_structure.bowtie_index')
                    self.logger.info("模块开始启动第{} bowtie_index程序".format(n))
                    bowtie_index.set_options({
                        'ref_fa': fa,
                        'bowtie_version': self.option("bowtie_version")
                    })

                    """绑定下一个将要运行的步骤"""
                    bowtie_index.on('end', self.finish_update, 'bowtie_index_{}'.format(n))
                    self.index_tools.append(bowtie_index)
                if len(self.index_tools) == 1:
                    self.bowtie_index.on('end', self.set_index_output, 'single_bowtie_index')
                    self.bowtie_index.on('end', self.set_step, {'end': self.step.bowtie_index})
                else:
                    self.on_rely(self.index_tools, self.set_index_output, 'index_tools')
                for tool in self.index_tools:
                    tool.run()

        if self.option('have_index') == 1:
            self.set_index_output()

    def set_index_output(self):
        self.logger.info("模块开始启动set index output程序")
        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)
        if not os.path.isdir(self.ref_single_seq_fa):
            os.mkdir(self.ref_single_seq_fa)
        if not os.path.isdir(self.ref_seq_fa_index):
            os.mkdir(self.ref_seq_fa_index)
        single_fa_dir = ""
        if self.option('have_single_seq') == 0:
            single_fa_dir = self.mapsplice_split_fa.output_dir
        if self.option('have_single_seq') == 1:
            single_fa_dir = self.option('single_fa_dir').path
        self.linkdir(single_fa_dir, os.path.basename(self.ref_single_seq_fa), r'^\S+.fa$')
        if self.option('have_index') == 1:
            index_dir = self.option('single_fa_index_dir')
            self.linkdir(index_dir, os.path.basename(self.ref_seq_fa_index), r'^\S+.ebwt$')
        if self.option('have_index') == 0:
            for tool in self.index_tools:
                self.linkdir(tool.output_dir, os.path.basename(self.ref_seq_fa_index), r'^\S+.ebwt$')
        self.check_fa_index_coherence()

    def check_fa_index_coherence(self):
        self.logger.info("模块开始检查fasta和index文件的一致性")
        index_file_set = set(os.listdir(self.ref_seq_fa_index))
        for fa_name in os.listdir(self.ref_single_seq_fa):
            m_fa = re.match(r'^(\S+\.fa)$', fa_name)
            if m_fa:
                seq_name = m_fa.group(1)
                self.logger.info("检查序列{}与其index的一致性".format(seq_name))
                ebwt1 = '{}.1.ebwt'.format(seq_name)
                ebwt2 = '{}.2.ebwt'.format(seq_name)
                ebwt3 = '{}.3.ebwt'.format(seq_name)
                ebwt4 = '{}.4.ebwt'.format(seq_name)
                rev_ebwt1 = '{}.rev.1.ebwt'.format(seq_name)
                rev_ebwt2 = '{}.rev.2.ebwt'.format(seq_name)
                condition = (ebwt1 in index_file_set) and (ebwt2 in index_file_set) and (ebwt3 in index_file_set) and (
                    ebwt4 in index_file_set) and (rev_ebwt1 in index_file_set) and (rev_ebwt2 in index_file_set)
                if not condition:
                    raise FileError('file {} in the fasta dir {} does not have legal(enough) index files')
            else:
                raise RuntimeError(
                    'the file name {} in path {} is an illegal single seq fasta file name!'.format(fa_name,
                                                                                                   self.ref_single_seq_fa))
        # self.step.add_steps('mapsplice_map')
        # self.step.mapsplice_map.start()
        # self.step.update()
        self.logger.info("模块结束检查fasta和index文件的一致性")
        self.mapsplice_map_run()

    def mapsplice_map_run(self):
        self.logger.info("模块开始运行mapslice map 程序")
        list_path = os.path.join(self.option("fastq_dir").path, "list.txt")
        self.logger.info("此次的样本文件路径为{}".format(list_path))
        file_sample = FileSampleFile()
        file_sample.set_path(list_path)
        # sample_list_dic = file_sample.get_sample_str()
        # sample_l_str = ",".join(
        #     os.path.join(self.option('fastq_dir').path, sample) for sample in sample_list_dic['left'])
        # sample_r_str = ",".join(
        #     os.path.join(self.option('fastq_dir').path, sample) for sample in sample_list_dic['right'])
        # sample_s_str = ",".join(
        #     os.path.join(self.option('fastq_dir').path, sample) for sample in sample_list_dic['single'])
        sample_list_dic = file_sample.get_list()
        map_options = {"ref_dir": self.ref_single_seq_fa,
                       "seq_type" : self.option('seq_type'),
                       "bowtie_index_dir": self.ref_seq_fa_index,
                       "threads": self.option('threads'), #"ref_gtf": self.option('ref_gtf'),
                       "segments_len": self.option('segments_len'),
                       "min_map_len": self.option('min_map_len'), "min_intron_len": self.option('min_intron_len'),
                       "max_intron_len": self.option('max_intron_len'),
                       "non_cano_double_anchor": self.option('non_cano_double_anchor'),
                       "non_cano_single_anchor": self.option('non_cano_single_anchor'),
                       "splice_mismatch_max_num": self.option('splice_mismatch_max_num'),
                       "max_append_mismatch": self.option('max_append_mismatch'),
                       "max_insertion_len": self.option('max_insertion_len'),
                       "max_deletion_len": self.option('max_deletion_len'),
                       "double_anchor_option": self.option('double_anchor_option'),
                       "single_anchor_option": self.option('single_anchor_option'),
                       "fusion_option": self.option('fusion_option'),
                       "fusion_non_canonical": self.option('fusion_non_canonical'),
                       "min_fusion_distance": self.option('min_fusion_distance')
                       }

        n = 0
        for sample in sample_list_dic.keys():
            n = n + 1
            if self.option('seq_type') == 'paired':
                if sample_list_dic[sample].has_key('r') and sample_list_dic[sample].has_key('l'):
                    map_options["_1_fq"] = os.path.join(self.option("fastq_dir").path, sample_list_dic[sample]['l'])
                    map_options["_2_fq"] = os.path.join(self.option("fastq_dir").path, sample_list_dic[sample]['r'])
                else:
                    raise KeyError("{} is not a normal paired fastq dic.".format(sample_list_dic[sample]))
                    # map_options["_1_fq"] = sample_l_str
                    # map_options["_2_fq"] = sample_r_str
            if self.option('seq_type') == 'single':
                map_options["_1_fq"] = os.path.join(self.option("fastq_dir").path, sample_list_dic[sample])
                # map_options["_1_fq"] = sample_s_str
            mapsplice_map = self.add_tool("ref_rna.gene_structure.mapsplice_map")
            self.step.add_steps('mapsplice_map_{}'.format(n))
            self.logger.info("此次map的参数为{}".format(map_options))
            mapsplice_map.set_options(map_options)
            step = getattr(self.step, 'mapsplice_map_{}'.format(n))
            step.start()
            """绑定下一个将要运行的步骤"""
            mapsplice_map.on('end', self.finish_update, 'mapsplice_map_{}'.format(n))
            self.mapsplice_map_tools.append(mapsplice_map)
        if len(self.mapsplice_map_tools) == 1:
            self.mapsplice_map.on('end', self.set_map_output, 'mapsplice_map')
        else:
            self.on_rely(self.mapsplice_map_tools, self.set_map_output, 'mapsplice_map')
        for tool in self.mapsplice_map_tools:
            tool.run()

    def link_ref(self):
        ref_fa = self.option('ref_fa').path
        self.ref_link = os.path.join(self.work_dir, os.path.basename(ref_fa))
        self.logger.info(self.ref_link)
        if os.path.exists(self.ref_link):
            os.remove(self.ref_link)
        os.symlink(ref_fa, self.ref_link)

    def set_map_output(self):
        self.logger.info('set map output')
        for map_tool in self.mapsplice_map_tools:
            self.linkdir(map_tool.output_dir, os.path.basename(map_tool.work_dir),r'^.*$')
        self.logger.info("set map output done")
        self.end()

    def linkdir(self, dirpath, dirname, pattern):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = [i for i in os.listdir(dirpath) if re.match(pattern, i)]
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        old_files = [os.path.join(dirpath, i) for i in allfiles]
        new_files = [os.path.join(newdir, i) for i in allfiles]
        for new_file in new_files:
            if os.path.exists(new_file):
                if os.path.isfile(new_file):
                    os.remove(new_file)
                else:
                    os.removedirs(new_file)
        for i in range(len(allfiles)):
            os.symlink(old_files[i], new_files[i])

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [r'.', '', 'mapsplice 运行结果']
        ]
        )
        result_dir.add_regexp_rules([
            ["map_result", '', 'mapsplice 运行结果'],
        ])
        super(MapspliceModule, self).end()
