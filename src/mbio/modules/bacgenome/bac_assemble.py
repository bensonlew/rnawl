# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict
class BacAssembleModule(Module):
    """
    微生物基因组组装
    author: gaohao
    last_modify: 2018.03.25
    """
    def __init__(self, work_id):
        super(BacAssembleModule, self).__init__(work_id)
        options = [
            {"name": "fq_dir", "type": "infile", "format": "sequence.fastq_dir"},  #
            {'name': 'sample_name', "type": "string"},  # 样本名
            {'name': 'seq_type', "type": "string"},  # 数据类型
            {'name': "sample_info", "type": "infile", "format": "sequence.barcode_info"},  # 样本信息表
            {"name": "scaffold", "type": "outfile", "format": "sequence.fasta"},  # 输出文件
            {"name": "genome_size", "type":"float", "default": 5000000},
            {"name": "depth", "type": "float", "default": 150},
            {'name': 'kmer', 'type': 'string', 'default': '21-47'},  # 质控软件选择
        ]
        self.add_option(options)
        self.step.add_steps('mp_cut','config_file','assemble','scaf_select','gapcloser')
        self.mp_cut = self.add_module('bacgenome.mp_cut_fq')
        self.config_file =self.add_tool('assemble.create_config')
        self.assemble = self.add_module('bacgenome.assemble_soap_denovo')
        self.scaf_select = self.add_tool('assemble.scaf_select')
        self.gapcloser = self.add_tool('assemble.gapcloser_scaf')
        self.sample_path = {}
        self.modules = []
        self.sample_info = {}
        self.file = self.work_dir + '/PE_list.txt'
        self.tk_file = self.work_dir + '/PE_list.tk.txt'
        self.seqtk_tools = []


    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('fq_dir').is_set:
            raise OptionError('必须输入clean data文件夹', code="21400501")
        if not self.option('sample_name'):
            raise OptionError('必须输入sample_name样品名称！', code="21400502")
        if not self.option('seq_type'):
            raise OptionError('必须输入序列文库类型名称！', code="21400503")
        else:
            if self.option('seq_type') not in ['PE','PE,MP','MP,PE']:
                raise OptionError('必须输入正确的文库类型名称！', code="21400504")
        if not self.option('sample_info').is_set:
            raise OptionError('必须输入sample_info文件！', code="21400505")


    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def run_mp_cut(self):
        opts = {
            "fq_dir": self.option('fq_dir'),
            "sample_info": self.option('sample_info'),
            "sample_name": self.option('sample_name'),
        }
        self.mp_cut.set_options(opts)
        self.mp_cut.run()
        self.step.mp_cut.finish()
        self.step.config_file.start()
        self.step.update()

    #抽取原始数据
    def run_seqtk(self):
        raw_new_map = self.get_seqtk_scale()
        for fq in raw_new_map:
            seqtk_fq = self.add_tool('bacgenome.seqtk')
            opts = {
                "fastq": fq,
                "outfastq" : raw_new_map[fq]['outfile'],
                "scale" : raw_new_map[fq]['scale']
            }
            seqtk_fq.set_options(opts)
            self.seqtk_tools.append(seqtk_fq)

    def run_config(self):
        opts = {
            "PE_list":self.tk_file,
            "MP_list":self.mp_cut.option('mp_list'),
            "sample_name": self.option('sample_name')
        }
        self.config_file.set_options(opts)
        self.config_file.run()
        self.step.config_file.finish()
        self.step.assemble.start()
        self.step.update()

    def run_assemble(self):
        opts = {
            "kmer": self.option('kmer'),
            "config": self.config_file.option('config_file'),
            "sample_name": self.option('sample_name'),
        }
        self.assemble.set_options(opts)
        self.assemble.run()
        self.step.assemble.finish()
        self.step.scaf_select.start()
        self.step.update()

    def run_select(self):
        opts = {
            "seq_dir": self.assemble.option('scafSeq'),
        }
        self.scaf_select.set_options(opts)
        self.scaf_select.run()
        self.step.scaf_select.finish()
        self.step.gapcloser.start()
        self.step.update()

    def run_gaploser(self):
        opts = {
            "seq_scaf": self.scaf_select.option('scf_seq'),
            "PE_list":self.file,
            "sample_name": self.option('sample_name'),
        }
        self.gapcloser.set_options(opts)
        self.gapcloser.on('end', self.set_output, 'gapcloser')
        self.gapcloser.run()
        self.step.gapcloser.finish()
        self.step.update()

    def run(self):
        """
        运行
        :return:
        """
        super(BacAssembleModule, self).run()
        self.get_info()
        self.get_list()
        time.sleep(2)
        self.run_seqtk()

        if self.option('seq_type') in ['PE']:
            self.gapcloser.on('end', self.end)
            self.on_rely(self.seqtk_tools, self.run_config)
            self.config_file.on('end', self.run_assemble)
            self.assemble.on('end', self.run_select)
            self.scaf_select.on('end', self.run_gaploser)
            #self.run_config()

        elif self.option('seq_type') in ['PE,MP','MP,PE']:
            self.gapcloser.on('end',self.end)
            self.on_rely(self.seqtk_tools, self.run_mp_cut)
            self.mp_cut.on('end',self.run_config)
            self.config_file.on('end', self.run_assemble)
            self.assemble.on('end', self.run_select)
            self.scaf_select.on('end', self.run_gaploser)
            #self.run_mp_cut()

        for seqtk_tool in self.seqtk_tools:
            seqtk_tool.run()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        newdir = os.path.join(self.work_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        allfiles = os.listdir(dirpath)
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
        if event["data"] == "gapcloser":
            self.linkdir(self.gapcloser.output_dir,self.output_dir + '/scf')
            self.option('scaffold').set_path(self.output_dir + '/scf/' + self.option('sample_name') + '.scaffold.fna')

    def end(self):
        super(BacAssembleModule, self).end()

    def get_info(self):
        """
        获得样本对应的路径信息，以及样本的
        :return:
        """
        with open(self.option('fq_dir').prop['path'] + '/list.txt')as fr:
            for line in fr:
                tmp = line.strip().split('\t')
                if tmp[1] in self.sample_path.keys():
                    self.sample_path[tmp[1]][tmp[2]] = tmp[0]
                else:
                    self.sample_path[tmp[1]] ={tmp[2]:tmp[0]}



    def get_list(self):
        """
        获得样本对应的路径信息，以及样本的
        :return:
        """
        self.logger.info(self.sample_path)
        with open(self.option("sample_info").prop['path'])as f:
            lines = f.readlines()
            for line in lines[0:]:
                tmp = line.strip().split('\t')
                self.sample_info[tmp[0]] = tmp[1] + '\t' + tmp[2]# 样本的插入片段长度
        base_dir = self.option("fq_dir").prop['path']
        with open(self.file, 'w') as f:
            for sample in self.sample_path:
                if re.search(r'PE', sample):
                    if sample in self.sample_info.keys():
                        if 's' not in self.sample_path[sample].keys():
                            f.write(sample + '\t' + base_dir + '/' + self.sample_path[sample]['l'] + '\t' + base_dir + '/' + self.sample_path[sample]['r'] + '\t-\t' + self.sample_info[sample] + '\n') #20200507
                        else:
                            f.write(sample + '\t' + base_dir + '/' + self.sample_path[sample]['l'] + '\t' + base_dir + '/' + self.sample_path[sample]['r'] + '\t' + base_dir + '/' + self.sample_path[sample]['s'] +'\t' + self.sample_info[sample] + '\n')
        f.close()

    #根据指定的基因组大小，深度，和read 长度 计算抽取的read 数
    def get_seqtk_scale(self):
        self.genome_size = self.option('genome_size')
        base_all = self.genome_size * self.option("depth")
        raw_new_map = {}
        dir_num = 0
        with open(self.file) as fr, open(self.tk_file, 'w') as fw:
            for line in fr:
                dir_num += 1
                spline = line.strip().split("\t")
                read_len = int(spline[5])
                read_num = int(base_all/read_len)
                raw_fq1_path = spline[1]
                raw_fq1 = os.path.basename(spline[1])
                out_dir = os.path.join(self.work_dir, str(dir_num))
                if not os.path.exists(out_dir):
                    os.mkdir(out_dir)
                spline[1] = os.path.join(out_dir, raw_fq1)
                raw_fq2_path = spline[2]
                raw_fq2 = os.path.basename(spline[2])
                spline[2] = os.path.join(out_dir, raw_fq2)
                raw_new_map[raw_fq1_path] = {'scale':3000000, 'outfile': spline[1]}
                raw_new_map[raw_fq2_path] = {'scale':3000000, 'outfile': spline[2]}
                fw.write('\t'.join(spline)+'\n')
        return raw_new_map








