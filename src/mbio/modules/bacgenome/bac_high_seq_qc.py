# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict
class BacHighSeqQcModule(Module):
    """
    Fastp质控
    author: wangzhaoyue
    last_modify: 2017.12.13
    """
    def __init__(self, work_id):
        super(BacHighSeqQcModule, self).__init__(work_id)
        options = [
            {"name": "sample_path", "type": "infile", "format": "sequence.file_sample"},
            # fastq路径list.txt文件，第一列路径，第二列样本名，第三列序列类型 l or r
            {'name': "sample_info", "type": "infile", "format": "sequence.barcode_info"},  # 样本信息表
            {'name': 'rm_single', 'type': 'bool', 'default': False},  # 是否舍弃质控后的single端
            {'name': 'phix_tool', 'type': 'string', 'default': "bwa", "choose": ["bwa", "bowtie"]},  # 去phix的方法
            {'name': 'flag', 'type': "string", "default": "4"},  # 提取没有比对上的reads,此处固定取值4，软件默认0
            {'name': 'readl', "type": "string"},  # 切除序列的阈值
            {'name': 'illuminaclip', 'type': "string", "default": "2:30:10"},  # 2:30:10
            {'name': 'leading', 'type': "string", "default": "3"},  # 切除首端碱基质量小于0的碱基或者N
            {'name': 'tailing', 'type': "string", "default": "3"},  # 切除末端碱基质量小于20的碱基或者N
            {'name': 'sliding_window', 'type': "string", "default": "4:15"},
            # 例50:20  Windows的size是50个碱基，其平均碱基质量小于20，则切除
            {'name': 'minlen', 'type': "string", "default": "36"},  # 最低reads长度
            {"name": "seqprep_quality", "type": "string", "default": '20'},
            {"name": "seqprep_length", "type": "string", "default": '25'},
            {"name": "adapter_a", "type": "string", "default": "AGATCGGAAGAGCACACGTC"},
            {"name": "adapter_b", "type": "string", "default": "AGATCGGAAGAGCGTCGTGT"},
            {"name": "sickle_quality", "type": "string", "default": '20'},
            {"name": "sickle_length", "type": "string", "default": '20'},
            {"name": "qual_type", "type": "string", "default": 'sanger'},
            {"name": "clean_list", "type": "outfile", "format": "bacgenome.list_file"}
        ]
        self.sample_path = defaultdict(list)
        self.read_info = {}
        self.sample_info = {}
        self.modules = []
        self.add_option(options)
    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('sample_path'):
            raise OptionError('必须输入样本文件夹对应的路径信息', code="21400601")
        row_num = len(open(self.option("sample_path").prop['path'], "r").readline().split())
        if row_num != 3:
            raise OptionError("PE序列list文件应该包括文件名、样本名和左右端说明三列", code="21400602")
        if not self.option('sample_info'):
            raise OptionError('必须输入样本信息文件，包括样本，文库类型、插入片段长度三列', code="21400603")
        return True
    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()
    def single_microbial_genome_qc_run(self):
        base_dir = os.path.dirname(self.option("sample_path").prop['path'])
        n = 0
        for sample in self.sample_path:
            self.single_microbial_genome_qc = self.add_module('bacgenome.bac_qc')
            self.step.add_steps('bac_qc{}'.format(n))
            opts = {
                "fq1": base_dir + '/' + self.sample_path[sample][0],
                "fq2": base_dir + '/' + self.sample_path[sample][1],
                'rm_single': self.option("rm_single"),
                "phix_tool": self.option("phix_tool"),
                "sample_name": sample,
                "readl":self.read_info[sample],
                "insert_size": self.sample_info[sample],
                "flag": self.option('flag'),
                "illuminaclip": self.option('illuminaclip'),
                "leading": self.option('leading'),
                "tailing": self.option('tailing'),
                "sliding_window": self.option('sliding_window'),
                "minlen": self.option('minlen'),
                "seqprep_quality": self.option('seqprep_quality'),
                "seqprep_length": self.option('seqprep_length'),
                "adapter_a": self.option('adapter_a'),
                "adapter_b": self.option('adapter_b'),
                "sickle_quality": self.option('sickle_quality'),
                "sickle_length": self.option('sickle_length'),
                "qual_type": self.option('qual_type'),
            }
            self.single_microbial_genome_qc.set_options(opts)
            step = getattr(self.step, 'bac_qc{}'.format(n))
            step.start()
            self.step.update()
            self.single_microbial_genome_qc.on('end', self.finish_update, 'bac_qc{}'.format(n))
            self.modules.append(self.single_microbial_genome_qc)
            n += 1
        self.logger.info(self.modules)
        self.on_rely(self.modules, self.set_output)
        self.step.update()
        for module in self.modules:
            module.run()

    def run(self):
        """
        运行
        :return:
        """
        super(BacHighSeqQcModule, self).run()
        self.get_info()
        time.sleep(2)
        self.single_microbial_genome_qc_run()

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

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for module in self.modules:
            self.logger.info(module)
            self.linkdir(module.output_dir + '/sickle', self.output_dir)
        path =self.work_dir + "/list.txt"
        time.sleep(2)
        with open(path,'w') as f:
            files = os.listdir(self.output_dir)
            for file in files:
                if re.search(r'.fq$',file):
                    tmp = file.split('.')
                    if tmp[-2] == '1':
                        f.write(file + '\t' + tmp[0] + '\t' + 'l' + '\n')
                    elif tmp[-2] == '2':
                        f.write(file + '\t' + tmp[0] + '\t' + 'r' + '\n')
                    else:
                        f.write(file + '\t' + tmp[0] + '\t' + 's' + '\n')
        f.close()
        if os.path.exists(self.output_dir + "/list.txt"):
            os.remove(self.output_dir + "/list.txt")
        os.link(path,self.output_dir + "/list.txt")
        self.option('clean_list').set_path(self.output_dir + "/list.txt")
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
        super(BacHighSeqQcModule, self).end()
    def get_info(self):
        """
        获得样本对应的路径信息，以及样本的
        :return:
        """
        with open(self.option("sample_info").prop['path'])as f:
            lines = f.readlines()
            for line in lines[0:]:
                tmp = line.strip().split('\t')
                self.sample_info[tmp[0]] = tmp[1]  # 样本的插入片段长度
                self.read_info[tmp[0]] =tmp[2]  #read长度
        with open(self.option('sample_path').prop['path'])as fr:
            for line in fr:
                # self.logger.info(self.sample_path)
                tmp = line.strip().split('\t')
                if tmp[1] in self.sample_info.keys():
                    if tmp[1] in self.sample_path.keys():
                        if tmp[2] == 'l':
                            self.sample_path[tmp[1]].insert(0, tmp[0])
                        else:
                            self.sample_path[tmp[1]].append(tmp[0])
                    else:
                        self.sample_path[tmp[1]].append(tmp[0])
                else:
                    self.set_error('需要质控的序列样本%s没有相关的样本信息，请核实！', variables=(tmp[1]), code="21400601")
            for key in self.sample_path.keys():
                if len(self.sample_path[key]) > 2:
                    self.set_error('需要质控的序列样本%s有重名，请改样本名或分开质控！', variables=(key), code="21400602")
                elif len(self.sample_info[key]) < 2:
                    self.set_error('样本%s对应的R1,R2序列不全,请核实！', variables=(key), code="21400603")


