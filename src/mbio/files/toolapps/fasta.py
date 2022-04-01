# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
from biocluster.iofile import File
from collections import defaultdict
import re
import subprocess
from biocluster.config import Config
import os
from biocluster.core.exceptions import FileError
from Bio import SeqIO


class FastaFile(File):
    """
    定义Fasta文件， 需安装seqstat工具软件
    """
    def __init__(self):
        super(FastaFile, self).__init__()
        self.seqstat_path = os.path.join(Config().SOFTWARE_DIR, "bioinfo/seq/biosquid_1.9g+cvs20050121/bin/seqstat")
        # fasta与gff相关的属性
        self._seq_ids = []
        self._seq_obj = []
        self._seq_id_len_dic = {}
        self.seq_type = ''  # 有 prot nucl两种类型

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(FastaFile, self).get_info()
        seqinfo = self.get_seq_info()
        self.set_property("file_format", seqinfo[0])
        self.set_property("seq_type", seqinfo[1])
        self.set_property("seq_number", seqinfo[2])
        self.set_property("bases", seqinfo[3])
        self.set_property("longest", seqinfo[4])
        self.set_property("shortest", seqinfo[5])
        newfile = self.check_fasta_unzip()
        self.set_property("fasta", newfile)

    def check_fasta_unzip(self):
        """
        检查fasta是否压缩，如果为压缩
        :return:
        """
        fasta_name = os.path.basename(self.prop['path'])
        self.fasta = self.prop['path']
        if re.search(r'.tar.gz', self.fasta):
            self.new_fasta = (self.fasta).strip(".tar.gz")
            unzip_cmd = "tar -zxvf " + self.fasta
            try:
                subprocess.check_output(unzip_cmd, shell=True)
                return self.new_fasta
            except subprocess.CalledProcessError:
                raise FileError("解压运行出错！", code="44400601")
        elif re.search(r'.gz', self.fasta):
            new_fasta_name = fasta_name.split(".gz")[0]
            dir_path = self.fasta.split("%s"%(fasta_name))[0]
            self.new_fasta = os.path.join(dir_path, new_fasta_name)
            unzip_cmd = "zcat " + self.fasta + " > " + self.new_fasta
            try:
                subprocess.check_output(unzip_cmd, shell=True)
                return self.new_fasta
            except subprocess.CalledProcessError:
                raise FileError("解压运行出错！", code="44400602")
        else:
            self.new_fasta = self.fasta
            return self.new_fasta

    def check(self):
        """
        检测文件是否满足要求,发生错误时应该触发FileError异常
        :return:
        """
        # print self.prop
        self.get_info()
        if super(FastaFile, self).check():
            if self.prop['file_format'] != 'FASTA':
                raise FileError("文件格式错误", code="44400603")
            if self.prop["seq_number"] < 1:
                raise FileError("应该至少含有一条序列", code="44400604")
        return True

    def ncbi_blast_tool_check(self):
        """
        供ncbi.blast Tool检查
        :return:
        """
        if self.check():
            if self.prop['seq_type'] not in {"DNA", "Protein"}:
                raise FileError("不支持此类型的Fasta进行blast比对", code="44400605")
        return True

    def check_fasta_type(self):
        tmp_file = self.prop['path']+'_tmp_check_type'
        with open(self.prop['path']) as fr, open(tmp_file,'w') as fw:
            fw.write('>check_type\n')
            check_type_seq = ''
            for line in fr:
                if line[0] == '>':
                    continue
                else:
                    check_type_seq += line
                    if len(check_type_seq) > 1000:
                        break
            fw.write(check_type_seq)

        try:
            subpro = subprocess.check_output(self.seqstat_path + " " + tmp_file, shell=True)
            result = subpro.split('\n')
            seq_type = re.split(r':\s+', result[6])[1]
            os.remove(tmp_file)
            return seq_type

        except subprocess.CalledProcessError:
            raise FileError("seqstat 运行出错！", code="44400606")

    def get_seq_info(self):
        """
        获取Fasta信息
        :return: (format,seq_type,seq_number,bases,longest,shortest)
        """
        try:
            subpro = subprocess.check_output(self.seqstat_path + " " + self.prop['path'], shell=True)
            result = subpro.split('\n')
            fformat = re.split(r':\s+', result[5])[1]
            #seq_type = re.split(r':\s+', result[6])[1]
            seq_type = self.check_fasta_type()
            seq_number = re.split(r':\s+', result[7])[1]
            bases = re.split(r':\s+', result[8])[1]
            shortest = re.split(r':\s+', result[9])[1]
            longest = re.split(r':\s+', result[10])[1]
            # print (fformat, seq_type, seq_number, bases, longest, shortest)
            return fformat, seq_type, seq_number, bases, longest, shortest
        except subprocess.CalledProcessError:
            raise FileError("seqstat 运行出错！", code="44400607")

    def get_all_seq_name(self):
        seq_name = defaultdict(int)
        for seq in SeqIO.parse(self.prop["path"], "fasta"):
            seq_name[seq.id] += 1
        dup_list = list()
        for k in seq_name.iterkeys():
            if seq_name[k] > 1:
                dup_list.append(k)
            if len(dup_list) > 0:
                str_ = "; ".join(dup_list)
                raise FileError("序列名:%s在输入的fasta文件里面重复", code="44400608")
        return seq_name

    def check_trinity(self):
        """
        检查是不是trinity生成的组装结果文件，检查序列名的开头是不是TRINITY_
        """
        with open(self.path, 'r') as f:
            check_num = 10  # 最小检查序列数
            for i in f:
                if check_num == 0:
                    break
                if i[0] == '>':
                    check_num -= 1
                    if i.startswith('>TRINITY_'):
                        pass
                    else:
                        raise FileError('文件不是trinity生成的文件', code="44400609")
        return True


    def get_contig_len(self):
        '''
        author： jinlinfang
        date：20170405
        :return:
        '''
        for seq in SeqIO.parse(self.prop["path"], "fasta"):
            self._seq_id_len_dic[seq.id] =int(len(seq.seq))
        return self._seq_id_len_dic

    def split_single_seq(self, output_dir):
        '''
        author： jinlinfang
        date：20170125
        实现将目标fasta文件分割为数个文件，
        每个文件只包含其一条序列，且文件名为seq_name.fa,所有文件放在output_dir里
        :param output_dir:
        :return:
        '''
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        if self.check():
            try:
                seq_records = SeqIO.parse(self.path, 'fasta')
                for seq_record in seq_records:
                    seq_seq = seq_record.seq
                    seq_name = seq_record.name
                    line = '>{}\n{}\n'.format(seq_name, seq_seq)
                    open(os.path.join(output_dir, seq_name + '.fa'), 'w').write(line)
            except Exception:
                raise FileError("get split fa to single seqs failed", code="44400610")

    def split(self, output, chunk=10000):
        """
        拆分Fasta文件成最大chunk大小的快
        :param output:  String 输出目录
        :param chunk:  int 块大小
        :return:
        """
        s, n = 1, 0
        wf = open("%s/%s.fa" % (output, s), 'w')
        with open(self.prop['path'], 'r') as f:
            while 1:
                line = f.readline()
                if not line:
                    wf.close()
                    break
                re_id = re.compile(r'^>(\S+)')
                m_id = re_id.match(line)
                if m_id is not None:
                    n += 1
                    if n == chunk + 1:
                        wf.close()
                        s += 1
                        n = 0
                        wf = open("%s/%s\.fa" % (output, s), 'w')
                wf.write(line)

    def merge_all_seq(self, file_name):
        # '''
        # 实现将目标fasta文件合并为一条序列
        # :param out_file:
        # '''
        if self.check():
            try:
                with open (file_name+'.fasta','w') as wf:
                    wf.write('>'+os.path.basename(file_name)+'\n')
                    for seq in SeqIO.parse(self.prop['path'], 'fasta'):
                        wf.write(str(seq.seq))
            except Exception:
                raise FileError("merge seqs failed", code="44400611")
