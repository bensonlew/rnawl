# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/3/10 18:08


import re, Bio, urllib2, regex, os
import subprocess
from biocluster.iofile import File
from collections import defaultdict
from biocluster.config import Config
from biocluster.core.exceptions import FileError

'''
gtf:gene transefer format
格式说明地址：http://mblab.wustl.edu/GTF22.html


'''


def check_seq_type(seq_type):
    return True


def dict_factory():
    return defaultdict(dict_factory)

def get_gene_list(gtf):
    gene_list = {}
    with open(gtf, "r") as file:
        for line in file:
            line = line.strip()
            tmp = line.split("\t")
            # print tmp[-1]
            m = re.match("transcript_id \"(.+)\";\sgene_id \"(.+?)\";", tmp[-1])
            if m:
                transcript_id = m.group(1)
                gene_id = m.group(2)
                if transcript_id not in gene_list.keys():
                    gene_list[transcript_id] = gene_id
    return gene_list


class GtfFile(File):
    # def __init__(self, fasta):  # fasta应为FastaFile对象
    def __init__(self):  # fasta应为FastaFile对象
        super(GtfFile, self).__init__()
        self._validate_gtf_tool = 'validate_gtf.pl'  # 此脚本
        # self._co_fasta = fasta
        self._contig_info = {}
        self._txpt_gene = {}
        self._txpt_len = {}
        self._cds_len = {}
        self.gtf2bed_path = Config().SOFTWARE_DIR + "/bioinfo/rna/scripts/gtf2bed.pl"
        # self._check_log_file = ''
        # self._structure_hierachy = dict_factory()

    def check(self):
        super(GtfFile, self).check()
        if self.prop["path"].endswith("gtf"):
            return True
        # self.__check_skechy()  # edited by shijin on 20170721 不适于第三列为gene的情况
        # self._check_chars()

    def check_format(self, fasta, so_file):

        pass
    
    def __check_skechy(self):
        '''

        粗略检查: 检查各个字段的字符是否符合规范
        1. tab分隔为9列
        2. 第九列 两个必须有的：gene_id value  ranscript_id value
        :return:
        '''
        for line in open(self.path):
            comment_m = regex.match(r'^#.+', line.strip())
            content_m = regex.match(
                r'^([^#]\S*?)\t+((\S+)\t+){7,7}((\btranscript_id\b|\bgene_id\b)\s+?\"(\S+?)\");.*((\btranscript_id\b|\bgene_id\b)\s+?\"(\S+?)\");(.*;)*$',
                line.strip())
            if content_m:
                if not {content_m.captures(5)[0], content_m.captures(8)[0]} == {'transcript_id', 'gene_id'}:
                    raise FileError('line error: %s 第9列必须有转录本id和基因id记录', variables = (line.strip()), code = "43702501")
                continue
            if not (comment_m or content_m):
                raise FileError(
                    'line %s is illegal in gtf file %s: it is not comment line(start with #) or tab-delimeted 9 colomuns line(the No9 line must contain gene_id txptid ) ', variables = (
                        line.strip(), self.path), code = "43702502")
    
    def check_in_detail(self, check_log_file):
        self.__check_skechy()
        # self.__check_hierachy()
        self.__check_gtf_bio_logic(check_log_file)
    
    def __check_gtf_bio_logic(self, log_file):
        '''
        此方法使用validate_gtf.pl文件检查gtf以下内容
        :return:
        '''
        if self._validate_gtf_tool:
            # tmp_out_txpt = os.path.join(os.path.dirname(self.path), os.path.basename(self.path) + '_tmp_txpt.gtf')
            validate_gtf_cmd = 'perl {} -fsm  {}'.format(self._validate_gtf_tool, self.path)
            open(log_file, 'wb').write(subprocess.check_output(validate_gtf_cmd, shell=True))
        else:
            raise FileError('gtf文件错误', code = "43702503")
    
    def _check_chars(self, merged=False):
        '''
        基本检查: 检查各个字段的字符是否符合规范
        1. tab分隔为9列
        2. 第九列 两个必须有的：gene_id value  ranscript_id value
        3. 每一列符合他们应有的规范
        :return:
        '''
        # gene_txpt_exon_dic = defaultdict(dict)
        for line in open(self.path):
            comment_m = regex.match(r'^#.+', line.strip())
            content_m = regex.match(
                r'^([^#]\S*?)\t+((\S+)\t+){7}((.*;)*((transcript_id|gene_id)\s+?\"(\S+?)\");.*((transcript_id|gene_id)\s+?\"(\S+?)\");(.*;)*)$',
                line.strip())
            if not (comment_m or content_m):
                raise FileError("%s", variables = (line), code = "43702504")
            if content_m:
                contig = content_m.captures(1)[0]
                seq_type = content_m.captures(2)[1].strip()
                start = content_m.captures(2)[2].strip()
                end = content_m.captures(2)[3].strip()
                frame = content_m.captures(2)[6].strip()
                strand = content_m.captures(2)[5].strip()
                contig_m = regex.match(r'^[\w.:^*$@!+?-|]+$', contig)  # contig的字符必须在[\w.:^*$@!+?-|]之内
                seq_type_m = check_seq_type(seq_type)  # seq_type必须在SO term集合之内
                start_m = regex.match(r'^\d+$', start)
                end_m = regex.match(r'^\d+$', end)
                frame_m = regex.match(r'^[\.120]$', frame)
                strand_m = regex.match(r'^[\.\?\-\+]$', strand)
                desc = content_m.captures(4)[0]

                if merged:
                    merged_m = regex.match(r'^.+?class_code "\w";$', desc)
                    if not merged_m:
                        raise FileError('illegal merged gtf', code = "43702505")
                if not (contig_m and seq_type_m and start_m and frame_m and end_m and strand_m):
                    raise FileError('line %s in gtf file %s is not legal.', variables = (line.strip(), self.path), code = "43702506")
    
    def __check_hierachy(self):
        '''
        包含__check_chars的功能
        :return:
        '''
        for line in open(self.path):
            comment_m = regex.match(r'^#.+', line.strip())
            content_m = regex.match(
                r'^([^#]\S*?)\t+((\S+)\t+){7}(.*;)*((transcript_id|gene_id)\s+?\"(\S+?)\");.*((transcript_id|gene_id)\s+?\"(\S+?)\");(.*;)*$',
                line.strip())

            if not (comment_m or content_m or regex.match(r'^$', line.strip())):
                raise FileError(
                    'line %s is illegal in gtf file %s: it is not comment line(start with #) or tab-delimeted 9 colomuns line(the No9 line must contain gene_id txptid ) ', variables = (
                        line.strip(), self.path), code = "43702507")

            if content_m:
                contig = content_m.captures(3)[0]
                seq_type = content_m.captures(3)[1]
                start = content_m.captures(3)[2]
                end = content_m.captures(3)[3]
                frame = content_m.captures(3)[5]
                strand = content_m.captures(3)[6]
                contig_m = regex.match(r'^[\w.:^*$@!+?-|]+$', contig)  # contig的字符必须在[\w.:^*$@!+?-|]之内
                seq_type_m = check_seq_type(seq_type)  #
                start_m = regex.match(r'^\d+$', start)
                end_m = regex.match(r'^\d+$', end)
                frame_m = regex.match(r'^[\.120]$', frame)
                strand_m = regex.match(r'^[\.\?\-\+]$', strand)
                if not (contig_m and seq_type_m and start_m and frame_m and end_m and strand_m):
                    raise FileError('line %s in gtf file %s is not legal.', variables = (line.strip(), self.path), code = "43702508")

                    # self._structure_hierachy[contig][]
    
    def check_gtf_for_merge(self):
        '''
        检查merged.gtf每一行的第九列
        :return:
        '''
        self._check_chars(merged=True)

    def to_isoform2unigene(self, t2g_file, g2t2p=None):
        '''
        转换为转录本基因对应关系表
        新基因返回最长转录本为yes已知基因最长已知转录本为yes
        已知基因的新转录本都为no
        return: 转录本\t基因\tyes_no\t长度\t蛋白
        '''
        t2g = self.get_txpt_gene_dic()
        t2l = self.get_txpt_len_dic()
        for i in t2g.keys():
            if t2l.has_key(i):
                pass
            elif i in self._cds_len:
                # 使用cds长度代替exon长度
                t2l[i] = self._cds_len[i]
            else:
                t2l[i] = 0

        t2g2l =[(i, t2g[i], t2l[i]) for i in t2g.keys()]
        t2g2l_sort = sorted(t2g2l, key= lambda x:x[2], reverse=True)
        gene_list = []

        if g2t2p:
            with open(g2t2p, 'rb') as f1:
                t2p = [(i.strip().split("\t")[1], i.strip().split("\t")[2]) for i in f1.readlines() if len(i.strip().split("\t")) >= 3]
            t2p_dict = dict(t2p)

        with open(t2g_file, 'wb') as f:
            for tran, gene, length in t2g2l_sort:
                if g2t2p:
                    if t2p_dict.has_key(tran):
                        pep = t2p_dict[tran]
                    else:
                        pep = ""
                else:
                    pep = tran

                if gene.startswith("MSTRG") or gene.startswith("XLOC") or gene.startswith("XLOC"):
                    if gene in gene_list:
                        f.write("{}\t{}\t{}\t{}\t{}\n".format(tran, gene, "no", length, pep))
                    else:
                        f.write("{}\t{}\t{}\t{}\t{}\n".format(tran, gene, "yes", length, pep))
                        gene_list.append(gene)
                elif not (tran.startswith("MSTRG") or tran.startswith("XLOC") or tran.startswith("TCONS")):
                    if gene in gene_list:
                        f.write("{}\t{}\t{}\t{}\t{}\n".format(tran, gene, "no", length, pep))
                    else:
                        f.write("{}\t{}\t{}\t{}\t{}\n".format(tran, gene, "yes", length, pep))
                        gene_list.append(gene)
                else:
                    f.write("{}\t{}\t{}\t{}\t{}\n".format(tran, gene, "no", length, pep))


    def get_txpt_len_dic(self):
        '''
        获取转录本长度信息
        liubinxu
        返回字典key:转录本， value:长度
        '''
        for line in open(self.path):
            txpt_id = ''
            gene_id = ''
            content_m = regex.match(
                r'^([^#][^\t]*?)\t+(([^\t]+)\t+){7}(.*;)*\s*((transcript_id|gene_id)\s+?\"(\S+?)\");.*((transcript_id|gene_id)\s+?\"(\S+?)\");(.*;)*',
                line.strip())
            if content_m:
                if 'transcript_id' in content_m.captures(6):
                    txpt_id = content_m.captures(7)[0]
                else:
                    txpt_id = content_m.captures(10)[0]
            gtf_line = line.strip().split("\t")
            if len(gtf_line) > 8 and gtf_line[2] == 'exon':
                exon_len = abs( int(gtf_line[4]) - int(gtf_line[3]) ) + 1
                if self._txpt_len.has_key(txpt_id):
                    self._txpt_len[txpt_id] += exon_len
                else:
                    self._txpt_len[txpt_id] = exon_len
            if len(gtf_line) > 8 and gtf_line[2].lower() == 'cds':
                cds_len = abs( int(gtf_line[4]) - int(gtf_line[3]) ) + 1
                if self._cds_len.has_key(txpt_id):
                    self._cds_len[txpt_id] += cds_len
                else:
                    self._cds_len[txpt_id] = cds_len
            else:
                pass
        return self._txpt_len

    def to_bed(self):
        bed_path = os.path.split(self.prop['path'])[0]
        bed = os.path.join(bed_path, os.path.split(self.prop['path'])[1] + ".bed")
        cmd = "perl {} {} > {}".format(self.gtf2bed_path, self.prop["path"], bed)
        try:
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError:
            os.remove(bed)
            raise FileError("运行出错", code = "43702509")
        pass

    def gtf_tbi(self):
        pass

    def get_txpt_gene_dic(self):
        '''
        获取转录本对应的基因
        返回字典key:转录本， value:基因
        liubinxu
        '''
        for line in open(self.path):
            txpt_id = ''
            gene_id = ''
            content_m = regex.match(
                r'^([^#][^\t]*?)\t+(([^\t]+)\t+){7}(.*;)*\s*((transcript_id|gene_id)\s+?\"(\S+?)\");.*((transcript_id|gene_id)\s+?\"(\S+?)\");(.*;)*',
                line.strip())
            if content_m:
                if 'transcript_id' in content_m.captures(6):
                    txpt_id = content_m.captures(7)[0]
                    gene_id = content_m.captures(10)[0]
                else:
                    txpt_id = content_m.captures(10)[0]
                    gene_id = content_m.captures(7)[0]
            if txpt_id:
                self._txpt_gene[txpt_id] = gene_id
        return self._txpt_gene

    def gtf_patch(self, gtf):
        """
        将gffread的输出文件中，只有转录本id，没有gene_id的行，加入gene_id
        :param gtf: gffread的输出文件
        :return:
        """
        gene_list = {}
        with open(gtf, "r") as file:
            for line in file:
                line = line.strip()
                tmp = line.split("\t")
                m = re.match("transcript_id \"(.+)\";\sgene_id \"(.+?)\";", tmp[-1])
                if m:
                    transcript_id = m.group(1)
                    gene_id = m.group(2)
                    if transcript_id not in gene_list.keys():
                        gene_list[transcript_id] = gene_id
        temp_gtf = os.path.split(self.prop["path"])[0] + "/tmp.gtf"
        w = open(temp_gtf, "w")
        with open(self.prop["path"], "r") as file:
            for line in file:
                line = line.strip()
                tmp = line.split("\t")
                if tmp[-1].find("gene_id") == -1:
                    m = re.match("transcript_id \"(.+?)\";", tmp[-1])
                    if m:
                        t_id = m.group(1)
                        if gene_id in gene_list.keys():
                            gene_id = gene_list[t_id]
                        elif re.match("transcript:(.+)", t_id):
                            t_id =  re.match("transcript:(.+)", t_id).group(1)
                            if t_id in gene_list.keys():
                                gene_id = gene_list[t_id]
                        else:
                            gene_id = t_id
                        lst = tmp[-1].split(";")
                        new_list = []
                        for item in lst:
                            new_list.append(item)
                            if item.startswith("transcript_id"):
                                new_list.append(" gene_id \"" + gene_id +"\"")
                        tmp[-1] = ";".join(new_list)
                        new_line = "\t".join(tmp)
                        w.write(new_line + "\n")
                elif tmp[-1].find("transcript_id") == -1:
                    self.logger.info("there is no transcript id in {}".format(line))
                else:
                    w.write(line + "\n")
        w.close()


if __name__ == "__main__":
    gtf = GtfFile()
    gtf.set_path("/mnt/ilustre/users/sanger-dev/workspace/20170410/Single_assembly_module_tophat_stringtie_zebra/Assembly/assembly_newtranscripts/merged.gtf")
    gtf.gtf_patch("/mnt/ilustre/users/sanger-dev/workspace/20170210/Refrna_refrna_test_01/FilecheckRef/Danio_rerio.GRCz10.85.gff3.gtf")
