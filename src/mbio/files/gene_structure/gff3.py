# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/3/2 15:56
import re, os, regex, sys
#
# sys.path.append('/mnt/hgfs/F/code_lib/SangerBiocluster/src/biocluster')
# sys.path.append('/mnt/hgfs/F/code_lib/SangerBiocluster/src/biocluster/core')

import subprocess
from sequence_ontology import SequenceOntologyFile
from collections import defaultdict
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from biocluster.config import Config
from mbio.files.sequence.fasta import FastaFile

'''
检查gff标准：
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
seq 特征符合sequence ontology 规定 ：https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/master/releases/so-xp.owl/so-simple.obo 此文件会有更新 但url地址不变
sequence_ontology_file_url = 'https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/master/releases/so-xp.owl/so.obo'
'''


class Gff3File(File):
    def __init__(self):
        super(Gff3File, self).__init__()
        self._properties = {}
        self._contigs = set()
        self._pos_set = set()
        self._contig_boundary_dic = {}
        self._seq_type_set = set()
        self._strand_set = set()
        self._phase_set = set()
        self._column_number = 9
        self._contig_max_len_dic = dict()
        self._gtf = ''  # gff3的配套gtf路径
        self._gffread_path = ''
        self._seq_pos = dict()
        self.check_parent = False
        self._genes = {}
        self._mrnas = {}
        self._exons = {}
        self._cds = {}
        self._feature_tree = {}
        self._attrs_record_set = {}

    def check(self):
        super(Gff3File, self).check()
        if (not self.prop["path"].endswith("gff")) and (not self.prop["path"].endswith("gff3")):
            raise FileError("gff文件格式不正确", code="42200501")
            # self.check_logical()

    def check_format(self, fasta_file, so_file):
        # self.check_lines()
        self.parse_and_check_column(fasta_file, so_file)
        self.check_logical()
        return True

    def parse_and_check_column(self, fa_file, so_file):
        contig_cmd = 'grep -v \'^#\' %s |awk -F \'\\t\' \'{print $1}\' |uniq| sort  -n -r | uniq' % self.path
        self._contigs = set(
            [contig.strip() for contig in subprocess.check_output(contig_cmd, shell=True).strip().split('\n')])
        seq_type_cmd = 'grep -v \'^#\' %s |awk -F \'\\t\' \'{print $3}\' |uniq| sort  -nr | uniq' % self.path
        self._seq_type_set = set(
            [seq_type.strip() for seq_type in subprocess.check_output(seq_type_cmd, shell=True).strip().split('\n')])
        strand_cmd = 'grep -v \'^#\' %s |awk -F \'\\t\' \'{print $7}\' |uniq| sort  -n -r | uniq' % self.path
        self._strand_set = set(
            [strand.strip() for strand in subprocess.check_output(strand_cmd, shell=True).strip().split('\n')])
        phase_cmd = 'grep -v \'^#\' %s |awk -F \'\\t\' \'{print $8}\' |uniq| sort  -n -r | uniq' % self.path
        self._phase_set = set(
            [phase.strip() for phase in subprocess.check_output(phase_cmd, shell=True).strip().split('\n')])

        pos_cmd = 'grep -v \'^#\' %s |awk -F \'\\t\' \'{printf $1":"$4"-"$5"\\n"}\' |uniq |sort' % self.path
        self._pos_set = set(
            [pos.strip() for pos in subprocess.check_output(pos_cmd, shell=True).strip().split("\n")])

        attr_cmd = 'grep -v \'^#\' %s |awk -F \'\\t\' \'{print $9}\' |uniq |sort' % self.path
        self._attrs_record_set = set(
            [attrs.strip() for attrs in subprocess.check_output(attr_cmd, shell=True).strip().split('\n')])
        self.check_contigs()
        self.pos_check(fa_file)
        self.check_seq_types(so_file, 'sequence_feature')
        self.check_strand()
        self.check_phase()

    def check_logical(self):
        self.phase_logical_check()
        self.attrs_format_check()
        type_id_cmd = 'grep  \'^[^#].*ID=\' %s |awk -F \'\\t\' \'{printf $3"\t"$9"\\n"}\' | uniq' % self.path
        type_id_content = [record.strip() for record in
                           subprocess.check_output(type_id_cmd, shell=True).strip().split('\n')]
        self.type_id_check(type_id_content)

    def attrs_format_check(self):
        for attr_record in self._attrs_record_set:
            if not regex.match(r'^(.+?=.+?;)+.+?=.+?;*$', attr_record.strip()):
                raise FileError('column 9 must match the format: tag=value', code = "42200502")

    def type_id_check(self, content):
        type_id_dic = defaultdict(int)
        for line in content:
            m = regex.search(r'(\S+)\t.*?ID=([^;=\"\'\t,]+)?', line.strip())
            if m:
                type_id_dic[(m.group(1), m.group(2))] += 1
            else:
                continue
        for count in type_id_dic.values():
            if count > 1:
                raise FileError('ID must be Uniq among same seq types', code = "42200503")

    def parent_relation_check(self, content):
        type_relation_set = set()
        for record in content:
            m = regex.search(r'(\S+)\t.*?Parent=([^;=\"\'\t,:]+)?:([^;=\"\'\t,:]+)?;', record)
            type_relation_set.add((m.group(1), m.group(2)))
        for type_pair in type_relation_set:
            self.check_seq_types_relation(type_pair[0], type_pair[1])

    def check_contigs(self):
        for contig in self._contigs:
            if (not regex.match(r'^[a-zA-Z0-9\.:^*$@!+_?-|]+$', contig)) or contig.startswith('>'):
                raise FileError('contig 错误', code = "42200504")
            else:
                continue
        return True

    def pos_check(self, fa_path):
        self._contig_boundary_dic = dict.fromkeys(self._contigs, (0, 0))  # 注意gff文件的pos起点
        for pos in self._pos_set:
            [contig, start, end] = regex.split(r'[:-]', pos)
            if not (regex.match(r'^\d+$', start) and regex.match(r'^\d+$', end)):
                raise FileError('illegal pos format', code = "42200505")
            if int(start) > int(end):
                raise FileError('%s record illegal pos in gff file %s', variables = (pos, self.path), code = "42200506")
            end = int(end)
            if end > self._contig_boundary_dic[contig][1]:
                self._contig_boundary_dic[contig][1] = end
        self.set_fasta_file(fa_path)
        contig_len_dic = self._fasta.get_contig_len()
        if self._contigs != set(contig_len_dic.keys()):
            raise FileError('gff3 and fasta file contigs id set not agreed', code = "42200507")
        for contig in contig_len_dic.keys():
            if contig_len_dic[contig] < self._contig_boundary_dic[contig][0]:
                raise FileError('illogical seq length ', code = "42200508")

    def check_seq_types(self, target_term, so_file):
        (id_set, name_set) = self.get_so_set(so_file, target_term, 'is_a')
        illegal_term_set = self._seq_type_set - id_set.union(name_set)
        if illegal_term_set:
            raise FileError('illegal seq types', code = "42200509")
        return True

    def parse_directives(self):
        directives_cmd = 'grep  \'^##\' %s' % self.path
        directives_content = [record.strip() for record in
                              subprocess.check_output(directives_cmd, shell=True).strip().split('\n')]
        for record in directives_content:
            v_m = regex.match(r'^##gff-version\s+(\S+)$', record)
            if v_m:
                self._version = v_m.group(1)
                continue
            seq_region_m = regex.search(r'^##sequence-region\s+(\S+\s+)+$', record.strip() + ' ')
            if seq_region_m:
                self._contigs_declar = {seq_region_m.captures(1)[0].strip(): (
                    int(seq_region_m.captures(1)[1].strip()), int(seq_region_m.captures(1)[2].strip()))}
                continue

    def check_lines(self):
        for line in open(self.path):
            if regex.match(r'^#', line) or regex.match(r'^$', line.strip()):
                continue
            else:
                if not regex.match(r'^[^#]\S+\t(.+?\t){7}(.+?=.+?;)*(.+?=.+?)*', line.strip()):
                    raise FileError('有不合格的行', code = "42200510")

    def set_so_file(self, value):
        if not os.path.isfile(value):
            raise FileError('so file does not exist', code = "42200511")
        self._so_file = SequenceOntologyFile()
        self._so_file.set_path(value)
        return self._so_file

    def set_gtf_file(self, value):
        self._gtf = value
        pass

    def set_gffread_path(self, value):
        self._gffread_path = value

    def set_gtf2bed_path(self, value):
        self._gtf2bed_path = value

    def set_fasta_file(self, fa):
        self._fasta = FastaFile()
        self._fasta.set_path(fa)

    def get_so_set(self, so_f, target_term_id, relation):
        if not regex.match(r'^SO:\d+$', target_term_id):
            raise FileError('illegal input so id', code = "42200512")
        if not os.path.isfile(so_f):
            raise FileError('so file does not exist', code = "42200513")
        so_file = SequenceOntologyFile()
        so_file.set_path(so_f)
        so_file.parse()
        return so_file.findAll(target_term_id, relation)

    def check_strand(self):
        for strand in self._strand_set:
            if not regex.match(r'^[\.\?\-\+]$', strand):
                raise FileError('illegal strand value', code = "42200514")

    def check_phase(self):
        for phase in self._phase_set:
            if not regex.match(r'^[\.120]$', phase):
                raise FileError('illegal phase value', code = "42200515")

    def phase_logical_check(self):
        cds_phase_cmd = 'grep -v \'^#\' %s |awk -F \'\\t\' \'$3~/CDS|cds/{print $8}\' |uniq | sort |uniq' % self.path
        cds_phase_set = set(
            [phase.strip() for phase in subprocess.check_output(cds_phase_cmd, shell=True).strip().split('\n')])
        if cds_phase_set & {'0', '1', '2', '.'} != cds_phase_set:
            raise FileError('illogical phase value: %s', variables = (cds_phase_set), code = "42200516")

    def get_genbank_assembly_id(self):
        if self._parse_status:
            for item in self._build_info_dic.keys():
                info_macth = regex.match(r'NCBI:(\S+)', self._build_info_dic[item])
                if regex.match(r'.*genome-build-accession.*', item) and info_macth:
                    self._genbank_assembly_id = info_macth.group(1)
            return self._genbank_assembly_id

    def to_gtf(self):
        temp_gtf = os.path.join(os.path.dirname(self._gtf), os.path.basename(self._gtf).split('.')[0]) + '_temp.gtf'
        to_gtf_cmd = '%s %s -T -O -o %s  ' % (self._gffread_path, self.path, temp_gtf)
        try:
            subprocess.check_output(to_gtf_cmd, shell=True)
        except:
            os.remove(temp_gtf)
            raise FileError("运行出错", code = "42200517")
        gtf = open(self._gtf, 'wb')
        for line in open(temp_gtf):
            newline = regex.sub(r'"(\S+?):(\S+?)";', '"\g<2>";', line)
            tmp = newline.strip().split("\t")[-1]
            if "gene_id" not in tmp:
                m = re.match("transcript_id \"(.+)\";", tmp)
                if m:
                    gene_id = self.get_parent(m.group(1))
                    newline = newline.strip() + " gene_id \"" + gene_id + "\";\n"
            gtf.write(newline)
        gtf.close()

    def gtf_to_bed(self, gtf_path):
        """
        gtf格式转bed格式
        """
        bed_path = os.path.split(gtf_path)[0]
        bed = os.path.join(bed_path, os.path.split(gtf_path)[1] + ".bed")
        cmd = "python {} -i {} -o {}".format(self._gtf2bed_path, gtf_path, bed)
        try:
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError:
            # os.remove(bed)
            raise FileError("运行出错", code = "42200518")
        return True

    def get_parent(self, transcript):
        if not self.check_parent:
            self.get_transcript_dict()
        gene_id = self._genes[transcript]
        return gene_id

    def get_transcript_dict(self):
        if not self.check_parent:
            self.transcripts = []
            with open(self.prop["path"], "r") as file:
                for line in file:
                    tmp = line.strip().split("\t")[-1]
                    m = re.match("ID=transcript:(.+);Parent=gene:(.+?);", tmp)
                    if m:
                        if m.group(1) not in self.transcripts:
                            self.transcripts.append(m.group(1))
                            self._genes[m.group(1)] = m.group(2)
            self.check_parent = True

    def split_gff_by_name(self, output_dir):
        '''
        add by zouxuan 20180416
        通过gff第二列的location名称，切分gff文件，得到相应的文件
        :return: split_file
        '''
        split_file = {}
        with open(self.prop["path"], "r") as file:
            head = file.readline()
            for line in file:
                split_line = line.strip().split("\t")
                one = split_line[1].split("_")[0]
                if one in split_file:
                    with open(split_file[one], 'a+') as f:
                        f.write(line)
                else:
                    split_file[one] = output_dir + '/' + one + '.gff'
                    if os.path.exists(split_file[one]):
                        os.remove(split_file[one])
                    with open(split_file[one], 'a+') as f:
                        f.write(head)
                        f.write(line)
        return split_file


if __name__ == '__main__':
    '''
    1. trans gff3 to gtf example:
    gff3 = Gff3File()
    gff3.set_path('/mnt/hgfs/F/temp/Homo_sapiens.GRCh38.87.gff3')
    gff3.set_gtf_file('/mnt/hgfs/F/temp/Homo_sapiens.gtf')
    gff3.set_gffread_path('/home/linfang/app/cufflinks/gffread')
    gff3.to_gtf()

    2.

    gff3 = Gff3File()
    gff3.set_path('/mnt/ilustre/users/sanger-dev/workspace/20170210/Refrna_refrna_test_01/FilecheckRef/Danio_rerio.GRCz10.85.gff3')
    gff3.set_gtf_file('/mnt/ilustre/users/sanger-dev/x.gtf')
    gff3.set_gffread_path('/mnt/ilustre/users/sanger-dev/app/bioinfo/rna/cufflinks-2.2.1/gffread')
    gff3.to_gtf()
    # type_id_cmd = 'grep  \'^[^#].*ID=\' %s |awk -F \'\\t\' \'{printf $3"\t"$9"\\n"}\' | uniq' % gff3.path
    # type_id_content = [record.strip() for record in
    #                    subprocess.check_output(type_id_cmd, shell=True).strip().split('\n')]
    # gff3.type_id_check(type_id_content)
    '''
    gff3 = Gff3File()
    gff3.set_path(
        '/mnt/ilustre/users/sanger-dev/workspace/20170210/Refrna_refrna_test_01/FilecheckRef/Danio_rerio.GRCz10.85.gff3')
    fasta = ''
    so_file = ''
    # gff3.check_format(fasta_file=fasta, so_file=so_file)
    # gff3.check_logical()
    gff3.check()
