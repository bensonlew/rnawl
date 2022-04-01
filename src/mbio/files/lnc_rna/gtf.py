# -*- coding: utf-8 -*-
# __author__ = fiona

import re
import Bio
import urllib2
import regex
import os
import subprocess
from biocluster.iofile import File
from collections import defaultdict
from biocluster.config import Config
from biocluster.core.exceptions import FileError

def check_seq_type(seq_type):
    return True

def dict_factory():
    return defaultdict(dict_factory)

def get_gene_list(gtf):
    gene_list = {}
    with open(gtf) as file:
        for line in file:
            line = line.strip()
            tmp = line.split("\t")
            m = re.match("transcript_id \"(.+)\";\sgene_id \"(.+?)\";", tmp[-1])
            if m:
                transcript_id = m.group(1)
                gene_id = m.group(2)
                if transcript_id not in gene_list.keys():
                    gene_list[transcript_id] = gene_id
    return gene_list

class GtfFile(File):
    def __init__(self):
        super(GtfFile, self).__init__()
        self._validate_gtf_tool = 'validate_gtf.pl'
        self._contig_info = {}
        self._txpt_gene = {}
        self._txpt_len = {}
        self._cds_len = {}
        self.gtf2bed_path = os.path.join(Config().PACKAGE_DIR, 'lnc_rna/gtf2bed.pl')
        self.gtf2standard = os.path.join(Config().PACKAGE_DIR, 'lnc_rna/gtf_standard.pl')

    def check(self):
        super(GtfFile, self).check()
        if self.prop["path"].endswith("gtf"):
            return True

    def check_format(self, fasta, so_file):
        pass

    def __check_skechy(self):
        for line in open(self.path):
            comment_m = regex.match(r'^#.+', line.strip())
            content_m = regex.match(
                r'^([^#]\S*?)\t+((\S+)\t+){7,7}((\btranscript_id\b|\bgene_id\b)\s+?\"(\S+?)\");.*((\btranscript_id\b|\bgene_id\b)\s+?\"(\S+?)\");(.*;)*$',
                line.strip())
            if content_m:
                if not {content_m.captures(5)[0], content_m.captures(8)[0]} == {'transcript_id', 'gene_id'}:
                    raise FileError('ERROR: Column 9 must have transcription id and gene id')
                continue
            if not (comment_m or content_m):
                raise FileError('line {} in {} is illegal.'.format(line, self.path))

    def check_in_detail(self, check_log_file):
        self.__check_skechy()
        self.__check_gtf_bio_logic(check_log_file)

    def __check_gtf_bio_logic(self, log_file):
        if self._validate_gtf_tool:
            validate_gtf_cmd = 'perl {} -fsm  {}'.format(self._validate_gtf_tool, self.path)
            open(log_file, 'wb').write(subprocess.check_output(validate_gtf_cmd, shell=True))
        else:
            raise FileError('GTF file error')

    def _check_chars(self, merged=False):
        for line in open(self.path):
            comment_m = regex.match(r'^#.+', line.strip())
            content_m = regex.match(
                r'^([^#]\S*?)\t+((\S+)\t+){7}((.*;)*((transcript_id|gene_id)\s+?\"(\S+?)\");.*((transcript_id|gene_id)\s+?\"(\S+?)\");(.*;)*)$',
                line.strip())
            if not (comment_m or content_m):
                raise FileError("%s", variables=(line), code="43702504")
            if content_m:
                contig = content_m.captures(1)[0]
                seq_type = content_m.captures(2)[1].strip()
                start = content_m.captures(2)[2].strip()
                end = content_m.captures(2)[3].strip()
                frame = content_m.captures(2)[6].strip()
                strand = content_m.captures(2)[5].strip()
                contig_m = regex.match(r'^[\w.:^*$@!+?-|]+$', contig)
                seq_type_m = check_seq_type(seq_type)
                start_m = regex.match(r'^\d+$', start)
                end_m = regex.match(r'^\d+$', end)
                frame_m = regex.match(r'^[\.120]$', frame)
                strand_m = regex.match(r'^[\.\?\-\+]$', strand)
                desc = content_m.captures(4)[0]
                if merged:
                    merged_m = regex.match(r'^.+?class_code "\w";$', desc)
                    if not merged_m:
                        raise FileError('illegal merged gtf', code="43702505")
                if not (contig_m and seq_type_m and start_m and frame_m and end_m and strand_m):
                    raise FileError('line {} in {} is illegal.'.format(line, self.path))

    def __check_hierachy(self):
        for line in open(self.path):
            comment_m = regex.match(r'^#.+', line.strip())
            content_m = regex.match(
                r'^([^#]\S*?)\t+((\S+)\t+){7}(.*;)*((transcript_id|gene_id)\s+?\"(\S+?)\");.*((transcript_id|gene_id)\s+?\"(\S+?)\");(.*;)*$',
                line.strip())
            if not (comment_m or content_m or regex.match(r'^$', line.strip())):
                raise FileError('line {} in {} is illegal.'.format(line, self.path))
            if content_m:
                contig = content_m.captures(3)[0]
                seq_type = content_m.captures(3)[1]
                start = content_m.captures(3)[2]
                end = content_m.captures(3)[3]
                frame = content_m.captures(3)[5]
                strand = content_m.captures(3)[6]
                contig_m = regex.match(r'^[\w.:^*$@!+?-|]+$', contig)
                seq_type_m = check_seq_type(seq_type)  #
                start_m = regex.match(r'^\d+$', start)
                end_m = regex.match(r'^\d+$', end)
                frame_m = regex.match(r'^[\.120]$', frame)
                strand_m = regex.match(r'^[\.\?\-\+]$', strand)
                if not (contig_m and seq_type_m and start_m and frame_m and end_m and strand_m):
                    raise FileError('line {} in {} is illegal.'.format(line, self.path))

    def check_gtf_for_merge(self):
        self._check_chars(merged=True)

    def to_isoform2unigene(self, t2g_file, g2t2p=None):
        t2g = self.get_txpt_gene_dic()
        t2l = self.get_txpt_len_dic()
        for i in t2g.keys():
            if t2l.has_key(i):
                pass
            elif i in self._cds_len:
                t2l[i] = self._cds_len[i]
            else:
                t2l[i] = 0
        t2g2l = [(i, t2g[i], t2l[i]) for i in t2g.keys()]
        t2g2l_sort = sorted(t2g2l, key=lambda x: x[2], reverse=True)
        gene_list = []
        if g2t2p:
            with open(g2t2p, 'rb') as f1:
                t2p = [(i.strip().split("\t")[1], i.strip().split("\t")[2]) for i in f1.readlines() if
                       len(i.strip().split("\t")) >= 3]
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
        for line in open(self.path):
            txpt_id = ''
            gene_id = ''
            content_m = regex.match(
                r'^([^#][^\t]*?)\t+(([^\t]+)\t+){7}(.*;)*\s*((transcript_id|gene_id)\s+?\"(\S+?)\");.*((transcript_id|gene_id)\s+?\"(\S+?)\");(.*;)*$',
                line.strip())
            if content_m:
                if 'transcript_id' in content_m.captures(6):
                    txpt_id = content_m.captures(7)[0]
                else:
                    txpt_id = content_m.captures(10)[0]
            gtf_line = line.strip().split("\t")
            if len(gtf_line) > 8 and gtf_line[2] == 'exon':
                exon_len = abs(int(gtf_line[4]) - int(gtf_line[3])) + 1
                if self._txpt_len.has_key(txpt_id):
                    self._txpt_len[txpt_id] += exon_len
                else:
                    self._txpt_len[txpt_id] = exon_len
            if len(gtf_line) > 8 and gtf_line[2].lower() == 'cds':
                cds_len = abs(int(gtf_line[4]) - int(gtf_line[3])) + 1
                if self._cds_len.has_key(txpt_id):
                    self._cds_len[txpt_id] += cds_len
                else:
                    self._cds_len[txpt_id] = cds_len
            else:
                pass
        return self._txpt_len

    def get_txpt_pos_dic(self, type="T"):
        pos = dict()
        for line in open(self.path):
            txpt_id = ''
            gene_id = ''
            content_m = regex.match(
                r'^([^#][^\t]*?)\t+(([^\t]+)\t+){7}(.*;)*\s*((transcript_id|gene_id)\s+?\"(\S+?)\");.*((transcript_id|gene_id)\s+?\"(\S+?)\");(.*;)*$',
                line.strip())
            if content_m:
                if 'transcript_id' in content_m.captures(6):
                    txpt_id = content_m.captures(7)[0]
                    gene_id = content_m.captures(10)[0]
                else:
                    txpt_id = content_m.captures(10)[0]
                    gene_id = content_m.captures(7)[0]
            gtf_line = line.strip().split("\t")
            if type=="G":
                txpt_id = gene_id
            if len(gtf_line) > 8 and gtf_line[2] == 'exon':
                if txpt_id in pos:
                    if pos[txpt_id]["start"] > int(gtf_line[3]):
                        pos[txpt_id]["start"] = int(gtf_line[3])
                    if pos[txpt_id]["end"] < int(gtf_line[4]):
                        pos[txpt_id]["end"] = int(gtf_line[4])
                else:
                    pos[txpt_id] = {
                        "chr": gtf_line[0],
                        "start": int(gtf_line[3]),
                        "end": int(gtf_line[4])
                    }
                    pos[txpt_id]["end"] = int(gtf_line[4])
            if len(gtf_line) > 8 and gtf_line[2].lower() == 'cds':
                if txpt_id in pos:
                    if pos[txpt_id]["start"] > int(gtf_line[3]):
                        pos[txpt_id]["start"] = int(gtf_line[3])
                    if pos[txpt_id]["end"] < int(gtf_line[4]):
                        pos[txpt_id]["end"] = int(gtf_line[4])
                else:
                    pos[txpt_id] = {
                        "chr": gtf_line[0],
                        "start": int(gtf_line[3]),
                        "end": int(gtf_line[4])
                    }
                    pos[txpt_id]["chr"] = gtf_line[0]
                    pos[txpt_id]["start"] = int(gtf_line[3])
                    pos[txpt_id]["end"] = int(gtf_line[4])
        return pos


    def filter_by_trans_list(self, trans_list, out_file):
        with open(self.path, 'r') as gtf_in, open(out_file, 'w') as gtf_out:
            for line in gtf_in:
                txpt_id = ''
                gene_id = ''
                content_m = regex.match(
                    r'^([^#][^\t]*?)\t+(([^\t]+)\t+){7}(.*;)*\s*((transcript_id|gene_id)\s+?\"(\S+?)\");.*((transcript_id|gene_id)\s+?\"(\S+?)\");(.*;)*$',
                    line.strip())
                if content_m:
                    if 'transcript_id' in content_m.captures(6):
                        txpt_id = content_m.captures(7)[0]
                    else:
                        txpt_id = content_m.captures(10)[0]
                    if txpt_id in trans_list:
                        gtf_out.write(line)



    def to_bed(self, bed=None):
        cmd = "perl {} --i {} >  {}".format(self.gtf2standard , self.prop["path"], self.prop["path"] + ".std")
        try:
            subprocess.check_output(cmd, shell=True)
            print cmd
        except subprocess.CalledProcessError:
            raise FileError("Operation error")
        pass


        bed_path = os.path.split(self.prop['path'])[0]
        if not bed:
            bed = os.path.join(bed_path, os.path.split(self.prop['path'])[1] + ".bed")
        cmd = "perl {} {} {}".format(self.gtf2bed_path, self.prop["path"] + ".std", bed)
        try:
            subprocess.check_output(cmd, shell=True)
            print cmd
        except subprocess.CalledProcessError:
            os.remove(bed)
            raise FileError("Operation error")
        pass

    def merge_gtf(self, add_gtf, out_gtf):
        with open(self.prop["path"], "r") as file_1, open(add_gtf, "r") as file_2, open(out_gtf, "w") as file_3:
            for line in file_1:
                file_3.write(line)
            for line in file_2:
                if line.startswith("#"):
                    pass
                else:
                    file_3.write(line)

    def gtf_tbi(self):
        pass

    def get_txpt_gene_dic(self):
        for line in open(self.path):
            txpt_id = ''
            gene_id = ''
            content_m = regex.match(
                r'^([^#][^\t]*?)\t+(([^\t]+)\t+){7}(.*;)*\s*((transcript_id|gene_id)\s+?\"(\S+?)\");.*((transcript_id|gene_id)\s+?\"(\S+?)\");(.*;)*$',
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
                            t_id = re.match("transcript:(.+)", t_id).group(1)
                            if t_id in gene_list.keys():
                                gene_id = gene_list[t_id]
                        else:
                            gene_id = t_id
                        lst = tmp[-1].split(";")
                        new_list = []
                        for item in lst:
                            new_list.append(item)
                            if item.startswith("transcript_id"):
                                new_list.append(" gene_id \"" + gene_id + "\"")
                        tmp[-1] = ";".join(new_list)
                        new_line = "\t".join(tmp)
                        w.write(new_line + "\n")
                elif tmp[-1].find("transcript_id") == -1:
                    self.logger.info("there is no transcript id in {}".format(line))
                else:
                    w.write(line + "\n")
        w.close()
