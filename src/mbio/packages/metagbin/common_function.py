# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import pandas as pd
import os
import re
import shutil
from Bio import SeqIO
from Bio import Entrez
from biocluster.config import Config
from mbio.packages.taxon.accession2taxon import taxon as tx
from mbio.packages.annotation.mg_annotation.mg_taxon import mg_taxon


def link_dir(olddir, newdir):
    """
    hard link directory from olddir to newdir
    :param olddir:
    :param newdir:
    :return:
    """
    if not os.path.isdir(olddir):
        raise Exception("不存在路径: %s" % olddir)
    allfiles = os.listdir(olddir)
    if not os.path.exists(newdir):
        os.makedirs(newdir)
    newfiles = [os.path.join(newdir,i) for i in allfiles]
    for newfile in newfiles:
        if os.path.exists(newfile):
            if os.path.isfile(newfile):
                os.remove(newfile)
            elif os.path.isdir(newfile):
                shutil.rmtree(newfile)
    if len(allfiles) >= 1:
        for i in allfiles:
            if os.path.isfile(os.path.join(olddir, i)):
                os.link(os.path.join(olddir, i),os.path.join(newdir, i))
            elif os.path.isdir(os.path.join(olddir, i)):
                link_dir(os.path.join(olddir, i),os.path.join(newdir, i))
    else:
        raise Exception("结果文件夹为空:%s" % allfiles)

def link_tax_dir(olddir, newdir):
    if not os.path.isdir(olddir):
        raise Exception("不存在路径: %s" % olddir)
    allfiles = os.listdir(olddir)
    if not os.path.exists(newdir):
        os.makedirs(newdir)
    for file in allfiles:
        oldfile = os.path.join(olddir, file)
        newfile = os.path.join(newdir, file)
        if os.path.isdir(newfile):
            shutil.rmtree(newfile)
        if os.path.isfile(oldfile):
            link_file(oldfile, newfile)
        else:
            link_dir(oldfile, newfile)

def link_file(oldfile, newfile):
    """
    hard link file from oldfile to newfile
    :param oldfile:
    :param newfile:
    :return:
    """
    if not os.path.isfile(oldfile):
        raise Exception("不存在文件：%s" % oldfile)
    if os.path.exists(newfile):
        os.remove(newfile)
    os.link(oldfile, newfile)

def bin_rename(dir,type):
    """
    主要bin改名字
    :param dir:
    :return:
    """
    files=os.listdir(dir)
    if type in ['concoct']:
        for file in files:
            os.rename(dir + "/" + file,dir + "/" + "bin" + file)
    elif type in ['maxbin']:
        for file in files:
            if re.search(r'.fasta',file):
                name = str(int(file.split('.')[1]))
                os.rename(dir + "/" + file, dir + "/" + "bin" + name + '.fa')
            else:
                shutil.move(dir + "/" + file,dir + '/../' + file)
    elif type in ['metabat']:
        for file in files:
            name = file.split('.')[1]
            os.rename(dir + "/" + file, dir + "/" + "bin" + name + '.fa')
    elif type in ['all']:
        n =1
        for file in files:
            os.rename(dir + "/" + file, dir + "/" + "bin" + str(n) + '.fa')
            n +=1

def export_longest_seq(oldfile, newfile, type="fasta"):
    max_len = 1
    max_record = 0
    for record in SeqIO.parse(oldfile, type):
        if len(record) > max_len:
            max_len = len(record)
            max_record = record
    SeqIO.write(max_record, newfile, type)

def download_from_ncbi(query_id, outfile, db="nucleotide"):
    handle = Entrez.efetch(db=db, id=query_id,rettype="fasta", retmode="text")
    text = handle.read()[:-1]
    with open(outfile, "w") as file:
        file.write(text)

def get_amphora_tax(task_id, bin_name):
    db = Config().get_mongo_client(mtype="metagbin", ref=False)[Config().get_mongo_dbname("metagbin", ref=False)]
    collection = db.assembly
    coll_detail = db.assembly_detail
    rsult = collection.find_one({"task_id": task_id, "genome_id": bin_name})
    if rsult:
        result = coll_detail.find_one({"assemble_id": rsult["main_id"]})
        return result["taxon"]
    else:
        return "unknown"

def get_ani_species(ani_ref_name):
    db_ref = Config().get_mongo_client(mtype="metagbin", ref=True)[Config().get_mongo_dbname("metagbin", ref=True)]
    collection = db_ref.bac_genome
    database_name = ani_ref_name.replace(".fna", "")
    result = collection.find_one({"database_name": database_name})
    if result:
        accession = result["accession"]
        isolation = accession.split("[")[-1].rstrip("]")
        return isolation
    else:
        return "unknown"

def get_pocp_genus(task_id, pocp_ref_name):
    db_ref = Config().get_mongo_client(mtype="metagbin", ref=True)[Config().get_mongo_dbname("metagbin", ref=True)]
    genome_coll = db_ref.bac_genome
    database_name = pocp_ref_name.replace(".faa", "").replace(".fna", "")
    result = genome_coll.find_one({"database_name": database_name})
    if result:
        tax_id = result["tax_id"]
    else:
        return get_amphora_tax(task_id, database_name)
    mg_taxons = mg_taxon()
    taxons = tx()
    taxons.taxid = int(tax_id)
    tracks = taxons.get_track()
    if not tracks:
        print "不存在的tax_id:%s" % tax_id
        return "unknown"
    taxons.track.reverse()
    track = ['{}{{{}}}'.format(i.spname, i.rank) for i in taxons.track]
    if not track:
        return "unknown"
    taxon_name = {'d': '', 'k': '', 'p': '', 'c': '', 'o': '', 'f': '', 'g': '', 's': ''}
    species_name = []
    for taxon in track:
        m = re.match(r"(.+){(.+)}$", taxon)
        m1 = m.group(1)
        m2 = m.group(2)
        taxon_name = mg_taxons.level_name(m1, m2, taxon_name)
    species_name = mg_taxons.add_unclassifed(taxon_name, species_name)
    print species_name
    return ";".join(species_name)
    # for names in species_name:
    #     if names.startswith("g__"):
    #         return names
    # return "unknown"

class Fasta(object):
    def __init__(self, path):
        super(Fasta, self).__init__()
        self.path = path
        if not os.path.isfile(path):
            raise Exception("没有找到文件%s" % path)
        self.base_sum = 0
        self.length = []
        self.no_c = 0
        self.no_g = 0
        self.no_a = 0
        self.no_t = 0
        self.no_n = 0
        self.n50 = 0
        self.longest = 0
        self.longest_len = 0

    def parse(self):
        for record in SeqIO.parse(open(self.path), "fasta"):
            seq = record.seq.lower()
            self.base_sum += len(seq)
            if len(seq) > self.longest_len:
                self.longest_len = len(seq)
                self.longest = record
            self.length.append(len(seq))
            self.no_c += seq.count("c")
            self.no_g += seq.count("g")
            self.no_a += seq.count("a")
            self.no_t += seq.count("t")
            self.no_n += seq.count("n")
        self.length.sort()
        self.length.reverse()

    def calculate_nxx(self, value=50): # example: value=50 if calculate N50
        pos = self.base_sum * float(value) / 100.0
        seq_sum = 0
        nxx = 0
        for seq_len in self.length:
            seq_sum += seq_len
            if pos <= seq_sum:
               nxx = seq_len
               break
        return nxx

    @property
    def seq_num(self):
        return len(self.length)

    @property
    def min_len(self):
        return min(self.length)

    @property
    def max_len(self):
        return max(self.length)
    @property
    def c_percent(self):
        c_percent = float(self.no_c) * 100 / self.base_sum
        return c_percent

    @property
    def g_percent(self):
        g_percent = float(self.no_g) * 100 / self.base_sum
        return g_percent

    @property
    def a_percent(self):
        a_percent = float(self.no_a) * 100 / self.base_sum
        return a_percent

    @property
    def t_percent(self):
        t_percent = float(self.no_t) * 100 / self.base_sum
        return t_percent

    @property
    def n_percent(self):
        n_percent = float(self.no_n) * 100 / self.base_sum
        return n_percent

    def export_longest_seq(self, outfile, type="fasta"):
        SeqIO.write(self.longest, outfile, type)

    def check_n50(self, lowest):
        n50 = self.calculate_nxx(50)
        if n50 >= lowest:
            return True