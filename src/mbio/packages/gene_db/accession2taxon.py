# -*- coding: utf-8 -*-
# __author__ = 'bingxu.liu'
# last_modified:20190808
# 改脚本参考ete3.ncbi_taxonomy.ncbiquery

import subprocess
from biocluster.config import Config
import sqlite3
import sys
import os
import tarfile
import math
import pickle
import argparse
import datetime
from ete3.coretype.tree import Tree

class taxon(object):
    """ncbi物种分类中一个节点分类的对象，提供一些查询方法和信息存储"""
    def __init__(self):
        self.acc2tax = Config().SOFTWARE_DIR + '/database/Annotation/latest_sqlite2019/accession2taxid.db'
        self.dbfile = Config().SOFTWARE_DIR + '/database/Annotation/latest_sqlite2019/taxonomy.db'
        self._connect()

    def _connect(self):
        self.db = sqlite3.connect(self.dbfile)
        self.acc2tax_db = sqlite3.connect(self.acc2tax)

    def get_accession_ranks(self, accession_ids, desired_ranks):
        '''
        根据物种accession_id列表，返回taxon分类列表
        '''
        ranks = []
        for accession_id in accession_ids:
            taxid = self.get_taxid_from_accid(accession_id)
            if taxid:
                ranks.append(self.get_desired_ranks(taxid, desired_ranks))
            else:
                ranks.append("unknown")
        return ranks

    def get_name_from_taxon(self, taxon_ids):
        '''
        根据taxid列表，返回species name列表
        '''
        species_names = []
        for taxon_id in taxon_ids:
            result = self.db.execute('SELECT spname FROM species WHERE taxid="%s"' %taxon_id)
            species_name  = result.fetchone()
            if species_name:
                species_names.append(species_name[0])
            else:
                species_names.append("unknown")
        return species_names

    def get_taxon_from_name(self, name):
        '''
        根据 species name 返回taxon_id
        '''
        result = self.db.execute('SELECT taxid FROM species WHERE spname="%s"' %name)
        taxon_id  = result.fetchone()
        if taxon_id:
            return str(taxon_id[0])
        else:
            return ""

    def get_track_from_name(self, name):
        '''
        根据 species name 返回taxon_id
        '''
        result = self.db.execute('SELECT taxid FROM species WHERE spname="%s"' %name)
        track  = result.fetchone()
        if track:
            return str(track[0])
        else:
            return None

    def get_taxid_from_accid(self, accession_id):
        '''
        根据NCBI accession_id获得taxon_id
        '''
        result = self.acc2tax_db.execute('SELECT taxid FROM prot_accession2taxid WHERE accv="%s"' %accession_id)
        taxid = result.fetchone()
        if not taxid:
            print "{} taxid not found in prot_accession2taxid".format(taxid)
            return 0
        else:
            taxid = taxid[0]
            return taxid

    def get_chilren_from_list(self, taxon_id_list):
        """
        根据列表获取所有子物种分类的taxon id
        """
        all_children = []
        for taxon_id in taxon_id_list:
            all_children.append(taxon_id)
            self.get_chilren_from_taxon(taxon_id, all_children)
        return set(all_children)

    def get_chilren_from_taxon(self, taxon_id, all_children):
        """
        根据taxon id获取所有子物种分类的taxon_id
        """
        cmd = "select taxid, parent FROM species WHERE parent IN (%s);" %taxon_id
        result = self.db.execute(cmd)
        for tax, spname in result.fetchall():
            all_children.append(tax)
            self.get_chilren_from_taxon(tax, all_children)


    def get_desired_ranks(self, taxid, desired_ranks=None):
        lineage = self.get_lineage(taxid)
        if lineage:
            lineage2ranks = self.get_rank(lineage)
            ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())

            desired_taxon = [ranks2lineage.get(rank, 0) for rank in desired_ranks]
            names = self.get_name_from_taxon(desired_taxon)
            return ";".join([names[n] + "{" + desired_ranks[n] + "}" for n,i in enumerate(desired_ranks)])
        else:
            return ";".join([ "unknown" + "{" + desired_ranks[n] + "}" for n,i in enumerate(desired_ranks)])


    def get_linkage_name(self, taxid):
        lineage = self.get_lineage(taxid)
        names = self.get_name_from_taxon(lineage)
        return ";".join("{}({})".format(name, taxon) for taxon,name in zip(lineage, names))



    def get_rank(self, taxids):
        'return a dictionary converting a list of taxids into their corresponding NCBI taxonomy rank'
        all_ids = set(taxids)
        all_ids.discard(None)
        all_ids.discard("")
        query = ','.join(['"%s"' %v for v in all_ids])
        cmd = "select taxid, rank FROM species WHERE taxid IN (%s);" %query
        result = self.db.execute(cmd)
        id2rank = {}
        for tax, spname in result.fetchall():
            id2rank[tax] = spname
        return id2rank

    def get_lineage(self, taxid):
        """Given a valid taxid number, return its corresponding lineage track as a
        hierarchically sorted list of parent taxids.
        """
        if not taxid:
            return None
        result = self.db.execute('SELECT track FROM species WHERE taxid=%s' %taxid)
        raw_track = result.fetchone()
        if not raw_track:
            #perhaps is an obsolete taxid
            _, merged_conversion = self._translate_merged([taxid])
            if taxid in merged_conversion:
                result = self.db.execute('SELECT track FROM species WHERE taxid=%s' %merged_conversion[taxid])
                raw_track = result.fetchone()
            # if not raise error
            if not raw_track:
                #raw_track = ["1"]
                print ValueError("%s taxid not found" %taxid)
            else:
                print "taxid {} was translated into {}".format(taxid, merged_conversion[taxid])
        if not raw_track:
            return 0
        else:
            track = list(map(int, raw_track[0].split(",")))
            return list(reversed(track))

    def update_taxonomy_database(self, taxdump_file=None):
        """Updates the ncbi taxonomy database by downloading and parsing the latest
        taxdump.tar.gz file from the NCBI FTP site (via HTTP).

        :param None taxdump_file: an alternative location of the taxdump.tax.gz file.
        """
        if not taxdump_file:
            self.update_db(self.dbfile)
        else:
            self.update_db(self.dbfile, taxdump_file)

    def update_db(self, dbfile, targz_file=None):
        basepath = os.path.split(dbfile)[0]
        if basepath and not os.path.exists(basepath):
            os.mkdir(basepath)

        if not targz_file:
            try:
                from urllib import urlretrieve
            except ImportError:
                from urllib.request import urlretrieve

            print('Downloading taxdump.tar.gz from NCBI FTP site (via HTTP)...')
            urlretrieve("http://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz", "taxdump.tar.gz")
            print('Done. Parsing...')
            targz_file = "taxdump.tar.gz"

        tar = tarfile.open(targz_file, 'r')
        t, synonyms = self.load_ncbi_tree_from_dump(tar)
        prepostorder = [int(node.name) for post, node in t.iter_prepostorder()]
        pickle.dump(prepostorder, open(dbfile+'.traverse.pkl', "wb"), 2)

        print("Updating database: %s ..." %dbfile)
        self.generate_table(t)
        open("syn.tab", "w").write('\n'.join(["%s\t%s" %(v[0],v[1]) for v in synonyms]))
        with open("merged.tab", "w") as merged:
            for line in tar.extractfile("merged.dmp"):
                line = str(line.decode())
                out_line = '\t'.join([_f.strip() for _f in line.split('|')[:2]])
                merged.write(out_line+'\n')
        try:
            self.upload_data(dbfile)
        except:
            raise
        else:
            os.system("rm syn.tab merged.tab taxa.tab")
            # remove only downloaded taxdump file
            if not targz_file:
                os.system("rm taxdump.tar.gz")

    def generate_table(self, t):
        OUT = open("taxa.tab", "w")
        for j, n in enumerate(t.traverse()):
            if j%1000 == 0:
                print("\r", j, "generating entries...")
            temp_node = n
            track = []
            while temp_node:
                track.append(temp_node.name)
                temp_node = temp_node.up
            if n.up:
                OUT.write('\t'.join([n.name, n.up.name, n.taxname, getattr(n, "common_name", ""), n.rank, ','.join(track)]) + "\n")
            else:
                OUT.write('\t'.join([n.name, "", n.taxname, getattr(n, "common_name", ""), n.rank, ','.join(track)]) + "\n")
        OUT.close()


    def load_ncbi_tree_from_dump(self, tar):
        # Download: http://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
        parent2child = {}
        name2node = {}
        node2taxname = {}
        synonyms = set()
        name2rank = {}
        node2common = {}
        print("Loading node names...")
        for line in tar.extractfile("names.dmp"):
            line = str(line.decode())
            fields =  [_f.strip() for _f in line.split("|")]
            nodename = fields[0]
            name_type = fields[3].lower()
            taxname = fields[1]
            if name_type == "scientific name":
                node2taxname[nodename] = taxname
            if name_type == "genbank common name":
                node2common[nodename] = taxname
            elif name_type in set(["synonym", "equivalent name", "genbank equivalent name",
                                   "anamorph", "genbank synonym", "genbank anamorph", "teleomorph"]):
                synonyms.add( (nodename, taxname) )
        print(len(node2taxname), "names loaded.")
        print(len(synonyms), "synonyms loaded.")

        print("Loading nodes...")
        for line in tar.extractfile("nodes.dmp"):
            line = str(line.decode())
            fields =  line.split("|")
            nodename = fields[0].strip()
            parentname = fields[1].strip()
            n = Tree()
            n.name = nodename
            n.taxname = node2taxname[nodename]
            if nodename in node2common:
                n.common_name = node2common[nodename]
            n.rank = fields[2].strip()
            parent2child[nodename] = parentname
            name2node[nodename] = n
        print(len(name2node), "nodes loaded.")

        print("Linking nodes...")
        for node in name2node:
           if node == "1":
               t = name2node[node]
           else:
               parent = parent2child[node]
               parent_node = name2node[parent]
               parent_node.add_child(name2node[node])
        print("Tree is loaded.")
        return t, synonyms

    def upload_data(self, dbfile):
        print('Uploading to', dbfile)
        basepath = os.path.split(dbfile)[0]
        if basepath and not os.path.exists(basepath):
            os.mkdir(basepath)

        create_cmd = """
        DROP TABLE IF EXISTS stats;
        DROP TABLE IF EXISTS species;
        DROP TABLE IF EXISTS synonym;
        DROP TABLE IF EXISTS merged;
        CREATE TABLE stats (version INT PRIMARY KEY);
        CREATE TABLE species (taxid INT PRIMARY KEY, parent INT, spname VARCHAR(50) COLLATE NOCASE, common VARCHAR(50) COLLATE NOCASE, rank VARCHAR(50), track TEXT);
        CREATE TABLE synonym (taxid INT,spname VARCHAR(50) COLLATE NOCASE, PRIMARY KEY (spname, taxid));
        CREATE TABLE merged (taxid_old INT, taxid_new INT);
        CREATE INDEX spname1 ON species (spname COLLATE NOCASE);
        CREATE INDEX spname2 ON synonym (spname COLLATE NOCASE);
        """
        for cmd in create_cmd.split(';'):
            self.db.execute(cmd)
        print()

        self.db.execute("INSERT INTO stats (version) VALUES (%d);" %int(datetime.datetime.now().strftime('%y%m%d')))
        self.db.commit()

        for i, line in enumerate(open("syn.tab")):
            if i%5000 == 0 :
                print('\rInserting synonyms:     % 6d' %i)
                sys.stderr.flush()
            taxid, spname = line.strip('\n').split('\t')
            self.db.execute("INSERT INTO synonym (taxid, spname) VALUES (?, ?);", (taxid, spname))
        print()
        self.db.commit()
        for i, line in enumerate(open("merged.tab")):
            if i%5000 == 0 :
                print('\rInserting taxid merges: % 6d' %i)
                sys.stderr.flush()
            taxid_old, taxid_new = line.strip('\n').split('\t')
            self.db.execute("INSERT INTO merged (taxid_old, taxid_new) VALUES (?, ?);", (taxid_old, taxid_new))
        print()
        self.db.commit()
        for i, line in enumerate(open("taxa.tab")):
            if i%5000 == 0 :
                print('\rInserting taxids:      % 6d' %i)
                sys.stderr.flush()
            taxid, parentid, spname, common, rank, lineage = line.strip('\n').split('\t')
            self.db.execute("INSERT INTO species (taxid, parent, spname, common, rank, track) VALUES (?, ?, ?, ?, ?, ?);", (taxid, parentid, spname, common, rank, lineage))
        print()
        self.db.commit()

    def _translate_merged(self, all_taxids):
        conv_all_taxids = set((list(map(int, all_taxids))))
        cmd = 'select taxid_old, taxid_new FROM merged WHERE taxid_old IN (%s)' %','.join(map(str, all_taxids))

        result = self.db.execute(cmd)
        conversion = {}
        for old, new in result.fetchall():
            conv_all_taxids.discard(int(old))
            conv_all_taxids.add(int(new))
            conversion[int(old)] = int(new)
        return conv_all_taxids, conversion

if __name__ == '__main__':  # for test
    a = taxon()
    parser = argparse.ArgumentParser()
    parser.add_argument('-update', type=str, required=False, help = 'update databases by dump_file', default = "")
    parser.add_argument('-accs', type=str, required=False, help = 'find accession taxons', default = "")
    parser.add_argument('-taxs', type=str, required=False, help = 'find taxons names', default = "")
    parser.add_argument('-acc_file', type=str, required=False, help = 'find accession taxons from file', default = "")
    parser.add_argument('-track', type=str, required=False, help = 'find track by names', default = "")
    parser.add_argument('-linkage', type=str, required=False, help = 'find taxon linkage names', default = "")

    args = parser.parse_args()
    print args.update
    if args.update:
        if os.path.exists(args.update):
            a.update_taxonomy_database(args.update)
        else:
            a.update_taxonomy_database()
    if args.accs:
        accs = args.accs.split(",")
        desired_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        ranks = a.get_accession_ranks(accs, desired_ranks)
        print ranks
    if args.taxs:
        taxs = args.taxs.split(",")
        names = a.get_name_from_taxon(taxs)
        print names
    if args.linkage:
        names = a.get_linkage_name(args.linkage)
        print names
    if args.acc_file:
        with open(args.acc_file, "r") as f:
            for line in f.readlines():
                acc = line.strip()
                taxon = a.get_taxid_from_accid(acc)
                print str(acc) + "\t" + str(taxon) + "\n"

    if args.track:
        taxon_id = a.get_track_from_name(args.track)

        if taxon_id:
            desired_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
            a = a.get_desired_ranks(taxon_id, desired_ranks)
            print a
