# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
'''
该脚本用于下载整理ENSEMBLE来源的数据库
'''
import re
import os
import sys
import ftplib
import subprocess
import argparse
import random
import requests
import datetime
import json

class Ensemble(object):

    def __init__(self):
        self.host = "ftp://ftp.ensemblgenomes.org"
        self.fa = ""
        self.release_version = "current"
        self.db_class = "plant"
        self.species_name = "unknown"
        self.version = ""
        self.genome_version = ""
        self.biomart_type = "type3"
        self.biomart_abr = ""
        self.rel_dir = "./"
        self.download_dir = "./"
        self.exon_fa = ""
        self.accession = ""
        self.taxon = 0
        self.json = dict()
        self.cmd_record = list()

    def run_cmd(self, cmd):
        print cmd + "\n"
        self.cmd_record.append(cmd)
        c = os.system(cmd)
        return c


    def ftp_download(self, file_ftp, file_out):
        '''
        根据ftp地址下载数据
        '''
        print("{}/{} download start".format(self.host, file_ftp))

        try:
            f=ftplib.FTP(self.host)
        except ftplib.error_perm:
            print('Can not contect"{}" '.format(self.host))
            return
        print('Connect "{}" successfully!'.format(self.host))
        try:
            f.login('anonymous','1')
        except ftplib.error_perm:
            print("Fail to login in !")
            f.quit()
            return
        try:
            f.retrbinary('RETR %s' % file_ftp, open(file_out,'wb').write)
            print("{} download sucessfully".format(file_ftp))
        except:
            print('{} 文件无法打开'.format(file_ftp))
            return

    def ftp_download2(self, file_ftp, file_out):
        '''
        根据ftp地址下载数据
        '''
        file_ftp = "ftp://" + self.host + "/" + file_ftp
        self.trys = 0
        if os.path.exists(file_out + ".finished"):
            print("{} download finished jump this step".format(file_ftp))
            return

        print("{} download start".format(file_ftp))

        base_name = os.path.basename(file_ftp)
        if os.path.exists(base_name):
            os.remove(base_name)
        os.system("wget {}".format(file_ftp))
        cmd = self.run_cmd("gunzip -f {}".format(base_name))
        if cmd != 0:
            print "error code is {}".format(cmd)
            self.ftp_re_download2(file_ftp, file_out)

        os.link(base_name[:-3], file_out[:-3])
        if os.path.exists(base_name[:-3]):
            os.remove(base_name[:-3])
        self.run_cmd("touch {}.finished".format(file_out))

    def ftp_re_download2(self, file_ftp, file_out):
        '''
        根据ftp地址下载数据
        '''
        print("{} 结果错误重新下载".format(file_ftp))
        self.trys += 1
        if self.trys > 25:
            raise("下载次数过多多，请检查路径是否正确" )

        self.run_cmd("wget -c {}".format(file_out, file_ftp))
        base_name = os.path.basename(file_ftp)
        cmd = self.run_cmd("gunzip -f {}".format(base_name))

        if cmd != 0:
            print "error code is {}".format(cmd)
            self.ftp_re_download2(file_ftp, file_out)


    def get_ftp_file(self):
        '''
        获取dna/gtf/cds/pep 文件在ftp中的相关地址，并下载
        '''
        if self.db_class == "vertebrates":
            db_class_dir = ""
        else:
            db_class_dir = self.db_class
        self.ftp_fa = "/".join(["pub",
                                db_class_dir,
                                self.release_version,
                                "fasta",
                                self.species_name.lower(),
                                "dna",
                                self.species_name + "." + self.genome_version + ".dna.toplevel.fa.gz"
        ])

        self.ftp_download2(self.ftp_fa, self.download_dir + self.species_name + ".dna.toplevel.fa.gz")
                # os.system("gunzip -f {}".format(self.download_dir + self.species_name + ".dna.toplevel.fa.gz"))
                # os.link(self.download_dir + self.species_name + ".dna.toplevel.fa",
                #        self.rel_dir + 'dna/' + self.species_name + ".dna.toplevel.fa")
        self.json.update({
            "dna_fa": self.rel_dir + 'dna/' + self.species_name + ".dna.toplevel.fa"
        })

        release_num = "." +  self.release_version.split("-")[-1]
        if self.genome_version.endswith(release_num):
            release_num = ""

        self.ftp_gtf = "/".join(["pub",
                               self.release_version,
                                 db_class_dir,
                               "gtf",
                               str(self.species_name).lower(),
                               self.species_name + "." + self.genome_version + release_num +  ".gtf.gz"
                               ])
        self.ftp_download2(self.ftp_gtf,
                           self.download_dir + self.species_name + "." + self.genome_version  + release_num + ".gtf.gz")

        self.json.update({
            "gtf": self.rel_dir + 'gtf/' + self.species_name + "." + self.genome_version + release_num +  ".gtf"
        })

        self.ftp_pep = "/".join(["pub",
                                 db_class_dir,
                               self.release_version,
                               "fasta",
                               self.species_name.lower(),
                               "pep",
                               self.species_name + "." + self.genome_version + ".pep.all.fa.gz"
                               ])

        self.ftp_download2(self.ftp_pep, self.download_dir + self.species_name + "." + self.genome_version + ".pep.all.fa.gz")

        self.json.update({
            "pep": self.rel_dir + 'cds/' + self.species_name + "." + self.genome_version + ".pep.all.fa"
        })

        self.ftp_cds = "/".join(["pub",
                                 db_class_dir,
                               self.release_version,
                               "fasta",
                               self.species_name.lower(),
                               "cds",
                               self.species_name + "." + self.genome_version + ".cds.all.fa.gz"
                               ])
        self.ftp_download2(self.ftp_cds, self.download_dir + self.species_name + "." + self.genome_version + ".cds.all.fa.gz")
        #os.link(self.download_dir + self.species_name + "." + self.genome_version + ".cds.all.fa",
        #            self.rel_dir + 'cds/' + self.species_name + "." + self.genome_version + ".cds.all.fa")
        self.json.update({
            "cds": self.rel_dir + 'cds/' + self.species_name + "." + self.genome_version + ".cds.all.fa"
        })

    def get_biomart_go_file(self):
        '''
        获取biomart GO注释文件
        '''
        path = sys.path[0]
        abr = self.biomart_abr + '_eg_gene'
        mart = self.db_class + '_mart'
        bio_class = self.db_class

        if self.db_class == "vertebrates" or self.db_class == "":
            mart = "default"
            abr = self.biomart_abr +  "_gene_ensembl"
            bio_class = "asia"
        get1 = "sh {} {} {} {}".format(os.path.join(path, "get_go_info.sh"), abr, mart, bio_class)
        outfile = abr + "_go.txt"

        if os.path.exists(self.rel_dir + 'get_go' + '.finished'):
            print "go 下载成功， 跳过此步"
        else:
            while 1:
                f = self.run_cmd(get1)
                if f != 130:
                    break
            with open(outfile, 'r') as out:
                lines = out.readlines()
                if lines[-1].strip() == "[success]":
                    pass
                else:
                    print lines
                    print "go 注释没有下载完整"
            self.run_cmd("mv {} {}".format(outfile, self.download_dir))
            self.run_cmd("sed '$d' -i {}".format(
                self.download_dir + outfile
            ))
        self.json.update({
            "go": self.rel_dir + 'GO/' + outfile
        })

    def get_biomart_enterz_file(self):
        '''
        获取biomart enterz注释文件
        '''
        path = sys.path[0]
        abr = self.biomart_abr + "_eg_gene"
        mart = self.db_class + '_mart'
        bio_class = self.db_class

        if self.db_class == "vertebrates" or self.db_class == "":
            mart = "default"
            abr = self.biomart_abr +  "_gene_ensembl"
            bio_class = "asia"
        get1 = "sh {} {} {} {}".format(os.path.join(path, "get_enterz_info.sh"), abr, mart, bio_class)
        outfile = abr + "_entrez.txt"
        if os.path.exists(self.rel_dir + 'get_enterz' + '.finished'):
            print "enterz 下载成功， 跳过此步"
        else:
            while 1:
                f = self.run_cmd(get1)
                if f != 130:
                    break
            with open(outfile, 'r') as out:
                lines = out.readlines()
                if lines[-1].strip() == "[success]":
                    pass
                else:
                    print lines
            self.run_cmd("mv {} {}".format(outfile, self.download_dir))
            self.run_cmd("sed '$d' -i {}".format(
                self.download_dir + outfile,
            ))
        self.json.update({
            "ensemble2entrez": self.rel_dir + 'NCBI/' + outfile
        })
    def get_acc_taxon(self):
        '''
        获取taxon ID accession 使用ensemble rest api
        '''
        species_name = self.species_name
        # .replace("_", " ")
        # server = "http://rest.ensemblgenomes.org"
        ext = "/info/assembly/{}?".format(species_name)
        # if self.db_class == "vertebrates":
        server = "http://rest.ensembl.org"
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

        if not r.ok:
            ext = "/info/assembly/{}?".format(species_name.replace("_", " "))
            r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
            if not r.ok:
                r.raise_for_status()
        else:
            decoded = r.json()
            with open(self.download_dir + 'assembly.json', 'wb') as f:
                json.dump(r.json(), f, indent=2)
            self.accession = decoded['assembly_accession']

        if len(species_name.split("_")) >2:
            species_name = "_".join(species_name.split("_")[:2])
        ext = "/taxonomy/classification/{}?".format(species_name)
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

        if not r.ok:
            ext = "/taxonomy/classification/{}?".format(species_name.replace("_", " "))
            r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
            if not r.ok:
                r.raise_for_status()
        else:
            self.taxon = int(r.json()[0]['id'])
            decoded = r.json()
            with open(self.download_dir + 'taxonomy.json', 'wb') as f:
                json.dump(r.json(), f, indent=2)

    def get_biomart_file(self):
        '''
        获取biomart文件
        '''
        path = sys.path[0]
        abr = self.biomart_abr + "_eg_gene"
        mart = self.db_class + '_mart'
        bio_class = self.db_class

        if self.db_class == "vertebrates" or self.db_class == "":
            mart = "default"
            abr = self.biomart_abr +  "_gene_ensembl"
            bio_class = "asia"
        get1 = "sh {} {} {} {}".format(os.path.join(path, "get_gene_info1.sh"), abr, mart, bio_class)
        get2 = "sh {} {} {} {}".format(os.path.join(path, "get_gene_info2.sh"), abr, mart, bio_class)
        get3 = "sh {} {} {} {}".format(os.path.join(path, "get_gene_info3.sh"), abr, mart, bio_class)
        outfile = abr + "_gene.txt"

        if os.path.exists(self.rel_dir + 'get_biomart' + '.finished'):
            with open(self.rel_dir + 'get_biomart' + '.finished', 'rb') as f:
                self.biomart_type = f.read()
            print "bio 下载成功， 跳过此步"
        else:
            while 1:
                f = self.run_cmd(get1)
                if f != 130:
                    break

            print get1
            # 按照type1 type2 type3的顺序依次获得biomart的注释
            with open(outfile, 'r') as out:
                lines = out.readlines()
                if lines[-1].strip() == "[success]":
                    self.biomart_type = "type1"
                    print "biomart download sucess"
                elif "Attribute" in lines[-1]:
                    print lines
                    print "try type2"
                    while 1:
                        f = self.run_cmd(get2)
                        if f != 130:
                            break
                    with open(outfile, 'r') as out2:
                        lines2 = out2.readlines()
                        if lines2[-1].strip() == "[success]":
                            self.biomart_type = "type2"
                            print "biomart download sucess"
                        elif "Attribute" in lines[-1]:
                            print lines2
                            print "try type3"
                            while 1:
                                f = self.run_cmd(get3)
                                if f != 130:
                                    break
                            with open(outfile, 'r') as out3:
                                lines3 = out3.readlines()
                                if lines3[-1].strip() == "[success]":
                                    self.biomart_type = "type3"
                                    print "biomart download sucess"
                                else:
                                    print "can't  get biomart file"
                        else:
                            print lines
                            pass
                else:
                    print lines
                    pass
            self.run_cmd("mv {} {}".format(outfile, self.download_dir))
            self.run_cmd("sed '$d' -i {}".format(
                self.download_dir + outfile,
            ))

            type2line = {"type1":"1,2,7", "type2":"1,2,5", "type3":"1,2,3"}
            self.run_cmd("cut -f {} {} > {}".format(
                type2line[self.biomart_type],
                self.download_dir + outfile,
                self.download_dir + 'g2t2p'
            ))

        self.json.update({
            "bio_mart_annot": self.rel_dir + 'biomart/' + outfile,
            "biomart_abr": abr,
            "biomart_gene_annotype": self.biomart_type,
            "g2t2p": self.rel_dir + 'Annotation_v2/g2t2p'
        })

    '''
    def get_index(self):
        #
        # 调用tools建索引，统计数据， 提取转录本序列
        #
        data = {
            "id": "genome_index" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v2.genome_index",
            "instant": False,
            "options": dict(
                fasta=os.getcwd() + '/' + self.json["dna_fa"],
                gtf=os.getcwd() + '/' + self.json["gtf"],
                index=os.getcwd() + '/' + ".".join(self.json["dna_fa"].split(".")[:-1]) + "_index",
                transcript=os.getcwd() + '/' + self.rel_dir + "dna/transcript.fa",
                stat=os.getcwd() + '/' + self.rel_dir + "gtf/genome_stat.xls"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        if os.path.exists(self.rel_dir + 'get_index' + '.finished'):
            print "index 运行成功， 跳过此步"
        else:
            wf.run()
        with open(self.json["dna_fa"] + ".fai") as f:
            length = sum([int(x.split("\t")[1]) for x in f.readlines()])/1000000
        self.json.update({
            "dna_index": ".".join(self.json["dna_fa"].split(".")[:-1]) + "_index",
            "gene_stat": self.rel_dir + "gtf/genome_stat.xls",
            "transcript": self.rel_dir + "dna/transcript.fa",
            "size": str(length),
        })
    '''


    def db2taxon(self, db_class):
        db2taxon = {
            "protists": "Protist",
            "vertebrates": "Animal",
            "metazoa": "Animal",
            "fungi": "Fungi",
            "plants": "Plant",
        }
        return db2taxon[db_class]

    def get_json(self):
        self.json.update({
            "taxon": self.db2taxon(self.db_class) ,
            "ensemble_class": self.db_class,
            "ensemble_web": "http://{}.ensembl.org/{}/Info/Index".format(self.db_class.replace("vertebrates", "www"), self.species_name),
            "ensemble_release": "Ensemble_" + self.release_version.replace("-", "_"),
            "name": self.species_name,
            "web_dir": self.fa_ftp,
            "common_name": "",
            "classification": "",
            "taxon_id": self.taxon,
            "assembly": self.genome_version,
            "accession": self.accession,
            "genebuild_method": "",
            "kegg_genome_abr": self.db2taxon(self.db_class).lower(),
            "kegg_genome_name": "",
            "ncbi_ensemble_tax": "",
            "kegg_use": "",
            "index": self.species_name + "." + self.genome_version,
            "kegg": "",
            "cog": "",
            "anno_path_v2": self.rel_dir + "/Annotation_v2",
        })
        with open(self.rel_dir + self.species_name + '.json', 'wb') as f:
            json.dump(self.json, f, indent=2)
        return self.json


    '''
    def get_orf(self):
        test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_ref_rna_v2'
        data = {
            "id": "annot_orf" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "module",
            "name": "ref_rna_v2.annot_orfpfam",
            "instant": False,
            "options": dict(
                pep =  os.getcwd() + "/" +  self.json['pep'],
                lines = 5000,
                gtf = os.getcwd() + "/" + self.json['gtf'],
                g2t2p = os.getcwd() + "/" + self.json['g2t2p']
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        if os.path.exists(self.rel_dir + 'get_orf' + '.finished'):
            print "orf 运行成功， 跳过此步"
        else:
            wf.run()
            CopyFile().linkdir(wf.output_dir, self.rel_dir + "Annotation_v2/annot_orfpfam")
    '''

    '''
    def get_annotation(self):
        db2nr = {
            "protists": "protist",
            "vertebrates": "metazoa",
            "metazoa": "metazoa",
            "fungi": "fungi",
            "plants": "viridiplantae",
        }
        data = {
            "id": "annot_db_" + self.species_name + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "module",
            "name": "ref_rna_v2.annot_mapdb",
            "instant": False,
            "options": dict(
                query = os.getcwd() + "/" + self.json['transcript'],
                method = "diamond",
                nr_db = db2nr[self.db_class],
                known_go = os.getcwd() + "/" + self.json['go'],
                lines = 20000
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        if os.path.exists(self.rel_dir + 'get_annotation' + '.finished'):
            print "annotation 运行成功， 跳过此步"
        else:
            wf.run()
            CopyFile().linkdir(wf.output_dir, self.rel_dir + "Annotation_v2/annot_mapdb")
    '''

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-release_version', type=str, required=False,
                        help = 'release-39 like ensemble version', default = "")
    parser.add_argument('-db_class', type=str, required=False, default = "",
                        help = 'protists,vertebrates,fungi,plants,metazoa\n default mean "" mean vertebrates')
    parser.add_argument('-species_name', type=str, required=False, default = "",
                        help = 'Callithrix_jacchus')
    parser.add_argument('-genome_version', type=str, required=False, default = "",
                        help = 'C_jacchus3.2.1')
    parser.add_argument('-host', type=str, required=False, default = "ftp.ensembl.org",
                        help = 'ftp.ensembl.org')
    parser.add_argument('-fa_ftp', type=str, required=False, default = "",
                        help = 'ftp://ftp.ensembl.org/pub/release-99/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa.gz\n' +
                        'ftp://ftp.ensemblgenomes.org/pub/release-46/fungi/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz\n' +
                        'be carefull should not add / in dir')

    args = parser.parse_args()
    ensemble = Ensemble()

    if args.fa_ftp:
        if args.fa_ftp.startswith("ftp://ftp") and args.fa_ftp.endswith("toplevel.fa.gz"):
            pass
        else:
            raise("fa_ftp should start with ftp://ftp and end with toplevel.fa.gz")

        dirs = args.fa_ftp.split("/")
        ensemble.release_version = dirs[4]
        if len(dirs) == 9:
            ensemble.db_class = "vertebrates"
        elif len(dirs) == 10:
            if dirs[4].startswith("release"):
                ensemble.db_class = dirs[5]
                ensemble.release_version = dirs[4]
            else:
                ensemble.db_class = dirs[4]
                ensemble.release_version = dirs[5]
            # ensemble.release_version = dirs[5]

        ensemble.host = dirs[2]


        names = dirs[-1].split('.dna')[0].split(".")
        ensemble.species_name = names[0]
        ensemble.genome_version = ".".join(names[1:])
        ensemble.fa_ftp = args.fa_ftp

    else:
        ensemble.release_version = args.release_version
        ensemble.genome_version = args.genome_version
        ensemble.db_class = args.db_class
        db_class_dir  = ensemble.db_class
        if ensemble.db_class == "vertebrates":
            ensemble.db_class = "vertebrates"
            db_class_dir
        ensemble.species_name = args.species_name
        ensemble.host = args.host

        ensemble.fa_ftp = "ftp://{}/pub/release-89/{}/fasta/{}/dna/{}.{}.dna.toplevel.fa.gz".format(args.host, db_class_dir,  args.species_name.lower(), args.species_name, args.genome_version)

    spes = ensemble.species_name.split("_")
    a = spes[0][0]
    b = "".join(spes[1:])
    ensemble.biomart_abr = (a+b).lower()

    rel_dir = ensemble.db_class + '/' + ensemble.species_name + '/Ensemble_' + ensemble.release_version.replace("-", "_") + "/"
    # for genome_dirs in ['Annotation', 'Annotation_v2', 'biomart', 'cds', 'COG', 'dna', 'GO', 'gtf', 'KEGG', 'NCBI', 'download']:
    for genome_dirs in ['download']:
        if os.path.exists(rel_dir + genome_dirs):
            pass
        else:
            os.makedirs(rel_dir + genome_dirs)
    ensemble.rel_dir  = rel_dir
    ensemble.download_dir = rel_dir + 'download/'
    ensemble.get_acc_taxon()

    if os.path.exists("{}.finished".format(ensemble.rel_dir + 'get_ftp')):
        pass
    else:
        ensemble.get_ftp_file()
        os.system("touch {}.finished".format(ensemble.rel_dir + 'get_ftp'))
    if os.path.exists("{}.finished".format(ensemble.rel_dir + 'get_biomart')):
        pass
    else:
        ensemble.get_biomart_file()
        os.system("touch {}.finished".format(ensemble.rel_dir + 'get_biomart'))

    with open(ensemble.rel_dir + 'get_biomart.finished', 'wb') as f:
        f.write(ensemble.biomart_type)
    if os.path.exists("{}.finished".format(ensemble.rel_dir + 'get_go')):
        pass
    else:
        ensemble.get_biomart_go_file()
        os.system("touch {}.finished".format(ensemble.rel_dir + 'get_go'))
    if os.path.exists("{}.finished".format(ensemble.rel_dir + 'get_enterz')):
        pass
    else:
        ensemble.get_biomart_enterz_file()
        os.system("touch {}.finished".format(ensemble.rel_dir + 'get_enterz'))

    with open(ensemble.rel_dir + 'down_load.cmds.sh', 'aw') as f:
        for cmd in ensemble.cmd_record:
            f.write("{}\n".format(cmd))

    '''
    ensemble.get_index()
    os.system("touch {}.finished".format(ensemble.rel_dir + 'get_index'))
    ensemble.get_annotation()
    os.system("touch {}.finished".format(ensemble.rel_dir + 'get_annotation'))
    ensemble.get_orf()
    os.system("touch {}.finished".format(ensemble.rel_dir + 'get_orf'))
    '''
    ensemble.get_json()
