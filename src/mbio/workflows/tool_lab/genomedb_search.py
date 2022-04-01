# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import group_detail_sort
from bson import ObjectId
import datetime
import shutil
import re,os
import time
import json
import xlsxwriter
from biocluster.workflow import Workflow
from biocluster.config import Config
import urllib2

class GenomedbSearchWorkflow(Workflow):
    """
    基因集venn图
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenomedbSearchWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "task_name", "type": "string", "default": ""},
            {"name": "file_id", "type": "string"},
            {"name": "task_type", "type": "int"},
            {"name": "submit_location", "type": "string"},
            {"name": "species", "type": "string", "default": ""},
            {"name": "ensembl", "type": "string", "default": "true"},
            {"name": "ncbi", "type": "string", "default": "true"},
            {"name": "local_only", "type": "string", "default": "true"},
            {"name": "main_id", "type": "string", "default": "5e78861fd041a47fdffffb4c"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.name_list = list()
        self.name2_list = list()
        self.search_type = "species_name"

    def get_search_type(self):
        try:
            int(self.option("species"))
            self.search_type = "taxon"
        except:
            pass


    def search_ncbi_online(self):
        results = list()
        esearch = self.config.SOFTWARE_DIR + '/bioinfo/biodb/edirect/esearch'
        efetch = self.config.SOFTWARE_DIR + '/bioinfo/biodb/edirect/efetch'
        perl_path = self.config.SOFTWARE_DIR + '/program/perl/perls/perl-5.24.0/bin'
        os.environ['PATH'] = perl_path  + ":" + os.environ['PATH']
        if self.search_type == "taxon":
            cmd = "{} -db assembly -query \"{} [TXID]\"  | {}  -format docsum  -json ".format(esearch, self.option("species"), efetch)
        else:
            cmd = "{} -db assembly -query \"{} [ORGN]\"  | {}  -format docsum  -json ".format(esearch, self.option("species"), efetch)

        self.logger.info(cmd)
        try:
            search_json = os.popen(cmd)

        except:
            self.logger.info("NCBI联网搜索失败 {}".format(search_json.readlines()))
            return results

        '''
        if return_code != 0:
            self.logger.info("NCBI联网搜索失败")
            return []
        '''
        try:
            summarys = json.load(search_json)
            with open(self.work_dir + '/ncbi_down.json', 'w') as f:
                f.write(json.dumps(summarys, indent=4))
        except:
            self.logger.info("NCBI联网搜索失败")
            return results

        for summary in summarys['DocumentSummarySet']['DocumentSummary']:
            n_dict = {
                "taxonomy_id": summary.get("Taxid", ""),
                "assembly": summary.get("AssemblyName", ""),
                "assembly_accession": summary.get("AssemblyAccession", ""),
                "assembly_level": summary.get("AssemblyStatus", ""),
                "size": summary.get("size", ""),
                "chr_num": summary.get("chr_num", ""),
                "ftp": str(summary.get("FtpPath_GenBank", "")),
                "url": "",
                "other_info": ";".join(["{}:{}".format(k,summary.get(k, ""))
                                        for k in ["Organism",
                                                  "AssemblyClass",
                                                  "SubmitterOrganization"
                                                  ]])
            }

            if 'Meta' in summary and 'Stats' in summary['Meta'] and 'Stat' in summary['Meta']['Stats']:
                for stat in summary['Meta']['Stats']['Stat']:
                    if stat['category'] == 'chromosome_count':
                        n_dict["chr_num"] = n_dict["chr_num"] + "chr:{}".format(stat['content'])
                    if stat['category'] == 'scaffold_count' and stat['sequence_tag'] == 'all':
                        n_dict["chr_num"] = n_dict["chr_num"] + "scaffold:{}".format(stat['content'])
                    if stat['category'] == 'total_length' and stat['sequence_tag'] == 'all':
                        n_dict["size"] = stat['content']
            if n_dict['ftp'].startswith("ftp"):
                results.append(n_dict)
        return results

    def search_ensembl_online(self):
        results = list()
        for name in self.name_list:
            try:
                url = "http://rest.ensembl.org/info/genomes/taxonomy/{}?content-type=application/json".format(name.replace(" ", "_"))
                self.logger.info("ensembl url is {}".format(url))
                data_json = urllib2.urlopen(url)
                e_dict = json.load(urllib2.urlopen(url))

                if e_dict.get("division") == "EnsemblVertebrates":
                    host = "ftp://ftp.ensembl.org/pub"
                    release = "release-100"
                    url_host = "http://asia.ensembl.org/{}/Info/Index"
                else:
                    host = "ftp://ftp.ensemblgenomes.org/pub/"
                    release = "release-46"
                    spe_class = e_dict.get("division").split('Ensembl')[1].lower()
                    release += "/" + spe_class
                    url_host = "http://" + spe_class + ".ensembl.org/{}/Info/Index"
                name_path = e_dict.get("scientific_name").replace(" ", "_").lower()
                ftp_path = "/".join([host,
                                     release,
                                     "fasta",
                                     name_path,
                                     "dna"])

                url_path = url_host.format(e_dict.get("scientific_name").replace(" ", "_"))


                self.name2_list.append(e_dict.get("scientific_name"))
                r_dict = {
                    "taxonomy_id": e_dict.get("taxon_id", ""),
                    "assembly": e_dict.get("assembly_name", ""),
                    "assembly_accession": e_dict.get("assembly_accession", ""),
                    "assembly_level": e_dict.get("assembly_level", ""),
                    "size": float(e_dict.get("size", ""))/1000000,
                    "chr_num": e_dict.get("chr_num", ""),
                    "ftp": ftp_path,
                    "url": url_path,
                    "other_info": ";".join(["{}:{}".format(k,e_dict.get(k, "")) for k in ['scientific_name', 'genebuild']])
                }
                results.append(r_dict)
            except Exception as e:
                self.logger.info("ensemble 没有查询到该物种 {}".format(e))
                pass
        return results

    def search_online(self, results):
        if self.option("ensembl") not in ["yes", "true"]:
            results.extend(self.search_ensembl_online())
        if self.option("ncbi") not in ["yes", "true"]:
            results.extend(self.search_ncbi_online())
        return results


    def run(self):
        self.start_listener()
        self.fire("start")
        relations = []
        self.get_search_type()

        self.taxon_db = self.api.api("gene_db.taxon_db")
        self.genome_db = self.api.api("gene_db.genome_db")

        self.taxon_list = []
        if self.option("species") != "":
            taxons = self.taxon_db.search_taxon(self.option("species"))
            self.taxon_list = list(set(taxons.values()))
            self.name_list = list(set(taxons.keys()))

        if self.search_type == "taxon":
            self.taxon_list = [self.option("species")]
            self.name_list = [self.taxon_db.search_taxon_by_id(self.option("species"))]

        self.logger.info("taxon is {}".format(self.taxon_list))
        results = list()

        if len(self.taxon_list) > 0:
            if self.option("ensembl") in ["yes", "true"]:
                ens_results = self.genome_db.search_ensembl(self.taxon_list)

                for e_dict in ens_results:
                    self.logger.info("ensembl {}".format(e_dict))
                    if e_dict.get("division") == "ENSEMBL":
                        host = "ftp://ftp.ensembl.org/pub"
                        release = "release-100"
                        url_host = "http://asia.ensembl.org/{}/Info/Index"
                    else:
                        host = "ftp://ftp.ensemblgenomes.org/pub/"
                        release = "release-46"
                        spe_class = e_dict.get("division").split('Ensembl')[1].lower()
                        release += "/" + spe_class
                        url_host = "http://" + spe_class + ".ensembl.org/{}/Info/Index"
                    name_path = e_dict.get("scientific_name").replace(" ", "_").lower()
                    ftp_path = "/".join([host,
                                         release,
                                         "fasta",
                                         name_path,
                                         "dna"])

                    url_path = url_host.format(e_dict.get("scientific_name").replace(" ", "_"))


                    self.name2_list.append(e_dict.get("scientific_name"))
                    r_dict = {
                        "taxonomy_id": e_dict.get("taxon_id", ""),
                        "assembly": e_dict.get("ensembl_assembly", ""),
                        "assembly_accession": e_dict.get("GCA_000951035.1", ""),
                        "assembly_level": e_dict.get("Full genebuild", ""),
                        "size": e_dict.get("size", ""),
                        "chr_num": e_dict.get("chr_num", ""),
                        "ftp": ftp_path,
                        "url": url_path,
                        "other_info": ";".join(["{}:{}".format(k,e_dict.get(k, "")) for k in ['scientific_name', 'common_name']])
                    }
                    results.append(r_dict)
            if self.option("ncbi") in ["yes", "true"]:
                ncbi_results = self.genome_db.search_ncbi(self.taxon_list)

                for n_dict in ncbi_results:
                    # self.logger.info("ncbi {}".format(n_dict))
                    acc_id = n_dict.get("assembly_accession")
                    if acc_id:
                        url_path = "https://www.ncbi.nlm.nih.gov/genome/?term={}".format(acc_id)

                        ftp_path = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/{}/{}/{}/".format(acc_id[4:7], acc_id[7:10], acc_id[10:13])
                    else:
                        url_path = ""
                        ftp_path = ""

                    r_dict = {
                        "taxonomy_id": n_dict.get("taxid", ""),
                        "assembly": n_dict.get("center", ""),
                        "assembly_accession": n_dict.get("assembly_accession", ""),
                        "assembly_level": n_dict.get("status", ""),
                        "size": n_dict.get("size__mb_", ""),
                        "chr_num": n_dict.get("scaffolds", ""),
                        "ftp": ftp_path,
                        "url": url_path,
                        "other_info": ";".join(["{}:{}".format(k,n_dict.get(k, "")) for k in ['_organism_name', 'center']])
                    }
                    results.append(r_dict)
        # if self.option("local_only") not in ["yes", "true"]:
        results = self.search_online(results)
        print results
        self.set_db(results)
        self.end(results)


    def end(self, results):
        out_file = self.output_dir + '/genome_search_result.xls'
        fields = ["taxonomy_id","assembly","assembly_accession","assembly_level","size","chr_num","ftp","url","other_info"]
        with open(out_file, 'w') as f:
            f.write("\t".join(fields) + "\n")
            for result in results:
                f.write("\t".join([result.get(field, '') for field in fields]) + "\n")
        super(GenomedbSearchWorkflow, self).end()


    def set_db(self, results):
        genome_api = self.api.api("tool_lab.genomedb_search")
        genome_api.add_genome_detail(self.option("main_id"), results)
        table_dict = {
            "column": [
                {"field": "taxonomy_id", "title": "taxonomy_id", "filter": "false", "sort": "false", "type": "string"},
                {"field": "assembly", "title": "assembly", "filter": "false", "sort": "false", "type": "string"},
                {"field": "assembly_accession", "title": "assembly_accession", "filter": "false", "sort": "false", "type": "string"},
                {"field": "assembly_level", "title": "assembly_level", "filter": "false", "sort": "false", "type": "string"},
                {"field": "size", "title": "size", "filter": "false", "sort": "false", "type": "string"},
                {"field": "chr_num", "title": "chr_num", "filter": "false", "sort": "false", "type": "string"},
                {"field": "ftp", "title": "ftp", "filter": "false", "sort": "false", "type": "string"},
                {"field": "url", "title": "url", "filter": "false", "sort": "false", "type": "string"},
                {"field": "other_info", "title": "other_info", "filter": "false", "sort": "false", "type": "string"}
            ],
            "condition": {}
        }
        table_info = json.dumps(table_dict, sort_keys=True, separators=(',', ':'))

        genome_api.update_db_record('genomedb_search',
                                    query_dict={"main_id": ObjectId(self.option("main_id"))},
                                    update_dict={'status': 'end',
                                                 'table_data': table_info})
