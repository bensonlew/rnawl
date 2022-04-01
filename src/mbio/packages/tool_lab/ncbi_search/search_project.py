# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import os
import json
from mbio.packages.ref_rna_v2.functions import tryforgood
from biocluster.config import Config
import requests
import re
import time


class SearchProject(object):
    """
    This script is used to search projects, use "SRP", "ERP", "DRP", "GSE" or keyword.
    project_fields = ["BioProject", "Title", "SampleNumbers", "GseAccession", "SraAccesion", "Summary"]
    accession examples: PRJNA128733, GSE22260, SRP002628
    """

    def __init__(self, accession=None, query=None, start=0, stop=20):
        perl_path = Config().SOFTWARE_DIR + '/program/perl/perls/perl-5.24.0/bin'
        os.environ['PATH'] = perl_path + ":" + os.environ['PATH']
        self.accession = accession
        self.query = query
        self.start = start
        self.stop = stop
        self.program = {
            'esearch': Config().SOFTWARE_DIR + '/bioinfo/biodb/edirect/esearch',
            'efetch': Config().SOFTWARE_DIR + '/bioinfo/biodb/edirect/efetch',
            'esummary': Config().SOFTWARE_DIR + '/bioinfo/biodb/edirect/esummary'
        }

    def run(self):
        if self.accession:
            if self.accession[0:3] in ["SRP", "DRP", "ERP"] or self.accession.startswith("PRJNA") or self.accession.startswith("PRJDB") or self.accession.startswith("PRJEB"):
                result = self.search_sra(self.accession)
                if result:
                    self.parse_sra_result(result, accession=self.accession)
                else:
                    raise Exception("search sra failed!")
            elif self.accession[0:3] in ["GSE"]:
                result = self.search_gds(self.accession)
                if result:
                    self.parse_gds_result(result)
                else:
                    raise Exception("search gds failed!")
            elif self.accession[0:3].startswith("ERP") or self.accession.startswith("PRJEB"):
                self.search_ebi(self.accession)
            else:
                raise Exception("{} not supported right now!".format(self.accession))
        if self.query:
            result = self.search_keyword(self.query, self.start, self.stop)
            self.parse_gds_result(result)

    def search_ebi(self, accession):
        url = "https://www.ebi.ac.uk/ena/data/view/" + accession + "&display=xml"
        study_accession = ""
        project_accession = ""
        title = ""
        name = ""
        try:
            xml = requests.get(url)
            pattern = re.compile('>(.*)<')
            if accession.startswith("ERP"):
                for line in xml:
                    if re.match("<PRIMARY_ID>(.*)</PRIMARY_ID>", line):
                        study_accession = pattern.findall(line)[0]
                    if re.match("<SECONDARY_ID>(.*)</SECONDARY_ID>", line):
                        project_accession = pattern.findall(line)[0]
                    if re.match("<CENTER_PROJECT_NAME>(.*)</CENTER_PROJECT_NAME>", line):
                        name = pattern.findall(line)[0]
                    if re.match("<STUDY_TITLE>(.*)</STUDY_TITLE>", line):
                        title = pattern.findall(line)[0]
            else:
                for line in xml:
                    if re.match("<PRIMARY_ID>(.*)</PRIMARY_ID>", line):
                        project_accession = pattern.findall(line)[0]
                    if re.match("<SECONDARY_ID>(.*)</SECONDARY_ID>", line):
                        study_accession = pattern.findall(line)[0]
                    if re.match("<NAME>(.*)</NAME>", line):
                        name = pattern.findall(line)[0]
                    if re.match("<TITLE>(.*)</TITLE>", line):
                        title = pattern.findall(line)[0]
        except:
            raise Exception("search ebi online failed")
        url = "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=" + accession + "&result=read_run&fields=run_accession"
        try:
            r = requests.get(url)
            n_samples = len(r.text.encode('utf-8').strip().split("\n")) - 1
        except:
            raise Exception("search ebi online failed")
        with open('metadata.txt', 'w') as w:
            header = ["Project", "Study", "Name", "Title", "SampleNumber"]
            w.write("\t".join(header) + "\n")
            w.write("\t".join([project_accession, study_accession, name, title, n_samples]) + "\n")

    def search_gds(self, accession):
        cmd = "{} -db gds -query {} | {}  -format docsum -json".format(self.program['esearch'], accession,
                                                                          self.program['efetch'])
        print(cmd)
        try:
            search_json = os.popen(cmd)
        except:
            print("search gds online failed: {}".format(search_json.readlines()))
            time.sleep(15)
            return False
        try:
            results = json.load(search_json)
            with open('gds_info.json', 'w') as f:
                f.write(json.dumps(results, indent=4))
            return results
        except:
            print("search gds online failed")
            time.sleep(15)
            return False

    def search_sra(self, accession):
        cmd = "{} -db sra -query {} | {}  -format docsum -json".format(self.program['esearch'], accession,
                                                                       self.program['efetch'])
        print(cmd)
        try:
            search_json = os.popen(cmd)
        except:
            print("search sra online failed: {}".format(search_json.readlines()))
            time.sleep(15)
            return False
        try:
            results = json.load(search_json)
            with open('sra_info.json', 'w') as f:
                f.write(json.dumps(results, indent=4))
            return results
        except:
            print("search sra online failed")
            time.sleep(15)
            return False

    @tryforgood
    def search_keyword(self, query, start, stop):
        cmd = "{} -db gds -query \"{}\" | {}  -format docsum -json -start {} -stop {}".format(self.program['esearch'],
                                                                                            query,
                                                                                            self.program['efetch'],
                                                                                            start, stop)
        print(cmd)
        try:
            search_json = os.popen(cmd)
        except:
            print("search gds online failed: {}".format(search_json.readlines()))
            time.sleep(15)
            return False
        try:
            results = json.load(search_json)
            with open('gds_info.json', 'w') as f:
                f.write(json.dumps(results, indent=4))
            return results
        except:
            print("search gds online failed")
            time.sleep(15)
            return False

    def parse_sra_result(self, result, accession=None):
        result_list = list()
        if isinstance(result['DocumentSummarySet']['DocumentSummary'], dict):
            expxml = result['DocumentSummarySet']['DocumentSummary']['ExpXml']
            runs = result['DocumentSummarySet']['DocumentSummary']['Runs']['Run']
            updatedate = result['DocumentSummarySet']['DocumentSummary']['UpdateDate']
            total_runs = expxml['Summary']['Statistics']['total_runs']
            if 'LIBRARY_STRATEGY' in expxml['Library_descriptor']:
                if isinstance(expxml['Library_descriptor']['LIBRARY_STRATEGY'], dict):
                    if 'content' in expxml['Library_descriptor']['LIBRARY_STRATEGY']:
                        library_strategy = expxml['Library_descriptor']['LIBRARY_STRATEGY']['content']
                    else:
                        library_strategy = ''
                else:
                    library_strategy = expxml['Library_descriptor']['LIBRARY_STRATEGY']
            else:
                library_strategy = ''
            if 'LIBRARY_SOURCE' in expxml['Library_descriptor']:
                if isinstance(expxml['Library_descriptor']['LIBRARY_SOURCE'], dict):
                    if 'content' in expxml['Library_descriptor']['LIBRARY_SOURCE']:
                        library_source = expxml['Library_descriptor']['LIBRARY_SOURCE']['content']
                    else:
                        library_source = ''
                else:
                    library_source = expxml['Library_descriptor']['LIBRARY_SOURCE']
            else:
                library_source = ''
            if 'LIBRARY_LAYOUT' in expxml['Library_descriptor']:
                if isinstance(expxml['Library_descriptor']['LIBRARY_LAYOUT'], dict):
                    if 'content' in expxml['Library_descriptor']['LIBRARY_LAYOUT']:
                        library_layout = expxml['Library_descriptor']['LIBRARY_LAYOUT']['content']
                    elif 'PAIRED' in expxml['Library_descriptor']['LIBRARY_LAYOUT']:
                        library_layout = 'PAIRED'
                    elif 'SINGLE' in expxml['Library_descriptor']['LIBRARY_LAYOUT']:
                        library_layout = 'SINGLE'
                    else:
                        library_layout = ''
                else:
                    library_layout = expxml['Library_descriptor']['LIBRARY_LAYOUT']
            else:
                library_layout = ''
            if 'LIBRARY_SELECTION' in expxml['Library_descriptor']:
                if isinstance(expxml['Library_descriptor']['LIBRARY_SELECTION'], dict):
                    if 'content' in expxml['Library_descriptor']['LIBRARY_SELECTION']:
                        library_selection = expxml['Library_descriptor']['LIBRARY_SELECTION']['content']
                    else:
                        library_selection = ''
                else:
                    library_selection = expxml['Library_descriptor']['LIBRARY_SELECTION']
            else:
                library_selection = ''
            if 'LIBRARY_NAME' in expxml['Library_descriptor']:
                if isinstance(expxml['Library_descriptor']['LIBRARY_NAME'], dict):
                    if 'content' in expxml['Library_descriptor']['LIBRARY_NAME']:
                        library_name = expxml['Library_descriptor']['LIBRARY_NAME']['content']
                    else:
                        library_name = ''
                else:
                    library_name = expxml['Library_descriptor']['LIBRARY_NAME']
            else:
                library_name = ''
            data = {
                "scientific_name": expxml['Organism']['ScientificName'] if 'ScientificName' in expxml[
                    'Organism'] else "",
                "taxonomy_id": expxml['Organism']['taxid'] if 'taxid' in expxml['Organism'] else "",
                "platform": expxml['Summary']['Platform']['instrument_model'] if 'instrument_model' in
                                                                                 expxml['Summary']['Platform'] else "",
                "study_acc": expxml['Study']['acc'] if 'acc' in expxml['Study'] else '',
                "study_name": expxml['Study']['name'] if 'name' in expxml['Study'] else '',
                "library_name": library_name,
                "library_strategy": library_strategy,
                "library_source": library_source,
                "library_layout": library_layout,
                "library_selection": library_selection,
                "biosample": expxml['Biosample'] if 'Biosample' in expxml else '',
                "bioproject": expxml['Bioproject'] if 'Bioproject' in expxml else '',
                "submitter": expxml['Submitter']['acc'] if 'Submitter' in expxml else '',
                "update_date": updatedate
            }
            if total_runs and int(total_runs) >= 2:
                for eachrun in runs:
                    acc = eachrun['acc']
                    if self.accession:
                        if self.accession[0:3] in ["SRR", "DRR", "ERR"]:
                            if self.accession != acc:
                                continue
                    sra_ftp = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/" + acc[
                                                                                                    0:3] + "/" + acc[
                                                                                                                 0:6] + "/" + acc[
                                                                                                                              0:9] + "/" + acc + ".sra"
                    url = "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=" + acc + "&result=read_run&fields=run_accession,fastq_ftp"
                    try:
                        r = requests.get(url)
                        fastq_ftp = r.text.encode('utf-8').split("\n")[1].split("\t")[1]
                    except:
                        fastq_ftp = ""
                    data.update({
                        "acc": acc,
                        "total_bases": eachrun['total_bases'],
                        "total_reads": eachrun['total_spots'],
                        "fastq_ftp": fastq_ftp,
                        "sra_ftp": sra_ftp,
                    })
                    run_data = data.copy()
                    result_list.append(run_data)
            elif total_runs and int(total_runs) == 1:
                acc = runs['acc']
                if self.accession:
                    if self.accession[0:3] in ["SRR", "DRR", "ERR"]:
                        if self.accession != acc:
                            return
                sra_ftp = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/" + acc[
                                                                                                0:3] + "/" + acc[
                                                                                                             0:6] + "/" + acc[
                                                                                                                          0:9] + "/" + acc + ".sra"
                url = "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=" + acc + "&result=read_run&fields=run_accession,fastq_ftp"
                try:
                    r = requests.get(url)
                    fastq_ftp = r.text.encode('utf-8').split("\n")[1].split("\t")[1]
                except:
                    fastq_ftp = ""
                data.update({
                    "acc": acc,
                    "total_bases": runs['total_bases'],
                    "total_reads": runs['total_spots'],
                    "fastq_ftp": fastq_ftp,
                    "sra_ftp": sra_ftp,
                })
                run_data = data.copy()
                result_list.append(run_data)
            else:
                return
        elif isinstance(result['DocumentSummarySet']['DocumentSummary'], list):
            for summary in result['DocumentSummarySet']['DocumentSummary']:
                if 'ExpXml' in summary:
                    updatedate = summary['UpdateDate']
                    expxml = summary['ExpXml']
                    if 'LIBRARY_STRATEGY' in expxml['Library_descriptor']:
                        if isinstance(expxml['Library_descriptor']['LIBRARY_STRATEGY'], dict):
                            if 'content' in expxml['Library_descriptor']['LIBRARY_STRATEGY']:
                                library_strategy = expxml['Library_descriptor']['LIBRARY_STRATEGY']['content']
                            else:
                                library_strategy = ''
                        else:
                            library_strategy = expxml['Library_descriptor']['LIBRARY_STRATEGY']
                    else:
                        library_strategy = ''
                    if 'LIBRARY_SOURCE' in expxml['Library_descriptor']:
                        if isinstance(expxml['Library_descriptor']['LIBRARY_SOURCE'], dict):
                            if 'content' in expxml['Library_descriptor']['LIBRARY_SOURCE']:
                                library_source = expxml['Library_descriptor']['LIBRARY_SOURCE']['content']
                            else:
                                library_source = ''
                        else:
                            library_source = expxml['Library_descriptor']['LIBRARY_SOURCE']
                    else:
                        library_source = ''
                    if 'LIBRARY_LAYOUT' in expxml['Library_descriptor']:
                        if isinstance(expxml['Library_descriptor']['LIBRARY_LAYOUT'], dict):
                            if 'content' in expxml['Library_descriptor']['LIBRARY_LAYOUT']:
                                library_layout = expxml['Library_descriptor']['LIBRARY_LAYOUT']['content']
                            elif 'PAIRED' in expxml['Library_descriptor']['LIBRARY_LAYOUT']:
                                library_layout = 'PAIRED'
                            elif 'SINGLE' in expxml['Library_descriptor']['LIBRARY_LAYOUT']:
                                library_layout = 'SINGLE'
                            else:
                                library_layout = ''
                        else:
                            library_layout = expxml['Library_descriptor']['LIBRARY_LAYOUT']
                    else:
                        library_layout = ''
                    if 'LIBRARY_SELECTION' in expxml['Library_descriptor']:
                        if isinstance(expxml['Library_descriptor']['LIBRARY_SELECTION'], dict):
                            if 'content' in expxml['Library_descriptor']['LIBRARY_SELECTION']:
                                library_selection = expxml['Library_descriptor']['LIBRARY_SELECTION']['content']
                            else:
                                library_selection = ''
                        else:
                            library_selection = expxml['Library_descriptor']['LIBRARY_SELECTION']
                    else:
                        library_selection = ''
                    if 'LIBRARY_NAME' in expxml['Library_descriptor']:
                        if isinstance(expxml['Library_descriptor']['LIBRARY_NAME'], dict):
                            if 'content' in expxml['Library_descriptor']['LIBRARY_NAME']:
                                library_name = expxml['Library_descriptor']['LIBRARY_NAME']['content']
                            else:
                                library_name = ''
                        else:
                            library_name = expxml['Library_descriptor']['LIBRARY_NAME']
                    else:
                        library_name = ''
                    data = {
                        "scientific_name": expxml['Organism']['ScientificName'] if 'ScientificName' in expxml[
                            'Organism'] else '',
                        "taxonomy_id": expxml['Organism']['taxid'] if 'taxid' in expxml['Organism'] else '',
                        "platform": expxml['Summary']['Platform']['instrument_model'] if 'instrument_model' in
                                                                                         expxml['Summary'][
                                                                                             'Platform'] else '',
                        "study_acc": expxml['Study']['acc'] if 'acc' in expxml['Study'] else '',
                        "study_name": expxml['Study']['name'] if 'name' in expxml['Study'] else '',
                        "library_name": library_name,
                        "library_strategy": library_strategy,
                        "library_source": library_source,
                        "library_layout": library_layout,
                        "library_selection": library_selection,
                        "biosample": expxml['Biosample'] if 'Biosample' in expxml else '',
                        "bioproject": expxml['Bioproject'] if 'Bioproject' in expxml else '',
                        "submitter": expxml['Submitter']['acc'] if 'Submitter' in expxml else '',
                        "update_date": updatedate
                    }
                    total_runs = expxml['Summary']['Statistics']['total_runs']
                    if total_runs and int(total_runs) >= 2:
                        for eachrun in summary['Runs']['Run']:
                            acc = eachrun['acc']
                            if self.accession:
                                if self.accession[0:3] in ["SRR", "DRR", "ERR"]:
                                    if self.accession != acc:
                                        continue
                            sra_ftp = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/" + acc[
                                                                                                            0:3] + "/" + acc[
                                                                                                                         0:6] + "/" + acc[
                                                                                                                                      0:9] + "/" + acc + ".sra"
                            url = "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=" + acc + "&result=read_run&fields=run_accession,fastq_ftp"
                            try:
                                r = requests.get(url)
                                fastq_ftp = r.text.encode('utf-8').split("\n")[1].split("\t")[1]
                            except:
                                fastq_ftp = ""
                            data.update({
                                "acc": acc,
                                "total_bases": eachrun['total_bases'],
                                "total_reads": eachrun['total_spots'],
                                "fastq_ftp": fastq_ftp,
                                "sra_ftp": sra_ftp,
                            })
                            run_data = data.copy()
                            result_list.append(run_data)
                    elif total_runs and int(total_runs) == 1:
                        run = summary['Runs']['Run']
                        acc = run['acc']
                        if self.accession:
                            if self.accession[0:3] in ["SRR", "DRR", "ERR"]:
                                if self.accession != acc:
                                    continue
                        sra_ftp = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/" + acc[
                                                                                                        0:3] + "/" + acc[
                                                                                                                     0:6] + "/" + acc[
                                                                                                                                  0:9] + "/" + acc + ".sra"
                        url = "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=" + acc + "&result=read_run&fields=run_accession,fastq_ftp"
                        try:
                            r = requests.get(url)
                            fastq_ftp = r.text.encode('utf-8').split("\n")[1].split("\t")[1]
                        except:
                            fastq_ftp = ""
                        data.update({
                            "acc": acc,
                            "total_bases": run['total_bases'],
                            "total_reads": run['total_spots'],
                            "fastq_ftp": fastq_ftp,
                            "sra_ftp": sra_ftp,
                        })
                        run_data = data.copy()
                        result_list.append(run_data)
                    else:
                        continue
        if accession:
            output = '{}_metadata.txt'.format(accession)
        else:
            output = 'metadata.txt'
        with open(output, 'w') as w:
            header = ["Accession", "UpdateDate", "TotalReads", "TotalBases", "FastqURL", "SraURL", "ScientificName", "Taxid",
                      "Platform", "LibraryName", "LibraryStrategy", "LibrarySource", "LibraryLayout",
                      "LibrarySelection", "StudyAcc", "StudyName", "Bioproject", "Biosample", "Submitter"]
            w.write("\t".join(header) + "\n")
            for sample in result_list:
                content = [sample["acc"], sample["update_date"], sample["total_reads"], sample["total_bases"], sample["fastq_ftp"],
                           sample["sra_ftp"], sample["scientific_name"], sample["taxonomy_id"], sample["platform"],
                           sample["library_name"], sample["library_strategy"], sample["library_source"],
                           sample["library_layout"], sample["library_selection"], sample["study_acc"],
                           sample["study_name"], sample["bioproject"], sample["biosample"], sample["submitter"]]
                w.write("\t".join(content) + "\n")

    def parse_gds_result(self, result):
        result_list = list()
        if isinstance(result['DocumentSummarySet']['DocumentSummary'], dict):
            summary = result['DocumentSummarySet']['DocumentSummary']
            data = {
                "BioProject": summary['BioProject'],
                "title": summary['title'],
                "summary": summary['summary'],
                "n_samples": summary['n_samples'],
                "gse_accession": "GSE" + str(summary['GSE']),
            }
            if "TargetObject" in summary['ExtRelations']['ExtRelation']:
                data.update({"sra_accession": summary['ExtRelations']['ExtRelation']['TargetObject']})
            result_list.append(data)
        elif isinstance(result['DocumentSummarySet']['DocumentSummary'], list):
            for summary in result['DocumentSummarySet']['DocumentSummary']:
                if summary['entryType'] != "GSE":
                    continue
                data = {
                    "BioProject": summary['BioProject'],
                    "title": summary['title'],
                    "summary": summary['summary'],
                    "n_samples": summary['n_samples'],
                    "gse_accession": "GSE" + str(summary['GSE']),
                }
                if 'ExtRelation' in summary['ExtRelations']:
                    if "TargetObject" in summary['ExtRelations']['ExtRelation']:
                        data.update({"sra_accession": summary['ExtRelations']['ExtRelation']['TargetObject']})
                else:
                    data.update({"sra_accession": ""})
                result_list.append(data)
        if self.accession:
            result = '{}_metadata.txt'.format(self.accession)
        else:
            result = 'metadata.txt'
        with open(result, 'w') as w:
            header = ["BioProject", "Title", "SampleNumber", "GseAccession", "SraAccesion", "Summary"]
            w.write("\t".join(header) + "\n")
            for doc in result_list:
                content = [doc["BioProject"], doc["title"], doc["n_samples"], doc["gse_accession"], doc["sra_accession"], doc["summary"]]
                w.write("\t".join(content) + "\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='This script is used to search projects, use "SRP", "ERP", "DRP", "GSE" or keyword.')
    parser.add_argument('-accession', type=str, default=None, help='accession examples: PRJNA128733, GSE22260, SRP002628')
    parser.add_argument('-keyword', type=str, default=None, help='Query string used to search.')
    parser.add_argument('-start', type=int, default=0, help='First record to fetch, used when -keyword is set.')
    parser.add_argument('-stop', type=int, default=20, help='Last record to fetch, used when -keyword is set.')

    args = parser.parse_args()

    if args.accession and args.keyword:
        raise Exception("accession and keyword cannot be searched at the same time.")
    if not (args.accession or args.keyword):
        raise Exception("accession and keyword cannot be empty at the same time.")
    search = SearchProject(args.accession, args.keyword, args.start, args.stop)
    search.run()
